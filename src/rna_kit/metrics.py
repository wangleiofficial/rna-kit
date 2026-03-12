from __future__ import annotations

import copy
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from Bio.PDB import PDBIO, Superimposer

from .alignment import StructureAlignment, infer_structure_alignment
from .exceptions import MetricCalculationError, SequenceMismatchError, ToolNotAvailableError
from .mc_annotate import MCAnnotateRunner
from .secondary_structure import (
    SecondaryStructureComparisonResult,
    SecondaryStructureResult,
    calculate_secondary_structure_for_structure,
    compare_secondary_structures,
)
from .structures import PDBStructure


@dataclass(frozen=True)
class RMSDResult:
    rmsd: float
    pvalue: float


@dataclass(frozen=True)
class InteractionNetworkResult:
    rmsd: float
    deformation_index: float
    inf_all: float
    inf_wc: float
    inf_nwc: float
    inf_stack: float


@dataclass(frozen=True)
class LDDTResult:
    lddt: float
    evaluated_atoms: int
    evaluated_pairs: int
    inclusion_radius: float
    per_residue: tuple["ResidueAssessment", ...] | None = None


@dataclass(frozen=True)
class ResidueAssessment:
    native_chain: str
    native_pos: int
    native_nt: str
    prediction_chain: str
    prediction_pos: int
    prediction_nt: str
    matched_atoms: int
    scored_atoms: int
    lddt: float | None
    local_rmsd: float | None
    mean_absolute_error: float | None
    max_absolute_error: float | None


@dataclass(frozen=True)
class AssessmentResult:
    rmsd: float
    pvalue: float
    deformation_index: float
    inf_all: float
    inf_wc: float
    inf_nwc: float
    inf_stack: float
    lddt: float
    lddt_evaluated_atoms: int
    lddt_evaluated_pairs: int
    secondary_structure_precision: float | None = None
    secondary_structure_recall: float | None = None
    secondary_structure_f1: float | None = None
    secondary_structure_jaccard: float | None = None
    secondary_structure: SecondaryStructureComparisonResult | None = None
    per_residue: tuple[ResidueAssessment, ...] | None = None


@dataclass(frozen=True)
class PreparedStructurePair:
    native: PDBStructure
    prediction: PDBStructure
    native_index: str | None
    prediction_index: str | None
    alignment: StructureAlignment | None
    used_sidecar_index: bool


def erf(z: float) -> float:
    t = 1.0 / (1.0 + 0.5 * abs(z))
    ans = 1 - t * math.exp(
        -z * z
        - 1.26551223
        + t
        * (
            1.00002368
            + t
            * (
                0.37409196
                + t
                * (
                    0.09678418
                    + t
                    * (
                        -0.18628806
                        + t
                        * (
                            0.27886807
                            + t
                            * (
                                -1.13520398
                                + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277))
                            )
                        )
                    )
                )
            )
        )
    )
    return ans if z >= 0.0 else -ans


class PDBComparer:
    BACKBONE_ATOMS = ["C1'", "C2'", "C3'", "C4'", "C5'", "O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "P"]
    HEAVY_ATOMS = ["C2", "C4", "C5", "C6", "C8", "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"]
    ALL_ATOMS = BACKBONE_ATOMS + HEAVY_ATOMS
    LDDT_THRESHOLDS = (0.5, 1.0, 2.0, 4.0)

    def rmsd(self, src_struct: PDBStructure, trg_struct: PDBStructure, fit_pdb: str | Path | None = None) -> float:
        src_atoms, trg_atoms = self._get_atoms_struct(self.ALL_ATOMS, src_struct.res_sequence(), trg_struct.res_sequence())
        superimposer = Superimposer()
        superimposer.set_atoms(src_atoms, trg_atoms)

        fit_structure = copy.deepcopy(trg_struct.struct)
        superimposer.apply(fit_structure.get_atoms())

        if fit_pdb is not None:
            io = PDBIO()
            io.set_structure(fit_structure)
            io.save(str(fit_pdb))

        return float(superimposer.rms)

    def pvalue(self, rmsd_value: float, residue_count: int, param: str) -> float:
        if param == "+":
            a, b = 5.1, 15.8
        elif param == "-":
            a, b = 6.4, 12.7
        else:
            raise ValueError(f"Wrong p-value parameter '{param}'. Expected '+' or '-'.")
        expected_rmsd = a * (residue_count**0.41) - b
        z_score = (rmsd_value - expected_rmsd) / 1.8
        return (1.0 + erf(z_score / (2**0.5))) / 2.0

    def inf(
        self,
        src_struct: PDBStructure,
        trg_struct: PDBStructure,
        interaction_type: str = "ALL",
        annotator: MCAnnotateRunner | None = None,
    ) -> float:
        annotator = annotator or MCAnnotateRunner()
        src_interactions = self._select_interactions(annotator.indexed_interactions(src_struct), interaction_type)
        trg_interactions = self._select_interactions(annotator.indexed_interactions(trg_struct), interaction_type)

        true_positives = len(set(src_interactions) & set(trg_interactions))
        false_negatives = len(set(src_interactions) - set(trg_interactions))
        false_positives = len(set(trg_interactions) - set(src_interactions))

        if true_positives == 0 and (false_positives == 0 or false_negatives == 0):
            return -1.0

        precision = true_positives / float(true_positives + false_positives)
        sensitivity = true_positives / float(true_positives + false_negatives)
        return math.sqrt(precision * sensitivity)

    def lddt(
        self,
        reference_struct: PDBStructure,
        model_struct: PDBStructure,
        inclusion_radius: float = 15.0,
        thresholds: tuple[float, ...] = LDDT_THRESHOLDS,
        include_per_residue: bool = False,
    ) -> LDDTResult:
        residue_pairs = self._get_residue_atom_pairs(reference_struct, model_struct)
        reference_atoms = [atom for pair in residue_pairs for atom in pair[2]]
        model_atoms = [atom for pair in residue_pairs for atom in pair[3]]

        if not reference_atoms:
            raise MetricCalculationError("No matched atoms available for lDDT calculation.")

        atom_scores, scored_atoms, total_pairs, total_score = self._lddt_atom_scores(
            reference_atoms,
            model_atoms,
            inclusion_radius=inclusion_radius,
            thresholds=thresholds,
        )
        if scored_atoms == 0:
            raise MetricCalculationError("No eligible atom pairs within the lDDT inclusion radius.")

        return LDDTResult(
            lddt=total_score / scored_atoms,
            evaluated_atoms=scored_atoms,
            evaluated_pairs=total_pairs,
            inclusion_radius=inclusion_radius,
            per_residue=self.per_residue_report(
                reference_struct,
                model_struct,
                inclusion_radius=inclusion_radius,
                thresholds=thresholds,
                residue_pairs=residue_pairs,
                atom_scores=atom_scores,
            )
            if include_per_residue
            else None,
        )

    def per_residue_report(
        self,
        reference_struct: PDBStructure,
        model_struct: PDBStructure,
        inclusion_radius: float = 15.0,
        thresholds: tuple[float, ...] = LDDT_THRESHOLDS,
        residue_pairs=None,
        atom_scores: list[float | None] | None = None,
    ) -> tuple[ResidueAssessment, ...]:
        residue_pairs = residue_pairs or self._get_residue_atom_pairs(reference_struct, model_struct)
        reference_atoms = [atom for pair in residue_pairs for atom in pair[2]]
        model_atoms = [atom for pair in residue_pairs for atom in pair[3]]
        if not reference_atoms:
            raise MetricCalculationError("No matched atoms available for per-residue reporting.")

        if atom_scores is None:
            atom_scores, _, _, _ = self._lddt_atom_scores(
                reference_atoms,
                model_atoms,
                inclusion_radius=inclusion_radius,
                thresholds=thresholds,
            )

        rotation, translation = self._fit_rotran(reference_atoms, model_atoms)
        reference_coords = [atom.get_coord() for atom in reference_atoms]
        fitted_model_coords = [atom.get_coord().dot(rotation) + translation for atom in model_atoms]

        reports: list[ResidueAssessment] = []
        atom_offset = 0
        for reference_record, model_record, residue_reference_atoms, residue_model_atoms in residue_pairs:
            atom_count = len(residue_reference_atoms)
            residue_reference_coords = reference_coords[atom_offset : atom_offset + atom_count]
            residue_model_coords = fitted_model_coords[atom_offset : atom_offset + atom_count]
            residue_atom_scores = [
                score
                for score in atom_scores[atom_offset : atom_offset + atom_count]
                if score is not None
            ]

            distances = [
                math.dist(reference_coord, model_coord)
                for reference_coord, model_coord in zip(residue_reference_coords, residue_model_coords)
            ]
            reports.append(
                ResidueAssessment(
                    native_chain=reference_record.chain,
                    native_pos=reference_record.pos,
                    native_nt=reference_record.nt,
                    prediction_chain=model_record.chain,
                    prediction_pos=model_record.pos,
                    prediction_nt=model_record.nt,
                    matched_atoms=atom_count,
                    scored_atoms=len(residue_atom_scores),
                    lddt=(sum(residue_atom_scores) / len(residue_atom_scores)) if residue_atom_scores else None,
                    local_rmsd=(sum(distance * distance for distance in distances) / atom_count) ** 0.5
                    if atom_count
                    else None,
                    mean_absolute_error=(sum(distances) / atom_count) if atom_count else None,
                    max_absolute_error=max(distances) if distances else None,
                )
            )
            atom_offset += atom_count

        return tuple(reports)

    def mcq(self, model_file: str | Path, target_file: str | Path, jar_path: str | Path | None = None) -> float:
        jar = Path(jar_path) if jar_path else _default_jar_path("mcq.ws.client-0.0.1-SNAPSHOT-jar-with-dependencies.jar")
        return self._run_java_metric(
            [
                "-cp",
                str(jar),
                "pl.poznan.put.mcq.ws.client.Global",
                "-m",
                str(model_file),
                "-t",
                str(target_file),
            ],
            parser=lambda stdout: float(stdout.strip()),
            tool_name="MCQ",
        )

    def gdt(self, reference_file: str | Path, model_file: str | Path, jar_path: str | Path | None = None) -> float:
        jar = Path(jar_path) if jar_path else _default_jar_path("gdt.jar")
        return self._run_java_metric(
            ["-jar", str(jar), str(model_file), str(reference_file)],
            parser=_parse_gdt_output,
            tool_name="GDT",
        )

    def _run_java_metric(self, args: list[str], parser, tool_name: str) -> float:
        if shutil.which("java") is None:
            raise ToolNotAvailableError(f"{tool_name} requires a Java runtime.")
        try:
            result = subprocess.run(["java", *args], check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:
            raise ToolNotAvailableError(
                f"{tool_name} failed: {exc.stderr.strip() or exc.stdout.strip()}"
            ) from exc
        return parser(result.stdout)

    def _get_atoms_residue(self, atom_list: list[str], src_residue, trg_residue) -> tuple[list[object], list[object]]:
        src_atoms: list[object] = []
        trg_atoms: list[object] = []
        src_candidates = [atom for atom in src_residue if atom.get_name() in atom_list]
        trg_candidates = [atom for atom in trg_residue if atom.get_name() in atom_list]

        for src_atom in src_candidates:
            src_name = src_atom.get_full_id()[4][0]
            for trg_atom in trg_candidates:
                trg_name = trg_atom.get_full_id()[4][0]
                if src_name == trg_name:
                    src_atoms.append(src_atom)
                    trg_atoms.append(trg_atom)
                    break
        return src_atoms, trg_atoms

    def _get_atoms_struct(
        self, atom_list: list[str], src_residues: list[object], trg_residues: list[object]
    ) -> tuple[list[object], list[object]]:
        if len(src_residues) != len(trg_residues):
            raise ValueError("Different number of residues.")
        src_atoms: list[object] = []
        trg_atoms: list[object] = []
        for src_residue, trg_residue in zip(src_residues, trg_residues):
            src_chunk, trg_chunk = self._get_atoms_residue(atom_list, src_residue, trg_residue)
            src_atoms.extend(src_chunk)
            trg_atoms.extend(trg_chunk)
        return src_atoms, trg_atoms

    def _fit_rotran(self, reference_atoms: list[object], model_atoms: list[object]):
        superimposer = Superimposer()
        superimposer.set_atoms(reference_atoms, model_atoms)
        return superimposer.rotran

    def _get_residue_atom_pairs(
        self,
        reference_struct: PDBStructure,
        model_struct: PDBStructure,
    ) -> list[tuple[object, object, list[object], list[object]]]:
        if len(reference_struct.res_seq) != len(model_struct.res_seq):
            raise ValueError("Different number of residues.")

        residue_pairs: list[tuple[object, object, list[object], list[object]]] = []
        for (_, reference_record), (_, model_record) in zip(
            reference_struct.selected_records(),
            model_struct.selected_records(),
        ):
            reference_atoms, model_atoms = self._get_atoms_residue(
                self.ALL_ATOMS,
                reference_record.residue,
                model_record.residue,
            )
            residue_pairs.append((reference_record, model_record, reference_atoms, model_atoms))
        return residue_pairs

    def _lddt_atom_scores(
        self,
        reference_atoms: list[object],
        model_atoms: list[object],
        inclusion_radius: float,
        thresholds: tuple[float, ...],
    ) -> tuple[list[float | None], int, int, float]:
        reference_coords = [tuple(float(value) for value in atom.get_coord()) for atom in reference_atoms]
        model_coords = [tuple(float(value) for value in atom.get_coord()) for atom in model_atoms]
        atom_scores: list[float | None] = [None] * len(reference_atoms)
        scored_atoms = 0
        total_pairs = 0
        total_score = 0.0
        threshold_count = float(len(thresholds))

        for atom_index, reference_coord in enumerate(reference_coords):
            atom_pairs = 0
            atom_score = 0.0

            for neighbor_index, neighbor_coord in enumerate(reference_coords):
                if atom_index == neighbor_index:
                    continue

                reference_distance = math.dist(reference_coord, neighbor_coord)
                if reference_distance == 0.0 or reference_distance > inclusion_radius:
                    continue

                model_distance = math.dist(model_coords[atom_index], model_coords[neighbor_index])
                delta = abs(reference_distance - model_distance)
                atom_score += sum(delta <= threshold for threshold in thresholds) / threshold_count
                atom_pairs += 1

            if atom_pairs == 0:
                continue

            atom_scores[atom_index] = atom_score / atom_pairs
            scored_atoms += 1
            total_pairs += atom_pairs
            total_score += atom_scores[atom_index]

        return atom_scores, scored_atoms, total_pairs, total_score

    def _select_interactions(
        self, interactions: list[tuple[str, int, int, str]], interaction_type: str
    ) -> list[tuple[str, int, int, str]]:
        if interaction_type == "ALL":
            return interactions
        if interaction_type == "PAIR":
            return [item for item in interactions if item[0] in {"PAIR_2D", "PAIR_3D"}]
        if interaction_type in {"PAIR_2D", "PAIR_3D", "STACK"}:
            return [item for item in interactions if item[0] == interaction_type]
        raise ValueError(
            f"Wrong interaction type '{interaction_type}'. Expected 'ALL', 'PAIR', 'PAIR_2D', 'PAIR_3D' or 'STACK'."
        )


def calculate_rmsd(
    native_file: str | Path,
    native_index: str | Path | None,
    prediction_file: str | Path,
    prediction_index: str | Path | None,
    pvalue_mode: str = "-",
) -> RMSDResult:
    prepared = prepare_structure_pair(native_file, native_index, prediction_file, prediction_index)
    native, prediction = prepared.native, prepared.prediction
    comparer = PDBComparer()
    rmsd_value = comparer.rmsd(prediction, native)
    return RMSDResult(
        rmsd=rmsd_value,
        pvalue=comparer.pvalue(rmsd_value, len(prediction.raw_sequence()), pvalue_mode),
    )


def calculate_interaction_network_fidelity(
    native_file: str | Path,
    native_index: str | Path | None,
    prediction_file: str | Path,
    prediction_index: str | Path | None,
    annotator: MCAnnotateRunner | None = None,
) -> InteractionNetworkResult:
    prepared = prepare_structure_pair(native_file, native_index, prediction_file, prediction_index)
    native, prediction = prepared.native, prepared.prediction
    comparer = PDBComparer()
    rmsd_value = comparer.rmsd(prediction, native)
    inf_all = comparer.inf(prediction, native, interaction_type="ALL", annotator=annotator)
    return InteractionNetworkResult(
        rmsd=rmsd_value,
        deformation_index=rmsd_value / inf_all,
        inf_all=inf_all,
        inf_wc=comparer.inf(prediction, native, interaction_type="PAIR_2D", annotator=annotator),
        inf_nwc=comparer.inf(prediction, native, interaction_type="PAIR_3D", annotator=annotator),
        inf_stack=comparer.inf(prediction, native, interaction_type="STACK", annotator=annotator),
    )


def calculate_lddt(
    native_file: str | Path,
    native_index: str | Path | None,
    prediction_file: str | Path,
    prediction_index: str | Path | None,
    inclusion_radius: float = 15.0,
    include_per_residue: bool = False,
) -> LDDTResult:
    prepared = prepare_structure_pair(native_file, native_index, prediction_file, prediction_index)
    native, prediction = prepared.native, prepared.prediction
    comparer = PDBComparer()
    return comparer.lddt(
        native,
        prediction,
        inclusion_radius=inclusion_radius,
        include_per_residue=include_per_residue,
    )


def calculate_assessment(
    native_file: str | Path,
    native_index: str | Path | None,
    prediction_file: str | Path,
    prediction_index: str | Path | None,
    pvalue_mode: str = "-",
    annotator: MCAnnotateRunner | None = None,
    inclusion_radius: float = 15.0,
    include_per_residue: bool = False,
    include_secondary_structure: bool = False,
    secondary_structure_runner: MCAnnotateRunner | None = None,
) -> AssessmentResult:
    prepared = prepare_structure_pair(native_file, native_index, prediction_file, prediction_index)
    return calculate_assessment_from_prepared(
        prepared,
        pvalue_mode=pvalue_mode,
        annotator=annotator,
        inclusion_radius=inclusion_radius,
        include_per_residue=include_per_residue,
        include_secondary_structure=include_secondary_structure,
        secondary_structure_runner=secondary_structure_runner,
    )


def calculate_assessment_from_prepared(
    prepared: PreparedStructurePair,
    pvalue_mode: str,
    annotator: MCAnnotateRunner | None,
    inclusion_radius: float,
    include_per_residue: bool = False,
    include_secondary_structure: bool = False,
    secondary_structure_runner: MCAnnotateRunner | None = None,
) -> AssessmentResult:
    native, prediction = prepared.native, prepared.prediction
    comparer = PDBComparer()

    rmsd_value = comparer.rmsd(prediction, native)
    pvalue = comparer.pvalue(rmsd_value, len(prediction.raw_sequence()), pvalue_mode)
    inf_all = comparer.inf(prediction, native, interaction_type="ALL", annotator=annotator)
    lddt_result = comparer.lddt(
        native,
        prediction,
        inclusion_radius=inclusion_radius,
        include_per_residue=include_per_residue,
    )
    secondary_structure_result = (
        compare_secondary_structures(
            native,
            prediction,
            runner=secondary_structure_runner or annotator,
        )
        if include_secondary_structure
        else None
    )

    return AssessmentResult(
        rmsd=rmsd_value,
        pvalue=pvalue,
        deformation_index=rmsd_value / inf_all,
        inf_all=inf_all,
        inf_wc=comparer.inf(prediction, native, interaction_type="PAIR_2D", annotator=annotator),
        inf_nwc=comparer.inf(prediction, native, interaction_type="PAIR_3D", annotator=annotator),
        inf_stack=comparer.inf(prediction, native, interaction_type="STACK", annotator=annotator),
        lddt=lddt_result.lddt,
        lddt_evaluated_atoms=lddt_result.evaluated_atoms,
        lddt_evaluated_pairs=lddt_result.evaluated_pairs,
        secondary_structure_precision=(
            secondary_structure_result.precision if secondary_structure_result is not None else None
        ),
        secondary_structure_recall=(
            secondary_structure_result.recall if secondary_structure_result is not None else None
        ),
        secondary_structure_f1=secondary_structure_result.f1 if secondary_structure_result is not None else None,
        secondary_structure_jaccard=(
            secondary_structure_result.jaccard if secondary_structure_result is not None else None
        ),
        secondary_structure=secondary_structure_result,
        per_residue=lddt_result.per_residue,
    )


def calculate_secondary_structure(
    pdb_file: str | Path,
    index_name: str | Path | None,
    runner: MCAnnotateRunner | None = None,
) -> SecondaryStructureResult:
    resolved_index, _ = _resolve_index_path(pdb_file, index_name, allow_sidecar=False)
    structure = PDBStructure.from_file(pdb_file, index_name=resolved_index)
    return calculate_secondary_structure_for_structure(structure, runner=runner)


def calculate_secondary_structure_comparison(
    native_file: str | Path,
    native_index: str | Path | None,
    prediction_file: str | Path,
    prediction_index: str | Path | None,
    runner: MCAnnotateRunner | None = None,
) -> SecondaryStructureComparisonResult:
    prepared = prepare_structure_pair(
        native_file,
        native_index,
        prediction_file,
        prediction_index,
        resolve_sidecar_indices=False,
    )
    return compare_secondary_structures(prepared.native, prepared.prediction, runner=runner)


def prepare_structure_pair(
    native_file: str | Path,
    native_index: str | Path | None,
    prediction_file: str | Path,
    prediction_index: str | Path | None,
    resolve_sidecar_indices: bool = True,
) -> PreparedStructurePair:
    resolved_native_index, native_sidecar = _resolve_index_path(
        native_file,
        native_index,
        allow_sidecar=resolve_sidecar_indices,
    )
    resolved_prediction_index, prediction_sidecar = _resolve_index_path(
        prediction_file,
        prediction_index,
        allow_sidecar=resolve_sidecar_indices,
    )
    native = PDBStructure.from_file(native_file, index_name=resolved_native_index)
    prediction = PDBStructure.from_file(prediction_file, index_name=resolved_prediction_index)

    if prediction.raw_sequence() == native.raw_sequence():
        return PreparedStructurePair(
            native=native,
            prediction=prediction,
            native_index=native.index_spec(),
            prediction_index=prediction.index_spec(),
            alignment=None,
            used_sidecar_index=native_sidecar or prediction_sidecar,
        )

    if resolved_native_index is not None and resolved_prediction_index is not None:
        raise SequenceMismatchError(
            "Result sequence does not match the reference sequence: "
            f"reference='{native.raw_sequence()}', prediction='{prediction.raw_sequence()}'."
        )

    alignment = infer_structure_alignment(native, prediction)
    matched_native = native.with_selected_indices(list(alignment.native_indices))
    matched_prediction = prediction.with_selected_indices(list(alignment.prediction_indices))
    if matched_prediction.raw_sequence() != matched_native.raw_sequence():
        raise SequenceMismatchError(
            "The inferred residue mapping still produced mismatched sequences: "
            f"reference='{matched_native.raw_sequence()}', prediction='{matched_prediction.raw_sequence()}'."
        )

    return PreparedStructurePair(
        native=matched_native,
        prediction=matched_prediction,
        native_index=alignment.native_index_spec,
        prediction_index=alignment.prediction_index_spec,
        alignment=alignment,
        used_sidecar_index=native_sidecar or prediction_sidecar,
    )


def _default_jar_path(filename: str) -> Path:
    return Path(__file__).resolve().parents[2] / "third_party" / "lib" / filename


def _resolve_index_path(
    pdb_file: str | Path,
    index_name: str | Path | None,
    allow_sidecar: bool = True,
) -> tuple[Path | None, bool]:
    if index_name is not None:
        return Path(index_name), False

    if not allow_sidecar:
        return None, False

    sidecar = Path(pdb_file).with_suffix(".index")
    if sidecar.exists():
        return sidecar, True
    return None, False


def _parse_gdt_output(stdout: str) -> float:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if len(lines) < 2:
        raise ToolNotAvailableError("GDT output was shorter than expected.")
    value = lines[1].split(",")[-1]
    if value == "NaN":
        return 0.0
    return float(value)
