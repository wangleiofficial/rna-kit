from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

from .arena import ArenaRunner
from .alignment import infer_structure_alignment
from .exceptions import ManifestFormatError, RNAAssessmentError
from .mc_annotate import MCAnnotateRunner, clone_with_annotation_overrides, existing_annotation_path
from .molprobity import MolProbityRunner
from .metrics import (
    AssessmentResult,
    PreparedStructurePair,
    calculate_assessment_from_prepared,
    prepared_structure_pair_context,
)


@dataclass(frozen=True)
class BenchmarkJob:
    prediction: str
    native: str | None = None
    label: str | None = None
    native_index: str | None = None
    prediction_index: str | None = None
    native_sequence_hint: str | None = None
    prediction_sequence_hint: str | None = None
    native_annotation: str | None = None
    prediction_annotation: str | None = None


@dataclass(frozen=True)
class ChainMappingResult:
    native_chain: str
    prediction_chain: str
    matched_residues: int


@dataclass(frozen=True)
class BenchmarkEntry:
    native: str
    prediction: str
    label: str | None
    status: str
    native_index: str | None
    prediction_index: str | None
    matched_residues: int | None
    chain_mappings: tuple[ChainMappingResult, ...]
    used_sidecar_index: bool = False
    used_sequence_hints: bool = False
    used_normalized_inputs: bool = False
    used_repaired_inputs: bool = False
    metrics: AssessmentResult | None = None
    error: str | None = None


@dataclass(frozen=True)
class BenchmarkResult:
    reference: str | None
    total_predictions: int
    succeeded: int
    failed: int
    entries: tuple[BenchmarkEntry, ...]


def describe_prepared_pair(
    prepared: PreparedStructurePair,
    prediction_file: str | Path,
    native_file: str | Path,
    label: str | None = None,
) -> BenchmarkEntry:
    alignment = prepared.alignment or infer_structure_alignment(prepared.native, prepared.prediction)
    return BenchmarkEntry(
        native=str(native_file),
        prediction=str(prediction_file),
        label=label,
        status="ready",
        native_index=prepared.native_index,
        prediction_index=prepared.prediction_index,
        matched_residues=len(prepared.native.res_seq),
        chain_mappings=tuple(
            ChainMappingResult(
                native_chain=chain_alignment.native_chain,
                prediction_chain=chain_alignment.prediction_chain,
                matched_residues=chain_alignment.matched_residues,
            )
            for chain_alignment in (alignment.chain_alignments if alignment is not None else ())
        ),
        used_sidecar_index=prepared.used_sidecar_index,
        used_sequence_hints=prepared.used_sequence_hints,
        used_normalized_inputs=prepared.used_normalized_inputs,
        used_repaired_inputs=prepared.used_repaired_inputs,
        metrics=None,
        error=None,
    )


def build_benchmark_jobs(
    native_file: str | Path,
    predictions: list[str | Path],
    native_index: str | Path | None = None,
    prediction_indices: dict[str, str | Path | None] | None = None,
    native_sequence_hint: str | Path | None = None,
    prediction_sequence_hints: dict[str, str | Path | None] | None = None,
) -> list[BenchmarkJob]:
    return [
        BenchmarkJob(
            native=str(native_file),
            prediction=str(prediction),
            native_index=None if native_index is None else str(native_index),
            prediction_index=(
                None
                if prediction_indices is None or prediction_indices.get(str(prediction)) is None
                else str(prediction_indices[str(prediction)])
            ),
            native_sequence_hint=None if native_sequence_hint is None else str(native_sequence_hint),
            prediction_sequence_hint=(
                None
                if prediction_sequence_hints is None or prediction_sequence_hints.get(str(prediction)) is None
                else str(prediction_sequence_hints[str(prediction)])
            ),
        )
        for prediction in predictions
    ]


def load_benchmark_manifest(manifest_path: str | Path) -> list[BenchmarkJob]:
    path = Path(manifest_path)
    if not path.exists():
        raise ManifestFormatError(f"Manifest file '{path}' does not exist.")

    if path.suffix.lower() == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        entries = payload["entries"] if isinstance(payload, dict) and "entries" in payload else payload
        if not isinstance(entries, list):
            raise ManifestFormatError("JSON manifest must be a list or an object with an 'entries' list.")
        return [_job_from_mapping(entry, path.parent, row_number=index + 1) for index, entry in enumerate(entries)]

    if path.suffix.lower() == ".csv":
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise ManifestFormatError("CSV manifest must include a header row.")
            return [
                _job_from_mapping(row, path.parent, row_number=index + 2)
                for index, row in enumerate(reader)
                if any(value and value.strip() for value in row.values())
            ]

    raise ManifestFormatError("Manifest format must be JSON or CSV.")


def run_benchmark(
    native_file: str | Path | None = None,
    predictions: list[str | Path] | None = None,
    native_index: str | Path | None = None,
    prediction_indices: dict[str, str | Path | None] | None = None,
    native_sequence_hint: str | Path | None = None,
    pvalue_mode: str = "-",
    annotator: MCAnnotateRunner | None = None,
    inclusion_radius: float = 15.0,
    jobs: list[BenchmarkJob] | None = None,
    include_per_residue: bool = False,
    include_secondary_structure: bool = False,
    secondary_structure_runner: MCAnnotateRunner | None = None,
    include_molprobity: bool = False,
    molprobity_runner: MolProbityRunner | None = None,
    repair_missing_atoms: bool = False,
    repair_runner: ArenaRunner | None = None,
    arena_option: int = 5,
) -> BenchmarkResult:
    if jobs is None:
        if native_file is None:
            raise ManifestFormatError("A reference structure is required when benchmark jobs are not provided.")
        jobs = build_benchmark_jobs(
            native_file=native_file,
            predictions=predictions or [],
            native_index=native_index,
            prediction_indices=prediction_indices,
            native_sequence_hint=native_sequence_hint,
        )

    entries: list[BenchmarkEntry] = []
    succeeded = 0
    base_binary = annotator.binary_path if annotator is not None else None
    base_cache_dir = annotator.cache_dir if annotator is not None else None
    referenced_natives = {job.native or (None if native_file is None else str(native_file)) for job in jobs}

    for job in jobs:
        resolved_native = job.native or (None if native_file is None else str(native_file))
        if resolved_native is None:
            raise ManifestFormatError(
                f"Benchmark job for prediction '{job.prediction}' is missing a reference structure."
            )

        try:
            with prepared_structure_pair_context(
                resolved_native,
                job.native_index if job.native_index is not None else native_index,
                job.prediction,
                job.prediction_index,
                native_sequence_hint=job.native_sequence_hint,
                prediction_sequence_hint=job.prediction_sequence_hint,
                repair_missing_atoms=repair_missing_atoms,
                repair_runner=repair_runner,
                arena_option=arena_option,
            ) as prepared:
                annotation_aliases = _prepared_annotation_aliases(
                    prepared,
                    resolved_native,
                    job.prediction,
                    explicit_native_annotation=job.native_annotation,
                    explicit_prediction_annotation=job.prediction_annotation,
                    cache_dir=base_cache_dir,
                )
                job_annotator = MCAnnotateRunner(
                    binary_path=base_binary,
                    cache_dir=base_cache_dir,
                    annotation_overrides=_annotation_overrides(job, resolved_native),
                )
                job_annotator = clone_with_annotation_overrides(job_annotator, annotation_aliases) or job_annotator
                job_secondary_runner = (
                    clone_with_annotation_overrides(secondary_structure_runner, annotation_aliases)
                    if secondary_structure_runner is not None
                    else None
                )
                metrics = calculate_assessment_from_prepared(
                    prepared,
                    pvalue_mode=pvalue_mode,
                    annotator=job_annotator,
                    inclusion_radius=inclusion_radius,
                    include_per_residue=include_per_residue,
                    include_secondary_structure=include_secondary_structure,
                    secondary_structure_runner=job_secondary_runner,
                    include_molprobity=include_molprobity,
                    molprobity_runner=molprobity_runner,
                )
                ready_entry = describe_prepared_pair(
                    prepared,
                    prediction_file=job.prediction,
                    native_file=resolved_native,
                    label=job.label,
                )
            entries.append(
                BenchmarkEntry(
                    native=ready_entry.native,
                    prediction=ready_entry.prediction,
                    label=ready_entry.label,
                    status="ok",
                    native_index=ready_entry.native_index,
                    prediction_index=ready_entry.prediction_index,
                    matched_residues=ready_entry.matched_residues,
                    chain_mappings=ready_entry.chain_mappings,
                    used_sidecar_index=ready_entry.used_sidecar_index,
                    used_sequence_hints=ready_entry.used_sequence_hints,
                    used_normalized_inputs=ready_entry.used_normalized_inputs,
                    used_repaired_inputs=ready_entry.used_repaired_inputs,
                    metrics=metrics,
                    error=None,
                )
            )
            succeeded += 1
        except RNAAssessmentError as exc:
            entries.append(
                BenchmarkEntry(
                    native=resolved_native,
                    prediction=job.prediction,
                    label=job.label,
                    status="error",
                    native_index=None,
                    prediction_index=None,
                    matched_residues=None,
                    chain_mappings=(),
                    used_sidecar_index=False,
                    used_sequence_hints=False,
                    used_normalized_inputs=False,
                    used_repaired_inputs=repair_missing_atoms,
                    metrics=None,
                    error=str(exc),
                )
            )

    return BenchmarkResult(
        reference=str(native_file) if native_file is not None and len(referenced_natives) == 1 else None,
        total_predictions=len(jobs),
        succeeded=succeeded,
        failed=len(jobs) - succeeded,
        entries=tuple(entries),
    )


def _annotation_overrides(job: BenchmarkJob, native_file: str) -> dict[str, str] | None:
    overrides: dict[str, str] = {}
    if job.native_annotation:
        overrides[native_file] = job.native_annotation
    if job.prediction_annotation:
        overrides[job.prediction] = job.prediction_annotation
    return overrides or None


def _prepared_annotation_aliases(
    prepared: PreparedStructurePair,
    native_source: str | Path,
    prediction_source: str | Path,
    *,
    explicit_native_annotation: str | Path | None,
    explicit_prediction_annotation: str | Path | None,
    cache_dir: str | Path | None,
) -> dict[str, str]:
    overrides: dict[str, str] = {}
    native_annotation = existing_annotation_path(
        native_source,
        explicit_annotation=explicit_native_annotation,
        cache_dir=cache_dir,
    )
    if native_annotation is not None:
        overrides[prepared.native.pdb_file] = str(native_annotation)
    prediction_annotation = existing_annotation_path(
        prediction_source,
        explicit_annotation=explicit_prediction_annotation,
        cache_dir=cache_dir,
    )
    if prediction_annotation is not None:
        overrides[prepared.prediction.pdb_file] = str(prediction_annotation)
    return overrides


def _job_from_mapping(row, base_dir: Path, row_number: int) -> BenchmarkJob:
    if not isinstance(row, dict):
        raise ManifestFormatError(f"Manifest row {row_number} must be an object.")

    prediction = _get_required_field(row, row_number, "prediction", aliases=("model",))
    native = _get_optional_field(row, "native", aliases=("reference",))
    native_sequence_hint = _get_optional_field(row, "native_fasta", aliases=("native_sequence",))
    prediction_sequence_hint = _get_optional_field(row, "prediction_fasta", aliases=("prediction_sequence",))
    return BenchmarkJob(
        prediction=_resolve_manifest_path(base_dir, prediction),
        native=_resolve_manifest_path(base_dir, native) if native else None,
        label=_get_optional_field(row, "label", aliases=("name",)),
        native_index=_resolve_manifest_path(base_dir, _get_optional_field(row, "native_index")),
        prediction_index=_resolve_manifest_path(base_dir, _get_optional_field(row, "prediction_index")),
        native_sequence_hint=_resolve_manifest_hint(base_dir, native_sequence_hint),
        prediction_sequence_hint=_resolve_manifest_hint(base_dir, prediction_sequence_hint),
        native_annotation=_resolve_manifest_path(base_dir, _get_optional_field(row, "native_annotation")),
        prediction_annotation=_resolve_manifest_path(base_dir, _get_optional_field(row, "prediction_annotation")),
    )


def _get_required_field(row: dict, row_number: int, field: str, aliases: tuple[str, ...] = ()) -> str:
    value = _get_optional_field(row, field, aliases=aliases)
    if value is None:
        names = ", ".join([field, *aliases])
        raise ManifestFormatError(f"Manifest row {row_number} is missing required field '{names}'.")
    return value


def _get_optional_field(row: dict, field: str, aliases: tuple[str, ...] = ()) -> str | None:
    for key in (field, *aliases):
        if key not in row:
            continue
        value = row[key]
        if value is None:
            return None
        stripped = str(value).strip()
        if not stripped:
            return None
        return stripped
    return None


def _resolve_manifest_path(base_dir: Path, value: str | None) -> str | None:
    if value is None:
        return None
    path = Path(value)
    if not path.is_absolute():
        manifest_relative = (base_dir / path).resolve()
        if manifest_relative.exists():
            path = manifest_relative
        else:
            path = path.resolve()
    return str(path)


def _resolve_manifest_hint(base_dir: Path, value: str | None) -> str | None:
    if value is None:
        return None
    return _resolve_manifest_path(base_dir, value) if _looks_like_path(value) else value


def _looks_like_path(value: str | None) -> bool:
    if value is None:
        return False
    stripped = value.strip()
    if not stripped:
        return False
    return any(token in stripped for token in ("/", "\\", ".", "~")) or stripped.lower().endswith((".fa", ".fasta", ".txt"))
