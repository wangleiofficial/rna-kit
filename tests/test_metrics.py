from __future__ import annotations

from pathlib import Path

import pytest

import rna_kit.metrics as metrics_module
from rna_kit import (
    calculate_assessment,
    calculate_ermsd,
    calculate_interaction_network_fidelity,
    calculate_lddt,
    calculate_mcq,
    calculate_rmsd,
)
from rna_kit.metrics import AssessmentResult

from .conftest import DATA_DIR, write_mmcif_from_pdb


def test_rmsd_self_comparison_is_zero() -> None:
    result = calculate_rmsd(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
    )

    assert result.rmsd == pytest.approx(0.0, abs=1e-8)
    assert 0.0 <= result.pvalue <= 1.0


def test_inf_self_comparison_is_one() -> None:
    result = calculate_interaction_network_fidelity(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
    )

    assert result.rmsd == pytest.approx(0.0, abs=1e-8)
    assert result.inf_all == pytest.approx(1.0, abs=1e-8)
    assert result.inf_wc == pytest.approx(1.0, abs=1e-8)
    assert result.inf_nwc == pytest.approx(1.0, abs=1e-8)
    assert result.inf_stack == pytest.approx(1.0, abs=1e-8)


def test_ermsd_self_comparison_is_zero() -> None:
    result = calculate_ermsd(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
    )

    assert result.ermsd == pytest.approx(0.0, abs=1e-8)
    assert result.evaluated_residues == 60


def test_ermsd_matches_expected_cross_structure_value() -> None:
    result = calculate_ermsd(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
        DATA_DIR / "14_ChenPostExp_2.pdb",
        DATA_DIR / "14_ChenPostExp_2.index",
    )

    assert result.ermsd == pytest.approx(1.276258422684251, abs=1e-8)
    assert result.evaluated_residues == 60


def test_mcq_uses_extracted_prepared_subset(monkeypatch) -> None:
    observed: dict[str, Path | str | None] = {}

    def fake_mcq(self, model_file, target_file, jar_path=None):
        observed["model_file"] = Path(model_file)
        observed["target_file"] = Path(target_file)
        observed["jar_path"] = None if jar_path is None else str(jar_path)
        assert Path(model_file).exists()
        assert Path(target_file).exists()
        assert "ATOM" in Path(model_file).read_text(encoding="utf-8")
        assert "ATOM" in Path(target_file).read_text(encoding="utf-8")
        return 0.4321

    monkeypatch.setattr(metrics_module.PDBComparer, "mcq", fake_mcq)

    result = calculate_mcq(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
        DATA_DIR / "14_ChenPostExp_2.pdb",
        DATA_DIR / "14_ChenPostExp_2.index",
    )

    assert result.mcq == pytest.approx(0.4321, abs=1e-8)
    assert result.evaluated_residues == 60
    assert observed["jar_path"] is None
    assert observed["model_file"] != DATA_DIR / "14_ChenPostExp_2.pdb"
    assert observed["target_file"] != DATA_DIR / "14_solution_0.pdb"


def test_cross_structure_metrics_are_finite() -> None:
    result = calculate_interaction_network_fidelity(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
        DATA_DIR / "14_ChenPostExp_2.pdb",
        DATA_DIR / "14_ChenPostExp_2.index",
    )

    assert result.rmsd == pytest.approx(7.751173243045826)
    assert result.deformation_index == pytest.approx(10.643784178530252)
    assert result.inf_all == pytest.approx(0.7282347248904991)
    assert result.inf_wc == pytest.approx(0.9375)
    assert result.inf_nwc == pytest.approx(0.25)
    assert result.inf_stack == pytest.approx(0.7082882469748285)


def test_lddt_self_comparison_is_one() -> None:
    result = calculate_lddt(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
        include_per_residue=True,
    )

    assert result.lddt == pytest.approx(1.0, abs=1e-8)
    assert result.evaluated_atoms > 0
    assert result.evaluated_pairs > 0
    assert result.per_residue is not None
    assert len(result.per_residue) == 60
    assert all(item.lddt == pytest.approx(1.0, abs=1e-8) for item in result.per_residue)
    assert all(item.local_rmsd == pytest.approx(0.0, abs=1e-8) for item in result.per_residue)


def test_assessment_returns_combined_metrics() -> None:
    result = calculate_assessment(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
        DATA_DIR / "14_ChenPostExp_2.pdb",
        DATA_DIR / "14_ChenPostExp_2.index",
        include_per_residue=True,
    )

    assert result.rmsd == pytest.approx(7.751173243045826)
    assert result.pvalue == pytest.approx(7.327471962526033e-15)
    assert result.inf_all == pytest.approx(0.7282347248904991)
    assert result.lddt == pytest.approx(0.6126129382795444)
    assert result.lddt_evaluated_atoms == 1287
    assert result.lddt_evaluated_pairs == 338908
    assert result.ermsd == pytest.approx(1.276258422684251, abs=1e-8)
    assert result.ermsd_evaluated_residues == 60
    assert result.per_residue is not None
    assert len(result.per_residue) == 60
    assert result.per_residue[0].native_chain == "A"
    assert result.per_residue[0].prediction_chain == "U"
    assert result.per_residue[0].matched_atoms > 0
    assert 0.0 <= result.per_residue[0].lddt <= 1.0
    assert result.per_residue[0].local_rmsd is not None


def test_assessment_can_include_mcq(monkeypatch) -> None:
    def fake_calculate_mcq_from_prepared(prepared, jar_path=None):
        assert jar_path == "mcq.jar"
        return metrics_module.MCQResult(mcq=0.4321, evaluated_residues=len(prepared.native.res_seq))

    monkeypatch.setattr(metrics_module, "calculate_mcq_from_prepared", fake_calculate_mcq_from_prepared)

    result = calculate_assessment(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
        DATA_DIR / "14_ChenPostExp_2.pdb",
        DATA_DIR / "14_ChenPostExp_2.index",
        include_mcq=True,
        mcq_jar_path="mcq.jar",
    )

    assert result.mcq == pytest.approx(0.4321, abs=1e-8)
    assert result.mcq_evaluated_residues == 60


def test_lddt_supports_mmcif_input(tmp_path) -> None:
    native_cif = write_mmcif_from_pdb(DATA_DIR / "14_solution_0.pdb", tmp_path / "native.cif")
    prediction_cif = write_mmcif_from_pdb(DATA_DIR / "14_ChenPostExp_2.pdb", tmp_path / "prediction.cif")

    result = calculate_lddt(
        native_cif,
        DATA_DIR / "14_solution_0.index",
        prediction_cif,
        DATA_DIR / "14_ChenPostExp_2.index",
        include_per_residue=True,
    )

    assert 0.0 <= result.lddt <= 1.0
    assert result.per_residue is not None
    assert len(result.per_residue) == 60


def test_assessment_auto_normalizes_when_preparation_fails(monkeypatch, tmp_path: Path) -> None:
    native_input = tmp_path / "native_raw.pdb"
    prediction_input = tmp_path / "prediction_raw.pdb"
    native_input.write_text("RAW\n", encoding="utf-8")
    prediction_input.write_text("RAW\n", encoding="utf-8")

    real_prepare = metrics_module.prepare_structure_pair
    prepare_calls: list[tuple[Path, Path]] = []

    def fake_prepare_structure_pair(
        native_file,
        native_index,
        prediction_file,
        prediction_index,
        **kwargs,
    ):
        native_path = Path(native_file)
        prediction_path = Path(prediction_file)
        prepare_calls.append((native_path, prediction_path))
        if native_path == native_input and prediction_path == prediction_input:
            raise metrics_module.MetricCalculationError("raw structure preparation failed")
        return real_prepare(
            DATA_DIR / "14_solution_0.pdb",
            DATA_DIR / "14_solution_0.index",
            DATA_DIR / "14_solution_0.pdb",
            DATA_DIR / "14_solution_0.index",
            resolve_sidecar_indices=False,
        )

    class FakeNormalizer:
        def normalize_or_raise(self, finput, foutput):
            Path(foutput).write_text((DATA_DIR / "14_solution_0.pdb").read_text(encoding="utf-8"), encoding="utf-8")
            return Path(foutput)

    observed: dict[str, bool] = {}

    def fake_calculate_assessment_from_prepared(prepared, **kwargs):
        observed["used_normalized_inputs"] = prepared.used_normalized_inputs
        return AssessmentResult(
            rmsd=0.0,
            pvalue=0.0,
            deformation_index=0.0,
            inf_all=1.0,
            inf_wc=1.0,
            inf_nwc=1.0,
            inf_stack=1.0,
            lddt=1.0,
            lddt_evaluated_atoms=1,
            lddt_evaluated_pairs=1,
        )

    monkeypatch.setattr(metrics_module, "prepare_structure_pair", fake_prepare_structure_pair)
    monkeypatch.setattr(
        metrics_module,
        "calculate_assessment_from_prepared",
        fake_calculate_assessment_from_prepared,
    )
    monkeypatch.setattr(
        metrics_module,
        "PDBNormalizer",
        type("FakeNormalizerFactory", (), {"from_defaults": staticmethod(lambda: FakeNormalizer())}),
    )

    result = calculate_assessment(
        native_input,
        None,
        prediction_input,
        None,
    )

    assert result.lddt == pytest.approx(1.0, abs=1e-8)
    assert observed["used_normalized_inputs"] is True
    assert len(prepare_calls) == 2
    assert prepare_calls[0] == (native_input, prediction_input)
    assert prepare_calls[1][0] != native_input
    assert prepare_calls[1][1] != prediction_input
