from __future__ import annotations

import pytest

from rna_kit import calculate_assessment, calculate_interaction_network_fidelity, calculate_lddt, calculate_rmsd

from .conftest import DATA_DIR


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
    assert result.per_residue is not None
    assert len(result.per_residue) == 60
    assert result.per_residue[0].native_chain == "A"
    assert result.per_residue[0].prediction_chain == "U"
    assert result.per_residue[0].matched_atoms > 0
    assert 0.0 <= result.per_residue[0].lddt <= 1.0
    assert result.per_residue[0].local_rmsd is not None
