from __future__ import annotations

from pathlib import Path

import pytest

from rna_kit import calculate_assessment, calculate_secondary_structure
from rna_kit.mc_annotate import MCAnnotateResult, RawInteraction
from rna_kit.secondary_structure import calculate_secondary_structure_for_structure, compare_secondary_structures
from rna_kit.secondary_structure_web import (
    render_secondary_structure_comparison_html,
    render_secondary_structure_html,
    write_secondary_structure_comparison_html,
)
from rna_kit.structures import PDBStructure

from .conftest import DATA_DIR


class FakeMCAnnotateRunner:
    def __init__(self, annotations: dict[str | Path, MCAnnotateResult]) -> None:
        self.annotations = {
            str(Path(pdb_file).resolve()): annotation for pdb_file, annotation in annotations.items()
        }

    def load(self, pdb_file: str | Path) -> MCAnnotateResult:
        return self.annotations[str(Path(pdb_file).resolve())]


def test_secondary_structure_uses_selected_residues_and_builds_dot_bracket() -> None:
    structure = PDBStructure.from_file(DATA_DIR / "14_solution_0.pdb", DATA_DIR / "14_solution_0.index")
    runner = FakeMCAnnotateRunner(
        {
            DATA_DIR / "14_solution_0.pdb": _build_annotation(
                structure,
                chain_id="A",
                pair_positions=[(1, 60), (2, 59)],
            )
        }
    )

    result = calculate_secondary_structure_for_structure(structure, runner=runner)

    assert result.selected_residues == 60
    assert result.base_pair_count == 2
    assert len(result.dot_bracket) == 60
    assert result.dot_bracket[:2] == "(("
    assert result.dot_bracket.count("(") == 2
    assert result.dot_bracket.count(")") == 2
    assert result.base_pairs[0].rank_1 == 0
    assert result.base_pairs[0].pos_2 == 60
    assert result.base_pairs[0].classification == "cisWW"


def test_secondary_structure_wrapper_does_not_use_sidecar_index_by_default() -> None:
    structure = PDBStructure.from_file(DATA_DIR / "14_solution_0.pdb")
    runner = FakeMCAnnotateRunner(
        {
            DATA_DIR / "14_solution_0.pdb": _build_annotation(
                structure,
                chain_id="A",
                pair_positions=[(1, 61)],
            )
        }
    )

    result = calculate_secondary_structure(
        DATA_DIR / "14_solution_0.pdb",
        None,
        runner=runner,
    )

    assert result.selected_residues == 122
    assert result.base_pair_count == 1
    assert result.base_pairs[0].chain_1 == "A"


def test_secondary_structure_wrapper_uses_explicit_index_when_requested() -> None:
    structure = PDBStructure.from_file(DATA_DIR / "14_solution_0.pdb", DATA_DIR / "14_solution_0.index")
    runner = FakeMCAnnotateRunner(
        {
            DATA_DIR / "14_solution_0.pdb": _build_annotation(
                structure,
                chain_id="A",
                pair_positions=[(1, 60)],
            )
        }
    )

    result = calculate_secondary_structure(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
        runner=runner,
    )

    assert result.selected_residues == 60
    assert result.base_pair_count == 1


def test_secondary_structure_comparison_matches_ranks_across_different_chain_ids() -> None:
    native_structure = PDBStructure.from_file(DATA_DIR / "14_solution_0.pdb", DATA_DIR / "14_solution_0.index")
    prediction_structure = PDBStructure.from_file(
        DATA_DIR / "14_ChenPostExp_2.pdb",
        DATA_DIR / "14_ChenPostExp_2.index",
    )
    runner = FakeMCAnnotateRunner(
        {
            DATA_DIR / "14_solution_0.pdb": _build_annotation(
                native_structure,
                chain_id="A",
                pair_positions=[(1, 60), (2, 59)],
            ),
            DATA_DIR / "14_ChenPostExp_2.pdb": _build_annotation(
                prediction_structure,
                chain_id="U",
                pair_positions=[(1, 60)],
            ),
        }
    )

    result = compare_secondary_structures(native_structure, prediction_structure, runner=runner)

    assert result.true_positives == 1
    assert result.false_positives == 0
    assert result.false_negatives == 1
    assert result.precision == pytest.approx(1.0)
    assert result.recall == pytest.approx(0.5)
    assert result.f1 == pytest.approx(2.0 / 3.0)
    assert result.native.base_pair_count == 2
    assert result.prediction.base_pair_count == 1
    assert len(result.true_positive_pairs) == 1
    assert len(result.false_positive_pairs) == 0
    assert len(result.false_negative_pairs) == 1
    assert result.true_positive_pairs[0].status == "tp"
    assert result.false_negative_pairs[0].status == "fn"


def test_assessment_can_include_secondary_structure_metrics() -> None:
    native_structure = PDBStructure.from_file(DATA_DIR / "14_solution_0.pdb", DATA_DIR / "14_solution_0.index")
    prediction_structure = PDBStructure.from_file(
        DATA_DIR / "14_ChenPostExp_2.pdb",
        DATA_DIR / "14_ChenPostExp_2.index",
    )
    runner = FakeMCAnnotateRunner(
        {
            DATA_DIR / "14_solution_0.pdb": _build_annotation(
                native_structure,
                chain_id="A",
                pair_positions=[(1, 60), (2, 59)],
            ),
            DATA_DIR / "14_ChenPostExp_2.pdb": _build_annotation(
                prediction_structure,
                chain_id="U",
                pair_positions=[(1, 60)],
            ),
        }
    )

    result = calculate_assessment(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
        DATA_DIR / "14_ChenPostExp_2.pdb",
        DATA_DIR / "14_ChenPostExp_2.index",
        include_secondary_structure=True,
        secondary_structure_runner=runner,
    )

    assert result.secondary_structure_precision == pytest.approx(1.0)
    assert result.secondary_structure_recall == pytest.approx(0.5)
    assert result.secondary_structure_f1 == pytest.approx(2.0 / 3.0)
    assert result.secondary_structure_jaccard == pytest.approx(0.5)
    assert result.secondary_structure is not None
    assert result.secondary_structure.native.base_pair_count == 2


def test_secondary_structure_html_renderer_returns_web_markup() -> None:
    structure = PDBStructure.from_file(DATA_DIR / "14_solution_0.pdb", DATA_DIR / "14_solution_0.index")
    runner = FakeMCAnnotateRunner(
        {
            DATA_DIR / "14_solution_0.pdb": _build_annotation(
                structure,
                chain_id="A",
                pair_positions=[(1, 60), (2, 59)],
            )
        }
    )

    result = calculate_secondary_structure_for_structure(structure, runner=runner)
    html_output = render_secondary_structure_html(result, title="Native structure")

    assert html_output.startswith("<!DOCTYPE html>")
    assert "Native structure" in html_output
    assert "fornac" in html_output
    assert "FornaC" in html_output
    assert "FornaC 1.1.8" in html_output
    assert "MC-Annotate cisWW" in html_output


def test_secondary_structure_comparison_html_writer_creates_file(tmp_path: Path) -> None:
    native_structure = PDBStructure.from_file(DATA_DIR / "14_solution_0.pdb", DATA_DIR / "14_solution_0.index")
    prediction_structure = PDBStructure.from_file(
        DATA_DIR / "14_ChenPostExp_2.pdb",
        DATA_DIR / "14_ChenPostExp_2.index",
    )
    runner = FakeMCAnnotateRunner(
        {
            DATA_DIR / "14_solution_0.pdb": _build_annotation(
                native_structure,
                chain_id="A",
                pair_positions=[(1, 60), (2, 59)],
            ),
            DATA_DIR / "14_ChenPostExp_2.pdb": _build_annotation(
                prediction_structure,
                chain_id="U",
                pair_positions=[(1, 60)],
            ),
        }
    )

    result = compare_secondary_structures(native_structure, prediction_structure, runner=runner)
    output_path = write_secondary_structure_comparison_html(result, tmp_path / "comparison.html", title="SS compare")

    content = output_path.read_text(encoding="utf-8")
    assert content.startswith("<!DOCTYPE html>")
    assert "SS compare" in content
    assert "False Negatives" in content
    assert "FornaC" in content
    assert "FornaC 1.1.8" in content


def _build_annotation(
    structure: PDBStructure,
    chain_id: str,
    pair_positions: list[tuple[int, int]],
) -> MCAnnotateResult:
    records = {record.pos: record for _, record in structure.selected_records()}
    interactions = []
    for pos_1, pos_2 in pair_positions:
        record_1 = records[pos_1]
        record_2 = records[pos_2]
        interactions.append(
            RawInteraction(
                type="PAIR_2D",
                chain_a=chain_id,
                pos_a=pos_1,
                nt_a=record_1.nt,
                chain_b=chain_id,
                pos_b=pos_2,
                nt_b=record_2.nt,
                extra_1="WW",
                extra_2="cis",
                extra_3="",
            )
        )
    return MCAnnotateResult(residues=[], interactions=interactions)
