from __future__ import annotations

import json
from pathlib import Path

from rna_kit import (
    build_assessment_report_document,
    build_secondary_structure_report_document,
    calculate_assessment_from_prepared,
    calculate_lddt,
    prepare_structure_pair,
    write_assessment_html_report,
    write_lddt_html_report,
    write_report_json,
    write_secondary_structure_html_report,
)
from rna_kit.secondary_structure import compare_secondary_structures

from .test_secondary_structure import FakeMCAnnotateRunner, _build_annotation
from .conftest import DATA_DIR


def test_assessment_report_writes_json_and_html(tmp_path: Path) -> None:
    prepared = prepare_structure_pair(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
        DATA_DIR / "14_ChenPostExp_2.pdb",
        DATA_DIR / "14_ChenPostExp_2.index",
    )
    runner = FakeMCAnnotateRunner(
        {
            DATA_DIR / "14_solution_0.pdb": _build_annotation(
                prepared.native,
                chain_id="A",
                pair_positions=[(1, 60), (2, 59)],
            ),
            DATA_DIR / "14_ChenPostExp_2.pdb": _build_annotation(
                prepared.prediction,
                chain_id="U",
                pair_positions=[(1, 60)],
            ),
        }
    )
    assessment = calculate_assessment_from_prepared(
        prepared,
        pvalue_mode="-",
        annotator=None,
        inclusion_radius=15.0,
        include_per_residue=True,
        include_secondary_structure=True,
        secondary_structure_runner=runner,
    )
    document = build_assessment_report_document(
        prepared,
        assessment,
        reference=DATA_DIR / "14_solution_0.pdb",
        prediction=DATA_DIR / "14_ChenPostExp_2.pdb",
        warnings=("example warning",),
        artifacts={"secondary_structure_html": "secondary.html"},
    )

    json_path = write_report_json(document, tmp_path / "assessment.json")
    html_path = write_assessment_html_report(document, tmp_path / "assessment.html")

    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert payload["metadata"]["schema_version"] == "1.0.0"
    assert payload["metrics"]["secondary_structure_f1"] == assessment.secondary_structure_f1
    html_content = html_path.read_text(encoding="utf-8")
    assert "RNA Kit Assessment Report" in html_content
    assert "Secondary Structure" in html_content
    assert "Per-residue lDDT" in html_content
    assert "Residue Heatmap" in html_content
    assert "example warning" in html_content
    assert "FornaC" in html_content
    assert "FornaC 1.1.8" in html_content


def test_secondary_structure_report_writes_json_and_html(tmp_path: Path) -> None:
    prepared = prepare_structure_pair(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
        DATA_DIR / "14_ChenPostExp_2.pdb",
        DATA_DIR / "14_ChenPostExp_2.index",
    )
    runner = FakeMCAnnotateRunner(
        {
            DATA_DIR / "14_solution_0.pdb": _build_annotation(
                prepared.native,
                chain_id="A",
                pair_positions=[(1, 60), (2, 59)],
            ),
            DATA_DIR / "14_ChenPostExp_2.pdb": _build_annotation(
                prepared.prediction,
                chain_id="U",
                pair_positions=[(1, 60)],
            ),
        }
    )
    comparison = compare_secondary_structures(prepared.native, prepared.prediction, runner=runner)
    document = build_secondary_structure_report_document(
        comparison,
        reference=DATA_DIR / "14_solution_0.pdb",
        prediction=DATA_DIR / "14_ChenPostExp_2.pdb",
    )

    json_path = write_report_json(document, tmp_path / "secondary.json")
    html_path = write_secondary_structure_html_report(document, tmp_path / "secondary.html")

    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert payload["comparison"]["f1"] == comparison.f1
    html_content = html_path.read_text(encoding="utf-8")
    assert "RNA Kit Secondary Structure Comparison" in html_content
    assert "False Negatives" in html_content
    assert "FornaC" in html_content
    assert "FornaC 1.1.8" in html_content


def test_lddt_html_report_writes_per_residue_visualization(tmp_path: Path) -> None:
    result = calculate_lddt(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
        DATA_DIR / "14_ChenPostExp_2.pdb",
        DATA_DIR / "14_ChenPostExp_2.index",
        include_per_residue=True,
    )

    html_path = write_lddt_html_report(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_ChenPostExp_2.pdb",
        result,
        tmp_path / "lddt.html",
    )

    html_content = html_path.read_text(encoding="utf-8")
    assert "RNA Kit lDDT Report" in html_content
    assert "Per-residue lDDT" in html_content
    assert "Residue Heatmap" in html_content
    assert "Local RMSD" in html_content
