from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

from rna_kit.structures import PDBStructure

from .conftest import DATA_DIR, PROJECT_ROOT


def test_rmsd_cli_returns_json() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "rmsd",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_solution_0.index"),
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_solution_0.index"),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["rmsd"] == pytest.approx(0.0, abs=1e-8)


def test_ermsd_cli_returns_json() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "ermsd",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_ChenPostExp_2.pdb"),
            "--native-index",
            str(DATA_DIR / "14_solution_0.index"),
            "--prediction-index",
            str(DATA_DIR / "14_ChenPostExp_2.index"),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["ermsd"] == pytest.approx(1.276258422684251, abs=1e-8)
    assert payload["evaluated_residues"] == 60


def test_normalize_cli_writes_file(tmp_path: Path) -> None:
    output_path = tmp_path / "normalized.pdb"
    subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "normalize",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(output_path),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    assert output_path.exists()


def test_assess_cli_returns_combined_metrics() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "assess",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_ChenPostExp_2.pdb"),
            "--native-index",
            str(DATA_DIR / "14_solution_0.index"),
            "--prediction-index",
            str(DATA_DIR / "14_ChenPostExp_2.index"),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["rmsd"] == pytest.approx(7.751173243045826)
    assert payload["lddt"] == pytest.approx(0.6126129382795444)
    assert payload["ermsd"] == pytest.approx(1.276258422684251, abs=1e-8)


def test_assess_cli_uses_sidecar_indices_when_flags_are_omitted() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "assess",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_ChenPostExp_2.pdb"),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["rmsd"] == pytest.approx(7.751173243045826)
    assert payload["inf_all"] == pytest.approx(0.7282347248904991)


def test_lddt_cli_returns_json() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "lddt",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_solution_0.pdb"),
            "--native-index",
            str(DATA_DIR / "14_solution_0.index"),
            "--prediction-index",
            str(DATA_DIR / "14_solution_0.index"),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["lddt"] == pytest.approx(1.0, abs=1e-8)


def test_lddt_cli_can_write_html_visualization(tmp_path: Path) -> None:
    html_path = tmp_path / "lddt.html"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "lddt",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_ChenPostExp_2.pdb"),
            "--native-index",
            str(DATA_DIR / "14_solution_0.index"),
            "--prediction-index",
            str(DATA_DIR / "14_ChenPostExp_2.index"),
            "--html",
            str(html_path),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["html_output"] == str(html_path)
    assert payload["per_residue"] is None
    html_content = html_path.read_text(encoding="utf-8")
    assert "RNA Kit lDDT Report" in html_content
    assert "Per-residue lDDT" in html_content
    assert "Residue Heatmap" in html_content


def test_assess_cli_can_emit_per_residue_report() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "assess",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_ChenPostExp_2.pdb"),
            "--native-index",
            str(DATA_DIR / "14_solution_0.index"),
            "--prediction-index",
            str(DATA_DIR / "14_ChenPostExp_2.index"),
            "--per-residue",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert len(payload["per_residue"]) == 60
    assert payload["per_residue"][0]["native_chain"] == "A"
    assert payload["per_residue"][0]["prediction_chain"] == "U"


def test_assess_cli_can_repair_missing_atoms_with_arena(tmp_path: Path) -> None:
    arena_path = _write_fake_arena_script(tmp_path)

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "assess",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_ChenPostExp_2.pdb"),
            "--native-index",
            str(DATA_DIR / "14_solution_0.index"),
            "--prediction-index",
            str(DATA_DIR / "14_ChenPostExp_2.index"),
            "--repair-missing-atoms",
            "--arena",
            str(arena_path),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["used_repaired_inputs"] is True
    assert payload["rmsd"] == pytest.approx(7.751173243045826)


def test_map_cli_reports_inferred_chain_mapping(tmp_path: Path) -> None:
    reference_path = tmp_path / "reference.pdb"
    prediction_path = tmp_path / "prediction.pdb"
    reference_path.write_text((DATA_DIR / "14_ChenPostExp_2.pdb").read_text(encoding="utf-8"), encoding="utf-8")
    prediction_path.write_text(
        _rewrite_chain_ids((DATA_DIR / "14_ChenPostExp_2.pdb").read_text(encoding="utf-8"), "Z"),
        encoding="utf-8",
    )

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "map",
            str(reference_path),
            str(prediction_path),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["matched_residues"] > 0
    assert payload["chain_mappings"][0]["native_chain"] == "U"
    assert payload["chain_mappings"][0]["prediction_chain"] == "Z"


def test_benchmark_cli_returns_multiple_results() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "benchmark",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_ChenPostExp_2.pdb"),
            "--native-index",
            str(DATA_DIR / "14_solution_0.index"),
            "--sort-by",
            "rmsd",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["total_predictions"] == 2
    assert payload["succeeded"] == 2
    assert payload["entries"][0]["metrics"]["rmsd"] == pytest.approx(0.0, abs=1e-8)


def test_benchmark_cli_can_repair_missing_atoms_with_arena(tmp_path: Path) -> None:
    arena_path = _write_fake_arena_script(tmp_path)

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "benchmark",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_solution_0.pdb"),
            "--native-index",
            str(DATA_DIR / "14_solution_0.index"),
            "--repair-missing-atoms",
            "--arena",
            str(arena_path),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["entries"][0]["used_repaired_inputs"] is True
    assert payload["entries"][0]["metrics"]["rmsd"] == pytest.approx(0.0, abs=1e-8)


def test_benchmark_cli_supports_json_manifest_and_per_residue(tmp_path: Path) -> None:
    manifest_path = tmp_path / "benchmark.json"
    manifest_path.write_text(
        json.dumps(
            [
                {
                    "label": "self",
                    "native": str(DATA_DIR / "14_solution_0.pdb"),
                    "native_index": str(DATA_DIR / "14_solution_0.index"),
                    "prediction": str(DATA_DIR / "14_solution_0.pdb"),
                },
                {
                    "label": "chen",
                    "native": str(DATA_DIR / "14_solution_0.pdb"),
                    "native_index": str(DATA_DIR / "14_solution_0.index"),
                    "prediction": str(DATA_DIR / "14_ChenPostExp_2.pdb"),
                },
            ]
        ),
        encoding="utf-8",
    )

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "benchmark",
            "--manifest",
            str(manifest_path),
            "--per-residue",
            "--sort-by",
            "rmsd",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["total_predictions"] == 2
    assert payload["entries"][0]["label"] == "self"
    assert len(payload["entries"][1]["metrics"]["per_residue"]) == 60


def test_benchmark_cli_supports_csv_manifest(tmp_path: Path) -> None:
    manifest_path = tmp_path / "benchmark.csv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["label", "native", "native_index", "prediction", "prediction_index"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "label": "self",
                "native": str(DATA_DIR / "14_solution_0.pdb"),
                "native_index": str(DATA_DIR / "14_solution_0.index"),
                "prediction": str(DATA_DIR / "14_solution_0.pdb"),
                "prediction_index": str(DATA_DIR / "14_solution_0.index"),
            }
        )

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "benchmark",
            "--manifest",
            str(manifest_path),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["total_predictions"] == 1
    assert payload["entries"][0]["label"] == "self"
    assert payload["entries"][0]["metrics"]["rmsd"] == pytest.approx(0.0, abs=1e-8)


def test_benchmark_cli_can_write_html_dashboard(tmp_path: Path) -> None:
    html_report = tmp_path / "benchmark.html"
    molprobity_path = _write_fake_molprobity_script(tmp_path)

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "benchmark",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_ChenPostExp_2.pdb"),
            "--native-index",
            str(DATA_DIR / "14_solution_0.index"),
            "--include-molprobity",
            "--molprobity",
            str(molprobity_path),
            "--html-report",
            str(html_report),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["html_report_output"] == str(html_report)
    detail_dir = Path(payload["detail_reports_dir"])
    assert detail_dir.exists()
    assert any(detail_dir.glob("*.html"))
    content = html_report.read_text(encoding="utf-8")
    assert "RNA Kit Benchmark Dashboard" in content
    assert "Best clashscore" in content
    assert "5.4200" in content
    assert "_reports/" in content


def test_secondary_structure_cli_returns_json(tmp_path: Path) -> None:
    annotation_path = _write_mcout_file(
        tmp_path / "14_solution_0.pdb.mcout",
        PDBStructure.from_file(DATA_DIR / "14_solution_0.pdb"),
        chain_id="A",
        pair_positions=[(1, 61), (2, 60)],
    )

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "secondary-structure",
            str(DATA_DIR / "14_solution_0.pdb"),
            "--annotation",
            str(annotation_path),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["base_pair_count"] == 2
    assert payload["selected_residues"] == 122
    assert payload["dot_bracket"][:2] == "(("


def test_assess_cli_can_include_secondary_structure_metrics(tmp_path: Path) -> None:
    native_annotation = _write_mcout_file(
        tmp_path / "14_solution_0.pdb.mcout",
        PDBStructure.from_file(DATA_DIR / "14_solution_0.pdb", DATA_DIR / "14_solution_0.index"),
        chain_id="A",
        pair_positions=[(1, 60), (2, 59)],
    )
    prediction_annotation = _write_mcout_file(
        tmp_path / "14_ChenPostExp_2.pdb.mcout",
        PDBStructure.from_file(DATA_DIR / "14_ChenPostExp_2.pdb", DATA_DIR / "14_ChenPostExp_2.index"),
        chain_id="U",
        pair_positions=[(1, 60)],
    )

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "assess",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_ChenPostExp_2.pdb"),
            "--native-index",
            str(DATA_DIR / "14_solution_0.index"),
            "--prediction-index",
            str(DATA_DIR / "14_ChenPostExp_2.index"),
            "--secondary-structure",
            "--native-annotation",
            str(native_annotation),
            "--prediction-annotation",
            str(prediction_annotation),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["secondary_structure_precision"] == pytest.approx(1.0)
    assert payload["secondary_structure_recall"] == pytest.approx(0.5)
    assert payload["secondary_structure_f1"] == pytest.approx(2.0 / 3.0)
    assert payload["secondary_structure"]["native"]["base_pair_count"] == 2


def test_secondary_compare_cli_can_write_html_visualization(tmp_path: Path) -> None:
    native_annotation = _write_mcout_file(
        tmp_path / "14_solution_0.pdb.mcout",
        PDBStructure.from_file(DATA_DIR / "14_solution_0.pdb", DATA_DIR / "14_solution_0.index"),
        chain_id="A",
        pair_positions=[(1, 60), (2, 59)],
    )
    prediction_annotation = _write_mcout_file(
        tmp_path / "14_ChenPostExp_2.pdb.mcout",
        PDBStructure.from_file(DATA_DIR / "14_ChenPostExp_2.pdb", DATA_DIR / "14_ChenPostExp_2.index"),
        chain_id="U",
        pair_positions=[(1, 60)],
    )
    html_path = tmp_path / "secondary_compare.html"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "secondary-compare",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_ChenPostExp_2.pdb"),
            "--native-index",
            str(DATA_DIR / "14_solution_0.index"),
            "--prediction-index",
            str(DATA_DIR / "14_ChenPostExp_2.index"),
            "--native-annotation",
            str(native_annotation),
            "--prediction-annotation",
            str(prediction_annotation),
            "--html",
            str(html_path),
            "--title",
            "Secondary compare",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["html_output"] == str(html_path)
    content = html_path.read_text(encoding="utf-8")
    assert content.startswith("<!DOCTYPE html>")
    assert "Secondary compare" in content
    assert "FornaC" in content
    assert "FornaC 1.1.8" in content


def test_tools_cli_reports_bundled_tools() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "tools",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    tools = {item["key"]: item for item in payload["tools"]}
    assert "arena" in tools
    assert tools["arena"]["supports_auto_download"] is True
    assert tools["cssr"]["available"] is True
    assert tools["mc_annotate"]["available"] is True
    assert tools["us_align"]["available"] is True


def test_repair_cli_writes_output(tmp_path: Path) -> None:
    arena_path = _write_fake_arena_script(tmp_path)
    output_path = tmp_path / "repaired.pdb"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "repair",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(output_path),
            "--arena",
            str(arena_path),
            "--arena-option",
            "5",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert output_path.exists()
    assert payload["output_file"] == str(output_path)
    assert payload["option"] == 5
    assert payload["used_auto_build"] is False


def test_map_cli_can_use_prediction_fasta_hint(tmp_path: Path) -> None:
    native_path = tmp_path / "native.pdb"
    prediction_path = tmp_path / "prediction_bad_names.pdb"
    fasta_path = tmp_path / "prediction.fasta"
    native_text = (DATA_DIR / "14_ChenPostExp_2.pdb").read_text(encoding="utf-8")
    native_path.write_text(native_text, encoding="utf-8")
    prediction_path.write_text(_rewrite_residue_names(native_text, "MOD"), encoding="utf-8")
    fasta_path.write_text(">U\n" + PDBStructure.from_file(native_path).raw_sequence() + "\n", encoding="utf-8")

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "map",
            str(native_path),
            str(prediction_path),
            "--prediction-fasta",
            str(fasta_path),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["used_sequence_hints"] is True
    assert payload["matched_residues"] > 0


def test_molprobity_cli_returns_json(tmp_path: Path) -> None:
    molprobity_path = _write_fake_molprobity_script(tmp_path)

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "molprobity",
            str(DATA_DIR / "14_solution_0.pdb"),
            "--molprobity",
            str(molprobity_path),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["clashscore"] == pytest.approx(5.42)
    assert payload["molprobity_score"] == pytest.approx(2.11)


def test_assess_cli_can_include_molprobity_metrics(tmp_path: Path) -> None:
    molprobity_path = _write_fake_molprobity_script(tmp_path)

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "assess",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_ChenPostExp_2.pdb"),
            "--native-index",
            str(DATA_DIR / "14_solution_0.index"),
            "--prediction-index",
            str(DATA_DIR / "14_ChenPostExp_2.index"),
            "--include-molprobity",
            "--molprobity",
            str(molprobity_path),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["molprobity"]["clashscore"] == pytest.approx(5.42)
    assert payload["molprobity"]["molprobity_score"] == pytest.approx(2.11)


def test_us_align_cli_can_write_html_viewer(tmp_path: Path) -> None:
    us_align_path = _write_fake_usalign_script(tmp_path)
    html_path = tmp_path / "us_align.html"
    output_dir = tmp_path / "us_align_assets"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "us-align",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_ChenPostExp_2.pdb"),
            "--us-align",
            str(us_align_path),
            "--output-dir",
            str(output_dir),
            "--html",
            str(html_path),
            "--title",
            "US-align Viewer",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["tm_score_reference"] == pytest.approx(0.41234)
    assert payload["tm_score_prediction"] == pytest.approx(0.29691)
    assert payload["html_output"] == str(html_path)
    assert Path(payload["superposed_prediction_output"]).exists()
    assert Path(payload["reference_structure_output"]).exists()
    assert "3Dmol-min.js" in html_path.read_text(encoding="utf-8")
    assert "US-align Viewer" in html_path.read_text(encoding="utf-8")
    assert "Prediction lDDT" in html_path.read_text(encoding="utf-8")
    assert "const residueScores = [" in html_path.read_text(encoding="utf-8")


def test_secondary_compare_cli_can_write_json_and_html_reports(tmp_path: Path) -> None:
    native_annotation = _write_mcout_file(
        tmp_path / "14_solution_0.pdb.mcout",
        PDBStructure.from_file(DATA_DIR / "14_solution_0.pdb", DATA_DIR / "14_solution_0.index"),
        chain_id="A",
        pair_positions=[(1, 60), (2, 59)],
    )
    prediction_annotation = _write_mcout_file(
        tmp_path / "14_ChenPostExp_2.pdb.mcout",
        PDBStructure.from_file(DATA_DIR / "14_ChenPostExp_2.pdb", DATA_DIR / "14_ChenPostExp_2.index"),
        chain_id="U",
        pair_positions=[(1, 60)],
    )
    json_report = tmp_path / "secondary_compare.json"
    html_report = tmp_path / "secondary_compare.html"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "secondary-compare",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_ChenPostExp_2.pdb"),
            "--native-index",
            str(DATA_DIR / "14_solution_0.index"),
            "--prediction-index",
            str(DATA_DIR / "14_ChenPostExp_2.index"),
            "--native-annotation",
            str(native_annotation),
            "--prediction-annotation",
            str(prediction_annotation),
            "--json-report",
            str(json_report),
            "--html-report",
            str(html_report),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["json_report_output"] == str(json_report)
    assert payload["html_report_output"] == str(html_report)
    report_payload = json.loads(json_report.read_text(encoding="utf-8"))
    assert report_payload["metadata"]["schema_version"] == "1.0.0"
    assert "False Negatives" in html_report.read_text(encoding="utf-8")
    assert "FornaC" in html_report.read_text(encoding="utf-8")
    assert "FornaC 1.1.8" in html_report.read_text(encoding="utf-8")


def test_assess_cli_html_report_includes_per_residue_visualization(tmp_path: Path) -> None:
    html_report = tmp_path / "assessment.html"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_kit",
            "assess",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_ChenPostExp_2.pdb"),
            "--native-index",
            str(DATA_DIR / "14_solution_0.index"),
            "--prediction-index",
            str(DATA_DIR / "14_ChenPostExp_2.index"),
            "--html-report",
            str(html_report),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["html_report_output"] == str(html_report)
    assert payload["per_residue"] is None
    html_content = html_report.read_text(encoding="utf-8")
    assert "Per-residue lDDT" in html_content
    assert "Residue Heatmap" in html_content


def _rewrite_chain_ids(pdb_text: str, chain_id: str) -> str:
    rows = []
    for row in pdb_text.splitlines():
        if row.startswith(("ATOM  ", "HETATM")):
            rows.append(f"{row[:21]}{chain_id}{row[22:]}")
        else:
            rows.append(row)
    return "\n".join(rows) + "\n"


def _rewrite_residue_names(pdb_text: str, residue_name: str) -> str:
    rows = []
    formatted = f"{residue_name:>3}"[:3]
    for row in pdb_text.splitlines():
        if row.startswith(("ATOM  ", "HETATM")):
            rows.append(f"{row[:17]}{formatted}{row[20:]}")
        else:
            rows.append(row)
    return "\n".join(rows) + "\n"


def _write_fake_usalign_script(tmp_path: Path) -> Path:
    script_path = tmp_path / "fake_usalign.py"
    script_path.write_text(
        "\n".join(
            [
                "#!/usr/bin/env python3",
                "import sys",
                "from pathlib import Path",
                "args = sys.argv[1:]",
                "prediction = Path(args[0])",
                "reference = Path(args[1])",
                "prefix = None",
                "for index, arg in enumerate(args):",
                "    if arg == '-o':",
                "        prefix = Path(args[index + 1])",
                "        break",
                "if prefix is None:",
                "    raise SystemExit(2)",
                "prefix.parent.mkdir(parents=True, exist_ok=True)",
                f"PDB_TEXT = {json.dumps(_build_fake_usalign_pdb())}",
                "prefix.with_suffix('.pdb').write_text(PDB_TEXT, encoding='utf-8')",
                "print(' ********************************************************************')",
                "print(' * US-align (Version 20241108)                                      *')",
                "print(' ********************************************************************')",
                "print()",
                "print(f'Name of Structure_1: {prediction}:A (to be superimposed onto Structure_2)')",
                "print(f'Name of Structure_2: {reference}:U')",
                "print('Length of Structure_1: 61 residues')",
                "print('Length of Structure_2: 61 residues')",
                "print()",
                "print('Aligned length= 36, RMSD=   3.72, Seq_ID=n_identical/n_aligned= 0.889')",
                "print('TM-score= 0.29691 (normalized by length of Structure_1: L=61, d0=2.17)')",
                "print('TM-score= 0.41234 (normalized by length of Structure_2: L=61, d0=2.17)')",
                "print('(You should use TM-score normalized by length of the reference structure)')",
                "print()",
                "print('(\":\" denotes residue pairs of d < 5.0 Angstrom, \".\" denotes other aligned residues)')",
                "print('AACGUU')",
                "print('::::..')",
                "print('AACGAA')",
                "print()",
                "print('#Total CPU time is  0.01 seconds')",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    script_path.chmod(0o755)
    return script_path


def _write_fake_molprobity_script(tmp_path: Path) -> Path:
    script_path = tmp_path / "fake_molprobity.py"
    script_path.write_text(
        "\n".join(
            [
                "#!/usr/bin/env python3",
                "print('All-atom clashscore = 5.42')",
                "print('MolProbity score = 2.11')",
                "print('Bad bonds = 1')",
                "print('Bad angles = 3')",
                "print('Pucker outliers = 2')",
                "print('Suite outliers = 4')",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    script_path.chmod(0o755)
    return script_path


def _write_fake_arena_script(tmp_path: Path) -> Path:
    script_path = tmp_path / "fake_arena.py"
    script_path.write_text(
        "\n".join(
            [
                "#!/usr/bin/env python3",
                "import sys",
                "from pathlib import Path",
                "input_path = Path(sys.argv[1])",
                "output_path = Path(sys.argv[2])",
                "output_path.write_text(input_path.read_text(encoding='utf-8'), encoding='utf-8')",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    script_path.chmod(0o755)
    return script_path


def _write_mcout_file(
    output_path: Path,
    structure: PDBStructure,
    chain_id: str,
    pair_positions: list[tuple[int, int]],
) -> Path:
    output_path.write_text(
        _build_mcout_output(structure, chain_id=chain_id, pair_positions=pair_positions),
        encoding="utf-8",
    )
    return output_path


def _build_mcout_output(
    structure: PDBStructure,
    chain_id: str,
    pair_positions: list[tuple[int, int]],
) -> str:
    records = {record.pos: record for _, record in structure.selected_records()}
    rows = [
        "",
        "Base-pairs ------------------------------------------------------",
    ]
    for pos_1, pos_2 in pair_positions:
        record_1 = records[pos_1]
        record_2 = records[pos_2]
        rows.append(
            f"{chain_id}{pos_1}-{chain_id}{pos_2} : "
            f"{record_1.nt}-{record_2.nt} Ww/Ww pairing antiparallel cis XIX "
        )
    return "\n".join(rows) + "\n"


def _build_fake_usalign_pdb() -> str:
    return (
        "ATOM      1  P     G A   1      10.000  10.000  10.000  1.00 20.00           P  \n"
        "ATOM      2  O5'   G A   1      11.000  10.000  10.000  1.00 20.00           O  \n"
        "TER\n"
        "ATOM      3  P     G B   1      20.000  20.000  20.000  1.00 20.00           P  \n"
        "ATOM      4  O5'   G B   1      21.000  20.000  20.000  1.00 20.00           O  \n"
        "TER\n"
        "END\n"
    )
