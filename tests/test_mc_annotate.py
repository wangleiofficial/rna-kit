from __future__ import annotations

import pytest

from rna_kit import MCAnnotateRunner, ToolNotAvailableError

from .conftest import DATA_DIR


def test_precomputed_annotation_is_used_without_binary() -> None:
    runner = MCAnnotateRunner(binary_path="/definitely/missing")
    result = runner.load(DATA_DIR / "14_solution_0.pdb")

    assert len(result.residues) > 0
    assert len(result.interactions) > 0


def test_missing_annotation_raises_without_available_binary(tmp_path) -> None:
    pdb_path = tmp_path / "model.pdb"
    pdb_path.write_text((DATA_DIR / "14_solution_0.pdb").read_text(encoding="utf-8"), encoding="utf-8")

    runner = MCAnnotateRunner(binary_path=tmp_path / "missing-binary")
    with pytest.raises(ToolNotAvailableError):
        runner.load(pdb_path)


def test_annotation_override_path_is_supported(tmp_path) -> None:
    pdb_path = tmp_path / "model.pdb"
    annotation_path = tmp_path / "custom_output.mcout"
    pdb_path.write_text((DATA_DIR / "14_solution_0.pdb").read_text(encoding="utf-8"), encoding="utf-8")
    annotation_path.write_text((DATA_DIR / "14_solution_0.pdb.mcout").read_text(encoding="utf-8"), encoding="utf-8")

    runner = MCAnnotateRunner(
        binary_path=tmp_path / "missing-binary",
        annotation_overrides={pdb_path: annotation_path},
    )
    result = runner.load(pdb_path)

    assert len(result.interactions) > 0
