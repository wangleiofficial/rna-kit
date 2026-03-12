from __future__ import annotations

from pathlib import Path

from rna_kit import PDBNormalizer, normalize_structure

from .conftest import DATA_DIR


def test_normalize_structure_writes_output(tmp_path: Path) -> None:
    output_path = tmp_path / "normalized.pdb"
    normalize_structure(DATA_DIR / "14_solution_0.pdb", output_path)

    content = output_path.read_text(encoding="utf-8").splitlines()
    assert output_path.exists()
    assert content[0].startswith("ATOM")
    assert content[-1] == "TER"


def test_normalizer_rejects_unknown_residue(tmp_path: Path) -> None:
    pdb_path = tmp_path / "bad.pdb"
    pdb_path.write_text(
        "ATOM      1  O5' XXX A   1      1.000   1.000   1.000  1.00 20.00           O\n",
        encoding="utf-8",
    )

    normalizer = PDBNormalizer.from_defaults()
    assert normalizer.normalize_file(pdb_path, tmp_path / "unused.pdb") is False
    assert any("Unknown residue name" in error for error in normalizer.errors)
