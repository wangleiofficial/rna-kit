from __future__ import annotations

from pathlib import Path

from rna_kit import PDBStructure, infer_structure_alignment, prepare_structure_pair

from .conftest import DATA_DIR


def test_prepare_pair_uses_sidecar_index_when_omitted() -> None:
    prepared = prepare_structure_pair(
        DATA_DIR / "14_solution_0.pdb",
        None,
        DATA_DIR / "14_ChenPostExp_2.pdb",
        None,
    )

    assert prepared.used_sidecar_index is True
    assert prepared.native_index == "A:1:31,A:33:29"
    assert prepared.prediction_index == "U:1:31,U:33:29"


def test_infer_alignment_handles_chain_name_mismatch(tmp_path: Path) -> None:
    reference_path = tmp_path / "reference.pdb"
    prediction_path = tmp_path / "prediction.pdb"
    reference_path.write_text((DATA_DIR / "14_ChenPostExp_2.pdb").read_text(encoding="utf-8"), encoding="utf-8")
    prediction_path.write_text(
        _rewrite_chain_ids((DATA_DIR / "14_ChenPostExp_2.pdb").read_text(encoding="utf-8"), "Z"),
        encoding="utf-8",
    )

    reference = PDBStructure.from_file(reference_path)
    prediction = PDBStructure.from_file(prediction_path)
    alignment = infer_structure_alignment(reference, prediction)

    assert alignment.matched_residues == len(reference.raw_sequence())
    assert alignment.chain_alignments[0].native_chain == "U"
    assert alignment.chain_alignments[0].prediction_chain == "Z"


def _rewrite_chain_ids(pdb_text: str, chain_id: str) -> str:
    rows = []
    for row in pdb_text.splitlines():
        if row.startswith(("ATOM  ", "HETATM")):
            rows.append(f"{row[:21]}{chain_id}{row[22:]}")
        else:
            rows.append(row)
    return "\n".join(rows) + "\n"
