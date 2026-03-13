from __future__ import annotations

from pathlib import Path

import pytest

from rna_kit import PDBStructure, RNAAssessmentError, infer_structure_alignment, prepare_structure_pair

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


def test_prepare_pair_can_use_sequence_hint_when_residue_names_are_uninformative(tmp_path: Path) -> None:
    native_path = tmp_path / "native.pdb"
    prediction_path = tmp_path / "prediction_bad_names.pdb"
    fasta_path = tmp_path / "prediction.fasta"
    native_text = (DATA_DIR / "14_ChenPostExp_2.pdb").read_text(encoding="utf-8")
    native_path.write_text(native_text, encoding="utf-8")
    prediction_path.write_text(_rewrite_residue_names(native_text, "MOD"), encoding="utf-8")
    fasta_path.write_text(">U\n" + PDBStructure.from_file(native_path).raw_sequence() + "\n", encoding="utf-8")

    with pytest.raises(RNAAssessmentError):
        prepare_structure_pair(
            native_path,
            None,
            prediction_path,
            None,
            resolve_sidecar_indices=False,
        )

    prepared = prepare_structure_pair(
        native_path,
        None,
        prediction_path,
        None,
        resolve_sidecar_indices=False,
        prediction_sequence_hint=fasta_path,
    )

    assert prepared.used_sequence_hints is True
    assert prepared.native_index is not None
    assert prepared.prediction_index is not None
    assert prepared.alignment is not None
    assert len(prepared.native.res_seq) == len(prepared.prediction.res_seq)


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
