from __future__ import annotations

from importlib.resources import files


def load_residue_mapping() -> dict[str, str]:
    return _load_mapping("residues.list")


def load_atom_mapping() -> dict[str, str]:
    return _load_mapping("atoms.list")


def _load_mapping(filename: str) -> dict[str, str]:
    rows = files("rna_kit.data").joinpath(filename).read_text(encoding="utf-8").splitlines()
    mapping: dict[str, str] = {}
    for row in rows:
        stripped = row.strip()
        if not stripped or stripped.startswith("#"):
            continue
        source, target = stripped.split()
        mapping[source] = target
    return mapping
