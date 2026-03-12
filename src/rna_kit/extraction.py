from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from Bio.PDB import PDBIO, PDBParser, Select

from .exceptions import InvalidResidueRangeError


@dataclass(frozen=True)
class ResidueRange:
    chain: str
    start: int
    count: int


def parse_residue_ranges(specification: str) -> list[ResidueRange]:
    ranges: list[ResidueRange] = []
    for chunk in specification.split(","):
        part = chunk.strip()
        if not part:
            continue
        pieces = part.split(":")
        if len(pieces) != 3:
            raise InvalidResidueRangeError(
                f"Invalid residue range '{part}'. Expected 'chain:start:count'."
            )
        chain, start, count = pieces
        ranges.append(ResidueRange(chain=chain, start=int(start), count=int(count)))
    if not ranges:
        raise InvalidResidueRangeError("At least one residue range is required.")
    return ranges


def extract_pdb(
    input_path: str | Path,
    residue_ranges: str | Iterable[ResidueRange],
    output_path: str | Path,
) -> Path:
    ranges = (
        parse_residue_ranges(residue_ranges)
        if isinstance(residue_ranges, str)
        else list(residue_ranges)
    )
    selected_keys = {
        f"{residue_range.chain}|{residue_id}"
        for residue_range in ranges
        for residue_id in range(residue_range.start, residue_range.start + residue_range.count)
    }

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("input", str(input_path))

    class ResidueSelection(Select):
        def accept_residue(self, residue) -> int:
            key = f"{residue.parent.id.strip()}|{residue.get_id()[1]}"
            return int(key in selected_keys)

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output), ResidueSelection())
    return output


def extract_PDB(input_path: str | Path, residue_ranges: str, output_path: str | Path) -> Path:
    return extract_pdb(input_path, residue_ranges, output_path)
