from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO, PDBParser

from .resources import load_residue_mapping


_RESIDUE_MAPPING = load_residue_mapping()


@dataclass(frozen=True)
class ResidueRecord:
    chain: str
    pos: int
    nt: str
    residue: object

    def key(self) -> str:
        return f"{self.chain}:{self.pos}"


class PDBStructure:
    def __init__(self) -> None:
        self._pdb_file: Path | None = None
        self._structure = None
        self._residues: list[ResidueRecord] = []
        self._selected_indices: list[int] = []
        self._residue_index: dict[str, list[int | None]] = {}
        self._interaction_cache: dict[str, list[tuple[str, int, int, str]]] = {}

    @classmethod
    def from_file(cls, pdb_file: str | Path, index_name: str | Path | None = None) -> "PDBStructure":
        structure = cls()
        structure.load(pdb_file, index_name=index_name)
        return structure

    def load(self, pdb_file: str | Path, index_name: str | Path | None = None) -> bool:
        self._pdb_file = Path(pdb_file)
        self._structure = load_structure_file(self._pdb_file)
        self._residues = []
        self._selected_indices = []
        self._residue_index = {}
        self._interaction_cache = {}

        model = self._structure[0]
        for index, chain in enumerate(model.child_list):
            del index
            for residue in chain.child_list:
                record = ResidueRecord(
                    chain=chain.id,
                    pos=residue.id[1],
                    nt=_canonical_residue_name(residue.resname.strip()),
                    residue=residue,
                )
                self._residues.append(record)
                self._residue_index[record.key()] = [len(self._residues) - 1, None]

        if index_name is None:
            self._select_all()
        else:
            self._load_index(Path(index_name))
        return True

    @property
    def pdb_file(self) -> str:
        if self._pdb_file is None:
            raise ValueError("Structure has not been loaded.")
        return str(self._pdb_file)

    @property
    def struct(self):
        if self._structure is None:
            raise ValueError("Structure has not been loaded.")
        return self._structure

    @property
    def res_seq(self) -> list[int]:
        return list(self._selected_indices)

    @property
    def res_list(self) -> list[ResidueRecord]:
        return list(self._residues)

    def raw_sequence(self) -> str:
        return "".join(self._residues[index].nt for index in self._selected_indices)

    def res_sequence(self) -> list[object]:
        return [self._residues[index].residue for index in self._selected_indices]

    def selected_records(self) -> list[tuple[int, ResidueRecord]]:
        return [(index, self._residues[index]) for index in self._selected_indices]

    def chain_records(self) -> dict[str, list[tuple[int, ResidueRecord]]]:
        chains: dict[str, list[tuple[int, ResidueRecord]]] = {}
        for index, record in self.selected_records():
            chains.setdefault(record.chain, []).append((index, record))
        return chains

    def with_selected_indices(self, selected_indices: list[int]) -> "PDBStructure":
        clone = object.__new__(PDBStructure)
        clone._pdb_file = self._pdb_file
        clone._structure = self._structure
        clone._residues = self._residues
        clone._selected_indices = list(selected_indices)
        clone._interaction_cache = {}
        clone._residue_index = {
            key: [value[0], None] for key, value in self._residue_index.items()
        }
        for rank, absolute_index in enumerate(clone._selected_indices):
            clone._residue_index[self._residues[absolute_index].key()][1] = rank
        return clone

    def index_spec(self) -> str:
        ranges: list[str] = []
        start_record: ResidueRecord | None = None
        last_record: ResidueRecord | None = None
        count = 0

        for _, record in self.selected_records():
            if (
                start_record is None
                or last_record is None
                or record.chain != last_record.chain
                or record.pos != last_record.pos + 1
            ):
                if start_record is not None:
                    ranges.append(f"{start_record.chain}:{start_record.pos}:{count}")
                start_record = record
                count = 1
            else:
                count += 1
            last_record = record

        if start_record is not None:
            ranges.append(f"{start_record.chain}:{start_record.pos}:{count}")
        return ",".join(ranges)

    def rank_of(self, chain: str, pos: int) -> int | None:
        entry = self._residue_index.get(f"{chain}:{pos}")
        return None if entry is None else entry[1]

    def cached_interactions(self, cache_key: str) -> list[tuple[str, int, int, str]] | None:
        return self._interaction_cache.get(cache_key)

    def set_cached_interactions(
        self, cache_key: str, interactions: list[tuple[str, int, int, str]]
    ) -> None:
        self._interaction_cache[cache_key] = interactions

    def _select_all(self) -> None:
        for index, residue in enumerate(self._residues):
            self._selected_indices.append(index)
            self._residue_index[residue.key()][1] = len(self._selected_indices) - 1

    def _load_index(self, index_file: Path) -> None:
        entries: list[tuple[str, int, int]] = []
        for row in index_file.read_text(encoding="utf-8").splitlines():
            stripped = row.strip()
            if not stripped or stripped.startswith("#"):
                continue
            for chunk in stripped.split(","):
                chain, pos, count = chunk.split(":")
                entries.append((chain, int(pos), int(count)))

        self._selected_indices = []
        for chain, pos, count in entries:
            start_index = self._absolute_index(chain, pos)
            if start_index is None:
                raise KeyError(f"Bad index key: '{chain}:{pos}'.")
            for absolute_index in range(start_index, start_index + count):
                if absolute_index >= len(self._residues):
                    raise IndexError(f"Index count overflows structure for '{chain}:{pos}:{count}'.")
                record = self._residues[absolute_index]
                if record.chain != chain:
                    raise IndexError(
                        f"Index '{chain}:{pos}:{count}' crosses into chain '{record.chain}'."
                    )
                self._selected_indices.append(absolute_index)
                self._residue_index[record.key()][1] = len(self._selected_indices) - 1

    def _absolute_index(self, chain: str, pos: int) -> int | None:
        entry = self._residue_index.get(f"{chain}:{pos}")
        return None if entry is None else int(entry[0])


def load_structure_file(pdb_file: str | Path):
    path = Path(pdb_file)
    parser = _structure_parser_for(path)
    return parser.get_structure("structure", str(path))


def write_structure_as_pdb(source_file: str | Path, output_file: str | Path):
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(load_structure_file(source_file))
    io.save(str(output_path))
    return output_path


def prepare_external_structure_input(
    source_file: str | Path,
    output_file: str | Path,
):
    source_path = Path(source_file)
    output_path = Path(output_file)
    if _is_mmcif_path(source_path):
        return write_structure_as_pdb(source_path, output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_bytes(source_path.read_bytes())
    return output_path


def _structure_parser_for(path: Path):
    if _is_mmcif_path(path):
        return MMCIFParser(QUIET=True)
    return PDBParser(QUIET=True)


def _is_mmcif_path(path: Path) -> bool:
    return path.suffix.lower() in {".cif", ".mmcif"}


def _canonical_residue_name(name: str) -> str:
    stripped = name.strip()
    mapped = _RESIDUE_MAPPING.get(stripped)
    if mapped and mapped != "-":
        return mapped
    upper = stripped.upper()
    if upper in {"A", "C", "G", "U"}:
        return upper
    if upper == "T":
        return "U"
    for candidate in (upper[:1], upper[-1:]):
        if candidate in {"A", "C", "G", "U"}:
            return candidate
        if candidate == "T":
            return "U"
    return "N"
