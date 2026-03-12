from __future__ import annotations

from pathlib import Path

from .exceptions import NormalizationError
from .resources import load_atom_mapping, load_residue_mapping


class PDBNormalizer:
    def __init__(
        self,
        residue_mapping: dict[str, str] | None = None,
        atom_mapping: dict[str, str] | None = None,
    ) -> None:
        self._residue_mapping = residue_mapping or load_residue_mapping()
        self._atom_mapping = atom_mapping or load_atom_mapping()
        self.errors: list[str] = []
        self._row_count = 0
        self._chain_found = False
        self._in_atom = False
        self._in_model = False

    @classmethod
    def from_defaults(cls) -> "PDBNormalizer":
        return cls()

    def parse(self, finput: str | Path, foutput: str | Path) -> bool:
        return self.normalize_file(finput, foutput)

    def normalize_file(self, finput: str | Path, foutput: str | Path) -> bool:
        self.errors = []
        self._row_count = 0
        self._chain_found = False
        self._in_atom = False
        self._in_model = False

        input_path = Path(finput)
        output_path = Path(foutput)
        output_rows: list[str] = []

        for raw_row in input_path.read_text(encoding="utf-8").splitlines():
            self._row_count += 1
            row = raw_row.rstrip("\n\r")
            rec_name = row[:6]

            if rec_name == "MODEL ":
                normalized = self._parse_model()
            elif rec_name == "ENDMDL":
                normalized = self._parse_endmdl()
            elif rec_name[:3] == "TER":
                normalized = self._parse_ter()
            elif rec_name in {"ATOM  ", "HETATM"}:
                normalized = self._parse_atom(row)
            else:
                continue

            if normalized:
                output_rows.append(normalized)

        if self._in_atom:
            output_rows.append("TER")

        if self.errors:
            return False

        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("\n".join(output_rows) + "\n", encoding="utf-8")
        return True

    def normalize_or_raise(self, finput: str | Path, foutput: str | Path) -> Path:
        output_path = Path(foutput)
        if not self.normalize_file(finput, output_path):
            raise NormalizationError("; ".join(self.errors))
        return output_path

    def _parse_model(self) -> str:
        if self._in_model:
            self._record_error("'ENDMDL' not found before the next MODEL record.")
        if self._in_atom:
            self._record_error("Missing 'MODEL' before atom declarations.")
        self._in_model = True
        return ""

    def _parse_endmdl(self) -> str:
        self._in_model = False
        self._in_atom = False
        return "ENDMDL"

    def _parse_ter(self) -> str:
        result = "TER" if self._in_atom else ""
        self._in_atom = False
        return result

    def _parse_atom(self, row: str) -> str:
        serial = row[6:11].strip()
        name = row[12:16].strip()
        alt_loc = row[16:17] or " "
        residue_name = row[17:20].strip()
        chain_id = row[21:22] or " "
        residue_seq = row[22:26].strip()
        insertion_code = row[26:27] or " "
        x = row[30:38].strip()
        y = row[38:46].strip()
        z = row[46:54].strip()
        temp_factor = row[60:66].strip() or "0.00"
        element = row[76:78].strip()
        charge = row[78:80].strip()

        residue_name = self._normalize_residue_name(residue_name)
        if residue_name is None:
            return ""

        atom_name = self._normalize_atom_name(name, residue_name)
        if atom_name is None:
            return ""

        chain_id = self._normalize_chain_id(chain_id)
        if chain_id is None:
            return ""

        self._in_atom = True
        formatted_name = f"{atom_name:<4}"
        formatted_element = element or atom_name[0]
        occupancy = "1.00"

        return (
            f"ATOM  {serial:>5} {formatted_name}{alt_loc}{residue_name:>3} {chain_id}"
            f"{int(residue_seq):>4}{insertion_code}   {x:>8}{y:>8}{z:>8}{occupancy:>6}"
            f"{temp_factor:>6}          {formatted_element:>2}{charge:>2}"
        )

    def _normalize_residue_name(self, residue_name: str) -> str | None:
        normalized = self._residue_mapping.get(residue_name)
        if normalized is None:
            self._record_error(f"Unknown residue name: '{residue_name}'.")
            return None
        if normalized == "-":
            return None
        return normalized

    def _normalize_atom_name(self, atom_name: str, residue_name: str) -> str | None:
        normalized = self._atom_mapping.get(atom_name)
        if normalized is None:
            self._record_error(f"Unknown atom name: '{atom_name}' in residue '{residue_name}'.")
            return None
        if normalized == "-":
            return None
        return normalized

    def _normalize_chain_id(self, chain_id: str) -> str | None:
        if chain_id != " ":
            self._chain_found = True
            return chain_id
        if self._chain_found:
            self._record_error("One of the chains is missing a chain identifier.")
            return None
        return "A"

    def _record_error(self, message: str) -> None:
        self.errors.append(f"Line {self._row_count}: {message}")
