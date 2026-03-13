from __future__ import annotations

from pathlib import Path

from .exceptions import SequenceHintError
from .structures import PDBStructure


def load_sequence_hints(
    hint: str | Path | None,
    structure: PDBStructure,
    *,
    label: str,
) -> dict[str, str] | None:
    if hint is None:
        return None

    text = _read_hint_text(hint)
    records = _parse_fasta(text) if text.lstrip().startswith(">") else [(None, _normalize_sequence(text))]
    if not records:
        raise SequenceHintError(f"{label} sequence hint did not contain any sequence records.")

    chains = list(structure.chain_records().items())
    if not chains:
        raise SequenceHintError(f"{label} structure does not contain any selected residues.")

    if len(records) == 1:
        return _assign_single_record(records[0][1], chains, label=label)
    return _assign_multiple_records(records, chains, label=label)


def _read_hint_text(hint: str | Path) -> str:
    path = Path(hint)
    if path.exists():
        return path.read_text(encoding="utf-8")
    return str(hint)


def _parse_fasta(text: str) -> list[tuple[str | None, str]]:
    records: list[tuple[str | None, str]] = []
    current_header: str | None = None
    current_chunks: list[str] = []
    for row in text.splitlines():
        stripped = row.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if current_header is not None:
                records.append((current_header, _normalize_sequence("".join(current_chunks))))
            current_header = stripped[1:].strip() or None
            current_chunks = []
            continue
        current_chunks.append(stripped)
    if current_header is not None:
        records.append((current_header, _normalize_sequence("".join(current_chunks))))
    return [(header, sequence) for header, sequence in records if sequence]


def _assign_single_record(
    sequence: str,
    chains: list[tuple[str, list[tuple[int, object]]]],
    *,
    label: str,
) -> dict[str, str]:
    if len(chains) == 1:
        chain_name, records = chains[0]
        _validate_sequence_length(sequence, len(records), label=label, chain=chain_name)
        return {chain_name: sequence}

    total_length = sum(len(records) for _, records in chains)
    if len(sequence) != total_length:
        raise SequenceHintError(
            f"{label} sequence hint length {len(sequence)} does not match the total selected residue count "
            f"{total_length} across chains {', '.join(chain for chain, _ in chains)}."
        )

    assigned: dict[str, str] = {}
    offset = 0
    for chain_name, records in chains:
        next_offset = offset + len(records)
        assigned[chain_name] = sequence[offset:next_offset]
        offset = next_offset
    return assigned


def _assign_multiple_records(
    records: list[tuple[str | None, str]],
    chains: list[tuple[str, list[tuple[int, object]]]],
    *,
    label: str,
) -> dict[str, str]:
    chain_names = [chain_name for chain_name, _ in chains]
    by_header = {_header_key(header): sequence for header, sequence in records if _header_key(header) is not None}

    if all(chain_name in by_header for chain_name in chain_names):
        assigned = {chain_name: by_header[chain_name] for chain_name in chain_names}
    elif len(records) == len(chains):
        assigned = {
            chain_name: sequence
            for (chain_name, _), (_, sequence) in zip(chains, records)
        }
    else:
        raise SequenceHintError(
            f"{label} FASTA records could not be mapped to structure chains {', '.join(chain_names)}. "
            "Use FASTA record identifiers that match chain IDs or provide one combined sequence."
        )

    for chain_name, chain_records in chains:
        _validate_sequence_length(assigned[chain_name], len(chain_records), label=label, chain=chain_name)
    return assigned


def _header_key(header: str | None) -> str | None:
    if not header:
        return None
    first_token = header.split()[0]
    if first_token.lower().startswith("chain:"):
        return first_token.split(":", 1)[1]
    return first_token


def _normalize_sequence(sequence: str) -> str:
    cleaned = "".join(char for char in sequence.upper() if char.isalpha())
    cleaned = cleaned.replace("T", "U")
    return cleaned


def _validate_sequence_length(sequence: str, residue_count: int, *, label: str, chain: str) -> None:
    if len(sequence) != residue_count:
        raise SequenceHintError(
            f"{label} sequence hint for chain '{chain}' has length {len(sequence)}, "
            f"but the selected structure residues for that chain have length {residue_count}."
        )
