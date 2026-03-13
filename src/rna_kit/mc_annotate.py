from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

from .exceptions import ToolExecutionError, ToolResolutionError
from .structures import prepare_external_structure_input
from .tools import default_tool_registry


@dataclass(frozen=True)
class RawInteraction:
    type: str
    chain_a: str
    pos_a: int
    nt_a: str
    chain_b: str
    pos_b: int
    nt_b: str
    extra_1: str
    extra_2: str
    extra_3: str


@dataclass(frozen=True)
class MCAnnotateResult:
    residues: list[tuple[str, str, str]]
    interactions: list[RawInteraction]


class MCAnnotateRunner:
    _PAIR_PATTERN = re.compile(
        r"^([A-Z]|'[0-9]'|)(\d+)-([A-Z]|'[0-9]'|)(\d+) : (\w+)-(\w+) ([\w']+)/([\w']+)"
        r"(?:.*)pairing( (parallel|antiparallel) (cis|trans))"
    )
    _STACK_PATTERN = re.compile(
        r"^([A-Z]|'[0-9]'|)(\d+)-([A-Z]|'[0-9]'|)(\d+) :.*(inward|upward|downward|outward).*"
    )

    def __init__(
        self,
        binary_path: str | Path | None = None,
        cache_dir: str | Path | None = None,
        annotation_overrides: dict[str | Path, str | Path] | None = None,
    ):
        self.binary_path = Path(binary_path) if binary_path else None
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self._cache: dict[Path, MCAnnotateResult] = {}
        self.annotation_overrides = {
            _normalize_path(pdb_file): Path(annotation_file)
            for pdb_file, annotation_file in (annotation_overrides or {}).items()
        }

    def cache_key(self, pdb_file: str | Path) -> str:
        return f"mc_annotate::{self.annotation_path_for(pdb_file)}"

    def annotation_path_for(self, pdb_file: str | Path) -> Path:
        pdb_path = _normalize_path(pdb_file)
        overridden = self.annotation_overrides.get(pdb_path)
        if overridden is not None:
            return overridden
        cache_dir = self.cache_dir or pdb_path.parent
        return cache_dir / f"{pdb_path.name}.mcout"

    def resolve_binary(self) -> Path | None:
        registry = default_tool_registry()
        status = registry.status("mc_annotate", override=self.binary_path)
        return None if status.binary_path is None else Path(status.binary_path)

    def load(self, pdb_file: str | Path) -> MCAnnotateResult:
        annotation_path = self.annotation_path_for(pdb_file)
        cached = self._cache.get(annotation_path)
        if cached is not None:
            return cached

        if not annotation_path.exists():
            self._generate_annotation(pdb_file, annotation_path)

        result = self.parse(annotation_path)
        self._cache[annotation_path] = result
        return result

    def indexed_interactions(self, structure) -> list[tuple[str, int, int, str]]:
        cache_key = self.cache_key(structure.pdb_file)
        cached = structure.cached_interactions(cache_key)
        if cached is not None:
            return cached

        annotation = self.load(structure.pdb_file)
        indexed: list[tuple[str, int, int, str]] = []
        for interaction in annotation.interactions:
            rank_a = structure.rank_of(interaction.chain_a, interaction.pos_a)
            rank_b = structure.rank_of(interaction.chain_b, interaction.pos_b)
            if rank_a is None or rank_b is None:
                continue

            if interaction.type == "STACK":
                extra = interaction.extra_1
            else:
                extra = f"{interaction.extra_1}{interaction.extra_2}"
            indexed.append((interaction.type, min(rank_a, rank_b), max(rank_a, rank_b), extra))

        structure.set_cached_interactions(cache_key, indexed)
        return indexed

    def parse(self, annotation_path: str | Path) -> MCAnnotateResult:
        residues: list[tuple[str, str, str]] = []
        interactions: list[RawInteraction] = []
        state = "out"
        model_count = 0

        for raw_line in Path(annotation_path).read_text(encoding="utf-8").splitlines():
            line = raw_line.strip()

            if line.startswith("Residue conformations"):
                if model_count == 0:
                    state = "residue"
                    model_count += 1
                    continue
                break

            if line.startswith("Base-pairs"):
                state = "pair"
                continue

            if line.startswith("Adjacent stackings") or line.startswith("Non-Adjacent stackings"):
                state = "stack"
                continue

            if line.endswith("----------"):
                state = "out"
                continue

            interaction: RawInteraction | None = None
            if state == "residue":
                parts = line.split()
                if len(parts) == 5:
                    residues.append((parts[0][0], parts[0][1:], parts[2]))
            elif state == "pair":
                match = self._PAIR_PATTERN.match(line)
                if match:
                    interaction = self._convert_pair(match.groups())
            elif state == "stack":
                match = self._STACK_PATTERN.match(line)
                if match:
                    interaction = self._convert_stack(match.groups())

            if interaction is not None:
                interactions.append(interaction)

        return MCAnnotateResult(residues=residues, interactions=interactions)

    def _generate_annotation(self, pdb_file: str | Path, annotation_path: Path) -> None:
        binary = self.resolve_binary()
        if binary is None:
            raise ToolResolutionError(
                "MC-Annotate is not available. Provide a precomputed '.mcout' file or set "
                "'RNA_KIT_MC_ANNOTATE'."
            )
        annotation_path.parent.mkdir(parents=True, exist_ok=True)
        source_path = Path(pdb_file)
        try:
            if source_path.suffix.lower() in {".cif", ".mmcif"}:
                with TemporaryDirectory(prefix="rna-kit-mc-annotate-") as temp_dir:
                    prepared_input = prepare_external_structure_input(
                        source_path,
                        Path(temp_dir) / f"{source_path.stem}.pdb",
                    )
                    result = subprocess.run(
                        [str(binary), str(prepared_input)],
                        check=True,
                        capture_output=True,
                        text=True,
                    )
            else:
                result = subprocess.run(
                    [str(binary), str(source_path)],
                    check=True,
                    capture_output=True,
                    text=True,
                )
        except OSError as exc:
            raise ToolExecutionError(
                f"Failed to execute MC-Annotate at '{binary}'. This repository currently bundles "
                "a Linux binary, so macOS users must provide a compatible executable or reuse "
                "precomputed '.mcout' files."
            ) from exc
        except subprocess.CalledProcessError as exc:
            raise ToolExecutionError(
                f"MC-Annotate failed for '{source_path}': {exc.stderr.strip() or exc.stdout.strip()}"
            ) from exc
        annotation_path.write_text(result.stdout, encoding="utf-8")

    def _convert_pair(self, groups: tuple[str, ...]) -> RawInteraction | None:
        int_a = groups[6][0].upper()
        int_b = groups[7][0].upper()

        if int_a not in {"W", "H", "S"} or int_b not in {"W", "H", "S"}:
            return None

        chain_a = groups[0].replace("'", "")
        pos_a = int(groups[1])
        nt_a = groups[4]

        chain_b = groups[2].replace("'", "")
        pos_b = int(groups[3])
        nt_b = groups[5]

        interaction_type = groups[10].lower()
        pair_name = "PAIR_2D" if f"{interaction_type}{int_a}{int_b}" == "cisWW" else "PAIR_3D"

        if ((chain_a == chain_b) and (pos_a < pos_b)) or (chain_a < chain_b):
            return RawInteraction(
                type=pair_name,
                chain_a=chain_a,
                pos_a=pos_a,
                nt_a=nt_a,
                chain_b=chain_b,
                pos_b=pos_b,
                nt_b=nt_b,
                extra_1=f"{int_a}{int_b}",
                extra_2=interaction_type,
                extra_3="",
            )
        return RawInteraction(
            type=pair_name,
            chain_a=chain_b,
            pos_a=pos_b,
            nt_a=nt_b,
            chain_b=chain_a,
            pos_b=pos_a,
            nt_b=nt_a,
            extra_1=f"{int_a}{int_b}",
            extra_2=interaction_type,
            extra_3="",
        )

    def _convert_stack(self, groups: tuple[str, ...]) -> RawInteraction:
        return RawInteraction(
            type="STACK",
            chain_a=groups[0].replace("'", ""),
            pos_a=int(groups[1]),
            nt_a="",
            chain_b=groups[2].replace("'", ""),
            pos_b=int(groups[3]),
            nt_b="",
            extra_1=groups[4],
            extra_2="",
            extra_3="",
        )


class MCAnnotate(MCAnnotateRunner):
    """Compatibility alias for the legacy class name."""


def existing_annotation_path(
    pdb_file: str | Path,
    *,
    explicit_annotation: str | Path | None = None,
    cache_dir: str | Path | None = None,
) -> Path | None:
    if explicit_annotation is not None:
        path = Path(explicit_annotation)
        return path if path.exists() else None

    pdb_path = _normalize_path(pdb_file)
    candidate_dir = Path(cache_dir) if cache_dir is not None else pdb_path.parent
    candidate = candidate_dir / f"{pdb_path.name}.mcout"
    if candidate.exists():
        return candidate
    fallback = pdb_path.parent / f"{pdb_path.name}.mcout"
    if fallback.exists():
        return fallback
    return None


def clone_with_annotation_overrides(
    runner: MCAnnotateRunner | None,
    overrides: dict[str | Path, str | Path],
) -> MCAnnotateRunner | None:
    if runner is None:
        return None
    if not hasattr(runner, "annotation_overrides"):
        return runner
    merged = dict(runner.annotation_overrides)
    merged.update({_normalize_path(key): Path(value) for key, value in overrides.items()})
    return MCAnnotateRunner(
        binary_path=getattr(runner, "binary_path", None),
        cache_dir=getattr(runner, "cache_dir", None),
        annotation_overrides=merged or None,
    )


def _normalize_path(path: str | Path) -> Path:
    return Path(path).expanduser().resolve()
