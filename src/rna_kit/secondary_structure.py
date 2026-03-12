from __future__ import annotations

import html
import platform
import re
import stat
import subprocess
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from .exceptions import MetricCalculationError, ReportGenerationError, ToolExecutionError, ToolResolutionError
from .mc_annotate import MCAnnotateResult, MCAnnotateRunner, RawInteraction
from .structures import PDBStructure
from .tools import default_tool_registry


@dataclass(frozen=True)
class CSSRBasePairRecord:
    chain_1: str | None
    pos_1: int
    nt_1: str
    chain_2: str | None
    pos_2: int
    nt_2: str
    pair: str
    classification: str
    saenger: str | None
    leontis_westhof: str | None
    dssr: str | None


@dataclass(frozen=True)
class CSSRAnnotation:
    base_pairs: tuple[CSSRBasePairRecord, ...]


@dataclass(frozen=True)
class SecondaryStructureBasePair:
    rank_1: int
    rank_2: int
    chain_1: str
    pos_1: int
    nt_1: str
    chain_2: str
    pos_2: int
    nt_2: str
    pair: str
    classification: str
    saenger: str | None
    leontis_westhof: str | None
    dssr: str | None


@dataclass(frozen=True)
class SecondaryStructureResult:
    sequence: str
    dot_bracket: str
    selected_residues: int
    base_pair_count: int
    base_pairs: tuple[SecondaryStructureBasePair, ...]


@dataclass(frozen=True)
class SecondaryStructureComparisonResult:
    native: SecondaryStructureResult
    prediction: SecondaryStructureResult
    true_positives: int
    false_positives: int
    false_negatives: int
    precision: float
    recall: float
    f1: float
    jaccard: float
    true_positive_pairs: tuple["SecondaryStructureComparisonPair", ...] = ()
    false_positive_pairs: tuple["SecondaryStructureComparisonPair", ...] = ()
    false_negative_pairs: tuple["SecondaryStructureComparisonPair", ...] = ()
    pair_details: tuple["SecondaryStructureComparisonPair", ...] = ()


@dataclass(frozen=True)
class SecondaryStructureComparisonPair:
    status: str
    rank_1: int
    rank_2: int
    native_pair: SecondaryStructureBasePair | None
    prediction_pair: SecondaryStructureBasePair | None


@dataclass(frozen=True)
class SecondaryStructureVisualizationResult:
    output: str
    format: str


class CSSRRunner:
    _PAIR_LINE = re.compile(
        r"^\s*(?P<count>\d+)\s+"
        r"(?P<nt1>\S+)\s+"
        r"(?P<nt2>\S+)\s+"
        r"(?P<pair>[A-Z]-[A-Z])\s+"
        r"(?P<classification>\S+)"
        r"(?:\s+(?P<saenger>\S+)\s+(?P<lw>\S+)\s+(?P<dssr>\S+))?\s*$"
    )
    _NT_TOKEN = re.compile(
        r"^(?:(?P<chain>[^.]+)\.)?(?P<nt>[A-Za-z])(?P<pos>-?\d+)(?:\^(?P<icode>.+))?$"
    )
    _DOWNLOAD_BASE_URL = "https://raw.githubusercontent.com/pylelab/CSSR/master/exe"

    def __init__(
        self,
        binary_path: str | Path | None = None,
        cache_dir: str | Path | None = None,
        download_dir: str | Path | None = None,
        auto_download: bool = True,
        annotation_overrides: dict[str | Path, str | Path] | None = None,
    ) -> None:
        self.binary_path = Path(binary_path) if binary_path else None
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.download_dir = Path(download_dir) if download_dir else None
        self.auto_download = auto_download
        self._cache: dict[Path, CSSRAnnotation] = {}
        self.annotation_overrides = {
            _normalize_path(pdb_file): Path(annotation_file)
            for pdb_file, annotation_file in (annotation_overrides or {}).items()
        }

    def annotation_path_for(self, pdb_file: str | Path) -> Path:
        pdb_path = _normalize_path(pdb_file)
        overridden = self.annotation_overrides.get(pdb_path)
        if overridden is not None:
            return overridden
        cache_dir = self.cache_dir or pdb_path.parent
        return cache_dir / f"{pdb_path.name}.cssr"

    def resolve_binary(self) -> Path:
        registry = default_tool_registry()
        status = registry.status("cssr", override=self.binary_path)
        if status.binary_path is not None:
            return Path(status.binary_path)

        platform_binary = _platform_binary_name()
        cache_candidate = self._binary_cache_dir() / platform_binary
        if cache_candidate.is_file():
            return cache_candidate

        if not self.auto_download:
            raise ToolResolutionError(
                "CSSR is not available. Set 'RNA_KIT_CSSR', pass '--cssr', "
                "or enable automatic download."
            )

        return self._download_binary(platform_binary)

    def load(self, pdb_file: str | Path) -> CSSRAnnotation:
        annotation_path = self.annotation_path_for(pdb_file)
        cached = self._cache.get(annotation_path)
        if cached is not None:
            return cached

        if not annotation_path.exists():
            self._generate_annotation(pdb_file, annotation_path)

        result = self.parse(annotation_path)
        self._cache[annotation_path] = result
        return result

    def parse(self, annotation_path: str | Path) -> CSSRAnnotation:
        rows = Path(annotation_path).read_text(encoding="utf-8").splitlines()
        base_pairs: list[CSSRBasePairRecord] = []
        in_table = False
        found_table = False

        for row in rows:
            stripped = row.strip()
            if stripped.startswith("List of ") and stripped.endswith(" base pairs"):
                in_table = True
                found_table = True
                continue
            if not in_table:
                continue
            if not stripped or stripped.startswith("*") or stripped.startswith("nt1"):
                continue

            match = self._PAIR_LINE.match(row)
            if not match:
                break

            nt_1 = self._parse_nt_token(match.group("nt1"))
            nt_2 = self._parse_nt_token(match.group("nt2"))
            base_pairs.append(
                CSSRBasePairRecord(
                    chain_1=nt_1[0],
                    pos_1=nt_1[2],
                    nt_1=nt_1[1],
                    chain_2=nt_2[0],
                    pos_2=nt_2[2],
                    nt_2=nt_2[1],
                    pair=match.group("pair"),
                    classification=match.group("classification"),
                    saenger=match.group("saenger"),
                    leontis_westhof=match.group("lw"),
                    dssr=match.group("dssr"),
                )
            )

        if not found_table:
            raise MetricCalculationError(f"CSSR output '{annotation_path}' did not contain a base-pair table.")

        return CSSRAnnotation(base_pairs=tuple(base_pairs))

    def _generate_annotation(self, pdb_file: str | Path, annotation_path: Path) -> None:
        binary = self.resolve_binary()
        annotation_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            subprocess.run(
                [str(binary), str(pdb_file), str(annotation_path), "-o", "2"],
                check=True,
                capture_output=True,
                text=True,
            )
        except OSError as exc:
            raise ToolExecutionError(f"Failed to execute CSSR at '{binary}'.") from exc
        except subprocess.CalledProcessError as exc:
            raise ToolExecutionError(
                f"CSSR failed for '{pdb_file}': {exc.stderr.strip() or exc.stdout.strip()}"
            ) from exc

    def _download_binary(self, binary_name: str) -> Path:
        target_dir = self._binary_cache_dir()
        target_dir.mkdir(parents=True, exist_ok=True)
        target_path = target_dir / binary_name
        url = f"{self._DOWNLOAD_BASE_URL}/{binary_name}"
        temp_path = target_path.with_name(f"{target_path.name}.tmp")

        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                temp_path.write_bytes(response.read())
        except OSError as exc:
            raise ToolResolutionError(
                f"Unable to download CSSR from '{url}'. Set 'RNA_KIT_CSSR' or pass '--cssr'."
            ) from exc

        temp_path.chmod(temp_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        temp_path.replace(target_path)
        return target_path

    def _binary_cache_dir(self) -> Path:
        if self.download_dir is not None:
            return self.download_dir
        return Path.home() / ".cache" / "rna-kit" / "bin"

    def _parse_nt_token(self, token: str) -> tuple[str | None, str, int]:
        match = self._NT_TOKEN.match(token)
        if match is None:
            raise MetricCalculationError(f"Unrecognized CSSR nucleotide token '{token}'.")
        return match.group("chain"), match.group("nt").upper(), int(match.group("pos"))


def calculate_secondary_structure_for_structure(
    structure: PDBStructure,
    runner: MCAnnotateRunner | None = None,
) -> SecondaryStructureResult:
    runner = runner or MCAnnotateRunner()
    annotation = runner.load(structure.pdb_file)
    base_pairs = _base_pairs_from_mc_annotate(structure, annotation)
    base_pairs.sort(key=lambda item: (item.rank_1, item.rank_2))
    return SecondaryStructureResult(
        sequence=structure.raw_sequence(),
        dot_bracket=_pairs_to_dot_bracket(
            len(structure.res_seq),
            [(pair.rank_1, pair.rank_2) for pair in base_pairs],
        ),
        selected_residues=len(structure.res_seq),
        base_pair_count=len(base_pairs),
        base_pairs=tuple(base_pairs),
    )


def compare_secondary_structures(
    native_structure: PDBStructure,
    prediction_structure: PDBStructure,
    runner: MCAnnotateRunner | None = None,
) -> SecondaryStructureComparisonResult:
    runner = runner or MCAnnotateRunner()
    native_result = calculate_secondary_structure_for_structure(native_structure, runner=runner)
    prediction_result = calculate_secondary_structure_for_structure(prediction_structure, runner=runner)

    native_pair_map = {(pair.rank_1, pair.rank_2): pair for pair in native_result.base_pairs}
    prediction_pair_map = {(pair.rank_1, pair.rank_2): pair for pair in prediction_result.base_pairs}
    native_pairs = set(native_pair_map)
    prediction_pairs = set(prediction_pair_map)
    true_positives = len(native_pairs & prediction_pairs)
    false_positives = len(prediction_pairs - native_pairs)
    false_negatives = len(native_pairs - prediction_pairs)
    true_positive_details = tuple(
        SecondaryStructureComparisonPair(
            status="tp",
            rank_1=rank_1,
            rank_2=rank_2,
            native_pair=native_pair_map[(rank_1, rank_2)],
            prediction_pair=prediction_pair_map[(rank_1, rank_2)],
        )
        for rank_1, rank_2 in sorted(native_pairs & prediction_pairs)
    )
    false_positive_details = tuple(
        SecondaryStructureComparisonPair(
            status="fp",
            rank_1=rank_1,
            rank_2=rank_2,
            native_pair=None,
            prediction_pair=prediction_pair_map[(rank_1, rank_2)],
        )
        for rank_1, rank_2 in sorted(prediction_pairs - native_pairs)
    )
    false_negative_details = tuple(
        SecondaryStructureComparisonPair(
            status="fn",
            rank_1=rank_1,
            rank_2=rank_2,
            native_pair=native_pair_map[(rank_1, rank_2)],
            prediction_pair=None,
        )
        for rank_1, rank_2 in sorted(native_pairs - prediction_pairs)
    )

    if not native_pairs and not prediction_pairs:
        precision = 1.0
        recall = 1.0
        f1 = 1.0
        jaccard = 1.0
    else:
        precision = true_positives / len(prediction_pairs) if prediction_pairs else 0.0
        recall = true_positives / len(native_pairs) if native_pairs else 0.0
        f1 = (2.0 * precision * recall / (precision + recall)) if precision + recall else 0.0
        union = native_pairs | prediction_pairs
        jaccard = true_positives / len(union) if union else 1.0

    return SecondaryStructureComparisonResult(
        native=native_result,
        prediction=prediction_result,
        true_positives=true_positives,
        false_positives=false_positives,
        false_negatives=false_negatives,
        precision=precision,
        recall=recall,
        f1=f1,
        jaccard=jaccard,
        true_positive_pairs=true_positive_details,
        false_positive_pairs=false_positive_details,
        false_negative_pairs=false_negative_details,
        pair_details=tuple(
            sorted(
                (*true_positive_details, *false_positive_details, *false_negative_details),
                key=lambda item: (item.rank_1, item.rank_2, item.status),
            )
        ),
    )


def _base_pairs_from_mc_annotate(
    structure: PDBStructure,
    annotation: MCAnnotateResult,
) -> list[SecondaryStructureBasePair]:
    base_pairs: list[SecondaryStructureBasePair] = []
    seen_pairs: set[tuple[int, int]] = set()

    for interaction in annotation.interactions:
        if interaction.type != "PAIR_2D":
            continue

        rank_1 = structure.rank_of(interaction.chain_a, interaction.pos_a)
        rank_2 = structure.rank_of(interaction.chain_b, interaction.pos_b)
        if rank_1 is None or rank_2 is None or rank_1 == rank_2:
            continue

        normalized = _normalize_mc_annotate_pair(interaction, rank_1, rank_2)
        pair_key = (normalized.rank_1, normalized.rank_2)
        if pair_key in seen_pairs:
            continue
        seen_pairs.add(pair_key)
        base_pairs.append(normalized)

    return base_pairs


def _normalize_mc_annotate_pair(
    interaction: RawInteraction,
    rank_1: int,
    rank_2: int,
) -> SecondaryStructureBasePair:
    if rank_1 <= rank_2:
        return SecondaryStructureBasePair(
            rank_1=rank_1,
            rank_2=rank_2,
            chain_1=interaction.chain_a,
            pos_1=interaction.pos_a,
            nt_1=interaction.nt_a,
            chain_2=interaction.chain_b,
            pos_2=interaction.pos_b,
            nt_2=interaction.nt_b,
            pair=f"{interaction.nt_a}-{interaction.nt_b}",
            classification="cisWW",
            saenger=None,
            leontis_westhof="cWW",
            dssr=None,
        )
    return SecondaryStructureBasePair(
        rank_1=rank_2,
        rank_2=rank_1,
        chain_1=interaction.chain_b,
        pos_1=interaction.pos_b,
        nt_1=interaction.nt_b,
        chain_2=interaction.chain_a,
        pos_2=interaction.pos_a,
        nt_2=interaction.nt_a,
        pair=f"{interaction.nt_a}-{interaction.nt_b}",
        classification="cisWW",
        saenger=None,
        leontis_westhof="cWW",
        dssr=None,
    )


def render_secondary_structure_svg(
    result: SecondaryStructureResult,
    title: str | None = None,
) -> str:
    width = _svg_width(result.selected_residues)
    left_margin = 72
    right_margin = 36
    top_margin = 64
    track_y = 170
    footer_y = track_y + 112
    height = footer_y + 70
    positions = _x_positions(result.selected_residues, width, left_margin, right_margin)

    parts = [
        _svg_header(width, height),
        _svg_style_block(),
        f'<text class="title" x="{width / 2:.1f}" y="34" text-anchor="middle">{html.escape(title or "Secondary Structure")}</text>',
        f'<text class="subtitle" x="{width / 2:.1f}" y="56" text-anchor="middle">base pairs={result.base_pair_count} residues={result.selected_residues}</text>',
        _render_track(
            result,
            positions=positions,
            baseline_y=track_y,
            title="Structure",
            pair_colors={(pair.rank_1, pair.rank_2): "#2f6db3" for pair in result.base_pairs},
            pair_labels={
                (pair.rank_1, pair.rank_2): _format_base_pair_label(pair, status="pair")
                for pair in result.base_pairs
            },
            title_y=96,
        ),
    ]
    parts.extend(_render_sequence_block(result.sequence, result.dot_bracket, left_margin, footer_y))
    parts.append("</svg>")
    return "\n".join(parts)


def render_secondary_structure_comparison_svg(
    result: SecondaryStructureComparisonResult,
    title: str | None = None,
) -> str:
    residue_count = max(result.native.selected_residues, result.prediction.selected_residues)
    width = _svg_width(residue_count)
    left_margin = 72
    right_margin = 36
    top_margin = 84
    native_y = 170
    prediction_y = 360
    footer_y = prediction_y + 118
    height = footer_y + 96
    positions = _x_positions(residue_count, width, left_margin, right_margin)

    native_pairs = {(pair.rank_1, pair.rank_2) for pair in result.native.base_pairs}
    prediction_pairs = {(pair.rank_1, pair.rank_2) for pair in result.prediction.base_pairs}
    shared_pairs = native_pairs & prediction_pairs

    native_colors = {
        pair: "#208b3a" if pair in shared_pairs else "#2f6db3"
        for pair in native_pairs
    }
    prediction_colors = {
        pair: "#208b3a" if pair in shared_pairs else "#d97706"
        for pair in prediction_pairs
    }
    native_labels = {
        (detail.rank_1, detail.rank_2): _format_comparison_pair_label(detail)
        for detail in (*result.true_positive_pairs, *result.false_negative_pairs)
    }
    prediction_labels = {
        (detail.rank_1, detail.rank_2): _format_comparison_pair_label(detail)
        for detail in (*result.true_positive_pairs, *result.false_positive_pairs)
    }

    parts = [
        _svg_header(width, height),
        _svg_style_block(),
        f'<text class="title" x="{width / 2:.1f}" y="34" text-anchor="middle">{html.escape(title or "Secondary Structure Comparison")}</text>',
        (
            f'<text class="subtitle" x="{width / 2:.1f}" y="58" text-anchor="middle">'
            f'precision={result.precision:.3f} recall={result.recall:.3f} '
            f'F1={result.f1:.3f} jaccard={result.jaccard:.3f}'
            "</text>"
        ),
        _render_legend(width),
        _render_track(
            result.native,
            positions=positions,
            baseline_y=native_y,
            title="Native",
            pair_colors=native_colors,
            pair_labels=native_labels,
            title_y=112,
        ),
        _render_track(
            result.prediction,
            positions=positions,
            baseline_y=prediction_y,
            title="Prediction",
            pair_colors=prediction_colors,
            pair_labels=prediction_labels,
            title_y=302,
        ),
    ]
    parts.extend(
        _render_sequence_block(
            result.native.sequence,
            result.native.dot_bracket,
            left_margin,
            footer_y,
            label="Native",
        )
    )
    parts.extend(
        _render_sequence_block(
            result.prediction.sequence,
            result.prediction.dot_bracket,
            left_margin,
            footer_y + 32,
            label="Prediction",
        )
    )
    parts.append("</svg>")
    return "\n".join(parts)


def write_secondary_structure_svg(
    result: SecondaryStructureResult,
    output_file: str | Path,
    title: str | None = None,
) -> Path:
    path = Path(output_file)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(render_secondary_structure_svg(result, title=title), encoding="utf-8")
    except OSError as exc:
        raise ReportGenerationError(f"Failed to write secondary-structure SVG to '{path}'.") from exc
    return path


def write_secondary_structure_comparison_svg(
    result: SecondaryStructureComparisonResult,
    output_file: str | Path,
    title: str | None = None,
) -> Path:
    path = Path(output_file)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(render_secondary_structure_comparison_svg(result, title=title), encoding="utf-8")
    except OSError as exc:
        raise ReportGenerationError(f"Failed to write secondary-structure comparison SVG to '{path}'.") from exc
    return path


def _normalize_pair_orientation(
    pair: CSSRBasePairRecord,
    chain_1: str,
    chain_2: str,
    rank_1: int,
    rank_2: int,
) -> SecondaryStructureBasePair:
    if rank_1 <= rank_2:
        return SecondaryStructureBasePair(
            rank_1=rank_1,
            rank_2=rank_2,
            chain_1=chain_1,
            pos_1=pair.pos_1,
            nt_1=pair.nt_1,
            chain_2=chain_2,
            pos_2=pair.pos_2,
            nt_2=pair.nt_2,
            pair=pair.pair,
            classification=pair.classification,
            saenger=pair.saenger,
            leontis_westhof=pair.leontis_westhof,
            dssr=pair.dssr,
        )
    return SecondaryStructureBasePair(
        rank_1=rank_2,
        rank_2=rank_1,
        chain_1=chain_2,
        pos_1=pair.pos_2,
        nt_1=pair.nt_2,
        chain_2=chain_1,
        pos_2=pair.pos_1,
        nt_2=pair.nt_1,
        pair=pair.pair,
        classification=pair.classification,
        saenger=pair.saenger,
        leontis_westhof=pair.leontis_westhof,
        dssr=pair.dssr,
    )


def _pairs_to_dot_bracket(length: int, pairs: list[tuple[int, int]]) -> str:
    characters = ["."] * length
    levels: list[list[tuple[int, int]]] = []

    for pair in sorted(pairs):
        level = _select_level(pair, levels)
        open_char, close_char = _dot_bracket_symbols(level)
        characters[pair[0]] = open_char
        characters[pair[1]] = close_char

    return "".join(characters)


def _select_level(pair: tuple[int, int], levels: list[list[tuple[int, int]]]) -> int:
    for index, intervals in enumerate(levels):
        if all(not _pairs_cross(pair, existing) for existing in intervals):
            intervals.append(pair)
            return index
    levels.append([pair])
    return len(levels) - 1


def _pairs_cross(left: tuple[int, int], right: tuple[int, int]) -> bool:
    return (left[0] < right[0] < left[1] < right[1]) or (right[0] < left[0] < right[1] < left[1])


def _dot_bracket_symbols(level: int) -> tuple[str, str]:
    symbols = [("(", ")"), ("[", "]"), ("{", "}"), ("<", ">")]
    if level < len(symbols):
        return symbols[level]

    alpha_index = level - len(symbols)
    if alpha_index >= 26:
        raise MetricCalculationError("The secondary structure contains more pseudoknot levels than supported.")
    return chr(ord("A") + alpha_index), chr(ord("a") + alpha_index)


def _default_chain_id(structure: PDBStructure) -> str | None:
    chain_ids = {record.chain for record in structure.res_list}
    if len(chain_ids) == 1:
        return next(iter(chain_ids))
    return None


def _platform_binary_name() -> str:
    system = platform.system().lower()
    machine = platform.machine().lower()
    if system == "darwin":
        return "CSSR-Mac64"
    if system == "linux" and machine in {"x86_64", "amd64"}:
        return "CSSR-Linux64"
    if system == "windows" and machine in {"x86_64", "amd64"}:
        return "CSSR-Win64.exe"
    raise ToolResolutionError(
        f"CSSR does not publish a precompiled binary for '{platform.system()} {platform.machine()}'. "
        "Provide a compatible executable with '--cssr' or 'RNA_KIT_CSSR'."
    )

def _normalize_path(path: str | Path) -> Path:
    return Path(path).expanduser().resolve()


def _svg_width(residue_count: int) -> int:
    return max(880, min(2200, 120 + residue_count * 14))


def _svg_header(width: int, height: int) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img">'
    )


def _svg_style_block() -> str:
    return """
<style>
  .title { font: 700 22px Helvetica, Arial, sans-serif; fill: #18212f; }
  .subtitle { font: 400 13px Helvetica, Arial, sans-serif; fill: #42526b; }
  .track-title { font: 700 14px Helvetica, Arial, sans-serif; fill: #18212f; }
  .axis { stroke: #7b8798; stroke-width: 1.5; }
  .tick { stroke: #b9c2cf; stroke-width: 1; }
  .tick-label { font: 12px Helvetica, Arial, sans-serif; fill: #5c6675; }
  .legend { font: 12px Helvetica, Arial, sans-serif; fill: #42526b; }
  .seq-label { font: 700 12px Helvetica, Arial, sans-serif; fill: #18212f; }
  .mono { font: 12px 'SFMono-Regular', Menlo, Consolas, monospace; fill: #253041; }
</style>
""".strip()


def _x_positions(residue_count: int, width: int, left_margin: int, right_margin: int) -> list[float]:
    if residue_count <= 1:
        return [float(left_margin)]
    usable = width - left_margin - right_margin
    return [left_margin + usable * index / (residue_count - 1) for index in range(residue_count)]


def _render_track(
    result: SecondaryStructureResult,
    positions: list[float],
    baseline_y: int,
    title: str,
    pair_colors: dict[tuple[int, int], str],
    pair_labels: dict[tuple[int, int], str],
    title_y: int,
) -> str:
    parts = [f'<text class="track-title" x="{positions[0]:.1f}" y="{title_y}">{html.escape(title)}</text>']
    parts.append(
        f'<line class="axis" x1="{positions[0]:.1f}" y1="{baseline_y}" x2="{positions[-1]:.1f}" y2="{baseline_y}" />'
    )
    parts.extend(_render_ticks(positions, baseline_y))
    for pair in result.base_pairs:
        pair_key = (pair.rank_1, pair.rank_2)
        parts.append(
            _render_arc(
                positions[pair.rank_1],
                positions[pair.rank_2],
                baseline_y,
                pair_colors[pair_key],
                pair_labels.get(pair_key),
            )
        )
    return "\n".join(parts)


def _render_arc(x1: float, x2: float, baseline_y: int, color: str, label: str | None) -> str:
    span = abs(x2 - x1)
    height = max(24.0, min(132.0, span * 0.42))
    mid_x = (x1 + x2) / 2.0
    control_y = baseline_y - height
    title = f"<title>{html.escape(label)}</title>" if label else ""
    return (
        "<path "
        f'd="M {x1:.2f} {baseline_y} Q {mid_x:.2f} {control_y:.2f} {x2:.2f} {baseline_y}" '
        f'fill="none" stroke="{color}" stroke-width="2.2" stroke-linecap="round" opacity="0.95">'
        f"{title}</path>"
    )


def _render_ticks(positions: list[float], baseline_y: int) -> list[str]:
    if not positions:
        return []
    tick_step = max(1, len(positions) // 12)
    tick_rows: list[str] = []
    for index, position in enumerate(positions, start=1):
        if index != 1 and index != len(positions) and (index - 1) % tick_step != 0:
            continue
        tick_rows.append(
            f'<line class="tick" x1="{position:.2f}" y1="{baseline_y - 6}" x2="{position:.2f}" y2="{baseline_y + 6}" />'
        )
        tick_rows.append(
            f'<text class="tick-label" x="{position:.2f}" y="{baseline_y + 24}" text-anchor="middle">{index}</text>'
        )
    return tick_rows


def _render_sequence_block(
    sequence: str,
    dot_bracket: str,
    left_margin: int,
    y: int,
    label: str | None = None,
) -> list[str]:
    prefix = f"{label}: " if label else ""
    return [
        f'<text class="seq-label" x="{left_margin}" y="{y}">{html.escape(prefix + "sequence")}</text>',
        f'<text class="mono" x="{left_margin + 76}" y="{y}">{html.escape(sequence)}</text>',
        f'<text class="seq-label" x="{left_margin}" y="{y + 16}">{html.escape(prefix + "dot")}</text>',
        f'<text class="mono" x="{left_margin + 76}" y="{y + 16}">{html.escape(dot_bracket)}</text>',
    ]


def _render_legend(width: int) -> str:
    items = [
        ("#208b3a", "TP / shared base pair"),
        ("#2f6db3", "FN / native only"),
        ("#d97706", "FP / prediction only"),
    ]
    x = width - 320
    parts = []
    for index, (color, label) in enumerate(items):
        y = 84 + index * 18
        parts.append(f'<line x1="{x}" y1="{y}" x2="{x + 20}" y2="{y}" stroke="{color}" stroke-width="3" />')
        parts.append(f'<text class="legend" x="{x + 28}" y="{y + 4}">{html.escape(label)}</text>')
    return "\n".join(parts)


def _format_base_pair_label(pair: SecondaryStructureBasePair, status: str) -> str:
    return (
        f"{status.upper()}: {pair.chain_1}:{pair.pos_1}({pair.nt_1}) - "
        f"{pair.chain_2}:{pair.pos_2}({pair.nt_2}) {pair.classification}"
    )


def _format_comparison_pair_label(detail: SecondaryStructureComparisonPair) -> str:
    pair = detail.native_pair or detail.prediction_pair
    if pair is None:
        return detail.status.upper()
    return _format_base_pair_label(pair, detail.status)
