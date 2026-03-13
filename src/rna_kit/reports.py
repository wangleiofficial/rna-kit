from __future__ import annotations

import html
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from .benchmark import BenchmarkEntry, BenchmarkResult, describe_prepared_pair
from .exceptions import ReportGenerationError, SchemaValidationError
from .metrics import AssessmentResult, LDDTResult, PreparedStructurePair, ResidueAssessment
from .secondary_structure import (
    SecondaryStructureComparisonPair,
    SecondaryStructureComparisonResult,
)
from .secondary_structure_web import (
    render_secondary_structure_assets,
    render_secondary_structure_comparison_component,
)
from .tools import ToolStatus

SCHEMA_VERSION = "1.0.0"


@dataclass(frozen=True)
class ReportMetadata:
    report_type: str
    schema_version: str
    library_version: str
    generated_at: str


@dataclass(frozen=True)
class AssessmentReportDocument:
    metadata: ReportMetadata
    reference: str
    prediction: str
    mapping: BenchmarkEntry
    metrics: AssessmentResult
    tool_statuses: tuple[ToolStatus, ...]
    warnings: tuple[str, ...] = ()
    artifacts: dict[str, str] | None = None


@dataclass(frozen=True)
class SecondaryStructureReportDocument:
    metadata: ReportMetadata
    reference: str
    prediction: str
    comparison: SecondaryStructureComparisonResult
    tool_statuses: tuple[ToolStatus, ...]
    warnings: tuple[str, ...] = ()
    artifacts: dict[str, str] | None = None


@dataclass(frozen=True)
class BenchmarkReportDocument:
    metadata: ReportMetadata
    result: BenchmarkResult
    tool_statuses: tuple[ToolStatus, ...]
    warnings: tuple[str, ...] = ()
    artifacts: dict[str, str] | None = None


def build_assessment_report_document(
    prepared: PreparedStructurePair,
    assessment: AssessmentResult,
    reference: str | Path,
    prediction: str | Path,
    tool_statuses: tuple[ToolStatus, ...] = (),
    warnings: tuple[str, ...] = (),
    artifacts: dict[str, str] | None = None,
) -> AssessmentReportDocument:
    return AssessmentReportDocument(
        metadata=_metadata("assessment"),
        reference=str(reference),
        prediction=str(prediction),
        mapping=describe_prepared_pair(prepared, prediction_file=prediction, native_file=reference),
        metrics=assessment,
        tool_statuses=tool_statuses,
        warnings=warnings,
        artifacts=artifacts,
    )


def build_secondary_structure_report_document(
    comparison: SecondaryStructureComparisonResult,
    reference: str | Path,
    prediction: str | Path,
    tool_statuses: tuple[ToolStatus, ...] = (),
    warnings: tuple[str, ...] = (),
    artifacts: dict[str, str] | None = None,
) -> SecondaryStructureReportDocument:
    return SecondaryStructureReportDocument(
        metadata=_metadata("secondary_structure_comparison"),
        reference=str(reference),
        prediction=str(prediction),
        comparison=comparison,
        tool_statuses=tool_statuses,
        warnings=warnings,
        artifacts=artifacts,
    )


def build_benchmark_report_document(
    result: BenchmarkResult,
    tool_statuses: tuple[ToolStatus, ...] = (),
    warnings: tuple[str, ...] = (),
    artifacts: dict[str, str] | None = None,
) -> BenchmarkReportDocument:
    return BenchmarkReportDocument(
        metadata=_metadata("benchmark"),
        result=result,
        tool_statuses=tool_statuses,
        warnings=warnings,
        artifacts=artifacts,
    )


def write_report_json(
    document: AssessmentReportDocument | SecondaryStructureReportDocument | BenchmarkReportDocument,
    output_file: str | Path,
) -> Path:
    _validate_document(document)
    path = Path(output_file)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(document), indent=2), encoding="utf-8")
    except OSError as exc:
        raise ReportGenerationError(f"Failed to write JSON report to '{path}'.") from exc
    return path


def write_assessment_html_report(
    document: AssessmentReportDocument,
    output_file: str | Path,
) -> Path:
    _validate_document(document)
    html_report = _render_assessment_html(document)
    return _write_html(output_file, html_report)


def write_secondary_structure_html_report(
    document: SecondaryStructureReportDocument,
    output_file: str | Path,
) -> Path:
    _validate_document(document)
    html_report = _render_secondary_structure_html(document)
    return _write_html(output_file, html_report)


def write_benchmark_html_report(
    document: BenchmarkReportDocument,
    output_file: str | Path,
) -> Path:
    _validate_document(document)
    path = Path(output_file)
    detail_dir = path.parent / f"{path.stem}_reports"
    detail_links: dict[int, str] = {}
    for index, entry in enumerate(document.result.entries, start=1):
        if entry.metrics is None:
            continue
        filename = f"{index:03d}-{_slugify(entry.label or Path(entry.prediction).stem)}.html"
        detail_path = detail_dir / filename
        detail_links[index] = str(Path(detail_dir.name) / filename)
        detail_html = _render_benchmark_entry_html(
            entry,
            metadata=document.metadata,
            tool_statuses=document.tool_statuses,
            dashboard_href=str(Path("..") / path.name),
        )
        _write_html(detail_path, detail_html)
    html_report = _render_benchmark_html(document, detail_links=detail_links)
    return _write_html(path, html_report)


def write_lddt_html_report(
    reference: str | Path,
    prediction: str | Path,
    result: LDDTResult,
    output_file: str | Path,
    title: str = "RNA Kit lDDT Report",
) -> Path:
    html_report = _render_lddt_html(reference, prediction, result, title)
    return _write_html(output_file, html_report)


def _render_assessment_html(
    document: AssessmentReportDocument,
) -> str:
    metrics = document.metrics
    summary_rows = _assessment_summary_rows(metrics)

    sections = [
        _html_header(
            "RNA Kit Assessment Report",
            document.metadata,
            include_secondary_assets=metrics.secondary_structure is not None,
        ),
        _html_identity_block(document.reference, document.prediction),
        _html_table("Summary Metrics", ("Metric", "Value"), summary_rows),
        _mapping_table(document.mapping),
        _tool_table(document.tool_statuses),
        _warnings_block(document.warnings),
    ]
    if metrics.per_residue is not None:
        sections.append(_per_residue_lddt_section(metrics.per_residue))
    if metrics.secondary_structure is not None:
        sections.append(_secondary_structure_section(metrics.secondary_structure))
    sections.append(_artifacts_block(document.artifacts))
    sections.append(_html_footer())
    return "\n".join(section for section in sections if section)


def _render_secondary_structure_html(
    document: SecondaryStructureReportDocument,
) -> str:
    comparison = document.comparison
    summary_rows = [
        ("Precision", f"{comparison.precision:.4f}"),
        ("Recall", f"{comparison.recall:.4f}"),
        ("F1", f"{comparison.f1:.4f}"),
        ("Jaccard", f"{comparison.jaccard:.4f}"),
        ("TP", str(comparison.true_positives)),
        ("FP", str(comparison.false_positives)),
        ("FN", str(comparison.false_negatives)),
    ]
    sections = [
        _html_header(
            "RNA Kit Secondary Structure Comparison",
            document.metadata,
            include_secondary_assets=True,
        ),
        _html_identity_block(document.reference, document.prediction),
        _html_table("Secondary-Structure Metrics", ("Metric", "Value"), summary_rows),
        _tool_table(document.tool_statuses),
        _warnings_block(document.warnings),
        _secondary_structure_section(comparison),
        _artifacts_block(document.artifacts),
        _html_footer(),
    ]
    return "\n".join(section for section in sections if section)


def _render_benchmark_html(
    document: BenchmarkReportDocument,
    *,
    detail_links: dict[int, str] | None = None,
) -> str:
    result = document.result
    succeeded = [entry for entry in result.entries if entry.metrics is not None]
    mean_rmsd = _mean([entry.metrics.rmsd for entry in succeeded])
    mean_lddt = _mean([entry.metrics.lddt for entry in succeeded])
    best_rmsd = min((entry.metrics.rmsd for entry in succeeded), default=None)
    best_lddt = max((entry.metrics.lddt for entry in succeeded), default=None)
    best_ss = max(
        (
            entry.metrics.secondary_structure_f1
            for entry in succeeded
            if entry.metrics.secondary_structure_f1 is not None
        ),
        default=None,
    )
    best_clash = min(
        (
            entry.metrics.molprobity.clashscore
            for entry in succeeded
            if entry.metrics.molprobity is not None and entry.metrics.molprobity.clashscore is not None
        ),
        default=None,
    )
    summary_cards = "".join(
        (
            _metric_card("Predictions", str(result.total_predictions)),
            _metric_card("Succeeded", str(result.succeeded)),
            _metric_card("Failed", str(result.failed)),
            _metric_card("Mean RMSD", _format_float(mean_rmsd)),
            _metric_card("Mean lDDT", _format_float(mean_lddt)),
            _metric_card("Best RMSD", _format_float(best_rmsd)),
            _metric_card("Best lDDT", _format_float(best_lddt)),
        )
        + ((_metric_card("Best SS F1", _format_float(best_ss)),) if best_ss is not None else ())
        + ((_metric_card("Best clashscore", _format_float(best_clash)),) if best_clash is not None else ())
    )

    rows = "".join(
        _benchmark_row(entry, detail_link=None if detail_links is None else detail_links.get(index))
        for index, entry in enumerate(result.entries, start=1)
    )
    sections = [
        _html_header("RNA Kit Benchmark Dashboard", document.metadata),
        (
            "<section><h2>Benchmark Summary</h2>"
            f"<p><strong>Reference:</strong> <code>{html.escape(result.reference or 'mixed')}</code></p>"
            f"<div class='metric-grid'>{summary_cards}</div>"
            "</section>"
        ),
        _tool_table(document.tool_statuses),
        _warnings_block(document.warnings),
        (
            "<section><h2>Results</h2><table>"
            "<thead><tr>"
            "<th>Label</th><th>Prediction</th><th>Status</th><th>Matched</th>"
            "<th>RMSD</th><th>INF_ALL</th><th>lDDT</th><th>SS F1</th>"
            "<th>Clashscore</th><th>MolProbity score</th><th>Error</th>"
            "</tr></thead>"
            f"<tbody>{rows}</tbody></table></section>"
        ),
        _artifacts_block(document.artifacts),
        _html_footer(),
    ]
    return "\n".join(section for section in sections if section)


def _render_benchmark_entry_html(
    entry: BenchmarkEntry,
    *,
    metadata: ReportMetadata,
    tool_statuses: tuple[ToolStatus, ...],
    dashboard_href: str,
) -> str:
    if entry.metrics is None:
        raise ReportGenerationError("Benchmark detail pages require a successful benchmark entry.")

    metrics = entry.metrics
    sections = [
        _html_header(
            f"RNA Kit Benchmark Entry: {entry.label or Path(entry.prediction).name}",
            metadata,
            include_secondary_assets=metrics.secondary_structure is not None,
        ),
        (
            "<section><h2>Navigation</h2>"
            f"<p><a href=\"{html.escape(dashboard_href)}\">Back to benchmark dashboard</a></p>"
            "</section>"
        ),
        _html_identity_block(entry.native, entry.prediction),
        _html_table("Summary Metrics", ("Metric", "Value"), _assessment_summary_rows(metrics)),
        _mapping_table(entry),
        _tool_table(tool_statuses),
    ]
    if metrics.per_residue is not None:
        sections.append(_per_residue_lddt_section(metrics.per_residue))
    if metrics.secondary_structure is not None:
        sections.append(_secondary_structure_section(metrics.secondary_structure))
    sections.append(_html_footer())
    return "\n".join(section for section in sections if section)


def _render_lddt_html(
    reference: str | Path,
    prediction: str | Path,
    result: LDDTResult,
    title: str,
) -> str:
    summary_rows = [
        ("lDDT", f"{result.lddt:.4f}"),
        ("Evaluated atoms", str(result.evaluated_atoms)),
        ("Evaluated pairs", str(result.evaluated_pairs)),
        ("Inclusion radius", f"{result.inclusion_radius:.1f} Å"),
    ]
    sections = [
        _html_header(title, _metadata("lddt")),
        _html_identity_block(str(reference), str(prediction)),
        _html_table("lDDT Summary", ("Metric", "Value"), summary_rows),
    ]
    if result.per_residue is not None:
        sections.append(_per_residue_lddt_section(result.per_residue))
    sections.append(_html_footer())
    return "\n".join(section for section in sections if section)


def _secondary_structure_section(
    comparison: SecondaryStructureComparisonResult,
) -> str:
    component = render_secondary_structure_comparison_component(
        comparison,
        title="Secondary Structure Comparison",
        component_id="secondary-structure-report-viewer",
    )
    return f"<section><h2>Secondary Structure</h2>{component}</section>"


def _per_residue_lddt_section(per_residue: tuple[ResidueAssessment, ...]) -> str:
    scored = [item for item in per_residue if item.lddt is not None]
    if not scored:
        return "<section><h2>Per-residue lDDT</h2><p>No per-residue lDDT values were available.</p></section>"

    mean_lddt = sum(item.lddt or 0.0 for item in scored) / len(scored)
    min_item = min(scored, key=lambda item: item.lddt if item.lddt is not None else 1.0)
    max_item = max(scored, key=lambda item: item.lddt if item.lddt is not None else -1.0)
    max_local_rmsd = max((item.local_rmsd or 0.0) for item in per_residue) or 1.0

    summary_cards = "".join(
        (
            _metric_card("Scored residues", str(len(scored))),
            _metric_card("Mean lDDT", f"{mean_lddt:.4f}"),
            _metric_card("Best residue", _residue_chip_label(max_item)),
            _metric_card("Worst residue", _residue_chip_label(min_item)),
        )
    )
    chips = "".join(_residue_chip(item) for item in per_residue)
    bars = "".join(_error_bar(item, max_local_rmsd) for item in per_residue)
    rows = "".join(_per_residue_row(item) for item in per_residue)

    return (
        "<section><h2>Per-residue lDDT</h2>"
        f"<div class='metric-grid'>{summary_cards}</div>"
        "<h3>Residue Heatmap</h3>"
        f"<div class='residue-strip'>{chips}</div>"
        "<h3>Local RMSD</h3>"
        f"<div class='bar-list'>{bars}</div>"
        "<h3>Residue Table</h3>"
        "<table><thead><tr>"
        "<th>Reference</th><th>Prediction</th><th>lDDT</th><th>Local RMSD</th>"
        "<th>Mean Abs Error</th><th>Max Abs Error</th><th>Matched atoms</th><th>Scored atoms</th>"
        "</tr></thead>"
        f"<tbody>{rows}</tbody></table>"
        "</section>"
    )


def _pair_list_block(title: str, pairs: tuple[SecondaryStructureComparisonPair, ...]) -> str:
    if not pairs:
        return f"<h3>{html.escape(title)}</h3><p>None</p>"
    items = "".join(f"<li>{html.escape(_pair_detail_line(pair))}</li>" for pair in pairs)
    return f"<h3>{html.escape(title)}</h3><ul>{items}</ul>"


def _pair_detail_line(pair: SecondaryStructureComparisonPair) -> str:
    primary = pair.native_pair or pair.prediction_pair
    if primary is None:
        return pair.status.upper()
    return (
        f"{pair.status.upper()} {primary.chain_1}:{primary.pos_1}({primary.nt_1}) - "
        f"{primary.chain_2}:{primary.pos_2}({primary.nt_2}) {primary.classification}"
    )


def _html_header(
    title: str,
    metadata: ReportMetadata,
    *,
    include_secondary_assets: bool = False,
) -> str:
    extra_assets = render_secondary_structure_assets() if include_secondary_assets else ""
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{html.escape(title)}</title>
  {extra_assets}
  <style>
    :root {{
      --bg: #f4f1ea;
      --panel: #fffdfa;
      --ink: #1f2933;
      --muted: #5b6778;
      --line: #d7d1c5;
      --accent: #2f6db3;
      --accent-2: #208b3a;
      --accent-3: #d97706;
    }}
    body {{
      margin: 0;
      font: 16px/1.5 Georgia, 'Iowan Old Style', serif;
      color: var(--ink);
      background:
        radial-gradient(circle at top left, rgba(47,109,179,0.08), transparent 28%),
        linear-gradient(180deg, #f7f4ee 0%, var(--bg) 100%);
    }}
    main {{
      max-width: 1240px;
      margin: 0 auto;
      padding: 36px 24px 64px;
    }}
    h1, h2, h3 {{
      font-family: 'Avenir Next', 'Trebuchet MS', sans-serif;
      letter-spacing: 0.02em;
      margin: 0 0 14px;
    }}
    section {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 18px;
      padding: 20px 22px;
      margin: 18px 0;
      box-shadow: 0 12px 32px rgba(31,41,51,0.06);
    }}
    .meta {{
      color: var(--muted);
      font-size: 14px;
      margin-top: 6px;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
    }}
    th, td {{
      text-align: left;
      padding: 8px 10px;
      border-bottom: 1px solid var(--line);
      vertical-align: top;
    }}
    th {{
      font-family: 'Avenir Next', 'Trebuchet MS', sans-serif;
      font-size: 13px;
      color: var(--muted);
      text-transform: uppercase;
      letter-spacing: 0.06em;
    }}
    ul {{
      margin: 8px 0 0;
      padding-left: 20px;
    }}
    code {{
      font-family: 'SFMono-Regular', Menlo, Consolas, monospace;
      font-size: 0.92em;
    }}
    .hero {{
      padding: 10px 0 8px;
    }}
    .metric-grid {{
      display: grid;
      gap: 12px;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      margin-bottom: 14px;
    }}
    .metric-card {{
      border: 1px solid var(--line);
      border-radius: 16px;
      background: linear-gradient(180deg, rgba(255,255,255,0.96), rgba(248,245,238,0.92));
      padding: 14px 16px;
    }}
    .metric-card span {{
      display: block;
      color: var(--muted);
      font-size: 0.85rem;
      margin-bottom: 4px;
    }}
    .metric-card strong {{
      font-family: 'Avenir Next', 'Trebuchet MS', sans-serif;
      font-size: 1.05rem;
    }}
    .residue-strip {{
      display: grid;
      gap: 10px;
      grid-template-columns: repeat(auto-fit, minmax(88px, 1fr));
      margin: 8px 0 18px;
    }}
    .residue-chip {{
      border: 1px solid color-mix(in srgb, var(--tile-color) 45%, var(--line));
      border-radius: 14px;
      background: linear-gradient(180deg, color-mix(in srgb, var(--tile-color) 26%, white), color-mix(in srgb, var(--tile-color) 64%, white));
      padding: 10px 10px 12px;
      min-height: 92px;
    }}
    .residue-chip-label {{
      display: block;
      font-family: 'Avenir Next', 'Trebuchet MS', sans-serif;
      font-size: 0.9rem;
      margin-bottom: 4px;
    }}
    .residue-chip small {{
      display: block;
      color: var(--muted);
    }}
    .residue-chip strong {{
      display: block;
      font-size: 1rem;
      margin-bottom: 4px;
    }}
    .bar-list {{
      display: grid;
      gap: 10px;
      margin: 8px 0 18px;
    }}
    .bar-row {{
      display: grid;
      grid-template-columns: minmax(110px, 140px) minmax(0, 1fr) 84px;
      gap: 12px;
      align-items: center;
    }}
    .bar-label {{
      font-family: 'Avenir Next', 'Trebuchet MS', sans-serif;
      font-size: 0.92rem;
    }}
    .bar-track {{
      height: 12px;
      border-radius: 999px;
      background: #ece6da;
      overflow: hidden;
      border: 1px solid var(--line);
    }}
    .bar-fill {{
      height: 100%;
      background: linear-gradient(90deg, #f59e0b, #d97706);
      border-radius: 999px;
    }}
    .bar-value {{
      text-align: right;
      color: var(--muted);
      font-variant-numeric: tabular-nums;
    }}
  </style>
</head>
<body>
<main>
  <div class="hero">
    <h1>{html.escape(title)}</h1>
    <div class="meta">schema={metadata.schema_version} library={metadata.library_version} generated={html.escape(metadata.generated_at)}</div>
  </div>"""


def _html_identity_block(reference: str, prediction: str) -> str:
    return (
        "<section><h2>Inputs</h2>"
        f"<p><strong>Reference:</strong> <code>{html.escape(reference)}</code><br />"
        f"<strong>Prediction:</strong> <code>{html.escape(prediction)}</code></p>"
        "</section>"
    )


def _html_table(title: str, headers: tuple[str, str], rows: list[tuple[str, str]]) -> str:
    body = "".join(
        f"<tr><td>{html.escape(label)}</td><td>{html.escape(value)}</td></tr>"
        for label, value in rows
    )
    return (
        f"<section><h2>{html.escape(title)}</h2><table>"
        f"<thead><tr><th>{html.escape(headers[0])}</th><th>{html.escape(headers[1])}</th></tr></thead>"
        f"<tbody>{body}</tbody></table></section>"
    )


def _mapping_table(mapping: BenchmarkEntry) -> str:
    rows = [
        ("Native index", mapping.native_index or "-"),
        ("Prediction index", mapping.prediction_index or "-"),
        ("Matched residues", "-" if mapping.matched_residues is None else str(mapping.matched_residues)),
        ("Used sidecar index", "yes" if mapping.used_sidecar_index else "no"),
        ("Used sequence hints", "yes" if mapping.used_sequence_hints else "no"),
        ("Used normalized inputs", "yes" if mapping.used_normalized_inputs else "no"),
        ("Used repaired inputs", "yes" if mapping.used_repaired_inputs else "no"),
    ]
    rows.extend(
        (
            f"Chain map {index}",
            f"{item.native_chain} -> {item.prediction_chain} ({item.matched_residues})",
        )
        for index, item in enumerate(mapping.chain_mappings, start=1)
    )
    return _html_table("Mapping", ("Field", "Value"), rows)


def _tool_table(statuses: tuple[ToolStatus, ...]) -> str:
    if not statuses:
        return ""
    body = "".join(
        "<tr>"
        f"<td>{html.escape(item.display_name)}</td>"
        f"<td>{'yes' if item.available else 'no'}</td>"
        f"<td>{html.escape(item.source or '-')}</td>"
        f"<td><code>{html.escape(item.binary_path or '-')}</code></td>"
        f"<td>{html.escape(item.notes or '-')}</td>"
        "</tr>"
        for item in statuses
    )
    return (
        "<section><h2>Tool Status</h2><table>"
        "<thead><tr><th>Tool</th><th>Available</th><th>Source</th><th>Binary</th><th>Notes</th></tr></thead>"
        f"<tbody>{body}</tbody></table></section>"
    )


def _warnings_block(warnings: tuple[str, ...]) -> str:
    if not warnings:
        return ""
    items = "".join(f"<li>{html.escape(item)}</li>" for item in warnings)
    return f"<section><h2>Warnings</h2><ul>{items}</ul></section>"


def _artifacts_block(artifacts: dict[str, str] | None) -> str:
    if not artifacts:
        return ""
    rows = [(key, value) for key, value in artifacts.items()]
    return _html_table("Artifacts", ("Kind", "Path"), rows)


def _html_footer() -> str:
    return "</main></body></html>"


def _metadata(report_type: str) -> ReportMetadata:
    return ReportMetadata(
        report_type=report_type,
        schema_version=SCHEMA_VERSION,
        library_version=_library_version(),
        generated_at=datetime.now(timezone.utc).isoformat(),
    )


def _library_version() -> str:
    try:
        return version("rna-kit")
    except PackageNotFoundError:
        return "0.0.0+local"


def _validate_document(
    document: AssessmentReportDocument | SecondaryStructureReportDocument | BenchmarkReportDocument,
) -> None:
    metadata = getattr(document, "metadata", None)
    if metadata is None or not metadata.schema_version:
        raise SchemaValidationError("Report document is missing metadata.schema_version.")
    if metadata.schema_version != SCHEMA_VERSION:
        raise SchemaValidationError(
            f"Report schema version '{metadata.schema_version}' does not match '{SCHEMA_VERSION}'."
        )


def _write_html(output_file: str | Path, content: str) -> Path:
    path = Path(output_file)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    except OSError as exc:
        raise ReportGenerationError(f"Failed to write HTML report to '{path}'.") from exc
    return path


def _metric_card(label: str, value: str) -> str:
    return (
        "<div class='metric-card'>"
        f"<span>{html.escape(label)}</span>"
        f"<strong>{html.escape(value)}</strong>"
        "</div>"
    )


def _residue_chip(item: ResidueAssessment) -> str:
    color = _lddt_color(item.lddt)
    lddt_value = "-" if item.lddt is None else f"{item.lddt:.3f}"
    mae_value = "-" if item.mean_absolute_error is None else f"{item.mean_absolute_error:.2f} Å"
    return (
        f"<div class='residue-chip' style='--tile-color: {html.escape(color)}'>"
        f"<span class='residue-chip-label'>{html.escape(_residue_chip_label(item))}</span>"
        f"<strong>{html.escape(lddt_value)}</strong>"
        f"<small>{html.escape(item.prediction_chain)}:{item.prediction_pos} {html.escape(item.prediction_nt)}</small>"
        f"<small>MAE {html.escape(mae_value)}</small>"
        "</div>"
    )


def _error_bar(item: ResidueAssessment, max_local_rmsd: float) -> str:
    label = _residue_chip_label(item)
    value = item.local_rmsd or 0.0
    ratio = max(0.0, min(1.0, value / max_local_rmsd))
    return (
        "<div class='bar-row'>"
        f"<div class='bar-label'>{html.escape(label)}</div>"
        f"<div class='bar-track'><div class='bar-fill' style='width: {ratio * 100:.2f}%'></div></div>"
        f"<div class='bar-value'>{value:.2f} Å</div>"
        "</div>"
    )


def _per_residue_row(item: ResidueAssessment) -> str:
    return (
        "<tr>"
        f"<td>{html.escape(f'{item.native_chain}:{item.native_pos} {item.native_nt}')}</td>"
        f"<td>{html.escape(f'{item.prediction_chain}:{item.prediction_pos} {item.prediction_nt}')}</td>"
        f"<td>{html.escape('-' if item.lddt is None else f'{item.lddt:.4f}')}</td>"
        f"<td>{html.escape('-' if item.local_rmsd is None else f'{item.local_rmsd:.4f}')}</td>"
        f"<td>{html.escape('-' if item.mean_absolute_error is None else f'{item.mean_absolute_error:.4f}')}</td>"
        f"<td>{html.escape('-' if item.max_absolute_error is None else f'{item.max_absolute_error:.4f}')}</td>"
        f"<td>{item.matched_atoms}</td>"
        f"<td>{item.scored_atoms}</td>"
        "</tr>"
    )


def _residue_chip_label(item: ResidueAssessment) -> str:
    return f"{item.native_chain}:{item.native_pos} {item.native_nt}"


def _lddt_color(score: float | None) -> str:
    if score is None:
        return "#d7dbe2"
    clamped = max(0.0, min(1.0, score))
    hue = 8 + clamped * 122
    lightness = 92 - clamped * 28
    return f"hsl({hue:.0f}, 68%, {lightness:.0f}%)"


def _benchmark_row(entry: BenchmarkEntry, detail_link: str | None = None) -> str:
    metrics = entry.metrics
    molprobity = None if metrics is None else metrics.molprobity
    label_text = entry.label or "-"
    if detail_link is not None and metrics is not None:
        label_cell = f"<a href=\"{html.escape(detail_link)}\">{html.escape(label_text)}</a>"
        prediction_cell = f"<a href=\"{html.escape(detail_link)}\"><code>{html.escape(entry.prediction)}</code></a>"
    else:
        label_cell = html.escape(label_text)
        prediction_cell = f"<code>{html.escape(entry.prediction)}</code>"
    return (
        "<tr>"
        f"<td>{label_cell}</td>"
        f"<td>{prediction_cell}</td>"
        f"<td>{html.escape(entry.status)}</td>"
        f"<td>{'-' if entry.matched_residues is None else entry.matched_residues}</td>"
        f"<td>{_format_float(None if metrics is None else metrics.rmsd)}</td>"
        f"<td>{_format_float(None if metrics is None else metrics.inf_all)}</td>"
        f"<td>{_format_float(None if metrics is None else metrics.lddt)}</td>"
        f"<td>{_format_float(None if metrics is None else metrics.secondary_structure_f1)}</td>"
        f"<td>{_format_float(None if molprobity is None else molprobity.clashscore)}</td>"
        f"<td>{_format_float(None if molprobity is None else molprobity.molprobity_score)}</td>"
        f"<td>{html.escape(entry.error or '-')}</td>"
        "</tr>"
    )


def _format_float(value: float | None) -> str:
    return "-" if value is None else f"{value:.4f}"


def _mean(values: list[float]) -> float | None:
    return None if not values else sum(values) / len(values)


def _assessment_summary_rows(metrics: AssessmentResult) -> list[tuple[str, str]]:
    rows = [
        ("RMSD", f"{metrics.rmsd:.4f}"),
        ("P-value", f"{metrics.pvalue:.3e}"),
        ("INF_ALL", f"{metrics.inf_all:.4f}"),
        ("INF_WC", f"{metrics.inf_wc:.4f}"),
        ("INF_NWC", f"{metrics.inf_nwc:.4f}"),
        ("INF_STACK", f"{metrics.inf_stack:.4f}"),
        ("lDDT", f"{metrics.lddt:.4f}"),
    ]
    if metrics.secondary_structure_f1 is not None:
        rows.extend(
            [
                ("SS Precision", f"{metrics.secondary_structure_precision:.4f}"),
                ("SS Recall", f"{metrics.secondary_structure_recall:.4f}"),
                ("SS F1", f"{metrics.secondary_structure_f1:.4f}"),
                ("SS Jaccard", f"{metrics.secondary_structure_jaccard:.4f}"),
            ]
        )
    if metrics.molprobity is not None:
        if metrics.molprobity.clashscore is not None:
            rows.append(("MolProbity clashscore", f"{metrics.molprobity.clashscore:.4f}"))
        if metrics.molprobity.molprobity_score is not None:
            rows.append(("MolProbity score", f"{metrics.molprobity.molprobity_score:.4f}"))
    return rows


def _slugify(value: str) -> str:
    cleaned = "".join(char.lower() if char.isalnum() else "-" for char in value.strip())
    cleaned = "-".join(part for part in cleaned.split("-") if part)
    return cleaned or "entry"
