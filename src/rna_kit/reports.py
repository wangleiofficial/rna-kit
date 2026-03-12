from __future__ import annotations

import html
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from .benchmark import BenchmarkEntry, describe_prepared_pair
from .exceptions import ReportGenerationError, SchemaValidationError
from .metrics import AssessmentResult, PreparedStructurePair
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


def write_report_json(document: AssessmentReportDocument | SecondaryStructureReportDocument, output_file: str | Path) -> Path:
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


def _render_assessment_html(
    document: AssessmentReportDocument,
) -> str:
    metrics = document.metrics
    summary_rows = [
        ("RMSD", f"{metrics.rmsd:.4f}"),
        ("P-value", f"{metrics.pvalue:.3e}"),
        ("INF_ALL", f"{metrics.inf_all:.4f}"),
        ("INF_WC", f"{metrics.inf_wc:.4f}"),
        ("INF_NWC", f"{metrics.inf_nwc:.4f}"),
        ("INF_STACK", f"{metrics.inf_stack:.4f}"),
        ("lDDT", f"{metrics.lddt:.4f}"),
    ]
    if metrics.secondary_structure_f1 is not None:
        summary_rows.extend(
            [
                ("SS Precision", f"{metrics.secondary_structure_precision:.4f}"),
                ("SS Recall", f"{metrics.secondary_structure_recall:.4f}"),
                ("SS F1", f"{metrics.secondary_structure_f1:.4f}"),
                ("SS Jaccard", f"{metrics.secondary_structure_jaccard:.4f}"),
            ]
        )

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


def _secondary_structure_section(
    comparison: SecondaryStructureComparisonResult,
) -> str:
    component = render_secondary_structure_comparison_component(
        comparison,
        title="Secondary Structure Comparison",
        component_id="secondary-structure-report-viewer",
    )
    return f"<section><h2>Secondary Structure</h2>{component}</section>"


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


def _validate_document(document: AssessmentReportDocument | SecondaryStructureReportDocument) -> None:
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
