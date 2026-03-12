from __future__ import annotations

import html
import json
from functools import lru_cache
from pathlib import Path

from .exceptions import ReportGenerationError
from .secondary_structure import (
    SecondaryStructureBasePair,
    SecondaryStructureComparisonPair,
    SecondaryStructureComparisonResult,
    SecondaryStructureResult,
)

FORNAC_VERSION = "1.1.8"
FORNAC_ASSET_DIR = Path(__file__).resolve().parents[2] / "third_party" / "web" / "fornac"
FORNAC_JS_FILE = FORNAC_ASSET_DIR / "fornac.js"
FORNAC_CSS_FILE = FORNAC_ASSET_DIR / "fornac.css"

PAIRED_COLOR = "#2f6db3"
TRUE_POSITIVE_COLOR = "#208b3a"
FALSE_POSITIVE_COLOR = "#d97706"
FALSE_NEGATIVE_COLOR = "#c2410c"
MIXED_ASSIGNMENT_COLOR = "#7c5e10"


def render_secondary_structure_assets() -> str:
    css = _asset_text(FORNAC_CSS_FILE)
    js = _asset_text(FORNAC_JS_FILE)
    return (
        f"<!-- FornaC {FORNAC_VERSION} -->\n"
        f"<style>\n{css}\n</style>\n"
        f"<script>\n{js}\n</script>"
    )


def render_secondary_structure_html(
    result: SecondaryStructureResult,
    title: str = "RNA Kit Secondary Structure",
) -> str:
    component = render_secondary_structure_component(
        result,
        title=title,
        component_id="secondary-structure-view",
    )
    return _document(title, component)


def render_secondary_structure_comparison_html(
    comparison: SecondaryStructureComparisonResult,
    title: str = "RNA Kit Secondary Structure Comparison",
) -> str:
    component = render_secondary_structure_comparison_component(
        comparison,
        title=title,
        component_id="secondary-structure-compare",
    )
    return _document(title, component)


def write_secondary_structure_html(
    result: SecondaryStructureResult,
    output_file: str | Path,
    title: str = "RNA Kit Secondary Structure",
) -> Path:
    path = Path(output_file)
    return _write_html(path, render_secondary_structure_html(result, title=title))


def write_secondary_structure_comparison_html(
    comparison: SecondaryStructureComparisonResult,
    output_file: str | Path,
    title: str = "RNA Kit Secondary Structure Comparison",
) -> Path:
    path = Path(output_file)
    return _write_html(path, render_secondary_structure_comparison_html(comparison, title=title))


def render_secondary_structure_component(
    result: SecondaryStructureResult,
    title: str,
    component_id: str,
) -> str:
    payload = {
        "title": title,
        "panels": [
            {
                "label": "Structure",
                "sequence": result.sequence,
                "dot_bracket": result.dot_bracket,
                "selected_residues": result.selected_residues,
                "base_pair_count": result.base_pair_count,
                "custom_colors": _color_text(_paired_residue_colors(result.base_pairs, PAIRED_COLOR)),
            }
        ],
    }
    return _component_markup(
        payload,
        component_id=component_id,
        summary_html=_single_summary(result),
        pair_lists="",
        legend_html=_legend(
            (
                ("Paired residues", PAIRED_COLOR),
            )
        ),
    )


def render_secondary_structure_comparison_component(
    comparison: SecondaryStructureComparisonResult,
    title: str,
    component_id: str,
) -> str:
    payload = {
        "title": title,
        "panels": [
            {
                "label": "Native",
                "sequence": comparison.native.sequence,
                "dot_bracket": comparison.native.dot_bracket,
                "selected_residues": comparison.native.selected_residues,
                "base_pair_count": comparison.native.base_pair_count,
                "custom_colors": _comparison_color_text(
                    comparison.true_positive_pairs,
                    comparison.false_negative_pairs,
                ),
            },
            {
                "label": "Prediction",
                "sequence": comparison.prediction.sequence,
                "dot_bracket": comparison.prediction.dot_bracket,
                "selected_residues": comparison.prediction.selected_residues,
                "base_pair_count": comparison.prediction.base_pair_count,
                "custom_colors": _comparison_color_text(
                    comparison.true_positive_pairs,
                    comparison.false_positive_pairs,
                ),
            },
        ],
    }
    return _component_markup(
        payload,
        component_id=component_id,
        summary_html=_comparison_summary(comparison),
        pair_lists=_comparison_pair_lists(comparison),
        legend_html=_legend(
            (
                ("TP residues", TRUE_POSITIVE_COLOR),
                ("FN residues", FALSE_NEGATIVE_COLOR),
                ("FP residues", FALSE_POSITIVE_COLOR),
                ("Mixed assignments", MIXED_ASSIGNMENT_COLOR),
            )
        ),
    )


def _component_markup(
    payload: dict,
    component_id: str,
    summary_html: str,
    pair_lists: str,
    legend_html: str,
) -> str:
    controls_id = f"{component_id}-controls"
    return f"""
<section class="fornac-web">
  <div class="fornac-web-head">
    <div>
      <h2>{html.escape(payload["title"])}</h2>
      <p class="fornac-web-note">MC-Annotate cisWW base pairs are rendered with FornaC.</p>
    </div>
    {summary_html}
  </div>
  <div class="fornac-web-toolbar">
    <div id="{html.escape(controls_id)}" class="fornac-web-controls">
      <label><input type="checkbox" data-control="numbering" checked /> Numbering</label>
      <label><input type="checkbox" data-control="outline" checked /> Node outline</label>
      <label><input type="checkbox" data-control="labels" checked /> Base labels</label>
      <label><input type="checkbox" data-control="links" checked /> Links</label>
      <label><input type="checkbox" data-control="pseudoknots" checked /> Pseudoknot links</label>
    </div>
    {legend_html}
  </div>
  <div id="{html.escape(component_id)}" class="fornac-web-viewer"></div>
  {pair_lists}
</section>
<script>
(() => {{
  const host = document.getElementById({json.dumps(component_id)});
  const controlsHost = document.getElementById({json.dumps(controls_id)});
  if (!host || !controlsHost) return;

  const payload = {json.dumps(payload)};
  const FornaCtor = (window.fornac && window.fornac.FornaContainer) || window.FornaContainer;
  if (!FornaCtor) {{
    host.innerHTML = '<div class="fornac-error">FornaC failed to initialize from bundled assets.</div>';
    return;
  }}

  const viewers = [];
  const toggleMethods = {{
    numbering: "displayNumbering",
    outline: "displayNodeOutline",
    labels: "displayNodeLabel",
    links: "displayLinks",
    pseudoknots: "displayPseudoknotLinks",
  }};

  function setToggle(container, controlName, enabled) {{
    const method = toggleMethods[controlName];
    if (!method || typeof container[method] !== "function") {{
      return;
    }}
    try {{
      container[method](enabled);
    }} catch (error) {{
      console.warn("FornaC toggle failed:", controlName, error);
    }}
  }}

  function refreshView(container) {{
    if (typeof container.setSize === "function") {{
      container.setSize();
    }}
  }}

  function applyControls() {{
    const states = Object.fromEntries(
      Array.from(controlsHost.querySelectorAll("input[data-control]")).map((input) => [input.dataset.control, input.checked])
    );
    viewers.forEach((item) => {{
      Object.entries(states).forEach(([controlName, enabled]) => setToggle(item.container, controlName, enabled));
      refreshView(item.container);
    }});
  }}

  function createPanel(panel, index) {{
    const card = document.createElement("section");
    card.className = "fornac-panel";

    const header = document.createElement("div");
    header.className = "fornac-panel-head";

    const heading = document.createElement("div");
    heading.innerHTML = `
      <h3>${{panel.label}}</h3>
      <p class="fornac-panel-meta">${{panel.selected_residues}} residues · ${{panel.base_pair_count}} base pairs</p>
    `;
    header.appendChild(heading);
    card.appendChild(header);

    const viewer = document.createElement("div");
    viewer.className = "fornac-viewer";
    viewer.id = `${{payload.title.replace(/[^a-zA-Z0-9]+/g, "-").toLowerCase() || "fornac"}}-${{index}}-${{Math.random().toString(36).slice(2, 8)}}`;
    card.appendChild(viewer);

    const sequenceBlock = document.createElement("pre");
    sequenceBlock.className = "fornac-sequence";
    sequenceBlock.textContent = `${{panel.sequence}}\\n${{panel.dot_bracket}}`;
    card.appendChild(sequenceBlock);

    host.appendChild(card);

    const container = new FornaCtor(`#${{viewer.id}}`, {{
      applyForce: false,
      allowPanningAndZooming: true,
      initialSize: [560, 360],
    }});
    container.addRNA(panel.dot_bracket, {{
      structure: panel.dot_bracket,
      sequence: panel.sequence,
    }});
    if (panel.custom_colors) {{
      container.addCustomColorsText(panel.custom_colors);
    }}
    refreshView(container);
    viewers.push({{ container }});

    if (typeof ResizeObserver !== "undefined") {{
      const observer = new ResizeObserver(() => refreshView(container));
      observer.observe(viewer);
    }} else {{
      window.addEventListener("resize", () => refreshView(container));
    }}
  }}

  payload.panels.forEach(createPanel);
  controlsHost.addEventListener("change", applyControls);
  applyControls();
}})();
</script>"""


def _single_summary(result: SecondaryStructureResult) -> str:
    return (
        '<div class="fornac-web-metrics">'
        f'<div><span>Residues</span><strong>{result.selected_residues}</strong></div>'
        f'<div><span>Base pairs</span><strong>{result.base_pair_count}</strong></div>'
        "</div>"
    )


def _comparison_summary(comparison: SecondaryStructureComparisonResult) -> str:
    return (
        '<div class="fornac-web-metrics">'
        f'<div><span>Precision</span><strong>{comparison.precision:.4f}</strong></div>'
        f'<div><span>Recall</span><strong>{comparison.recall:.4f}</strong></div>'
        f'<div><span>F1</span><strong>{comparison.f1:.4f}</strong></div>'
        f'<div><span>Jaccard</span><strong>{comparison.jaccard:.4f}</strong></div>'
        "</div>"
    )


def _comparison_pair_lists(comparison: SecondaryStructureComparisonResult) -> str:
    return (
        '<div class="fornac-pair-grid">'
        f'{_pair_list("True Positives", comparison.true_positive_pairs)}'
        f'{_pair_list("False Positives", comparison.false_positive_pairs)}'
        f'{_pair_list("False Negatives", comparison.false_negative_pairs)}'
        "</div>"
    )


def _pair_list(title: str, pairs: tuple[SecondaryStructureComparisonPair, ...]) -> str:
    if not pairs:
        body = "<li>None</li>"
    else:
        body = "".join(f"<li>{html.escape(_pair_label(pair))}</li>" for pair in pairs)
    return f'<div class="fornac-pair-list"><h3>{html.escape(title)}</h3><ul>{body}</ul></div>'


def _pair_label(pair: SecondaryStructureComparisonPair) -> str:
    primary = pair.native_pair or pair.prediction_pair
    if primary is None:
        return pair.status.upper()
    return (
        f"{primary.chain_1}:{primary.pos_1}({primary.nt_1}) - "
        f"{primary.chain_2}:{primary.pos_2}({primary.nt_2}) {primary.classification} [{pair.status.upper()}]"
    )


def _legend(items: tuple[tuple[str, str], ...]) -> str:
    chips = "".join(
        f'<span class="fornac-legend-chip"><span class="fornac-legend-swatch" style="background:{html.escape(color)}"></span>{html.escape(label)}</span>'
        for label, color in items
    )
    return f'<div class="fornac-legend">{chips}</div>'


def _paired_residue_colors(
    base_pairs: tuple[SecondaryStructureBasePair, ...],
    color: str,
) -> dict[int, str]:
    colors: dict[int, str] = {}
    for pair in base_pairs:
        colors[pair.rank_1] = color
        colors[pair.rank_2] = color
    return colors


def _comparison_color_text(
    primary_pairs: tuple[SecondaryStructureComparisonPair, ...],
    secondary_pairs: tuple[SecondaryStructureComparisonPair, ...],
) -> str | None:
    colors: dict[int, str] = {}
    _add_pair_colors(colors, primary_pairs, TRUE_POSITIVE_COLOR)
    fallback_color = FALSE_NEGATIVE_COLOR
    if secondary_pairs and secondary_pairs[0].status == "fp":
        fallback_color = FALSE_POSITIVE_COLOR
    _add_pair_colors(colors, secondary_pairs, fallback_color)
    return _color_text(colors)


def _add_pair_colors(
    colors: dict[int, str],
    pairs: tuple[SecondaryStructureComparisonPair, ...],
    color: str,
) -> None:
    for pair in pairs:
        _assign_color(colors, pair.rank_1, color)
        _assign_color(colors, pair.rank_2, color)


def _assign_color(colors: dict[int, str], rank: int, color: str) -> None:
    previous = colors.get(rank)
    if previous is None or previous == color:
        colors[rank] = color
        return
    colors[rank] = MIXED_ASSIGNMENT_COLOR


def _color_text(colors: dict[int, str]) -> str | None:
    if not colors:
        return None
    return "\n".join(f"{rank + 1}:{color}" for rank, color in sorted(colors.items()))


def _document(title: str, component: str) -> str:
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{html.escape(title)}</title>
  {render_secondary_structure_assets()}
  <style>
    * {{ box-sizing: border-box; }}
    :root {{
      --bg: #f4f1ea;
      --panel: rgba(255, 253, 250, 0.96);
      --ink: #1f2933;
      --muted: #52606d;
      --line: #d7d1c5;
    }}
    body {{
      margin: 0;
      font: 16px/1.5 Georgia, 'Iowan Old Style', serif;
      color: var(--ink);
      background:
        radial-gradient(circle at top left, rgba(47,109,179,0.1), transparent 28%),
        linear-gradient(180deg, #f7f4ee 0%, var(--bg) 100%);
    }}
    main {{
      max-width: 1320px;
      margin: 0 auto;
      padding: 32px 24px 56px;
    }}
    .fornac-web {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 20px;
      padding: 22px;
      box-shadow: 0 18px 44px rgba(31, 41, 51, 0.08);
    }}
    .fornac-web-head {{
      display: flex;
      justify-content: space-between;
      gap: 18px;
      align-items: start;
      flex-wrap: wrap;
      margin-bottom: 18px;
    }}
    .fornac-web h2, .fornac-web h3 {{
      margin: 0 0 10px;
      font-weight: 600;
    }}
    .fornac-web-note {{
      margin: 0;
      color: var(--muted);
    }}
    .fornac-web-toolbar {{
      display: flex;
      justify-content: space-between;
      align-items: start;
      gap: 16px;
      flex-wrap: wrap;
      margin-bottom: 18px;
    }}
    .fornac-web-controls {{
      display: flex;
      gap: 10px 14px;
      flex-wrap: wrap;
      align-items: center;
    }}
    .fornac-web-controls label {{
      display: inline-flex;
      align-items: center;
      gap: 8px;
      padding: 8px 12px;
      border: 1px solid var(--line);
      border-radius: 999px;
      background: rgba(255,255,255,0.7);
      font-size: 0.95rem;
    }}
    .fornac-legend {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px 12px;
      align-items: center;
    }}
    .fornac-legend-chip {{
      display: inline-flex;
      align-items: center;
      gap: 8px;
      color: var(--muted);
      font-size: 0.92rem;
    }}
    .fornac-legend-swatch {{
      width: 12px;
      height: 12px;
      border-radius: 999px;
      display: inline-block;
      box-shadow: inset 0 0 0 1px rgba(31,41,51,0.18);
    }}
    .fornac-web-viewer {{
      display: grid;
      gap: 18px;
      grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
    }}
    .fornac-panel {{
      border: 1px solid var(--line);
      border-radius: 16px;
      padding: 16px;
      background: rgba(255,255,255,0.62);
      overflow: hidden;
    }}
    .fornac-panel-meta {{
      margin: 0;
      color: var(--muted);
    }}
    .fornac-viewer {{
      min-height: 360px;
      border: 1px solid var(--line);
      border-radius: 14px;
      background: linear-gradient(180deg, #f7f8fb 0%, #eef2f7 100%);
      overflow: hidden;
      position: relative;
      isolation: isolate;
    }}
    .fornac-viewer svg {{
      width: 100%;
      height: 100%;
    }}
    .fornac-sequence {{
      margin: 14px 0 0;
      padding: 12px;
      overflow-x: auto;
      border-radius: 14px;
      background: #f8f5ee;
      border: 1px solid var(--line);
      font: 14px/1.45 'SFMono-Regular', Menlo, Monaco, Consolas, monospace;
    }}
    .fornac-web-metrics {{
      display: flex;
      gap: 12px;
      flex-wrap: wrap;
    }}
    .fornac-web-metrics div {{
      min-width: 112px;
      padding: 10px 12px;
      border: 1px solid var(--line);
      border-radius: 14px;
      background: rgba(255,255,255,0.7);
    }}
    .fornac-web-metrics span {{
      display: block;
      color: var(--muted);
      font-size: 0.86rem;
    }}
    .fornac-web-metrics strong {{
      font-size: 1.05rem;
    }}
    .fornac-pair-grid {{
      display: grid;
      gap: 14px;
      grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
      margin-top: 18px;
    }}
    .fornac-pair-list {{
      border: 1px solid var(--line);
      border-radius: 16px;
      background: rgba(255,255,255,0.6);
      padding: 14px 16px;
    }}
    .fornac-pair-list ul {{
      margin: 0;
      padding-left: 20px;
    }}
    .fornac-error {{
      padding: 18px;
      border: 1px solid rgba(194,65,12,0.28);
      border-radius: 16px;
      background: rgba(194,65,12,0.08);
      color: #8a2c0d;
    }}
    @media (max-width: 960px) {{
      main {{ padding: 20px 14px 40px; }}
    }}
  </style>
</head>
<body>
  <main>
    {component}
  </main>
</body>
</html>
"""


@lru_cache(maxsize=2)
def _asset_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ReportGenerationError(f"Failed to read bundled FornaC asset '{path}'.") from exc


def _write_html(path: Path, content: str) -> Path:
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    except OSError as exc:
        raise ReportGenerationError(f"Failed to write secondary-structure HTML to '{path}'.") from exc
    return path
