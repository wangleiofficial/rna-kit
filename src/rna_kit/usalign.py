from __future__ import annotations

import html
import json
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

from .exceptions import ReportGenerationError, ToolExecutionError, ToolResolutionError
from .tools import default_tool_registry


@dataclass(frozen=True)
class USAlignResult:
    reference: str
    prediction: str
    reference_name: str | None
    prediction_name: str | None
    reference_length: int
    prediction_length: int
    aligned_length: int
    rmsd: float
    sequence_identity: float
    tm_score_reference: float
    tm_score_prediction: float
    alignment_prediction: str | None
    alignment_markup: str | None
    alignment_reference: str | None
    binary_path: str
    output_directory: str | None = None
    reference_structure_output: str | None = None
    prediction_structure_output: str | None = None
    superposition_output: str | None = None
    superposed_prediction_output: str | None = None


class USAlignRunner:
    _NAME_PATTERN = re.compile(
        r"^Name of Structure_(?P<index>[12]):\s*(?P<name>.+?)(?:\s+\(to be superimposed onto Structure_2\))?$"
    )
    _LENGTH_PATTERN = re.compile(r"^Length of Structure_(?P<index>[12]):\s*(?P<length>\d+)\s+residues$")
    _SUMMARY_PATTERN = re.compile(
        r"^Aligned length=\s*(?P<aligned>\d+), RMSD=\s*(?P<rmsd>[0-9.]+), "
        r"Seq_ID=n_identical/n_aligned=\s*(?P<seqid>[0-9.]+)$"
    )
    _TM_PATTERN = re.compile(
        r"^TM-score=\s*(?P<score>[0-9.]+)\s+\(normalized by length of Structure_(?P<index>[12]):"
    )

    def __init__(self, binary_path: str | Path | None = None) -> None:
        self.binary_path = Path(binary_path) if binary_path else None

    def resolve_binary(self) -> Path:
        registry = default_tool_registry()
        status = registry.status("us_align", override=self.binary_path)
        if status.binary_path is not None:
            return Path(status.binary_path)
        raise ToolResolutionError(
            "US-align is not available. Set 'RNA_KIT_US_ALIGN', pass '--us-align', "
            "or use the bundled executable in 'third_party/bin'."
        )

    def align(
        self,
        reference_file: str | Path,
        prediction_file: str | Path,
        output_dir: str | Path | None = None,
    ) -> USAlignResult:
        if output_dir is None:
            with TemporaryDirectory(prefix="rna-kit-usalign-") as temp_dir:
                return self._align(
                    reference_file,
                    prediction_file,
                    work_dir=Path(temp_dir),
                    persist_artifacts=False,
                )

        work_dir = Path(output_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
        return self._align(
            reference_file,
            prediction_file,
            work_dir=work_dir,
            persist_artifacts=True,
        )

    def _align(
        self,
        reference_file: str | Path,
        prediction_file: str | Path,
        *,
        work_dir: Path,
        persist_artifacts: bool,
    ) -> USAlignResult:
        binary = self.resolve_binary()
        reference_path = Path(reference_file)
        prediction_path = Path(prediction_file)

        reference_input = _prepare_input_structure(reference_path, work_dir / "reference_input.pdb", persist_artifacts)
        prediction_input = _prepare_input_structure(prediction_path, work_dir / "prediction_input.pdb", persist_artifacts)
        output_prefix = work_dir / "sup"

        try:
            completed = subprocess.run(
                [
                    str(binary),
                    str(prediction_input),
                    str(reference_input),
                    "-mol",
                    "RNA",
                    "-o",
                    str(output_prefix),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
        except OSError as exc:
            raise ToolExecutionError(f"Failed to execute US-align at '{binary}'.") from exc
        except subprocess.CalledProcessError as exc:
            message = exc.stderr.strip() or exc.stdout.strip() or "unknown error"
            raise ToolExecutionError(
                f"US-align failed for prediction '{prediction_path}' against reference '{reference_path}': {message}"
            ) from exc

        superposition_output = output_prefix.with_suffix(".pdb")
        superposed_prediction_output = work_dir / "superposed_prediction.pdb"
        if superposition_output.exists():
            _extract_first_model(superposition_output, superposed_prediction_output)

        return _parse_output(
            completed.stdout,
            reference=reference_path,
            prediction=prediction_path,
            binary=binary,
            reference_input=reference_input if persist_artifacts else None,
            prediction_input=prediction_input if persist_artifacts else None,
            output_dir=work_dir if persist_artifacts else None,
            superposition_output=superposition_output if persist_artifacts and superposition_output.exists() else None,
            superposed_prediction_output=(
                superposed_prediction_output if persist_artifacts and superposed_prediction_output.exists() else None
            ),
        )


def calculate_us_align(
    reference_file: str | Path,
    prediction_file: str | Path,
    runner: USAlignRunner | None = None,
    output_dir: str | Path | None = None,
) -> USAlignResult:
    runner = runner or USAlignRunner()
    return runner.align(
        reference_file,
        prediction_file,
        output_dir=output_dir,
    )


def render_us_align_html(
    reference_file: str | Path,
    superposed_prediction_file: str | Path,
    result: USAlignResult,
    title: str = "RNA Kit Structure Alignment",
) -> str:
    reference_path = Path(reference_file)
    superposed_prediction_path = Path(superposed_prediction_file)
    try:
        reference_text = reference_path.read_text(encoding="utf-8")
        superposed_prediction_text = superposed_prediction_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ReportGenerationError("Failed to read aligned structures for HTML visualization.") from exc

    summary_rows = [
        ("Aligned length", str(result.aligned_length)),
        ("US-align RMSD", f"{result.rmsd:.4f}"),
        ("Seq_ID", f"{result.sequence_identity:.4f}"),
        ("TM-score (reference)", f"{result.tm_score_reference:.5f}"),
        ("TM-score (prediction)", f"{result.tm_score_prediction:.5f}"),
        ("Reference length", str(result.reference_length)),
        ("Prediction length", str(result.prediction_length)),
    ]
    summary_html = "".join(
        f"<tr><th>{html.escape(label)}</th><td>{html.escape(value)}</td></tr>" for label, value in summary_rows
    )
    alignment_block = ""
    if result.alignment_prediction and result.alignment_markup and result.alignment_reference:
        alignment_block = (
            "<section class='alignment'><h2>Sequence-Independent Alignment</h2><pre>"
            f"{html.escape(result.alignment_prediction)}\n"
            f"{html.escape(result.alignment_markup)}\n"
            f"{html.escape(result.alignment_reference)}"
            "</pre></section>"
        )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{html.escape(title)}</title>
  <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
  <style>
    * {{
      box-sizing: border-box;
    }}
    :root {{
      --bg: #f4f1ea;
      --panel: rgba(255, 253, 250, 0.96);
      --ink: #1f2933;
      --muted: #52606d;
      --line: #d7d1c5;
      --reference: #2f6db3;
      --prediction: #d97706;
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
      max-width: 1280px;
      margin: 0 auto;
      padding: 32px 24px 56px;
    }}
    h1, h2 {{
      margin: 0 0 12px;
      font-weight: 600;
    }}
    p {{
      margin: 0 0 12px;
      color: var(--muted);
    }}
    .panel {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 18px;
      padding: 18px 20px;
      box-shadow: 0 18px 44px rgba(31, 41, 51, 0.08);
      margin-bottom: 18px;
    }}
    .viewer-panel {{
      overflow: hidden;
    }}
    #viewer {{
      width: 100%;
      max-width: 100%;
      height: 720px;
      position: relative;
      border-radius: 16px;
      border: 1px solid var(--line);
      background: linear-gradient(180deg, #f7f8fb 0%, #eef2f7 100%);
      overflow: hidden;
      isolation: isolate;
    }}
    #viewer > div,
    #viewer canvas {{
      width: 100% !important;
      height: 100% !important;
      max-width: 100%;
      max-height: 100%;
      display: block;
    }}
    .grid {{
      display: grid;
      gap: 18px;
      grid-template-columns: minmax(0, 1.8fr) minmax(280px, 0.9fr);
      align-items: start;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.95rem;
    }}
    th, td {{
      text-align: left;
      padding: 8px 0;
      border-bottom: 1px solid var(--line);
      vertical-align: top;
    }}
    .legend {{
      display: flex;
      gap: 16px;
      flex-wrap: wrap;
      margin-top: 14px;
      font-size: 0.95rem;
    }}
    .swatch {{
      width: 12px;
      height: 12px;
      border-radius: 999px;
      display: inline-block;
      margin-right: 8px;
      vertical-align: middle;
    }}
    .alignment pre {{
      margin: 0;
      padding: 14px;
      overflow-x: auto;
      border-radius: 14px;
      background: #f8f5ee;
      border: 1px solid var(--line);
      font: 14px/1.45 'SFMono-Regular', Menlo, Monaco, Consolas, monospace;
    }}
    @media (max-width: 960px) {{
      main {{
        padding: 20px 14px 40px;
      }}
      .grid {{
        grid-template-columns: 1fr;
      }}
      #viewer {{
        height: 520px;
      }}
    }}
  </style>
</head>
<body>
  <main>
    <section class="panel">
      <h1>{html.escape(title)}</h1>
      <p>Reference: {html.escape(result.reference_name or result.reference)}</p>
      <p>Prediction: {html.escape(result.prediction_name or result.prediction)}</p>
    </section>
    <section class="grid">
      <div class="panel viewer-panel">
        <div id="viewer"></div>
        <div class="legend">
          <span><span class="swatch" style="background: var(--reference);"></span>Reference structure</span>
          <span><span class="swatch" style="background: var(--prediction);"></span>Prediction after US-align superposition</span>
        </div>
      </div>
      <aside class="panel">
        <h2>Alignment Summary</h2>
        <table>
          <tbody>
            {summary_html}
          </tbody>
        </table>
      </aside>
    </section>
    {alignment_block}
  </main>
  <script>
    const viewer = $3Dmol.createViewer("viewer", {{ backgroundColor: "white" }});
    const referencePdb = {json.dumps(reference_text)};
    const superposedPredictionPdb = {json.dumps(superposed_prediction_text)};
    viewer.addModel(referencePdb, "pdb");
    viewer.setStyle({{ model: 0 }}, {{
      cartoon: {{ color: "{_css_color('reference')}", opacity: 0.9 }},
      stick: {{ color: "{_css_color('reference')}", radius: 0.14 }}
    }});
    viewer.addModel(superposedPredictionPdb, "pdb");
    viewer.setStyle({{ model: 1 }}, {{
      cartoon: {{ color: "{_css_color('prediction')}", opacity: 0.45 }},
      stick: {{ color: "{_css_color('prediction')}", radius: 0.14 }}
    }});
    viewer.zoomTo();
    viewer.render();
  </script>
</body>
</html>
"""


def write_us_align_html(
    reference_file: str | Path,
    superposed_prediction_file: str | Path,
    output_file: str | Path,
    result: USAlignResult,
    title: str = "RNA Kit Structure Alignment",
) -> Path:
    path = Path(output_file)
    html_output = render_us_align_html(
        reference_file,
        superposed_prediction_file,
        result=result,
        title=title,
    )
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(html_output, encoding="utf-8")
    except OSError as exc:
        raise ReportGenerationError(f"Failed to write US-align HTML visualization to '{path}'.") from exc
    return path


def _prepare_input_structure(
    source: Path,
    target: Path,
    persist_artifacts: bool,
) -> Path:
    if persist_artifacts:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        return target
    return source


def _parse_output(
    stdout: str,
    *,
    reference: Path,
    prediction: Path,
    binary: Path,
    reference_input: Path | None,
    prediction_input: Path | None,
    output_dir: Path | None,
    superposition_output: Path | None,
    superposed_prediction_output: Path | None,
) -> USAlignResult:
    names: dict[int, str] = {}
    lengths: dict[int, int] = {}
    tm_scores: dict[int, float] = {}
    aligned_length: int | None = None
    rmsd: float | None = None
    sequence_identity: float | None = None

    lines = stdout.splitlines()
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue

        name_match = USAlignRunner._NAME_PATTERN.match(stripped)
        if name_match:
            names[int(name_match.group("index"))] = name_match.group("name")
            continue

        length_match = USAlignRunner._LENGTH_PATTERN.match(stripped)
        if length_match:
            lengths[int(length_match.group("index"))] = int(length_match.group("length"))
            continue

        summary_match = USAlignRunner._SUMMARY_PATTERN.match(stripped)
        if summary_match:
            aligned_length = int(summary_match.group("aligned"))
            rmsd = float(summary_match.group("rmsd"))
            sequence_identity = float(summary_match.group("seqid"))
            continue

        tm_match = USAlignRunner._TM_PATTERN.match(stripped)
        if tm_match:
            tm_scores[int(tm_match.group("index"))] = float(tm_match.group("score"))

    alignment_prediction = None
    alignment_markup = None
    alignment_reference = None
    for index, line in enumerate(lines):
        if line.startswith("(\":"):
            if index + 3 >= len(lines):
                break
            alignment_prediction = lines[index + 1].strip()
            alignment_markup = lines[index + 2].rstrip()
            alignment_reference = lines[index + 3].strip()
            break

    if (
        1 not in lengths
        or 2 not in lengths
        or aligned_length is None
        or rmsd is None
        or sequence_identity is None
        or 1 not in tm_scores
        or 2 not in tm_scores
    ):
        raise ToolExecutionError("US-align output was incomplete and could not be parsed.")

    return USAlignResult(
        reference=str(reference),
        prediction=str(prediction),
        reference_name=names.get(2),
        prediction_name=names.get(1),
        reference_length=lengths[2],
        prediction_length=lengths[1],
        aligned_length=aligned_length,
        rmsd=rmsd,
        sequence_identity=sequence_identity,
        tm_score_reference=tm_scores[2],
        tm_score_prediction=tm_scores[1],
        alignment_prediction=alignment_prediction,
        alignment_markup=alignment_markup,
        alignment_reference=alignment_reference,
        binary_path=str(binary),
        output_directory=None if output_dir is None else str(output_dir),
        reference_structure_output=None if reference_input is None else str(reference_input),
        prediction_structure_output=None if prediction_input is None else str(prediction_input),
        superposition_output=None if superposition_output is None else str(superposition_output),
        superposed_prediction_output=(
            None if superposed_prediction_output is None else str(superposed_prediction_output)
        ),
    )


def _extract_first_model(source: Path, output: Path) -> None:
    records: list[str] = []
    started = False

    for row in source.read_text(encoding="utf-8").splitlines():
        record = row[:6].strip()
        if record in {"ATOM", "HETATM"}:
            started = True
            records.append(row)
            continue
        if started and record == "TER":
            records.append(row)
            break

    if not records:
        raise ToolExecutionError(f"US-align produced '{source}', but no aligned ATOM records were found.")
    if records[-1][:3] != "TER":
        records.append("TER")
    records.append("END")

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(records) + "\n", encoding="utf-8")


def _css_color(name: str) -> str:
    if name == "reference":
        return "#2f6db3"
    return "#d97706"
