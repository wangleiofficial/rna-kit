# RNA Kit

<div class="hero">
  <p class="hero-badges">
    <img alt="Python" src="https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white">
    <img alt="Biopython" src="https://img.shields.io/badge/BioPython-1.83%2B-2E8B57">
    <img alt="CLI" src="https://img.shields.io/badge/CLI-rna--kit-1F6FEB">
    <img alt="RNA" src="https://img.shields.io/badge/RNA-Structure%20Evaluation-CB6D51">
    <img alt="HTML Reports" src="https://img.shields.io/badge/Reports-HTML%20%2B%20JSON-0F766E">
  </p>
  <p class="hero-lead">
    RNA Kit is a practical toolkit for RNA structure normalization, residue mapping,
    3D and 2D evaluation, benchmarking, and report generation.
  </p>
  <p class="hero-actions">
    <a class="md-button md-button--primary" href="#quick-start">Quick Start</a>
    <a class="md-button" href="platform-support/">Platform Support</a>
  </p>
</div>

## What It Covers

- `normalize` for cleaning input structures before downstream processing
- `repair` for filling missing RNA atoms with Arena
- `map` for residue correspondence inspection and sequence-guided alignment
- `assess` for one-shot evaluation with RMSD, eRMSD, INF, lDDT, and optional secondary-structure and MolProbity metrics
- `benchmark` for manifest-driven batch evaluation and HTML dashboards
- `us-align` for global RNA superposition plus browser-based 3D visualization

## Documentation Map

- [Getting Started](getting-started.md): installation, optional Phenix setup, and first commands
- [Core Concepts and Commands](core-concepts.md): what `normalize`, `map`, `assess`, `lddt`, `secondary-compare`, `us-align`, and `benchmark` actually do
- [Assessment and Mapping Reference](assessment-reference.md): `assess` output fields, `.index` behavior, sequence-guided mapping, and HTML report types
- [Workflows and Benchmarking](workflows-and-benchmarking.md): common workflows, manifest formats, and benchmark output structure
- [Tooling, API, and Project Layout](tooling-and-api.md): external tool resolution, Python API, release notes, acknowledgments, and repository structure

## Quick Start

Install the package online with `pip` directly from GitHub:

```bash
python -m pip install "git+https://github.com/wangleiofficial/rna-kit.git"
```

Run the main single-pair evaluation:

```bash
rna-kit assess \
  examples/data/14_solution_0.pdb \
  examples/data/14_ChenPostExp_2.pdb \
  --secondary-structure \
  --per-residue \
  --html-report examples/output/assessment.html
```

Build and preview the documentation locally:

```bash
git clone https://github.com/wangleiofficial/rna-kit.git
cd rna-kit
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[docs]'
mkdocs serve
```

The GitHub Pages workflow publishes the same site from the `docs/` directory.

## Main Workflows

<div class="grid cards" markdown>

-   `assess`

    ---

    The main evaluation pipeline for one native/prediction pair.

    Output:
    `RMSD`, `P-value`, `INF`, `lDDT`, optional secondary-structure metrics, optional MolProbity metrics, and HTML/JSON reports.

-   `benchmark`

    ---

    Batch evaluation for multiple predictions.

    Output:
    Summary table, benchmark dashboard, and one detail page per successful target.

-   `map`

    ---

    Inspection tool for residue correspondence.

    Use it when chain IDs differ, numbering is messy, or you want to verify what `assess` will compare.

-   `repair`

    ---

    Missing-atom repair with Arena.

    Use it directly, or enable it from `assess` / `benchmark` with `--repair-missing-atoms`.

</div>

## Core Commands

```bash
rna-kit assess native.pdb prediction.pdb
rna-kit ermsd native.pdb prediction.pdb
rna-kit assess native.pdb prediction.pdb --repair-missing-atoms
rna-kit lddt native.pdb prediction.pdb --html lddt.html
rna-kit secondary-compare native.pdb prediction.pdb --html secondary.html
rna-kit us-align native.pdb prediction.pdb --html us_align.html
rna-kit benchmark --manifest benchmark_manifest.json --html-report benchmark.html
rna-kit map native.pdb prediction.pdb --prediction-fasta prediction.fasta
rna-kit repair input.pdb repaired_output.pdb
rna-kit tools
```

## Sequence-Guided Mapping

Sequence hints are optional inputs for cases where structure residue names or numbering are unreliable.

They do **not** change coordinates or improve the model itself.
They only help RNA Kit decide which residues correspond across two structures.

Supported forms:

- FASTA file via `--native-fasta` or `--prediction-fasta`
- raw sequence string via `--native-sequence` or `--prediction-sequence`

Example:

```bash
rna-kit map native.pdb prediction_bad_names.pdb \
  --prediction-fasta examples/data/14_solution_0.fasta
```

## Reports and Visualization

RNA Kit produces browser-friendly outputs:

- per-residue lDDT HTML reports
- secondary-structure comparison pages
- `US-align` 3D structure pages
- benchmark dashboards with per-target detail pages

These reports are static files. You can open them directly in a browser or publish them from GitHub Pages.

## Platform Support

The library core runs on Linux, macOS, and Windows.
Bundled and external tools have different support levels.

See the full matrix here:

- [Platform Support](platform-support.md)

## Repository Links

- Repository: [github.com/wangleiofficial/rna-kit](https://github.com/wangleiofficial/rna-kit)
- README: [README.md](https://github.com/wangleiofficial/rna-kit/blob/master/README.md)
- Release notes template: [release-body.md](https://github.com/wangleiofficial/rna-kit/blob/master/.github/release-body.md)
