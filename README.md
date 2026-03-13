# rna-kit

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
![Biopython](https://img.shields.io/badge/BioPython-1.83%2B-2E8B57)
![CLI](https://img.shields.io/badge/CLI-rna--kit-1F6FEB)
![RNA](https://img.shields.io/badge/RNA-Structure%20Evaluation-CB6D51)
![HTML](https://img.shields.io/badge/Reports-HTML%20%2B%20JSON-0F766E)
[![Docs](https://img.shields.io/badge/Docs-GitHub%20Pages-1F883D?logo=github&logoColor=white)](https://wangleiofficial.github.io/rna-kit/)

**rna-kit** is a toolkit for RNA structure normalization, residue mapping, scoring, benchmarking, and HTML reporting.
It is designed for the common workflow: take a native RNA structure, take one or more predicted structures, and generate reproducible evaluation results without hand-editing scripts.

## Quick Navigation

- [What It Does](#what-it-does)
- [Platform Support](#platform-support)
- [Documentation Site](#documentation-site)
- [Installation](#installation)
- [Start Here](#start-here)
- [Core Concepts](#core-concepts)
- [Core Commands](#core-commands)
- [Important Commands](#what-each-important-command-does)
- [`repair`](#repair-command)
- [`normalize`](#normalize-command)
- [`map`](#map-command)
- [`assess`](#assess-command)
- [`lddt`](#lddt-command)
- [`secondary-compare`](#secondary-compare-command)
- [`us-align`](#us-align-command)
- [`benchmark`](#benchmark-command)
- [What `assess` Outputs](#assess-output-fields)
- [Do You Need an `.index` File?](#index-files)
- [Sequence-Guided Mapping](#sequence-guided-mapping)
- [HTML Outputs](#html-outputs)
- [Batch Benchmarking](#batch-benchmarking)
- [Optional Third-Party Tools](#optional-third-party-tools)
- [GitHub Releases](#github-releases)
- [Python API](#python-api)
- [Repository Layout](#repository-layout)
- [Testing](#testing)

<a id="what-it-does"></a>
## What It Does

`rna-kit` currently supports:

- RNA PDB normalization
- mmCIF input support for evaluation and external-tool workflows
- Automatic residue mapping between native and predicted structures
- FASTA/sequence-guided residue mapping when residue names or numbering are unreliable
- missing-atom repair through Arena
- RMSD and P-value
- INF / DI metrics through `MC-Annotate`
- all-atom `lDDT` in pure Python
- per-residue local error reporting
- MolProbity-compatible geometry validation
- RNA secondary-structure extraction and comparison
- `US-align` integration for global RNA superposition and TM-score
- HTML outputs for lDDT, benchmark dashboards, secondary structure, US-align, and combined assessment
- per-target benchmark detail pages linked from the batch dashboard
- batch benchmarking from direct inputs or CSV/JSON manifests
- GitHub release notes with a platform support matrix

<a id="platform-support"></a>
## Platform Support

### Core Runtime

![Linux](https://img.shields.io/badge/Linux-supported-2EA44F?logo=linux&logoColor=white)
![macOS](https://img.shields.io/badge/macOS-supported-000000?logo=apple&logoColor=white)
![Windows](https://img.shields.io/badge/Windows-supported-0078D4?logo=windows&logoColor=white)
![HTML Reports](https://img.shields.io/badge/HTML%20Reports-local%20browser%20support-0F766E?logo=html5&logoColor=white)

- The pure Python core works on Linux, macOS, and Windows
- `PDB` input is supported natively, and `mmCIF` is converted automatically when external tools require `PDB`
- HTML reports open locally in a browser without a running server

### Optional Tooling

![Arena](https://img.shields.io/badge/Arena-Linux%20%26%20macOS%20auto--build-1F6FEB)
![MC-Annotate](https://img.shields.io/badge/MC--Annotate-Linux%20bundled-8B5CF6)
![US-align](https://img.shields.io/badge/US--align-bundled%20on%20all%20platforms-C2410C)
![MolProbity](https://img.shields.io/badge/MolProbity%20%2F%20Phenix-external%20install-4D7C0F)

- `Arena`: auto-builds from source on Linux and macOS; on Windows, provide your own executable
- `MC-Annotate`: bundled and tested on Linux; on macOS and Windows, use precomputed `.mcout` files or provide your own binary
- `US-align`: bundled for Linux, macOS, and Windows
- `MolProbity / Phenix`: not bundled; supported on any platform where the command-line tools are installed

Detailed notes:
- [docs/platform-support.md](docs/platform-support.md)

<a id="documentation-site"></a>
## Documentation Site

The project now includes a dedicated documentation site built with `MkDocs Material`.

- source pages live in [docs/](docs/)
- the site configuration is in [mkdocs.yml](mkdocs.yml)
- GitHub Pages deployment is handled by [.github/workflows/pages.yml](.github/workflows/pages.yml)
- the published site URL is `https://wangleiofficial.github.io/rna-kit/`
- the workflow attempts to enable **GitHub Actions** deployment automatically on first run
- if your repository or organization blocks automatic Pages enablement, confirm **Settings > Pages > Build and deployment > Source = GitHub Actions**

To preview the site locally:

```bash
python -m pip install -e '.[docs]'
mkdocs serve
```

<a id="installation"></a>
## Installation

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
```

Main CLI:

```bash
rna-kit --help
```

### Optional: Install Phenix / MolProbity

`rna-kit` can run without Phenix, but `MolProbity`-related metrics require an external tool.
The most practical setup is to install **Phenix**, which provides command-line validation tools such as `phenix.clashscore` and, depending on the installation, `phenix.molprobity`.

Recommended setup:

1. Download Phenix from the official page:
   `https://phenix-online.org/download`
2. Install it.
   A typical non-interactive Linux or macOS shell install looks like:

```bash
bash phenix-<version>-<platform>.sh -b -p "$HOME/apps/phenix-<version>"
```

3. Load the Phenix shell environment:

```bash
source "$HOME/apps/phenix-<version>/phenix_env.sh"
```

4. Confirm the command is visible:

```bash
which phenix.clashscore
phenix.clashscore --help
```

5. Check that `rna-kit` can see it:

```bash
rna-kit tools
```

If you prefer, you can point `rna-kit` directly at the executable:

```bash
rna-kit molprobity model.pdb --molprobity "$(which phenix.clashscore)"
```

Notes:

- `rna-kit` checks `--molprobity`, then `RNA_KIT_MOLPROBITY`, then `PATH`
- a standalone MolProbity installation may also work if it exposes a compatible `molprobity` or `clashscore` command
- in practice, the most reliable path is Phenix

Official references:

- Phenix download: `https://phenix-online.org/download`
- Phenix installation and environment setup: `https://phenix-online.org/documentation/install-setup-run.html`
- MolProbity in Phenix: `https://www.phenix-online.org/documentation/reference/molprobity_tool.html`

<a id="start-here"></a>
## Start Here

If you only want one command to evaluate a prediction against a reference, use `assess`.

```bash
rna-kit assess \
  examples/data/14_solution_0.pdb \
  examples/data/14_ChenPostExp_2.pdb \
  --secondary-structure \
  --include-molprobity \
  --per-residue \
  --html-report examples/output/assessment.html
```

This command gives you:

- global 3D metrics
- optional secondary-structure metrics
- optional MolProbity geometry metrics
- per-residue lDDT and local error values
- an HTML report

For most users, this is enough. You do not need to run `map` first, and you usually do not need to run `normalize` first.
`assess` already tries automatic residue mapping internally, and if raw input preparation fails it will retry with temporary normalized copies of the input structures.

If you want only local quality visualization:

```bash
rna-kit lddt \
  examples/data/14_solution_0.pdb \
  examples/data/14_ChenPostExp_2.pdb \
  --html examples/output/lddt_report.html
```

If you want only secondary-structure comparison:

```bash
rna-kit secondary-compare \
  examples/data/14_solution_0.pdb \
  examples/data/14_ChenPostExp_2.pdb \
  --html examples/output/secondary_compare_fornac.html
```

If you want only 3D structural superposition:

```bash
rna-kit us-align \
  examples/data/14_solution_0.pdb \
  examples/data/14_ChenPostExp_2.pdb \
  --html examples/output/us_align.html
```

If you want a benchmark dashboard:

```bash
rna-kit benchmark \
  --manifest benchmark_manifest.json \
  --secondary-structure \
  --html-report examples/output/benchmark.html
```

If a structure is missing RNA atoms and you want to repair it before evaluation:

```bash
rna-kit repair \
  examples/data/14_ChenPostExp_2.pdb \
  examples/output/14_ChenPostExp_2.repaired.pdb
```

<a id="core-concepts"></a>
## Core Concepts

Before using the commands, it helps to separate four different tasks:

1. `normalize`: clean one structure file so that downstream tools can read it more consistently
2. `map`: figure out which residues in the prediction correspond to which residues in the reference
3. `assess`: compute scores once the residue correspondence is known
4. `benchmark`: run `assess` repeatedly for many predictions and summarize the results

`normalize` and `map` do not score a model by themselves. They prepare the inputs so the scoring commands behave predictably.

<a id="core-commands"></a>
## Core Commands

These are the commands most users need:

```bash
rna-kit assess native.pdb prediction.pdb
rna-kit repair input.pdb repaired_output.pdb
rna-kit lddt native.pdb prediction.pdb --html out.html
rna-kit secondary-compare native.pdb prediction.pdb --html out.html
rna-kit us-align native.pdb prediction.pdb --html out.html
rna-kit benchmark --manifest benchmark_manifest.json
rna-kit map native.pdb prediction.pdb
rna-kit molprobity prediction.pdb
rna-kit tools
```

Supporting commands:

```bash
rna-kit normalize input.pdb output.pdb
rna-kit extract input.pdb A:1:20,A:25:10 output.pdb
rna-kit rmsd native.pdb native.index prediction.pdb prediction.index
rna-kit inf native.pdb native.index prediction.pdb prediction.index
```

<a id="what-each-important-command-does"></a>
## What Each Important Command Does

<a id="repair-command"></a>
### `repair`

```bash
rna-kit repair input.pdb repaired_output.pdb
```

Purpose:

- repair missing RNA atoms with Arena
- write a new structure file that can be used for later evaluation
- keep the original input file unchanged

Use it when:

- a prediction is missing backbone or base atoms
- you want a repaired structure before running `assess`, `lddt`, or `secondary-compare`
- you want explicit control over whether repaired coordinates are used in evaluation

Notes:

- Arena mode `5` is the default in `rna-kit`
- mode `5` fills missing atoms without moving existing input atoms
- if Arena is not installed, `rna-kit` can auto-build it from the official source repository on Linux and macOS

<a id="normalize-command"></a>
### `normalize`

```bash
rna-kit normalize input.pdb output.pdb
```

Purpose:

- clean a structure file before evaluation
- standardize residue and atom naming where possible
- make downstream parsing more stable

Use it when:

- the input comes from different modeling tools or databases
- residue names or atom names are inconsistent
- you want a normalized file for reproducible preprocessing

Do not think of `normalize` as a scoring command. It prepares one structure file. It does not compare two structures.

<a id="map-command"></a>
### `map`

```bash
rna-kit map native.pdb prediction.pdb
```

Purpose:

- infer which residues in the native and prediction should be compared
- report chain mapping and generated index ranges
- let you inspect the comparison setup before computing metrics

Use it when:

- chain names differ between the two files
- residue numbering is messy
- you want to verify what `assess` will actually compare
- a score looks suspicious and you want to debug the mapping first

`map` does not compute RMSD, INF, or lDDT. It only answers: "what is being compared to what?"

<a id="assess-command"></a>
### `assess`

```bash
rna-kit assess native.pdb prediction.pdb
```

Purpose:

- run the main evaluation pipeline for one reference/prediction pair
- combine global 3D metrics, local quality, and optional secondary-structure and MolProbity metrics

Use it when:

- you want one command that gives the main evaluation result
- you want JSON output or an HTML report for a single prediction

What it already does for you:

- tries explicit indices if you provide them
- otherwise looks for sidecar `.index` files
- otherwise tries automatic residue mapping
- retries with temporary normalized inputs if structure preparation fails on the raw files

<a id="lddt-command"></a>
### `lddt`

```bash
rna-kit lddt native.pdb prediction.pdb --html lddt.html
```

Purpose:

- focus on local structural correctness
- report global lDDT plus per-residue local quality when requested

Use it when:

- you care more about local neighborhood quality than a single global RMSD
- you want per-nucleotide HTML visualization

<a id="secondary-compare-command"></a>
### `secondary-compare`

```bash
rna-kit secondary-compare native.pdb prediction.pdb --html secondary.html
```

Purpose:

- compare base-pairing patterns between native and prediction
- report precision, recall, F1, and Jaccard for RNA secondary structure

Use it when:

- you want to know whether stems and base pairs were recovered
- 2D structure is more important than atom-level 3D fit

<a id="us-align-command"></a>
### `us-align`

```bash
rna-kit us-align native.pdb prediction.pdb --html us_align.html
```

Purpose:

- run global structural superposition with `US-align`
- generate an interactive 3D web view

Use it when:

- you want a whole-structure alignment and TM-score style summary
- you want a browser view of the aligned structures

Unlike `assess`, `US-align` performs its own structural alignment and does not use `.index` files.

<a id="benchmark-command"></a>
### `benchmark`

```bash
rna-kit benchmark --manifest benchmark_manifest.json --html-report benchmark.html
```

Purpose:

- evaluate many predictions in one run
- generate a batch dashboard plus one detail page per successful entry

Use it when:

- you are comparing many models or many methods
- you want a leaderboard-like summary with linked detailed reports

<a id="assess-output-fields"></a>
## What `assess` Outputs

`rna-kit assess native.pdb prediction.pdb` returns JSON. The core fields are:

| Field | Meaning |
| --- | --- |
| `rmsd` | all-atom RMSD after superposition |
| `pvalue` | RMSD significance estimate |
| `deformation_index` | `RMSD / INF_ALL` |
| `inf_all` | overall interaction network fidelity |
| `inf_wc` | Watson-Crick interaction fidelity |
| `inf_nwc` | non-Watson-Crick interaction fidelity |
| `inf_stack` | stacking interaction fidelity |
| `lddt` | global all-atom lDDT |
| `lddt_evaluated_atoms` | number of atoms scored in lDDT |
| `lddt_evaluated_pairs` | number of local atom pairs scored |

With `--secondary-structure`, these fields are added:

| Field | Meaning |
| --- | --- |
| `secondary_structure_precision` | base-pair precision |
| `secondary_structure_recall` | base-pair recall |
| `secondary_structure_f1` | base-pair F1 |
| `secondary_structure_jaccard` | base-pair Jaccard score |

With `--per-residue`, `per_residue` is added. Each residue record contains:

- `native_chain`, `native_pos`, `native_nt`
- `prediction_chain`, `prediction_pos`, `prediction_nt`
- `matched_atoms`, `scored_atoms`
- `lddt`
- `local_rmsd`
- `mean_absolute_error`
- `max_absolute_error`

<a id="index-files"></a>
## Do You Need an `.index` File?

Usually, no.

For the high-level commands:

- `assess`
- `lddt`
- `map`
- `benchmark`

`rna-kit` behaves like this:

1. if you explicitly pass `--native-index` or `--prediction-index`, it uses them
2. otherwise it looks for sidecar files such as `model.index`
3. if no `.index` file exists, it tries automatic residue mapping

That means most users can start with:

```bash
rna-kit assess native.pdb prediction.pdb
```

Use an `.index` file only when you need strict control over the evaluation region, for example:

- evaluate only one RNA chain
- skip unresolved residues
- compare a discontinuous subset
- reproduce a benchmark with a fixed residue range

Current low-level commands `rmsd` and `inf` still require explicit index files.

Input format note:

- `assess`, `lddt`, `secondary-compare`, `benchmark`, `molprobity`, and `us-align` accept both `PDB` and `mmCIF`
- bundled external tools still run on `PDB`, so `rna-kit` converts `mmCIF` inputs automatically when needed

To inspect what `rna-kit` will compare before scoring, use:

```bash
rna-kit map native.pdb prediction.pdb
```

This prints the inferred `native_index`, `prediction_index`, matched residue count, and chain mapping.

<a id="sequence-guided-mapping"></a>
## Sequence-Guided Mapping

When residue names are noisy, non-canonical, or numbering is unreliable, you can provide sequence hints.

This is useful for cases such as:

- residue names rewritten as `UNK`, `MOD`, or other placeholders
- structures converted from pipelines that preserve coordinates but lose clean nucleotide labels
- messy residue numbering where sequence identity is more trustworthy than residue names

Commands that support sequence hints:

- `assess`
- `lddt`
- `secondary-compare`
- `map`
- `benchmark`

What a sequence hint does:

- it does not change the structure coordinates
- it does not improve the prediction itself
- it only helps `rna-kit` decide which residues correspond across two structures

In other words, sequence hints are only for mapping. They are useful when residue names in the structure file are not trustworthy enough for automatic matching.

You can pass either a FASTA file or a raw sequence string:

```bash
rna-kit map native.pdb prediction.pdb --prediction-fasta prediction.fasta
rna-kit assess native.pdb prediction.pdb --native-fasta native.fasta --prediction-fasta prediction.fasta
rna-kit lddt native.pdb prediction.pdb --prediction-sequence ACGUACGU...
```

Hint rules:

- a single-chain structure can use a single FASTA record
- a multi-chain structure can use one combined sequence whose total length matches all selected residues
- a multi-record FASTA can use record identifiers that match chain IDs, such as `>A` or `>chain:A`

`map` output includes `used_sequence_hints: true` when a sequence hint was applied.

Example intuition:

- if your prediction PDB has residues renamed to something uninformative like `MOD`
- but you still know the real RNA sequence
- then `--prediction-fasta prediction.fasta` tells `rna-kit` how to align residues by sequence instead of by unreliable residue labels

<a id="html-outputs"></a>
## HTML Outputs

### lDDT report

`lddt --html out.html` writes a standalone HTML report with:

- global lDDT summary
- per-nucleotide lDDT heatmap
- local RMSD bar chart
- residue table with local error values

Example:

- [examples/output/lddt_report.html](examples/output/lddt_report.html)

### Combined assessment report

`assess --html-report out.html` writes a combined HTML report containing:

- input and mapping summary
- global RMSD / INF / lDDT metrics
- optional MolProbity clashscore and geometry metrics
- optional secondary-structure summary
- per-residue lDDT visualization when available

### Benchmark dashboard

`benchmark --html-report out.html` writes a dashboard for many predictions.

- benchmark-level summary cards
- sortable batch output in CLI JSON
- per-model table with RMSD / INF / lDDT / SS F1 / MolProbity values
- one HTML detail page per successful entry
- links from the summary table into each detail page
- failed job reporting in the same page

### Secondary-structure web view

`secondary-structure --html out.html` and `secondary-compare --html out.html` write browser-based views for RNA secondary structure.

- single-structure dot-bracket inspection
- native vs prediction comparison
- mismatch highlighting for base pairs

Example:

- [examples/output/secondary_compare_fornac.html](examples/output/secondary_compare_fornac.html)

### US-align web view

`us-align --html out.html` writes a standalone 3D structure viewer powered by `3Dmol.js`.

- interactive 3D rotation and zoom
- aligned native and prediction structures
- US-align summary metrics in the page
- prediction structure colored by per-residue lDDT when residue mapping can be inferred

Example:

- [examples/output/us_align.html](examples/output/us_align.html)

## Typical Workflows

### 1. Evaluate one prediction against one native structure

```bash
rna-kit assess native.pdb prediction.pdb --html-report assessment.html
```

### 2. Inspect only residue mapping

```bash
rna-kit map native.pdb prediction.pdb
```

### 3. Generate a local quality report

```bash
rna-kit lddt native.pdb prediction.pdb --html lddt.html
```

### 4. Compare secondary structure

```bash
rna-kit secondary-compare native.pdb prediction.pdb --html secondary.html
```

### 5. Align whole structures with US-align

```bash
rna-kit us-align native.pdb prediction.pdb --html us_align.html
```

### 6. Run a batch benchmark

```bash
rna-kit benchmark \
  --manifest benchmark_manifest.json \
  --secondary-structure \
  --include-molprobity \
  --per-residue \
  --html-report benchmark.html
```

<a id="batch-benchmarking"></a>
## Batch Benchmarking

`benchmark` supports:

- direct input paths
- glob expansion
- JSON manifest
- CSV manifest
- optional per-entry sequence hints
- HTML dashboard plus linked per-entry detail pages

Supported manifest fields:

- `prediction` or `model`
- `native` or `reference`
- `label`
- `native_index`
- `prediction_index`
- `native_fasta` or `native_sequence`
- `prediction_fasta` or `prediction_sequence`
- `native_annotation`
- `prediction_annotation`

JSON example:

```json
[
  {
    "label": "method_a",
    "native": "examples/data/14_solution_0.pdb",
    "native_index": "examples/data/14_solution_0.index",
    "prediction": "examples/data/14_ChenPostExp_2.pdb"
  },
  {
    "label": "method_b_with_sequence_hint",
    "native": "examples/data/14_solution_0.pdb",
    "prediction": "examples/data/14_ChenPostExp_2.pdb",
    "prediction_fasta": "examples/data/14_ChenPostExp_2.fasta"
  }
]
```

CSV example:

```csv
label,native,native_index,prediction,prediction_index,prediction_fasta
method_a,examples/data/14_solution_0.pdb,examples/data/14_solution_0.index,examples/data/14_ChenPostExp_2.pdb,,
method_b_with_sequence_hint,examples/data/14_solution_0.pdb,,examples/data/14_ChenPostExp_2.pdb,,examples/data/14_ChenPostExp_2.fasta
```

If you write:

```bash
rna-kit benchmark --manifest benchmark_manifest.json --html-report benchmark.html
```

`rna-kit` generates:

- `benchmark.html` for the batch summary
- `benchmark_reports/` beside it, containing one detail page per successful entry

<a id="optional-third-party-tools"></a>
## Optional Third-Party Tools

The core package depends only on `biopython`, but some workflows use external tools.

### Arena

Used for:

- repairing missing RNA atoms
- generating a repaired structure file before evaluation

Resolution order:

1. `--arena`
2. `RNA_KIT_ARENA`
3. `PATH`
4. cached auto-built binary under `~/.cache/rna-kit/bin`

Current behavior:

- `rna-kit` does not bundle Arena binaries in the repository
- if no executable is available, `rna-kit repair` can auto-build Arena from the official source repository on Linux and macOS
- Windows users should provide a compiled Arena executable explicitly

Official source:

- `https://github.com/pylelab/Arena`

### MC-Annotate

Used for:

- `INF`
- `DI`
- secondary-structure extraction

Resolution order:

1. explicit annotation files such as `--annotation`, `--native-annotation`, `--prediction-annotation`
2. sidecar `.mcout`
3. `--mc-annotate`
4. `RNA_KIT_MC_ANNOTATE`
5. `PATH`
6. bundled `third_party/bin/MC-Annotate`

### MolProbity

Used for:

- clashscore
- MolProbity score
- bond / angle / RNA geometry outlier summaries when the selected binary reports them

Resolution order:

1. `--molprobity`
2. `RNA_KIT_MOLPROBITY`
3. `PATH`

### US-align

Used for:

- whole-structure superposition
- TM-score
- aligned structure export

Resolution order:

1. `--us-align`
2. `RNA_KIT_US_ALIGN`
3. `PATH`
4. bundled `third_party/bin/USalign*`

Check current tool availability with:

```bash
rna-kit tools
```

<a id="github-releases"></a>
## GitHub Releases

Tagged releases can publish:

- source and wheel artifacts
- release notes
- a platform support matrix

The repository includes:

- [docs/platform-support.md](docs/platform-support.md)
- [.github/release-body.md](.github/release-body.md)
- [.github/workflows/release.yml](.github/workflows/release.yml)

If you create a tag like `v0.1.0`, the release workflow can attach the current platform notes to the GitHub Release page.

<a id="python-api"></a>
## Python API

```python
from pathlib import Path

from rna_kit import calculate_assessment, calculate_us_align, normalize_structure

data_dir = Path("examples/data")

normalize_structure(
    data_dir / "14_solution_0.pdb",
    Path("normalized.pdb"),
)

assessment = calculate_assessment(
    data_dir / "14_solution_0.pdb",
    None,
    data_dir / "14_ChenPostExp_2.pdb",
    None,
    include_per_residue=True,
    include_secondary_structure=True,
    include_molprobity=False,
)

print(assessment.rmsd)
print(assessment.lddt)
print(assessment.secondary_structure_f1)

alignment = calculate_us_align(
    data_dir / "14_solution_0.pdb",
    data_dir / "14_ChenPostExp_2.pdb",
)

print(alignment.tm_score_reference)
print(alignment.tm_score_prediction)
```

For a larger example script, see [examples/basic_usage.py](examples/basic_usage.py).

<a id="repository-layout"></a>
## Repository Layout

```text
src/rna_kit/           main implementation
src/RNA_normalizer/    legacy compatibility layer
examples/data/         example inputs and precomputed annotations
examples/output/       generated HTML outputs
examples/golden/       notes for golden outputs
third_party/bin/       bundled executables
third_party/web/       bundled frontend assets
third_party/lib/       optional JAR dependencies
tests/                 automated test suite
```

## Acknowledgments

This project is maintained as **rna-kit**.
Parts of the implementation were adapted from the earlier **RNA_assessment** codebase.

<a id="testing"></a>
## Testing

Run the test suite with:

```bash
python -m pytest -q
```
