# rna-kit

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
![Biopython](https://img.shields.io/badge/BioPython-1.83%2B-2E8B57)
![CLI](https://img.shields.io/badge/CLI-rna--kit-1F6FEB)
![RNA](https://img.shields.io/badge/RNA-Structure%20Evaluation-CB6D51)
![HTML](https://img.shields.io/badge/Reports-HTML%20%2B%20JSON-0F766E)

**rna-kit** is a toolkit for RNA structure normalization, residue mapping, scoring, benchmarking, and HTML reporting.
It is designed for the common workflow: take a native RNA structure, take one or more predicted structures, and generate reproducible evaluation results without hand-editing scripts.

## What It Does

`rna-kit` currently supports:

- RNA PDB normalization
- mmCIF input support for evaluation and external-tool workflows
- Automatic residue mapping between native and predicted structures
- RMSD and P-value
- INF / DI metrics through `MC-Annotate`
- all-atom `lDDT` in pure Python
- per-residue local error reporting
- MolProbity-compatible geometry validation
- RNA secondary-structure extraction and comparison
- `US-align` integration for global RNA superposition and TM-score
- HTML outputs for lDDT, benchmark dashboards, secondary structure, US-align, and combined assessment
- batch benchmarking from direct inputs or CSV/JSON manifests

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

## Core Commands

These are the commands most users need:

```bash
rna-kit assess native.pdb prediction.pdb
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

## Batch Benchmarking

`benchmark` supports:

- direct input paths
- glob expansion
- JSON manifest
- CSV manifest

Supported manifest fields:

- `prediction` or `model`
- `native` or `reference`
- `label`
- `native_index`
- `prediction_index`
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
  }
]
```

CSV example:

```csv
label,native,native_index,prediction,prediction_index
method_a,examples/data/14_solution_0.pdb,examples/data/14_solution_0.index,examples/data/14_ChenPostExp_2.pdb,
```

## Optional Third-Party Tools

The core package depends only on `biopython`, but some workflows use external tools.

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

## Testing

Run the test suite with:

```bash
python -m pytest -q
```
