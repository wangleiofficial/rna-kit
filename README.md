# rna-kit

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
![Biopython](https://img.shields.io/badge/BioPython-1.83%2B-2E8B57)
![CLI](https://img.shields.io/badge/CLI-rna--kit-1F6FEB)
![Secondary Structure](https://img.shields.io/badge/Secondary%20Structure-MC--Annotate%20%2B%20FornaC-CB6D51)
![Alignment](https://img.shields.io/badge/Alignment-US--align-6F42C1)

**rna-kit** is a practical toolkit for RNA structure normalization, alignment, evaluation, and reporting.
It packages the parts that are usually scattered across scripts and third-party tools into one reproducible workflow.

## 🧬 What rna-kit provides

- PDB normalization for RNA structures
- RMSD and P-value calculation on indexed or auto-mapped residues
- INF and DI metrics from `MC-Annotate`
- Pure Python all-atom `lDDT`
- RNA secondary-structure extraction and comparison
- `US-align` integration for global RNA superposition and TM-score
- HTML reports and browser-based visualization
- Batch benchmarking from direct inputs or CSV/JSON manifests

## ✨ Highlights

- **One CLI** for normalization, mapping, scoring, benchmarking, and reporting
- **Minimal manual setup** for bundled examples
- **Automatic residue mapping** when `.index` files are not provided
- **Web visualization** for secondary structure and 3D superposition outputs

## 🚀 Installation

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
```

The main CLI is now:

```bash
rna-kit --help
```

## ⚡ Quick Start

Run a combined 3D assessment:

```bash
rna-kit assess \
  examples/data/14_solution_0.pdb \
  examples/data/14_ChenPostExp_2.pdb \
  --per-residue \
  --secondary-structure
```

Compare RNA secondary structure and export a web view:

```bash
rna-kit secondary-compare \
  examples/data/14_solution_0.pdb \
  examples/data/14_ChenPostExp_2.pdb \
  --html examples/output/secondary_compare_fornac.html \
  --title "Native vs Prediction"
```

Align two RNA structures with `US-align` and export an interactive HTML viewer:

```bash
rna-kit us-align \
  examples/data/14_solution_0.pdb \
  examples/data/14_ChenPostExp_2.pdb \
  --html examples/output/us_align.html
```

If a matching `*.index` file exists next to the input PDB, `rna-kit` will use it automatically.
If not, the tool falls back to chain-aware residue mapping.

## 🛠 CLI Commands

### Core workflows

```bash
rna-kit normalize INPUT.pdb OUTPUT.pdb
rna-kit extract INPUT.pdb RESIDUES OUTPUT.pdb
rna-kit map native.pdb prediction.pdb
rna-kit assess native.pdb prediction.pdb --secondary-structure --per-residue
rna-kit benchmark native.pdb prediction_1.pdb prediction_2.pdb
```

### Individual metrics

```bash
rna-kit rmsd native.pdb native.index prediction.pdb prediction.index
rna-kit inf native.pdb native.index prediction.pdb prediction.index
rna-kit lddt native.pdb prediction.pdb --native-index native.index --prediction-index prediction.index
rna-kit secondary-structure model.pdb --html out.html
rna-kit secondary-compare native.pdb prediction.pdb --html out.html
rna-kit us-align native.pdb prediction.pdb --html out.html
rna-kit tools
```

### Reports

```bash
rna-kit assess native.pdb prediction.pdb \
  --secondary-structure \
  --secondary-structure-html examples/output/secondary.html \
  --json-report examples/output/assessment.json \
  --html-report examples/output/assessment.html
```

## 📊 Batch Benchmarking

`benchmark` supports direct inputs, glob expansion, and CSV/JSON manifests.

Example:

```bash
rna-kit benchmark \
  --manifest benchmark_manifest.json \
  --secondary-structure \
  --sort-by secondary_structure_f1 \
  --descending \
  --per-residue
```

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

## 🧪 Python API

```python
from pathlib import Path

from rna_kit import (
    calculate_assessment,
    calculate_secondary_structure,
    calculate_us_align,
    normalize_structure,
    prepare_structure_pair,
)

data_dir = Path("examples/data")
output_dir = Path("examples/output")
output_dir.mkdir(exist_ok=True)

normalize_structure(
    data_dir / "14_solution_0.pdb",
    output_dir / "14_solution_0.normalized.pdb",
)

assessment = calculate_assessment(
    data_dir / "14_solution_0.pdb",
    None,
    data_dir / "14_ChenPostExp_2.pdb",
    None,
    include_per_residue=True,
    include_secondary_structure=True,
)
print(assessment.rmsd, assessment.lddt, assessment.secondary_structure_f1)

secondary = calculate_secondary_structure(data_dir / "14_solution_0.pdb", None)
print(secondary.dot_bracket)

alignment = calculate_us_align(
    data_dir / "14_solution_0.pdb",
    data_dir / "14_ChenPostExp_2.pdb",
)
print(alignment.tm_score_reference, alignment.tm_score_prediction)

prepared = prepare_structure_pair(
    data_dir / "14_solution_0.pdb",
    None,
    data_dir / "14_ChenPostExp_2.pdb",
    None,
)
print(prepared.native_index, prepared.prediction_index)
```

## 🖼 Visual Outputs

### Secondary structure

`secondary-structure` and `secondary-compare` export a **FornaC-based HTML viewer**.

- Single-structure view for dot-bracket inspection
- Native vs prediction comparison view
- Toggle controls for numbering, node labels, links, and pseudoknot links
- Color-coded residue states for shared or mismatched base-pair assignments

The FornaC frontend is bundled inside the repository, so the exported HTML can be opened locally without additional downloads.

### 3D superposition

`us-align --html out.html` exports a standalone HTML viewer powered by `3Dmol.js`.
The command also writes aligned structure files into a sibling `*_assets/` directory unless you specify `--output-dir`.

## 🔧 Third-Party Tool Resolution

The core library depends only on `biopython`, but several optional workflows use external tools.

### MC-Annotate

Used for:

- `INF`
- `DI`
- secondary-structure extraction

Resolution order:

1. `--annotation`, `--native-annotation`, `--prediction-annotation`
2. sidecar `.mcout` next to the input PDB
3. `--mc-annotate`
4. `RNA_KIT_MC_ANNOTATE`
5. `PATH`
6. bundled `third_party/bin/MC-Annotate`

### US-align

Resolution order:

1. `--us-align`
2. `RNA_KIT_US_ALIGN`
3. `PATH`
4. bundled `third_party/bin/USalign*`

### CSSR

CSSR remains registered in the tool system, but the current secondary-structure workflow uses `MC-Annotate`.

## 📁 Repository Layout

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

## 🙏 Acknowledgments

This project is now maintained and presented as **rna-kit**.
Parts of the implementation were adapted from the earlier **RNA_assessment** codebase, which helped shape the initial normalization and evaluation workflow.

## ✅ Testing

Run the test suite with:

```bash
pytest
```

The repository also includes GitHub Actions coverage for unit tests and real bundled binary integration checks.
