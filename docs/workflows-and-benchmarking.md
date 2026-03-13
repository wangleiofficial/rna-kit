# Workflows and Benchmarking

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

## Benchmark Dashboard

`benchmark --html-report out.html` writes a dashboard for many predictions.

- benchmark-level summary cards
- sortable batch output in CLI JSON
- per-model table with RMSD / eRMSD / optional MCQ / INF / lDDT / SS F1 / MolProbity values
- one HTML detail page per successful entry
- links from the summary table into each detail page
- failed job reporting in the same page
