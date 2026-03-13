# Core Concepts and Commands

## Core Concepts

Before using the commands, it helps to separate four different tasks:

1. `normalize`: clean one structure file so downstream tools can read it more consistently
2. `map`: figure out which residues in the prediction correspond to which residues in the reference
3. `assess`: compute scores once the residue correspondence is known
4. `benchmark`: run `assess` repeatedly for many predictions and summarize the results

`normalize` and `map` do not score a model by themselves. They prepare the inputs so the scoring commands behave
predictably.

## Core Commands

These are the commands most users need:

```bash
rna-kit assess native.pdb prediction.pdb
rna-kit ermsd native.pdb prediction.pdb
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

## What Each Important Command Does

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

`normalize` is not a scoring command. It prepares one structure file. It does not compare two structures.

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

### `ermsd`

```bash
rna-kit ermsd native.pdb prediction.pdb
```

Purpose:

- calculate RNA eRMSD using base-relative geometry
- compare nucleobase arrangement rather than only atom-level superposition

Use it when:

- you want an RNA-specific 3D similarity measure
- ordinary RMSD is too sensitive to global rigid-body differences
- you care about relative base organization

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
