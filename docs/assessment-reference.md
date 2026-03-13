# Assessment and Mapping Reference

## What `assess` Outputs

`rna-kit assess native.pdb prediction.pdb` returns JSON. The core fields are:

| Field | Meaning |
| --- | --- |
| `rmsd` | all-atom RMSD after superposition |
| `ermsd` | RNA eRMSD based on base-relative geometry |
| `mcq` | RNA MCQ torsion-based similarity metric, when `--include-mcq` is enabled |
| `pvalue` | RMSD significance estimate |
| `deformation_index` | `RMSD / INF_ALL` |
| `inf_all` | overall interaction network fidelity |
| `inf_wc` | Watson-Crick interaction fidelity |
| `inf_nwc` | non-Watson-Crick interaction fidelity |
| `inf_stack` | stacking interaction fidelity |
| `lddt` | global all-atom lDDT |
| `lddt_evaluated_atoms` | number of atoms scored in lDDT |
| `lddt_evaluated_pairs` | number of local atom pairs scored |
| `ermsd_evaluated_residues` | number of residues used for eRMSD |
| `mcq_evaluated_residues` | number of residues used for MCQ |

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

## HTML Outputs

### lDDT report

`lddt --html out.html` writes a standalone HTML report with:

- global lDDT summary
- per-nucleotide lDDT heatmap
- local RMSD bar chart
- residue table with local error values

Example:

- [examples/output/lddt_report.html](https://github.com/wangleiofficial/rna-kit/blob/master/examples/output/lddt_report.html)

### Combined assessment report

`assess --html-report out.html` writes a combined HTML report containing:

- input and mapping summary
- global RMSD / eRMSD / optional MCQ / INF / lDDT metrics
- optional MolProbity clashscore and geometry metrics
- optional secondary-structure summary
- per-residue lDDT visualization when available

### Secondary-structure web view

`secondary-structure --html out.html` and `secondary-compare --html out.html` write browser-based views for RNA secondary structure.

- single-structure dot-bracket inspection
- native vs prediction comparison
- mismatch highlighting for base pairs

Example:

- [examples/output/secondary_compare_fornac.html](https://github.com/wangleiofficial/rna-kit/blob/master/examples/output/secondary_compare_fornac.html)

### US-align web view

`us-align --html out.html` writes a standalone 3D structure viewer powered by `3Dmol.js`.

- interactive 3D rotation and zoom
- aligned native and prediction structures
- US-align summary metrics in the page
- prediction structure colored by per-residue lDDT when residue mapping can be inferred

Example:

- [examples/output/us_align.html](https://github.com/wangleiofficial/rna-kit/blob/master/examples/output/us_align.html)
