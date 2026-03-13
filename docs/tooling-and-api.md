# Tooling, API, and Project Layout

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

### MCQ

Used for:

- RNA backbone similarity scoring
- optional `assess` and `benchmark` outputs
- standalone `mcq` command

Resolution order:

1. `java` from `PATH`
2. `--mcq-jar` if provided
3. bundled `third_party/lib/mcq.ws.client-0.0.1-SNAPSHOT-jar-with-dependencies.jar`

Current behavior:

- `rna-kit` bundles the MCQ client JAR in the repository
- MCQ still requires a working Java runtime
- the `tools` command does not validate MCQ because it is a JAR-plus-runtime workflow rather than a single binary

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

## GitHub Releases

Tagged releases can publish:

- source and wheel artifacts
- release notes
- a platform support matrix

The repository includes:

- [docs/platform-support.md](platform-support.md)
- [.github/release-body.md](https://github.com/wangleiofficial/rna-kit/blob/master/.github/release-body.md)
- [.github/workflows/release.yml](https://github.com/wangleiofficial/rna-kit/blob/master/.github/workflows/release.yml)

If you create a tag like `v0.1.0`, the release workflow can attach the current platform notes to the GitHub Release page.

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

For a larger example script, see [examples/basic_usage.py](https://github.com/wangleiofficial/rna-kit/blob/master/examples/basic_usage.py).

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
