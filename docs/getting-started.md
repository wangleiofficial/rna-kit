# Getting Started

## Installation

For regular use, install `rna-kit` online with `pip` directly from GitHub:

```bash
python -m pip install "git+https://github.com/wangleiofficial/rna-kit.git"
```

This gives you the CLI immediately:

```bash
rna-kit --help
```

`pip install rna-kit` is not documented here because the package is not currently published on PyPI.

If you want a local development checkout instead:

```bash
git clone https://github.com/wangleiofficial/rna-kit.git
cd rna-kit
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
```

If you also want to build the documentation site locally:

```bash
python -m pip install -e '.[docs]'
```

## Optional: Install Phenix / MolProbity

`rna-kit` can run without Phenix, but `MolProbity`-related metrics require an external tool.
The most practical setup is to install **Phenix**, which provides command-line validation tools such as
`phenix.clashscore` and, depending on the installation, `phenix.molprobity`.

Recommended setup:

1. Download Phenix from the official page: `https://phenix-online.org/download`
2. Install it. A typical non-interactive Linux or macOS shell install looks like:

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
- optional MCQ when `--include-mcq` is enabled
- optional secondary-structure metrics
- optional MolProbity geometry metrics
- per-residue lDDT and local error values
- an HTML report

For most users, this is enough. You do not need to run `map` first, and you usually do not need to run
`normalize` first. `assess` already tries automatic residue mapping internally, and if raw input preparation
fails it retries with temporary normalized copies of the input structures.

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
