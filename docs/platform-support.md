# Platform Support

## Core Library

| Component | Linux | macOS | Windows | Notes |
| --- | --- | --- | --- | --- |
| Python package | yes | yes | yes | Pure Python core depends on `biopython` |
| PDB input | yes | yes | yes | Native support |
| mmCIF input | yes | yes | yes | Converted to PDB automatically when required by external tools |
| HTML reports | yes | yes | yes | Open locally in a browser |

## Bundled Tools

| Tool | Linux | macOS | Windows | Notes |
| --- | --- | --- | --- | --- |
| Arena | auto-build from source | auto-build from source | manual executable recommended | `rna-kit` can build Arena from the official repository on first use of `repair` |
| MC-Annotate | bundled and tested | not bundled for native execution | not bundled | Use precomputed `.mcout` or provide your own executable on non-Linux platforms |
| US-align | bundled | bundled | bundled | `rna-kit` resolves platform-specific binaries automatically |
| CSSR | bundled | bundled | bundled | Currently optional; secondary-structure workflow uses MC-Annotate by default |

## External Tools

| Tool | Linux | macOS | Windows | Notes |
| --- | --- | --- | --- | --- |
| Phenix / MolProbity | supported if installed | supported if installed | supported if installed | `rna-kit` does not bundle Phenix; set `RNA_KIT_MOLPROBITY` or use `--molprobity` |

## CI Coverage

- Unit tests run on `ubuntu-latest` and `macos-latest`
- Real bundled-tool checks currently run in CI for:
  - `US-align`
  - bundled `MC-Annotate` on Linux only

## Practical Notes

- For the widest compatibility, use the core Python metrics plus precomputed annotation files
- For missing-atom repair, `rna-kit repair` can auto-build Arena on Linux and macOS if local build tools are installed
- For geometry validation, install Phenix and expose `phenix.clashscore` or `phenix.molprobity`
- For reproducible releases, prefer opening the GitHub Release page and reading the platform notes attached to that version
