This directory documents the repository's bundled golden artifacts.

Current golden assets live in `examples/data/`:

- `*.mcout`: precomputed MC-Annotate outputs for INF regression tests
- `*.cssr`: precomputed CSSR outputs for secondary-structure regression tests

These files are intentionally versioned so the test suite and example commands
can exercise stable, known-good annotations without requiring users to run
third-party binaries first.
