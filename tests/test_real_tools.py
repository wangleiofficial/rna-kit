from __future__ import annotations

import os
import platform
from pathlib import Path

import pytest

from rna_kit import calculate_secondary_structure, calculate_us_align
from rna_kit.mc_annotate import MCAnnotateRunner

from .conftest import DATA_DIR, PROJECT_ROOT


pytestmark = pytest.mark.skipif(
    os.environ.get("RNA_KIT_RUN_REAL_TOOLS") != "1",
    reason="Real bundled-binary checks are only enabled in CI.",
)


def test_secondary_structure_uses_precomputed_mc_annotate_output() -> None:
    result = calculate_secondary_structure(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
    )

    assert result.base_pair_count > 0
    assert result.dot_bracket


def test_bundled_us_align_binary_runs(tmp_path: Path) -> None:
    result = calculate_us_align(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_ChenPostExp_2.pdb",
        output_dir=tmp_path / "us_align",
    )

    assert result.aligned_length > 0
    assert result.tm_score_reference > 0.0
    assert result.superposed_prediction_output is not None
    assert Path(result.superposed_prediction_output).exists()


@pytest.mark.skipif(platform.system() != "Linux", reason="Bundled MC-Annotate only runs on Linux.")
def test_bundled_mc_annotate_binary_runs() -> None:
    runner = MCAnnotateRunner(binary_path=PROJECT_ROOT / "third_party" / "bin" / "MC-Annotate")
    result = runner.load(DATA_DIR / "14_solution_0.pdb")

    assert len(result.interactions) > 0
