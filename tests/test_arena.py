from __future__ import annotations

from pathlib import Path

import pytest

from rna_kit import ArenaRunner, repair_missing_atoms
from rna_kit.arena import _arena_cache_binary_name
from rna_kit.exceptions import ToolResolutionError

from .conftest import DATA_DIR, write_mmcif_from_pdb


def test_arena_runner_repairs_structure_with_explicit_binary(tmp_path: Path) -> None:
    binary_path = _write_fake_arena_script(tmp_path)
    output_path = tmp_path / "repaired.pdb"

    result = repair_missing_atoms(
        DATA_DIR / "14_solution_0.pdb",
        output_path,
        runner=ArenaRunner(binary_path=binary_path),
        option=5,
    )

    assert output_path.exists()
    assert result.output_file == str(output_path)
    assert result.option == 5
    assert result.used_auto_build is False
    assert (tmp_path / "seen_option.txt").read_text(encoding="utf-8").strip() == "5"


def test_arena_runner_converts_mmcif_input_to_pdb(tmp_path: Path) -> None:
    mmcif_path = write_mmcif_from_pdb(DATA_DIR / "14_solution_0.pdb", tmp_path / "model.cif")
    binary_path = _write_fake_arena_script(tmp_path)

    repair_missing_atoms(
        mmcif_path,
        tmp_path / "repaired.pdb",
        runner=ArenaRunner(binary_path=binary_path),
    )

    seen_path = (tmp_path / "seen_input.txt").read_text(encoding="utf-8").strip()
    assert seen_path.endswith(".pdb")


def test_arena_runner_can_auto_build_when_binary_is_missing(monkeypatch, tmp_path: Path) -> None:
    runner = ArenaRunner(build_dir=tmp_path)
    built_binary = tmp_path / _arena_cache_binary_name()
    built_binary.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    built_binary.chmod(0o755)

    monkeypatch.setattr(runner, "_build_binary", lambda: built_binary)

    resolved, used_auto_build = runner.resolve_binary()

    assert resolved == built_binary
    assert used_auto_build is True


def test_arena_runner_does_not_fallback_when_override_is_missing(tmp_path: Path) -> None:
    runner = ArenaRunner(binary_path=tmp_path / "missing-arena")

    with pytest.raises(ToolResolutionError):
        runner.resolve_binary()


def _write_fake_arena_script(tmp_path: Path) -> Path:
    script_path = tmp_path / "fake_arena.py"
    seen_input = tmp_path / "seen_input.txt"
    seen_option = tmp_path / "seen_option.txt"
    script_path.write_text(
        "\n".join(
            [
                "#!/usr/bin/env python3",
                "import sys",
                "from pathlib import Path",
                "input_path = Path(sys.argv[1])",
                "output_path = Path(sys.argv[2])",
                "option = sys.argv[3]",
                f"Path({str(seen_input)!r}).write_text(str(input_path), encoding='utf-8')",
                f"Path({str(seen_option)!r}).write_text(str(option), encoding='utf-8')",
                "output_path.write_text(input_path.read_text(encoding='utf-8'), encoding='utf-8')",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    script_path.chmod(0o755)
    return script_path
