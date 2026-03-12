from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from rna_kit import calculate_rmsd

from .conftest import PROJECT_ROOT


def test_rna_kit_import_alias_exports_public_api() -> None:
    assert callable(calculate_rmsd)


def test_rna_kit_module_entrypoint_works() -> None:
    env = os.environ.copy()
    src_path = str(Path(PROJECT_ROOT) / "src")
    env["PYTHONPATH"] = src_path if not env.get("PYTHONPATH") else f"{src_path}:{env['PYTHONPATH']}"
    result = subprocess.run(
        [sys.executable, "-m", "rna_kit", "--help"],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )

    assert "rna-kit" in result.stdout
