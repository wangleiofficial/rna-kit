from __future__ import annotations

from pathlib import Path

import pytest

from rna_kit import default_tool_registry
from rna_kit.exceptions import ToolResolutionError

from .conftest import PROJECT_ROOT


def test_tool_registry_reports_bundled_cssr_and_mc_annotate() -> None:
    registry = default_tool_registry()
    cssr = registry.status("cssr")
    mc_annotate = registry.status("mc_annotate")
    us_align = registry.status("us_align")

    assert cssr.available
    assert cssr.source == "bundled"
    assert cssr.binary_path is not None
    assert Path(cssr.binary_path).exists()
    assert mc_annotate.available
    assert mc_annotate.binary_path is not None
    assert us_align.available
    assert us_align.source == "bundled"
    assert us_align.binary_path is not None
    assert Path(us_align.binary_path).exists()


def test_tool_registry_raises_for_missing_required_binary(tmp_path: Path) -> None:
    registry = default_tool_registry()
    missing = tmp_path / "missing-tool"
    with pytest.raises(ToolResolutionError):
        registry.require_binary("molprobity", override=missing)
