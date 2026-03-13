from __future__ import annotations

import os
import platform
import shutil
from dataclasses import dataclass
from pathlib import Path

from .exceptions import ToolResolutionError


@dataclass(frozen=True)
class ToolPluginSpec:
    key: str
    display_name: str
    category: str
    env_var: str | None
    path_names: tuple[str, ...]
    bundled_paths: tuple[str, ...]
    supports_auto_download: bool = False
    notes: str | None = None


@dataclass(frozen=True)
class ToolStatus:
    key: str
    display_name: str
    category: str
    available: bool
    source: str | None
    binary_path: str | None
    supports_auto_download: bool
    notes: str | None = None


class ToolRegistry:
    def __init__(self, specs: tuple[ToolPluginSpec, ...]) -> None:
        self._specs = {spec.key: spec for spec in specs}

    def spec(self, key: str) -> ToolPluginSpec:
        return self._specs[key]

    def list_statuses(self, overrides: dict[str, str | Path | None] | None = None) -> tuple[ToolStatus, ...]:
        overrides = overrides or {}
        return tuple(self.status(key, override=overrides.get(key)) for key in sorted(self._specs))

    def status(self, key: str, override: str | Path | None = None) -> ToolStatus:
        spec = self.spec(key)
        binary, source = self._resolve_binary(key, override)
        return ToolStatus(
            key=spec.key,
            display_name=spec.display_name,
            category=spec.category,
            available=binary is not None,
            source=source,
            binary_path=None if binary is None else str(binary),
            supports_auto_download=spec.supports_auto_download,
            notes=spec.notes,
        )

    def require_binary(self, key: str, override: str | Path | None = None) -> Path:
        binary, _ = self._resolve_binary(key, override)
        if binary is None:
            spec = self.spec(key)
            env_names = _env_var_candidates(spec.env_var)
            env_hint = ""
            if env_names:
                env_hint = " or set " + " / ".join(f"'{name}'" for name in env_names)
            raise ToolResolutionError(
                f"{spec.display_name} is not available. Provide a compatible executable{env_hint}."
            )
        return binary

    def _resolve_binary(
        self,
        key: str,
        override: str | Path | None = None,
    ) -> tuple[Path | None, str | None]:
        spec = self.spec(key)
        if override is not None:
            candidate = Path(override)
            if candidate.is_file():
                return candidate, "override"
            return None, "override"

        for env_var in _env_var_candidates(spec.env_var):
            if env_var in os.environ:
                candidate = Path(os.environ[env_var])
                if candidate.is_file():
                    return candidate, "env"

        for path_name in spec.path_names:
            resolved = shutil.which(path_name)
            if resolved:
                return Path(resolved), "path"

        repository_root = Path(__file__).resolve().parents[2]
        for relative_path in spec.bundled_paths:
            candidate = repository_root / relative_path
            if candidate.is_file():
                return candidate, "bundled"

        return None, None


def default_tool_registry() -> ToolRegistry:
    arena_cached = (str(_arena_cache_binary_path()),)
    cssr_bundled = ("third_party/bin/CSSR", f"third_party/bin/{_cssr_platform_binary_name()}")
    us_align_bundled = ("third_party/bin/USalign", f"third_party/bin/{_us_align_platform_binary_name()}")
    specs = (
        ToolPluginSpec(
            key="arena",
            display_name="Arena",
            category="structure_repair",
            env_var="RNA_KIT_ARENA",
            path_names=("Arena",),
            bundled_paths=arena_cached,
            supports_auto_download=True,
            notes="Missing-atom repair for RNA structures; can auto-build from the official source repository.",
        ),
        ToolPluginSpec(
            key="cssr",
            display_name="CSSR",
            category="secondary_structure",
            env_var="RNA_KIT_CSSR",
            path_names=("CSSR",),
            bundled_paths=cssr_bundled,
            supports_auto_download=True,
            notes="Canonical RNA secondary-structure assignment.",
        ),
        ToolPluginSpec(
            key="mc_annotate",
            display_name="MC-Annotate",
            category="interaction_annotation",
            env_var="RNA_KIT_MC_ANNOTATE",
            path_names=("MC-Annotate", "MC-Annotate.pl"),
            bundled_paths=("third_party/bin/MC-Annotate", "MC-Annotate"),
            notes="Interaction network annotation for INF metrics.",
        ),
        ToolPluginSpec(
            key="molprobity",
            display_name="MolProbity",
            category="geometry_validation",
            env_var="RNA_KIT_MOLPROBITY",
            path_names=("phenix.molprobity", "phenix.clashscore", "molprobity", "clashscore"),
            bundled_paths=(),
            notes="Geometry validation and clash analysis.",
        ),
        ToolPluginSpec(
            key="us_align",
            display_name="US-align",
            category="structure_alignment",
            env_var="RNA_KIT_US_ALIGN",
            path_names=("USalign", "US-align"),
            bundled_paths=us_align_bundled,
            notes="Global RNA structure superposition and TM-score metrics.",
        ),
    )
    return ToolRegistry(specs)


def _env_var_candidates(env_var: str | None) -> tuple[str, ...]:
    if env_var is None:
        return ()
    return (env_var,)


def _cssr_platform_binary_name() -> str:
    system = platform.system().lower()
    machine = platform.machine().lower()
    if system == "darwin":
        return "CSSR-Mac64"
    if system == "linux" and machine in {"x86_64", "amd64"}:
        return "CSSR-Linux64"
    if system == "windows" and machine in {"x86_64", "amd64"}:
        return "CSSR-Win64.exe"
    return "CSSR"


def _us_align_platform_binary_name() -> str:
    system = platform.system().lower()
    machine = platform.machine().lower()
    if system == "darwin":
        return "USalign-Mac64"
    if system == "linux" and machine in {"x86_64", "amd64"}:
        return "USalign-Linux64"
    if system == "windows" and machine in {"x86_64", "amd64"}:
        return "USalign-Win64.exe"
    return "USalign"


def _arena_cache_binary_path() -> Path:
    system = platform.system().lower()
    machine = platform.machine().lower().replace(" ", "_")
    return Path.home() / ".cache" / "rna-kit" / "bin" / f"Arena-{system}-{machine}"
