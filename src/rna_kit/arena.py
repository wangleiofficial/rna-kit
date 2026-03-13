from __future__ import annotations

import platform
import shutil
import stat
import subprocess
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

from .exceptions import ToolExecutionError, ToolResolutionError
from .structures import prepare_external_structure_input
from .tools import default_tool_registry


@dataclass(frozen=True)
class ArenaRepairResult:
    input_file: str
    output_file: str
    option: int
    binary_path: str
    used_auto_build: bool


class ArenaRunner:
    _ARCHIVE_URL = "https://codeload.github.com/pylelab/Arena/zip/refs/heads/main"

    def __init__(
        self,
        binary_path: str | Path | None = None,
        build_dir: str | Path | None = None,
        auto_build: bool = True,
    ) -> None:
        self.binary_path = Path(binary_path) if binary_path else None
        self.build_dir = Path(build_dir) if build_dir else None
        self.auto_build = auto_build

    def resolve_binary(self) -> tuple[Path, bool]:
        registry = default_tool_registry()
        if self.binary_path is not None:
            return registry.require_binary("arena", override=self.binary_path), False

        status = registry.status("arena", override=self.binary_path)
        if status.binary_path is not None:
            return Path(status.binary_path), False

        cache_candidate = self._binary_cache_dir() / _arena_cache_binary_name()
        if cache_candidate.is_file():
            return cache_candidate, True

        if not self.auto_build:
            raise ToolResolutionError(
                "Arena is not available. Set 'RNA_KIT_ARENA', pass '--arena', "
                "or enable automatic source build."
            )

        return self._build_binary(), True

    def repair(
        self,
        input_file: str | Path,
        output_file: str | Path,
        *,
        option: int = 5,
    ) -> ArenaRepairResult:
        if option not in {1, 2, 3, 4, 5}:
            raise ToolExecutionError("Arena option must be one of 1, 2, 3, 4, or 5.")

        binary, used_auto_build = self.resolve_binary()
        source_path = Path(input_file)
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            if source_path.suffix.lower() in {".cif", ".mmcif"}:
                with TemporaryDirectory(prefix="rna-kit-arena-") as temp_dir:
                    prepared_input = prepare_external_structure_input(
                        source_path,
                        Path(temp_dir) / f"{source_path.stem}.pdb",
                    )
                    self._run(binary, prepared_input, output_path, option)
            else:
                self._run(binary, source_path, output_path, option)
        except OSError as exc:
            raise ToolExecutionError(f"Failed to execute Arena at '{binary}'.") from exc

        return ArenaRepairResult(
            input_file=str(source_path),
            output_file=str(output_path),
            option=option,
            binary_path=str(binary),
            used_auto_build=used_auto_build,
        )

    def _run(
        self,
        binary: Path,
        input_path: Path,
        output_path: Path,
        option: int,
    ) -> None:
        try:
            subprocess.run(
                [str(binary), str(input_path), str(output_path), str(option)],
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as exc:
            message = exc.stderr.strip() or exc.stdout.strip() or "unknown error"
            raise ToolExecutionError(f"Arena failed for '{input_path}': {message}") from exc

    def _build_binary(self) -> Path:
        if platform.system().lower() == "windows":
            raise ToolResolutionError(
                "Automatic Arena source builds are currently supported on Linux and macOS only. "
                "Provide a compiled Arena executable with '--arena' or 'RNA_KIT_ARENA'."
            )

        if shutil.which("make") is None:
            raise ToolResolutionError(
                "Arena automatic source build requires 'make'. Install build tools or provide '--arena'."
            )

        target_dir = self._binary_cache_dir()
        target_dir.mkdir(parents=True, exist_ok=True)
        target_path = target_dir / _arena_cache_binary_name()

        with TemporaryDirectory(prefix="rna-kit-arena-build-") as temp_dir:
            archive_path = Path(temp_dir) / "Arena.zip"
            extract_dir = Path(temp_dir) / "source"
            try:
                with urllib.request.urlopen(self._ARCHIVE_URL, timeout=120) as response:
                    archive_path.write_bytes(response.read())
            except OSError as exc:
                raise ToolResolutionError(
                    "Unable to download Arena source from the official repository. "
                    "Set 'RNA_KIT_ARENA' or pass '--arena'."
                ) from exc

            with zipfile.ZipFile(archive_path) as archive:
                archive.extractall(extract_dir)

            source_roots = [item for item in extract_dir.iterdir() if item.is_dir()]
            if not source_roots:
                raise ToolResolutionError("Downloaded Arena source archive did not contain a source directory.")
            source_root = source_roots[0]

            try:
                subprocess.run(
                    ["make", "Arena"],
                    cwd=source_root,
                    check=True,
                    capture_output=True,
                    text=True,
                )
            except subprocess.CalledProcessError as exc:
                message = exc.stderr.strip() or exc.stdout.strip() or "unknown error"
                raise ToolExecutionError(
                    f"Failed to build Arena from source: {message}"
                ) from exc

            built_binary = source_root / "Arena"
            if not built_binary.is_file():
                raise ToolResolutionError("Arena source build completed without producing an executable.")

            temp_binary = target_path.with_name(f"{target_path.name}.tmp")
            shutil.copy2(built_binary, temp_binary)
            temp_binary.chmod(temp_binary.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
            temp_binary.replace(target_path)

        return target_path

    def _binary_cache_dir(self) -> Path:
        if self.build_dir is not None:
            return self.build_dir
        return Path.home() / ".cache" / "rna-kit" / "bin"


def repair_missing_atoms(
    input_file: str | Path,
    output_file: str | Path,
    *,
    runner: ArenaRunner | None = None,
    option: int = 5,
) -> ArenaRepairResult:
    runner = runner or ArenaRunner()
    return runner.repair(input_file, output_file, option=option)


def _arena_cache_binary_name() -> str:
    system = platform.system().lower()
    machine = platform.machine().lower().replace(" ", "_")
    return f"Arena-{system}-{machine}"
