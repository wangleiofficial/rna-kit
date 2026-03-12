from __future__ import annotations

from pathlib import Path

from .alignment import ChainAlignment, StructureAlignment, StructureMatcher, infer_structure_alignment
from .benchmark import (
    BenchmarkEntry,
    BenchmarkJob,
    BenchmarkResult,
    ChainMappingResult,
    build_benchmark_jobs,
    describe_prepared_pair,
    load_benchmark_manifest,
    run_benchmark,
)
from .extraction import extract_pdb
from .metrics import (
    AssessmentResult,
    InteractionNetworkResult,
    LDDTResult,
    PreparedStructurePair,
    ResidueAssessment,
    RMSDResult,
    calculate_assessment,
    calculate_assessment_from_prepared,
    calculate_secondary_structure,
    calculate_secondary_structure_comparison,
    calculate_interaction_network_fidelity,
    calculate_lddt,
    calculate_rmsd,
    prepare_structure_pair,
)
from .normalization import PDBNormalizer
from .reports import (
    AssessmentReportDocument,
    ReportMetadata,
    SCHEMA_VERSION,
    SecondaryStructureReportDocument,
    build_assessment_report_document,
    build_secondary_structure_report_document,
    write_assessment_html_report,
    write_report_json,
    write_secondary_structure_html_report,
)
from .secondary_structure import (
    CSSRAnnotation,
    CSSRBasePairRecord,
    CSSRRunner,
    SecondaryStructureBasePair,
    SecondaryStructureComparisonResult,
    SecondaryStructureResult,
    SecondaryStructureVisualizationResult,
    render_secondary_structure_comparison_svg,
    render_secondary_structure_svg,
    write_secondary_structure_comparison_svg,
    write_secondary_structure_svg,
)
from .secondary_structure_web import (
    render_secondary_structure_comparison_component,
    render_secondary_structure_comparison_html,
    render_secondary_structure_html,
    render_secondary_structure_component,
    write_secondary_structure_comparison_html,
    write_secondary_structure_html,
)
from .tools import ToolPluginSpec, ToolRegistry, ToolStatus, default_tool_registry
from .usalign import USAlignResult, USAlignRunner, calculate_us_align, render_us_align_html, write_us_align_html


def normalize_structure(
    input_file: str | Path,
    output_file: str | Path,
    normalizer: PDBNormalizer | None = None,
) -> Path:
    normalizer = normalizer or PDBNormalizer.from_defaults()
    return normalizer.normalize_or_raise(input_file, output_file)


def extract_structure(
    input_file: str | Path,
    residue_ranges: str,
    output_file: str | Path,
) -> Path:
    return extract_pdb(input_file, residue_ranges, output_file)


__all__ = [
    "AssessmentResult",
    "BenchmarkEntry",
    "BenchmarkJob",
    "BenchmarkResult",
    "ChainAlignment",
    "ChainMappingResult",
    "InteractionNetworkResult",
    "LDDTResult",
    "PreparedStructurePair",
    "ResidueAssessment",
    "RMSDResult",
    "CSSRAnnotation",
    "CSSRBasePairRecord",
    "CSSRRunner",
    "AssessmentReportDocument",
    "ReportMetadata",
    "SCHEMA_VERSION",
    "SecondaryStructureBasePair",
    "SecondaryStructureComparisonResult",
    "SecondaryStructureResult",
    "SecondaryStructureReportDocument",
    "SecondaryStructureVisualizationResult",
    "StructureAlignment",
    "StructureMatcher",
    "ToolPluginSpec",
    "ToolRegistry",
    "ToolStatus",
    "USAlignResult",
    "USAlignRunner",
    "build_benchmark_jobs",
    "build_assessment_report_document",
    "build_secondary_structure_report_document",
    "calculate_assessment",
    "calculate_assessment_from_prepared",
    "calculate_secondary_structure",
    "calculate_secondary_structure_comparison",
    "calculate_interaction_network_fidelity",
    "calculate_lddt",
    "calculate_rmsd",
    "calculate_us_align",
    "describe_prepared_pair",
    "extract_structure",
    "default_tool_registry",
    "infer_structure_alignment",
    "load_benchmark_manifest",
    "normalize_structure",
    "prepare_structure_pair",
    "render_secondary_structure_comparison_svg",
    "render_secondary_structure_comparison_component",
    "render_secondary_structure_comparison_html",
    "render_secondary_structure_svg",
    "render_secondary_structure_component",
    "render_secondary_structure_html",
    "render_us_align_html",
    "run_benchmark",
    "write_assessment_html_report",
    "write_report_json",
    "write_secondary_structure_html_report",
    "write_secondary_structure_comparison_html",
    "write_secondary_structure_comparison_svg",
    "write_secondary_structure_html",
    "write_secondary_structure_svg",
    "write_us_align_html",
]
