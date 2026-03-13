from __future__ import annotations

import argparse
import glob
import json
from dataclasses import asdict
from pathlib import Path

from .api import extract_structure, normalize_structure
from .arena import ArenaRunner, repair_missing_atoms
from .benchmark import build_benchmark_jobs, describe_prepared_pair, load_benchmark_manifest, run_benchmark
from .exceptions import ManifestFormatError, RNAAssessmentError, ReportGenerationError
from .mc_annotate import MCAnnotateRunner, clone_with_annotation_overrides, existing_annotation_path
from .metrics import (
    PreparedStructurePair,
    calculate_assessment,
    calculate_assessment_from_prepared,
    calculate_ermsd,
    calculate_mcq,
    calculate_secondary_structure,
    calculate_secondary_structure_comparison,
    calculate_interaction_network_fidelity,
    calculate_lddt,
    calculate_lddt_from_prepared,
    calculate_rmsd,
    prepare_structure_pair,
    prepared_structure_pair_context,
)
from .molprobity import MolProbityRunner, calculate_molprobity
from .reports import (
    build_benchmark_report_document,
    build_assessment_report_document,
    build_secondary_structure_report_document,
    write_benchmark_html_report,
    write_assessment_html_report,
    write_lddt_html_report,
    write_report_json,
    write_secondary_structure_html_report,
)
from .tools import default_tool_registry
from .secondary_structure_web import (
    write_secondary_structure_comparison_html,
    write_secondary_structure_html,
)
from .usalign import USAlignRunner, calculate_us_align, write_us_align_html


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rna-kit")
    subparsers = parser.add_subparsers(dest="command", required=True)

    normalize_parser = subparsers.add_parser("normalize", help="Normalize a PDB file.")
    normalize_parser.add_argument("input")
    normalize_parser.add_argument("output")

    repair_parser = subparsers.add_parser("repair", help="Repair missing RNA atoms with Arena.")
    repair_parser.add_argument("input")
    repair_parser.add_argument("output")
    repair_parser.add_argument("--arena")
    repair_parser.add_argument(
        "--arena-option",
        type=int,
        default=5,
        help="Arena repair mode. Use 5 to fill missing atoms without moving existing input atoms.",
    )

    extract_parser = subparsers.add_parser("extract", help="Extract indexed residues from a PDB file.")
    extract_parser.add_argument("input")
    extract_parser.add_argument("residues")
    extract_parser.add_argument("output")

    rmsd_parser = subparsers.add_parser("rmsd", help="Calculate RMSD and p-value.")
    rmsd_parser.add_argument("native")
    rmsd_parser.add_argument("native_index")
    rmsd_parser.add_argument("prediction")
    rmsd_parser.add_argument("prediction_index")
    rmsd_parser.add_argument("--pvalue-mode", default="-", choices=["+", "-"])

    ermsd_parser = subparsers.add_parser("ermsd", help="Calculate RNA eRMSD.")
    _add_structure_pair_arguments(ermsd_parser)
    ermsd_parser.add_argument("--cutoff", type=float, default=2.4)

    mcq_parser = subparsers.add_parser("mcq", help="Calculate RNA MCQ.")
    _add_structure_pair_arguments(mcq_parser)
    mcq_parser.add_argument("--mcq-jar")

    inf_parser = subparsers.add_parser("inf", help="Calculate interaction network fidelity metrics.")
    inf_parser.add_argument("native")
    inf_parser.add_argument("native_index")
    inf_parser.add_argument("prediction")
    inf_parser.add_argument("prediction_index")
    inf_parser.add_argument("--mc-annotate")
    inf_parser.add_argument("--annotation-cache-dir")
    inf_parser.add_argument("--native-annotation")
    inf_parser.add_argument("--prediction-annotation")

    lddt_parser = subparsers.add_parser("lddt", help="Calculate all-atom lDDT.")
    _add_structure_pair_arguments(lddt_parser)
    lddt_parser.add_argument("--inclusion-radius", type=float, default=15.0)
    lddt_parser.add_argument("--per-residue", action="store_true")
    lddt_parser.add_argument("--html")
    lddt_parser.add_argument("--title")

    tools_parser = subparsers.add_parser(
        "tools",
        help="Show availability of optional third-party tools.",
    )
    tools_parser.add_argument("--arena")
    tools_parser.add_argument("--cssr")
    tools_parser.add_argument("--mc-annotate")
    tools_parser.add_argument("--molprobity")
    tools_parser.add_argument("--us-align")

    molprobity_parser = subparsers.add_parser(
        "molprobity",
        help="Validate RNA geometry with MolProbity-compatible tools.",
    )
    molprobity_parser.add_argument("input")
    molprobity_parser.add_argument("--molprobity")

    us_align_parser = subparsers.add_parser(
        "us-align",
        help="Align RNA structures with US-align and optionally export an interactive HTML viewer.",
    )
    us_align_parser.add_argument("native")
    us_align_parser.add_argument("prediction")
    us_align_parser.add_argument("--us-align")
    us_align_parser.add_argument("--output-dir")
    us_align_parser.add_argument("--html")
    us_align_parser.add_argument("--title")

    secondary_parser = subparsers.add_parser(
        "secondary-structure",
        help="Parse RNA secondary structure from a 3D structure with MC-Annotate.",
    )
    secondary_parser.add_argument("input")
    secondary_parser.add_argument("--index")
    secondary_parser.add_argument("--annotation")
    secondary_parser.add_argument("--mc-annotate")
    secondary_parser.add_argument("--annotation-cache-dir")
    secondary_parser.add_argument("--html")
    secondary_parser.add_argument("--title")

    secondary_compare_parser = subparsers.add_parser(
        "secondary-compare",
        help="Compare RNA secondary structure for a reference/prediction pair with MC-Annotate.",
    )
    _add_structure_pair_arguments(secondary_compare_parser)
    secondary_compare_parser.add_argument("--mc-annotate")
    secondary_compare_parser.add_argument("--annotation-cache-dir")
    secondary_compare_parser.add_argument("--native-annotation")
    secondary_compare_parser.add_argument("--prediction-annotation")
    secondary_compare_parser.add_argument("--html")
    secondary_compare_parser.add_argument("--title")
    secondary_compare_parser.add_argument("--json-report")
    secondary_compare_parser.add_argument("--html-report")

    assess_parser = subparsers.add_parser(
        "assess",
        help="Calculate RMSD, P-value, INF and lDDT for a reference/prediction pair.",
    )
    _add_structure_pair_arguments(assess_parser)
    assess_parser.add_argument("--pvalue-mode", default="-", choices=["+", "-"])
    assess_parser.add_argument("--mc-annotate")
    assess_parser.add_argument("--annotation-cache-dir")
    assess_parser.add_argument("--native-annotation")
    assess_parser.add_argument("--prediction-annotation")
    assess_parser.add_argument("--inclusion-radius", type=float, default=15.0)
    assess_parser.add_argument("--per-residue", action="store_true")
    assess_parser.add_argument("--secondary-structure", action="store_true")
    assess_parser.add_argument("--include-molprobity", action="store_true")
    assess_parser.add_argument("--repair-missing-atoms", action="store_true")
    assess_parser.add_argument("--arena")
    assess_parser.add_argument("--arena-option", type=int, default=5)
    assess_parser.add_argument("--ermsd-cutoff", type=float, default=2.4)
    assess_parser.add_argument("--include-mcq", action="store_true")
    assess_parser.add_argument("--mcq-jar")
    assess_parser.add_argument("--secondary-structure-html")
    assess_parser.add_argument("--json-report")
    assess_parser.add_argument("--html-report")
    assess_parser.add_argument("--molprobity")

    map_parser = subparsers.add_parser(
        "map",
        help="Infer residue mapping and generated index specifications for a structure pair.",
    )
    _add_structure_pair_arguments(map_parser)

    benchmark_parser = subparsers.add_parser(
        "benchmark",
        help="Run batch assessment for one reference structure against many predictions.",
    )
    benchmark_parser.add_argument("native", nargs="?")
    benchmark_parser.add_argument("predictions", nargs="*")
    benchmark_parser.add_argument("--native-index")
    benchmark_parser.add_argument("--native-fasta")
    benchmark_parser.add_argument("--native-sequence")
    benchmark_parser.add_argument("--prediction-glob")
    benchmark_parser.add_argument("--manifest")
    benchmark_parser.add_argument("--pvalue-mode", default="-", choices=["+", "-"])
    benchmark_parser.add_argument("--mc-annotate")
    benchmark_parser.add_argument("--annotation-cache-dir")
    benchmark_parser.add_argument("--native-annotation")
    benchmark_parser.add_argument("--inclusion-radius", type=float, default=15.0)
    benchmark_parser.add_argument("--per-residue", action="store_true")
    benchmark_parser.add_argument("--secondary-structure", action="store_true")
    benchmark_parser.add_argument("--include-molprobity", action="store_true")
    benchmark_parser.add_argument("--repair-missing-atoms", action="store_true")
    benchmark_parser.add_argument("--arena")
    benchmark_parser.add_argument("--arena-option", type=int, default=5)
    benchmark_parser.add_argument("--ermsd-cutoff", type=float, default=2.4)
    benchmark_parser.add_argument("--include-mcq", action="store_true")
    benchmark_parser.add_argument("--mcq-jar")
    benchmark_parser.add_argument("--molprobity")
    benchmark_parser.add_argument("--json-report")
    benchmark_parser.add_argument("--html-report")
    benchmark_parser.add_argument(
        "--sort-by",
        choices=[
            "input",
            "rmsd",
            "ermsd",
            "mcq",
            "pvalue",
            "inf_all",
            "lddt",
            "secondary_structure_f1",
            "molprobity_clashscore",
        ],
        default="input",
    )
    benchmark_parser.add_argument("--descending", action="store_true")

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        if args.command == "normalize":
            path = normalize_structure(args.input, args.output)
            print(json.dumps({"output": str(path)}, indent=2))
            return 0

        if args.command == "repair":
            result = repair_missing_atoms(
                args.input,
                args.output,
                runner=ArenaRunner(binary_path=args.arena),
                option=args.arena_option,
            )
            print(json.dumps(asdict(result), indent=2))
            return 0

        if args.command == "extract":
            path = extract_structure(args.input, args.residues, args.output)
            print(json.dumps({"output": str(path)}, indent=2))
            return 0

        if args.command == "tools":
            registry = default_tool_registry()
            statuses = registry.list_statuses(
                overrides={
                    "arena": getattr(args, "arena", None),
                    "cssr": getattr(args, "cssr", None),
                    "mc_annotate": getattr(args, "mc_annotate", None),
                    "molprobity": getattr(args, "molprobity", None),
                    "us_align": getattr(args, "us_align", None),
                }
            )
            print(json.dumps({"tools": [asdict(item) for item in statuses]}, indent=2))
            return 0

        if args.command == "molprobity":
            result = calculate_molprobity(
                args.input,
                runner=MolProbityRunner(binary_path=args.molprobity),
            )
            print(json.dumps(asdict(result), indent=2))
            return 0

        if args.command == "us-align":
            result = calculate_us_align(
                args.native,
                args.prediction,
                runner=_build_us_align_runner(args),
                output_dir=_resolve_us_align_output_dir(args),
            )
            payload = asdict(result)
            if args.html:
                per_residue = None
                try:
                    with prepared_structure_pair_context(
                        args.native,
                        None,
                        args.prediction,
                        None,
                        resolve_sidecar_indices=False,
                    ) as prepared:
                        lddt_result = calculate_lddt_from_prepared(
                            prepared,
                            include_per_residue=True,
                        )
                        per_residue = lddt_result.per_residue
                except RNAAssessmentError:
                    per_residue = None
                if result.superposed_prediction_output is None:
                    raise ReportGenerationError("US-align HTML output requires a persistent superposed structure.")
                reference_view = result.reference_structure_output or args.native
                html_path = write_us_align_html(
                    reference_view,
                    result.superposed_prediction_output,
                    args.html,
                    result=result,
                    title=args.title or "RNA Kit Structure Alignment",
                    per_residue=per_residue,
                )
                payload["html_output"] = str(html_path)
            print(json.dumps(payload, indent=2))
            return 0

        if args.command == "rmsd":
            result = calculate_rmsd(
                args.native,
                args.native_index,
                args.prediction,
                args.prediction_index,
                pvalue_mode=args.pvalue_mode,
            )
            print(json.dumps(asdict(result), indent=2))
            return 0

        if args.command == "ermsd":
            result = calculate_ermsd(
                args.native,
                args.native_index,
                args.prediction,
                args.prediction_index,
                cutoff=args.cutoff,
                native_sequence_hint=_sequence_hint_argument(args, "native"),
                prediction_sequence_hint=_sequence_hint_argument(args, "prediction"),
            )
            print(json.dumps(asdict(result), indent=2))
            return 0

        if args.command == "mcq":
            result = calculate_mcq(
                args.native,
                args.native_index,
                args.prediction,
                args.prediction_index,
                jar_path=args.mcq_jar,
                native_sequence_hint=_sequence_hint_argument(args, "native"),
                prediction_sequence_hint=_sequence_hint_argument(args, "prediction"),
            )
            print(json.dumps(asdict(result), indent=2))
            return 0

        if args.command == "lddt":
            result = calculate_lddt(
                args.native,
                args.native_index,
                args.prediction,
                args.prediction_index,
                inclusion_radius=args.inclusion_radius,
                include_per_residue=args.per_residue or bool(args.html),
                native_sequence_hint=_sequence_hint_argument(args, "native"),
                prediction_sequence_hint=_sequence_hint_argument(args, "prediction"),
            )
            payload = asdict(result)
            if not args.per_residue:
                payload["per_residue"] = None
            if args.html:
                html_path = write_lddt_html_report(
                    args.native,
                    args.prediction,
                    result,
                    args.html,
                    title=args.title or "RNA Kit lDDT Report",
                )
                payload["html_output"] = str(html_path)
            print(json.dumps(payload, indent=2))
            return 0

        if args.command == "secondary-structure":
            result = calculate_secondary_structure(
                args.input,
                args.index,
                runner=_build_secondary_structure_runner(args, input_file=args.input),
            )
            payload = asdict(result)
            if args.html:
                path = write_secondary_structure_html(result, args.html, title=args.title or "RNA Kit Secondary Structure")
                payload["html_output"] = str(path)
            print(json.dumps(payload, indent=2))
            return 0

        if args.command == "secondary-compare":
            result = calculate_secondary_structure_comparison(
                args.native,
                args.native_index,
                args.prediction,
                args.prediction_index,
                runner=_build_secondary_structure_runner(
                    args,
                    native_file=args.native,
                    prediction_file=args.prediction,
                ),
                native_sequence_hint=_sequence_hint_argument(args, "native"),
                prediction_sequence_hint=_sequence_hint_argument(args, "prediction"),
            )
            payload = asdict(result)
            if args.html:
                path = write_secondary_structure_comparison_html(
                    result,
                    args.html,
                    title=args.title or "RNA Kit Secondary Structure Comparison",
                )
                payload["html_output"] = str(path)
            if args.json_report or args.html_report:
                document = build_secondary_structure_report_document(
                    result,
                    reference=args.native,
                    prediction=args.prediction,
                    tool_statuses=_tool_statuses(args, include_mc_annotate=True),
                    artifacts={"html": str(path)} if args.html else None,
                )
                if args.json_report:
                    json_path = write_report_json(document, args.json_report)
                    payload["json_report_output"] = str(json_path)
                if args.html_report:
                    html_path = write_secondary_structure_html_report(document, args.html_report)
                    payload["html_report_output"] = str(html_path)
            print(json.dumps(payload, indent=2))
            return 0

        if args.command == "map":
            prepared = prepare_structure_pair(
                args.native,
                args.native_index,
                args.prediction,
                args.prediction_index,
                native_sequence_hint=_sequence_hint_argument(args, "native"),
                prediction_sequence_hint=_sequence_hint_argument(args, "prediction"),
            )
            entry = describe_prepared_pair(prepared, args.prediction, args.native)
            payload = asdict(entry)
            payload["used_sidecar_index"] = prepared.used_sidecar_index
            payload["used_sequence_hints"] = prepared.used_sequence_hints
            print(json.dumps(payload, indent=2))
            return 0

        annotator = _build_annotator(args)
        secondary_structure_runner = (
            _build_secondary_structure_runner(args, native_file=args.native, prediction_file=args.prediction)
            if getattr(args, "secondary_structure", False)
            else None
        )
        if args.command == "inf":
            result = calculate_interaction_network_fidelity(
                args.native,
                args.native_index,
                args.prediction,
                args.prediction_index,
                annotator=annotator,
                native_sequence_hint=_sequence_hint_argument(args, "native"),
                prediction_sequence_hint=_sequence_hint_argument(args, "prediction"),
            )
            print(json.dumps(asdict(result), indent=2))
            return 0

        if args.command == "benchmark":
            jobs = _collect_benchmark_jobs(
                native=args.native,
                native_index=args.native_index,
                native_sequence_hint=_sequence_hint_argument(args, "native"),
                predictions=args.predictions,
                prediction_glob=args.prediction_glob,
                manifest=args.manifest,
                native_annotation=args.native_annotation,
            )
            if not jobs:
                parser.exit(2, "benchmark requires at least one prediction path, --prediction-glob, or --manifest.\n")
            result = run_benchmark(
                native_file=args.native,
                pvalue_mode=args.pvalue_mode,
                annotator=annotator,
                inclusion_radius=args.inclusion_radius,
                jobs=jobs,
                include_per_residue=args.per_residue,
                include_secondary_structure=args.secondary_structure,
                secondary_structure_runner=secondary_structure_runner,
                include_molprobity=args.include_molprobity,
                molprobity_runner=MolProbityRunner(binary_path=args.molprobity) if args.include_molprobity else None,
                repair_missing_atoms=args.repair_missing_atoms,
                repair_runner=_build_arena_runner(args) if args.repair_missing_atoms else None,
                arena_option=args.arena_option,
                ermsd_cutoff=args.ermsd_cutoff,
                include_mcq=args.include_mcq,
                mcq_jar_path=args.mcq_jar,
            )
            result = _sort_benchmark_result(result, args.sort_by, args.descending)
            payload = asdict(result)
            if args.json_report or args.html_report:
                document = build_benchmark_report_document(
                    result,
                    tool_statuses=_tool_statuses(
                        args,
                        include_arena=args.repair_missing_atoms,
                        include_mc_annotate=True,
                        include_molprobity=args.include_molprobity,
                    ),
                )
                if args.json_report:
                    json_path = write_report_json(document, args.json_report)
                    payload["json_report_output"] = str(json_path)
                if args.html_report:
                    html_path = write_benchmark_html_report(document, args.html_report)
                    payload["html_report_output"] = str(html_path)
                    payload["detail_reports_dir"] = str(html_path.parent / f"{html_path.stem}_reports")
            print(json.dumps(payload, indent=2))
            return 0

        with prepared_structure_pair_context(
            args.native,
            args.native_index,
            args.prediction,
            args.prediction_index,
            native_sequence_hint=_sequence_hint_argument(args, "native"),
            prediction_sequence_hint=_sequence_hint_argument(args, "prediction"),
            repair_missing_atoms=args.repair_missing_atoms,
            repair_runner=_build_arena_runner(args) if args.repair_missing_atoms else None,
            arena_option=args.arena_option,
        ) as prepared:
            annotator = _adapt_annotator_for_prepared(args, annotator, prepared)
            secondary_structure_runner = _adapt_annotator_for_prepared(args, secondary_structure_runner, prepared)
            result = calculate_assessment_from_prepared(
                prepared,
                pvalue_mode=args.pvalue_mode,
                annotator=annotator,
                inclusion_radius=args.inclusion_radius,
                include_per_residue=args.per_residue or bool(args.html_report),
                include_secondary_structure=args.secondary_structure,
                secondary_structure_runner=secondary_structure_runner,
                include_molprobity=args.include_molprobity,
                molprobity_runner=MolProbityRunner(binary_path=args.molprobity) if args.include_molprobity else None,
                ermsd_cutoff=args.ermsd_cutoff,
                include_mcq=args.include_mcq,
                mcq_jar_path=args.mcq_jar,
            )
            payload = asdict(result)
            payload["used_sidecar_index"] = prepared.used_sidecar_index
            payload["used_sequence_hints"] = prepared.used_sequence_hints
            payload["used_normalized_inputs"] = prepared.used_normalized_inputs
            payload["used_repaired_inputs"] = prepared.used_repaired_inputs
            if not args.per_residue:
                payload["per_residue"] = None
            artifacts = {}
            if args.secondary_structure_html and result.secondary_structure is not None:
                html_path = write_secondary_structure_comparison_html(
                    result.secondary_structure,
                    args.secondary_structure_html,
                    title="RNA Kit Secondary Structure Comparison",
                )
                artifacts["secondary_structure_html"] = str(html_path)
                payload["secondary_structure_html_output"] = str(html_path)
            if args.json_report or args.html_report:
                document = build_assessment_report_document(
                    prepared,
                    result,
                    reference=args.native,
                    prediction=args.prediction,
                    tool_statuses=_tool_statuses(
                        args,
                        include_arena=args.repair_missing_atoms,
                        include_mc_annotate=True,
                        include_molprobity=args.include_molprobity,
                    ),
                    artifacts=artifacts or None,
                )
                if args.json_report:
                    json_path = write_report_json(document, args.json_report)
                    payload["json_report_output"] = str(json_path)
                if args.html_report:
                    html_path = write_assessment_html_report(document, args.html_report)
                    payload["html_report_output"] = str(html_path)
        print(json.dumps(payload, indent=2))
        return 0
    except RNAAssessmentError as exc:
        parser.exit(1, f"{exc}\n")


def _add_structure_pair_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("native")
    parser.add_argument("prediction")
    parser.add_argument("--native-index")
    parser.add_argument("--prediction-index")
    parser.add_argument("--native-fasta")
    parser.add_argument("--prediction-fasta")
    parser.add_argument("--native-sequence")
    parser.add_argument("--prediction-sequence")


def _build_annotator(args: argparse.Namespace) -> MCAnnotateRunner:
    annotation_overrides = {}
    if getattr(args, "native_annotation", None):
        annotation_overrides[args.native] = args.native_annotation
    if getattr(args, "prediction_annotation", None):
        annotation_overrides[args.prediction] = args.prediction_annotation

    return MCAnnotateRunner(
        binary_path=getattr(args, "mc_annotate", None),
        cache_dir=getattr(args, "annotation_cache_dir", None),
        annotation_overrides=annotation_overrides or None,
    )


def _build_secondary_structure_runner(
    args: argparse.Namespace,
    *,
    input_file: str | None = None,
    native_file: str | None = None,
    prediction_file: str | None = None,
) -> MCAnnotateRunner:
    annotation_overrides = {}
    annotation = getattr(args, "annotation", None)
    if input_file is not None and annotation is not None:
        annotation_overrides[input_file] = annotation
    if native_file is not None and getattr(args, "native_annotation", None):
        annotation_overrides[native_file] = getattr(args, "native_annotation")
    if prediction_file is not None and getattr(args, "prediction_annotation", None):
        annotation_overrides[prediction_file] = getattr(args, "prediction_annotation")
    return MCAnnotateRunner(
        binary_path=getattr(args, "mc_annotate", None),
        cache_dir=getattr(args, "annotation_cache_dir", None),
        annotation_overrides=annotation_overrides or None,
    )


def _build_us_align_runner(args: argparse.Namespace) -> USAlignRunner:
    return USAlignRunner(binary_path=getattr(args, "us_align", None))


def _build_arena_runner(args: argparse.Namespace) -> ArenaRunner:
    return ArenaRunner(binary_path=getattr(args, "arena", None))


def _tool_statuses(
    args: argparse.Namespace,
    *,
    include_arena: bool = False,
    include_mc_annotate: bool = False,
    include_molprobity: bool = False,
):
    registry = default_tool_registry()
    requested = []
    if include_arena:
        requested.append(("arena", getattr(args, "arena", None)))
    if include_mc_annotate:
        requested.append(("mc_annotate", getattr(args, "mc_annotate", None)))
    if include_molprobity:
        requested.append(("molprobity", getattr(args, "molprobity", None)))
    return tuple(registry.status(key, override=override) for key, override in requested)


def _adapt_annotator_for_prepared(
    args: argparse.Namespace,
    runner: MCAnnotateRunner | None,
    prepared: PreparedStructurePair,
) -> MCAnnotateRunner | None:
    if runner is None:
        return None
    aliases = _annotation_aliases_for_prepared(args, prepared, cache_dir=getattr(args, "annotation_cache_dir", None))
    return clone_with_annotation_overrides(runner, aliases) or runner


def _annotation_aliases_for_prepared(
    args: argparse.Namespace,
    prepared: PreparedStructurePair,
    *,
    cache_dir: str | Path | None,
) -> dict[str, str]:
    overrides: dict[str, str] = {}
    native_annotation = existing_annotation_path(
        args.native,
        explicit_annotation=getattr(args, "native_annotation", None),
        cache_dir=cache_dir,
    )
    if native_annotation is not None:
        overrides[prepared.native.pdb_file] = str(native_annotation)
    prediction_annotation = existing_annotation_path(
        args.prediction,
        explicit_annotation=getattr(args, "prediction_annotation", None),
        cache_dir=cache_dir,
    )
    if prediction_annotation is not None:
        overrides[prepared.prediction.pdb_file] = str(prediction_annotation)
    return overrides


def _resolve_us_align_output_dir(args: argparse.Namespace) -> Path | None:
    output_dir = getattr(args, "output_dir", None)
    if output_dir:
        return Path(output_dir)
    html_output = getattr(args, "html", None)
    if not html_output:
        return None
    html_path = Path(html_output)
    base = html_path.with_suffix("") if html_path.suffix else html_path
    return base.parent / f"{base.name}_assets"


def _collect_prediction_paths(predictions: list[str], prediction_glob: str | None) -> list[str]:
    collected = list(predictions)
    if prediction_glob:
        collected.extend(sorted(glob.glob(prediction_glob)))
    deduplicated: list[str] = []
    seen: set[str] = set()
    for item in collected:
        if item in seen:
            continue
        seen.add(item)
        deduplicated.append(item)
    return deduplicated


def _collect_benchmark_jobs(
    native: str | None,
    native_index: str | None,
    native_sequence_hint: str | None,
    predictions: list[str],
    prediction_glob: str | None,
    manifest: str | None,
    native_annotation: str | None,
):
    jobs = []
    if manifest:
        jobs.extend(load_benchmark_manifest(manifest))

    prediction_paths = _collect_prediction_paths(predictions, prediction_glob)
    if prediction_paths:
        if native is None:
            raise ManifestFormatError("benchmark path mode requires a reference structure.")
        jobs.extend(
            build_benchmark_jobs(
                native_file=native,
                predictions=prediction_paths,
                native_index=native_index,
                native_sequence_hint=native_sequence_hint,
            )
        )

    if native_annotation is not None:
        jobs = [
            job.__class__(
                prediction=job.prediction,
                native=job.native,
                label=job.label,
                native_index=job.native_index,
                prediction_index=job.prediction_index,
                native_sequence_hint=job.native_sequence_hint,
                prediction_sequence_hint=job.prediction_sequence_hint,
                native_annotation=native_annotation if job.native_annotation is None else job.native_annotation,
                prediction_annotation=job.prediction_annotation,
            )
            for job in jobs
        ]

    if native is not None:
        jobs = [
            job.__class__(
                prediction=job.prediction,
                native=native if job.native is None else job.native,
                label=job.label,
                native_index=native_index if job.native_index is None else job.native_index,
                prediction_index=job.prediction_index,
                native_sequence_hint=(
                    native_sequence_hint if job.native_sequence_hint is None else job.native_sequence_hint
                ),
                prediction_sequence_hint=job.prediction_sequence_hint,
                native_annotation=job.native_annotation,
                prediction_annotation=job.prediction_annotation,
            )
            for job in jobs
        ]
    return jobs


def _sequence_hint_argument(args: argparse.Namespace, prefix: str) -> str | None:
    sequence_value = getattr(args, f"{prefix}_sequence", None)
    if sequence_value is not None:
        return sequence_value
    return getattr(args, f"{prefix}_fasta", None)


def _sort_benchmark_result(result, sort_by: str, descending: bool):
    if sort_by == "input":
        return result

    def metric_value(entry):
        if entry.metrics is None:
            return (1, 0.0)
        if sort_by == "molprobity_clashscore":
            if entry.metrics.molprobity is None or entry.metrics.molprobity.clashscore is None:
                return (1, 0.0)
            value = entry.metrics.molprobity.clashscore
        else:
            value = getattr(entry.metrics, sort_by)
            if value is None:
                return (1, 0.0)
        sort_value = -value if descending else value
        return (0, sort_value)

    succeeded = [entry for entry in result.entries if entry.metrics is not None]
    failed = [entry for entry in result.entries if entry.metrics is None]
    succeeded.sort(key=metric_value)
    return result.__class__(
        reference=result.reference,
        total_predictions=result.total_predictions,
        succeeded=result.succeeded,
        failed=result.failed,
        entries=tuple(succeeded + failed),
    )


if __name__ == "__main__":
    raise SystemExit(main())
