from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from rna_kit import (
    MCAnnotateRunner,
    USAlignRunner,
    calculate_assessment_from_prepared,
    calculate_us_align,
    describe_prepared_pair,
    normalize_structure,
    prepare_structure_pair,
    write_us_align_html,
)


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parent
    data_dir = root / "data"

    parser = argparse.ArgumentParser(
        description="Evaluate an RNA prediction against a reference structure.",
    )
    parser.add_argument(
        "reference",
        nargs="?",
        default=data_dir / "14_solution_0.pdb",
        type=Path,
        help="Reference PDB file.",
    )
    parser.add_argument(
        "prediction",
        nargs="?",
        default=data_dir / "14_ChenPostExp_2.pdb",
        type=Path,
        help="Predicted PDB file.",
    )
    parser.add_argument("--reference-index", type=Path)
    parser.add_argument("--prediction-index", type=Path)
    parser.add_argument("--reference-annotation", type=Path)
    parser.add_argument("--prediction-annotation", type=Path)
    parser.add_argument("--mc-annotate", type=Path)
    parser.add_argument("--annotation-cache-dir", type=Path)
    parser.add_argument("--normalize-reference", action="store_true")
    parser.add_argument("--inclusion-radius", type=float, default=15.0)
    parser.add_argument("--per-residue", action="store_true")
    parser.add_argument("--secondary-structure", action="store_true")
    parser.add_argument("--include-us-align", action="store_true")
    parser.add_argument("--us-align", type=Path)
    parser.add_argument("--us-align-output-dir", type=Path)
    parser.add_argument("--us-align-html", type=Path)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parent
    output_dir = root / "output"
    output_dir.mkdir(exist_ok=True)

    original_reference_path = args.reference
    reference_path = original_reference_path
    if args.normalize_reference:
        reference_path = normalize_structure(
            original_reference_path,
            output_dir / f"{original_reference_path.stem}.normalized.pdb",
        )

    annotation_overrides = {}
    if args.reference_annotation:
        annotation_overrides[reference_path] = args.reference_annotation
    elif args.normalize_reference:
        original_annotation = original_reference_path.with_name(f"{original_reference_path.name}.mcout")
        if original_annotation.exists():
            annotation_overrides[reference_path] = original_annotation
    if args.prediction_annotation:
        annotation_overrides[args.prediction] = args.prediction_annotation

    annotator = MCAnnotateRunner(
        binary_path=args.mc_annotate,
        cache_dir=args.annotation_cache_dir,
        annotation_overrides=annotation_overrides or None,
    )
    prepared = prepare_structure_pair(
        reference_path,
        args.reference_index,
        args.prediction,
        args.prediction_index,
    )
    result = calculate_assessment_from_prepared(
        prepared,
        annotator=annotator,
        pvalue_mode="-",
        inclusion_radius=args.inclusion_radius,
        include_per_residue=args.per_residue,
        include_secondary_structure=args.secondary_structure,
        secondary_structure_runner=annotator if args.secondary_structure else None,
    )
    description = describe_prepared_pair(prepared, args.prediction, reference_path)
    include_us_align = args.include_us_align or args.us_align_output_dir is not None or args.us_align_html is not None
    us_align_result = None
    us_align_html_output = None
    if include_us_align:
        us_align_output_dir = args.us_align_output_dir
        if us_align_output_dir is None and args.us_align_html is not None:
            us_align_output_dir = output_dir / "us_align_assets"
        us_align_result = calculate_us_align(
            reference_path,
            args.prediction,
            runner=USAlignRunner(binary_path=args.us_align),
            output_dir=us_align_output_dir,
        )
        if args.us_align_html is not None:
            if us_align_result.superposed_prediction_output is None:
                raise RuntimeError("US-align HTML output requires a persistent output directory.")
            reference_view = us_align_result.reference_structure_output or str(reference_path)
            html_path = write_us_align_html(
                reference_view,
                us_align_result.superposed_prediction_output,
                args.us_align_html,
                result=us_align_result,
            )
            us_align_html_output = str(html_path)

    print(
        json.dumps(
            {
                "reference": str(reference_path),
                "prediction": str(args.prediction),
                "native_index": description.native_index,
                "prediction_index": description.prediction_index,
                "matched_residues": description.matched_residues,
                "chain_mappings": asdict(description)["chain_mappings"],
                "metrics": asdict(result),
                "us_align": None if us_align_result is None else asdict(us_align_result),
                "us_align_html_output": us_align_html_output,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
