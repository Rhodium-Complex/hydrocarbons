"""Main module for generating hydrocarbon structures and their SMILES representations."""
import argparse
from pathlib import Path

import export_structures_pdf
import export_structures_svg
import generation_pipeline

MIN_CARBON = 2
MAX_CARBON = 9
EXPORT_NONE = "none"
EXPORT_PDF = "pdf"
EXPORT_SVG = "svg"


def main(
    min_carbon=MIN_CARBON,
    max_carbon=MAX_CARBON,
    workers=None,
    include_smiles=True,
    include_stereo=False,
    export_format=EXPORT_NONE,
    output=export_structures_pdf.DEFAULT_PDF_OUTPUT,
    output_prefix=export_structures_svg.DEFAULT_SVG_OUTPUT_PREFIX,
    include_tetrahedral_stereo=False,
):
    """Generate hydrocarbon structures and optionally return their SMILES strings."""
    if export_format != EXPORT_NONE:
        groups = generation_pipeline.run_export_smiles_groups(
            min_carbon=min_carbon,
            max_carbon=max_carbon,
            workers=workers,
            include_stereo=include_stereo,
            include_tetrahedral_stereo=include_tetrahedral_stereo,
        )
        try:
            if export_format == EXPORT_PDF:
                export_structures_pdf.export_formula_smiles_pdf(groups, output)
                print(f"wrote {output}", flush=True)
                return output
            output_paths = export_structures_svg.export_formula_smiles_svg_pages(
                groups,
                output_prefix,
            )
        except RuntimeError as exc:
            raise SystemExit(str(exc)) from exc
        for output_path in output_paths:
            print(f"wrote {output_path}", flush=True)
        return output_paths

    if include_smiles:
        return generation_pipeline.run_generation_smiles_groups(
            min_carbon=min_carbon,
            max_carbon=max_carbon,
            workers=workers,
            include_stereo=include_stereo,
            include_tetrahedral_stereo=include_tetrahedral_stereo,
            log_step=generation_pipeline.print_step_result,
        )
    generation_pipeline.run_generation(
        min_carbon=min_carbon,
        max_carbon=max_carbon,
        workers=workers,
        log_step=generation_pipeline.print_step_result,
    )
    return None


def parse_args():
    """Parse command-line arguments for the hydrocarbon generation script."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-carbon", type=int, default=MIN_CARBON)
    parser.add_argument("--max-carbon", type=int, default=MAX_CARBON)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--no-smiles", action="store_true")
    parser.add_argument("--include-stereo", action="store_true")
    parser.add_argument("--include-tetrahedral-stereo", action="store_true")
    parser.add_argument(
        "--export",
        choices=(EXPORT_NONE, EXPORT_PDF, EXPORT_SVG),
        default=EXPORT_NONE,
        help="Export generated structures instead of returning SMILES groups.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=export_structures_pdf.DEFAULT_PDF_OUTPUT,
        help="PDF output path used with --export pdf.",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=export_structures_svg.DEFAULT_SVG_OUTPUT_PREFIX,
        help="SVG output prefix used with --export svg.",
    )
    args = parser.parse_args()
    if args.no_smiles and args.export != EXPORT_NONE:
        parser.error("--no-smiles cannot be combined with --export")
    return args


if __name__ == "__main__":
    args = parse_args()
    main(
        min_carbon=args.min_carbon,
        max_carbon=args.max_carbon,
        workers=args.workers,
        include_smiles=not args.no_smiles,
        include_stereo=args.include_stereo,
        include_tetrahedral_stereo=args.include_tetrahedral_stereo,
        export_format=args.export,
        output=args.output,
        output_prefix=args.output_prefix,
    )
