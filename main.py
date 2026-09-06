"""Main module for generating hydrocarbon structures and their SMILES representations."""
import argparse
from pathlib import Path

import export_structures_pdf
import export_structures_svg
import generation_pipeline
from svg_progress import SvgProgressBar

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
    svg_workers=1,
):
    """Generate hydrocarbon structures and optionally return their SMILES strings."""
    if export_format == EXPORT_SVG:
        try:
            with SvgProgressBar(
                print_step=generation_pipeline.print_step_result,
            ) as progress, export_structures_svg.SvgPageWriter(
                output_prefix,
                on_progress=progress.update,
                workers=svg_workers,
            ) as writer:
                generation_pipeline.stream_export_smiles_groups(
                    min_carbon=min_carbon,
                    max_carbon=max_carbon,
                    workers=workers,
                    include_stereo=include_stereo,
                    include_tetrahedral_stereo=include_tetrahedral_stereo,
                    consume_group=writer.add_group,
                    log_step=progress.log_step,
                )
                return writer.finish()
        except RuntimeError as exc:
            raise SystemExit(str(exc)) from exc

    if export_format == EXPORT_PDF:
        groups = generation_pipeline.run_export_smiles_groups(
            min_carbon=min_carbon,
            max_carbon=max_carbon,
            workers=workers,
            include_stereo=include_stereo,
            include_tetrahedral_stereo=include_tetrahedral_stereo,
        )
        try:
            export_structures_pdf.export_formula_smiles_pdf(groups, output)
        except RuntimeError as exc:
            raise SystemExit(str(exc)) from exc
        print(f"wrote {output}", flush=True)
        return output

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
    parser.add_argument(
        "--svg-workers", type=int, default=1,
        help="SVG page-rendering processes, separate from --workers (default: 1).",
    )
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
    if args.svg_workers < 1:
        parser.error("--svg-workers must be at least 1")
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
        svg_workers=args.svg_workers,
    )
