"""Main module for generating hydrocarbon structures and their SMILES representations."""
import argparse

import generation_pipeline

MIN_CARBON = 2
MAX_CARBON = 9


def main(
    min_carbon=MIN_CARBON,
    max_carbon=MAX_CARBON,
    workers=None,
    include_smiles=True,
    include_stereo=False,
):
    """Generate hydrocarbon structures and optionally return their SMILES strings."""
    return generation_pipeline.run_generation(
        min_carbon=min_carbon,
        max_carbon=max_carbon,
        workers=workers,
        include_smiles=include_smiles,
        include_stereo=include_stereo,
        log_step=generation_pipeline.print_step_result,
    )


def parse_args():
    """Parse command-line arguments for the hydrocarbon generation script."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-carbon", type=int, default=MIN_CARBON)
    parser.add_argument("--max-carbon", type=int, default=MAX_CARBON)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--no-smiles", action="store_true")
    parser.add_argument("--include-stereo", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        min_carbon=args.min_carbon,
        max_carbon=args.max_carbon,
        workers=args.workers,
        include_smiles=not args.no_smiles,
        include_stereo=args.include_stereo,
    )


