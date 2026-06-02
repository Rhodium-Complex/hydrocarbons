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
):
    """Generate hydrocarbon structures and optionally return their SMILES strings."""
    return generation_pipeline.run_generation(
        min_carbon=min_carbon,
        max_carbon=max_carbon,
        workers=workers,
        include_smiles=include_smiles,
        log_step=lambda result: print(
            generation_pipeline.format_step_result(result),
            flush=True,
        ),
    )


def parse_args():
    """Parse command-line arguments for the hydrocarbon generation script."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-carbon", type=int, default=MIN_CARBON)
    parser.add_argument("--max-carbon", type=int, default=MAX_CARBON)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--no-smiles", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    all_smiles_results = main(
        min_carbon=args.min_carbon,
        max_carbon=args.max_carbon,
        workers=args.workers,
        include_smiles=not args.no_smiles,
    )


# --- Example Output Counts (C: H=...) ---
# 2 :1 1 1 0
# 3 :1 2 3 2 1
# 4 :2 5 9 11 7 3
# 5 :3 10 26 40 40 21 6
# 6 :5 25 77 159 217 185 85 19
# 7 :9 56 222 574 1029 1229 920 356 50
# 8 :18 139 652 2069 4656 7396 7950 5289 1804 204

# 8/28 21:50-9/1 22:33
# 2 :...1 ...1 ...1 ...
# 3 :...1 ...2 ...3 ...2 ...1
# 4 :...2 ...5 ...9 ...11 ...7 ...3
# 5 :...3 ...10 ...26 ...40 ...40 ...21 ...6
# 6 :...5 ...25 ...77 ...159 ...217 ...185 ...85 ...19
# 7 :...9 ...56 ...222 ...575 ...1031 ...1230 ...920 ...356 ...50
# 8 :...18 ...139 ...654 ...2082 ...4679 ...7437 ...7982 ...5308 ...1804
# 9 :...35 ...338 ...1902 ...7244 ...19983 ...40139 ...57771 ...56437
# 10 :...75 ...852 ...5568 ...24938 ...81909 ...201578 ..369067
