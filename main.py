"""Main module for generating hydrocarbon structures and their SMILES representations."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import itertools
import os
import time

import converter
import molecule_transformations
import structure_generator

MIN_CARBON = 2
MAX_CARBON = 7
SVG_GROUP_SIZE = 190
SVG_OUTPUT_DIR = "./image/"


def _count_structures(structure_groups):
    return sum(len(structures) for structures in structure_groups)


def main(
    min_carbon=MIN_CARBON,
    max_carbon=MAX_CARBON,
    workers=None,
    include_smiles=True,
):
    """Generate hydrocarbon structures and optionally return their SMILES strings."""
    with (
        ProcessPoolExecutor(max_workers=workers) as dehydro_executor,
        ProcessPoolExecutor(max_workers=workers) as structure_executor,
    ):
        all_smiles_results = ["N#N", "N#N", "N#N", "C"] if include_smiles else []

        for carbon_count in range(min_carbon, max_carbon + 1):
            current_carbon_structures = []
            for hydrogen_count in range(0, carbon_count * 2 + 3, 2)[::-1]:
                step_start = time.perf_counter()
                if include_smiles:
                    all_smiles_results.append("N#N")

                dehydro_start = time.perf_counter()
                future_dehydro = dehydro_executor.map(
                    molecule_transformations.unique_dehydro_mols,
                    current_carbon_structures,
                )
                current_carbon_structures = list(future_dehydro)
                dehydro_seconds = time.perf_counter() - dehydro_start

                build_start = time.perf_counter()
                future_structure = structure_executor.map(
                    structure_generator.build_structure,
                    structure_generator.build_carbon_hydrogen_combination(
                        carbon_count,
                        hydrogen_count,
                    ),
                )

                for structures in future_structure:
                    current_carbon_structures += structures
                build_seconds = time.perf_counter() - build_start

                structure_count = _count_structures(current_carbon_structures)
                smiles_seconds = 0.0
                output_count = structure_count

                if include_smiles:
                    smiles_start = time.perf_counter()
                    flattened_structures = itertools.chain.from_iterable(
                        current_carbon_structures
                    )
                    future_smiles = map(converter.mat2smiles, flattened_structures)
                    results_before_adding = len(all_smiles_results)
                    all_smiles_results += list(future_smiles)
                    output_count = len(all_smiles_results) - results_before_adding
                    smiles_seconds = time.perf_counter() - smiles_start

                total_seconds = time.perf_counter() - step_start
                print(
                    f"C={carbon_count} H={hydrogen_count} "
                    f"dehydro={dehydro_seconds:.4f}s "
                    f"build={build_seconds:.4f}s "
                    f"smiles={smiles_seconds:.4f}s "
                    f"count={output_count} "
                    f"total={total_seconds:.4f}s",
                    flush=True,
                )
    return all_smiles_results


def parse_args():
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
