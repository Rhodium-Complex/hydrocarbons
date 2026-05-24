"""Main module for generating hydrocarbon structures and their SMILES representations."""
from concurrent.futures import ProcessPoolExecutor
import itertools
import os

import converter
import molecule_transformations
import structure_generator

MIN_CARBON = 2
MAX_CARBON = 7
SVG_GROUP_SIZE = 190
SVG_OUTPUT_DIR = "./image/"


def main():
    """Generate hydrocarbon structures and return their SMILES strings."""
    with (
        ProcessPoolExecutor() as dehydro_executor,
        ProcessPoolExecutor() as structure_executor,
    ):
        all_smiles_results = [
            "N#N",
            "N#N",
            "N#N",
            "C",
        ]

        for carbon_count in range(MIN_CARBON, MAX_CARBON + 1):
            current_carbon_structures = []
            print(carbon_count, ":", end="")
            for hydrogen_count in range(0, carbon_count * 2 + 3, 2)[::-1]:
                all_smiles_results.append("N#N")

                future_dehydro = dehydro_executor.map(
                    molecule_transformations.unique_dehydro_mols,
                    current_carbon_structures,
                )
                current_carbon_structures = list(future_dehydro)
                print(".", end="")

                future_structure = structure_executor.map(
                    structure_generator.build_structure,
                    structure_generator.build_carbon_hydrogen_combination(
                        carbon_count,
                        hydrogen_count,
                    ),
                )
                print(".", end="")

                for structures in future_structure:
                    current_carbon_structures += structures

                flattened_structures = itertools.chain.from_iterable(
                    current_carbon_structures
                )
                future_smiles = map(converter.mat2smiles, flattened_structures)
                results_before_adding = len(all_smiles_results)
                all_smiles_results += list(future_smiles)
                print(len(all_smiles_results) - results_before_adding, end=" ")

            print("")
    return all_smiles_results


if __name__ == "__main__":
    all_smiles_results = main()


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
