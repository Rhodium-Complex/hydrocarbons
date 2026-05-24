"""Benchmark generation and uniqueness filtering of hydrocarbon structures."""
import argparse
from pathlib import Path
import sys
import time
import tracemalloc

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import graph_utils
import molecule
import molecule_transformations
import numpy as np
import structure_generator


def benchmark_generation(max_carbon: int) -> None:
    """Measure generated structure counts, runtime, and peak memory."""
    for carbon_count in range(2, max_carbon + 1):
        total_inputs = 0
        total_structures = 0
        tracemalloc.start()
        start = time.perf_counter()
        for hydrogen_count in range(0, carbon_count * 2 + 3, 2)[::-1]:
            combinations = list(
                structure_generator.build_carbon_hydrogen_combination(
                    carbon_count,
                    hydrogen_count,
                )
            )
            total_inputs += len(combinations)
            for combination in combinations:
                total_structures += len(
                    structure_generator.build_structure(combination)
                )
        elapsed = time.perf_counter() - start
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        print(
            f"C={carbon_count} inputs={total_inputs} "
            f"unique_outputs={total_structures} seconds={elapsed:.4f} "
            f"peak_bytes={peak}"
        )


def representative_molecules(
    carbon_count: int,
    hydrogen_count: int,
    min_candidates: int,
) -> list:
    """Return a representative molecule set for uniqueness benchmarks."""
    for combination in structure_generator.build_carbon_hydrogen_combination(
        carbon_count,
        hydrogen_count,
    ):
        carbon_type_list = np.asarray(combination)
        candidate_mols = [
            graph_utils.canonicalize(candidate)
            for candidate in structure_generator.create_single_bonds_map(
                carbon_type_list
            )
            if graph_utils.is_connected_graph(candidate)
        ]
        unique_candidates = np.unique(candidate_mols, axis=0)
        if len(unique_candidates) >= min_candidates:
            return [molecule.Molecule(candidate) for candidate in unique_candidates]
    raise RuntimeError("No representative molecule set found")


def benchmark_unique_mols(
    carbon_count: int,
    hydrogen_count: int,
    min_candidates: int,
    repeats: int,
) -> None:
    """Measure unique molecule filtering for a representative input set."""
    molecules = representative_molecules(carbon_count, hydrogen_count, min_candidates)
    for repeat_index in range(repeats):
        tracemalloc.start()
        start = time.perf_counter()
        unique = list(molecule_transformations.unique_mols(molecules))
        elapsed = time.perf_counter() - start
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        print(
            f"repeat={repeat_index + 1} input_molecules={len(molecules)} "
            f"unique_outputs={len(unique)} seconds={elapsed:.4f} "
            f"peak_bytes={peak}"
        )


def main() -> None:
    """Parse arguments and run benchmarks."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-carbon", type=int, default=5)
    parser.add_argument("--unique-carbon", type=int, default=6)
    parser.add_argument("--unique-hydrogen", type=int, default=8)
    parser.add_argument("--min-candidates", type=int, default=2)
    parser.add_argument("--unique-repeats", type=int, default=3)
    args = parser.parse_args()

    print("generation")
    benchmark_generation(args.max_carbon)
    print("unique_mols")
    benchmark_unique_mols(
        args.unique_carbon,
        args.unique_hydrogen,
        args.min_candidates,
        args.unique_repeats,
    )


if __name__ == "__main__":
    main()
