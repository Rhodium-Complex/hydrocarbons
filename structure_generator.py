"""Module for generating hydrocarbon structures from carbon valence patterns."""
from collections.abc import Generator

import numpy as np

import deduplication
import graph_utils
import molecule


def build_carbon_hydrogen_combination(
    c_num: int,
    h_num: int,
) -> Generator[np.ndarray, None, None]:
    """Generate possible carbon valence patterns for a formula."""

    def generate_combinations(
        c_num: int,
        h_num: int,
        current_combination: list,
    ) -> Generator[np.ndarray, None, None]:
        if c_num == 0 and h_num == 0:
            yield np.array(current_combination)
            return
        if c_num < 0 or h_num < 0:
            return

        current_hydrogen_count = 1
        if len(current_combination) > 0:
            current_hydrogen_count = current_combination[0]

        for current_hydrogen in range(current_hydrogen_count, 5):
            new_combination = [current_hydrogen] + current_combination
            yield from generate_combinations(
                c_num - 1,
                h_num - (4 - current_hydrogen),
                new_combination,
            )

    yield from generate_combinations(c_num, h_num, [])


def create_single_bonds_map(
    single_bond_degrees: np.ndarray,
) -> Generator[np.ndarray, None, None]:
    """Generate connected single-bond adjacency matrices for a degree pattern."""

    def generate_child_nodes(
        adjacency_matrix,
        remaining_degrees,
        equivalent_degree_flags,
        candidate_start_index=0,
    ) -> Generator[np.ndarray, None, None]:
        """Recursively add valid single bonds to the adjacency matrix."""
        atom_count = len(remaining_degrees)
        unprocessed = np.where(remaining_degrees > 0)[0]
        if len(unprocessed) == 0:
            if graph_utils.is_connected_graph(adjacency_matrix):
                yield adjacency_matrix
            return
        next_node_index = unprocessed[0]

        if next_node_index == atom_count:
            return

        for bond_candidate_index in range(candidate_start_index, atom_count):
            # Reject invalid or duplicate edge choices before copying state.
            if bond_candidate_index == next_node_index:
                continue
            if adjacency_matrix[bond_candidate_index][next_node_index] == 1:
                continue
            if remaining_degrees[bond_candidate_index] == 0:
                continue
            # Equal adjacent degree entries represent interchangeable labels at
            # this point in the search. Once the immediate next equivalent
            # choice has been considered, later equivalent choices would only
            # relabel the same partial matrix and can be skipped.
            if (
                equivalent_degree_flags[bond_candidate_index]
                and bond_candidate_index - next_node_index > 1
            ):
                continue

            # Copy only after all pruning checks pass; this path is hot.
            updated_adj_matrix = adjacency_matrix.copy()
            updated_bond_counts = remaining_degrees.copy()
            current_equivalent_flags = equivalent_degree_flags.copy()

            updated_adj_matrix[next_node_index][bond_candidate_index] = 1
            updated_adj_matrix[bond_candidate_index][next_node_index] = 1
            updated_bond_counts[next_node_index] -= 1
            updated_bond_counts[bond_candidate_index] -= 1
            current_equivalent_flags[bond_candidate_index + 1] = False

            yield from generate_child_nodes(
                updated_adj_matrix,
                updated_bond_counts,
                current_equivalent_flags,
                bond_candidate_index
                if updated_bond_counts[next_node_index]
                else 0,
            )

    atom_count = len(single_bond_degrees)
    equivalent_degree_flags = (
        [False]
        + [
            single_bond_degrees[i + 1] == single_bond_degrees[i]
            for i in range(atom_count - 1)
        ]
        + [True]
    )
    equivalent_degree_flags[1] = False

    yield from generate_child_nodes(
        np.full((atom_count, atom_count), 0, dtype="i4"),
        single_bond_degrees,
        equivalent_degree_flags,
    )


def build_structure(input_structure):
    """Build unique molecule structures from a carbon valence pattern."""
    # Generate and canonicalize connected candidates before deduplication.
    candidate_mols = [
        graph_utils.canonicalize(candidate)
        for candidate in create_single_bonds_map(input_structure)
    ]
    # Drop duplicate labeled matrices first, then remove structural isomorphs.
    unique_mols_list = deduplication.unique_mols(
        molecule.Molecule(candidate)
        for candidate in np.unique(candidate_mols, axis=0)
    )
    return [[molecule_entity] for molecule_entity in unique_mols_list]
