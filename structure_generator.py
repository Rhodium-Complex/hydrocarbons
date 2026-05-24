"""Module for generating hydrocarbon structures from carbon valence patterns."""
from collections.abc import Generator

import numpy as np

import graph_utils
import molecule
import molecule_transformations


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
    carbon_type_list: list | np.ndarray,
) -> Generator[np.ndarray, None, None]:
    """Generate single-bond adjacency matrices for a carbon valence pattern."""

    def generate_child_nodes(
        adjacency_matrix,
        bond_degrees,
        visitable_nodes,
        start_index=0,
    ) -> Generator[np.ndarray, None, None]:
        """Recursively add valid single bonds to the adjacency matrix."""
        unprocessed = np.where(bond_degrees > 0)[0]
        if len(unprocessed) == 0:
            yield adjacency_matrix
            return
        next_node_index = unprocessed[0]

        if next_node_index == len(bond_degrees):
            return

        for bond_candidate_index in range(start_index, len(bond_degrees)):
            if bond_candidate_index == next_node_index:
                continue
            if adjacency_matrix[bond_candidate_index][next_node_index] == 1:
                continue
            if bond_degrees[bond_candidate_index] == 0:
                continue
            if (
                visitable_nodes[bond_candidate_index]
                and bond_candidate_index - next_node_index > 1
            ):
                continue

            updated_adj_matrix = adjacency_matrix.copy()
            updated_bond_counts = bond_degrees.copy()
            current_visit_flags = visitable_nodes.copy()

            updated_adj_matrix[next_node_index][bond_candidate_index] = 1
            updated_adj_matrix[bond_candidate_index][next_node_index] = 1
            updated_bond_counts[next_node_index] -= 1
            updated_bond_counts[bond_candidate_index] -= 1
            current_visit_flags[bond_candidate_index + 1] = False

            yield from generate_child_nodes(
                updated_adj_matrix,
                updated_bond_counts,
                current_visit_flags,
                bond_candidate_index
                if updated_bond_counts[next_node_index]
                else 0,
            )

    num_carbon_types = len(carbon_type_list)
    visited = (
        [False]
        + [
            carbon_type_list[i + 1] == carbon_type_list[i]
            for i in range(num_carbon_types - 1)
        ]
        + [True]
    )
    visited[1] = False

    yield from generate_child_nodes(
        np.full((num_carbon_types, num_carbon_types), 0, dtype="i4"),
        carbon_type_list,
        visited,
    )


def build_structure(input_structure):
    """Build unique molecule structures from a carbon valence pattern."""
    candidate_mols = [
        graph_utils.canonicalize(candidate)
        for candidate in create_single_bonds_map(input_structure)
        if graph_utils.is_connected_graph(candidate)
    ]
    unique_mols_list = molecule_transformations.unique_mols(
        molecule.Molecule(candidate)
        for candidate in np.unique(candidate_mols, axis=0)
    )
    return [[molecule_entity] for molecule_entity in unique_mols_list]
