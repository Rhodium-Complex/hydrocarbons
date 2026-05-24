import unittest

import converter
import graph_utils
import molecule
import molecule_transformations
import numpy as np
import structure_generator


def formula_counts_for_carbon(carbon_count):
    current_carbon_structures = []
    counts = []
    for hydrogen_count in range(0, carbon_count * 2 + 3, 2)[::-1]:
        current_carbon_structures = [
            molecule_transformations.unique_dehydro_mols(structures)
            for structures in current_carbon_structures
        ]
        for combination in structure_generator.build_carbon_hydrogen_combination(carbon_count, hydrogen_count):
            current_carbon_structures += structure_generator.build_structure(combination)
        counts.append((hydrogen_count, sum(len(structures) for structures in current_carbon_structures)))
    return counts


class GenerationCorrectnessTests(unittest.TestCase):
    def test_formula_counts_c2_to_c3(self):
        expected_counts = {
            2: [(6, 1), (4, 1), (2, 1), (0, 0)],
            3: [(8, 1), (6, 2), (4, 3), (2, 2), (0, 1)],
        }
        for carbon_count, expected_count in expected_counts.items():
            with self.subTest(carbon_count=carbon_count):
                self.assertEqual(formula_counts_for_carbon(carbon_count), expected_count)

    def test_single_bond_generation_has_no_self_loops(self):
        for combination in ([2, 2], [2, 2, 2]):
            with self.subTest(combination=combination):
                for bonds in structure_generator.create_single_bonds_map(np.array(combination)):
                    self.assertTrue(np.all(np.diag(bonds) == 0))

    def test_unique_mols_yields_each_input_at_most_once(self):
        combination = np.array([2, 1, 1])
        candidate_mols = [
            graph_utils.canonicalize(candidate)
            for candidate in structure_generator.create_single_bonds_map(combination)
            if graph_utils.is_connected_graph(candidate)
        ]
        unique_candidates = np.unique(candidate_mols, axis=0)
        unique = list(molecule_transformations.unique_mols(molecule.Molecule(candidate) for candidate in unique_candidates))

        self.assertEqual(len(unique_candidates), 1)
        self.assertEqual(len(unique), 1)

    def test_ring_smiles_conversion_terminates(self):
        cyclopropane = molecule.Molecule(
            np.array(
                [
                    [0, 1, 1],
                    [1, 0, 1],
                    [1, 1, 0],
                ]
            )
        )

        self.assertEqual(converter.mat2smiles(cyclopropane), "C1CC1")


if __name__ == "__main__":
    unittest.main()
