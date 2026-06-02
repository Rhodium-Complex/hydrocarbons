""" test cases to verify the correctness of the hydrocarbon generation pipeline. """
import unittest
from unittest import mock

import converter
import generation_pipeline
import graph_utils
import molecule
import molecule_transformations
import numpy as np
import structure_generator


def formula_counts_for_carbon(carbon_count):
    """Helper function to compute the counts of structures for a given carbon count."""
    current_carbon_structures = []
    counts = []
    for hydrogen_count in range(0, carbon_count * 2 + 3, 2)[::-1]:
        current_carbon_structures = [
            molecule_transformations.unique_dehydro_mols(structures)
            for structures in current_carbon_structures
        ]
        for combination in structure_generator.build_carbon_hydrogen_combination(
            carbon_count,
            hydrogen_count,
        ):
            generated_structures = structure_generator.build_structure(combination)
            current_carbon_structures += generated_structures
        counts.append(
            (
                hydrogen_count,
                sum(len(structures) for structures in current_carbon_structures),
            )
        )
    return counts


class GenerationCorrectnessTests(unittest.TestCase):
    """Unit tests to verify the correctness of the hydrocarbon generation pipeline."""
    def test_formula_counts_c2_to_c3(self):
        """Test that the counts of generated structures for C2 and C3 match expected values."""
        expected_counts = {
            2: [(6, 1), (4, 1), (2, 1), (0, 0)],
            3: [(8, 1), (6, 2), (4, 3), (2, 2), (0, 1)],
        }
        for carbon_count, expected_count in expected_counts.items():
            with self.subTest(carbon_count=carbon_count):
                self.assertEqual(
                    formula_counts_for_carbon(carbon_count),
                    expected_count,
                )

    def test_single_bond_generation_has_no_self_loops(self):
        """Test that the generated single-bond adjacency matrices do not contain self-loops."""
        for combination in ([2, 2], [2, 2, 2]):
            with self.subTest(combination=combination):
                for bonds in structure_generator.create_single_bonds_map(
                    np.array(combination)
                ):
                    self.assertTrue(np.all(np.diag(bonds) == 0))

    def test_unique_mols_yields_each_input_at_most_once(self):
        """Test that the unique_mols function returns each input molecule at most once."""
        combination = np.array([2, 1, 1])
        candidate_mols = [
            graph_utils.canonicalize(candidate)
            for candidate in structure_generator.create_single_bonds_map(combination)
            if graph_utils.is_connected_graph(candidate)
        ]
        unique_candidates = np.unique(candidate_mols, axis=0)
        unique = list(
            molecule_transformations.unique_mols(
                molecule.Molecule(candidate)
                for candidate in unique_candidates
            )
        )

        self.assertEqual(len(unique_candidates), 1)
        self.assertEqual(len(unique), 1)

    def test_ring_smiles_conversion_terminates(self):
        """Test that the SMILES conversion for a high-symmetry ring terminates."""
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

    def test_high_symmetry_tetravalent_c9_generation(self):
        """Test that the generation of structures with 9 carbons and 4 hydrogens terminates and 
        produces the expected count."""
        structures = structure_generator.build_structure(np.full(9, 4))

        self.assertEqual(sum(len(group) for group in structures), 16)

    def test_pipeline_step_counts_c2_to_c3_without_smiles(self):
        """Test that the generation pipeline produces the expected counts for C2 and C3 
        without generating SMILES strings."""
        observed_counts = []

        generation_pipeline.run_generation(
            min_carbon=2,
            max_carbon=3,
            workers=2,
            include_smiles=False,
            log_step=lambda result: observed_counts.append(
                (result.carbon_count, result.hydrogen_count, result.count)
            ),
        )

        self.assertEqual(
            observed_counts,
            [
                (2, 6, 1),
                (2, 4, 1),
                (2, 2, 1),
                (2, 0, 0),
                (3, 8, 1),
                (3, 6, 2),
                (3, 4, 3),
                (3, 2, 2),
                (3, 0, 1),
            ],
        )

    def test_pipeline_workers_one_uses_sequential_map(self):
        """Test that the generation pipeline does not use ProcessPoolExecutor 
        when workers is set to 1."""
        observed_counts = []

        with mock.patch("generation_pipeline.ProcessPoolExecutor") as executor:
            generation_pipeline.run_generation(
                min_carbon=2,
                max_carbon=2,
                workers=1,
                include_smiles=False,
                log_step=lambda result: observed_counts.append(
                    (result.carbon_count, result.hydrogen_count, result.count)
                ),
            )

        executor.assert_not_called()
        self.assertEqual(
            observed_counts,
            [
                (2, 6, 1),
                (2, 4, 1),
                (2, 2, 1),
                (2, 0, 0),
            ],
        )


if __name__ == "__main__":
    unittest.main()
