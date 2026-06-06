""" test cases to verify the correctness of the hydrocarbon generation pipeline. """
import unittest
from unittest import mock

import converter
import generation_pipeline
import graph_utils
import isomorphism
import molecule
import molecule_transformations
import numpy as np
import stereochemistry
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

    def test_isomorphism_automorphisms_for_symmetric_chain(self):
        """Test automorphism enumeration for a symmetric carbon chain."""
        butane_bonds = np.array(
            [
                [0, 1, 0, 0],
                [1, 0, 1, 0],
                [0, 1, 0, 1],
                [0, 0, 1, 0],
            ]
        )

        self.assertEqual(
            set(isomorphism.automorphisms(butane_bonds)),
            {
                (0, 1, 2, 3),
                (3, 2, 1, 0),
            },
        )

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

    def test_canonicalize_preserves_c11_h24_candidate_size(self):
        """Test that canonicalization does not drop vertices for C11H24 candidates."""
        for combination in structure_generator.build_carbon_hydrogen_combination(11, 24):
            for candidate in structure_generator.create_single_bonds_map(combination):
                if not graph_utils.is_connected_graph(candidate):
                    continue

                canonical_candidate = graph_utils.canonicalize(candidate)

                self.assertEqual(canonical_candidate.shape, candidate.shape)

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

    def test_pipeline_smiles_output_keeps_formula_separators(self):
        """Test that flat SMILES output keeps downstream formula group markers."""
        smiles_results = generation_pipeline.run_generation(
            min_carbon=2,
            max_carbon=2,
            workers=1,
            include_smiles=True,
        )

        self.assertEqual(smiles_results[:4], ["N#N", "N#N", "N#N", "C"])
        self.assertEqual(smiles_results[4], "N#N")
        self.assertIn("CC", smiles_results)

    def test_structured_smiles_groups_c2(self):
        """Test grouped SMILES output for B5 PDF export."""
        groups = generation_pipeline.run_generation_smiles_groups(
            min_carbon=2,
            max_carbon=2,
            workers=1,
        )

        self.assertEqual(
            [group.label for group in groups],
            ["CH4", "C2H6", "C2H4", "C2H2", "C2H0"],
        )
        smiles_by_label = {group.label: group.smiles for group in groups}
        self.assertEqual(smiles_by_label["CH4"], ["C"])
        self.assertIn("CC", smiles_by_label["C2H6"])
        self.assertIn("C=C", smiles_by_label["C2H4"])
        self.assertIn("C#C", smiles_by_label["C2H2"])
        self.assertEqual(smiles_by_label["C2H0"], [])

    def test_ethene_and_propene_have_no_ez_assignments(self):
        """Test that terminal/simple alkenes without two distinct ligand pairs are not E/Z."""
        ethene = molecule.Molecule(
            np.array(
                [
                    [0, 2],
                    [2, 0],
                ]
            )
        )
        propene = molecule.Molecule(
            np.array(
                [
                    [0, 2, 0],
                    [2, 0, 1],
                    [0, 1, 0],
                ]
            )
        )

        ethene_analysis = stereochemistry.analyze_ez(ethene)
        propene_analysis = stereochemistry.analyze_ez(propene)
        self.assertEqual(ethene_analysis.double_bonds, ())
        self.assertEqual(propene_analysis.double_bonds, ())
        self.assertEqual(
            ethene_analysis.assignments,
            (stereochemistry.EzAssignment(labels=()),),
        )

    def test_two_butene_expands_to_e_and_z_smiles(self):
        """Test that 2-butene is expanded into separate E/Z outputs."""
        two_butene = molecule.Molecule(
            np.array(
                [
                    [0, 1, 0, 0],
                    [1, 0, 2, 0],
                    [0, 2, 0, 1],
                    [0, 0, 1, 0],
                ]
            )
        )

        analysis = stereochemistry.analyze_ez(two_butene)
        smiles = [
            converter.mat2stereo_smiles(two_butene, assignment, analysis)
            for assignment in analysis.assignments
        ]

        self.assertEqual(len(analysis.assignments), 2)
        self.assertEqual(smiles, ["C/C=C/C", "C/C=C\\C"])

    def test_branched_alkene_expands_to_slash_stereo_smiles(self):
        """Test that branched E/Z alkenes produce drawable slash SMILES."""
        methyl_pentene = molecule.Molecule(
            np.array(
                [
                    [0, 1, 0, 0, 0, 0],
                    [1, 0, 2, 0, 0, 0],
                    [0, 2, 0, 1, 0, 0],
                    [0, 0, 1, 0, 1, 1],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                ]
            )
        )

        analysis = stereochemistry.analyze_ez(methyl_pentene)
        smiles = [
            converter.mat2stereo_smiles(methyl_pentene, assignment, analysis)
            for assignment in analysis.assignments
        ]

        self.assertEqual(smiles, ["C/C=C/C(C)(C)", "C/C=C\\C(C)(C)"])
        self.assertTrue(all("[E:" not in text and "[Z:" not in text for text in smiles))

    def test_equal_ligands_do_not_create_ez_double_bond(self):
        """Test that equal substituents on one alkene carbon prevent E/Z assignment."""
        equal_methyls = molecule.Molecule(
            np.array(
                [
                    [0, 2, 1, 1, 0],
                    [2, 0, 0, 0, 1],
                    [1, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0],
                ]
            )
        )

        self.assertEqual(stereochemistry.analyze_ez(equal_methyls).double_bonds, ())

    def test_multiple_double_bonds_expand_assignments(self):
        """Test that two independent E/Z double bonds produce four assignments."""
        octadiene = molecule.Molecule(
            np.array(
                [
                    [0, 1, 0, 0, 0, 0, 0, 0],
                    [1, 0, 2, 0, 0, 0, 0, 0],
                    [0, 2, 0, 1, 0, 0, 0, 0],
                    [0, 0, 1, 0, 1, 0, 0, 0],
                    [0, 0, 0, 1, 0, 2, 0, 0],
                    [0, 0, 0, 0, 2, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1, 0, 1],
                    [0, 0, 0, 0, 0, 0, 1, 0],
                ]
            )
        )

        analysis = stereochemistry.analyze_ez(octadiene)
        self.assertEqual(len(analysis.double_bonds), 2)
        self.assertEqual(len(analysis.assignments), 4)

    def test_cyclohexatriene_double_bonds_are_formal_ez_candidates(self):
        """Test that ring double bonds are treated as formal E/Z candidates."""
        cyclohexatriene = molecule.Molecule(
            np.array(
                [
                    [0, 2, 0, 0, 0, 1],
                    [2, 0, 1, 0, 0, 0],
                    [0, 1, 0, 2, 0, 0],
                    [0, 0, 2, 0, 1, 0],
                    [0, 0, 0, 1, 0, 2],
                    [1, 0, 0, 0, 2, 0],
                ]
            )
        )

        self.assertEqual(len(stereochemistry.analyze_ez(cyclohexatriene).double_bonds), 3)

    def test_generation_stereo_flag_preserves_default_outputs(self):
        """Test that stereo output is opt-in and default generation remains unchanged."""
        default_groups = generation_pipeline.run_generation_smiles_groups(
            min_carbon=4,
            max_carbon=4,
            workers=1,
            include_methane=False,
        )
        stereo_groups = generation_pipeline.run_generation_smiles_groups(
            min_carbon=4,
            max_carbon=4,
            workers=1,
            include_methane=False,
            include_stereo=True,
        )

        default_c4h8 = {
            group.label: group.smiles for group in default_groups
        }["C4H8"]
        stereo_c4h8 = {
            group.label: group.smiles for group in stereo_groups
        }["C4H8"]
        self.assertEqual(len(default_c4h8), 5)
        self.assertGreater(len(stereo_c4h8), len(default_c4h8))
        self.assertIn("C/C=C/C", stereo_c4h8)

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
