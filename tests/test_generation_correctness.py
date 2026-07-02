""" test cases to verify the correctness of the hydrocarbon generation pipeline. """
import importlib.util
import itertools
import unittest
from unittest import mock

import converter
import deduplication
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
        for combination in ([2, 2], [2, 2, 2], [3, 3, 2, 2, 1, 1, 1, 1]):
            with self.subTest(combination=combination):
                for bonds in structure_generator.create_single_bonds_map(
                    np.array(combination)
                ):
                    self.assertTrue(np.all(np.diag(bonds) == 0))

    def test_single_bond_generation_yields_connected_degree_matches(self):
        """Test that single-bond candidates are connected and match requested degrees."""
        for combination in ([2, 1, 1], [3, 3, 2, 2, 1, 1, 1, 1]):
            with self.subTest(combination=combination):
                requested_degrees = np.array(combination)
                for bonds in structure_generator.create_single_bonds_map(
                    requested_degrees
                ):
                    self.assertTrue(graph_utils.is_connected_graph(bonds))
                    np.testing.assert_array_equal(bonds.sum(axis=0), requested_degrees)

    def test_single_bond_generation_omits_disconnected_degree_patterns(self):
        """Test that degree patterns with no connected realization yield no candidates."""
        self.assertEqual(
            list(structure_generator.create_single_bonds_map(np.array([1, 1, 1, 1]))),
            [],
        )

    def test_unique_mols_yields_each_input_at_most_once(self):
        """Test that the unique_mols function returns each input molecule at most once."""
        combination = np.array([2, 1, 1])
        candidate_mols = [
            graph_utils.canonicalize(candidate)
            for candidate in structure_generator.create_single_bonds_map(combination)
        ]
        unique_candidates = np.unique(candidate_mols, axis=0)
        unique = list(
            deduplication.unique_mols(
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

        self.assertEqual(converter.mat2smiles_variants(cyclopropane), ["C1CC1"])

    def test_high_symmetry_tetravalent_c9_generation(self):
        """Test that the generation of structures with 9 carbons and 4 hydrogens terminates and 
        produces the expected count."""
        structures = structure_generator.build_structure(np.full(9, 4))

        self.assertEqual(sum(len(group) for group in structures), 16)

    def test_canonicalize_preserves_c11_h24_candidate_size(self):
        """Test that canonicalization does not drop vertices for C11H24 candidates."""
        for combination in structure_generator.build_carbon_hydrogen_combination(11, 24):
            for candidate in structure_generator.create_single_bonds_map(combination):
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

    def test_simple_carbon_hydrogen_ligands_skip_graph_refinement(self):
        """Test that locally distinct alkene ligands do not invoke WL refinement."""
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

        with mock.patch.object(
            stereochemistry,
            "_refined_atom_colors",
            wraps=stereochemistry._refined_atom_colors,
        ) as refined:
            analysis = stereochemistry.analyze_ez(two_butene)

        self.assertEqual(len(analysis.double_bonds), 1)
        refined.assert_not_called()

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

        self.assertEqual(smiles, ["C/C=C/C(C)C", "C/C=C\\C(C)C"])
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
        smiles = [
            converter.mat2stereo_smiles(octadiene, assignment, analysis)
            for assignment in analysis.assignments
        ]
        self.assertEqual(
            smiles,
            [
                "C/C=C/C/C=C/CC",
                "C/C=C/C/C=C\\CC",
                "C/C=C\\C/C=C/CC",
                "C/C=C\\C/C=C\\CC",
            ],
        )
        self.assertTrue(all("[E:" not in text and "[Z:" not in text for text in smiles))

    def test_branched_diene_keeps_standard_slash_stereo(self):
        """Test that branching does not force multiple E/Z labels into the suffix."""
        bonds = np.zeros((9, 9), dtype=int)
        for atom1, atom2, order in (
            (0, 1, 1),
            (1, 2, 2),
            (2, 3, 1),
            (3, 4, 1),
            (4, 5, 2),
            (5, 6, 1),
            (6, 7, 1),
            (3, 8, 1),
        ):
            bonds[atom1][atom2] = order
            bonds[atom2][atom1] = order
        branched_diene = molecule.Molecule(bonds)

        analysis = stereochemistry.analyze_ez(branched_diene)
        smiles = [
            converter.mat2stereo_smiles(branched_diene, assignment, analysis)
            for assignment in analysis.assignments
        ]

        self.assertEqual(len(analysis.double_bonds), 2)
        self.assertEqual(len(smiles), 4)
        self.assertTrue(all("/" in text or "\\" in text for text in smiles))
        self.assertTrue(all("[E:" not in text and "[Z:" not in text for text in smiles))

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

    def test_ring_and_exocyclic_ez_are_both_embedded(self):
        """Test that ring closures and external alkenes both carry slash stereo."""
        molecule_obj = molecule.Molecule(
            np.array(
                [
                    [0, 1, 0, 0, 0, 0],
                    [1, 0, 2, 0, 0, 0],
                    [0, 2, 0, 1, 0, 0],
                    [0, 0, 1, 0, 1, 1],
                    [0, 0, 0, 1, 0, 2],
                    [0, 0, 0, 1, 2, 0],
                ]
            )
        )
        analysis = stereochemistry.analyze_ez(molecule_obj)
        assignment = next(
            item
            for item in analysis.assignments
            if item.labels == ((1, 2, "Z"), (4, 5, "Z"))
        )

        smiles = converter.mat2stereo_smiles(molecule_obj, assignment, analysis)

        self.assertEqual(smiles, "C/C=C\\C/1/C=C1")
        self.assertNotIn("[Z:", smiles)

    def test_conflicting_ring_constraints_fall_back_to_non_stereo_smiles(self):
        """Test that inconsistent formal ring assignments remain visible without stereo."""
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

        variants = converter.mat2smiles_variants(cyclohexatriene, include_stereo=True)

        self.assertEqual(len(variants), 4)
        non_stereo = converter.mat2smiles_variants(cyclohexatriene)[0]
        self.assertEqual(variants[1], non_stereo)
        self.assertEqual(variants[3], non_stereo)
        self.assertTrue(all("[E:" not in text and "[Z:" not in text for text in variants))

    def test_ring_stereo_survives_atom_renumbering(self):
        """Test that ring slash constraints remain representable after renumbering."""
        bonds = np.array(
            [
                [0, 1, 0, 0, 0, 0],
                [1, 0, 2, 0, 0, 0],
                [0, 2, 0, 1, 0, 0],
                [0, 0, 1, 0, 1, 1],
                [0, 0, 0, 1, 0, 2],
                [0, 0, 0, 1, 2, 0],
            ]
        )
        permutation = [3, 5, 1, 4, 0, 2]
        original = molecule.Molecule(bonds)
        renumbered = molecule.Molecule(bonds[np.ix_(permutation, permutation)])

        original_variants = converter.mat2smiles_variants(original, include_stereo=True)
        renumbered_variants = converter.mat2smiles_variants(renumbered, include_stereo=True)

        self.assertEqual(len(original_variants), len(renumbered_variants))
        self.assertTrue(all(original_variants))
        self.assertTrue(all(renumbered_variants))

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

    def test_single_tetrahedral_center_expands_to_two_smiles(self):
        bonds = np.zeros((7, 7), dtype=int)
        for atom1, atom2 in ((0, 1), (0, 2), (2, 3), (0, 4), (4, 5), (5, 6)):
            bonds[atom1][atom2] = bonds[atom2][atom1] = 1
        molecule_obj = molecule.Molecule(bonds)

        analysis = stereochemistry.analyze_chiral(molecule_obj)
        smiles = converter.mat2smiles_variants(
            molecule_obj,
            include_tetrahedral_stereo=True,
        )

        self.assertEqual(len(analysis.tetrahedral_centers), 1)
        self.assertEqual(len(smiles), 2)
        self.assertTrue(any("[C@H]" in value for value in smiles))
        self.assertTrue(any("[C@@H]" in value for value in smiles))

    def test_equal_tetrahedral_ligands_remain_inactive(self):
        bonds = np.zeros((5, 5), dtype=int)
        for atom1, atom2 in ((0, 1), (0, 2), (0, 3), (3, 4)):
            bonds[atom1][atom2] = bonds[atom2][atom1] = 1
        molecule_obj = molecule.Molecule(bonds)

        analysis = stereochemistry.analyze_chiral(molecule_obj)
        smiles = converter.mat2smiles_variants(molecule_obj, False, True)

        self.assertEqual(len(analysis.tetrahedral_centers), 1)
        self.assertEqual(analysis.assignments[0].active_centers, frozenset())
        self.assertEqual(len(smiles), 1)
        self.assertNotIn("@", smiles[0])

    def test_odd_stabilizer_makes_tetrahedral_center_inactive(self):
        bonds = np.zeros((5, 5), dtype=int)
        for atom1, atom2 in ((0, 1), (0, 2), (0, 3), (3, 4)):
            bonds[atom1][atom2] = bonds[atom2][atom1] = 1
        analysis = stereochemistry.analyze_chiral(molecule.Molecule(bonds))
        assignment = analysis.assignments[0]

        active = stereochemistry.active_chiral_centers(
            bonds,
            analysis,
            assignment,
        )

        self.assertGreater(len(analysis.automorphisms), 1)
        self.assertNotIn(("T", analysis.tetrahedral_centers[0].atom), active)

    def test_ez_labels_filter_automorphisms_for_chiral_enumeration(self):
        bonds = np.zeros((8, 8), dtype=int)
        centers = (
            stereochemistry.TetrahedralCenter(0, (-1, 2, 3, 4), 1),
            stereochemistry.TetrahedralCenter(1, (-1, 5, 6, 7), 1),
        )
        identity = tuple(range(8))
        swap_halves = (1, 0, 5, 6, 7, 2, 3, 4)
        analysis = stereochemistry.ChiralAnalysis(
            tetrahedral_centers=centers,
            allene_centers=(),
            assignments=(),
            automorphisms=(identity, swap_halves),
        )
        same_labels = stereochemistry.EzAssignment(((2, 3, "E"), (5, 6, "E")))
        different_labels = stereochemistry.EzAssignment(((2, 3, "E"), (5, 6, "Z")))

        symmetric = stereochemistry.chiral_assignments_for_ez(
            bonds, analysis, same_labels
        )
        symmetry_broken = stereochemistry.chiral_assignments_for_ez(
            bonds, analysis, different_labels
        )

        self.assertEqual(len(symmetric), 3)
        self.assertEqual(len(symmetry_broken), 4)

    def test_chiral_analysis_is_invariant_under_atom_renumbering(self):
        bonds = np.zeros((8, 8), dtype=int)
        for atom1, atom2 in (
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (2, 6), (3, 7)
        ):
            bonds[atom1][atom2] = bonds[atom2][atom1] = 1
        permutation = [6, 3, 0, 7, 2, 5, 1, 4]

        original = stereochemistry.analyze_chiral(molecule.Molecule(bonds))
        renumbered = stereochemistry.analyze_chiral(
            molecule.Molecule(bonds[np.ix_(permutation, permutation)])
        )

        self.assertEqual(len(original.automorphisms), len(renumbered.automorphisms))
        self.assertEqual(len(original.assignments), len(renumbered.assignments))
        self.assertEqual(
            sorted(len(item.active_centers) for item in original.assignments),
            sorted(len(item.active_centers) for item in renumbered.assignments),
        )

    def test_configuration_dependent_center_adds_fourth_isomer(self):
        bonds = np.array(
            [
                [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                [1, 0, 1, 0, 1, 0, 0, 0, 0, 0],
                [0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 1, 1, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 1, 1, 0],
                [0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 1, 0, 0, 1],
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
            ]
        )
        molecule_obj = molecule.Molecule(bonds)

        analysis = stereochemistry.analyze_chiral(molecule_obj)
        smiles = converter.mat2smiles_variants(molecule_obj, False, True)

        self.assertEqual(len(analysis.tetrahedral_centers), 3)
        self.assertEqual(len(smiles), 4)
        self.assertEqual(sum(value.count("[C@") == 3 for value in smiles), 2)
        self.assertEqual(sum(value.count("[C@") == 2 for value in smiles), 2)

    def test_alkane_stereoisomer_counts_match_oeis_a000628(self):
        expected = (1, 1, 1, 1, 2, 3, 5, 11, 24, 55, 136, 345, 900)
        actual = [1]
        for carbon_count in range(1, 13):
            if carbon_count == 1:
                count = 1
            else:
                groups = []
                for combination in structure_generator.build_carbon_hydrogen_combination(
                    carbon_count,
                    2 * carbon_count + 2,
                ):
                    groups += structure_generator.build_structure(combination)
                count = sum(
                    len(converter.mat2smiles_variants(item, False, True))
                    for item in itertools.chain.from_iterable(groups)
                )
            actual.append(count)

        self.assertEqual(tuple(actual), expected)

    def test_symmetric_two_center_molecule_has_three_assignments(self):
        bonds = np.zeros((8, 8), dtype=int)
        for atom1, atom2 in (
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (2, 6), (3, 7)
        ):
            bonds[atom1][atom2] = bonds[atom2][atom1] = 1
        molecule_obj = molecule.Molecule(bonds)

        analysis = stereochemistry.analyze_chiral(molecule_obj)
        smiles = converter.mat2smiles_variants(
            molecule_obj,
            include_tetrahedral_stereo=True,
        )

        self.assertEqual(len(analysis.tetrahedral_centers), 2)
        self.assertEqual(len(analysis.assignments), 3)
        self.assertEqual(len(smiles), 3)

    def test_allene_like_center_expands_to_two_smiles(self):
        bonds = np.zeros((5, 5), dtype=int)
        for atom1, atom2, order in ((0, 1, 2), (1, 2, 2), (0, 3, 1), (2, 4, 1)):
            bonds[atom1][atom2] = bonds[atom2][atom1] = order
        molecule_obj = molecule.Molecule(bonds)

        analysis = stereochemistry.analyze_chiral(molecule_obj)
        smiles = converter.mat2smiles_variants(
            molecule_obj,
            include_tetrahedral_stereo=True,
        )

        self.assertEqual(len(analysis.allene_centers), 1)
        self.assertEqual(len(smiles), 2)
        self.assertTrue(any("@AL1" in value for value in smiles))
        self.assertTrue(any("@AL2" in value for value in smiles))

    def test_odd_cumulene_is_not_allene_like_stereogenic(self):
        bonds = np.zeros((6, 6), dtype=int)
        for atom1, atom2, order in (
            (0, 1, 2), (1, 2, 2), (2, 3, 2), (0, 4, 1), (3, 5, 1)
        ):
            bonds[atom1][atom2] = bonds[atom2][atom1] = order
        molecule_obj = molecule.Molecule(bonds)

        self.assertEqual(stereochemistry.analyze_chiral(molecule_obj).allene_centers, ())

    def test_extended_even_cumulene_is_allene_like(self):
        bonds = np.zeros((7, 7), dtype=int)
        for atom1, atom2, order in (
            (0, 1, 2), (1, 2, 2), (2, 3, 2), (3, 4, 2), (0, 5, 1), (4, 6, 1)
        ):
            bonds[atom1][atom2] = bonds[atom2][atom1] = order
        molecule_obj = molecule.Molecule(bonds)

        analysis = stereochemistry.analyze_chiral(molecule_obj)
        smiles = converter.mat2smiles_variants(
            molecule_obj,
            include_tetrahedral_stereo=True,
        )

        self.assertEqual(len(analysis.allene_centers), 1)
        self.assertEqual(len(smiles), 2)

    def test_ez_and_tetrahedral_configurations_form_product(self):
        bonds = np.zeros((8, 8), dtype=int)
        for atom1, atom2, order in (
            (0, 1, 1), (0, 2, 1), (2, 3, 1), (0, 4, 1),
            (4, 5, 1), (5, 6, 2), (6, 7, 1),
        ):
            bonds[atom1][atom2] = bonds[atom2][atom1] = order
        molecule_obj = molecule.Molecule(bonds)

        smiles = converter.mat2smiles_variants(
            molecule_obj,
            include_stereo=True,
            include_tetrahedral_stereo=True,
        )

        self.assertEqual(len(smiles), 4)
        self.assertTrue(all("@" in value for value in smiles))
        self.assertTrue(all("/" in value or "\\" in value for value in smiles))

    def test_tetrahedral_outputs_survive_atom_renumbering(self):
        bonds = np.zeros((7, 7), dtype=int)
        for atom1, atom2 in ((0, 1), (0, 2), (2, 3), (0, 4), (4, 5), (5, 6)):
            bonds[atom1][atom2] = bonds[atom2][atom1] = 1
        permutation = [4, 2, 6, 0, 5, 1, 3]
        original = molecule.Molecule(bonds)
        renumbered = molecule.Molecule(bonds[np.ix_(permutation, permutation)])

        original_smiles = converter.mat2smiles_variants(original, False, True)
        renumbered_smiles = converter.mat2smiles_variants(renumbered, False, True)

        self.assertEqual(len(original_smiles), 2)
        self.assertEqual(len(renumbered_smiles), 2)

    def test_ring_tetrahedral_centers_render_valid_ring_smiles(self):
        bonds = np.zeros((8, 8), dtype=int)
        for atom1, atom2 in (
            (0, 1), (1, 2), (2, 3), (3, 4),
            (4, 5), (5, 0), (0, 6), (1, 7),
        ):
            bonds[atom1][atom2] = bonds[atom2][atom1] = 1
        molecule_obj = molecule.Molecule(bonds)

        analysis = stereochemistry.analyze_chiral(molecule_obj)
        smiles = converter.mat2smiles_variants(molecule_obj, False, True)

        self.assertEqual(len(analysis.tetrahedral_centers), 2)
        self.assertEqual(len(smiles), 3)
        self.assertTrue(all("@" in value and "1" in value for value in smiles))

    @unittest.skipUnless(importlib.util.find_spec("rdkit"), "RDKit is not installed")
    def test_atom_centered_stereo_smiles_parse_with_rdkit(self):
        from rdkit import Chem

        tetra_bonds = np.zeros((7, 7), dtype=int)
        for atom1, atom2 in ((0, 1), (0, 2), (2, 3), (0, 4), (4, 5), (5, 6)):
            tetra_bonds[atom1][atom2] = tetra_bonds[atom2][atom1] = 1
        allene_bonds = np.zeros((5, 5), dtype=int)
        for atom1, atom2, order in ((0, 1, 2), (1, 2, 2), (0, 3, 1), (2, 4, 1)):
            allene_bonds[atom1][atom2] = allene_bonds[atom2][atom1] = order

        tetra_smiles = converter.mat2smiles_variants(
            molecule.Molecule(tetra_bonds), False, True
        )
        allene_smiles = converter.mat2smiles_variants(
            molecule.Molecule(allene_bonds), False, True
        )

        tetra_molecules = [Chem.MolFromSmiles(value) for value in tetra_smiles]
        self.assertTrue(all(item is not None for item in tetra_molecules))
        self.assertEqual(
            len({
                Chem.MolToSmiles(item, isomericSmiles=True)
                for item in tetra_molecules
                if item is not None
            }),
            2,
        )
        self.assertTrue(all(Chem.MolFromSmiles(value) is not None for value in allene_smiles))

    def test_tetrahedral_flag_is_independent_from_ez_flag(self):
        groups = generation_pipeline.run_generation_smiles_groups(
            min_carbon=7,
            max_carbon=7,
            workers=1,
            include_methane=False,
            include_tetrahedral_stereo=True,
        )

        values = [smiles for group in groups for smiles in group.smiles]
        self.assertTrue(any("@" in smiles for smiles in values))

    def test_pipeline_workers_one_uses_sequential_map(self):
        """Test that the generation pipeline does not use ProcessPoolExecutor 
        when workers is set to 1."""
        observed_counts = []

        with mock.patch("generation_pipeline.ProcessPoolExecutor") as executor:
            generation_pipeline.run_generation(
                min_carbon=2,
                max_carbon=2,
                workers=1,
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
