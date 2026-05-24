import itertools
import unittest

import graph_utils
import molecule
import molecule_transformations
import numpy as np
import structure_generator


def reference_unique_mols(molecule_matrix):
    unique_molecule_records = {}
    for molecule_obj in molecule_matrix:
        molecule_fingerprint = molecule_obj.fingerprint
        bond_fingerprints = graph_utils.morgan(molecule_obj.bonds)
        bond_fingerprints = [list(np.where(bond_fingerprints == i)[0]) for i in np.unique(bond_fingerprints)[::-1]]

        for bond_fingerprint_group in bond_fingerprints:
            for permutation in itertools.permutations(bond_fingerprint_group):
                canonical_permutation = tuple(permutation)
                if molecule_fingerprint not in unique_molecule_records:
                    unique_molecule_records[molecule_fingerprint] = [
                        molecule_obj.bonds[:, canonical_permutation][canonical_permutation, :]
                    ]
                    yield molecule_obj
                    continue

                for i in itertools.product(*[itertools.permutations(j) for j in bond_fingerprints]):
                    i = sum(i, ())
                    current_bond_structure = molecule_obj.bonds[:, i][i, :]
                    if any(np.array_equal(current_bond_structure, j) for j in unique_molecule_records[molecule_fingerprint]):
                        break
                else:
                    unique_molecule_records[molecule_fingerprint].append(
                        molecule_obj.bonds[:, canonical_permutation][canonical_permutation, :]
                    )
                    yield molecule_obj


def build_count_for_carbon(carbon_count):
    total_structures = 0
    for hydrogen_count in range(0, carbon_count * 2 + 3, 2)[::-1]:
        for combination in structure_generator.build_carbon_hydrogen_combination(carbon_count, hydrogen_count):
            total_structures += len(structure_generator.build_structure(combination))
    return total_structures


def representative_molecules():
    carbon_count = 5
    hydrogen_count = 6
    for combination in structure_generator.build_carbon_hydrogen_combination(carbon_count, hydrogen_count):
        candidate_mols = [
            graph_utils.canonicalize(candidate)
            for candidate in structure_generator.create_single_bonds_map(combination)
            if graph_utils.is_connected_graph(candidate)
        ]
        unique_candidates = np.unique(candidate_mols, axis=0)
        if len(unique_candidates) >= 5:
            return [molecule.Molecule(candidate) for candidate in unique_candidates]
    raise AssertionError("No representative molecules found")


class GenerationCorrectnessTests(unittest.TestCase):
    def test_build_structure_counts_c2_to_c5(self):
        expected_counts = {
            2: 2,
            3: 11,
            4: 110,
            5: 1072,
        }
        for carbon_count, expected_count in expected_counts.items():
            with self.subTest(carbon_count=carbon_count):
                self.assertEqual(build_count_for_carbon(carbon_count), expected_count)

    def test_unique_mols_matches_reference_order(self):
        molecules = representative_molecules()
        reference = list(reference_unique_mols(molecules))
        optimized = list(molecule_transformations.unique_mols(molecules))

        self.assertEqual(len(optimized), len(reference))
        for optimized_mol, reference_mol in zip(optimized, reference):
            self.assertTrue(np.array_equal(optimized_mol.bonds, reference_mol.bonds))


if __name__ == "__main__":
    unittest.main()
