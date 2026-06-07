"""Module for generating dehydrogenated molecular variants."""
from collections.abc import Iterable

import numpy as np

import deduplication
import molecule


def unique_dehydro_mols(
    molecules: Iterable[molecule.Molecule],
) -> list[molecule.Molecule]:
    """Return unique molecules after one possible dehydrogenation step."""
    max_bond_order = 3

    def generate_candidates(source_molecules: Iterable[molecule.Molecule]):
        seen_candidate_keys = set()
        for source_molecule in source_molecules:
            source_bonds = source_molecule.bonds
            atom_count = len(source_bonds)
            bond_orders = np.sum(source_bonds, axis=1)

            for atom1 in range(atom_count):
                if bond_orders[atom1] > max_bond_order:
                    continue
                for atom2 in range(atom1 + 1, atom_count):
                    bond_order = source_bonds[atom1][atom2]
                    if bond_order == 0:
                        continue
                    if bond_order == max_bond_order:
                        continue
                    if bond_orders[atom2] > max_bond_order:
                        continue

                    candidate_bonds = source_bonds.copy()
                    candidate_bonds[atom1][atom2] += 1
                    candidate_bonds[atom2][atom1] += 1
                    candidate_key = np.ascontiguousarray(
                        candidate_bonds,
                        dtype=np.uint8,
                    ).tobytes()
                    if candidate_key in seen_candidate_keys:
                        continue
                    seen_candidate_keys.add(candidate_key)
                    yield molecule.Molecule(candidate_bonds)

    return list(deduplication.unique_mols(generate_candidates(molecules)))
