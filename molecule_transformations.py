"""Module for generating dehydrogenated molecular variants."""
from collections.abc import Iterable

import numpy as np

import deduplication
import molecule


def unique_dehydro_mols(
    molecules: Iterable[molecule.Molecule],
) -> list[molecule.Molecule]:
    """Return unique molecules after one possible dehydrogenation step."""

    def dehydrogenate_once(source_molecule: molecule.Molecule):
        max_bond_order = 3
        bond_orders = np.sum(source_molecule.bonds, axis=1)

        for i in range(len(source_molecule)):
            if bond_orders[i] > max_bond_order:
                continue
            for j in range(i + 1, len(source_molecule)):
                if source_molecule.bonds[i][j] == 0:
                    continue
                if source_molecule.bonds[i][j] == max_bond_order:
                    continue
                if bond_orders[j] > max_bond_order:
                    continue

                candidate_bonds = source_molecule.bonds.copy()
                candidate_bonds[i][j] += 1
                candidate_bonds[j][i] += 1
                yield candidate_bonds

    def generate_candidates(source_molecules: Iterable[molecule.Molecule]):
        seen_candidate_keys = set()
        for source_molecule in source_molecules:
            for candidate_bonds in dehydrogenate_once(source_molecule):
                candidate_key = np.ascontiguousarray(
                    candidate_bonds,
                    dtype=np.uint8,
                ).tobytes()
                if candidate_key in seen_candidate_keys:
                    continue
                seen_candidate_keys.add(candidate_key)
                yield molecule.Molecule(candidate_bonds)

    return list(deduplication.unique_mols(generate_candidates(molecules)))
