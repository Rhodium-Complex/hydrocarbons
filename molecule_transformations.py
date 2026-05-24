"""Module for generating unique molecular structures and their dehydrogenated variants."""
from collections.abc import Generator, Iterable
import itertools

import numpy as np

import graph_utils
import isomorphism
import molecule


def unique_mols(
    molecule_matrix: Iterable[molecule.Molecule],
) -> Generator[molecule.Molecule, None, None]:
    """Yield structurally unique molecules from an iterable of molecules."""
    unique_molecule_records = {}
    seen_labeled_molecules = set()
    for molecule_obj in molecule_matrix:
        molecule_fingerprint = molecule_obj.fingerprint
        bonds_for_keys = np.ascontiguousarray(molecule_obj.bonds, dtype=np.uint8)
        labeled_key = (len(molecule_obj), bonds_for_keys.tobytes())
        if labeled_key in seen_labeled_molecules:
            continue
        seen_labeled_molecules.add(labeled_key)

        bond_fingerprints = graph_utils.morgan(molecule_obj.bonds)
        bond_fingerprints = [
            list(np.where(bond_fingerprints == fingerprint)[0])
            for fingerprint in np.unique(bond_fingerprints)[::-1]
        ]
        canonical_permutation = tuple(itertools.chain.from_iterable(bond_fingerprints))

        if molecule_fingerprint not in unique_molecule_records:
            canonical_key = isomorphism.permuted_matrix_key(
                bonds_for_keys,
                canonical_permutation,
            )[1]
            unique_molecule_records[molecule_fingerprint] = {
                len(molecule_obj): {canonical_key}
            }
            yield molecule_obj
            continue

        record_keys = unique_molecule_records[molecule_fingerprint].setdefault(
            len(molecule_obj),
            set(),
        )
        if isomorphism.has_permutation_match(
            bonds_for_keys,
            bond_fingerprints,
            record_keys,
        ):
            continue

        canonical_key = isomorphism.permuted_matrix_key(
            bonds_for_keys,
            canonical_permutation,
        )[1]
        record_keys.add(canonical_key)
        yield molecule_obj
    del unique_molecule_records


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

    return list(unique_mols(generate_candidates(molecules)))
