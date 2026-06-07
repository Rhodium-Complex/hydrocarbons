"""Structural deduplication helpers for molecular bond matrices."""
from collections.abc import Iterable
import itertools

import numpy as np

import graph_utils
import isomorphism
import molecule


def unique_mols(
    molecule_matrix: Iterable[molecule.Molecule],
) -> Iterable[molecule.Molecule]:
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
