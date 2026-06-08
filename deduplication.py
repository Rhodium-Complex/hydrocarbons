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
    take = np.take
    for molecule_obj in molecule_matrix:
        molecule_fingerprint = molecule_obj.fingerprint
        bonds_for_keys = np.ascontiguousarray(molecule_obj.bonds, dtype=np.uint8)
        atom_count = len(molecule_obj)
        labeled_key = (atom_count, bonds_for_keys.tobytes())
        if labeled_key in seen_labeled_molecules:
            continue
        seen_labeled_molecules.add(labeled_key)

        bond_fingerprints = graph_utils.morgan(molecule_obj.bonds)
        bond_fingerprints = [
            list(np.where(bond_fingerprints == fingerprint)[0])
            for fingerprint in np.unique(bond_fingerprints)[::-1]
        ]
        canonical_permutation = tuple(itertools.chain.from_iterable(bond_fingerprints))

        records_by_size = unique_molecule_records.get(molecule_fingerprint)
        if records_by_size is None:
            records_by_size = {}
            unique_molecule_records[molecule_fingerprint] = records_by_size

        record_keys = records_by_size.setdefault(atom_count, set())
        if record_keys:
            if isomorphism.has_permutation_match(
                bonds_for_keys,
                bond_fingerprints,
                record_keys,
            ):
                continue

        canonical_key = take(
            take(bonds_for_keys, canonical_permutation, axis=0),
            canonical_permutation,
            axis=1,
        ).tobytes()
        record_keys.add(canonical_key)
        yield molecule_obj
