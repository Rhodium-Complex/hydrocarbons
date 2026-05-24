"""Module for generating unique molecular structures and their dehydrogenated variants."""
from collections.abc import Generator, Iterable
import importlib
from importlib.util import find_spec
import itertools
from typing import Callable

import numpy as np

import graph_utils
import molecule

_rust_module = (
    importlib.import_module("hydrocarbon_rust")
    if find_spec("hydrocarbon_rust")
    else None
)
_rust_has_permutation_match: Callable | None = (
    getattr(_rust_module, "has_permutation_match", None)
    if _rust_module is not None
    else None
)


def _permuted_matrix_key(bonds: np.ndarray, permutation: tuple) -> tuple:
    permuted_bonds = np.take(np.take(bonds, permutation, axis=0), permutation, axis=1)
    return (len(permutation), permuted_bonds.tobytes())


def _has_permutation_match_python(
    matrix_bytes: bytes,
    num_atoms: int,
    bond_fingerprints: list,
    record_keys: tuple,
) -> bool:
    allowed_candidates_by_position = []
    for group in bond_fingerprints:
        group_candidates = tuple(group)
        allowed_candidates_by_position.extend([group_candidates] * len(group))

    permutation = [-1] * num_atoms
    used_candidates = [False] * num_atoms

    def compatible_records(
        record_position: int,
        matrix_index: int,
        candidate_records: tuple,
    ) -> tuple:
        matched_records = []
        matrix_row_offset = matrix_index * num_atoms
        for record_key in candidate_records:
            record_row_offset = record_position * num_atoms
            for previous_position in range(record_position):
                previous_matrix_index = permutation[previous_position]
                if (
                    matrix_bytes[matrix_row_offset + previous_matrix_index]
                    != record_key[record_row_offset + previous_position]
                ):
                    break
            else:
                matched_records.append(record_key)
        return tuple(matched_records)

    def search(record_position: int, candidate_records: tuple) -> bool:
        if record_position == num_atoms:
            return bool(candidate_records)

        for matrix_index in allowed_candidates_by_position[record_position]:
            if used_candidates[matrix_index]:
                continue

            next_candidate_records = compatible_records(
                record_position,
                matrix_index,
                candidate_records,
            )
            if not next_candidate_records:
                continue

            permutation[record_position] = matrix_index
            used_candidates[matrix_index] = True
            if search(record_position + 1, next_candidate_records):
                return True
            used_candidates[matrix_index] = False
            permutation[record_position] = -1
        return False

    return search(0, record_keys)


def _has_permutation_match(
    bonds: np.ndarray,
    bond_fingerprints: list,
    record_keys: set,
) -> bool:
    num_atoms = len(bonds)
    if not record_keys:
        return False

    matrix_bytes = bonds.tobytes()
    if _rust_has_permutation_match is not None:
        return _rust_has_permutation_match(
            matrix_bytes,
            num_atoms,
            bond_fingerprints,
            tuple(record_keys),
        )

    return _has_permutation_match_python(
        matrix_bytes,
        num_atoms,
        bond_fingerprints,
        tuple(record_keys),
    )


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
            canonical_key = _permuted_matrix_key(
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
        if _has_permutation_match(
            bonds_for_keys,
            bond_fingerprints,
            record_keys,
        ):
            continue

        canonical_key = _permuted_matrix_key(
            bonds_for_keys,
            canonical_permutation,
        )[1]
        record_keys.add(canonical_key)
        yield molecule_obj
    del unique_molecule_records


def unique_dehydro_mols(dehydro_stream) -> list:
    """Return unique molecules after one possible dehydrogenation step."""

    def dehydrogenase(mol):
        max_bond_order = 3
        bond_orders = np.sum(mol.bonds, axis=1)

        for i in range(len(mol)):
            if bond_orders[i] > max_bond_order:
                continue
            for j in range(i + 1, len(mol)):
                if mol.bonds[i][j] == 0:
                    continue
                if mol.bonds[i][j] == max_bond_order:
                    continue
                if bond_orders[j] > max_bond_order:
                    continue

                tmp = mol.bonds.copy()
                tmp[i][j] += 1
                tmp[j][i] += 1
                yield tmp

    def generate_molecules(molecule_groups):
        seen_candidate_keys = set()
        for row in molecule_groups:
            for candidate_bonds in dehydrogenase(row):
                candidate_key = np.ascontiguousarray(
                    candidate_bonds,
                    dtype=np.uint8,
                ).tobytes()
                if candidate_key in seen_candidate_keys:
                    continue
                seen_candidate_keys.add(candidate_key)
                yield molecule.Molecule(candidate_bonds)

    return list(unique_mols(generate_molecules(dehydro_stream)))
