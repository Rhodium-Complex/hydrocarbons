"""Isomorphism helpers for molecular bond matrices."""
import importlib
from importlib.util import find_spec
from typing import Callable

import numpy as np

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


def permuted_matrix_key(bonds: np.ndarray, permutation: tuple) -> tuple:
    """Return a size-tagged byte key for a permuted bond matrix."""
    permuted_bonds = np.take(np.take(bonds, permutation, axis=0), permutation, axis=1)
    return (len(permutation), permuted_bonds.tobytes())


def automorphisms(bonds: np.ndarray) -> list[tuple[int, ...]]:
    """Return graph automorphisms for a single bond matrix."""
    size = len(bonds)
    signatures = [
        (
            int(np.sum(bonds[index])),
            tuple(sorted(int(order) for order in bonds[index] if order > 0)),
        )
        for index in range(size)
    ]
    candidates_by_position = [
        tuple(
            candidate
            for candidate, signature in enumerate(signatures)
            if signature == signatures[position]
        )
        for position in range(size)
    ]
    search_order = sorted(
        range(size),
        key=lambda index: len(candidates_by_position[index]),
    )
    mapping = [-1] * size
    used = [False] * size
    results = []

    def is_compatible(source: int, target: int) -> bool:
        for previous_source, previous_target in enumerate(mapping):
            if previous_target == -1:
                continue
            if bonds[source][previous_source] != bonds[target][previous_target]:
                return False
        return True

    def search(depth: int) -> None:
        if depth == size:
            results.append(tuple(mapping))
            return
        source = search_order[depth]
        for target in candidates_by_position[source]:
            if used[target] or not is_compatible(source, target):
                continue
            mapping[source] = target
            used[target] = True
            search(depth + 1)
            used[target] = False
            mapping[source] = -1

    search(0)
    return results


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


def has_permutation_match(
    bonds: np.ndarray,
    bond_fingerprints: list,
    record_keys: set,
) -> bool:
    """Return whether bonds matches any representative under allowed permutations."""
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
