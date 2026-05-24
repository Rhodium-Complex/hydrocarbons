from collections.abc import Generator, Iterable
import itertools

import numpy as np

import graph_utils
import molecule

try:
    from hydrocarbon_rust import has_permutation_match as _rust_has_permutation_match
except ImportError:
    _rust_has_permutation_match = None


def _permuted_matrix_key(bonds: np.ndarray, permutation: tuple) -> tuple:
    permuted_bonds = np.take(np.take(bonds, permutation, axis=0), permutation, axis=1)
    return (len(permutation), permuted_bonds.tobytes())


def _matching_prefix_records(
    matrix_bytes: bytes,
    num_atoms: int,
    permutation: tuple,
    record_keys: tuple,
) -> tuple:
    matched_records = []
    for record_key in record_keys:
        for row_position, row_index in enumerate(permutation):
            matrix_row_offset = row_index * num_atoms
            record_row_offset = row_position * num_atoms
            for column_position, column_index in enumerate(permutation):
                matrix_value = matrix_bytes[matrix_row_offset + column_index]
                record_value = record_key[record_row_offset + column_position]
                if matrix_value != record_value:
                    break
            else:
                continue
            break
        else:
            matched_records.append(record_key)
    return tuple(matched_records)


def _has_permutation_match_python(
    matrix_bytes: bytes,
    num_atoms: int,
    permutation_groups: tuple,
    record_keys: tuple,
) -> bool:
    def search(
        group_index: int,
        permutation_prefix: tuple,
        candidate_records: tuple,
    ) -> bool:
        if group_index == len(permutation_groups):
            return bool(candidate_records)

        for group_permutation in permutation_groups[group_index]:
            next_permutation = permutation_prefix + group_permutation
            next_candidate_records = _matching_prefix_records(
                matrix_bytes,
                num_atoms,
                next_permutation,
                candidate_records,
            )
            if next_candidate_records and search(
                group_index + 1,
                next_permutation,
                next_candidate_records,
            ):
                return True
        return False

    return search(0, (), record_keys)


def _has_permutation_match(
    bonds: np.ndarray,
    bond_fingerprints: list,
    permutation_groups: tuple,
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
        permutation_groups,
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
        permutation_groups = tuple(
            tuple(itertools.permutations(group)) for group in bond_fingerprints
        )
        if _has_permutation_match(
            bonds_for_keys,
            bond_fingerprints,
            permutation_groups,
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
                yield molecule.Molecule(tmp)

    def generate_molecules(molecule_groups):
        for row in molecule_groups:
            yield from dehydrogenase(row)

    return list(unique_mols(generate_molecules(dehydro_stream)))
