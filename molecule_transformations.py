import graph_utils,molecule
import itertools
import numpy as np
from typing import Generator

try:
    from hydrocarbon_rust import has_permutation_match as _rust_has_permutation_match
except ImportError:
    _rust_has_permutation_match = None


def _permuted_matrix_key(bonds: np.ndarray, permutation: tuple) -> tuple:
    permuted_bonds = np.take(np.take(bonds, permutation, axis=0), permutation, axis=1)
    return (len(permutation), permuted_bonds.tobytes())


def _matching_prefix_records(matrix_bytes: bytes, num_atoms: int, permutation: tuple, record_keys: tuple) -> tuple:
    matched_records = []
    for record_key in record_keys:
        for row_position, row_index in enumerate(permutation):
            matrix_row_offset = row_index * num_atoms
            record_row_offset = row_position * num_atoms
            for column_position, column_index in enumerate(permutation):
                if matrix_bytes[matrix_row_offset + column_index] != record_key[record_row_offset + column_position]:
                    break
            else:
                continue
            break
        else:
            matched_records.append(record_key)
    return tuple(matched_records)


def _has_permutation_match_python(matrix_bytes: bytes, num_atoms: int, permutation_groups: tuple, record_keys: tuple) -> bool:
    def search(group_index: int, permutation_prefix: tuple, candidate_records: tuple) -> bool:
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
            if next_candidate_records and search(group_index + 1, next_permutation, next_candidate_records):
                return True
        return False

    return search(0, (), record_keys)


def _has_permutation_match(bonds: np.ndarray, bond_fingerprints: list, permutation_groups: tuple, record_keys: set) -> bool:
    num_atoms = len(bonds)
    full_record_keys = tuple(key for key_size, key in record_keys if key_size == num_atoms)
    if not full_record_keys:
        return False

    matrix_bytes = bonds.tobytes()
    if _rust_has_permutation_match is not None:
        return _rust_has_permutation_match(matrix_bytes, num_atoms, bond_fingerprints, full_record_keys)

    return _has_permutation_match_python(matrix_bytes, num_atoms, permutation_groups, full_record_keys)

def unique_mols(molecule_matrix: list) -> Generator[molecule.Molecule, None, None]:
    """
    分子のリストから一意な分子構造のみを抽出します。
    
    このメソッドは、グラフ同型判定を使用して、異なる構造表現を持つが
    化学的に同一の分子を除去します。Morganアルゴリズムを使用して
    頂点ラベリングを行い、それに基づいて同型判定を行います。
    
    Args:
        molecule_list: 分子オブジェクトのリスト
        
    Yields:
        一意な分子構造
    """
    unique_molecule_records = {}
    for molecule_obj  in molecule_matrix:
        molecule_fingerprint = molecule_obj .fingerprint
        bonds_for_keys = np.ascontiguousarray(molecule_obj .bonds, dtype=np.uint8)
        bond_fingerprints = graph_utils.morgan(molecule_obj .bonds)
        bond_fingerprints = [list(np.where(bond_fingerprints == i)[0]) for i in np.unique(bond_fingerprints)[::-1]]
        permutation_groups = tuple(tuple(itertools.permutations(group)) for group in bond_fingerprints)
        
        for group_permutations in permutation_groups:
            for permutation in group_permutations:
                canonical_permutation = tuple(permutation)
                canonical_key = _permuted_matrix_key(bonds_for_keys, canonical_permutation)
                if molecule_fingerprint not in unique_molecule_records:
                    unique_molecule_records[molecule_fingerprint] = {canonical_key}
                    yield molecule_obj 
                    continue
                
                if not _has_permutation_match(
                    bonds_for_keys,
                    bond_fingerprints,
                    permutation_groups,
                    unique_molecule_records[molecule_fingerprint],
                ):
                    unique_molecule_records[molecule_fingerprint].add(canonical_key)
                    yield molecule_obj 
    del unique_molecule_records

def unique_dehydro_mols(dehydro_stream)-> list:
    """
    脱水素化された分子のリストから一意な分子構造のみを抽出します。
    このメソッドは、グラフ同型判定を使用して、異なる構造表現を持つが
    化学的に同一の分子を除去します。Morganアルゴリズムを使用して
    頂点ラベリングを行い、それに基づいて同型判定を行います。
    Args:
        dehydro_stream: 脱水素化された分子のリスト
    Returns:
        list: 一意な分子構造のリスト（二次元リスト形式）
    """
    
    def dehydrogenase(mol):
        MAX_BOND_ORDER = 3
        
        # 結合次数の増加（脱水素化）を試みる
        # 化学的には2つの原子間の結合電子対が増えることを表現
        for i in range(len(mol)):
            # 原子iの結合の合計が最大結合次数を超える場合はスキップ
            if sum(mol.bonds[i]) > MAX_BOND_ORDER: continue
            for j in range(i + 1, len(mol)):
                if mol.bonds[i][j] == 0: continue
                if mol.bonds[i][j] == MAX_BOND_ORDER: continue
                if sum(mol.bonds[j]) > MAX_BOND_ORDER: continue
                
                tmp = mol.bonds.copy()
                tmp[i][j] +=1
                tmp[j][i] +=1
                yield molecule.Molecule(tmp)
                
    def generate_molecules(j):
        for row in j:
            yield from dehydrogenase(row)
    
    return list(unique_mols(generate_molecules(dehydro_stream)))
