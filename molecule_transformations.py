import graph_utils,molecule
import itertools
import numpy as np
from typing import Generator

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
        bond_fingerprints = graph_utils.morgan(molecule_obj .bonds)
        bond_fingerprints = [list(np.where(bond_fingerprints == i)[0]) for i in np.unique(bond_fingerprints)[::-1]]
        
        for bond_fingerprint_group in bond_fingerprints:
            for permutation in itertools.permutations(bond_fingerprint_group):
                canonical_permutation = sum(permutation, ())
                if molecule_fingerprint not in unique_molecule_records:
                    unique_molecule_records[molecule_fingerprint] = [molecule_obj .bonds[:,canonical_permutation][canonical_permutation,:]]
                    yield molecule_obj 
                    continue
                
                for i in itertools.product(*[ itertools.permutations(j) for j in bond_fingerprints]):
                    i = sum(i,())
                    current_bond_structure = molecule_obj .bonds[:,i][i,:]
                    if any(np.array_equal(current_bond_structure, j) for j in unique_molecule_records[molecule_fingerprint]): break
                else:
                    unique_molecule_records[molecule_fingerprint].append(molecule_obj .bonds[:,canonical_permutation][canonical_permutation,:])
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