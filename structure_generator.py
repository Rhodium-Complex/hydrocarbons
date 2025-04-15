import numpy as np
import graph_utils, molecule, molecule_transformations
from typing import Generator

def build_carbon_hydrogen_combination(c_num: int, h_num: int) -> Generator[np.ndarray, None, None]:
    """
    与えられた炭素数と水素数から、可能な組み合わせを生成します。
    
    各炭素は最大4つの結合を形成でき、未使用の結合は水素と結合すると仮定します。
    
    Args:
        c_num: 炭素の数
        h_num: 水素の数
    Returns:
        generator: 炭素と水素の組み合わせの生成器
    """
    def generate_combinations(c_num: int, h_num: int, current_combination: list) -> Generator[np.ndarray, None, None]:
        if c_num == 0 and h_num == 0:
            yield np.array(current_combination)
            return
        if c_num < 0 or h_num < 0:
            return

        current_hydrogen_count = 1
        if len(current_combination) > 0:
            current_hydrogen_count = current_combination[0]
            
        # 隣接する水素の数が1、2、3、4のそれぞれのケースを処理
        for current_hydrogen in range(current_hydrogen_count, 5):
            new_combination = [current_hydrogen] + current_combination
            yield from generate_combinations(c_num - 1, h_num - (4 - current_hydrogen), new_combination)

    yield from generate_combinations(c_num, h_num, [])

def create_single_bonds_map(carbon_type_list: list) -> Generator[np.ndarray, None, None]:
    """
    単結合の隣接行列を生成する。
    Args:
        carbon_type_list: 炭素の種類のリスト
    Returns:
        generator: 単結合の隣接行列の生成器
    
    """
    
    def generate_child_nodes(adjacency_matrix, bond_degrees, visitable_nodes, start_index=0) -> Generator[np.ndarray, None, None]:
        """
        子ノードを生成するための再帰関数。
        
        Args:
            adj_matrix: 現在の隣接行列
            bond_degrees: 現在の結合次数
            is_unvisitable: 現在の訪問可能なノードのリスト
            start_index: 探索を開始するインデックス
        Yields:
            adj_matrix: 更新された隣接行列
        """
        unprocessed = np.where(bond_degrees > 0)[0]  # 未処理の最初のindexを取得。
        if len(unprocessed) == 0:         # 全点が次数を満足している場合、探索終了。
            yield adjacency_matrix
            return
        next_node_index = unprocessed[0]

        # 未処理の点が最後の1つだった場合、探索失敗。
        if next_node_index == len(bond_degrees): return

        for bond_candidate_index in range(start_index, len(bond_degrees)):
            if adjacency_matrix[bond_candidate_index][next_node_index] == 1 or bond_degrees[bond_candidate_index] == 0:
                continue
            
            # 候補が探索済みもしくは未探索のメチン、メチレン、メタンの最初でなければスキップ
            # これは同型の分子構造を抑制するために必要な制約です
            if visitable_nodes[bond_candidate_index] and bond_candidate_index-next_node_index>1: continue

            (updated_adj_matrix, updated_bond_counts, current_visit_flags) = (adjacency_matrix.copy(), bond_degrees.copy(), visitable_nodes.copy())

            (updated_adj_matrix[next_node_index][bond_candidate_index], updated_adj_matrix[bond_candidate_index][next_node_index]) = (1, 1)
            (updated_bond_counts[next_node_index], updated_bond_counts[bond_candidate_index]) = (updated_bond_counts[next_node_index] - 1, updated_bond_counts[bond_candidate_index] - 1)
            current_visit_flags[bond_candidate_index + 1] = False

            yield from generate_child_nodes(
                updated_adj_matrix, 
                updated_bond_counts, 
                current_visit_flags, 
                bond_candidate_index if updated_bond_counts[next_node_index] else 0)
    
    Num_carbon_types = len(carbon_type_list)
    visited = (
        [False]
        + [carbon_type_list[i + 1] == carbon_type_list[i] for i in range(Num_carbon_types - 1)]
        + [True]
    )
    visited[1] = False
    
    yield from generate_child_nodes(
        np.full((Num_carbon_types, Num_carbon_types), 0, dtype="i4"),
        carbon_type_list,
        visited
    )

def build_structure(input_structure):
    """
    炭素タイプのリストから可能な分子構造を構築します。
    
    Args:
        input_structure: 炭素の種類（水素結合数）を表すnumpy配列
        
    Returns:
        list: 一意な分子構造のリスト（二次元リスト形式）
    """
    candidate_mols = [graph_utils.canonicalize(i) for i in create_single_bonds_map(input_structure) if graph_utils.is_connected_graph(i)]
    unique_mols_list = molecule_transformations.unique_mols(molecule.Molecule(i) for i in np.unique(candidate_mols, axis=0))
    return [[molecule_entity] for molecule_entity in unique_mols_list]