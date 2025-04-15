import numpy as np
from collections import deque

MAX_PRIORITY = 11
PRIORITY_THRESHOLD = 10

def is_connected_graph(adj_matrix: np.ndarray) -> bool:
    """
    グラフが連結かどうかを確認します。
    
    Args:
        adj_matrix: グラフの隣接行列（正方行列であること）
        
    Returns:
        bool: グラフが連結であればTrue、そうでなければFalse
    """
    visited = np.zeros(len(adj_matrix), dtype=bool)
    visited[0] = True
    for _ in range(len(adj_matrix) - 1):
        visited |= np.dot(adj_matrix, visited).astype(bool)
        if visited.all():
            return True
    return False

def morgan(adjacency_matrix: np.ndarray) -> np.ndarray:
    """
    Morgan アルゴリズムを使用して頂点に識別ラベルを割り当てます。
    
    Args:
        adjacency_matrix: グラフの隣接行列
        
    Returns:
        numpy配列: 各頂点のMorganラベル
    """
    
    # マジックナンバーを定数として定義
    BASE = 4  # 基数として使用(結合は4を超えないため)
    
    connect_map = adjacency_matrix
    flag = 0
    
    # 初期ラベルは各頂点の接続数に基づく
    node_labels = BASE ** np.sum(connect_map, axis = 1) # 4は基数として使用
    unique_node_labels = np.unique(node_labels)
    while len(unique_node_labels) > flag:
        flag = len(unique_node_labels)
        for i in range(len(unique_node_labels)):
            node_labels[node_labels == unique_node_labels[i]] = BASE ** (flag-i)
        node_labels = connect_map @ node_labels
        unique_node_labels = np.unique(node_labels)
    return node_labels

def assign_priorities(canonical_labels : np.ndarray, adjacency_matrix: np.ndarray, max_priority: int = 11) -> np.ndarray:
    """
    ノードに優先順位を割り当てます。

    Args:
        cano_num: Morganアルゴリズムで生成されたノードラベル
        adjacency_matrix: グラフの隣接行列
        max_priority: 最大優先順位

    Returns:
        np.ndarray: 優先順位が割り当てられたノードラベル
    """
    node_priorities = np.zeros(len(canonical_labels ))
    weighted_matrix = (adjacency_matrix > 0) * canonical_labels  + adjacency_matrix
    node_priorities[canonical_labels .argmax()] = max_priority

    def priority_generator():
        for i in range(max_priority)[::-1]:
            yield i

    priority_gen = priority_generator()
    stack = deque([i for i in weighted_matrix[canonical_labels .argmax()].argsort() if weighted_matrix[canonical_labels .argmax()][i]])
    while stack:
        i = stack.popleft()
        if node_priorities[i] == 0:
            node_priorities[i] = next(priority_gen)
            stack.extend([j for j in weighted_matrix[i].argsort() if weighted_matrix[i][j]])

    return node_priorities * 10

def build_main_branch(canonical_labels: np.ndarray, binary_connection: np.ndarray, adjacency_matrix: np.ndarray, priority_threshold: int = 10) -> list:
    """
    主枝を構築します。

    Args:
        canonical_labels: 優先順位が割り当てられたノードラベル
        binary_connection: 二値化された接続情報
        adjacency_matrix: グラフの隣接行列
        priority_threshold: 優先度の閾値

    Returns:
        list: 主枝を構成するノードの順序
    """
    current_node = canonical_labels.argmax()
    node_order = [current_node]
    canonical_labels[current_node] = 0
    edge_values = binary_connection[current_node] * canonical_labels + adjacency_matrix[current_node]

    while max(edge_values) > priority_threshold:
        current_node = edge_values.argmax()
        node_order = [current_node] + node_order
        canonical_labels[current_node] = 0
        edge_values = binary_connection[current_node] * canonical_labels + adjacency_matrix[current_node]

    return node_order

def build_side_branches(node_order: list, canonical_labels: np.ndarray, binary_connection: np.ndarray, adjacency_matrix: np.ndarray, priority_threshold: int = 10) -> list:
    """
    側枝を構築します。

    Args:
        node_order: 主枝を構成するノードの順序
        canonical_labels: 優先順位が割り当てられたノードラベル
        binary_connection: 二値化された接続情報
        adjacency_matrix: グラフの隣接行列
        priority_threshold: 優先度の閾値

    Returns:
        list: 主枝と側枝を含む全ノードの順序
    """
    counter = len(node_order) - 1
    edge_values = binary_connection[node_order[counter]] * canonical_labels + adjacency_matrix[node_order[counter]]
    if max(edge_values) <= priority_threshold:
        counter = 0

    while max(canonical_labels) > priority_threshold:
        while max(edge_values) <= priority_threshold:
            if counter >= len(node_order):
                break
            edge_values = binary_connection[node_order[counter]] * canonical_labels + adjacency_matrix[node_order[counter]]
            counter += 1

        current_node = edge_values.argmax()
        node_order = node_order + [current_node]
        canonical_labels[current_node] = 0
        edge_values = binary_connection[current_node] * canonical_labels
        counter += 1

        if max(edge_values) == 0:
            counter = 0

    return node_order

def canonicalize(adjacency_matrix: np.ndarray) -> np.ndarray:
    """
    グラフの正規化表現を生成します。
    同型なグラフが同一の表現を持つように変換します。

    Args:
        adjacency_matrix: グラフの隣接行列

    Returns:
        numpy配列: 正規化された隣接行列
    """
    binary_connection = (adjacency_matrix > 0)
    canonical_labels = morgan(adjacency_matrix)

    # ノードに優先順位を割り当てる
    canonical_labels = assign_priorities(canonical_labels, adjacency_matrix, MAX_PRIORITY)

    # 主枝を構築
    node_order = build_main_branch(canonical_labels, binary_connection, adjacency_matrix, PRIORITY_THRESHOLD)

    # 側枝を構築
    node_order = build_side_branches(node_order, canonical_labels, binary_connection, adjacency_matrix, PRIORITY_THRESHOLD)

    return np.array(adjacency_matrix[:, node_order][node_order, :])