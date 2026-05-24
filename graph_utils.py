from collections import deque

import numpy as np

MAX_PRIORITY = 11
PRIORITY_THRESHOLD = 10


def is_connected_graph(adj_matrix: np.ndarray) -> bool:
    """Return whether an adjacency matrix describes a connected graph."""
    visited = np.zeros(len(adj_matrix), dtype=bool)
    visited[0] = True
    for _ in range(len(adj_matrix) - 1):
        visited |= np.dot(adj_matrix, visited).astype(bool)
        if visited.all():
            return True
    return False


def morgan(adjacency_matrix: np.ndarray) -> np.ndarray:
    """Assign Morgan-style labels to graph vertices."""
    base = 4

    connect_map = adjacency_matrix
    flag = 0

    node_labels = base ** np.sum(connect_map, axis=1)
    unique_node_labels = np.unique(node_labels)
    while len(unique_node_labels) > flag:
        flag = len(unique_node_labels)
        for i in range(len(unique_node_labels)):
            node_labels[node_labels == unique_node_labels[i]] = base ** (flag - i)
        node_labels = connect_map @ node_labels
        unique_node_labels = np.unique(node_labels)
    return node_labels


def assign_priorities(
    canonical_labels: np.ndarray,
    adjacency_matrix: np.ndarray,
    max_priority: int = 11,
) -> np.ndarray:
    """Assign traversal priorities from canonical labels and bond orders."""
    node_priorities = np.zeros(len(canonical_labels))
    weighted_matrix = (adjacency_matrix > 0) * canonical_labels + adjacency_matrix
    root_index = canonical_labels.argmax()
    node_priorities[root_index] = max_priority

    def priority_generator():
        for i in range(max_priority)[::-1]:
            yield i

    priority_gen = priority_generator()
    stack = deque(
        [
            i
            for i in weighted_matrix[root_index].argsort()
            if weighted_matrix[root_index][i]
        ]
    )
    while stack:
        i = stack.popleft()
        if node_priorities[i] == 0:
            node_priorities[i] = next(priority_gen)
            stack.extend(
                j for j in weighted_matrix[i].argsort() if weighted_matrix[i][j]
            )

    return node_priorities * 10


def build_main_branch(
    canonical_labels: np.ndarray,
    binary_connection: np.ndarray,
    adjacency_matrix: np.ndarray,
    priority_threshold: int = 10,
) -> list:
    """Build the ordered main branch for canonicalization."""
    current_node = canonical_labels.argmax()
    node_order = [current_node]
    canonical_labels[current_node] = 0
    edge_values = (
        binary_connection[current_node] * canonical_labels
        + adjacency_matrix[current_node]
    )

    while max(edge_values) > priority_threshold:
        current_node = edge_values.argmax()
        node_order = [current_node] + node_order
        canonical_labels[current_node] = 0
        edge_values = (
            binary_connection[current_node] * canonical_labels
            + adjacency_matrix[current_node]
        )

    return node_order


def build_side_branches(
    node_order: list,
    canonical_labels: np.ndarray,
    binary_connection: np.ndarray,
    adjacency_matrix: np.ndarray,
    priority_threshold: int = 10,
) -> list:
    """Append side branches to the canonical node order."""
    counter = len(node_order) - 1
    edge_values = (
        binary_connection[node_order[counter]] * canonical_labels
        + adjacency_matrix[node_order[counter]]
    )
    if max(edge_values) <= priority_threshold:
        counter = 0

    while max(canonical_labels) > priority_threshold:
        while max(edge_values) <= priority_threshold:
            if counter >= len(node_order):
                break
            edge_values = (
                binary_connection[node_order[counter]] * canonical_labels
                + adjacency_matrix[node_order[counter]]
            )
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
    """Return a canonical adjacency matrix for an isomorphic graph class."""
    binary_connection = adjacency_matrix > 0
    canonical_labels = morgan(adjacency_matrix)

    canonical_labels = assign_priorities(
        canonical_labels,
        adjacency_matrix,
        MAX_PRIORITY,
    )
    node_order = build_main_branch(
        canonical_labels,
        binary_connection,
        adjacency_matrix,
        PRIORITY_THRESHOLD,
    )
    node_order = build_side_branches(
        node_order,
        canonical_labels,
        binary_connection,
        adjacency_matrix,
        PRIORITY_THRESHOLD,
    )

    return np.array(adjacency_matrix[:, node_order][node_order, :])
