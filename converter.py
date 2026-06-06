"""Module for converting molecule bond matrices into SMILES-like strings."""
import stereochemistry


def mat2smiles(mat):
    """Convert a molecule bond matrix into a compact SMILES-like string."""
    bond_matrix = mat.bonds
    bond_symbols = ["", "", "=", "#"]
    smiles_parts = [""] * len(bond_matrix[0])
    counter = 1

    for i in range(len(bond_matrix[0]))[::-1]:
        flag = False
        smiles_parts[i] = "C" + smiles_parts[i]
        for j in range(i)[::-1]:
            bond = bond_matrix[i][j]
            if bond == 0:
                continue
            if flag:
                smiles_parts[i] += bond_symbols[bond] + "{" + str(counter) + "}"
                smiles_parts[j] = "{" + str(counter) + "}" + smiles_parts[j]
                counter += 1
            else:
                smiles_parts[i] = bond_symbols[bond] + smiles_parts[i]
            flag = True
        tmp = ""
        flag = False
        for j in range(i, len(bond_matrix[0]))[::-1]:
            if bond_matrix[i][j] == 0:
                continue
            if smiles_parts[j] == "":
                continue
            tmp = "(" + smiles_parts[j] + ")" + tmp if flag else smiles_parts[j]
            smiles_parts[j] = ""
            flag = True
        smiles_parts[i] += tmp
    smiles_parts = smiles_parts[0]

    counter = 1
    def sanit(num):
        if num < 10:
            return str(num)
        return "%" + str(num)

    while "{" in smiles_parts:
        i = smiles_parts.find("{")
        l = smiles_parts.find("}")
        ring_token = smiles_parts[i:l + 1]

        # 正規表現を使用してラベルを置換
        # 例: {1} を適切な SMILES 記法に変換
        smiles_parts = smiles_parts.replace(ring_token, sanit(counter), 2)
        counter = counter + 1
    return smiles_parts


def _simple_path_order(bond_matrix):
    binary_bonds = bond_matrix > 0
    degrees = binary_bonds.sum(axis=1)
    if len(bond_matrix) == 1:
        return [0]
    if degrees.max() > 2:
        return None
    endpoints = [int(index) for index, degree in enumerate(degrees) if degree == 1]
    if len(endpoints) != 2:
        return None

    order = []
    previous = -1
    current = endpoints[0]
    while True:
        order.append(current)
        neighbors = [
            int(neighbor)
            for neighbor in binary_bonds[current].nonzero()[0]
            if int(neighbor) != previous
        ]
        if not neighbors:
            break
        previous, current = current, neighbors[0]
        if current in order:
            return None

    if len(order) != len(bond_matrix):
        return None
    return order


def _path_stereo_smiles(mat, assignment, analysis):
    if len(assignment.labels) != 1:
        return None

    bond_matrix = mat.bonds
    order = _simple_path_order(bond_matrix)
    if order is None:
        return None

    atom1, atom2, label = assignment.labels[0]
    double_bond = analysis.double_bond_by_edge.get(
        (min(atom1, atom2), max(atom1, atom2))
    )
    if double_bond is None:
        return None
    high_ligands = {double_bond.high_ligand1, double_bond.high_ligand2}
    if stereochemistry.HYDROGEN_LIGAND in high_ligands:
        return None

    try:
        double_index = next(
            index
            for index in range(len(order) - 1)
            if {order[index], order[index + 1]} == {atom1, atom2}
        )
    except StopIteration:
        return None

    if double_index == 0 or double_index + 2 >= len(order):
        return None
    if order[double_index - 1] not in high_ligands:
        return None
    if order[double_index + 2] not in high_ligands:
        return None

    left_slash = "/"
    right_slash = "/" if label == "E" else "\\"
    parts = ["C"]
    for index in range(len(order) - 1):
        left = order[index]
        right = order[index + 1]
        bond = bond_matrix[left][right]
        if bond == 1 and index == double_index - 1:
            bond_symbol = left_slash
        elif bond == 1 and index == double_index + 1:
            bond_symbol = right_slash
        elif bond == 1:
            bond_symbol = ""
        elif bond == 2:
            bond_symbol = "="
        elif bond == 3:
            bond_symbol = "#"
        else:
            return None
        parts.append(bond_symbol + "C")
    return "".join(parts)


def mat2stereo_smiles(mat, assignment, analysis=None):
    """Convert a molecule and E/Z assignment to a stereo-aware SMILES-like string."""
    if not assignment.labels:
        return mat2smiles(mat)
    if analysis is None:
        analysis = stereochemistry.analyze_ez(mat)

    path_smiles = _path_stereo_smiles(mat, assignment, analysis)
    if path_smiles is not None:
        return path_smiles

    assignment_suffix = "".join(
        f" [{label}:{atom1 + 1}-{atom2 + 1}]"
        for atom1, atom2, label in assignment.labels
    )
    return mat2smiles(mat) + assignment_suffix
