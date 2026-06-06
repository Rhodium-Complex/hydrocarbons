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


def _tree_stereo_smiles(mat, assignment, analysis):
    if len(assignment.labels) != 1:
        return None

    bond_matrix = mat.bonds
    if int((bond_matrix > 0).sum() // 2) != len(bond_matrix) - 1:
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
    if (
        bond_matrix[double_bond.high_ligand1][double_bond.atom1] != 1
        or bond_matrix[double_bond.high_ligand2][double_bond.atom2] != 1
    ):
        return None

    main_path = (
        double_bond.high_ligand1,
        double_bond.atom1,
        double_bond.atom2,
        double_bond.high_ligand2,
    )
    main_next = {
        main_path[index]: main_path[index + 1]
        for index in range(len(main_path) - 1)
    }
    main_edges = {
        frozenset((main_path[index], main_path[index + 1]))
        for index in range(len(main_path) - 1)
    }
    stereo_single_bonds = {
        frozenset((main_path[0], main_path[1])): "/",
        frozenset((main_path[2], main_path[3])): "/" if label == "E" else "\\",
    }

    def bond_symbol(left, right):
        edge = frozenset((left, right))
        if edge in stereo_single_bonds:
            return stereo_single_bonds[edge]
        bond = bond_matrix[left][right]
        if bond == 1:
            return ""
        elif bond == 2:
            return "="
        elif bond == 3:
            return "#"
        return None

    def emit_atom(atom, parent=None):
        parts = ["C"]
        next_main_atom = main_next.get(atom)
        branch_neighbors = [
            int(neighbor)
            for neighbor in bond_matrix[atom].nonzero()[0]
            if int(neighbor) != parent
            and int(neighbor) != next_main_atom
            and frozenset((atom, int(neighbor))) not in main_edges
        ]
        for neighbor in sorted(branch_neighbors):
            symbol = bond_symbol(atom, neighbor)
            if symbol is None:
                return None
            branch = emit_atom(neighbor, atom)
            if branch is None:
                return None
            parts.append(f"({symbol}{branch})")

        if next_main_atom is not None:
            symbol = bond_symbol(atom, next_main_atom)
            if symbol is None:
                return None
            child = emit_atom(next_main_atom, atom)
            if child is None:
                return None
            parts.append(symbol + child)
        return "".join(parts)

    return emit_atom(main_path[0])


def mat2stereo_smiles(mat, assignment, analysis=None):
    """Convert a molecule and E/Z assignment to a stereo-aware SMILES-like string."""
    if not assignment.labels:
        return mat2smiles(mat)
    if analysis is None:
        analysis = stereochemistry.analyze_ez(mat)

    tree_smiles = _tree_stereo_smiles(mat, assignment, analysis)
    if tree_smiles is not None:
        return tree_smiles

    assignment_suffix = "".join(
        f" [{label}:{atom1 + 1}-{atom2 + 1}]"
        for atom1, atom2, label in assignment.labels
    )
    return mat2smiles(mat) + assignment_suffix
