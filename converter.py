"""Module for converting molecule bond matrices into SMILES-like strings."""
import stereochemistry

BOND_SYMBOLS = {
    1: "",
    2: "=",
    3: "#",
}


def mat2smiles(mat):
    """Convert a molecule bond matrix into a compact SMILES-like string."""
    bond_matrix = mat.bonds
    smiles_parts = [""] * len(bond_matrix[0])
    ring_placeholder = 1

    for atom in range(len(bond_matrix[0]))[::-1]:
        has_parent_bond = False
        smiles_parts[atom] = "C" + smiles_parts[atom]
        for neighbor in range(atom)[::-1]:
            bond = bond_matrix[atom][neighbor]
            if bond == 0:
                continue
            if has_parent_bond:
                token = "{" + str(ring_placeholder) + "}"
                smiles_parts[atom] += BOND_SYMBOLS[bond] + token
                smiles_parts[neighbor] = token + smiles_parts[neighbor]
                ring_placeholder += 1
            else:
                smiles_parts[atom] = BOND_SYMBOLS[bond] + smiles_parts[atom]
            has_parent_bond = True

        child_smiles = ""
        has_child = False
        for neighbor in range(atom, len(bond_matrix[0]))[::-1]:
            if bond_matrix[atom][neighbor] == 0:
                continue
            if smiles_parts[neighbor] == "":
                continue
            if has_child:
                child_smiles = "(" + smiles_parts[neighbor] + ")" + child_smiles
            else:
                child_smiles = smiles_parts[neighbor]
            smiles_parts[neighbor] = ""
            has_child = True
        smiles_parts[atom] += child_smiles

    smiles = smiles_parts[0]

    counter = 1

    def ring_number(num):
        if num < 10:
            return str(num)
        return "%" + str(num)

    while "{" in smiles:
        start = smiles.find("{")
        end = smiles.find("}")
        ring_token = smiles[start:end + 1]
        smiles = smiles.replace(ring_token, ring_number(counter), 2)
        counter = counter + 1
    return smiles


def _stereo_main_path(bond_matrix, assignment, analysis):
    single_label = assignment.single_label
    if single_label is None:
        return None

    double_bond = analysis.double_bond_for_assignment(assignment)
    if double_bond is None:
        return None
    if not double_bond.has_carbon_high_ligands:
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
    _atom1, _atom2, label = single_label
    return main_path, label


def _acyclic_stereo_smiles(mat, assignment, analysis):
    bond_matrix = mat.bonds
    if int((bond_matrix > 0).sum() // 2) != len(bond_matrix) - 1:
        return None

    path_and_label = _stereo_main_path(bond_matrix, assignment, analysis)
    if path_and_label is None:
        return None
    main_path, label = path_and_label

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
        return BOND_SYMBOLS.get(bond_matrix[left][right])

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

    acyclic_smiles = _acyclic_stereo_smiles(mat, assignment, analysis)
    if acyclic_smiles is not None:
        return acyclic_smiles

    assignment_suffix = "".join(
        f" [{label}:{atom1 + 1}-{atom2 + 1}]"
        for atom1, atom2, label in assignment.labels
    )
    return mat2smiles(mat) + assignment_suffix


def mat2smiles_variants(mat, include_stereo: bool = False) -> list[str]:
    """Return one or more SMILES outputs for a molecule."""
    if not include_stereo:
        return [mat2smiles(mat)]
    analysis = stereochemistry.analyze_ez(mat)
    return [
        mat2stereo_smiles(mat, assignment, analysis)
        for assignment in analysis.assignments
    ]
