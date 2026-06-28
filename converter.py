"""Module for converting molecule bond matrices into SMILES-like strings."""
import itertools

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


def _acyclic_stereo_smiles(mat, assignment, analysis):
    """Encode every compatible E/Z label on an acyclic carbon skeleton."""
    bond_matrix = mat.bonds
    atom_count = len(bond_matrix)
    if int((bond_matrix > 0).sum() // 2) != atom_count - 1:
        return None

    double_bond_by_edge = {
        (double_bond.atom1, double_bond.atom2): double_bond
        for double_bond in analysis.double_bonds
    }
    stereo_data = []
    preferred_next = {}
    for atom1, atom2, label in assignment.labels:
        double_bond = double_bond_by_edge.get((min(atom1, atom2), max(atom1, atom2)))
        if double_bond is None or not double_bond.has_carbon_high_ligands:
            return None
        high1 = double_bond.high_ligand1
        high2 = double_bond.high_ligand2
        if bond_matrix[high1][double_bond.atom1] != 1 or bond_matrix[double_bond.atom2][high2] != 1:
            return None
        stereo_data.append((double_bond, label))
        for left, right in zip(
            (high1, double_bond.atom1, double_bond.atom2),
            (double_bond.atom1, double_bond.atom2, high2),
        ):
            preferred_next.setdefault(left, right)

    first_double_bond = stereo_data[0][0]
    root = first_double_bond.high_ligand1
    parent = {root: None}
    children = {atom: [] for atom in range(atom_count)}
    stack = [root]
    while stack:
        atom = stack.pop()
        for neighbor in sorted(int(value) for value in bond_matrix[atom].nonzero()[0]):
            if neighbor == parent[atom]:
                continue
            if neighbor in parent:
                return None
            parent[neighbor] = atom
            children[atom].append(neighbor)
            stack.append(neighbor)
    if len(parent) != atom_count:
        return None

    def emitted_reversed(left, right):
        """Return whether an edge is emitted opposite to left -> right."""
        if parent.get(right) == left:
            return False
        if parent.get(left) == right:
            return True
        raise ValueError("stereo ligand edge is not in the spanning tree")

    constraints = []
    marked_edges = set()
    for double_bond, label in stereo_data:
        ligand_edge1 = frozenset((double_bond.high_ligand1, double_bond.atom1))
        ligand_edge2 = frozenset((double_bond.atom2, double_bond.high_ligand2))
        reversed1 = emitted_reversed(double_bond.high_ligand1, double_bond.atom1)
        reversed2 = emitted_reversed(double_bond.atom2, double_bond.high_ligand2)
        differs = (label == "Z") ^ reversed1 ^ reversed2
        constraints.append((ligand_edge1, ligand_edge2, differs))
        marked_edges.update((ligand_edge1, ligand_edge2))

    marks = {}
    adjacency = {edge: [] for edge in marked_edges}
    for edge1, edge2, differs in constraints:
        adjacency[edge1].append((edge2, differs))
        adjacency[edge2].append((edge1, differs))
    for start in sorted(marked_edges, key=lambda edge: tuple(sorted(edge))):
        if start in marks:
            continue
        marks[start] = False
        pending = [start]
        while pending:
            edge = pending.pop()
            for other, differs in adjacency[edge]:
                expected = marks[edge] ^ differs
                if other in marks and marks[other] != expected:
                    return None
                if other not in marks:
                    marks[other] = expected
                    pending.append(other)

    def bond_symbol(left, right):
        edge = frozenset((left, right))
        if edge in marks:
            return "\\" if marks[edge] else "/"
        return BOND_SYMBOLS.get(bond_matrix[left][right])

    def emit_atom(atom):
        atom_children = children[atom]
        continuation = preferred_next.get(atom)
        if continuation not in atom_children:
            continuation = atom_children[-1] if atom_children else None
        branches = [child for child in atom_children if child != continuation]
        parts = ["C"]
        for child in branches:
            symbol = bond_symbol(atom, child)
            child_text = emit_atom(child)
            if symbol is None or child_text is None:
                return None
            parts.append(f"({symbol}{child_text})")
        if continuation is not None:
            symbol = bond_symbol(atom, continuation)
            child_text = emit_atom(continuation)
            if symbol is None or child_text is None:
                return None
            parts.append(symbol + child_text)
        return "".join(parts)

    return emit_atom(root)


def mat2stereo_smiles(mat, assignment, analysis=None):
    """Convert a molecule and E/Z assignment to a stereo-aware SMILES-like string."""
    if not assignment.labels:
        return mat2smiles(mat)
    if analysis is None:
        analysis = stereochemistry.analyze_ez(mat)

    acyclic_smiles = _acyclic_stereo_smiles(mat, assignment, analysis)
    if acyclic_smiles is not None:
        return acyclic_smiles

    # Preserve as much standard slash stereochemistry as possible.  A label is
    # left in the extension only when no compatible standard rendering that
    # includes it can be produced by the current renderer.
    for encoded_count in range(len(assignment.labels) - 1, 0, -1):
        for encoded_labels in itertools.combinations(assignment.labels, encoded_count):
            partial_assignment = stereochemistry.EzAssignment(labels=encoded_labels)
            partial_smiles = _acyclic_stereo_smiles(
                mat,
                partial_assignment,
                analysis,
            )
            if partial_smiles is None:
                continue
            encoded = set(encoded_labels)
            unresolved_suffix = "".join(
                f" [{label}:{atom1 + 1}-{atom2 + 1}]"
                for atom1, atom2, label in assignment.labels
                if (atom1, atom2, label) not in encoded
            )
            return partial_smiles + unresolved_suffix

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
