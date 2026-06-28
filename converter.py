"""Convert hydrocarbon bond matrices into SMILES-like strings."""
import itertools

import stereochemistry


BOND_SYMBOLS = {1: "", 2: "=", 3: "#"}


def _mat2smiles(mat, directional_bonds=None):
    """Render a bond matrix, optionally marking selected single bonds."""
    bonds = mat.bonds
    directional_bonds = directional_bonds or {}
    parts = [""] * len(bonds)
    ring_placeholder = 1

    def bond_symbol(atom1, atom2, order):
        return directional_bonds.get(
            frozenset((atom1, atom2)),
            BOND_SYMBOLS[order],
        )

    for atom in range(len(bonds) - 1, -1, -1):
        has_parent = False
        parts[atom] = "C" + parts[atom]
        for neighbor in range(atom - 1, -1, -1):
            order = bonds[atom][neighbor]
            if order == 0:
                continue
            if has_parent:
                token = "{" + str(ring_placeholder) + "}"
                parts[atom] += bond_symbol(atom, neighbor, order) + token
                parts[neighbor] = token + parts[neighbor]
                ring_placeholder += 1
            else:
                parts[atom] = bond_symbol(neighbor, atom, order) + parts[atom]
            has_parent = True

        children = []
        for neighbor in range(len(bonds) - 1, atom - 1, -1):
            if bonds[atom][neighbor] != 0 and parts[neighbor]:
                children.append(parts[neighbor])
                parts[neighbor] = ""
        if children:
            branches = "".join(f"({child})" for child in reversed(children[1:]))
            parts[atom] += branches + children[0]

    smiles = parts[0]
    ring_number = 1
    while "{" in smiles:
        start = smiles.find("{")
        end = smiles.find("}", start)
        token = smiles[start : end + 1]
        replacement = str(ring_number) if ring_number < 10 else f"%{ring_number}"
        smiles = smiles.replace(token, replacement, 2)
        ring_number += 1
    return smiles


def mat2smiles(mat):
    """Convert a molecule bond matrix into a compact SMILES-like string."""
    return _mat2smiles(mat)


def _stereo_data(mat, labels, analysis):
    """Resolve assignment labels to double bonds usable by slash SMILES."""
    by_edge = {
        (double_bond.atom1, double_bond.atom2): double_bond
        for double_bond in analysis.double_bonds
    }
    resolved = []
    for atom1, atom2, label in labels:
        double_bond = by_edge.get((min(atom1, atom2), max(atom1, atom2)))
        if double_bond is None or not double_bond.has_carbon_high_ligands:
            return None
        if (
            mat.bonds[double_bond.high_ligand1][double_bond.atom1] != 1
            or mat.bonds[double_bond.atom2][double_bond.high_ligand2] != 1
        ):
            return None
        resolved.append((double_bond, label))
    return resolved


def _solve_directional_bonds(stereo_data, edge_reversed, allowed_edges=None):
    """Solve relative slash directions for a collection of E/Z constraints."""
    adjacency = {}
    for double_bond, label in stereo_data:
        edge1 = frozenset((double_bond.high_ligand1, double_bond.atom1))
        edge2 = frozenset((double_bond.atom2, double_bond.high_ligand2))
        if allowed_edges is not None and (
            edge1 not in allowed_edges or edge2 not in allowed_edges
        ):
            return None
        reversed1 = edge_reversed(double_bond.high_ligand1, double_bond.atom1)
        reversed2 = edge_reversed(double_bond.atom2, double_bond.high_ligand2)
        if reversed1 is None or reversed2 is None:
            return None
        differs = (label == "Z") ^ reversed1 ^ reversed2
        adjacency.setdefault(edge1, []).append((edge2, differs))
        adjacency.setdefault(edge2, []).append((edge1, differs))

    marks = {}
    for start in sorted(adjacency, key=lambda edge: tuple(sorted(edge))):
        if start in marks:
            continue
        marks[start] = False
        pending = [start]
        while pending:
            edge = pending.pop()
            for other, differs in adjacency[edge]:
                expected = marks[edge] ^ differs
                if other in marks:
                    if marks[other] != expected:
                        return None
                else:
                    marks[other] = expected
                    pending.append(other)
    return {edge: "\\" if mark else "/" for edge, mark in marks.items()}


def _acyclic_stereo_smiles(mat, assignment, analysis):
    """Render slash SMILES using a traversal rooted at the first E/Z bond."""
    bonds = mat.bonds
    atom_count = len(bonds)
    if int((bonds > 0).sum() // 2) != atom_count - 1:
        return None
    stereo_data = _stereo_data(mat, assignment.labels, analysis)
    if not stereo_data:
        return None

    preferred_next = {}
    for double_bond, _label in stereo_data:
        path = (
            double_bond.high_ligand1,
            double_bond.atom1,
            double_bond.atom2,
            double_bond.high_ligand2,
        )
        for left, right in zip(path, path[1:]):
            preferred_next.setdefault(left, right)

    root = stereo_data[0][0].high_ligand1
    parent = {root: None}
    children = {atom: [] for atom in range(atom_count)}
    pending = [root]
    while pending:
        atom = pending.pop()
        for neighbor in sorted(int(value) for value in bonds[atom].nonzero()[0]):
            if neighbor == parent[atom]:
                continue
            if neighbor in parent:
                return None
            parent[neighbor] = atom
            children[atom].append(neighbor)
            pending.append(neighbor)
    if len(parent) != atom_count:
        return None

    def edge_reversed(left, right):
        if parent.get(right) == left:
            return False
        if parent.get(left) == right:
            return True
        return None

    directional_bonds = _solve_directional_bonds(stereo_data, edge_reversed)
    if directional_bonds is None:
        return None

    def emit(atom):
        continuation = preferred_next.get(atom)
        if continuation not in children[atom]:
            continuation = children[atom][-1] if children[atom] else None
        ordered_children = [
            child for child in children[atom] if child != continuation
        ]
        if continuation is not None:
            ordered_children.append(continuation)

        result = ["C"]
        for child in ordered_children:
            symbol = directional_bonds.get(
                frozenset((atom, child)),
                BOND_SYMBOLS.get(bonds[atom][child]),
            )
            child_smiles = emit(child)
            if symbol is None or child_smiles is None:
                return None
            text = symbol + child_smiles
            result.append(text if child == continuation else f"({text})")
        return "".join(result)

    return emit(root)


def _base_tree_stereo_smiles(mat, assignment, analysis):
    """Mark E/Z ligand bonds used by the normal SMILES spanning tree."""
    stereo_data = _stereo_data(mat, assignment.labels, analysis)
    if not stereo_data:
        return None

    tree_edges = set()
    for atom in range(len(mat.bonds)):
        for neighbor in range(atom - 1, -1, -1):
            if mat.bonds[atom][neighbor] != 0:
                tree_edges.add(frozenset((neighbor, atom)))
                break

    directional_bonds = _solve_directional_bonds(
        stereo_data,
        edge_reversed=lambda left, right: left > right,
        allowed_edges=tree_edges,
    )
    if directional_bonds is None:
        return None
    return _mat2smiles(mat, directional_bonds)


def _standard_stereo_smiles(mat, assignment, analysis):
    return _acyclic_stereo_smiles(mat, assignment, analysis) or _base_tree_stereo_smiles(
        mat,
        assignment,
        analysis,
    )


def _assignment_suffix(labels):
    return "".join(
        f" [{label}:{atom1 + 1}-{atom2 + 1}]"
        for atom1, atom2, label in labels
    )


def mat2stereo_smiles(mat, assignment, analysis=None):
    """Convert a molecule and E/Z assignment to a stereo-aware SMILES-like string."""
    if not assignment.labels:
        return mat2smiles(mat)
    analysis = analysis or stereochemistry.analyze_ez(mat)

    labels = assignment.labels
    for encoded_count in range(len(labels), 0, -1):
        for encoded_labels in itertools.combinations(labels, encoded_count):
            partial = stereochemistry.EzAssignment(labels=encoded_labels)
            smiles = _standard_stereo_smiles(mat, partial, analysis)
            if smiles is None:
                continue
            encoded = set(encoded_labels)
            unresolved = tuple(label for label in labels if label not in encoded)
            return smiles + _assignment_suffix(unresolved)
    return mat2smiles(mat) + _assignment_suffix(labels)


def mat2smiles_variants(mat, include_stereo: bool = False) -> list[str]:
    """Return one or more SMILES outputs for a molecule."""
    if not include_stereo:
        return [mat2smiles(mat)]
    analysis = stereochemistry.analyze_ez(mat)
    return [
        mat2stereo_smiles(mat, assignment, analysis)
        for assignment in analysis.assignments
    ]
