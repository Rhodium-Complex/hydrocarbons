"""炭化水素の結合行列をSMILES文字列へ変換する。"""

import stereochemistry
from generation_output import CumuleneStereo, StructureVariant


BOND_SYMBOLS = {1: "", 2: "=", 3: "#"}


# 通常SMILESとE/Z方向制約
def _render_smiles(mat, directional_bonds=None):
    """結合行列を描画し、必要なら指定単結合へ方向記号を付ける。"""
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


def _stereo_data(mat, labels, analysis):
    """配置ラベルをスラッシュSMILESで表現可能な二重結合へ対応付ける。"""
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


def _solve_directional_bonds(stereo_data, edge_reversed):
    """複数のE/Z制約を同時に満たす相対的なスラッシュ方向を解く。"""
    adjacency = {}
    for double_bond, label in stereo_data:
        edge1 = frozenset((double_bond.high_ligand1, double_bond.atom1))
        edge2 = frozenset((double_bond.atom2, double_bond.high_ligand2))
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
    """最初のE/Z結合を起点とする走査で非環状の立体SMILESを描画する。"""
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


def _stereo_traversal(bonds, root):
    """環状立体SMILESの描画に使う決定的な深さ優先走査を構築する。"""
    parent = {root: None}
    children = {atom: [] for atom in range(len(bonds))}
    preorder = []

    def visit(atom):
        preorder.append(atom)
        for neighbor in sorted(int(value) for value in bonds[atom].nonzero()[0]):
            if neighbor == parent[atom] or neighbor in parent:
                continue
            parent[neighbor] = atom
            children[atom].append(neighbor)
            visit(neighbor)

    visit(root)
    if len(parent) != len(bonds):
        return None

    tree_edges = {
        frozenset((atom, parent_atom))
        for atom, parent_atom in parent.items()
        if parent_atom is not None
    }
    rank = {atom: index for index, atom in enumerate(preorder)}
    closure_edges = []
    for atom1, row in enumerate(bonds):
        for atom2, bond_order in enumerate(row[atom1 + 1 :], start=atom1 + 1):
            edge = frozenset((atom1, atom2))
            if bond_order and edge not in tree_edges:
                closure_edges.append(edge)
    return parent, children, rank, tree_edges, closure_edges


# 環を含むSMILES走査と原子中心立体表記
def _render_cyclic_stereo(mat, traversal, directional_bonds, atom_tokens=None):
    """DFS木と方向を考慮した環閉鎖辺を標準SMILESとして描画する。"""
    bonds = mat.bonds
    atom_tokens = atom_tokens or {}
    parent, children, rank, _tree_edges, closure_edges = traversal
    root = next(atom for atom, parent_atom in parent.items() if parent_atom is None)
    closures = {atom: [] for atom in range(len(bonds))}
    for ring_number, edge in enumerate(closure_edges, start=1):
        atom1, atom2 = sorted(edge, key=lambda atom: rank[atom])
        token = str(ring_number) if ring_number < 10 else f"%{ring_number}"
        symbol = directional_bonds.get(edge, BOND_SYMBOLS[bonds[atom1][atom2]])
        closures[atom1].append(symbol + token)
        closures[atom2].append(token)

    def emit(atom):
        result = [atom_tokens.get(atom, "C"), *closures[atom]]
        continuation = children[atom][-1] if children[atom] else None
        for child in children[atom]:
            edge = frozenset((atom, child))
            symbol = directional_bonds.get(edge, BOND_SYMBOLS[bonds[atom][child]])
            text = symbol + emit(child)
            result.append(text if child == continuation else f"({text})")
        return "".join(result)

    return emit(root)


def _edge_reversed(traversal, left, right, bonds=None):
    """走査上の辺がleftからrightとは逆向きに出力されるかを返す。"""
    parent, _children, rank, tree_edges, _closure_edges = traversal
    if bonds is not None and bonds[left][right] == 0:
        return None
    edge = frozenset((left, right))
    if edge in tree_edges:
        return parent.get(right) != left
    return not (rank[left] < rank[right])


def _lexical_neighbors(atom, traversal, hydrogen_count=0):
    """原子中心のSMILES立体表記で使われる字句順に隣接原子を返す。"""
    parent, children, _, _tree_edges, closure_edges = traversal
    ordered = [stereochemistry.HYDROGEN_LIGAND] * hydrogen_count
    if parent[atom] is not None:
        ordered.append(parent[atom])
    ring_neighbors = []
    for ring_number, edge in enumerate(closure_edges, start=1):
        if atom in edge:
            other = next(value for value in edge if value != atom)
            ring_neighbors.append((ring_number, other))
    ordered.extend(other for _number, other in sorted(ring_neighbors))
    ordered.extend(children[atom])
    return tuple(ordered)


def _chiral_atom_tokens(
    chiral_analysis,
    chiral_assignment,
    traversal,
    active_centers=None,
):
    """抽象的な配置ビットを走査順に対応した原子トークンへ変換する。"""
    labels = {(kind, atom): bit for kind, atom, bit in chiral_assignment.labels}
    active_centers = active_centers or chiral_assignment.active_centers
    tokens = {}
    for center in chiral_analysis.tetrahedral_centers:
        if ("T", center.atom) not in active_centers:
            continue
        bit = labels.get(("T", center.atom))
        if bit is None:
            continue
        lexical = _lexical_neighbors(
            center.atom,
            traversal,
            center.hydrogen_count,
        )
        if set(lexical) != set(center.ligands) or len(lexical) != 4:
            continue
        parity = stereochemistry._permutation_is_odd(center.ligands, lexical)
        marker = "@@" if bit ^ int(parity) else "@"
        hydrogen = "H" if center.hydrogen_count else ""
        tokens[center.atom] = f"[C{marker}{hydrogen}]"

    for center in chiral_analysis.allene_centers:
        if ("A", center.center_atom) not in active_centers:
            continue
        bit = labels.get(("A", center.center_atom))
        if bit is None:
            continue
        terminal1, terminal2 = center.path[0], center.path[-1]
        lexical1 = tuple(
            value
            for value in _lexical_neighbors(terminal1, traversal, int(
                stereochemistry.HYDROGEN_LIGAND in center.ligands1
            ))
            if value != center.path[1]
        )
        lexical2 = tuple(
            value
            for value in _lexical_neighbors(terminal2, traversal, int(
                stereochemistry.HYDROGEN_LIGAND in center.ligands2
            ))
            if value != center.path[-2]
        )
        if set(lexical1) != set(center.ligands1) or set(lexical2) != set(
            center.ligands2
        ):
            continue
        parity = stereochemistry._permutation_is_odd(center.ligands1, lexical1) ^ (
            stereochemistry._permutation_is_odd(center.ligands2, lexical2)
        )
        marker = "AL2" if bit ^ int(parity) else "AL1"
        tokens[center.center_atom] = f"[C@{marker}]"
    return tokens


def _combined_stereo_smiles(
    mat,
    chiral_assignment,
    chiral_analysis,
    ez_assignment=None,
    ez_analysis=None,
    traversal=None,
):
    """原子中心立体と表現可能なE/Z制約をまとめて描画する。"""
    traversal = traversal or _stereo_traversal(mat.bonds, 0)
    if traversal is None:
        return _render_smiles(mat)
    directional_bonds = {}
    if ez_assignment is not None and ez_assignment.labels:
        stereo_data = _stereo_data(mat, ez_assignment.labels, ez_analysis)
        if stereo_data:
            solved = _solve_directional_bonds(
                stereo_data,
                lambda left, right: _edge_reversed(traversal, left, right),
            )
            if solved is not None:
                directional_bonds = solved
    active_centers = (
        stereochemistry.active_chiral_centers(
            mat.bonds,
            chiral_analysis,
            chiral_assignment,
            ez_assignment,
        )
        if ez_assignment is not None and ez_assignment.labels
        else chiral_assignment.active_centers
    )
    atom_tokens = _chiral_atom_tokens(
        chiral_analysis,
        chiral_assignment,
        traversal,
        active_centers,
    )
    return _render_cyclic_stereo(
        mat,
        traversal,
        directional_bonds,
        atom_tokens,
    )


# 公開変換経路で使用するE/Z描画とフォールバック
def _cyclic_stereo_smiles(mat, assignment, analysis, traversal=None):
    """環閉鎖辺を含むアルケンの相対配置制約をすべて描画する。"""
    stereo_data = _stereo_data(mat, assignment.labels, analysis)
    if not stereo_data:
        return None

    root = stereo_data[0][0].high_ligand1
    traversal = traversal or _stereo_traversal(mat.bonds, root)
    if traversal is None:
        return None
    directional_bonds = _solve_directional_bonds(
        stereo_data,
        edge_reversed=lambda left, right: _edge_reversed(
            traversal, left, right, mat.bonds
        ),
    )
    if directional_bonds is None:
        return None
    return _render_cyclic_stereo(mat, traversal, directional_bonds)


def _standard_stereo_smiles(mat, assignment, analysis, cyclic_traversal=None):
    return _acyclic_stereo_smiles(mat, assignment, analysis) or _cyclic_stereo_smiles(
        mat,
        assignment,
        analysis,
        cyclic_traversal,
    )


def mat2stereo_smiles(mat, assignment, analysis=None):
    """アルケンの相対配置1件を標準的なスラッシュSMILESへ変換する。

    E/Zは二つの相対配置を区別する内部ラベルであり、正式なCIP命名を
    保証しない。制約を矛盾なく表現できない場合は立体情報なしのSMILESを返す。
    """
    if not assignment.labels:
        return _render_smiles(mat)
    analysis = analysis or stereochemistry.analyze_ez(mat)

    return _standard_stereo_smiles(mat, assignment, analysis) or _render_smiles(mat)


def mat2smiles_variants(
    mat,
    include_stereo: bool = False,
    include_tetrahedral_stereo: bool = False,
) -> list[str]:
    """1分子について指定された立体配置を含むSMILES候補を返す。"""
    if not include_stereo and not include_tetrahedral_stereo:
        return [_render_smiles(mat)]
    if include_tetrahedral_stereo:
        chiral_analysis = stereochemistry.analyze_chiral(mat)
        if not chiral_analysis.tetrahedral_centers and not chiral_analysis.allene_centers:
            return mat2smiles_variants(mat, include_stereo, False)
        ez_analysis = stereochemistry.analyze_ez(mat) if include_stereo else None
        ez_assignments = (
            ez_analysis.assignments
            if ez_analysis is not None
            else (stereochemistry.EzAssignment(labels=()),)
        )
        variants = []
        seen = set()
        traversal = _stereo_traversal(mat.bonds, 0)
        for ez_assignment in ez_assignments:
            chiral_assignments = stereochemistry.chiral_assignments_for_ez(
                mat.bonds,
                chiral_analysis,
                ez_assignment,
            )
            for chiral_assignment in chiral_assignments:
                smiles = _combined_stereo_smiles(
                    mat,
                    chiral_assignment,
                    chiral_analysis,
                    ez_assignment,
                    ez_analysis,
                    traversal,
                )
                if smiles in seen:
                    continue
                seen.add(smiles)
                variants.append(smiles)
        return variants
    analysis = stereochemistry.analyze_ez(mat)
    if not analysis.double_bonds:
        return [_render_smiles(mat)]

    bonds = mat.bonds
    is_cyclic = int((bonds > 0).sum() // 2) != len(bonds) - 1
    cyclic_traversal = None
    if is_cyclic:
        root = analysis.double_bonds[0].high_ligand1
        cyclic_traversal = _stereo_traversal(bonds, root)

    variants = []
    non_stereo_smiles = None
    for assignment in analysis.assignments:
        smiles = _standard_stereo_smiles(
            mat,
            assignment,
            analysis,
            cyclic_traversal,
        )
        if smiles is None:
            non_stereo_smiles = non_stereo_smiles or _render_smiles(mat)
            smiles = non_stereo_smiles
        variants.append(smiles)

    seen = set()
    for index, smiles in enumerate(variants):
        if smiles in seen:
            non_stereo_smiles = non_stereo_smiles or _render_smiles(mat)
            variants[index] = non_stereo_smiles
        else:
            seen.add(smiles)
    return variants


def mat2structure_variants(
    mat,
    include_stereo: bool = False,
    include_tetrahedral_stereo: bool = False,
):
    """Return directly drawable stereochemical structure variants."""
    base = mat2smiles_variants(mat, include_stereo, include_tetrahedral_stereo)
    bonds = tuple(tuple(int(value) for value in row) for row in mat.bonds)
    if not include_stereo:
        return [StructureVariant(value) for value in base]
    analysis = stereochemistry.analyze_cumulene_ez(mat)
    if not analysis.centers:
        return [StructureVariant(value) for value in base]
    center_by_edge = {
        (center.path[0], center.path[-1]): center for center in analysis.centers
    }
    result = []
    for value in base:
        for assignment in analysis.assignments:
            decorations = []
            for atom1, atom2, label in assignment.labels:
                center = center_by_edge[(atom1, atom2)]
                decorations.append(
                    CumuleneStereo(
                        center.path,
                        center.high_ligand1,
                        center.high_ligand2,
                        label,
                    )
                )
            result.append(StructureVariant(value, bonds, tuple(decorations)))
    return result
