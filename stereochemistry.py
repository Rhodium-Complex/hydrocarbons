"""炭化水素の結合行列から立体配置を列挙する補助処理。"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import itertools

import numpy as np

import isomorphism
import molecule


HYDROGEN_LIGAND = -1
EzLabel = tuple[int, int, str]
ChiralLabel = tuple[str, int, int]
EzAction = dict[tuple[int, int], tuple[tuple[int, int], bool]]


@dataclass(frozen=True)
class EzDoubleBond:
    """両端に区別可能な配位子を持つ二重結合。"""

    atom1: int
    atom2: int
    high_ligand1: int
    high_ligand2: int
    other_ligand1: int | None = None
    other_ligand2: int | None = None

    @property
    def has_carbon_high_ligands(self) -> bool:
        """両側の代表配位子が明示的な炭素原子かを返す。"""
        return (
            self.high_ligand1 != HYDROGEN_LIGAND
            and self.high_ligand2 != HYDROGEN_LIGAND
        )


@dataclass(frozen=True)
class EzAssignment:
    """立体配置を持つ二重結合へ割り当てたE/Zラベル群。"""

    labels: tuple[EzLabel, ...]

    @property
    def single_label(self) -> EzLabel | None:
        """割当て対象が1結合だけの場合に、そのE/Zラベルを返す。"""
        if len(self.labels) != 1:
            return None
        return self.labels[0]

    @property
    def single_edge(self) -> tuple[int, int] | None:
        """割当て対象が1結合だけの場合に、正規化した辺を返す。"""
        single_label = self.single_label
        if single_label is None:
            return None
        atom1, atom2, _label = single_label
        return (min(atom1, atom2), max(atom1, atom2))


@dataclass(frozen=True)
class EzAnalysis:
    """1分子のE/Z候補二重結合と対称性を除いた配置割当て。"""

    double_bonds: tuple[EzDoubleBond, ...]
    assignments: tuple[EzAssignment, ...]
    conditional_double_bonds: tuple[EzDoubleBond, ...] = ()

    @property
    def all_double_bonds(self) -> tuple[EzDoubleBond, ...]:
        """Return static and configuration-dependent E/Z candidates."""
        return self.double_bonds + self.conditional_double_bonds

    def double_bond_for_assignment(
        self,
        assignment: EzAssignment,
    ) -> EzDoubleBond | None:
        """単一のE/Z割当てが指す二重結合を返す。"""
        edge = assignment.single_edge
        if edge is None:
            return None
        for double_bond in self.all_double_bonds:
            if edge == (double_bond.atom1, double_bond.atom2):
                return double_bond
        return None


@dataclass(frozen=True)
class CumuleneEzCenter:
    """Odd cumulene chain with distinguishable ligand pairs at both ends."""

    path: tuple[int, ...]
    high_ligand1: int
    high_ligand2: int


@dataclass(frozen=True)
class CumuleneEzAnalysis:
    centers: tuple[CumuleneEzCenter, ...]
    assignments: tuple[EzAssignment, ...]


@dataclass(frozen=True)
class TetrahedralCenter:
    """四面体炭素と、決定的に定めた4配位子の基準順序。"""

    atom: int
    ligands: tuple[int, int, int, int]
    hydrogen_count: int


@dataclass(frozen=True)
class AlleneCenter:
    """両末端に区別可能な配位子対を持つ偶数クムレン鎖。"""

    center_atom: int
    path: tuple[int, ...]
    ligands1: tuple[int, int]
    ligands2: tuple[int, int]


@dataclass(frozen=True)
class ChiralAssignment:
    """四面体中心とallene-like中心の二値配置割当て。"""

    labels: tuple[ChiralLabel, ...]
    active_centers: frozenset[tuple[str, int]] = frozenset()


@dataclass(frozen=True)
class ChiralAnalysis:
    """原子中心の立体要素と、対称性を除いた配置割当て。"""

    tetrahedral_centers: tuple[TetrahedralCenter, ...]
    allene_centers: tuple[AlleneCenter, ...]
    assignments: tuple[ChiralAssignment, ...]
    automorphisms: tuple[tuple[int, ...], ...] = ()


# 原子・配位子の同一性判定
def _compressed_color_ids(signatures: list[tuple]) -> list[int]:
    """シグネチャへ決定的で連続した整数の色番号を割り当てる。"""
    color_by_signature = {
        signature: index + 1
        for index, signature in enumerate(sorted(set(signatures)))
    }
    return [color_by_signature[signature] for signature in signatures]


def _refined_atom_colors(blocked: np.ndarray) -> list[int]:
    """Weisfeiler-Lehman型の反復細分化により原子色を求める。"""
    hydrogens = molecule.implicit_hydrogens(blocked)
    neighbors = [
        tuple(int(neighbor) for neighbor in np.where(blocked[index] > 0)[0])
        for index in range(len(blocked))
    ]
    colors = _compressed_color_ids(
        [
            (
                int(hydrogens[index]),
                len(neighbors[index]),
                int(np.sum(blocked[index])),
            )
            for index in range(len(blocked))
        ]
    )

    for _ in range(max(1, len(blocked) * 2)):
        next_colors = _compressed_color_ids(
            [
                (
                    colors[index],
                    tuple(
                        sorted(
                            (
                                int(blocked[index][neighbor]),
                                colors[neighbor],
                            )
                            for neighbor in neighbors[index]
                        )
                    ),
                )
                for index in range(len(blocked))
            ]
        )
        if next_colors == colors:
            break
        colors = next_colors
    return colors


def _ligands_for_atom(
    bonds: np.ndarray,
    atom: int,
    double_bond_partner: int,
    hydrogens: np.ndarray,
) -> list[int]:
    ligands = [
        int(neighbor)
        for neighbor in np.where(bonds[atom] > 0)[0]
        if int(neighbor) != double_bond_partner
    ]
    ligands.extend([HYDROGEN_LIGAND] * int(hydrogens[atom]))
    return ligands


def _high_priority_ligand(
    colors: list[int],
    ligands: list[int],
) -> int | None:
    if len(ligands) < 2:
        return None
    ranked = sorted(
        (
            0 if ligand == HYDROGEN_LIGAND else colors[ligand] + 1,
            ligand,
        )
        for ligand in ligands
    )
    if ranked[-1][0] == ranked[-2][0]:
        return None
    return ranked[-1][1]


def _locally_resolved_high_ligand(ligands: list[int]) -> tuple[bool, int | None]:
    """Resolve ligand identity without graph refinement when it is unambiguous.

    Hydrocarbon alkene ends with one explicit carbon and one implicit hydrogen
    always select the carbon.  Missing or all-hydrogen ligand sets can also be
    rejected locally.  Two explicit carbon ligands still require graph-based
    equivalence refinement.
    """
    if len(ligands) < 2:
        return True, None
    carbon_ligands = [
        ligand for ligand in ligands if ligand != HYDROGEN_LIGAND
    ]
    if not carbon_ligands:
        return True, None
    if len(ligands) == 2 and len(carbon_ligands) == 1:
        return True, carbon_ligands[0]
    return False, None


def _other_ligand(ligands: list[int], high_ligand: int) -> int | None:
    """Return the other ligand in a two-ligand alkene end."""
    for ligand in ligands:
        if ligand != high_ligand:
            return ligand
    return None


# E/Z候補の検出と対称性を考慮した配置列挙
def _find_ez_double_bonds_for_bonds(
    bonds: np.ndarray,
) -> tuple[tuple[EzDoubleBond, ...], tuple[EzDoubleBond, ...]]:
    hydrogens = molecule.implicit_hydrogens(bonds)
    cumulene_edges = {
        frozenset((left, right))
        for path in _double_bond_paths(bonds)
        if len(path) > 2
        for left, right in zip(path, path[1:])
    }
    double_bonds = []
    conditional_double_bonds = []
    for atom1, bond_row in enumerate(bonds):
        for atom2, bond_order in enumerate(bond_row[atom1 + 1 :], start=atom1 + 1):
            if bond_order != 2:
                continue
            if frozenset((atom1, atom2)) in cumulene_edges:
                continue

            ligands1 = _ligands_for_atom(bonds, atom1, atom2, hydrogens)
            ligands2 = _ligands_for_atom(bonds, atom2, atom1, hydrogens)
            resolved1, high1 = _locally_resolved_high_ligand(ligands1)
            resolved2, high2 = _locally_resolved_high_ligand(ligands2)
            if (resolved1 and high1 is None) or (resolved2 and high2 is None):
                continue

            if resolved1 and resolved2:
                assert high1 is not None and high2 is not None
                double_bonds.append(
                    EzDoubleBond(
                        atom1=atom1,
                        atom2=atom2,
                        high_ligand1=high1,
                        high_ligand2=high2,
                        other_ligand1=_other_ligand(ligands1, high1),
                        other_ligand2=_other_ligand(ligands2, high2),
                    )
                )
                continue

            blocked = bonds.copy()
            blocked[atom1][atom2] = 0
            blocked[atom2][atom1] = 0
            colors = _refined_atom_colors(blocked)
            if not resolved1:
                high1 = _high_priority_ligand(colors, ligands1)
            if not resolved2:
                high2 = _high_priority_ligand(colors, ligands2)
            if high1 is not None and high2 is not None:
                double_bonds.append(
                    EzDoubleBond(
                        atom1=atom1,
                        atom2=atom2,
                        high_ligand1=high1,
                        high_ligand2=high2,
                        other_ligand1=_other_ligand(ligands1, high1),
                        other_ligand2=_other_ligand(ligands2, high2),
                    )
                )
                continue

            # A same-colour pair may become distinguishable after another
            # E/Z assignment.  Keep a deterministic carbon reference for
            # slash-SMILES rendering, but do not make it a static candidate.
            carbon_ligands1 = [value for value in ligands1 if value != HYDROGEN_LIGAND]
            carbon_ligands2 = [value for value in ligands2 if value != HYDROGEN_LIGAND]
            if not carbon_ligands1 or not carbon_ligands2:
                continue
            reference1 = high1 if high1 is not None else min(carbon_ligands1)
            reference2 = high2 if high2 is not None else min(carbon_ligands2)
            conditional_double_bonds.append(
                EzDoubleBond(
                    atom1=atom1,
                    atom2=atom2,
                    high_ligand1=reference1,
                    high_ligand2=reference2,
                    other_ligand1=_other_ligand(ligands1, reference1),
                    other_ligand2=_other_ligand(ligands2, reference2),
                )
            )
    return tuple(double_bonds), tuple(conditional_double_bonds)


def _enumerate_assignments(
    bonds: np.ndarray,
    double_bonds: tuple[EzDoubleBond, ...],
    automorphisms: tuple[tuple[int, ...], ...] | None = None,
) -> tuple[EzAssignment, ...]:
    if not double_bonds:
        return (EzAssignment(labels=()),)
    if len(double_bonds) == 1:
        double_bond = double_bonds[0]
        return tuple(
            EzAssignment(
                labels=((double_bond.atom1, double_bond.atom2, label),)
            )
            for label in ("E", "Z")
        )

    if automorphisms is None:
        colors = _refined_atom_colors(bonds)
        automorphisms = tuple(
            [tuple(range(len(bonds)))]
            if len(set(colors)) == len(colors)
            else isomorphism.automorphisms(bonds)
        )
    by_edge = _ez_bonds_by_edge(double_bonds)
    actions = tuple(
        _compile_ez_action(by_edge, automorphism)
        for automorphism in automorphisms
    )

    assignments = []
    seen_keys = set()
    for labels in itertools.product(("E", "Z"), repeat=len(double_bonds)):
        assignment = EzAssignment(
            labels=tuple(
                (
                    double_bond.atom1,
                    double_bond.atom2,
                    label,
                )
                for double_bond, label in zip(double_bonds, labels)
            )
        )
        canonical_key = min(
            _map_ez_labels(assignment.labels, action)
            for action in actions
        )
        if canonical_key in seen_keys:
            continue
        seen_keys.add(canonical_key)
        assignments.append(assignment)
    return tuple(assignments)


def _ez_bonds_by_edge(
    double_bonds: tuple[EzDoubleBond, ...],
) -> dict[tuple[int, int], EzDoubleBond]:
    return {
        (double_bond.atom1, double_bond.atom2): double_bond
        for double_bond in double_bonds
    }


def _compile_ez_action(
    by_edge: dict[tuple[int, int], EzDoubleBond],
    automorphism: tuple[int, ...],
) -> EzAction:
    """Compile one automorphism into E/Z edge moves and bit reversals."""
    action = {}
    for edge, source in by_edge.items():
        mapped_atom1 = automorphism[source.atom1]
        mapped_atom2 = automorphism[source.atom2]
        mapped_edge = (min(mapped_atom1, mapped_atom2), max(mapped_atom1, mapped_atom2))
        target = by_edge[mapped_edge]
        target_ref1 = target.high_ligand1 if target.atom1 == mapped_atom1 else target.high_ligand2
        target_ref2 = target.high_ligand1 if target.atom1 == mapped_atom2 else target.high_ligand2
        action[edge] = (
            mapped_edge,
            (automorphism[source.high_ligand1] != target_ref1)
            ^ (automorphism[source.high_ligand2] != target_ref2),
        )
    return action


def _map_ez_labels(labels: tuple[EzLabel, ...], action: EzAction) -> tuple[EzLabel, ...]:
    """Map relative E/Z labels with a precompiled automorphism action."""
    mapped_labels = []
    for atom1, atom2, label in labels:
        mapped_edge, flipped = action[(min(atom1, atom2), max(atom1, atom2))]
        mapped_label = "Z" if flipped and label == "E" else "E" if flipped else label
        mapped_labels.append(
            (*mapped_edge, mapped_label)
        )
    return tuple(sorted(mapped_labels))


def _ez_action_for_labels(
    labels: tuple[EzLabel, ...],
    automorphism: tuple[int, ...],
    by_edge: dict[tuple[int, int], EzDoubleBond] | None,
) -> EzAction:
    if by_edge is not None:
        return _compile_ez_action(by_edge, automorphism)
    return {
        (min(atom1, atom2), max(atom1, atom2)): (
            (min(automorphism[atom1], automorphism[atom2]), max(automorphism[atom1], automorphism[atom2])),
            False,
        )
        for atom1, atom2, _label in labels
    }


def _conditional_ez_is_active(
    double_bond: EzDoubleBond,
    labels: tuple[EzLabel, ...],
    automorphisms: tuple[tuple[int, ...], ...],
    actions: tuple[EzAction, ...],
) -> bool:
    """Return whether existing stereo labels distinguish both alkene ends."""
    raw_labels = tuple(sorted(labels))

    def side_is_distinct(atom: int, reference: int, other: int | None) -> bool:
        if other is None:
            return True
        for automorphism, action in zip(automorphisms, actions):
            if automorphism[atom] != atom or automorphism[reference] != other:
                continue
            if _map_ez_labels(raw_labels, action) == raw_labels:
                return False
        return True

    return side_is_distinct(
        double_bond.atom1, double_bond.high_ligand1, double_bond.other_ligand1
    ) and side_is_distinct(
        double_bond.atom2, double_bond.high_ligand2, double_bond.other_ligand2
    )


def _expand_conditional_ez_assignments(
    static_assignments: tuple[EzAssignment, ...],
    conditional_double_bonds: tuple[EzDoubleBond, ...],
    automorphisms: tuple[tuple[int, ...], ...],
    all_double_bonds: tuple[EzDoubleBond, ...],
) -> tuple[EzAssignment, ...]:
    by_edge = _ez_bonds_by_edge(all_double_bonds)
    actions = tuple(
        _compile_ez_action(by_edge, automorphism)
        for automorphism in automorphisms
    )
    completed = []
    seen_states = set()

    def canonical(labels: tuple[EzLabel, ...]) -> tuple[EzLabel, ...]:
        return min(_map_ez_labels(labels, action) for action in actions)

    def visit(labels: tuple[EzLabel, ...]) -> None:
        labels = canonical(labels)
        if labels in seen_states:
            return
        seen_states.add(labels)
        assigned_edges = {(atom1, atom2) for atom1, atom2, _label in labels}
        enabled = [
            double_bond
            for double_bond in conditional_double_bonds
            if (double_bond.atom1, double_bond.atom2) not in assigned_edges
            and _conditional_ez_is_active(double_bond, labels, automorphisms, actions)
        ]
        if not enabled:
            completed.append(EzAssignment(labels))
            return
        double_bond = enabled[0]
        for label in ("E", "Z"):
            visit(tuple(sorted((*labels, (double_bond.atom1, double_bond.atom2, label)))))

    for assignment in static_assignments:
        visit(assignment.labels)
    return tuple(completed)


def analyze_ez(
    molecule_obj,
    automorphisms: tuple[tuple[int, ...], ...] | None = None,
) -> EzAnalysis:
    """1分子についてE/Z候補と配置割当てを解析する。"""
    double_bonds, conditional_double_bonds = _find_ez_double_bonds_for_bonds(
        molecule_obj.bonds
    )
    if conditional_double_bonds and double_bonds:
        all_double_bonds = double_bonds + conditional_double_bonds
        automorphisms = automorphisms or tuple(isomorphism.automorphisms(molecule_obj.bonds))
        assignments = _expand_conditional_ez_assignments(
            _enumerate_assignments(molecule_obj.bonds, double_bonds, automorphisms),
            conditional_double_bonds,
            automorphisms,
            all_double_bonds,
        )
    else:
        # No prior E/Z decoration can activate a conditional center.  Retain
        # the previous WL fast path for the overwhelmingly common case.
        assignments = _enumerate_assignments(molecule_obj.bonds, double_bonds)
    return EzAnalysis(
        double_bonds=double_bonds,
        assignments=assignments,
        conditional_double_bonds=conditional_double_bonds,
    )


def analyze_cumulene_ez(molecule_obj) -> CumuleneEzAnalysis:
    """Return symmetry-unique E/Z assignments for odd cumulene chains."""
    bonds = molecule_obj.bonds
    hydrogens = molecule.implicit_hydrogens(bonds)
    centers = []
    for raw_path in _double_bond_paths(bonds):
        double_bond_count = len(raw_path) - 1
        if double_bond_count < 3 or double_bond_count % 2 == 0:
            continue
        path = raw_path if raw_path[0] < raw_path[-1] else tuple(reversed(raw_path))
        ligands1 = _terminal_ligands(bonds, path[0], path[1], hydrogens)
        ligands2 = _terminal_ligands(bonds, path[-1], path[-2], hydrogens)
        if len(ligands1) != 2 or len(ligands2) != 2:
            continue
        blocked = bonds.copy()
        for left, right in zip(path, path[1:]):
            blocked[left][right] = blocked[right][left] = 0
        colors = _refined_atom_colors(blocked)
        high1 = _high_priority_ligand(colors, ligands1)
        high2 = _high_priority_ligand(colors, ligands2)
        if high1 is None or high2 is None:
            continue
        centers.append(CumuleneEzCenter(path, high1, high2))

    pseudo_bonds = tuple(
        EzDoubleBond(center.path[0], center.path[-1], center.high_ligand1, center.high_ligand2)
        for center in centers
    )
    return CumuleneEzAnalysis(
        centers=tuple(centers),
        assignments=_enumerate_assignments(bonds, pseudo_bonds),
    )


# 四面体中心・allene-like中心の検出
def _ligand_order(colors: list[int], ligands: list[int]) -> tuple[int, ...] | None:
    signatures = [
        (0, ligand) if ligand == HYDROGEN_LIGAND else (colors[ligand] + 1, ligand)
        for ligand in ligands
    ]
    if len({signature[0] for signature in signatures}) != len(signatures):
        return None
    return tuple(ligand for _color, ligand in sorted(signatures))


def _find_tetrahedral_centers(bonds: np.ndarray) -> tuple[TetrahedralCenter, ...]:
    hydrogens = molecule.implicit_hydrogens(bonds)
    centers = []
    for atom in range(len(bonds)):
        neighbors = [int(value) for value in np.where(bonds[atom] > 0)[0]]
        hydrogen_count = int(hydrogens[atom])
        if hydrogen_count > 1 or len(neighbors) + hydrogen_count != 4:
            continue
        if any(bonds[atom][neighbor] != 1 for neighbor in neighbors):
            continue
        ligands = neighbors + [HYDROGEN_LIGAND] * hydrogen_count
        ordered = tuple(
            sorted(
                ligands,
                key=lambda ligand: (
                    ligand != HYDROGEN_LIGAND,
                    ligand,
                ),
            )
        )
        centers.append(
            TetrahedralCenter(atom, ordered, hydrogen_count)  # type: ignore[arg-type]
        )
    return tuple(centers)


def _double_bond_paths(bonds: np.ndarray) -> list[tuple[int, ...]]:
    adjacency = {
        atom: [int(n) for n in np.where(bonds[atom] == 2)[0]]
        for atom in range(len(bonds))
    }
    endpoints = sorted(atom for atom, values in adjacency.items() if len(values) == 1)
    visited_edges = set()
    paths = []
    for start in endpoints:
        first_edge = frozenset((start, adjacency[start][0]))
        if first_edge in visited_edges:
            continue
        path = [start]
        previous = None
        current = start
        while True:
            candidates = [n for n in adjacency[current] if n != previous]
            if not candidates:
                break
            next_atom = candidates[0]
            edge = frozenset((current, next_atom))
            if edge in visited_edges:
                break
            visited_edges.add(edge)
            path.append(next_atom)
            previous, current = current, next_atom
            if len(adjacency[current]) != 2:
                break
        if len(path) > 2:
            paths.append(tuple(path))
    return paths


def _terminal_ligands(
    bonds: np.ndarray,
    terminal: int,
    chain_neighbor: int,
    hydrogens: np.ndarray,
) -> list[int]:
    ligands = [
        int(value)
        for value in np.where(bonds[terminal] > 0)[0]
        if int(value) != chain_neighbor
    ]
    ligands.extend([HYDROGEN_LIGAND] * int(hydrogens[terminal]))
    return ligands


def _find_allene_centers(bonds: np.ndarray) -> tuple[AlleneCenter, ...]:
    hydrogens = molecule.implicit_hydrogens(bonds)
    centers = []
    for raw_path in _double_bond_paths(bonds):
        if (len(raw_path) - 1) % 2 or len(raw_path) < 3:
            continue
        path = raw_path if raw_path[0] < raw_path[-1] else tuple(reversed(raw_path))
        ligands1 = _terminal_ligands(bonds, path[0], path[1], hydrogens)
        ligands2 = _terminal_ligands(bonds, path[-1], path[-2], hydrogens)
        if len(ligands1) != 2 or len(ligands2) != 2:
            continue
        blocked = bonds.copy()
        for left, right in zip(path, path[1:]):
            blocked[left][right] = blocked[right][left] = 0
        colors = _refined_atom_colors(blocked)
        ordered1 = _ligand_order(colors, ligands1)
        ordered2 = _ligand_order(colors, ligands2)
        if ordered1 is None or ordered2 is None:
            continue
        centers.append(
            AlleneCenter(
                center_atom=path[len(path) // 2],
                path=path,
                ligands1=ordered1,  # type: ignore[arg-type]
                ligands2=ordered2,  # type: ignore[arg-type]
            )
        )
    return tuple(centers)


@lru_cache(maxsize=65536)
def _permutation_is_odd(source: tuple[int, ...], target: tuple[int, ...]) -> bool:
    permutation = [target.index(value) for value in source]
    inversions = 0
    for left, left_value in enumerate(permutation):
        for right_value in permutation[left + 1 :]:
            inversions += left_value > right_value
    return bool(inversions % 2)


def _map_ligand(ligand: int, automorphism: tuple[int, ...]) -> int:
    return ligand if ligand == HYDROGEN_LIGAND else automorphism[ligand]


# 自己同型写像による原子中心配置の変換と列挙
def _compile_chiral_action(elements, tetra_by_atom, allene_by_center, automorphism):
    """1つの自己同型写像を移動先中心と配置反転の有無へ変換する。"""
    action = []
    for kind, atom in elements:
        mapped_atom = automorphism[atom]
        if kind == "T":
            source = tetra_by_atom[atom]
            target = tetra_by_atom[mapped_atom]
            mapped_ligands = tuple(
                _map_ligand(value, automorphism) for value in source.ligands
            )
            parity = _permutation_is_odd(mapped_ligands, target.ligands)
        else:
            source = allene_by_center[atom]
            target = allene_by_center[mapped_atom]
            mapped_end = automorphism[source.path[0]]
            target_pairs = (
                (target.ligands1, target.ligands2)
                if mapped_end == target.path[0]
                else (target.ligands2, target.ligands1)
            )
            mapped_pairs = (
                tuple(_map_ligand(value, automorphism) for value in source.ligands1),
                tuple(_map_ligand(value, automorphism) for value in source.ligands2),
            )
            parity = _permutation_is_odd(mapped_pairs[0], target_pairs[0]) ^ (
                _permutation_is_odd(mapped_pairs[1], target_pairs[1])
            )
        action.append((kind, mapped_atom, int(parity)))
    return tuple(action)


def _map_chiral_labels(labels, action):
    return tuple(
        sorted(
            (kind, mapped_atom, bit ^ parity)
            for (_kind, _atom, bit), (kind, mapped_atom, parity) in zip(labels, action)
        )
    )


def _enumerate_chiral_assignments(
    bonds: np.ndarray,
    tetrahedral_centers: tuple[TetrahedralCenter, ...],
    allene_centers: tuple[AlleneCenter, ...],
    automorphisms: tuple[tuple[int, ...], ...] | None = None,
) -> tuple[ChiralAssignment, ...]:
    elements = [*(('T', center.atom) for center in tetrahedral_centers), *(
        ('A', center.center_atom) for center in allene_centers
    )]
    if not elements:
        return (ChiralAssignment(labels=()),)
    tetra_by_atom = {center.atom: center for center in tetrahedral_centers}
    allene_by_center = {center.center_atom: center for center in allene_centers}
    automorphisms = automorphisms or tuple(isomorphism.automorphisms(bonds))

    actions = tuple(
        _compile_chiral_action(
            elements, tetra_by_atom, allene_by_center, automorphism
        )
        for automorphism in automorphisms
    )

    def mapped_key(bits, action):
        return tuple(
            sorted(
                (kind, mapped_atom, bit ^ parity)
                for bit, (kind, mapped_atom, parity) in zip(bits, action)
            )
        )

    def unique_bit_patterns():
        seen_by_depth = [set() for _ in range(len(elements) + 1)]

        def visit(bits):
            depth = len(bits)
            if depth == len(elements):
                yield bits
                return
            prefix_ids = set(elements[: depth + 1])
            for bit in (0, 1):
                next_bits = bits + (bit,)
                keys = []
                for action in actions:
                    key = mapped_key(next_bits, action)
                    if {label[:2] for label in key} == prefix_ids:
                        keys.append(key)
                canonical_prefix = min(keys)
                if canonical_prefix in seen_by_depth[depth + 1]:
                    continue
                seen_by_depth[depth + 1].add(canonical_prefix)
                yield from visit(next_bits)

        yield from visit(())

    assignments = []
    seen = set()
    for bits in unique_bit_patterns():
        keys = [mapped_key(bits, action) for action in actions]
        canonical = min(keys)
        if canonical in seen:
            continue
        seen.add(canonical)
        raw_key = tuple(
            sorted(
                (kind, atom, bit)
                for (kind, atom), bit in zip(elements, bits)
            )
        )
        active_centers = {
            ("A", center.center_atom) for center in allene_centers
        }
        for center in tetrahedral_centers:
            other_labels = tuple(
                label
                for label in raw_key
                if label[:2] != ("T", center.atom)
            )
            center_index = elements.index(("T", center.atom))
            has_odd_stabilizer = any(
                action[center_index][1] == center.atom
                and bool(action[center_index][2])
                and tuple(
                    label
                    for label in mapped_key(bits, action)
                    if label[:2] != ("T", center.atom)
                ) == other_labels
                for action in actions
            )
            if not has_odd_stabilizer:
                active_centers.add(("T", center.atom))
        assignments.append(
            ChiralAssignment(
                labels=tuple(
                    (kind, atom, bit) for (kind, atom), bit in zip(elements, bits)
                ),
                active_centers=frozenset(active_centers),
            )
        )
    return tuple(assignments)


def analyze_chiral(molecule_obj) -> ChiralAnalysis:
    """四面体・allene-like立体要素と配置割当てを解析する。"""
    bonds = molecule_obj.bonds
    tetrahedral = _find_tetrahedral_centers(bonds)
    allenes = _find_allene_centers(bonds)
    automorphisms = (
        tuple(isomorphism.automorphisms(bonds))
        if tetrahedral or allenes
        else ()
    )
    return ChiralAnalysis(
        tetrahedral_centers=tetrahedral,
        allene_centers=allenes,
        assignments=_enumerate_chiral_assignments(
            bonds,
            tetrahedral,
            allenes,
            automorphisms,
        ),
        automorphisms=automorphisms,
    )


def active_chiral_centers(
    bonds: np.ndarray,
    analysis: ChiralAnalysis,
    assignment: ChiralAssignment,
    ez_assignment: EzAssignment | None = None,
    ez_analysis: EzAnalysis | None = None,
) -> frozenset[tuple[str, int]]:
    """選択した全立体配置を考慮した後に有効となる中心を返す。"""
    tetra_by_atom = {center.atom: center for center in analysis.tetrahedral_centers}
    allene_by_center = {center.center_atom: center for center in analysis.allene_centers}
    raw_labels = tuple(sorted(assignment.labels))
    raw_ez = tuple(sorted(ez_assignment.labels)) if ez_assignment is not None else ()
    ez_bonds = (
        _ez_bonds_by_edge(ez_analysis.all_double_bonds)
        if ez_analysis is not None
        else None
    )

    automorphisms = analysis.automorphisms or tuple(isomorphism.automorphisms(bonds))
    ez_actions = tuple(
        _ez_action_for_labels(raw_ez, automorphism, ez_bonds)
        for automorphism in automorphisms
    )
    elements = tuple(label[:2] for label in assignment.labels)
    active = {("A", center.center_atom) for center in analysis.allene_centers}
    for center in analysis.tetrahedral_centers:
        other_labels = tuple(
            label for label in raw_labels if label[:2] != ("T", center.atom)
        )
        has_odd_stabilizer = False
        for automorphism, ez_action in zip(automorphisms, ez_actions):
            if automorphism[center.atom] != center.atom:
                continue
            action = _compile_chiral_action(
                elements, tetra_by_atom, allene_by_center, automorphism
            )
            center_action = action[elements.index(("T", center.atom))]
            if not center_action[2]:
                continue
            mapped_other = tuple(
                label
                for label in _map_chiral_labels(assignment.labels, action)
                if label[:2] != ("T", center.atom)
            )
            if (
                mapped_other == other_labels
                and _map_ez_labels(raw_ez, ez_action) == raw_ez
            ):
                has_odd_stabilizer = True
                break
        if not has_odd_stabilizer:
            active.add(("T", center.atom))
    return frozenset(active)


def chiral_assignments_for_ez(
    bonds: np.ndarray,
    analysis: ChiralAnalysis,
    ez_assignment: EzAssignment,
    ez_analysis: EzAnalysis | None = None,
) -> tuple[ChiralAssignment, ...]:
    """E/Zラベルを保存する自己同型写像の下で原子中心配置を列挙する。"""
    if not ez_assignment.labels:
        return analysis.assignments
    raw_ez = tuple(sorted(ez_assignment.labels))
    ez_bonds = (
        _ez_bonds_by_edge(ez_analysis.all_double_bonds)
        if ez_analysis is not None
        else None
    )
    ez_actions = tuple(
        _ez_action_for_labels(raw_ez, automorphism, ez_bonds)
        for automorphism in analysis.automorphisms
    )

    preserving = tuple(
        automorphism
        for automorphism, ez_action in zip(analysis.automorphisms, ez_actions)
        if _map_ez_labels(raw_ez, ez_action) == raw_ez
    )
    return _enumerate_chiral_assignments(
        bonds,
        analysis.tetrahedral_centers,
        analysis.allene_centers,
        preserving,
    )
