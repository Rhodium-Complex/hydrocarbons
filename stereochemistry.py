"""Late-stage E/Z stereochemistry helpers for hydrocarbon bond matrices."""
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


@dataclass(frozen=True)
class EzDoubleBond:
    """A double bond whose two ends have distinguishable ligands."""

    atom1: int
    atom2: int
    high_ligand1: int
    high_ligand2: int

    @property
    def has_carbon_high_ligands(self) -> bool:
        """Return whether both high-priority ligands are explicit carbons."""
        return (
            self.high_ligand1 != HYDROGEN_LIGAND
            and self.high_ligand2 != HYDROGEN_LIGAND
        )


@dataclass(frozen=True)
class EzAssignment:
    """E/Z labels assigned to stereogenic double bonds."""

    labels: tuple[EzLabel, ...]

    @property
    def single_label(self) -> EzLabel | None:
        """Return the only E/Z label when this assignment targets one bond."""
        if len(self.labels) != 1:
            return None
        return self.labels[0]

    @property
    def single_edge(self) -> tuple[int, int] | None:
        """Return the normalized edge for a single-bond assignment."""
        single_label = self.single_label
        if single_label is None:
            return None
        atom1, atom2, _label = single_label
        return (min(atom1, atom2), max(atom1, atom2))


@dataclass(frozen=True)
class EzAnalysis:
    """E/Z-capable double bonds and unique assignments for one molecule."""

    double_bonds: tuple[EzDoubleBond, ...]
    assignments: tuple[EzAssignment, ...]

    def double_bond_for_assignment(
        self,
        assignment: EzAssignment,
    ) -> EzDoubleBond | None:
        """Return the double bond described by a single-bond E/Z assignment."""
        edge = assignment.single_edge
        if edge is None:
            return None
        for double_bond in self.double_bonds:
            if edge == (double_bond.atom1, double_bond.atom2):
                return double_bond
        return None


@dataclass(frozen=True)
class TetrahedralCenter:
    """A tetrahedral carbon and its deterministic four-ligand reference order."""

    atom: int
    ligands: tuple[int, int, int, int]
    hydrogen_count: int


@dataclass(frozen=True)
class AlleneCenter:
    """An even cumulene chain with distinguishable ligand pairs at both ends."""

    center_atom: int
    path: tuple[int, ...]
    ligands1: tuple[int, int]
    ligands2: tuple[int, int]


@dataclass(frozen=True)
class ChiralAssignment:
    """Binary configurations for tetrahedral and allene-like centers."""

    labels: tuple[ChiralLabel, ...]
    active_centers: frozenset[tuple[str, int]] = frozenset()


@dataclass(frozen=True)
class ChiralAnalysis:
    """Atom-centered stereogenic elements and symmetry-unique assignments."""

    tetrahedral_centers: tuple[TetrahedralCenter, ...]
    allene_centers: tuple[AlleneCenter, ...]
    assignments: tuple[ChiralAssignment, ...]
    automorphisms: tuple[tuple[int, ...], ...] = ()


def _compressed_color_ids(signatures: list[tuple]) -> list[int]:
    """Return deterministic dense integer color ids for signatures."""
    color_by_signature = {
        signature: index + 1
        for index, signature in enumerate(sorted(set(signatures)))
    }
    return [color_by_signature[signature] for signature in signatures]


def _refined_atom_colors(blocked: np.ndarray) -> list[int]:
    """Return integer atom colors from Weisfeiler-Lehman style refinement."""
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


def _find_ez_double_bonds_for_bonds(bonds: np.ndarray) -> tuple[EzDoubleBond, ...]:
    hydrogens = molecule.implicit_hydrogens(bonds)
    double_bonds = []
    for atom1, bond_row in enumerate(bonds):
        for atom2, bond_order in enumerate(bond_row[atom1 + 1 :], start=atom1 + 1):
            if bond_order != 2:
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
            if high1 is None or high2 is None:
                continue

            double_bonds.append(
                EzDoubleBond(
                    atom1=atom1,
                    atom2=atom2,
                    high_ligand1=high1,
                    high_ligand2=high2,
                )
            )
    return tuple(double_bonds)


def _enumerate_assignments(
    bonds: np.ndarray,
    double_bonds: tuple[EzDoubleBond, ...],
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

    colors = _refined_atom_colors(bonds)
    automorphisms = (
        [tuple(range(len(bonds)))]
        if len(set(colors)) == len(colors)
        else isomorphism.automorphisms(bonds)
    )

    def assignment_key(
        assignment: EzAssignment,
        automorphism: tuple[int, ...],
    ) -> tuple[EzLabel, ...]:
        mapped_labels = []
        for atom1, atom2, label in assignment.labels:
            mapped_atom1 = automorphism[atom1]
            mapped_atom2 = automorphism[atom2]
            mapped_labels.append(
                (
                    min(mapped_atom1, mapped_atom2),
                    max(mapped_atom1, mapped_atom2),
                    label,
                )
            )
        return tuple(sorted(mapped_labels))

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
            assignment_key(assignment, automorphism)
            for automorphism in automorphisms
        )
        if canonical_key in seen_keys:
            continue
        seen_keys.add(canonical_key)
        assignments.append(assignment)
    return tuple(assignments)


def analyze_ez(molecule_obj) -> EzAnalysis:
    """Return the full E/Z analysis for one molecule."""
    double_bonds = _find_ez_double_bonds_for_bonds(molecule_obj.bonds)
    return EzAnalysis(
        double_bonds=double_bonds,
        assignments=_enumerate_assignments(molecule_obj.bonds, double_bonds),
    )


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

    def compile_action(automorphism):
        action = []
        for kind, atom in elements:
            mapped_atom = automorphism[atom]
            parity = False
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
                    tuple(_map_ligand(v, automorphism) for v in source.ligands1),
                    tuple(_map_ligand(v, automorphism) for v in source.ligands2),
                )
                parity = _permutation_is_odd(mapped_pairs[0], target_pairs[0]) ^ (
                    _permutation_is_odd(mapped_pairs[1], target_pairs[1])
                )
            action.append((kind, mapped_atom, int(parity)))
        return tuple(action)

    actions = tuple(compile_action(automorphism) for automorphism in automorphisms)

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
    """Return tetrahedral and allene-like stereogenic elements and assignments."""
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
) -> frozenset[tuple[str, int]]:
    """Return centers active after all selected stereo decorations are applied."""
    tetra_by_atom = {center.atom: center for center in analysis.tetrahedral_centers}
    allene_by_center = {center.center_atom: center for center in analysis.allene_centers}
    raw_labels = tuple(sorted(assignment.labels))
    raw_ez = tuple(sorted(ez_assignment.labels)) if ez_assignment is not None else ()

    def map_labels(automorphism):
        mapped = []
        for kind, atom, bit in assignment.labels:
            mapped_atom = automorphism[atom]
            parity = False
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
                    tuple(_map_ligand(v, automorphism) for v in source.ligands1),
                    tuple(_map_ligand(v, automorphism) for v in source.ligands2),
                )
                parity = _permutation_is_odd(mapped_pairs[0], target_pairs[0]) ^ (
                    _permutation_is_odd(mapped_pairs[1], target_pairs[1])
                )
            mapped.append((kind, mapped_atom, bit ^ int(parity)))
        return tuple(sorted(mapped))

    def map_ez(automorphism):
        return tuple(
            sorted(
                (
                    min(automorphism[atom1], automorphism[atom2]),
                    max(automorphism[atom1], automorphism[atom2]),
                    label,
                )
                for atom1, atom2, label in raw_ez
            )
        )

    automorphisms = analysis.automorphisms or tuple(isomorphism.automorphisms(bonds))
    active = {("A", center.center_atom) for center in analysis.allene_centers}
    for center in analysis.tetrahedral_centers:
        other_labels = tuple(
            label for label in raw_labels if label[:2] != ("T", center.atom)
        )
        has_odd_stabilizer = False
        for automorphism in automorphisms:
            if automorphism[center.atom] != center.atom:
                continue
            target = tetra_by_atom[center.atom]
            mapped_ligands = tuple(
                _map_ligand(value, automorphism) for value in center.ligands
            )
            if not _permutation_is_odd(mapped_ligands, target.ligands):
                continue
            mapped_other = tuple(
                label
                for label in map_labels(automorphism)
                if label[:2] != ("T", center.atom)
            )
            if mapped_other == other_labels and map_ez(automorphism) == raw_ez:
                has_odd_stabilizer = True
                break
        if not has_odd_stabilizer:
            active.add(("T", center.atom))
    return frozenset(active)


def chiral_assignments_for_ez(
    bonds: np.ndarray,
    analysis: ChiralAnalysis,
    ez_assignment: EzAssignment,
) -> tuple[ChiralAssignment, ...]:
    """Enumerate chiral configurations under automorphisms preserving E/Z labels."""
    if not ez_assignment.labels:
        return analysis.assignments
    raw_ez = tuple(sorted(ez_assignment.labels))

    def mapped_ez(automorphism):
        return tuple(
            sorted(
                (
                    min(automorphism[atom1], automorphism[atom2]),
                    max(automorphism[atom1], automorphism[atom2]),
                    label,
                )
                for atom1, atom2, label in raw_ez
            )
        )

    preserving = tuple(
        automorphism
        for automorphism in analysis.automorphisms
        if mapped_ez(automorphism) == raw_ez
    )
    return _enumerate_chiral_assignments(
        bonds,
        analysis.tetrahedral_centers,
        analysis.allene_centers,
        preserving,
    )
