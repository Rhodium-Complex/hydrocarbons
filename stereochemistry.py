"""Late-stage E/Z stereochemistry helpers for hydrocarbon bond matrices."""
from __future__ import annotations

from dataclasses import dataclass
import itertools

import numpy as np

import isomorphism
import molecule


HYDROGEN_LIGAND = -1
EzLabel = tuple[int, int, str]


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


def _find_ez_double_bonds_for_bonds(bonds: np.ndarray) -> tuple[EzDoubleBond, ...]:
    hydrogens = molecule.implicit_hydrogens(bonds)
    double_bonds = []
    for atom1, bond_row in enumerate(bonds):
        for atom2, bond_order in enumerate(bond_row[atom1 + 1 :], start=atom1 + 1):
            if bond_order != 2:
                continue

            blocked = bonds.copy()
            blocked[atom1][atom2] = 0
            blocked[atom2][atom1] = 0
            colors = _refined_atom_colors(blocked)
            ligands1 = _ligands_for_atom(bonds, atom1, atom2, hydrogens)
            ligands2 = _ligands_for_atom(bonds, atom2, atom1, hydrogens)
            high1 = _high_priority_ligand(colors, ligands1)
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
