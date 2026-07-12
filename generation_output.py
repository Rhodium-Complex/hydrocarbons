"""Output data structures for generated hydrocarbon SMILES."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CumuleneStereo:
    """Drawable relative configuration for one odd cumulene chain."""

    path: tuple[int, ...]
    high_ligand1: int
    high_ligand2: int
    configuration: str


@dataclass(frozen=True)
class StructureVariant:
    """One stereochemical structure and the data required to draw it."""

    smiles: str
    bonds: tuple[tuple[int, ...], ...] | None = None
    cumulenes: tuple[CumuleneStereo, ...] = ()


@dataclass(frozen=True)
class FormulaSmilesGroup:
    """SMILES strings grouped by molecular formula."""

    label: str
    carbon_count: int
    hydrogen_count: int
    variants: list[StructureVariant]


def format_formula_label(carbon_count: int, hydrogen_count: int) -> str:
    """Return a compact hydrocarbon formula label."""
    if carbon_count == 1:
        return f"CH{hydrogen_count}"
    return f"C{carbon_count}H{hydrogen_count}"
