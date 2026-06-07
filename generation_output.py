"""Output data structures for generated hydrocarbon SMILES."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FormulaSmilesGroup:
    """SMILES strings grouped by molecular formula."""

    label: str
    carbon_count: int
    hydrogen_count: int
    smiles: list[str]


def format_formula_label(carbon_count: int, hydrogen_count: int) -> str:
    """Return a compact hydrocarbon formula label."""
    if carbon_count == 1:
        return f"CH{hydrogen_count}"
    return f"C{carbon_count}H{hydrogen_count}"
