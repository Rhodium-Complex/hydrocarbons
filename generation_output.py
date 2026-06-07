"""Output data structures and adapters for generated hydrocarbon SMILES."""
from __future__ import annotations

from dataclasses import dataclass

FORMULA_SEPARATOR = "N#N"


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


def flatten_legacy_smiles_groups(
    groups: list[FormulaSmilesGroup],
) -> list[str]:
    """Return the historical flat SMILES list with formula separators."""
    results = [FORMULA_SEPARATOR, FORMULA_SEPARATOR, FORMULA_SEPARATOR]
    for group in groups:
        if group.carbon_count == 1 and group.hydrogen_count == 4:
            results.extend(group.smiles)
            continue
        results.append(FORMULA_SEPARATOR)
        results.extend(group.smiles)
    return results
