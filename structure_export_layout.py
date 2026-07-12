"""Shared layout helpers for hydrocarbon structure exports."""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
import importlib
import logging
import math
import sys

import generation_output

MM_TO_POINTS = 72.0 / 25.4
PAGE_SIZE = (142 * MM_TO_POINTS, 219 * MM_TO_POINTS)
GRID_COLUMNS = 16
GRID_ROWS = 25
PAGE_MARGIN_MM = 0
CELL_PADDING_MM = 0
PAGE_MARGIN = PAGE_MARGIN_MM * MM_TO_POINTS
CELL_PADDING = CELL_PADDING_MM * MM_TO_POINTS
GRID_CAPACITY = GRID_COLUMNS * GRID_ROWS
STRUCTURE_BOND_LINE_WIDTH = 1.0


@dataclass(frozen=True)
class StructureCell:
    """One printable cell in the B5 structure grid."""

    kind: str
    text: str
    variant: generation_output.StructureVariant | None = None


@dataclass(frozen=True)
class CellPosition:
    """Printable grid position for diagnostics."""

    page: int
    row: int
    column: int
    formula: str


class _WarningCaptureHandler(logging.Handler):
    def __init__(self) -> None:
        super().__init__(level=logging.WARNING)
        self.messages: list[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.messages.append(record.getMessage())


def build_structure_cells(
    groups: list[generation_output.FormulaSmilesGroup],
) -> list[StructureCell]:
    """Flatten formula groups into header and structure cells."""
    cells = []
    for group in groups:
        cells.append(StructureCell(kind="formula", text=group.label))
        cells.extend(
            StructureCell(kind="smiles", text=variant.smiles, variant=variant)
            for variant in group.variants
        )
    return cells


def build_formula_lookup(cells: list[StructureCell]) -> list[str]:
    """Return the current formula label for each flattened cell."""
    labels = []
    current_label = ""
    for cell in cells:
        if cell.kind == "formula":
            current_label = cell.text
        labels.append(current_label)
    return labels


def cell_position_for_index(
    cell_index: int,
    formula_label: str,
) -> CellPosition:
    """Return one-based page, row, and column for a flattened cell index."""
    page_cell_index = cell_index % GRID_CAPACITY
    return CellPosition(
        page=cell_index // GRID_CAPACITY + 1,
        row=page_cell_index // GRID_COLUMNS + 1,
        column=page_cell_index % GRID_COLUMNS + 1,
        formula=formula_label,
    )


def page_count_for_cell_count(cell_count: int) -> int:
    """Return the number of B5 pages needed for printable cells."""
    return max(1, (cell_count + GRID_CAPACITY - 1) // GRID_CAPACITY)


@contextmanager
def capture_rdkit_warnings():
    """Capture RDKit warning records emitted while drawing one cell."""
    rd_base = importlib.import_module("rdkit.rdBase")
    rd_base.LogToPythonLogger()
    logger = logging.getLogger("rdkit")
    handler = _WarningCaptureHandler()
    logger.addHandler(handler)
    try:
        yield handler.messages
    finally:
        logger.removeHandler(handler)


def report_cell_warnings(
    warnings: list[str],
    position: CellPosition,
    smiles: str,
) -> None:
    """Print RDKit draw warnings with enough grid context to find the cell."""
    for warning in warnings:
        print(
            "RDKit warning at "
            f"page={position.page} row={position.row} column={position.column} "
            f"formula={position.formula} smiles={smiles}: {warning}",
            file=sys.stderr,
            flush=True,
        )


def prepare_variant_molecule(variant, chem, depictor):
    """Build a molecule and apply the requested odd-cumulene geometry."""
    if variant.bonds is None:
        return chem.MolFromSmiles(variant.smiles)
    bonds = variant.bonds
    editable = chem.RWMol()
    for _ in bonds:
        editable.AddAtom(chem.Atom(6))
    types = {1: chem.BondType.SINGLE, 2: chem.BondType.DOUBLE, 3: chem.BondType.TRIPLE}
    for left, row in enumerate(bonds):
        for right, order in enumerate(row[left + 1:], left + 1):
            if order:
                editable.AddBond(left, right, types[order])
    mol = editable.GetMol()
    chem.SanitizeMol(mol)
    depictor.Compute2DCoords(mol)
    conformer = mol.GetConformer()

    def explicit_ligand(terminal, neighbor, high):
        if high >= 0:
            return high, 1
        values = [i for i, order in enumerate(bonds[terminal]) if order and i != neighbor]
        return (values[0], -1) if values else (None, -1)

    def side(first, last, ligand):
        a, b, p = (conformer.GetAtomPosition(i) for i in (first, last, ligand))
        return (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x)

    def component(start, blocked):
        result, pending = set(), [start]
        while pending:
            atom = pending.pop()
            if atom in result or atom in blocked:
                continue
            result.add(atom)
            pending.extend(i for i, order in enumerate(bonds[atom]) if order)
        return result

    for stereo in variant.cumulenes:
        first, last = stereo.path[0], stereo.path[-1]
        ligand1, factor1 = explicit_ligand(first, stereo.path[1], stereo.high_ligand1)
        ligand2, factor2 = explicit_ligand(last, stereo.path[-2], stereo.high_ligand2)
        if ligand1 is None or ligand2 is None:
            continue
        same = side(first, last, ligand1) * factor1 * side(first, last, ligand2) * factor2 > 0
        if same == (stereo.configuration == "Z"):
            continue
        origin, target = conformer.GetAtomPosition(first), conformer.GetAtomPosition(last)
        dx, dy = target.x - origin.x, target.y - origin.y
        length = math.hypot(dx, dy)
        if not length:
            continue
        ux, uy = dx / length, dy / length
        for atom in component(ligand2, set(stereo.path)):
            point = conformer.GetAtomPosition(atom)
            rx, ry = point.x - origin.x, point.y - origin.y
            parallel, perpendicular = rx * ux + ry * uy, rx * -uy + ry * ux
            point.x = origin.x + parallel * ux + perpendicular * uy
            point.y = origin.y + parallel * uy - perpendicular * ux
            conformer.SetAtomPosition(atom, point)
    return mol
