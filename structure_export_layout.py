"""Shared layout helpers for hydrocarbon structure exports."""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
import importlib
import logging
import sys

import generation_pipeline

MM_TO_POINTS = 72.0 / 25.4
PAGE_SIZE = (176 * MM_TO_POINTS, 250 * MM_TO_POINTS)
GRID_COLUMNS = 10
GRID_ROWS = 16
PAGE_MARGIN_MM = 6
CELL_PADDING_MM = 1
PAGE_MARGIN = PAGE_MARGIN_MM * MM_TO_POINTS
CELL_PADDING = CELL_PADDING_MM * MM_TO_POINTS
GRID_CAPACITY = GRID_COLUMNS * GRID_ROWS
STRUCTURE_BOND_LINE_WIDTH = 1.0


@dataclass(frozen=True)
class StructureCell:
    """One printable cell in the B5 structure grid."""

    kind: str
    text: str


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
    groups: list[generation_pipeline.FormulaSmilesGroup],
) -> list[StructureCell]:
    """Flatten formula groups into header and structure cells."""
    cells = []
    for group in groups:
        cells.append(StructureCell(kind="formula", text=group.label))
        cells.extend(StructureCell(kind="smiles", text=smiles) for smiles in group.smiles)
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
