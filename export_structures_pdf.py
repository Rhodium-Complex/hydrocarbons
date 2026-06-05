"""Export generated hydrocarbon structures to a B5 PDF grid."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import dataclass
from io import BytesIO
import importlib
import logging
from pathlib import Path
import sys
import textwrap
from typing import Any

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
class PdfCell:
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


def build_pdf_cells(
    groups: list[generation_pipeline.FormulaSmilesGroup],
) -> list[PdfCell]:
    """Flatten formula groups into header and structure cells."""
    cells = []
    for group in groups:
        cells.append(PdfCell(kind="formula", text=group.label))
        cells.extend(PdfCell(kind="smiles", text=smiles) for smiles in group.smiles)
    return cells


def build_formula_lookup(cells: list[PdfCell]) -> list[str]:
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


def _load_pdf_dependencies():
    try:
        chem = importlib.import_module("rdkit.Chem")
        draw2d = importlib.import_module("rdkit.Chem.Draw.rdMolDraw2D")
        colors = importlib.import_module("reportlab.lib.colors")
        image_reader = importlib.import_module("reportlab.lib.utils").ImageReader
        canvas = importlib.import_module("reportlab.pdfgen.canvas")
    except ImportError as exc:
        raise RuntimeError(
            "PDF export requires optional dependencies. "
            "Install them with: pip install rdkit reportlab pillow"
        ) from exc
    return chem, draw2d, colors, image_reader, canvas


def _smiles_to_png_bytes(
    smiles: str,
    image_size: int,
    chem: Any,
    draw2d: Any,
) -> bytes | None:
    mol = chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    drawer = draw2d.MolDraw2DCairo(image_size, image_size)
    options = drawer.drawOptions()
    options.bondLineWidth = STRUCTURE_BOND_LINE_WIDTH
    draw2d.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _draw_centered_text(pdf, text: str, x: float, y: float, width: float, height: float):
    pdf.drawCentredString(x + width / 2, y + height / 2 - 4, text)


def _draw_fallback_text(pdf, smiles: str, x: float, y: float, width: float, height: float):
    pdf.setFont("Helvetica", 5)
    lines = textwrap.wrap(smiles, width=18)[:4]
    lines.append("parse failed")
    line_height = 6
    total_height = len(lines) * line_height
    start_y = y + (height + total_height) / 2 - line_height
    for index, line in enumerate(lines):
        pdf.drawCentredString(x + width / 2, start_y - index * line_height, line)


def export_formula_smiles_pdf(
    groups: list[generation_pipeline.FormulaSmilesGroup],
    output_path: str | Path,
) -> None:
    """Write formula-grouped SMILES structures to a B5 portrait PDF."""
    chem, draw2d, colors, image_reader, canvas = _load_pdf_dependencies()
    cells = build_pdf_cells(groups)
    formula_lookup = build_formula_lookup(cells)
    output_path = Path(output_path)

    page_width, page_height = PAGE_SIZE
    cell_width = (page_width - 2 * PAGE_MARGIN) / GRID_COLUMNS
    cell_height = (page_height - 2 * PAGE_MARGIN) / GRID_ROWS
    image_points = min(cell_width, cell_height) - 2 * CELL_PADDING
    image_pixels = max(64, int(image_points * 2))

    pdf = canvas.Canvas(str(output_path), pagesize=PAGE_SIZE)
    for cell_index, cell in enumerate(cells):
        if cell_index and cell_index % GRID_CAPACITY == 0:
            pdf.showPage()

        position = cell_position_for_index(cell_index, formula_lookup[cell_index])
        row = position.row - 1
        column = position.column - 1
        x = PAGE_MARGIN + column * cell_width
        y = page_height - PAGE_MARGIN - (row + 1) * cell_height

        if cell.kind == "formula":
            pdf.setFont("Helvetica-Bold", 10)
            pdf.setFillColor(colors.black)
            _draw_centered_text(pdf, cell.text, x, y, cell_width, cell_height)
            continue

        with capture_rdkit_warnings() as warnings:
            png_bytes = _smiles_to_png_bytes(cell.text, image_pixels, chem, draw2d)
        report_cell_warnings(warnings, position, cell.text)
        if png_bytes is None:
            pdf.setFillColor(colors.black)
            _draw_fallback_text(pdf, cell.text, x, y, cell_width, cell_height)
            continue

        image = image_reader(BytesIO(png_bytes))
        image_x = x + (cell_width - image_points) / 2
        image_y = y + (cell_height - image_points) / 2
        pdf.drawImage(
            image,
            image_x,
            image_y,
            width=image_points,
            height=image_points,
            preserveAspectRatio=True,
            anchor="c",
            mask="auto",
        )

    pdf.save()


def parse_args():
    """Parse command line arguments for B5 PDF export."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-carbon", type=int, default=2)
    parser.add_argument("--max-carbon", type=int, default=8)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--output", type=Path, default=Path("structures_b5.pdf"))
    return parser.parse_args()


def main() -> None:
    """Generate structures and export them to a B5 PDF."""
    args = parse_args()
    groups = generation_pipeline.run_generation_smiles_groups(
        min_carbon=args.min_carbon,
        max_carbon=args.max_carbon,
        workers=args.workers,
        include_methane=True,
        log_step=lambda result: print(
            generation_pipeline.format_step_result(result),
            flush=True,
        ),
    )
    try:
        export_formula_smiles_pdf(groups, args.output)
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from exc
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
