"""Export generated hydrocarbon structures to a B5 PDF grid."""
from __future__ import annotations

from io import BytesIO
import importlib
from pathlib import Path
import textwrap
from typing import Any

import generation_output
import structure_export_layout as layout

DEFAULT_OUTPUT_DIR = Path("outputs")
DEFAULT_PDF_OUTPUT = DEFAULT_OUTPUT_DIR / "structures_b5.pdf"


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
            "Install them with: pip install rdkit reportlab"
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
    options.bondLineWidth = layout.STRUCTURE_BOND_LINE_WIDTH
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
    groups: list[generation_output.FormulaSmilesGroup],
    output_path: str | Path,
    report_warnings: bool = True,
) -> None:
    """Write formula-grouped SMILES structures to a B5 portrait PDF."""
    chem, draw2d, colors, image_reader, canvas = _load_pdf_dependencies()
    cells = layout.build_structure_cells(groups)
    formula_lookup = layout.build_formula_lookup(cells)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    page_width, page_height = layout.PAGE_SIZE
    cell_width = (page_width - 2 * layout.PAGE_MARGIN) / layout.GRID_COLUMNS
    cell_height = (page_height - 2 * layout.PAGE_MARGIN) / layout.GRID_ROWS
    image_points = min(cell_width, cell_height) - 2 * layout.CELL_PADDING
    image_pixels = max(64, int(image_points * 2))

    pdf = canvas.Canvas(str(output_path), pagesize=layout.PAGE_SIZE)
    for cell_index, cell in enumerate(cells):
        if cell_index and cell_index % layout.GRID_CAPACITY == 0:
            pdf.showPage()

        position = layout.cell_position_for_index(cell_index, formula_lookup[cell_index])
        row = position.row - 1
        column = position.column - 1
        x = layout.PAGE_MARGIN + column * cell_width
        y = page_height - layout.PAGE_MARGIN - (row + 1) * cell_height

        if cell.kind == "formula":
            pdf.setFont("Helvetica-Bold", 10)
            pdf.setFillColor(colors.black)
            _draw_centered_text(pdf, cell.text, x, y, cell_width, cell_height)
            continue

        if report_warnings:
            with layout.capture_rdkit_warnings() as warnings:
                png_bytes = _smiles_to_png_bytes(cell.text, image_pixels, chem, draw2d)
            layout.report_cell_warnings(warnings, position, cell.text)
        else:
            rd_base = importlib.import_module("rdkit.rdBase")
            with rd_base.BlockLogs():
                png_bytes = _smiles_to_png_bytes(cell.text, image_pixels, chem, draw2d)
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
