"""Export generated hydrocarbon structures to editable B5 SVG pages."""
from __future__ import annotations

import argparse
import html
import importlib
import re
from pathlib import Path
import textwrap
from typing import Any

import export_structures_pdf
import generation_pipeline

PAGE_WIDTH_MM = 176
PAGE_HEIGHT_MM = 250
SVG_NAMESPACE = "http://www.w3.org/2000/svg"


def _load_svg_dependencies():
    try:
        chem = importlib.import_module("rdkit.Chem")
        draw2d = importlib.import_module("rdkit.Chem.Draw.rdMolDraw2D")
    except ImportError as exc:
        raise RuntimeError(
            "SVG export requires RDKit. Install it with: pip install rdkit"
        ) from exc
    return chem, draw2d


def _strip_svg_wrapper(svg_text: str) -> str:
    svg_text = re.sub(r"^\s*<\?xml[^>]*>\s*", "", svg_text)
    svg_text = re.sub(r"^\s*<svg[^>]*>", "", svg_text)
    svg_text = re.sub(r"</svg>\s*$", "", svg_text)
    return svg_text.strip()


def _smiles_to_svg_fragment(
    smiles: str,
    image_size: int,
    chem: Any,
    draw2d: Any,
) -> str | None:
    mol = chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    drawer = draw2d.MolDraw2DSVG(image_size, image_size)
    options = drawer.drawOptions()
    options.bondLineWidth = export_structures_pdf.STRUCTURE_BOND_LINE_WIDTH
    draw2d.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return _strip_svg_wrapper(drawer.GetDrawingText())


def _svg_text(text: str, x: float, y: float, size: float, weight: str = "normal") -> str:
    escaped = html.escape(text)
    return (
        f'<text x="{x:.3f}" y="{y:.3f}" text-anchor="middle" '
        f'dominant-baseline="middle" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size:.3f}" font-weight="{weight}">{escaped}</text>'
    )


def _fallback_text(smiles: str, x: float, y: float, width: float, height: float) -> str:
    lines = textwrap.wrap(smiles, width=18)[:4]
    lines.append("parse failed")
    line_height = 6.0
    start_y = y + height / 2 - (len(lines) - 1) * line_height / 2
    return "\n".join(
        _svg_text(line, x + width / 2, start_y + index * line_height, 5)
        for index, line in enumerate(lines)
    )


def _output_path_for_page(output_prefix: str | Path, page_number: int) -> Path:
    output_prefix = Path(output_prefix)
    stem = output_prefix.stem if output_prefix.suffix else output_prefix.name
    parent = output_prefix.parent
    return parent / f"{stem}_page_{page_number:03d}.svg"


def export_formula_smiles_svg_pages(
    groups: list[generation_pipeline.FormulaSmilesGroup],
    output_prefix: str | Path,
) -> list[Path]:
    """Write formula-grouped SMILES structures to editable B5 SVG pages."""
    chem, draw2d = _load_svg_dependencies()
    cells = export_structures_pdf.build_pdf_cells(groups)
    formula_lookup = export_structures_pdf.build_formula_lookup(cells)
    page_width, page_height = export_structures_pdf.PAGE_SIZE
    cell_width = (
        page_width - 2 * export_structures_pdf.PAGE_MARGIN
    ) / export_structures_pdf.GRID_COLUMNS
    cell_height = (
        page_height - 2 * export_structures_pdf.PAGE_MARGIN
    ) / export_structures_pdf.GRID_ROWS
    image_points = min(cell_width, cell_height) - 2 * export_structures_pdf.CELL_PADDING
    image_pixels = max(64, int(image_points * 2))
    scale = image_points / image_pixels
    page_count = export_structures_pdf.page_count_for_cell_count(len(cells))
    output_paths = []

    for page_index in range(page_count):
        page_cells = cells[
            page_index
            * export_structures_pdf.GRID_CAPACITY : (page_index + 1)
            * export_structures_pdf.GRID_CAPACITY
        ]
        elements = [
            (
                f'<svg xmlns="{SVG_NAMESPACE}" width="{PAGE_WIDTH_MM}mm" '
                f'height="{PAGE_HEIGHT_MM}mm" viewBox="0 0 {page_width:.3f} '
                f'{page_height:.3f}">'
            )
        ]

        for page_cell_index, cell in enumerate(page_cells):
            cell_index = page_index * export_structures_pdf.GRID_CAPACITY + page_cell_index
            position = export_structures_pdf.cell_position_for_index(
                cell_index,
                formula_lookup[cell_index],
            )
            row = page_cell_index // export_structures_pdf.GRID_COLUMNS
            column = page_cell_index % export_structures_pdf.GRID_COLUMNS
            x = export_structures_pdf.PAGE_MARGIN + column * cell_width
            y = export_structures_pdf.PAGE_MARGIN + row * cell_height

            if cell.kind == "formula":
                elements.append(
                    _svg_text(
                        cell.text,
                        x + cell_width / 2,
                        y + cell_height / 2,
                        10,
                        "bold",
                    )
                )
                continue

            with export_structures_pdf.capture_rdkit_warnings() as warnings:
                fragment = _smiles_to_svg_fragment(cell.text, image_pixels, chem, draw2d)
            export_structures_pdf.report_cell_warnings(warnings, position, cell.text)
            if fragment is None:
                elements.append(_fallback_text(cell.text, x, y, cell_width, cell_height))
                continue

            image_x = x + (cell_width - image_points) / 2
            image_y = y + (cell_height - image_points) / 2
            elements.append(
                f'<g transform="translate({image_x:.3f} {image_y:.3f}) '
                f'scale({scale:.6f})">{fragment}</g>'
            )

        elements.append("</svg>")
        output_path = _output_path_for_page(output_prefix, page_index + 1)
        output_path.write_text("\n".join(elements), encoding="utf-8")
        output_paths.append(output_path)

    return output_paths


def parse_args():
    """Parse command line arguments for B5 SVG export."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-carbon", type=int, default=2)
    parser.add_argument("--max-carbon", type=int, default=8)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--output-prefix", type=Path, default=Path("structures_b5"))
    return parser.parse_args()


def main() -> None:
    """Generate structures and export them to B5 SVG pages."""
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
        output_paths = export_formula_smiles_svg_pages(groups, args.output_prefix)
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from exc
    for output_path in output_paths:
        print(f"wrote {output_path}", flush=True)


if __name__ == "__main__":
    main()
