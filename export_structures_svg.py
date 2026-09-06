"""Export generated hydrocarbon structures to editable B5 SVG pages."""
from __future__ import annotations

from collections.abc import Iterable
from collections import deque
from concurrent.futures import ProcessPoolExecutor
import html
import importlib
import multiprocessing
import re
from pathlib import Path
import textwrap

import generation_output
import structure_export_layout as layout

PAGE_WIDTH_MM = 142
PAGE_HEIGHT_MM = 219
SVG_NAMESPACE = "http://www.w3.org/2000/svg"
DEFAULT_OUTPUT_DIR = Path("outputs")
DEFAULT_SVG_OUTPUT_PREFIX = DEFAULT_OUTPUT_DIR / "structures_b5"
_worker_dependencies = None


def _initialize_svg_worker():
    global _worker_dependencies
    _worker_dependencies = _load_svg_dependencies()


def _write_svg_page_task(*args):
    return _write_svg_page(*args, dependencies=_worker_dependencies)


def _load_svg_dependencies():
    try:
        chem = importlib.import_module("rdkit.Chem")
        draw2d = importlib.import_module("rdkit.Chem.Draw.rdMolDraw2D")
        depictor = importlib.import_module("rdkit.Chem.rdDepictor")
    except ImportError as exc:
        raise RuntimeError(
            "SVG export requires RDKit. Install it with: pip install rdkit"
        ) from exc
    return chem, draw2d, depictor


def _strip_svg_wrapper(svg_text: str) -> str:
    svg_text = re.sub(r"^\s*<\?xml[^>]*>\s*", "", svg_text)
    svg_text = re.sub(r"^\s*<svg[^>]*>", "", svg_text)
    svg_text = re.sub(r"</svg>\s*$", "", svg_text)
    return svg_text.strip()


def _variant_to_svg_fragment(variant, image_size, chem, draw2d, depictor):
    mol = layout.prepare_variant_molecule(variant, chem, depictor)
    if mol is None:
        return None
    drawer = draw2d.MolDraw2DSVG(image_size, image_size)
    drawer.drawOptions().bondLineWidth = layout.STRUCTURE_BOND_LINE_WIDTH
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return _strip_svg_wrapper(drawer.GetDrawingText())


def _svg_text(text: str, x: float, y: float, size: float, weight: str = "normal") -> str:
    escaped = html.escape(text)
    return (
        f'<text x="{x:.3f}" y="{y:.3f}" text-anchor="middle" '
        f'dominant-baseline="middle" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size:.3f}" font-weight="{weight}">{escaped}</text>'
    )


def _fallback_text(smiles: str, x: float, y: float, width: float, height: float,
                   reason: str = "parse failed") -> str:
    lines = textwrap.wrap(smiles, width=18)[:4]
    lines.append("draw failed" if reason != "parse failed" else reason)
    line_height = 6.0
    start_y = y + height / 2 - (len(lines) - 1) * line_height / 2
    description = f"<desc>{html.escape(smiles + ': ' + reason)}</desc>\n"
    return description + "\n".join(
        _svg_text(line, x + width / 2, start_y + index * line_height, 5)
        for index, line in enumerate(lines)
    )


def _output_path_for_page(output_prefix: str | Path, page_number: int) -> Path:
    output_prefix = Path(output_prefix)
    stem = output_prefix.stem if output_prefix.suffix else output_prefix.name
    parent = output_prefix.parent
    parent.mkdir(parents=True, exist_ok=True)
    return parent / f"{stem}_page_{page_number:03d}.svg"


class SvgPageWriter:
    """Write full pages with a bounded queue; use as a context manager."""

    def __init__(self, output_prefix, report_warnings=True, on_page_written=None,
                 workers=1, on_progress=None):
        if workers < 1:
            raise ValueError("SVG workers must be at least 1")
        self._dependencies = _load_svg_dependencies()
        self._output_prefix = output_prefix
        self._report_warnings = report_warnings
        self._on_page_written = on_page_written
        self._on_progress = on_progress
        self._cells = []
        self._formula_labels = []
        self._current_formula = ""
        self._output_paths = []
        self._finished = False
        self._workers = workers
        self._executor = None
        self._pending = deque()
        self._submitted_pages = 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        if exc_type is None:
            self.finish()
        else:
            self.close()

    def close(self):
        """Release queued pages and stop workers, including on generation errors."""
        for future in self._pending:
            future.cancel()
        self._pending.clear()
        self._cells.clear()
        self._formula_labels.clear()
        self._finished = True
        if self._executor is not None:
            self._executor.shutdown(wait=True, cancel_futures=True)
            self._executor = None

    def add_group(self, group: generation_output.FormulaSmilesGroup) -> None:
        if self._finished:
            raise RuntimeError("Cannot add structures after SVG export has finished")
        self._current_formula = group.label
        self._add_cell(layout.StructureCell(kind="formula", text=group.label))
        for variant in group.variants:
            self._add_cell(layout.StructureCell(
                kind="smiles", text=variant.smiles, variant=variant,
            ))

    def _add_cell(self, cell):
        self._cells.append(cell)
        self._formula_labels.append(self._current_formula)
        if len(self._cells) == layout.GRID_CAPACITY:
            self._flush_page()

    def _flush_page(self):
        if self._workers > 1:
            if self._executor is None:
                # SVG submissions can happen while the generation pool is live.
                # Spawn avoids forking its threads and initialized RDKit state.
                self._executor = ProcessPoolExecutor(
                    max_workers=self._workers,
                    mp_context=multiprocessing.get_context("spawn"),
                    initializer=_initialize_svg_worker,
                )
            if len(self._pending) >= self._workers * 2:
                self._collect_page()
        args = (
            self._cells, self._formula_labels, self._submitted_pages + 1,
            self._output_prefix, self._report_warnings,
        )
        self._submitted_pages += 1
        self._report_progress()
        if self._executor is None:
            self._record_path(_write_svg_page(*args, dependencies=self._dependencies))
        else:
            self._pending.append(self._executor.submit(_write_svg_page_task, *args))
        # Submitted lists must remain intact until the executor serializes them.
        self._cells = []
        self._formula_labels = []
        while self._pending and self._pending[0].done():
            self._collect_page()

    def _collect_page(self):
        self._record_path(self._pending.popleft().result())

    def _record_path(self, path):
        self._output_paths.append(path)
        self._report_progress()
        if self._on_page_written is not None:
            self._on_page_written(path)

    def _report_progress(self):
        if self._on_progress is not None:
            self._on_progress(len(self._output_paths), self._submitted_pages)

    def finish(self) -> list[Path]:
        """Write the last partial page once, or one empty page for empty input."""
        if self._finished:
            return self._output_paths
        try:
            if self._cells or not self._submitted_pages:
                self._flush_page()
            while self._pending:
                self._collect_page()
            return self._output_paths
        finally:
            self.close()


def export_formula_smiles_svg_pages(
    groups: Iterable[generation_output.FormulaSmilesGroup],
    output_prefix: str | Path,
    report_warnings: bool = True,
    workers: int = 1,
) -> list[Path]:
    """Write iterable formula groups without flattening the entire export."""
    with SvgPageWriter(output_prefix, report_warnings=report_warnings, workers=workers) as writer:
        for group in groups:
            writer.add_group(group)
        return writer.finish()


def _write_svg_page(
    page_cells, formula_labels, page_number, output_prefix,
    report_warnings, dependencies,
) -> Path:
    chem, draw2d, depictor = dependencies
    page_width, page_height = layout.PAGE_SIZE
    cell_width = (
        page_width - 2 * layout.PAGE_MARGIN
    ) / layout.GRID_COLUMNS
    cell_height = (
        page_height - 2 * layout.PAGE_MARGIN
    ) / layout.GRID_ROWS
    image_points = min(cell_width, cell_height) - 2 * layout.CELL_PADDING
    image_pixels = max(64, int(image_points * 2))
    scale = image_points / image_pixels
    elements = [
        (
            f'<svg xmlns="{SVG_NAMESPACE}" width="{PAGE_WIDTH_MM}mm" '
            f'height="{PAGE_HEIGHT_MM}mm" viewBox="0 0 {page_width:.3f} '
            f'{page_height:.3f}">'
        )
    ]

    for page_cell_index, cell in enumerate(page_cells):
        cell_index = (page_number - 1) * layout.GRID_CAPACITY + page_cell_index
        position = layout.cell_position_for_index(
            cell_index,
            formula_labels[page_cell_index],
        )
        row = page_cell_index // layout.GRID_COLUMNS
        column = page_cell_index % layout.GRID_COLUMNS
        x = layout.PAGE_MARGIN + column * cell_width
        y = layout.PAGE_MARGIN + row * cell_height

        if cell.kind == "formula":
            elements.append(
                _svg_text(
                    cell.text,
                    x + cell_width / 2,
                    y + cell_height / 2,
                    8,
                    "bold",
                )
            )
            continue

        failure_reason = "parse failed"
        try:
            if report_warnings:
                with layout.capture_rdkit_warnings() as warnings:
                    fragment = _variant_to_svg_fragment(
                        cell.variant, image_pixels, chem, draw2d, depictor
                    )
                layout.report_cell_warnings(warnings, position, cell.text)
            else:
                rd_base = importlib.import_module("rdkit.rdBase")
                with rd_base.BlockLogs():
                    fragment = _variant_to_svg_fragment(
                        cell.variant, image_pixels, chem, draw2d, depictor
                    )
        except (RuntimeError, ValueError) as exc:
            # A bad depiction must not discard this page or stop enumeration.
            # File I/O and process failures remain fatal.
            fragment = None
            failure_reason = str(exc)
            if report_warnings:
                layout.report_cell_warnings(
                    [f"SVG draw failed: {exc}"], position, cell.text,
                )
        if fragment is None:
            elements.append(_fallback_text(
                cell.text, x, y, cell_width, cell_height, failure_reason,
            ))
            continue

        image_x = x + (cell_width - image_points) / 2
        image_y = y + (cell_height - image_points) / 2
        elements.append(
            f'<g transform="translate({image_x:.3f} {image_y:.3f}) '
            f'scale({scale:.6f})">{fragment}</g>'
        )

    elements.append("</svg>")
    output_path = _output_path_for_page(output_prefix, page_number)
    output_path.write_text("\n".join(elements), encoding="utf-8")
    return output_path
