import importlib.util
import tempfile
import unittest
from pathlib import Path

import export_structures_pdf
import generation_output
import structure_export_layout as layout


class PdfExportTests(unittest.TestCase):
    def test_pdf_cells_include_formula_headers(self):
        groups = [
            generation_output.FormulaSmilesGroup("CH4", 1, 4, ["C"]),
            generation_output.FormulaSmilesGroup("C2H0", 2, 0, []),
        ]

        cells = layout.build_structure_cells(groups)

        self.assertEqual(
            cells,
            [
                layout.StructureCell("formula", "CH4"),
                layout.StructureCell("smiles", "C"),
                layout.StructureCell("formula", "C2H0"),
            ],
        )

    def test_grid_capacity_and_page_count(self):
        self.assertAlmostEqual(layout.PAGE_SIZE[0], 176 * 72 / 25.4)
        self.assertAlmostEqual(layout.PAGE_SIZE[1], 250 * 72 / 25.4)
        self.assertEqual(layout.GRID_COLUMNS, 10)
        self.assertEqual(layout.GRID_ROWS, 16)
        self.assertEqual(layout.GRID_CAPACITY, 160)
        self.assertEqual(layout.page_count_for_cell_count(0), 1)
        self.assertEqual(layout.page_count_for_cell_count(160), 1)
        self.assertEqual(layout.page_count_for_cell_count(161), 2)

    @unittest.skipUnless(
        importlib.util.find_spec("rdkit")
        and importlib.util.find_spec("reportlab")
        and importlib.util.find_spec("PIL"),
        "optional PDF dependencies are not installed",
    )
    def test_pdf_export_smoke(self):
        groups = [
            generation_output.FormulaSmilesGroup("CH4", 1, 4, ["C"]),
            generation_output.FormulaSmilesGroup("bad", 0, 0, ["not-a-smiles"]),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "structures.pdf"
            export_structures_pdf.export_formula_smiles_pdf(
                groups,
                output_path,
                report_warnings=False,
            )

            self.assertTrue(output_path.exists())
            self.assertGreater(output_path.stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
