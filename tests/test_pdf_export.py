import importlib.util
import tempfile
import unittest
from pathlib import Path

import export_structures_pdf
import generation_pipeline


class PdfExportTests(unittest.TestCase):
    def test_pdf_cells_include_formula_headers(self):
        groups = [
            generation_pipeline.FormulaSmilesGroup("CH4", 1, 4, ["C"]),
            generation_pipeline.FormulaSmilesGroup("C2H0", 2, 0, []),
        ]

        cells = export_structures_pdf.build_pdf_cells(groups)

        self.assertEqual(
            cells,
            [
                export_structures_pdf.PdfCell("formula", "CH4"),
                export_structures_pdf.PdfCell("smiles", "C"),
                export_structures_pdf.PdfCell("formula", "C2H0"),
            ],
        )

    def test_grid_capacity_and_page_count(self):
        self.assertAlmostEqual(export_structures_pdf.PAGE_SIZE[0], 176 * 72 / 25.4)
        self.assertAlmostEqual(export_structures_pdf.PAGE_SIZE[1], 250 * 72 / 25.4)
        self.assertEqual(export_structures_pdf.GRID_COLUMNS, 10)
        self.assertEqual(export_structures_pdf.GRID_ROWS, 16)
        self.assertEqual(export_structures_pdf.GRID_CAPACITY, 160)
        self.assertEqual(export_structures_pdf.page_count_for_cell_count(0), 1)
        self.assertEqual(export_structures_pdf.page_count_for_cell_count(160), 1)
        self.assertEqual(export_structures_pdf.page_count_for_cell_count(161), 2)

    @unittest.skipUnless(
        importlib.util.find_spec("rdkit")
        and importlib.util.find_spec("reportlab")
        and importlib.util.find_spec("PIL"),
        "optional PDF dependencies are not installed",
    )
    def test_pdf_export_smoke(self):
        groups = [
            generation_pipeline.FormulaSmilesGroup("CH4", 1, 4, ["C"]),
            generation_pipeline.FormulaSmilesGroup("bad", 0, 0, ["not-a-smiles"]),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "structures.pdf"
            export_structures_pdf.export_formula_smiles_pdf(groups, output_path)

            self.assertTrue(output_path.exists())
            self.assertGreater(output_path.stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
