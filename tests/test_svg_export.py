import re
import tempfile
import unittest
from pathlib import Path

import export_structures_svg
import generation_pipeline


class SvgExportTests(unittest.TestCase):
    def test_output_path_for_page(self):
        self.assertEqual(
            export_structures_svg._output_path_for_page("structures_b5.svg", 2),
            Path("structures_b5_page_002.svg"),
        )

    def test_svg_export_smoke(self):
        groups = [
            generation_pipeline.FormulaSmilesGroup("CH4", 1, 4, ["C"]),
            generation_pipeline.FormulaSmilesGroup("bad", 0, 0, ["not-a-smiles"]),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            output_prefix = Path(tmpdir) / "structures"
            output_paths = export_structures_svg.export_formula_smiles_svg_pages(
                groups,
                output_prefix,
                report_warnings=False,
            )

            self.assertEqual(len(output_paths), 1)
            svg_text = output_paths[0].read_text(encoding="utf-8")
            self.assertIn("<svg", svg_text)
            self.assertIn("CH4", svg_text)
            self.assertIn("parse failed", svg_text)
            self.assertRegex(svg_text, re.compile(r"<path|<line|<text"))


if __name__ == "__main__":
    unittest.main()
