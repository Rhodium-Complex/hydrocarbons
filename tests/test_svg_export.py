import re
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import export_structures_svg
import generation_output
import generation_pipeline
import main


class SvgExportTests(unittest.TestCase):
    def test_zero_length_vector_error_keeps_cell_and_continues_page(self):
        groups = [generation_output.FormulaSmilesGroup("C2H4", 2, 4, [
            generation_output.StructureVariant("C=C"),
            generation_output.StructureVariant("CC"),
        ])]
        for error_type in (RuntimeError, ValueError):
            with self.subTest(error_type=error_type), tempfile.TemporaryDirectory() as tmpdir:
                with mock.patch(
                    "export_structures_svg._variant_to_svg_fragment",
                    side_effect=[error_type("Cannot normalize a zero length vector"),
                                 '<path id="next-molecule" />'],
                ):
                    paths = export_structures_svg.export_formula_smiles_svg_pages(
                        groups, Path(tmpdir) / "structures", report_warnings=False,
                    )
                text = paths[0].read_text()
                self.assertIn("draw failed", text)
                self.assertIn("C=C: Cannot normalize a zero length vector", text)
                self.assertIn('id="next-molecule"', text)

    def test_progress_reports_submitted_and_completed_pages(self):
        with tempfile.TemporaryDirectory() as tmpdir, mock.patch(
            "structure_export_layout.GRID_CAPACITY", 1,
        ):
            progress = mock.Mock()
            with export_structures_svg.SvgPageWriter(
                Path(tmpdir) / "structures", on_progress=progress,
            ) as writer:
                writer.add_group(generation_output.FormulaSmilesGroup("empty", 2, 0, []))
            self.assertEqual(progress.call_args_list, [mock.call(0, 1), mock.call(1, 1)])

    def test_parallel_spawn_pages_match_serial_pages(self):
        variants = [generation_output.StructureVariant("C=C") for _ in range(401)]
        variants[-1] = generation_output.StructureVariant("not-a-smiles")
        groups = [generation_output.FormulaSmilesGroup("C2H4", 2, 4, variants)]
        with tempfile.TemporaryDirectory() as tmpdir:
            before = export_structures_svg.export_formula_smiles_svg_pages(
                groups, Path(tmpdir) / "serial", report_warnings=False,
            )
            after = export_structures_svg.export_formula_smiles_svg_pages(
                iter(groups), Path(tmpdir) / "parallel", report_warnings=False, workers=2,
            )
            self.assertEqual(len(after), 2)
            self.assertEqual(
                [path.read_text() for path in before],
                [path.read_text() for path in after],
            )

    def test_parallel_queue_is_bounded_and_submitted_cells_are_not_cleared(self):
        outstanding = peak = 0

        def submit(function, cells, labels, number, prefix, warnings):
            nonlocal outstanding, peak
            outstanding += 1
            peak = max(peak, outstanding)
            future = mock.Mock()
            future.done.return_value = False

            def result():
                nonlocal outstanding
                self.assertEqual(len(cells), 1)
                self.assertEqual(labels, [str(number)])
                outstanding -= 1
                return Path(f"page_{number}.svg")

            future.result.side_effect = result
            return future

        with mock.patch("export_structures_svg.ProcessPoolExecutor") as pool, mock.patch(
            "structure_export_layout.GRID_CAPACITY", 1,
        ):
            pool.return_value.submit.side_effect = submit
            written = []
            with export_structures_svg.SvgPageWriter(
                "unused", workers=2, on_page_written=written.append,
            ) as writer:
                for number in range(1, 10):
                    writer.add_group(generation_output.FormulaSmilesGroup(str(number), 2, 0, []))
                paths = writer.finish()
            self.assertEqual(peak, 4)
            self.assertEqual(outstanding, 0)
            self.assertEqual(paths, [Path(f"page_{i}.svg") for i in range(1, 10)])
            self.assertEqual(written, paths)
            pool.return_value.shutdown.assert_called_once_with(wait=True, cancel_futures=True)

    def test_parallel_write_failure_propagates_and_stops_workers(self):
        with mock.patch("export_structures_svg.ProcessPoolExecutor") as pool:
            future = pool.return_value.submit.return_value
            future.done.return_value = False
            future.result.side_effect = OSError("disk full")
            with self.assertRaisesRegex(OSError, "disk full"):
                with export_structures_svg.SvgPageWriter("unused", workers=2) as writer:
                    writer.add_group(generation_output.FormulaSmilesGroup("empty", 2, 0, []))
            pool.return_value.shutdown.assert_called_once_with(wait=True, cancel_futures=True)

    def test_full_pages_are_written_and_partial_page_carries_between_formulas(self):
        variant = generation_output.StructureVariant
        group = generation_output.FormulaSmilesGroup
        with tempfile.TemporaryDirectory() as tmpdir, mock.patch(
            "structure_export_layout.GRID_CAPACITY", 4,
        ), mock.patch(
            "export_structures_svg._variant_to_svg_fragment", return_value="<path />",
        ), mock.patch("structure_export_layout.report_cell_warnings") as warnings:
            written = []
            writer = export_structures_svg.SvgPageWriter(
                Path(tmpdir) / "structures", on_page_written=written.append,
            )
            writer.add_group(group("A", 2, 6, [variant("CC"), variant("C=C")]))
            self.assertEqual(written, [])
            writer.add_group(group("B", 2, 2, [variant("C#C")]))
            self.assertEqual(len(written), 1)
            self.assertTrue(written[0].exists())
            writer.add_group(group("empty", 2, 0, []))
            paths = writer.finish()
            self.assertEqual(len(paths), 2)
            self.assertIn(">A</text>", paths[0].read_text())
            self.assertIn(">B</text>", paths[0].read_text())
            self.assertNotIn(">B</text>", paths[1].read_text())
            self.assertIn(">empty</text>", paths[1].read_text())
            # The molecule continuing onto page two still has formula B.
            self.assertEqual(warnings.call_args_list[-1].args[1].formula, "B")
            self.assertEqual(warnings.call_args_list[-1].args[1].page, 2)
            self.assertEqual(writer.finish(), paths)
            self.assertEqual(len(written), 2)

    def test_exact_page_and_empty_input_have_no_extra_page(self):
        with tempfile.TemporaryDirectory() as tmpdir, mock.patch(
            "structure_export_layout.GRID_CAPACITY", 1,
        ):
            writer = export_structures_svg.SvgPageWriter(Path(tmpdir) / "full")
            writer.add_group(generation_output.FormulaSmilesGroup("empty", 2, 0, []))
            self.assertEqual(len(writer.finish()), 1)
            paths = export_structures_svg.export_formula_smiles_svg_pages(
                iter(()), Path(tmpdir) / "empty",
            )
            self.assertEqual(len(paths), 1)
            self.assertIn("</svg>", paths[0].read_text())

    def test_streamed_groups_match_collected_output_and_release_lists(self):
        for stereo in (False, True):
            with self.subTest(stereo=stereo):
                expected = generation_pipeline.run_generation_smiles_groups(
                    2, 4, workers=1,
                    include_stereo=stereo, include_tetrahedral_stereo=stereo,
                )
                observed, borrowed = [], []

                def consume(group):
                    self.assertTrue(all(not items for items in borrowed))
                    observed.append(generation_output.FormulaSmilesGroup(
                        group.label, group.carbon_count, group.hydrogen_count,
                        list(group.variants),
                    ))
                    borrowed.append(group.variants)

                with mock.patch("generation_pipeline.print_step_result") as log:
                    generation_pipeline.stream_export_smiles_groups(
                        2, 4, consume, workers=1,
                        include_stereo=stereo, include_tetrahedral_stereo=stereo,
                    )
                self.assertEqual(observed, expected)
                self.assertTrue(all(not items for items in borrowed))
                self.assertEqual(
                    sum(call.args[0].count for call in log.call_args_list),
                    sum(len(group.variants) for group in expected[1:]),
                )

    def test_svg_cli_uses_streaming_generation(self):
        with tempfile.TemporaryDirectory() as tmpdir, mock.patch(
            "generation_pipeline.run_export_smiles_groups",
        ) as collect, mock.patch("structure_export_layout.GRID_CAPACITY", 4), mock.patch(
            "builtins.print",
        ):
            paths = main.main(
                min_carbon=2, max_carbon=2, workers=1,
                export_format="svg", output_prefix=Path(tmpdir) / "structures",
            )
            collect.assert_not_called()
            self.assertEqual(len(paths), 3)
            self.assertTrue(all(path.exists() for path in paths))

    def test_output_path_for_page(self):
        self.assertEqual(
            export_structures_svg._output_path_for_page("structures_b5.svg", 2),
            Path("structures_b5_page_002.svg"),
        )

    def test_svg_export_smoke(self):
        groups = [
            generation_output.FormulaSmilesGroup("CH4", 1, 4, [generation_output.StructureVariant("C")]),
            generation_output.FormulaSmilesGroup("bad", 0, 0, [generation_output.StructureVariant("not-a-smiles")]),
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
