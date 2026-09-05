"""Regression coverage for fused conversion and bounded process dispatch."""
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import hashlib
import json
import multiprocessing
import unittest
from unittest import mock

import numpy as np

import generation_pipeline as pipeline


def output_digest(groups):
    values = [
        (group.carbon_count, group.hydrogen_count,
         [variant.smiles for variant in group.variants])
        for group in groups
    ]
    return hashlib.sha256(json.dumps(values).encode()).hexdigest()


class ParallelPipelineTests(unittest.TestCase):
    # Captured from the original, unfused pipeline, C2..C6 plus methane.
    ORIGINAL_DIGEST = "37ae3f8f8c60f1249cf384934fc1d40429305af2c7c6eb2259819e0cc778b67b"

    def test_fused_output_matches_original_including_order(self):
        rows = []
        with mock.patch("generation_pipeline.ProcessPoolExecutor") as executor:
            groups = pipeline.run_generation_smiles_groups(
                2, 6, workers=1, log_step=rows.append,
            )
        executor.assert_not_called()
        self.assertEqual(output_digest(groups), self.ORIGINAL_DIGEST)
        for row in rows:
            self.assertTrue(row.smiles_fused)
            self.assertEqual(row.smiles_seconds, 0)
            self.assertIn("smiles=fused", pipeline.format_step_result(row))
            self.assertEqual(row.total_seconds, row.dehydro_seconds + row.build_seconds)

    def test_spawn_matches_original_and_uses_one_generation_pool(self):
        factory = partial(ProcessPoolExecutor, mp_context=multiprocessing.get_context("spawn"))
        with mock.patch("generation_pipeline.ProcessPoolExecutor", side_effect=factory) as pools:
            groups = pipeline.run_generation_smiles_groups(2, 6, workers=2)
        self.assertEqual(pools.call_count, 1)
        self.assertEqual(output_digest(groups), self.ORIGINAL_DIGEST)

    def test_fusion_preserves_groups_through_every_dehydrogenation(self):
        plain = list(pipeline._iter_formula_structure_steps(6, 6, workers=1))
        fused = list(pipeline._iter_formula_structure_steps(6, 6, workers=1, fuse_smiles=True))
        for expected, actual in zip(plain, fused, strict=True):
            self.assertEqual(expected.hydrogen_count, actual.hydrogen_count)
            self.assertEqual(len(expected.structures), len(actual.structures))
            for left, right in zip(expected.structures, actual.structures, strict=True):
                self.assertEqual(len(left), len(right))
                for a, b in zip(left, right, strict=True):
                    np.testing.assert_array_equal(a.bonds, b.bonds)
            self.assertEqual(
                actual.outputs,
                pipeline._structure_outputs(expected.structures, False),
            )

    def test_no_smiles_skips_conversion_and_empty_formula_stays_empty(self):
        rows = []
        with mock.patch("generation_pipeline.converter.mat2structure_variants") as convert:
            pipeline.run_generation(2, 2, workers=1, log_step=rows.append)
        convert.assert_not_called()
        self.assertEqual([row.count for row in rows], [1, 1, 1, 0])
        self.assertTrue(all(not row.smiles_fused for row in rows))
        groups = pipeline.run_generation_smiles_groups(2, 2, workers=1)
        self.assertEqual(groups[-1].variants, [])
        self.assertEqual(groups[1].variants[0].smiles, "CC")

    def test_stereo_keeps_separate_conversion_under_spawn(self):
        factory = partial(ProcessPoolExecutor, mp_context=multiprocessing.get_context("spawn"))
        for ez, tetra in ((True, False), (False, True), (True, True)):
            with self.subTest(ez=ez, tetra=tetra):
                options = dict(include_stereo=ez, include_tetrahedral_stereo=tetra)
                expected = pipeline.run_generation_smiles_groups(7, 7, workers=1, **options)
                rows = []
                with mock.patch("generation_pipeline.ProcessPoolExecutor", side_effect=factory):
                    actual = pipeline.run_generation_smiles_groups(
                        7, 7, workers=2, log_step=rows.append, **options,
                    )
                self.assertEqual(actual, expected)
                self.assertTrue(all(not row.smiles_fused for row in rows))
                self.assertGreater(sum(row.smiles_seconds for row in rows), 0)

    def test_dispatch_bounds_uncollected_batches_and_keeps_order(self):
        uncollected = 0
        peak = 0

        def submit(function, *args):
            nonlocal uncollected, peak
            uncollected += 1
            peak = max(peak, uncollected)
            future = mock.Mock()

            def result():
                nonlocal uncollected
                uncollected -= 1
                return function(*args)

            future.result.side_effect = result
            return future

        executor = mock.Mock()
        executor.submit.side_effect = submit
        result = list(pipeline._map_generation_task(
            executor, abs, range(-19, 0), chunksize=3, workers=2,
        ))
        self.assertEqual(result, list(range(19, 0, -1)))
        self.assertEqual(peak, 4)
        self.assertEqual(uncollected, 0)

    def test_dispatch_failure_cancels_pending_batches(self):
        futures = [mock.Mock() for _ in range(4)]
        futures[0].result.side_effect = ValueError("worker failed")
        executor = mock.Mock()
        executor.submit.side_effect = futures
        with self.assertRaisesRegex(ValueError, "worker failed"):
            list(pipeline._map_generation_task(executor, abs, range(10), workers=2))
        for future in futures[1:]:
            future.cancel.assert_called_once()


if __name__ == "__main__":
    unittest.main()
