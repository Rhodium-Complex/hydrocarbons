"""Benchmark delayed E/Z stereochemistry expansion."""
import argparse
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import generation_pipeline


def benchmark_stereo(
    min_carbon: int,
    max_carbon: int,
    workers: int | None,
) -> None:
    """Print timing and output counts for stereo generation."""
    start = time.perf_counter()
    groups = generation_pipeline.run_generation_smiles_groups(
        min_carbon=min_carbon,
        max_carbon=max_carbon,
        workers=workers,
        include_methane=False,
        include_stereo=True,
    )
    elapsed = time.perf_counter() - start
    output_count = sum(len(group.variants) for group in groups)
    print(
        f"min_carbon={min_carbon} max_carbon={max_carbon} "
        f"workers={workers} stereo_outputs={output_count} "
        f"seconds={elapsed:.4f}",
        flush=True,
    )
    for group in groups:
        print(f"{group.label} {len(group.variants)}", flush=True)


def parse_args():
    """Parse benchmark arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-carbon", type=int, default=6)
    parser.add_argument("--max-carbon", type=int, default=6)
    parser.add_argument("--workers", type=int, default=1)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    benchmark_stereo(
        min_carbon=args.min_carbon,
        max_carbon=args.max_carbon,
        workers=args.workers,
    )
