"""Benchmark the full generation pipeline by C/H step."""
import argparse
from pathlib import Path
import sys
import time
import tracemalloc

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import generation_pipeline
try:
    from benchmarks.process_memory import ProcessTreeMemory
except ModuleNotFoundError as exc:
    if exc.name != "psutil":
        raise
    raise SystemExit('Install benchmark dependencies: python -m pip install ".[benchmark]"') from exc


def benchmark_pipeline(
    min_carbon: int,
    max_carbon: int,
    workers: int | None,
    include_smiles: bool,
    trace_memory: bool,
) -> None:
    """Run the pipeline and print step timings in the standard log format."""

    def log_step(result: generation_pipeline.GenerationStepResult) -> None:
        memory_suffix = ""
        if trace_memory:
            _, peak = tracemalloc.get_traced_memory()
            memory_suffix = f" peak_bytes={peak}"
        print(generation_pipeline.format_step_result(result) + memory_suffix, flush=True)

    if trace_memory:
        tracemalloc.start()
    try:
        with ProcessTreeMemory() as memory:
            start = time.perf_counter()
            if include_smiles:
                outputs = generation_pipeline.run_generation_smiles_groups(
                    min_carbon=min_carbon,
                    max_carbon=max_carbon,
                    workers=workers,
                    log_step=log_step,
                )
            else:
                outputs = None
                generation_pipeline.run_generation(
                    min_carbon=min_carbon,
                    max_carbon=max_carbon,
                    workers=workers,
                    log_step=log_step,
                )
            wall_seconds = time.perf_counter() - start
        print(
            f"wall_seconds={wall_seconds:.4f} "
            f"peak_process_tree_rss_bytes={memory.peak_bytes}",
            flush=True,
        )
        del outputs
    finally:
        if trace_memory:
            tracemalloc.stop()


def main() -> None:
    """Parse arguments and run the pipeline benchmark."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-carbon", type=int, default=2)
    parser.add_argument("--max-carbon", type=int, default=8)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--include-smiles", action="store_true")
    parser.add_argument("--trace-memory", action="store_true")
    args = parser.parse_args()

    benchmark_pipeline(
        min_carbon=args.min_carbon,
        max_carbon=args.max_carbon,
        workers=args.workers,
        include_smiles=args.include_smiles,
        trace_memory=args.trace_memory,
    )


if __name__ == "__main__":
    main()
