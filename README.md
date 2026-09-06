# Molecular Generator Program
This Python program is designed to enumerate molecules composed of carbon and hydrogen. It allows users to input a specific molecular formula and checks whether the given formula is composed of carbon and hydrogen.

## Novelty
This implementation aims to fill the gap in the comprehensive enumeration of hydrocarbons. Despite the significant interest in hydrocarbon enumeration, there's a lack of comprehensive databases or programs addressing this, making this one of the initial implementations in this domain.

## Usage
Run the generator from the project root:

```powershell
.\.venv\Scripts\python.exe -u main.py --max-carbon 8 --workers 4
```

For faster enumeration-only runs, skip SMILES conversion:

```powershell
.\.venv\Scripts\python.exe -u main.py --max-carbon 8 --workers 4 --no-smiles
```

Useful options:

- `--min-carbon`: first carbon count to generate.
- `--max-carbon`: last carbon count to generate.
- `--workers`: number of worker processes in the shared generation pool;
  `1` runs sequentially. Ordinary SMILES conversion runs in the same tasks.
- `--svg-workers`: SVG page-rendering processes, separate from generation
  workers; defaults to `1` (sequential rendering).
- `--no-smiles`: skip SMILES conversion and report structure counts only.
- `--include-stereo`: expand formal E/Z stereoisomers during final SMILES output.
- `--include-tetrahedral-stereo`: expand tetrahedral-carbon and allene-like
  configurations using `@`/`@@` and `@AL1`/`@AL2` SMILES notation.

For E/Z-aware SMILES output, enable stereo expansion:

```powershell
.\.venv\Scripts\python.exe -u main.py --max-carbon 8 --workers 4 --include-stereo
```

Atom-centered stereochemistry can be enabled independently or together with E/Z:

```powershell
.\.venv\Scripts\python.exe -u main.py --max-carbon 8 --workers 4 --include-tetrahedral-stereo
.\.venv\Scripts\python.exe -u main.py --max-carbon 8 --workers 4 --include-stereo --include-tetrahedral-stereo
```

The binary configurations are enumeration labels, not guaranteed CIP R/S or
Ra/Sa names. Tetrahedral carbon and even cumulene (allene-like) stereochemistry
are supported; helicity and conformational stereochemistry are not.

Tetrahedral enumeration includes configuration-dependent centers: a CH or
quaternary carbon may become stereogenic only after neighboring configurations
are assigned. Symmetry is evaluated on the fully decorated molecular graph, so
inactive equal-ligand candidates do not receive an `@` marker. For acyclic
saturated hydrocarbons, counts are regression-tested against OEIS A000628.

Export generated structures to a B5 PDF grid:

```powershell
.\.venv\Scripts\python.exe -u main.py --max-carbon 8 --workers 4 --export pdf
```

Export editable B5 SVG pages:

```powershell
.\.venv\Scripts\python.exe -u main.py --max-carbon 8 --workers 4 --export svg
```

SVG export writes full pages as each molecular formula (carbon/hydrogen count)
is generated, instead of collecting all formulas first. A page contains 400
cells, including formula headings. An unfinished page carries over to the next
formula; the last partial page is written when generation finishes. The terminal
shows a single SVG progress bar (completed / submitted pages). It is cleared
when the submitted pages finish, before each formula timing line, and when the
export exits. The denominator is pages submitted so far, not a forecast of the
entire enumeration. Redirected output contains no progress bar or per-file logs.

If RDKit raises a drawing/coordinate error such as
`Cannot normalize a zero length vector`, export continues with the next cell.
The failed cell shows its SMILES and `draw failed`; its SVG `<desc>` retains the
full SMILES and error message. Drawing warnings still report the page, row,
column and formula. Disk-write errors and worker-process failures still stop
the export.

To render multiple SVG pages concurrently:

```powershell
.\.venv\Scripts\python.exe -u main.py --max-carbon 10 --workers 12 --export svg --svg-workers 4
```

SVG workers render and save separate pages while generation can continue.
At most `2 * svg_workers` pages are submitted but not yet collected; generation
waits when that limit is reached. Page numbers and returned file order remain
stable even if files finish in a different order. The SVG pool adds processes
and memory usage alongside the generation pool, so its size is controlled
independently. Use `--svg-workers 1` for the previous synchronous behavior.

Exported SMILES lists are released after each formula. The current formula's
output and structures needed for the next dehydrogenation still occupy memory;
this does not put a fixed memory limit on structure generation. PDF export
continues to collect its input before writing.

By default, PDF and SVG files are written under `outputs/`, which is ignored by Git.
Export also accepts `--include-stereo`, `--output`, and `--output-prefix`.

## 構造マップ
`main.main()` から生成、SMILES 変換、PDF/SVG 出力までの大まかな呼び出し関係は
[ARCHITECTURE.md](ARCHITECTURE.md) にまとめています。

The PDF export requires the optional `pdf` dependencies.

```powershell
.\.venv\Scripts\python.exe -m pip install ".[pdf]"
```

The command prints one timing line per carbon/hydrogen step. Ordinary SMILES
conversion runs after deduplication inside each dehydrogenation/build task:

```text
C= 8 H= 6 dehydro=   0.4120s build=   0.2810s smiles=fused count=   7982 total=   0.6930s
```

With `smiles=fused`, `dehydro` and `build` include SMILES conversion. Stereo
conversion remains a separate stage with its own timing. With `--no-smiles`,
conversion is skipped:

```text
C= 8 H= 6 dehydro=   1.2144s build=   0.4415s smiles=   0.0000s count=   7982 total=   1.6560s
```

## Benchmarks
Install the optional benchmark dependency first:

```powershell
.\.venv\Scripts\python.exe -m pip install ".[benchmark]"
```

Benchmark the full generation pipeline:

```powershell
.\.venv\Scripts\python.exe -u benchmarks\benchmark_pipeline.py --max-carbon 8 --workers 4
```

Use `--include-smiles` to benchmark the fused ordinary SMILES path. For CPU
scaling comparisons, run C8 and C9 separately with `--min-carbon` equal to
`--max-carbon`, compare `--workers 1`, `4`, `12`, and `24`, and take the median
of three runs per setting.

Each run reports `wall_seconds`, including process pool startup/shutdown, and
`peak_process_tree_rss_bytes`, the peak sampled sum of parent/child RSS at
20 ms intervals. RSS counts shared pages in each process and can miss shorter
peaks; it is not unique physical memory usage.

Also include parent-process Python allocation tracking (adds profiling overhead):

```powershell
.\.venv\Scripts\python.exe -u benchmarks\benchmark_pipeline.py --max-carbon 8 --workers 4 --trace-memory
```

Benchmark lower-level generation primitives:

```powershell
.\.venv\Scripts\python.exe benchmarks\benchmark_generation.py --max-carbon 5
```

## Tests
Run the correctness tests from the project root:

```powershell
.\.venv\Scripts\python.exe -m unittest -v
```

## Language and Libraries Used
Python 3.11+ with NumPy.
The optional PDF export uses RDKit and ReportLab. The SVG export uses RDKit.
## Important Note
This program conducts a basic verification of the input molecule. It is not suitable for advanced chemical evaluations or in-depth structural analysis.

**Disclaimer**: The use of this program is for reference purposes only. It does not guarantee chemical accuracy.
