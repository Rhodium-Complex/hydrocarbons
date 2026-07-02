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
- `--workers`: number of worker processes for generation.
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

By default, PDF and SVG files are written under `outputs/`, which is ignored by Git.
Export also accepts `--include-stereo`, `--output`, and `--output-prefix`.

## 構造マップ
`main.main()` から生成、SMILES 変換、PDF/SVG 出力までの大まかな呼び出し関係は
[ARCHITECTURE.md](ARCHITECTURE.md) にまとめています。

The PDF export requires the optional `pdf` dependencies.

```powershell
.\.venv\Scripts\python.exe -m pip install ".[pdf]"
```

The command prints one timing line per carbon/hydrogen step:

```text
C= 8 H= 6 dehydro=   1.2144s build=   0.4415s smiles=   0.0000s count=   7982 total=   1.6560s
```

## Benchmarks
Benchmark the full generation pipeline:

```powershell
.\.venv\Scripts\python.exe -u benchmarks\benchmark_pipeline.py --max-carbon 8 --workers 4
```

Include memory tracking:

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
