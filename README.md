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

For E/Z-aware SMILES output, enable stereo expansion:

```powershell
.\.venv\Scripts\python.exe -u main.py --max-carbon 8 --workers 4 --include-stereo
```

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
The optional PDF export uses RDKit, ReportLab, and Pillow. The SVG export uses RDKit.
## Important Note
This program conducts a basic verification of the input molecule. It is not suitable for advanced chemical evaluations or in-depth structural analysis.

**Disclaimer**: The use of this program is for reference purposes only. It does not guarantee chemical accuracy.
