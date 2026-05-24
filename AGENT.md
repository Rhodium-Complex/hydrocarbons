# AGENT.md

## Project Goal

This repository focuses on extremely high-performance hydrocarbon structure generation and comparison.

Performance is a primary requirement.
Readable but slow abstractions are not preferred.

The repository is expected to evolve toward:
- low allocation
- cache-friendly layouts
- compact graph/state representations
- Rust acceleration for hot paths

---

# Optimization Priorities

Always prioritize in this order:

1. Correctness
2. Algorithmic improvement
3. Allocation reduction
4. Copy reduction
5. Memory locality
6. Data-oriented design
7. Numba
8. Cython
9. Rust/PyO3

Do NOT prioritize:
- abstraction purity
- object-oriented refactors
- excessive genericity
- unnecessary vectorization

---

# Profiling Rules

Before optimization:

- run profiling
- identify hot paths
- measure allocation behavior

Preferred tools:
- scalene
- cProfile
- line_profiler

Always include:
- before/after benchmark
- memory impact
- correctness validation

---

# Important Constraints

This codebase is NOT dense linear algebra.

Avoid assumptions that:
- NumPy vectorization is automatically faster
- broadcasting is always beneficial
- temporary arrays are acceptable

Frequent issues:
- ndarray copies
- fancy indexing allocations
- transpose copies
- temporary boolean arrays
- Python object churn

---

# Preferred Data Structures

Prefer:
- contiguous ndarray
- uint8 / uint16 / uint32
- bit-packed representations
- immutable compact states
- index indirection instead of physical reordering

Avoid:
- dtype=object
- nested Python lists
- set of tuples
- deepcopy
- repeated concatenate/vstack

---

# Numba Guidelines

Numba is preferred when:
- loops are numeric
- operations are local
- object-free execution is possible

Prefer:
- @numba.njit(cache=True)
- explicit ndarray inputs
- contiguous arrays

Avoid:
- pandas
- Python dict-heavy logic
- object arrays

---

# Rust Migration Strategy

Rust is expected for ultimate performance-critical paths.

Candidates:
- graph canonicalization
- state deduplication
- hashing
- combinatorial search
- adjacency comparison

Preferred stack:
- Rust
- PyO3
- maturin

Python should remain orchestration-focused.

---

# Correctness Rules

All optimizations must preserve:
- canonical structure uniqueness
- deterministic ordering
- graph equivalence behavior

Add tests before optimization whenever possible.

---

# Benchmark Rules

Benchmarks must:
- exclude JIT warmup
- use realistic workloads
- measure both CPU and memory
- include regression comparison

---

# Agent Behavior

Agents should:
- propose multiple optimization approaches
- explain tradeoffs
- avoid speculative rewrites
- avoid full rewrites unless justified by profiling

Never perform broad refactors without measurable benefit.