"""Pipeline orchestration for hydrocarbon structure generation."""
from concurrent.futures import ProcessPoolExecutor
from collections import deque
from contextlib import ExitStack, nullcontext
from dataclasses import dataclass
import itertools
import os
import time
from collections.abc import Callable, Iterable, Iterator
from typing import TypeVar

import converter
from generation_output import FormulaSmilesGroup, StructureVariant, format_formula_label
import molecule
import molecule_transformations
import structure_generator

MoleculeGroup = list[molecule.Molecule]
MoleculeGroups = list[MoleculeGroup]
MIN_PARALLEL_STEREO_STRUCTURES = 256
TaskInput = TypeVar("TaskInput")
TaskOutput = TypeVar("TaskOutput")


@dataclass(frozen=True)
class GenerationStepResult:
    """Timing and count information for one C/H generation step."""

    carbon_count: int
    hydrogen_count: int
    dehydro_seconds: float
    build_seconds: float
    smiles_seconds: float
    count: int
    total_seconds: float
    smiles_fused: bool = False


@dataclass(frozen=True)
class _FormulaStructureStep:
    """Generated structures and timings for one molecular formula."""

    carbon_count: int
    hydrogen_count: int
    structures: MoleculeGroups
    dehydro_seconds: float
    build_seconds: float
    outputs: list[StructureVariant] | None = None


def count_structures(structure_groups: MoleculeGroups) -> int:
    """Return the number of molecules stored in grouped structure lists."""
    return sum(len(structures) for structures in structure_groups)


def format_step_result(result: GenerationStepResult) -> str:
    """Format a generation step result using the existing CLI log format."""
    smiles_timing = "fused" if result.smiles_fused else f"{result.smiles_seconds:>9.4f}s"
    return (
        f"C={result.carbon_count:>2} H={result.hydrogen_count:>2} "
        f"dehydro={result.dehydro_seconds:>9.4f}s "
        f"build={result.build_seconds:>9.4f}s "
        f"smiles={smiles_timing} "
        f"count={result.count:>7} "
        f"total={result.total_seconds:>9.4f}s"
    )


def print_step_result(result: GenerationStepResult) -> None:
    """Print a generation step result immediately."""
    print(format_step_result(result), flush=True)


def _map_generation_task(
    executor: ProcessPoolExecutor | None,
    function: Callable[[TaskInput], TaskOutput],
    iterable: Iterable[TaskInput],
    chunksize: int = 1,
    workers: int | None = None,
) -> Iterator[TaskOutput]:
    """Map in input order with at most two uncollected batches per worker."""
    if executor is None:
        yield from map(function, iterable)
        return
    pending = deque()
    items = iter(iterable)
    limit = _worker_count(workers) * 2
    try:
        while True:
            while len(pending) < limit:
                batch = tuple(itertools.islice(items, chunksize))
                if not batch:
                    break
                pending.append(executor.submit(_run_task_batch, function, batch))
            if not pending:
                return
            yield from pending.popleft().result()
    finally:
        for future in pending:
            future.cancel()


def _run_task_batch(function, batch):
    return [function(item) for item in batch]


def _worker_count(workers: int | None) -> int:
    if workers is not None:
        return workers
    cpu_count = getattr(os, "process_cpu_count", os.cpu_count)() or 1
    return min(cpu_count, 61) if os.name == "nt" else cpu_count


def _dehydro_chunksize(group_count: int, workers: int | None) -> int:
    target_batches = _worker_count(workers) * 8
    return max(1, min(64, (group_count + target_batches - 1) // target_batches))


def _dehydro_task(task):
    """Convert only the survivors of this group's structural deduplication."""
    group, fuse_smiles = task
    structures = molecule_transformations.unique_dehydro_mols(group)
    outputs = _plain_outputs([structures]) if fuse_smiles else None
    return structures, outputs


def _build_task(task):
    """Return new skeleton groups and, optionally, their ordinary SMILES."""
    pattern, fuse_smiles = task
    structures = structure_generator.build_structure(pattern)
    outputs = _plain_outputs(structures) if fuse_smiles else None
    return structures, outputs


def _plain_outputs(structures: MoleculeGroups) -> list[StructureVariant]:
    return [
        variant
        for group in structures
        for item in group
        for variant in converter.mat2structure_variants(item)
    ]


def _executor_context(workers: int | None):
    if workers == 1:
        return nullcontext(None)
    return ProcessPoolExecutor(max_workers=workers)


def _iter_formula_structure_steps(
    min_carbon: int,
    max_carbon: int,
    workers: int | None = None,
    fuse_smiles: bool = False,
) -> Iterator[_FormulaStructureStep]:
    """Yield generated structures for each formula in the existing order."""
    with _executor_context(workers) as executor:
        for carbon_count in range(min_carbon, max_carbon + 1):
            current_carbon_structures: MoleculeGroups = []
            for hydrogen_count in range(0, carbon_count * 2 + 3, 2)[::-1]:
                dehydro_start = time.perf_counter()
                current_carbon_structures = [
                    structures for structures in current_carbon_structures if structures
                ]
                future_dehydro = _map_generation_task(
                    executor,
                    _dehydro_task,
                    ((group, fuse_smiles) for group in current_carbon_structures),
                    chunksize=_dehydro_chunksize(len(current_carbon_structures), workers),
                    workers=workers,
                )
                next_structures = []
                outputs = [] if fuse_smiles else None
                for structures, variants in future_dehydro:
                    next_structures.append(structures)
                    if outputs is not None:
                        outputs.extend(variants)
                current_carbon_structures = next_structures
                dehydro_seconds = time.perf_counter() - dehydro_start

                build_start = time.perf_counter()
                future_structure = _map_generation_task(
                    executor,
                    _build_task,
                    (
                        (pattern, fuse_smiles)
                        for pattern in structure_generator.build_carbon_hydrogen_combination(
                            carbon_count, hydrogen_count,
                        )
                    ),
                    workers=workers,
                )

                for structures, variants in future_structure:
                    current_carbon_structures += structures
                    if outputs is not None:
                        outputs.extend(variants)
                build_seconds = time.perf_counter() - build_start

                yield _FormulaStructureStep(
                    carbon_count=carbon_count,
                    hydrogen_count=hydrogen_count,
                    structures=current_carbon_structures,
                    dehydro_seconds=dehydro_seconds,
                    build_seconds=build_seconds,
                    outputs=outputs,
                )


def _structure_variants_task(task: tuple[molecule.Molecule, bool, bool]):
    molecule_obj, include_stereo, include_tetrahedral_stereo = task
    return converter.mat2structure_variants(
        molecule_obj, include_stereo, include_tetrahedral_stereo
    )


def _smiles_chunksize(structure_count: int, workers: int | None) -> int:
    """Choose a coarse chunksize to reduce multiprocessing dispatch overhead."""
    worker_count = workers if workers is not None else (os.cpu_count() or 1)
    if worker_count <= 1:
        return 1
    return max(1, structure_count // (worker_count * 8))


def _should_parallelize_stereo_smiles(
    structure_count: int,
    include_stereo: bool,
    include_tetrahedral_stereo: bool,
    workers: int | None,
) -> bool:
    """Return whether stereo SMILES conversion is large enough to parallelize."""
    return (
        (include_stereo or include_tetrahedral_stereo)
        and workers != 1
        and structure_count >= MIN_PARALLEL_STEREO_STRUCTURES
    )


def _structure_outputs(
    structure_groups: MoleculeGroups,
    include_stereo: bool,
    include_tetrahedral_stereo: bool = False,
    executor: ProcessPoolExecutor | None = None,
    workers: int | None = None,
):
    structures = list(itertools.chain.from_iterable(structure_groups))
    tasks = (
        (item, include_stereo, include_tetrahedral_stereo) for item in structures
    )
    if executor is None or len(structures) < MIN_PARALLEL_STEREO_STRUCTURES:
        nested = (_structure_variants_task(task) for task in tasks)
    else:
        nested = _map_generation_task(
            executor,
            _structure_variants_task,
            tasks,
            chunksize=_smiles_chunksize(len(structures), workers),
            workers=workers,
        )
    return list(itertools.chain.from_iterable(nested))


def _run_generation_pipeline(
    min_carbon: int,
    max_carbon: int,
    workers: int | None = None,
    include_smiles: bool = False,
    include_methane: bool = True,
    include_stereo: bool = False,
    include_tetrahedral_stereo: bool = False,
    log_step: Callable[[GenerationStepResult], None] | None = None,
) -> list[FormulaSmilesGroup] | None:
    """Run the shared generation loop, optionally collecting SMILES groups."""
    formula_groups = []
    if include_smiles and include_methane:
        formula_groups.append(
            FormulaSmilesGroup(
                label=format_formula_label(1, 4),
                carbon_count=1,
                hydrogen_count=4,
                variants=[StructureVariant("C")],
            )
        )

    with ExitStack() as stack:
        smiles_executor = None
        fuse_smiles = include_smiles and not (include_stereo or include_tetrahedral_stereo)
        for step in _iter_formula_structure_steps(
            min_carbon, max_carbon, workers, fuse_smiles=fuse_smiles,
        ):
            structure_count = count_structures(step.structures)
            smiles_seconds = 0.0
            output_count = structure_count

            if include_smiles and not fuse_smiles:
                smiles_start = time.perf_counter()
                if (
                    smiles_executor is None
                    and _should_parallelize_stereo_smiles(
                        structure_count,
                        include_stereo,
                        include_tetrahedral_stereo,
                        workers,
                    )
                ):
                    smiles_executor = stack.enter_context(_executor_context(workers))
                outputs = _structure_outputs(
                    step.structures,
                    include_stereo,
                    include_tetrahedral_stereo,
                    executor=smiles_executor,
                    workers=workers,
                )
                smiles_seconds = time.perf_counter() - smiles_start
            elif fuse_smiles:
                outputs = step.outputs

            if include_smiles:
                output_count = len(outputs)
                formula_groups.append(
                    FormulaSmilesGroup(
                        label=format_formula_label(
                            step.carbon_count,
                            step.hydrogen_count,
                        ),
                        carbon_count=step.carbon_count,
                        hydrogen_count=step.hydrogen_count,
                        variants=outputs,
                    )
                )

            if log_step is not None:
                log_step(
                    GenerationStepResult(
                        carbon_count=step.carbon_count,
                        hydrogen_count=step.hydrogen_count,
                        dehydro_seconds=step.dehydro_seconds,
                        build_seconds=step.build_seconds,
                        smiles_seconds=smiles_seconds,
                        count=output_count,
                        total_seconds=(
                            step.dehydro_seconds + step.build_seconds + smiles_seconds
                        ),
                        smiles_fused=fuse_smiles,
                    )
                )

    if include_smiles:
        return formula_groups
    return None


def run_generation(
    min_carbon: int,
    max_carbon: int,
    workers: int | None = None,
    log_step: Callable[[GenerationStepResult], None] | None = None,
) -> None:
    """Generate hydrocarbon structures and optionally log structure counts."""
    _run_generation_pipeline(
        min_carbon=min_carbon,
        max_carbon=max_carbon,
        workers=workers,
        include_smiles=False,
        log_step=log_step,
    )


def run_generation_smiles_groups(
    min_carbon: int,
    max_carbon: int,
    workers: int | None = None,
    include_methane: bool = True,
    include_stereo: bool = False,
    log_step: Callable[[GenerationStepResult], None] | None = None,
    include_tetrahedral_stereo: bool = False,
) -> list[FormulaSmilesGroup]:
    """Generate SMILES strings grouped by molecular formula."""
    formula_groups = _run_generation_pipeline(
        min_carbon=min_carbon,
        max_carbon=max_carbon,
        workers=workers,
        include_smiles=True,
        include_methane=include_methane,
        include_stereo=include_stereo,
        include_tetrahedral_stereo=include_tetrahedral_stereo,
        log_step=log_step,
    )
    if formula_groups is None:
        return []
    return formula_groups


def run_export_smiles_groups(
    min_carbon: int,
    max_carbon: int,
    workers: int | None = None,
    include_stereo: bool = False,
    include_tetrahedral_stereo: bool = False,
) -> list[FormulaSmilesGroup]:
    """Generate formula-grouped SMILES with standard export logging."""
    return run_generation_smiles_groups(
        min_carbon=min_carbon,
        max_carbon=max_carbon,
        workers=workers,
        include_methane=True,
        include_stereo=include_stereo,
        include_tetrahedral_stereo=include_tetrahedral_stereo,
        log_step=print_step_result,
    )
