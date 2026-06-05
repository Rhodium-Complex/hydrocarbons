"""Pipeline orchestration for hydrocarbon structure generation."""
from concurrent.futures import ProcessPoolExecutor
from contextlib import nullcontext
from dataclasses import dataclass
import itertools
import time
from collections.abc import Callable, Iterable, Iterator
from typing import TypeVar

import converter
import molecule
import molecule_transformations
import structure_generator

MoleculeGroup = list[molecule.Molecule]
MoleculeGroups = list[MoleculeGroup]
FORMULA_SEPARATOR = "N#N"
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


@dataclass(frozen=True)
class FormulaSmilesGroup:
    """SMILES strings grouped by molecular formula."""

    label: str
    carbon_count: int
    hydrogen_count: int
    smiles: list[str]


def count_structures(structure_groups: MoleculeGroups) -> int:
    """Return the number of molecules stored in grouped structure lists."""
    return sum(len(structures) for structures in structure_groups)


def format_step_result(result: GenerationStepResult) -> str:
    """Format a generation step result using the existing CLI log format."""
    return (
        f"C={result.carbon_count:>2} H={result.hydrogen_count:>2} "
        f"dehydro={result.dehydro_seconds:>9.4f}s "
        f"build={result.build_seconds:>9.4f}s "
        f"smiles={result.smiles_seconds:>9.4f}s "
        f"count={result.count:>7} "
        f"total={result.total_seconds:>9.4f}s"
    )


def _map_generation_task(
    executor: ProcessPoolExecutor | None,
    function: Callable[[TaskInput], TaskOutput],
    iterable: Iterable[TaskInput],
) -> Iterator[TaskOutput]:
    """Map a generation task, using direct iteration for sequential runs."""
    if executor is None:
        return map(function, iterable)
    return executor.map(function, iterable)


def _append_formula_separator(results: list[str]) -> None:
    """Append the external formula-group marker used by downstream exports."""
    results.append(FORMULA_SEPARATOR)


def format_formula_label(carbon_count: int, hydrogen_count: int) -> str:
    """Return a compact hydrocarbon formula label."""
    if carbon_count == 1:
        return f"CH{hydrogen_count}"
    return f"C{carbon_count}H{hydrogen_count}"


def _executor_context(workers: int | None):
    if workers == 1:
        return nullcontext(None)
    return ProcessPoolExecutor(max_workers=workers)


def run_generation(
    min_carbon: int,
    max_carbon: int,
    workers: int | None = None,
    include_smiles: bool = True,
    log_step: Callable[[GenerationStepResult], None] | None = None,
) -> list[str]:
    """Generate hydrocarbon structures and optionally return their SMILES strings."""
    with (
        _executor_context(workers) as dehydro_executor,
        _executor_context(workers) as structure_executor,
    ):
        all_smiles_results = (
            [FORMULA_SEPARATOR, FORMULA_SEPARATOR, FORMULA_SEPARATOR, "C"]
            if include_smiles
            else []
        )

        for carbon_count in range(min_carbon, max_carbon + 1):
            current_carbon_structures: MoleculeGroups = []
            for hydrogen_count in range(0, carbon_count * 2 + 3, 2)[::-1]:
                step_start = time.perf_counter()
                if include_smiles:
                    _append_formula_separator(all_smiles_results)

                dehydro_start = time.perf_counter()
                current_carbon_structures = [
                    structures for structures in current_carbon_structures if structures
                ]
                future_dehydro = _map_generation_task(
                    dehydro_executor,
                    molecule_transformations.unique_dehydro_mols,
                    current_carbon_structures,
                )
                current_carbon_structures = list(future_dehydro)
                dehydro_seconds = time.perf_counter() - dehydro_start

                build_start = time.perf_counter()
                future_structure = _map_generation_task(
                    structure_executor,
                    structure_generator.build_structure,
                    structure_generator.build_carbon_hydrogen_combination(
                        carbon_count,
                        hydrogen_count,
                    ),
                )

                for structures in future_structure:
                    current_carbon_structures += structures
                build_seconds = time.perf_counter() - build_start

                structure_count = count_structures(current_carbon_structures)
                smiles_seconds = 0.0
                output_count = structure_count

                if include_smiles:
                    smiles_start = time.perf_counter()
                    flattened_structures = itertools.chain.from_iterable(
                        current_carbon_structures
                    )
                    future_smiles = map(converter.mat2smiles, flattened_structures)
                    results_before_adding = len(all_smiles_results)
                    all_smiles_results += list(future_smiles)
                    output_count = len(all_smiles_results) - results_before_adding
                    smiles_seconds = time.perf_counter() - smiles_start

                step_result = GenerationStepResult(
                    carbon_count=carbon_count,
                    hydrogen_count=hydrogen_count,
                    dehydro_seconds=dehydro_seconds,
                    build_seconds=build_seconds,
                    smiles_seconds=smiles_seconds,
                    count=output_count,
                    total_seconds=time.perf_counter() - step_start,
                )
                if log_step is not None:
                    log_step(step_result)
    return all_smiles_results


def run_generation_smiles_groups(
    min_carbon: int,
    max_carbon: int,
    workers: int | None = None,
    include_methane: bool = True,
    log_step: Callable[[GenerationStepResult], None] | None = None,
) -> list[FormulaSmilesGroup]:
    """Generate SMILES strings grouped by molecular formula."""
    formula_groups = []
    if include_methane:
        formula_groups.append(
            FormulaSmilesGroup(
                label=format_formula_label(1, 4),
                carbon_count=1,
                hydrogen_count=4,
                smiles=["C"],
            )
        )

    with (
        _executor_context(workers) as dehydro_executor,
        _executor_context(workers) as structure_executor,
    ):
        for carbon_count in range(min_carbon, max_carbon + 1):
            current_carbon_structures: MoleculeGroups = []
            for hydrogen_count in range(0, carbon_count * 2 + 3, 2)[::-1]:
                step_start = time.perf_counter()

                dehydro_start = time.perf_counter()
                current_carbon_structures = [
                    structures for structures in current_carbon_structures if structures
                ]
                future_dehydro = _map_generation_task(
                    dehydro_executor,
                    molecule_transformations.unique_dehydro_mols,
                    current_carbon_structures,
                )
                current_carbon_structures = list(future_dehydro)
                dehydro_seconds = time.perf_counter() - dehydro_start

                build_start = time.perf_counter()
                future_structure = _map_generation_task(
                    structure_executor,
                    structure_generator.build_structure,
                    structure_generator.build_carbon_hydrogen_combination(
                        carbon_count,
                        hydrogen_count,
                    ),
                )

                for structures in future_structure:
                    current_carbon_structures += structures
                build_seconds = time.perf_counter() - build_start

                smiles_start = time.perf_counter()
                flattened_structures = itertools.chain.from_iterable(
                    current_carbon_structures
                )
                smiles = list(map(converter.mat2smiles, flattened_structures))
                smiles_seconds = time.perf_counter() - smiles_start

                formula_groups.append(
                    FormulaSmilesGroup(
                        label=format_formula_label(carbon_count, hydrogen_count),
                        carbon_count=carbon_count,
                        hydrogen_count=hydrogen_count,
                        smiles=smiles,
                    )
                )

                if log_step is not None:
                    log_step(
                        GenerationStepResult(
                            carbon_count=carbon_count,
                            hydrogen_count=hydrogen_count,
                            dehydro_seconds=dehydro_seconds,
                            build_seconds=build_seconds,
                            smiles_seconds=smiles_seconds,
                            count=len(smiles),
                            total_seconds=time.perf_counter() - step_start,
                        )
                    )
    return formula_groups
