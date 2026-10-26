"""Independent parameter sweeps and unit-labelled, flat result records."""

import csv
from dataclasses import asdict, replace
import math
from pathlib import Path
from typing import Iterable

from .model import System
from .solver import Result, solve


def linear_values(start: float, stop: float, points: int) -> list[float]:
    if not math.isfinite(start) or not math.isfinite(stop):
        raise ValueError("sweep endpoints must be finite")
    if not isinstance(points, int) or isinstance(points, bool) or points < 2:
        raise ValueError("sweep points must be an integer of at least two")
    values = [(1 - i / (points - 1)) * start + i / (points - 1) * stop
              for i in range(points)]
    values[0], values[-1] = start, stop
    return values


def result_record(system: System, result: Result) -> dict:
    record = asdict(system)
    for prefix, state in (("initial", result.initial), ("final", result.final)):
        record.update({f"{prefix}_{name}": value for name, value in asdict(state).items()})
        record[f"{prefix}_concentration_difference_mol_m3"] = state.concentration_difference_mol_m3
    record.update(
        status=result.status,
        converged=result.converged,
        iterations=result.iterations,
        bracket_width_m3=result.bracket_width_m3,
        transferred_to_left_m3=result.transferred_to_left_m3,
        volume_error_m3=result.volume_error_m3,
    )
    return record


def sweep_temperature(system: System, temperatures_k: Iterable[float], **solver_options) -> list[dict]:
    """Hold geometry, initial volumes, and solute amounts fixed."""
    systems = [replace(system, temperature_k=value) for value in temperatures_k]
    if not systems:
        raise ValueError("a sweep must contain at least one value")
    return [dict(sweep_temperature_k=item.temperature_k,
                 **result_record(item, solve(item, **solver_options))) for item in systems]


def sweep_concentration(
    system: System, differences_mol_m3: Iterable[float], **solver_options,
) -> list[dict]:
    """Vary initial RIGHT-minus-LEFT concentration at fixed arithmetic mean.

    Initial volumes, geometry, and temperature stay fixed. Solute amounts change
    between experiments, but each experiment conserves its own solute amounts.
    """
    mean = (system.solute_left_mol / system.volume_left_m3
            + system.solute_right_mol / system.volume_right_m3) / 2
    if not math.isfinite(mean):
        raise ValueError("mean initial concentration must be finite")
    cases = []
    for difference in differences_mol_m3:
        if not math.isfinite(difference) or abs(difference) > 2 * mean:
            raise ValueError("concentration difference must be finite and within +/- twice the mean")
        item = replace(system,
                       solute_left_mol=(mean - difference / 2) * system.volume_left_m3,
                       solute_right_mol=(mean + difference / 2) * system.volume_right_m3)
        cases.append((difference, item))
    if not cases:
        raise ValueError("a sweep must contain at least one value")
    return [dict(sweep_initial_right_minus_left_mol_m3=difference,
                 **result_record(item, solve(item, **solver_options)))
            for difference, item in cases]


def write_csv(path: str | Path, records: list[dict]) -> None:
    if not records:
        raise ValueError("cannot write an empty result table")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
