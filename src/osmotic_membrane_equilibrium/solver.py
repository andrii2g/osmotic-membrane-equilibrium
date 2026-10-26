"""Bracketed equilibrium solver, including dry-compartment limits."""

from dataclasses import dataclass
import math
from typing import Literal

from .model import State, System, evaluate


@dataclass(frozen=True)
class Result:
    initial: State
    final: State
    status: Literal["converged", "dry_boundary", "max_iterations", "precision_limit"]
    iterations: int
    bracket_width_m3: float

    @property
    def converged(self) -> bool:
        return self.status == "converged"

    @property
    def transferred_to_left_m3(self) -> float:
        return self.final.volume_left_m3 - self.initial.volume_left_m3

    @property
    def volume_error_m3(self) -> float:
        return (
            self.final.volume_left_m3 + self.final.volume_right_m3
            - self.initial.volume_left_m3 - self.initial.volume_right_m3
        )


def solve(
    system: System,
    *,
    pressure_tolerance_pa: float = 1e-6,
    volume_tolerance_m3: float = 1e-15,
    relative_volume_tolerance: float = 1e-12,
    max_iterations: int = 200,
) -> Result:
    """Find F = hydrostatic difference - osmotic difference = 0 by bisection.

    Success requires pressure tolerance AND volume bracket tolerance, unless
    an exactly zero residual is found. Boundary outcomes are not wet equilibria.
    """
    for name, value in (
        ("pressure_tolerance_pa", pressure_tolerance_pa),
        ("volume_tolerance_m3", volume_tolerance_m3),
        ("relative_volume_tolerance", relative_volume_tolerance),
    ):
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"{name} must be finite and nonnegative")
    if not isinstance(max_iterations, int) or isinstance(max_iterations, bool) or max_iterations < 1:
        raise ValueError("max_iterations must be a positive integer")

    initial = evaluate(system, system.volume_left_m3)
    total = system.total_volume_m3
    low, high = 0.0, total
    tolerance = volume_tolerance_m3 + relative_volume_tolerance * total
    if not math.isfinite(tolerance):
        raise ValueError("combined volume tolerance must be finite")

    if system.solute_left_mol == system.solute_right_mol == system.gravity_m_s2 == 0:
        return Result(initial, initial, "converged", 0, 0.0)
    if system.solute_left_mol == 0:
        boundary = evaluate(system, low)
        if boundary.residual_pa >= 0:
            return Result(initial, boundary, "dry_boundary", 0, 0.0)
    if system.solute_right_mol == 0:
        boundary = evaluate(system, high)
        if boundary.residual_pa <= 0:
            return Result(initial, boundary, "dry_boundary", 0, 0.0)

    best = initial
    for iteration in range(1, max_iterations + 1):
        middle = low + (high - low) / 2
        if middle == low or middle == high:
            return Result(initial, best, "precision_limit", iteration - 1, high - low)
        state = evaluate(system, middle)
        if abs(state.residual_pa) < abs(best.residual_pa):
            best = state
        if state.residual_pa == 0:
            return Result(initial, state, "converged", iteration, 0.0)
        if state.residual_pa < 0:
            low = middle
        else:
            high = middle
        if abs(state.residual_pa) <= pressure_tolerance_pa and high - low <= tolerance:
            return Result(initial, state, "converged", iteration, high - low)
    return Result(initial, best, "max_iterations", max_iterations, high - low)
