"""Ideal dilute osmosis between two rigid, equal-area compartments."""

from .model import GAS_CONSTANT, State, System, evaluate
from .solver import Result, solve

__all__ = ["GAS_CONSTANT", "State", "System", "evaluate", "Result", "solve"]
