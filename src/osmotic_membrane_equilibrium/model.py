"""SI-unit model: equal areas, common density, fixed solute, open free surfaces."""

from dataclasses import dataclass
import math

GAS_CONSTANT = 8.31446261815324  # J / (mol K)


@dataclass(frozen=True)
class System:
    """Amounts are moles of osmotically active particles, not formula units."""

    volume_left_m3: float = 0.001
    volume_right_m3: float = 0.001
    solute_left_mol: float = 0.0001
    solute_right_mol: float = 0.0002
    temperature_k: float = 298.15
    area_m2: float = 0.01
    density_kg_m3: float = 1000.0
    gravity_m_s2: float = 9.80665

    def __post_init__(self) -> None:
        for name, value in vars(self).items():
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            if name in {"solute_left_mol", "solute_right_mol", "gravity_m_s2"}:
                if value < 0:
                    raise ValueError(f"{name} must be nonnegative")
            elif value <= 0:
                raise ValueError(f"{name} must be positive")
        if not math.isfinite(self.total_volume_m3):
            raise ValueError("total volume must be finite")

    @property
    def total_volume_m3(self) -> float:
        return self.volume_left_m3 + self.volume_right_m3


@dataclass(frozen=True)
class State:
    volume_left_m3: float
    volume_right_m3: float
    height_left_m: float
    height_right_m: float
    concentration_left_mol_m3: float
    concentration_right_mol_m3: float
    osmotic_left_pa: float
    osmotic_right_pa: float
    hydrostatic_difference_pa: float
    osmotic_difference_pa: float
    residual_pa: float

    @property
    def concentration_difference_mol_m3(self) -> float:
        return self.concentration_left_mol_m3 - self.concentration_right_mol_m3


def evaluate(system: System, volume_left_m3: float) -> State:
    """Evaluate left-minus-right pressures at a conserved liquid distribution.

    A dry, solute-free side has concentration zero by its limiting value.
    A dry side containing solute is outside this model's domain.
    """
    if not math.isfinite(volume_left_m3):
        raise ValueError("left volume must be finite")
    total = system.total_volume_m3
    if not 0 <= volume_left_m3 <= total:
        raise ValueError("left volume must be between zero and total volume")
    volume_right = total - volume_left_m3
    concentrations = []
    for amount, volume in (
        (system.solute_left_mol, volume_left_m3),
        (system.solute_right_mol, volume_right),
    ):
        if volume == 0 and amount > 0:
            raise ValueError("a compartment containing solute cannot be dry")
        concentrations.append(amount / volume if volume else 0.0)
    c_left, c_right = concentrations
    h_left = volume_left_m3 / system.area_m2
    h_right = volume_right / system.area_m2
    rt = GAS_CONSTANT * system.temperature_k
    pi_left, pi_right = rt * c_left, rt * c_right
    hydrostatic = system.density_kg_m3 * system.gravity_m_s2 * (h_left - h_right)
    osmotic = pi_left - pi_right
    state = State(
        volume_left_m3, volume_right, h_left, h_right, c_left, c_right,
        pi_left, pi_right, hydrostatic, osmotic, hydrostatic - osmotic,
    )
    if not all(math.isfinite(value) for value in vars(state).values()):
        raise ValueError("derived state exceeds floating-point range; rescale inputs")
    return state
