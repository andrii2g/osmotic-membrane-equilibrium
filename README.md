# Osmotic Membrane Equilibrium

A small Python model of water redistribution between two rigid, equal-area
compartments separated by a membrane that passes water but retains solute.
Requires Python 3.10 or later; the runtime has no third-party dependencies.

All model inputs use SI units. Solute amounts count moles of dissolved particles.
For a non-dissociating solute, these are ordinary solute moles.

```python
from osmotic_membrane_equilibrium import System, solve

result = solve(System())
print(result.status)
print(result.final.height_left_m, result.final.height_right_m)
print(result.final.residual_pa)
```

The conserved total volume is `V = V_left_initial + V_right_initial`. Bisection
finds the final left volume `v` from:

```text
rho*g*(2*v - V)/A - R*T*(n_left/v - n_right/(V - v)) = 0
```

The residual is hydrostatic pressure difference minus osmotic pressure
difference, with both differences defined as left minus right. Its derivative is
positive, so a wet equilibrium is unique when it exists. Both compartments with
positive solute amounts guarantee an interior root.

The model assumes ideal dilute, well-mixed solutions at a common fixed
temperature, incompressible liquid with a common constant density, negligible
solute volume, equal atmospheric pressure, rigid compartments of unlimited
height, and a membrane at a common bottom reference level. It excludes membrane
deformation, solute transport, nonideal activities, and flow dynamics.

With a solute-free side, the wet equilibrium can lie beyond the available water.
The solver then reports `dry_boundary`, not convergence; the reported dry-side
concentration is its limiting value of zero. A nonzero residual can remain there,
and the wet-compartment model stops at that boundary.

Run the tests from a source checkout in PowerShell:

```powershell
$env:PYTHONPATH = "$PWD/src"
python -m unittest discover -s tests -v
```
