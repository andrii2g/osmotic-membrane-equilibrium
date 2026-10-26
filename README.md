# Osmotic Membrane Equilibrium

A small Python model of water redistribution between two rigid, equal-area
compartments separated by a membrane that passes water but retains solute.
Requires Python 3.10 or later; the runtime has no third-party dependencies.

Install and run from PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e .
.\.venv\Scripts\osmotic-equilibrium.exe solve
.\.venv\Scripts\python.exe -m osmotic_membrane_equilibrium solve --json
```

The default experiment starts with 1000 mL on each side, 0.1 mmol of solute
particles on the left and 0.2 mmol on the right, 100 cm² area per compartment,
and a temperature of 298.15 K. Customize the experiment with unit-labelled flags:

```powershell
python -m osmotic_membrane_equilibrium solve --left-volume-ml 750 --right-volume-ml 1250 --temperature-k 310
python -m osmotic_membrane_equilibrium solve --json --output artifacts/solution.csv
python -m osmotic_membrane_equilibrium sweep-temperature --start-k 278.15 --stop-k 338.15 --points 13 --output artifacts/temperature.csv
python -m osmotic_membrane_equilibrium sweep-concentration --start-mol-m3 -0.2 --stop-mol-m3 0.2 --points 21 --output artifacts/concentration.csv
```

Use the virtual environment's Python, or activate it, for these commands.
Each subcommand supports `--help`. CSV and JSON contain all initial and final
state values, individual osmotic pressures, signed pressure and concentration
differences, volume error, transfer, solver iterations, final bracket width, and
pressure residual. Column names include their SI units. One mol/m³ equals one
mmol/L; one mol/L equals 1000 mol/m³.

Temperature sweeps hold solute amounts, initial volumes, and geometry fixed.
Concentration sweeps hold the arithmetic mean of initial concentrations fixed
and vary **right minus left** concentration. This sweep coordinate has the
opposite sign to the left-minus-right state differences. Concentrations must
remain nonnegative, so the sweep difference cannot exceed twice the mean in
absolute value. Different sweep experiments can have different solute amounts;
each individual experiment conserves its own amounts.

Exit codes: `0` for convergence, `1` for solver exhaustion or floating-point
precision limits, `2` for invalid arguments or file errors, and `3` for a dry
boundary. Sweeps return `1` if any numerical failure occurs, otherwise `3` if
any experiment reaches a dry boundary. Result files include these statuses.

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

Run the tests from a source checkout in PowerShell (or omit `PYTHONPATH` after
installing the package):

```powershell
$env:PYTHONPATH = "$PWD/src"
python -m unittest discover -s tests -v
```
