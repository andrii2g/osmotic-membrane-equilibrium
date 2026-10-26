import math
import unittest
from dataclasses import replace

from osmotic_membrane_equilibrium import GAS_CONSTANT, System, evaluate, solve


class ModelTests(unittest.TestCase):
    def test_pressure_units(self):
        system = System()
        state = evaluate(system, system.volume_left_m3)
        self.assertAlmostEqual(state.height_left_m, 0.1)
        self.assertAlmostEqual(state.concentration_left_mol_m3, 0.1)
        self.assertAlmostEqual(state.osmotic_left_pa, 0.1 * GAS_CONSTANT * 298.15)
        self.assertAlmostEqual(state.hydrostatic_difference_pa, 0)

    def test_invalid_inputs(self):
        for field, value in (
            ("volume_left_m3", 0), ("volume_right_m3", -1),
            ("solute_left_mol", -1), ("solute_right_mol", math.nan),
            ("temperature_k", 0), ("area_m2", 0),
            ("density_kg_m3", -1), ("gravity_m_s2", -1),
            ("temperature_k", math.inf),
        ):
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                replace(System(), **{field: value})
        for volume in (-1, 1, math.nan, 0, System().total_volume_m3):
            with self.subTest(volume=volume), self.assertRaises(ValueError):
                evaluate(System(), volume)


class SolverTests(unittest.TestCase):
    def assert_balanced(self, system, result):
        self.assertTrue(result.converged, result.status)
        state = result.final
        self.assertLessEqual(abs(state.residual_pa), 1e-6)
        self.assertAlmostEqual(state.volume_left_m3 + state.volume_right_m3,
                               system.total_volume_m3, delta=1e-18)
        self.assertAlmostEqual(state.concentration_left_mol_m3 * state.volume_left_m3,
                               system.solute_left_mol, delta=1e-18)
        self.assertAlmostEqual(state.concentration_right_mol_m3 * state.volume_right_m3,
                               system.solute_right_mol, delta=1e-18)

    def test_default_conserves_volume_and_solute_and_balances_pressure(self):
        system = System()
        result = solve(system)
        self.assert_balanced(system, result)
        self.assertLess(result.transferred_to_left_m3, 0)
        self.assertGreater(result.final.concentration_right_mol_m3,
                           result.final.concentration_left_mol_m3)
        self.assertLess(result.final.concentration_right_mol_m3
                        - result.final.concentration_left_mol_m3, 0.1)
        self.assertLessEqual(result.bracket_width_m3, 3e-15)

    def test_symmetric_solutions(self):
        system = System(solute_right_mol=0.0001)
        result = solve(system)
        self.assert_balanced(system, result)
        self.assertEqual(result.transferred_to_left_m3, 0)

    def test_pure_water_equalizes_unequal_heights(self):
        system = System(volume_left_m3=0.0005, volume_right_m3=0.0015,
                        solute_left_mol=0, solute_right_mol=0)
        result = solve(system)
        self.assert_balanced(system, result)
        self.assertEqual(result.final.height_left_m, result.final.height_right_m)

    def test_zero_gravity_matches_analytical_solution(self):
        for left, right in ((1e-4, 2e-4), (1e-6, 8e-4), (3e-4, 1e-4)):
            system = System(solute_left_mol=left, solute_right_mol=right, gravity_m_s2=0)
            result = solve(system)
            self.assert_balanced(system, result)
            expected = system.total_volume_m3 * left / (left + right)
            self.assertAlmostEqual(result.final.volume_left_m3, expected, delta=3e-15)

    def test_swapping_compartments_swaps_results(self):
        system = System(volume_left_m3=0.0007, volume_right_m3=0.0013)
        swapped = replace(system, volume_left_m3=system.volume_right_m3,
                          volume_right_m3=system.volume_left_m3,
                          solute_left_mol=system.solute_right_mol,
                          solute_right_mol=system.solute_left_mol)
        a, b = solve(system), solve(swapped)
        self.assert_balanced(system, a)
        self.assert_balanced(swapped, b)
        self.assertAlmostEqual(a.final.volume_left_m3, b.final.volume_right_m3, delta=3e-15)
        self.assertAlmostEqual(a.transferred_to_left_m3, -b.transferred_to_left_m3, delta=3e-15)

    def test_solute_free_side_can_have_a_wet_equilibrium(self):
        system = System(solute_left_mol=0)
        self.assert_balanced(system, solve(system))

    def test_dry_boundary_is_reported_not_claimed_as_pressure_balance(self):
        for left, right in ((0, 0.01), (0.01, 0)):
            result = solve(System(solute_left_mol=left, solute_right_mol=right))
            self.assertEqual(result.status, "dry_boundary")
            self.assertFalse(result.converged)
            self.assertEqual(min(result.final.volume_left_m3, result.final.volume_right_m3), 0)
            self.assertGreater(abs(result.final.residual_pa), 1)

    def test_zero_gravity_pure_water_preserves_initial_distribution(self):
        system = System(volume_left_m3=0.0007, volume_right_m3=0.0013,
                        solute_left_mol=0, solute_right_mol=0, gravity_m_s2=0)
        result = solve(system)
        self.assert_balanced(system, result)
        self.assertEqual(result.transferred_to_left_m3, 0)

    def test_iteration_exhaustion(self):
        result = solve(System(), max_iterations=1)
        self.assertEqual(result.status, "max_iterations")
        self.assertFalse(result.converged)
        self.assertEqual(result.iterations, 1)

    def test_unattainable_tolerance_reports_precision_limit(self):
        result = solve(System(), pressure_tolerance_pa=0,
                       volume_tolerance_m3=0, relative_volume_tolerance=0)
        self.assertIn(result.status, ("precision_limit", "converged"))
        if result.converged:
            self.assertEqual(result.final.residual_pa, 0)

    def test_invalid_solver_options(self):
        for options in ({"max_iterations": 0}, {"max_iterations": 1.5},
                        {"pressure_tolerance_pa": -1},
                        {"volume_tolerance_m3": math.nan},
                        {"relative_volume_tolerance": math.inf}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                solve(System(), **options)


if __name__ == "__main__":
    unittest.main()
