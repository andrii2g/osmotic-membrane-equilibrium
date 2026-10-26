import csv
import io
import json
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from osmotic_membrane_equilibrium import System
from osmotic_membrane_equilibrium.cli import main
from osmotic_membrane_equilibrium.experiments import (
    linear_values, sweep_concentration, sweep_temperature, write_csv,
)


class ExperimentTests(unittest.TestCase):
    def test_temperature_sweep_holds_amounts_and_increases_height_difference(self):
        system = System()
        records = sweep_temperature(system, [278.15, 298.15, 338.15])
        heights = []
        for row in records:
            self.assertTrue(row["converged"])
            self.assertEqual(row["solute_left_mol"], system.solute_left_mol)
            self.assertEqual(row["solute_right_mol"], system.solute_right_mol)
            self.assertLess(abs(row["volume_error_m3"]), 1e-18)
            self.assertLess(abs(row["final_residual_pa"]), 1e-6)
            heights.append(row["final_height_right_m"] - row["final_height_left_m"])
        self.assertLess(heights[0], heights[1])
        self.assertLess(heights[1], heights[2])

    def test_concentration_sweep_holds_mean_and_changes_direction(self):
        records = sweep_concentration(System(), [-0.2, 0, 0.2])
        for delta, row in zip([-0.2, 0, 0.2], records):
            self.assertTrue(row["converged"])
            self.assertAlmostEqual(row["initial_concentration_left_mol_m3"]
                                   + row["initial_concentration_right_mol_m3"], 0.3)
            self.assertAlmostEqual(row["initial_concentration_right_mol_m3"]
                                   - row["initial_concentration_left_mol_m3"], delta)
        self.assertGreater(records[0]["transferred_to_left_m3"], 0)
        self.assertEqual(records[1]["transferred_to_left_m3"], 0)
        self.assertLess(records[2]["transferred_to_left_m3"], 0)

    def test_invalid_sweep_values(self):
        for function, values in ((sweep_temperature, []), (sweep_temperature, [0]),
                                 (sweep_concentration, []), (sweep_concentration, [0.31]),
                                 (sweep_concentration, [float("nan")])):
            with self.subTest(function=function.__name__, values=values), self.assertRaises(ValueError):
                function(System(), values)
        for values in ((0, 1, 1), (0, float("inf"), 2)):
            with self.assertRaises(ValueError):
                linear_values(*values)

    def test_csv_round_trip_has_metrics_and_exact_sweep_endpoints(self):
        values = linear_values(278.15, 338.15, 13)
        self.assertEqual((values[0], values[-1]), (278.15, 338.15))
        records = sweep_temperature(System(), values)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "nested" / "temperatures.csv"
            write_csv(path, records)
            with path.open(newline="", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 13)
        self.assertEqual(float(rows[0]["sweep_temperature_k"]), 278.15)
        self.assertIn("final_concentration_difference_mol_m3", rows[0])
        self.assertIn("final_osmotic_left_pa", rows[0])


class CliTests(unittest.TestCase):
    def run_main(self, args):
        output = io.StringIO()
        with redirect_stdout(output):
            status = main(args)
        return status, output.getvalue()

    def test_json_and_convenience_unit_conversion(self):
        status, output = self.run_main(["solve", "--left-volume-ml", "500", "--json"])
        row = json.loads(output)
        self.assertEqual(status, 0)
        self.assertEqual(row["volume_left_m3"], 0.0005)
        self.assertEqual(row["area_m2"], 0.01)
        self.assertEqual(row["solute_left_mol"], 0.0001)
        self.assertLess(abs(row["final_residual_pa"]), 1e-6)

    def test_human_report_and_sweep_json(self):
        status, output = self.run_main(["solve"])
        self.assertEqual(status, 0)
        self.assertIn("Initial", output)
        self.assertIn("residual:", output)
        status, output = self.run_main(["sweep-temperature", "--points", "3", "--json"])
        self.assertEqual(status, 0)
        self.assertEqual(len(json.loads(output)), 3)
        status, output = self.run_main(["sweep-concentration", "--points", "3", "--json"])
        self.assertEqual(status, 0)
        self.assertEqual(len(json.loads(output)), 3)

    def test_failure_and_boundary_exit_codes(self):
        status, _ = self.run_main(["solve", "--max-iterations", "1", "--json"])
        self.assertEqual(status, 1)
        status, output = self.run_main(["solve", "--left-solute-mmol", "0",
                                        "--right-solute-mmol", "10", "--json"])
        self.assertEqual(status, 3)
        self.assertEqual(json.loads(output)["status"], "dry_boundary")

    def test_invalid_inputs_have_clean_errors(self):
        with redirect_stderr(io.StringIO()) as output, self.assertRaises(SystemExit) as error:
            main(["solve", "--temperature-k", "0"])
        self.assertEqual(error.exception.code, 2)
        self.assertIn("must be positive", output.getvalue())
        self.assertNotIn("Traceback", output.getvalue())

    def test_module_entry_point(self):
        process = subprocess.run([sys.executable, "-m", "osmotic_membrane_equilibrium",
                                  "solve", "--json"], capture_output=True, text=True)
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertTrue(json.loads(process.stdout)["converged"])


if __name__ == "__main__":
    unittest.main()
