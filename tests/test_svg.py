import io
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
import tempfile
import unittest
from xml.etree import ElementTree as ET

from osmotic_membrane_equilibrium import System, solve
from osmotic_membrane_equilibrium.cli import main
from osmotic_membrane_equilibrium.svg import render_svg

NS = {"s": "http://www.w3.org/2000/svg"}


class SvgTests(unittest.TestCase):
    def test_water_rectangles_share_physical_scale(self):
        system = System()
        result = solve(system)
        root = ET.fromstring(render_svg(system, result))
        rectangles = [node for node in root.findall("s:rect", NS)
                      if node.get("class", "").startswith("water-")]
        self.assertEqual(len(rectangles), 4)
        scales = [float(node.get("height")) / float(node.get("data-height-m"))
                  for node in rectangles]
        self.assertLess(max(scales) - min(scales), 1e-4)
        self.assertAlmostEqual(float(rectangles[2].get("data-height-m")), result.final.height_left_m)
        self.assertEqual(len(root.findall("s:line[@class='membrane']", NS)), 2)
        self.assertEqual(len(root.findall("s:polyline", NS)), 2)
        self.assertIn("common height scale", " ".join(root.itertext()))

    def test_boundary_pure_water_and_failed_solver_reports_are_valid(self):
        for system, options in (
            (System(solute_left_mol=0, solute_right_mol=0.01), {}),
            (System(solute_left_mol=0.01, solute_right_mol=0), {}),
            (System(solute_left_mol=0, solute_right_mol=0), {}),
            (System(solute_left_mol=0, solute_right_mol=0, gravity_m_s2=0), {}),
            (System(), {"max_iterations": 1}),
        ):
            with self.subTest(system=system, options=options):
                result = solve(system, **options)
                root = ET.fromstring(render_svg(system, result))
                text = " ".join(root.itertext())
                self.assertIn(result.status.replace("_", " "), text)
                self.assertIn("not a time simulation", text)
                for rect in root.findall("s:rect", NS):
                    self.assertGreaterEqual(float(rect.get("height", 0)), 0)
                self.assertFalse(root.findall(".//s:script", NS))
                self.assertFalse(root.findall(".//s:image", NS))

    def test_cli_creates_svg_and_csv_in_nested_directory(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "reports" / "example.svg"
            csv_path = path.with_suffix(".csv")
            with redirect_stdout(io.StringIO()):
                code = main(["solve", "--svg", str(path), "--output", str(csv_path)])
            self.assertEqual(code, 0)
            self.assertTrue(csv_path.is_file())
            self.assertEqual(ET.parse(path).getroot().get("role"), "img")

    def test_cli_rejects_colliding_report_paths_before_writing(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "same.svg"
            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
                main(["solve", "--svg", str(path), "--output", str(path)])
            self.assertEqual(error.exception.code, 2)
            self.assertFalse(path.exists())


if __name__ == "__main__":
    unittest.main()
