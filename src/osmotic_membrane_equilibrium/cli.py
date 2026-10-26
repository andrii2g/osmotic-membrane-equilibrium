"""Command-line experiments; input convenience units, SI-unit result records."""

import argparse
import json
from pathlib import Path

from .experiments import (
    linear_values, result_record, sweep_concentration, sweep_temperature, write_csv,
)
from .model import System
from .solver import solve
from .svg import write_svg


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Ideal dilute osmosis in equal-area compartments.")
    commands = parser.add_subparsers(dest="command", required=True)
    for command, help_text in (
        ("solve", "Find the equilibrium water distribution."),
        ("sweep-temperature", "Vary temperature at fixed solute amounts and geometry."),
        ("sweep-concentration", "Vary initial right-minus-left concentration at fixed mean."),
    ):
        child = commands.add_parser(command, help=help_text)
        child.add_argument("--left-volume-ml", type=float, default=1000, help="Initial left liquid volume (default: 1000 mL).")
        child.add_argument("--right-volume-ml", type=float, default=1000, help="Initial right liquid volume (default: 1000 mL).")
        child.add_argument("--left-solute-mmol", type=float, default=0.1, help="Left solute particles (default: 0.1 mmol).")
        child.add_argument("--right-solute-mmol", type=float, default=0.2, help="Right solute particles (default: 0.2 mmol).")
        child.add_argument("--temperature-k", type=float, default=298.15)
        child.add_argument("--area-cm2", type=float, default=100, help="Area of EACH compartment (default: 100 cm^2).")
        child.add_argument("--density-kg-m3", type=float, default=1000)
        child.add_argument("--gravity-m-s2", type=float, default=9.80665)
        child.add_argument("--pressure-tolerance-pa", type=float, default=1e-6)
        child.add_argument("--volume-tolerance-m3", type=float, default=1e-15)
        child.add_argument("--relative-volume-tolerance", type=float, default=1e-12)
        child.add_argument("--max-iterations", type=int, default=200)
        child.add_argument("--json", action="store_true", help="Print SI-unit result records as JSON.")
        child.add_argument("--output", type=Path, help="Write all SI-unit metrics to CSV.")
        if command == "solve":
            child.add_argument("--svg", type=Path, help="Write a standalone compartment and pressure-balance report.")
        if command != "solve":
            child.add_argument("--points", type=int, default=13)
        if command == "sweep-temperature":
            child.add_argument("--start-k", type=float, default=278.15)
            child.add_argument("--stop-k", type=float, default=338.15)
        elif command == "sweep-concentration":
            child.add_argument("--start-mol-m3", type=float, default=0)
            child.add_argument("--stop-mol-m3", type=float, default=0.3)
    return parser


def _print_solution(record: dict) -> None:
    print(f"Status: {record['status']}")
    print("                         Initial          Final")
    for label, field, scale, unit in (
        ("Left height", "height_left_m", 100, "cm"),
        ("Right height", "height_right_m", 100, "cm"),
        ("Left concentration", "concentration_left_mol_m3", 1, "mol/m^3"),
        ("Right concentration", "concentration_right_mol_m3", 1, "mol/m^3"),
        ("Osmotic difference", "osmotic_difference_pa", 1, "Pa"),
        ("Hydrostatic difference", "hydrostatic_difference_pa", 1, "Pa"),
    ):
        print(f"{label:<23} {record['initial_' + field] * scale:>12.6g} "
              f"{record['final_' + field] * scale:>14.6g} {unit}")
    print(f"Transfer to left: {record['transferred_to_left_m3'] * 1e6:.6g} mL")
    print(f"Iterations: {record['iterations']}; residual: {record['final_residual_pa']:.6g} Pa")
    print(f"Volume error: {record['volume_error_m3']:.6g} m^3; "
          f"bracket width: {record['bracket_width_m3']:.6g} m^3")
    if record["status"] == "dry_boundary":
        print("The solute-free compartment reaches the dry boundary; no wet equilibrium exists.")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "solve" and args.svg and args.output:
            if args.svg.resolve() == args.output.resolve():
                raise ValueError("SVG and CSV output paths must be different")
        system = System(
            volume_left_m3=args.left_volume_ml * 1e-6,
            volume_right_m3=args.right_volume_ml * 1e-6,
            solute_left_mol=args.left_solute_mmol * 1e-3,
            solute_right_mol=args.right_solute_mmol * 1e-3,
            temperature_k=args.temperature_k,
            area_m2=args.area_cm2 * 1e-4,
            density_kg_m3=args.density_kg_m3,
            gravity_m_s2=args.gravity_m_s2,
        )
        options = dict(
            pressure_tolerance_pa=args.pressure_tolerance_pa,
            volume_tolerance_m3=args.volume_tolerance_m3,
            relative_volume_tolerance=args.relative_volume_tolerance,
            max_iterations=args.max_iterations,
        )
        if args.command == "solve":
            result = solve(system, **options)
            records = [result_record(system, result)]
            if args.svg:
                write_svg(args.svg, system, result)
        elif args.command == "sweep-temperature":
            records = sweep_temperature(system, linear_values(args.start_k, args.stop_k, args.points), **options)
        else:
            records = sweep_concentration(system, linear_values(args.start_mol_m3, args.stop_mol_m3, args.points), **options)
        if args.output:
            write_csv(args.output, records)
        if args.json:
            print(json.dumps(records[0] if args.command == "solve" else records, indent=2, allow_nan=False))
        elif args.command == "solve":
            _print_solution(records[0])
        else:
            print(f"{len(records)} experiments; "
                  f"{sum(row['converged'] for row in records)} converged; "
                  f"maximum absolute residual {max(abs(row['final_residual_pa']) for row in records):.6g} Pa")
        if args.output and not args.json:
            print(f"CSV: {args.output}")
        if args.command == "solve" and args.svg and not args.json:
            print(f"SVG: {args.svg}")
        if any(row["status"] in {"max_iterations", "precision_limit"} for row in records):
            return 1
        if any(row["status"] == "dry_boundary" for row in records):
            return 3
        return 0
    except (ValueError, OSError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
