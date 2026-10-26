"""Standalone SVG reports with a shared height scale and pressure-balance plot."""

import math
from pathlib import Path
from xml.etree import ElementTree as ET

from .model import State, System, evaluate
from .solver import Result

SVG = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG)
INK = "#183344"
MUTED = "#536c7a"
TEAL = "#087f8c"
PURPLE = "#8054bb"


def _node(parent, tag, text=None, **attributes):
    node = ET.SubElement(parent, f"{{{SVG}}}{tag}", {
        key.replace("_", "-"): str(value) for key, value in attributes.items()
    })
    if text is not None:
        node.text = str(text)
    return node


def _text(parent, x, y, text, size=14, color=INK, **attributes):
    return _node(parent, "text", text, x=x, y=y, font_size=size, fill=color, **attributes)


def _panel(parent, x, y, width, height):
    return _node(parent, "rect", x=x, y=y, width=width, height=height,
                 rx=18, fill="#ffffff", stroke="#dce7eb")


def _compartments(parent, x, title, state: State, system: System, scale: float):
    _panel(parent, x, 108, 554, 392)
    _text(parent, x + 24, 141, title, 20, font_weight="bold")
    bottom, top, width = 370, 174, 148
    # The submerged connector is a schematic, with negligible storage volume.
    _node(parent, "rect", x=x + 200, y=354, width=134, height=16, fill="#d0eef4")
    _node(parent, "path", d=f"M{x + 200} 354 H{x + 334} M{x + 200} 370 H{x + 334}",
          fill="none", stroke=INK, stroke_width=2)
    maximum_amount = max(system.solute_left_mol, system.solute_right_mol)
    for side, vessel_x, height, volume, amount, concentration, osmotic, color in (
        ("left", x + 52, state.height_left_m, state.volume_left_m3,
         system.solute_left_mol, state.concentration_left_mol_m3, state.osmotic_left_pa, "#afe0ec"),
        ("right", x + 334, state.height_right_m, state.volume_right_m3,
         system.solute_right_mol, state.concentration_right_mol_m3, state.osmotic_right_pa, "#d3c6ee"),
    ):
        pixels = height * scale
        surface = bottom - pixels
        _text(parent, vessel_x + width / 2, 163, side.upper(), 11, MUTED,
              text_anchor="middle", letter_spacing=2)
        _node(parent, "rect", x=vessel_x, y=f"{surface:.6f}", width=width,
              height=f"{pixels:.6f}", fill=color,
              **{"class": f"water-{side}", "data-height-m": repr(height)})
        count = round(20 * amount / maximum_amount) if maximum_amount else 0
        for i in range(count):
            dot_x = vessel_x + 18 + (i % 5) * 28
            dot_y = bottom - (12 + i // 5 * 22) / 95 * max(0, pixels - 4)
            _node(parent, "circle", cx=dot_x, cy=f"{dot_y:.4f}", r=2.6,
                  fill=TEAL if side == "left" else PURPLE, opacity=0.8)
        _node(parent, "path", d=f"M{vessel_x} {top} V{bottom} H{vessel_x + width} V{top}",
              fill="none", stroke=INK, stroke_width=2)
        if pixels:
            _node(parent, "line", x1=vessel_x, x2=vessel_x + width,
                  y1=surface, y2=surface, stroke=TEAL if side == "left" else PURPLE,
                  stroke_width=2)
        _text(parent, vessel_x + width / 2, surface - 9, f"{height * 100:.4g} cm", 15,
              text_anchor="middle", font_weight="bold")
        _text(parent, vessel_x + width / 2, 396,
              f"{volume * 1e6:.5g} mL · {amount * 1000:.4g} mmol", 13, text_anchor="middle")
        _text(parent, vessel_x + width / 2, 420,
              f"c = {concentration:.5g} mol/m³", 13, MUTED, text_anchor="middle")
        _text(parent, vessel_x + width / 2, 444,
              f"π = {osmotic:.5g} Pa", 13, MUTED, text_anchor="middle")
    _node(parent, "line", x1=x + 267, x2=x + 267, y1=344, y2=379,
          stroke=TEAL, stroke_width=5, stroke_dasharray="4 3", **{"class": "membrane"})
    _text(parent, x + 267, 396, "membrane", 10, MUTED, text_anchor="middle")
    _text(parent, x + 267, 413, "water only", 10, MUTED, text_anchor="middle")
    _text(parent, x + 24, 478,
          f"ΔP = {state.hydrostatic_difference_pa:.5g} Pa     Δπ = {state.osmotic_difference_pa:.5g} Pa",
          14)


def _pressure_plot(parent, system: System, result: Result):
    _panel(parent, 32, 526, 740, 394)
    _text(parent, 56, 560, "Where the pressures balance", 20, font_weight="bold")
    _text(parent, 56, 584, "Both differences are left minus right. The curves meet at wet equilibrium.", 13, MUTED)
    total = system.total_volume_m3
    start, final = result.initial.volume_left_m3, result.final.volume_left_m3
    margin = max(total * 0.08, abs(start - final) * 0.25)
    low, high = min(start, final) - margin, max(start, final) + margin
    if low <= 0:
        low = min(start, final) / 2 if system.solute_left_mol else 0
    if high >= total:
        high = total - (total - max(start, final)) / 2 if system.solute_right_mol else total
    if system.solute_left_mol:
        low = max(math.nextafter(0, total), low)
    if system.solute_right_mol:
        high = min(math.nextafter(total, 0), high)
    samples = [evaluate(system, low + (high - low) * i / 160) for i in range(161)]
    values = [0.0] + [value for state in samples for value in
                      (state.hydrostatic_difference_pa, state.osmotic_difference_pa)]
    y_low, y_high = min(values), max(values)
    padding = max((y_high - y_low) * 0.12, 1.0)
    y_low, y_high = y_low - padding, y_high + padding
    left, right, top, bottom = 106, 740, 615, 830

    def px(volume):
        return left + (volume - low) / (high - low) * (right - left)

    def py(pressure):
        return bottom - (pressure - y_low) / (y_high - y_low) * (bottom - top)

    for i in range(5):
        fraction = i / 4
        pressure = y_low + fraction * (y_high - y_low)
        y = py(pressure)
        _node(parent, "line", x1=left, x2=right, y1=y, y2=y, stroke="#e4ecef")
        _text(parent, left - 10, y + 4, f"{pressure:.4g}", 11, MUTED, text_anchor="end")
        volume = low + fraction * (high - low)
        x = px(volume)
        _node(parent, "line", x1=x, x2=x, y1=top, y2=bottom, stroke="#edf2f4")
        _text(parent, x, bottom + 21, f"{(volume - start) * 1e6:.4g}", 11, MUTED, text_anchor="middle")
    _text(parent, 56, 607, "Pa", 12, MUTED)
    _node(parent, "line", x1=left, x2=right, y1=py(0), y2=py(0), stroke="#a9bac3")
    for name, color in (("hydrostatic_difference_pa", TEAL), ("osmotic_difference_pa", PURPLE)):
        points = " ".join(f"{px(s.volume_left_m3):.3f},{py(getattr(s, name)):.3f}" for s in samples)
        _node(parent, "polyline", points=points, fill="none", stroke=color,
              stroke_width=2.8, **{"class": name})
    for volume, label, color in ((start, "initial", MUTED),
                                 (final, "equilibrium" if result.converged else "final", INK)):
        x = px(volume)
        _node(parent, "line", x1=x, x2=x, y1=top, y2=bottom,
              stroke=color, stroke_dasharray="4 4", opacity=0.6)
        _text(parent, x + 5, top + 14 if label == "initial" else top + 31,
              label, 11, color)
    for pressure, color in ((result.final.hydrostatic_difference_pa, TEAL),
                            (result.final.osmotic_difference_pa, PURPLE)):
        _node(parent, "circle", cx=px(final), cy=py(pressure), r=4.5,
              fill=color, stroke="white", stroke_width=1.5)
    _text(parent, (left + right) / 2, 878, "Water transferred to left (mL); negative means transfer to right", 13,
          MUTED, text_anchor="middle")
    for x, color, label in ((106, TEAL, "Hydrostatic ΔP"), (350, PURPLE, "Osmotic Δπ")):
        _node(parent, "line", x1=x, x2=x + 28, y1=902, y2=902, stroke=color, stroke_width=3)
        _text(parent, x + 37, 906, label, 13)


def render_svg(system: System, result: Result) -> str:
    root = ET.Element(f"{{{SVG}}}svg", {
        "width": "1200", "height": "1000", "viewBox": "0 0 1200 1000",
        "role": "img", "aria-labelledby": "report-title report-description",
        "font-family": "Segoe UI, DejaVu Sans, sans-serif",
    })
    _node(root, "title", "Osmotic membrane equilibrium report", id="report-title")
    _node(root, "desc", f"Initial and final water heights on a common scale. "
          f"Status: {result.status}. Final pressure residual {result.final.residual_pa:.6g} Pa. "
          "The pressure plot shows candidate distributions, not elapsed time.", id="report-description")
    _node(root, "rect", width=1200, height=1000, fill="#f1f6f8")
    _text(root, 32, 49, "Water moves. Solute stays.", 30, font_weight="bold")
    _text(root, 32, 79,
          f"Ideal dilute osmosis · {system.temperature_k:g} K · {system.area_m2 * 1e4:g} cm² per compartment · common height scale",
          15, MUTED)
    _node(root, "rect", x=974, y=24, width=194, height=37, rx=18,
          fill="#d9efe8" if result.converged else "#fbe6c8")
    _text(root, 1071, 48, result.status.replace("_", " "), 14,
          text_anchor="middle", font_weight="bold")
    maximum = max(result.initial.height_left_m, result.initial.height_right_m,
                  result.final.height_left_m, result.final.height_right_m)
    scale = 185 / maximum
    _compartments(root, 32, "01 / Initial distribution", result.initial, system, scale)
    final_title = "02 / Equilibrium" if result.converged else "02 / Final estimate"
    if result.status == "dry_boundary":
        final_title = "02 / Dry boundary — wet model stops"
    _compartments(root, 614, final_title, result.final, system, scale)
    _pressure_plot(root, system, result)
    _panel(root, 796, 526, 372, 394)
    _text(root, 820, 560, "Result", 20, font_weight="bold")
    transfer = result.transferred_to_left_m3 * 1e6
    direction = "to left" if transfer > 0 else "to right" if transfer < 0 else "no net transfer"
    for y, label, value in (
        (601, "Redistribution", f"{abs(transfer):.5g} mL · {direction}"),
        (660, "Final heights · left / right", f"{result.final.height_left_m * 100:.5g} / {result.final.height_right_m * 100:.5g} cm"),
        (719, "Pressure residual · ΔP − Δπ", f"{result.final.residual_pa:.4g} Pa"),
        (778, "Bisection iterations", str(result.iterations)),
        (837, "Total-volume error", f"{result.volume_error_m3:.4g} m³"),
    ):
        _text(root, 820, y, label, 13, MUTED)
        _text(root, 820, y + 26, value, 20, font_weight="bold")
    _text(root, 32, 952,
          "Idealized: dilute, well mixed, incompressible; impermeant solute; rigid compartments; common atmosphere and density.",
          13, MUTED)
    _text(root, 32, 977,
          "Dots symbolize fixed solute amounts. The connector has negligible volume. The pressure plot is not a time simulation.",
          13, MUTED)
    return ET.tostring(root, encoding="unicode", xml_declaration=True)


def write_svg(path: str | Path, system: System, result: Result) -> None:
    content = render_svg(system, result)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
