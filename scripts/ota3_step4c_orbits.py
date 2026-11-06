#!/usr/bin/env python3
"""Step 4c: general symmetry-orbit allocation for multiple groups.

Demonstrates:

* three differential pairs sharing one 1-D array;
* a 1:2 ratioed group (weighted/device-level common centroid);
* a mixed array: two pairs + one 1:2 group;
* 2-D 180-degree interleaving, including multi-unit devices.

The devices in this script are synthetic but use the same unit-cell geometry
as the real ota3 pair; this is an algorithm demo, not a netlist-accurate
layout.  The real ota3 top still uses the k=1 ABBA block.

Run from the project root:

    .venv/bin/python scripts/ota3_step4c_orbits.py
"""

from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cc.place.interleave import (  # noqa: E402
    Group,
    allocate_orbits_1d,
    allocate_orbits_2d,
    group_metrics,
    print_group_metrics,
)
from cc.gds.array_writer import generate_group_gds  # noqa: E402
from cc.gds.render import render_gds  # noqa: E402


def _synthetic_device(name: str) -> dict:
    # Same unit-cell geometry as ota3 M2/M5/M4/M3.
    return {
        "name": name,
        "nch": False,
        "model": "pch",
        "w": 2.5,
        "l": 5.0,
        "nf": 1,
        "multi": 1,
        "attr": [],
    }


def _write_gds(name: str, layout, out_dir: Path) -> None:
    devices = []
    seen = set()
    for row in (layout if layout and isinstance(layout[0], list) else [layout]):
        for device in row:
            if device is not None and device not in seen:
                seen.add(device)
                devices.append(device)
    name_to_id = {device: index + 1 for index, device in enumerate(devices)}

    def to_ids(row):
        return [0 if device is None else name_to_id[device] for device in row]

    matrix = (
        [to_ids(row) for row in layout]
        if layout and isinstance(layout[0], list)
        else [to_ids(layout)]
    )

    output = out_dir / f"{name}.gds"
    info = generate_group_gds(
        devices=[_synthetic_device(device) for device in devices],
        cc_matrix=matrix,
        output_path=output,
        array_name=name,
    )
    print(f"  GDS: {info['path']}")
    render_gds(output, out_dir / f"{name}.png", dpi=200)


def main() -> None:
    out_dir = ROOT / "build" / "ota3"
    out_dir.mkdir(parents=True, exist_ok=True)

    for filename in (
        "orbits_mixed_groups",
        "orbits_ratio_1to2",
        "orbits_2d_multiunit",
    ):
        for suffix in (".gds", ".png"):
            stale = out_dir / f"{filename}{suffix}"
            if stale.exists():
                stale.unlink()

    # ------------------------------------------------------------------
    # 1-D demo A: three pairs, one unit per device, 6 cells.
    # ------------------------------------------------------------------
    groups_a = [
        Group("A", ("A0", "A1"), (1, 1)),
        Group("B", ("B0", "B1"), (1, 1)),
        Group("C", ("C0", "C1"), (1, 1)),
    ]
    layout_a = allocate_orbits_1d(groups_a, 6)
    print("\n1-D three pairs layout:", layout_a)
    print_group_metrics("1-D three pairs", group_metrics([layout_a], groups_a))

    # ------------------------------------------------------------------
    # 1-D demo B: one 1:2 ratio group, 3 cells.
    # ------------------------------------------------------------------
    groups_b = [Group("R", ("R_unit", "R_twice"), (1, 2))]
    layout_b = allocate_orbits_1d(groups_b, 3)
    print("\n1-D 1:2 ratio layout:", layout_b)
    print_group_metrics("1-D 1:2 ratio", group_metrics([layout_b], groups_b))

    # ------------------------------------------------------------------
    # 1-D demo C: two pairs + one 1:2 ratio, 7 cells.
    # ------------------------------------------------------------------
    groups_c = [
        Group("A", ("A0", "A1"), (1, 1)),
        Group("B", ("B0", "B1"), (1, 1)),
        Group("R", ("R_unit", "R_twice"), (1, 2)),
    ]
    layout_c = allocate_orbits_1d(groups_c, 7)
    print("\n1-D mixed groups layout:", layout_c)
    print_group_metrics("1-D mixed groups", group_metrics([layout_c], groups_c))

    # ------------------------------------------------------------------
    # 2-D demo D: two pairs, 2x2 checkerboard.
    # ------------------------------------------------------------------
    groups_d = [
        Group("A", ("A0", "A1"), (1, 1)),
        Group("B", ("B0", "B1"), (1, 1)),
    ]
    layout_d = allocate_orbits_2d(groups_d, 2, 2)
    print("\n2-D two pairs layout:")
    for row in layout_d:
        print("  ", row)
    print_group_metrics("2-D two pairs", group_metrics(layout_d, groups_d))

    # ------------------------------------------------------------------
    # 2-D demo E: two pairs, two units per device, 2x4.
    # ------------------------------------------------------------------
    groups_e = [
        Group("A", ("A0", "A1"), (2, 2)),
        Group("B", ("B0", "B1"), (2, 2)),
    ]
    layout_e = allocate_orbits_2d(groups_e, 2, 4)
    print("\n2-D multi-unit layout:")
    for row in layout_e:
        print("  ", row)
    print_group_metrics("2-D multi-unit", group_metrics(layout_e, groups_e))

    # ------------------------------------------------------------------
    # GDS/PNG for representative cases.
    # ------------------------------------------------------------------
    print("\n=== GDS outputs ===")
    _write_gds("orbits_ratio_1to2", layout_b, out_dir)
    _write_gds("orbits_mixed_groups", layout_c, out_dir)
    _write_gds("orbits_2d_multiunit", layout_e, out_dir)


if __name__ == "__main__":
    main()
