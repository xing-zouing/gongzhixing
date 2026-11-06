#!/usr/bin/env python3
"""Step 7 prototype: global vertical-axis, device-level pair placement.

Model:

* one global vertical symmetry axis at x = 0;
* each matching pair is placed left/right of the axis at the same y;
* pairs are stacked in rows; every pair centroid lies on the global axis;
* x-gradient cancels by left/right mirroring;
* y-gradient cancels because both devices of a pair share the same y;
* unmatched devices fill side/row gaps (no symmetry constraint).

One row per pair for this first version.  Devices are one cell per netlist
device; ``multi`` is not expanded into unit instances yet.

Run from the project root:

    .venv/bin/python scripts/ota3_step7_global_axis.py
"""

from __future__ import annotations

from pathlib import Path
import sys

import gdspy

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cc.gds.device_cells import device_cell  # noqa: E402
from cc.gds.render import render_gds  # noqa: E402
from cc.io.netlist_parser import load_case  # noqa: E402
from cc.plan.ota3_groups import GROUPS  # noqa: E402


def main() -> None:
    out_dir = ROOT / "build" / "ota3" / "global_axis"
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in out_dir.glob("ota3_global_axis.*"):
        stale.unlink()

    devices, _ = load_case(
        ROOT / "examples" / "ota3.sp",
        ROOT / "examples" / "ota3.sym",
        subckt="OTA2_my",
    )

    pairs = []
    for spec in GROUPS:
        devs = [devices[name] for name in spec["devices"]]
        if len(devs) != 2:
            continue
        cells = []
        for dev in devs:
            cell, width, height = device_cell(dev)
            cells.append({"name": dev.name, "cell": cell, "w": width, "h": height})
        pairs.append({"group": spec["name"], "devices": cells})

    unmatched = []
    for name in ("M9", "M12"):
        cell, width, height = device_cell(devices[name])
        unmatched.append({"name": name, "cell": cell, "w": width, "h": height})

    gap = 2.0
    margin = 5.0

    # One pair per row, centered on x = 0.  Sort by height (tallest first)
    # to keep the outer rows easy to fill.
    pairs.sort(key=lambda p: max(d["h"] for d in p["devices"]), reverse=True)

    rows = []
    cursor_y = 0.0
    for pair in pairs:
        left, right = pair["devices"]
        half = max(left["w"], right["w"]) / 2.0 + gap / 2.0
        row_height = max(left["h"], right["h"])
        row_width = 2.0 * (half + max(left["w"], right["w"]) / 2.0)
        center_y = cursor_y + row_height / 2.0
        rows.append(
            {
                "pair": pair,
                "center_y": center_y,
                "half": half,
                "height": row_height,
                "width": row_width,
                "left": left,
                "right": right,
            }
        )
        cursor_y += row_height + gap

    top_width = max(row["width"] for row in rows)
    top_height = cursor_y - gap

    placements = []

    def place(device: dict, cx: float, cy: float) -> None:
        placements.append({**device, "cx": cx, "cy": cy})

    for row in rows:
        place(row["left"], -row["half"], row["center_y"])
        place(row["right"], +row["half"], row["center_y"])

    # Fill side gaps of the narrower rows with the unmatched devices.
    for device in unmatched:
        placed = False
        for row in sorted(rows, key=lambda r: r["width"]):
            side_space = (top_width - row["width"]) / 2.0
            if side_space >= device["w"] + gap:
                # Put one on the far left of the row; unmatched devices do
                # not need to be symmetric.
                cx = -top_width / 2.0 + device["w"] / 2.0 + gap
                place(device, cx, row["center_y"])
                row["width"] = max(row["width"], 2.0 * (abs(cx) + device["w"] / 2.0))
                placed = True
                break
        if not placed:
            # Fallback: append a dedicated row at the bottom.
            cy = cursor_y + device["h"] / 2.0
            place(device, 0.0, cy)
            cursor_y += device["h"] + gap
            top_height = cursor_y - gap

    # Recompute extents from real device boxes.
    min_x = min(p["cx"] - p["w"] / 2.0 for p in placements)
    max_x = max(p["cx"] + p["w"] / 2.0 for p in placements)
    min_y = min(p["cy"] - p["h"] / 2.0 for p in placements)
    max_y = max(p["cy"] + p["h"] / 2.0 for p in placements)

    top_width = (max_x - min_x) + 2 * margin
    top_height = (max_y - min_y) + 2 * margin
    shift_x = -min_x + margin
    shift_y = -min_y + margin

    top_lib = gdspy.GdsLibrary()
    top = gdspy.Cell("ota3_global_axis")
    top_lib.add(top)
    top.add(gdspy.Rectangle((0, 0), (top_width, top_height), layer=99, datatype=0))

    for item in placements:
        x = item["cx"] + shift_x - item["w"] / 2.0
        y = item["cy"] + shift_y - item["h"] / 2.0
        top.add(gdspy.CellReference(item["cell"], origin=(x, y)))
        top.add(gdspy.Label(item["name"], (x + 0.4, y + 0.4), layer=100))
        top_lib.add(item["cell"])

    gds_path = out_dir / "ota3_global_axis.gds"
    top_lib.write_gds(str(gds_path))
    info = render_gds(gds_path, out_dir / "ota3_global_axis.png", dpi=180)

    print("=== global vertical-axis placement ===")
    print(f"top size: {top_width:.2f} x {top_height:.2f} um "
          f"(area {top_width * top_height:.1f} um^2)")
    print(f"GDS: {gds_path}")
    print(f"PNG: {info['png']}")

    for row in rows:
        left = next(p for p in placements if p["name"] == row["left"]["name"])
        right = next(p for p in placements if p["name"] == row["right"]["name"])
        print(
            f"  {row['pair']['group']:28s} "
            f"x-sum={left['cx'] + right['cx']:+.3f} "
            f"y-diff={left['cy'] - right['cy']:+.3f} "
            f"y={left['cy']:.2f}"
        )


if __name__ == "__main__":
    main()
