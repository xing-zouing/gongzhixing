#!/usr/bin/env python3
"""Step 6 prototype: device-level global D2-symmetric placement.

The old flow placed each matching group as one block and then stacked blocks
along a common axis.  This prototype instead:

* splits every two-device matching group into two individual device cells;
* places each pair on the global x-axis or y-axis, symmetric about the center;
* uses the other axis for other pairs, forming a D2-symmetric cross;
* fills the remaining central gap with unmatched devices M9/M12.

The matched devices' centroids are D2-symmetric; unmatched devices are free
to fill gaps.  This is a placement prototype: each netlist device is one cell
(multiplicity is not expanded into separate unit instances yet).

Run from the project root:

    .venv/bin/python scripts/ota3_step6_global_d2.py
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
import sys

import gdspy

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cc.device_generation.Capacitor import Capacitor  # noqa: E402
from cc.device_generation.Mosfet import Mosfet  # noqa: E402
from cc.gds.render import render_gds  # noqa: E402
from cc.io.netlist_parser import Device, load_case  # noqa: E402
from cc.plan.ota3_groups import GROUPS  # noqa: E402


def _normalized(cell: gdspy.Cell, name: str) -> tuple[gdspy.Cell, float, float]:
    cell.flatten()
    bbox = cell.get_bounding_box()
    if bbox is None:
        raise ValueError(f"cell {name} has no geometry")
    ox, oy = bbox[0]
    width = float(bbox[1][0] - bbox[0][0])
    height = float(bbox[1][1] - bbox[0][1])
    normalized = gdspy.Cell(name, True)
    for (layer, datatype), polygons in cell.get_polygons(by_spec=True).items():
        for polygon in polygons:
            shifted = [(x - ox, y - oy) for x, y in polygon]
            normalized.add(gdspy.Polygon(shifted, layer=layer, datatype=datatype))
    return normalized, width, height


def _device_cell(dev: Device) -> tuple[gdspy.Cell, float, float]:
    if dev.kind == "mos":
        mos = Mosfet(
            nch=dev.nch,
            name=dev.name,
            w=dev.unit_width,
            l=dev.unit_length,
            nf=dev.nf,
            attr=dev.attr,
            spectre=False,
            pinConType=None,
            bulkCon=[0],
        )
        return _normalized(mos.cell, dev.name)
    if dev.kind == "cap":
        cap = Capacitor(
            name=dev.name,
            w=dev.unit_width,
            sp=float(dev.params.get("s", 0.07)),
            nf=dev.nr,
            l=dev.unit_length,
            m_bot=int(dev.params.get("stm", 3)),
            m_top=int(dev.params.get("spm", 5)),
            attr=["2t"] if "2t" in dev.model.lower() else [],
            f_tip=float(dev.params.get("ftip", 0.14)),
            flatten=True,
        )
        return _normalized(cap.cell, dev.name)
    raise NotImplementedError(f"unsupported device kind: {dev.kind}")


def main() -> None:
    out_dir = ROOT / "build" / "ota3" / "global_d2"
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in out_dir.glob("ota3_global_d2.*"):
        stale.unlink()

    devices, _ = load_case(ROOT / "examples" / "ota3.sp", ROOT / "examples" / "ota3.sym", subckt="OTA2_my")

    # Build one cell per device for each two-device matching group.
    pairs = []
    for spec in GROUPS:
        devs = [devices[name] for name in spec["devices"]]
        if len(devs) != 2:
            continue
        cells = []
        for dev in devs:
            cell, width, height = _device_cell(dev)
            cells.append({"name": dev.name, "cell": cell, "w": width, "h": height})
        pairs.append(
            {
                "group": spec["name"],
                "devices": cells,
                "max_w": max(c["w"] for c in cells),
                "max_h": max(c["h"] for c in cells),
            }
        )

    unmatched = []
    for name in ("M9", "M12"):
        dev = devices[name]
        cell, width, height = _device_cell(dev)
        unmatched.append({"name": name, "cell": cell, "w": width, "h": height})

    gap = 2.0
    n = len(pairs)
    best = None

    # Choose which pairs go on the horizontal x-axis vs the vertical y-axis.
    indices = range(n)
    for h_count in range(1, n):
        for h_set in combinations(indices, h_count):
            h_set = set(h_set)
            h_pairs = [pairs[i] for i in range(n) if i in h_set]
            v_pairs = [pairs[i] for i in range(n) if i not in h_set]

            h_half_height = max(p["max_h"] for p in h_pairs) / 2.0
            v_half_width = max(p["max_w"] for p in v_pairs) / 2.0
            central_gap = 2.0 * max(h_half_height, v_half_width) + gap

            h_half = central_gap / 2.0 + sum(p["max_w"] + gap for p in h_pairs)
            v_half = central_gap / 2.0 + sum(p["max_h"] + gap for p in v_pairs)
            width = 2.0 * h_half
            height = 2.0 * v_half
            aspect = width / height if height else 1.0
            score = width * height * (1.0 + 0.15 * abs(aspect - 1.0))
            if best is None or score < best[0]:
                best = (score, h_set, h_half, v_half, central_gap)

    _, h_set, h_half, v_half, central_gap = best
    h_pairs = [pairs[i] for i in range(n) if i in h_set]
    v_pairs = [pairs[i] for i in range(n) if i not in h_set]

    placements: list[dict] = []

    def place(device: dict, cx: float, cy: float) -> None:
        placements.append({**device, "cx": cx, "cy": cy})

    # Horizontal pairs: left/right, both on y=0.
    cursor = central_gap / 2.0
    for pair in h_pairs:
        distance = cursor + pair["max_w"] / 2.0
        left, right = pair["devices"]
        place(left, -distance, 0.0)
        place(right, +distance, 0.0)
        cursor = distance + pair["max_w"] / 2.0 + gap

    # Vertical pairs: bottom/top, both on x=0.
    cursor = central_gap / 2.0
    for pair in v_pairs:
        distance = cursor + pair["max_h"] / 2.0
        bottom, top = pair["devices"]
        place(bottom, 0.0, -distance)
        place(top, 0.0, +distance)
        cursor = distance + pair["max_h"] / 2.0 + gap

    # Fill the central gap with unmatched devices (row first, then column).
    if unmatched:
        total_w = sum(d["w"] for d in unmatched) + gap * (len(unmatched) - 1)
        total_h = sum(d["h"] for d in unmatched) + gap * (len(unmatched) - 1)
        if total_w <= central_gap:
            cursor_x = -total_w / 2.0
            for device in unmatched:
                place(device, cursor_x + device["w"] / 2.0, 0.0)
                cursor_x += device["w"] + gap
        elif total_h <= central_gap:
            cursor_y = -total_h / 2.0
            for device in unmatched:
                place(device, 0.0, cursor_y + device["h"] / 2.0)
                cursor_y += device["h"] + gap
        else:
            # Fallback: put them in the lower-right quadrant.
            cx = h_half + gap + max(d["w"] for d in unmatched) / 2.0
            cy = -v_half - gap - max(d["h"] for d in unmatched) / 2.0
            cursor_x = cx - sum(d["w"] for d in unmatched) / 2.0
            for device in unmatched:
                place(device, cursor_x + device["w"] / 2.0, cy)
                cursor_x += device["w"] + gap

    # Bounding box from real device extents, plus margin.
    margin = 5.0
    min_x = min(p["cx"] - p["w"] / 2.0 for p in placements)
    max_x = max(p["cx"] + p["w"] / 2.0 for p in placements)
    min_y = min(p["cy"] - p["h"] / 2.0 for p in placements)
    max_y = max(p["cy"] + p["h"] / 2.0 for p in placements)

    top_width = (max_x - min_x) + 2 * margin
    top_height = (max_y - min_y) + 2 * margin
    shift_x = -min_x + margin
    shift_y = -min_y + margin

    top_lib = gdspy.GdsLibrary()
    top = gdspy.Cell("ota3_global_d2")
    top_lib.add(top)
    top.add(gdspy.Rectangle((0, 0), (top_width, top_height), layer=99, datatype=0))

    for item in placements:
        x = item["cx"] + shift_x - item["w"] / 2.0
        y = item["cy"] + shift_y - item["h"] / 2.0
        top.add(gdspy.CellReference(item["cell"], origin=(x, y)))
        top.add(gdspy.Label(item["name"], (x + 0.4, y + 0.4), layer=100))
        top_lib.add(item["cell"])

    gds_path = out_dir / "ota3_global_d2.gds"
    top_lib.write_gds(str(gds_path))
    info = render_gds(gds_path, out_dir / "ota3_global_d2.png", dpi=180)

    print("=== global D2 placement ===")
    print("horizontal pairs:", [p["group"] for p in h_pairs])
    print("vertical pairs  :", [p["group"] for p in v_pairs])
    print(f"central gap     : {central_gap:.2f}")
    print(f"top size        : {top_width:.2f} x {top_height:.2f} um "
          f"(area {top_width * top_height:.1f} um^2)")
    print(f"GDS             : {gds_path}")
    print(f"PNG             : {info['png']}")

    # Report pair centroid symmetry.
    for pair in pairs:
        left, right = pair["devices"]
        left_pos = next(p for p in placements if p["name"] == left["name"])
        right_pos = next(p for p in placements if p["name"] == right["name"])
        sx = left_pos["cx"] + right_pos["cx"]
        sy = left_pos["cy"] + right_pos["cy"]
        print(f"  {pair['group']:28s} centroid sum=({sx:+.3f},{sy:+.3f})")


if __name__ == "__main__":
    main()
