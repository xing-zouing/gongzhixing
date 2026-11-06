#!/usr/bin/env python3
"""Step 10: unit-level placement for multi-unit matched devices.

M0/M7 have ``multi=4``.  Instead of folding the four units into one rigid
macro, this script places each unit independently:

* device M0 centroid at ``(-d, y)``; its four units at ``(-d ± dx, y ± dy)``;
* device M7 centroid at ``(+d, y)``; its four units at ``(+d ∓ dx, y ± dy)``
  with mirrored unit cells;
* other devices can occupy the space between the upper and lower unit rows;
* every other pair/single keeps the macro + 0/90-degree rotation model.

The SA objective remains area + overlap + HPWL + compactness + aspect.

Run from the project root:

    .venv/bin/python scripts/ota3_step10_unit_level.py
"""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path
import sys

import gdspy
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle as MplRectangle  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cc.gds.cell_ops import rotate_cell_90  # noqa: E402
from cc.gds.device_cells import (  # noqa: E402
    device_cell,
    device_macro,
    mirror_cell_horizontally,
)
from cc.gds.render import render_gds  # noqa: E402
from cc.io.netlist_parser import load_case  # noqa: E402
from cc.plan.ota3_groups import GROUPS  # noqa: E402

from scripts.ota3_step9_sa_refine import (  # noqa: E402
    Rect,
    _device_entry,
    _dims,
    _overlap_area,
)


def _unit_rects(pair, state):
    y, rotation, left_offsets, right_offsets = state
    rects = []
    device_of_rect = {}
    left_name = pair["left"]["name"]
    right_name = pair["right"]["name"]
    unit_w, unit_h = pair["unit_w"], pair["unit_h"]
    if rotation:
        unit_w, unit_h = unit_h, unit_w

    for name, offsets in ((left_name, left_offsets), (right_name, right_offsets)):
        for index, (ox, oy) in enumerate(offsets):
            cx = ox
            cy = y + oy
            rect_name = f"{name}#{index}"
            rects.append(
                Rect(rect_name, cx - unit_w / 2, cy - unit_h / 2,
                     cx + unit_w / 2, cy + unit_h / 2)
            )
            device_of_rect[rect_name] = name
    return rects, device_of_rect


def _checkerboard_offsets(
    unit_w: float,
    unit_h: float,
    count: int,
    gap: float,
) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
    """Initial interleaved checkerboard for the two devices.

    A 4x4 grid (for 8 units each) is split by checkerboard parity; each
    device gets 8 cells whose centroid is the grid center.  The offsets of
    each device are then shifted to sum to zero, preserving that centroid.
    """
    total = 2 * count
    columns = math.ceil(math.sqrt(total))
    rows = math.ceil(total / columns)
    step_x = unit_w + gap
    step_y = unit_h + gap

    left: list[tuple[float, float]] = []
    right: list[tuple[float, float]] = []
    for index in range(total):
        row, col = divmod(index, columns)
        ox = (col - (columns - 1) / 2.0) * step_x
        oy = (row - (rows - 1) / 2.0) * step_y
        if (row + col) % 2 == 0:
            left.append((ox, oy))
        else:
            right.append((ox, oy))

    def zero_mean(points):
        if not points:
            return []
        mx = sum(point[0] for point in points) / len(points)
        my = sum(point[1] for point in points) / len(points)
        return [(x - mx, y - my) for x, y in points]

    return zero_mean(left), zero_mean(right)


def _band_offsets(
    unit_w: float,
    unit_h: float,
    count: int,
    gap: float,
    central_gap: float = 12.0,
) -> list[tuple[float, float]]:
    """Two columns split into upper/lower groups with a central band."""
    columns = 2
    rows = math.ceil(count / columns)
    half_rows = rows // 2
    step_x = unit_w + gap
    step_y = unit_h + gap
    offsets = []
    for index in range(count):
        row, col = divmod(index, columns)
        ox = (col - (columns - 1) / 2.0) * step_x
        if row < half_rows:
            oy = central_gap / 2.0 + (half_rows - 1 - row + 0.5) * step_y
        else:
            oy = -central_gap / 2.0 - (row - half_rows + 0.5) * step_y
        offsets.append((ox, oy))
    mx = sum(point[0] for point in offsets) / len(offsets)
    my = sum(point[1] for point in offsets) / len(offsets)
    return [(x - mx, y - my) for x, y in offsets]


def _build_rects(unit_state, macro_state, single_state, unit_pairs, macro_pairs, singles):
    rects = []
    device_of_rect = {}
    for pair in unit_pairs:
        pair_rects, mapping = _unit_rects(pair, unit_state[pair["name"]])
        rects.extend(pair_rects)
        device_of_rect.update(mapping)
    for pair in macro_pairs:
        d, y, rotation = macro_state[pair["name"]]
        left, right = pair["left"], pair["right"]
        lw, lh = _dims(left, rotation)
        rw, rh = _dims(right, rotation)
        rects.append(Rect(left["name"], -d - lw / 2, y - lh / 2,
                          -d + lw / 2, y + lh / 2))
        rects.append(Rect(right["name"], d - rw / 2, y - rh / 2,
                          d + rw / 2, y + rh / 2))
        device_of_rect[left["name"]] = left["name"]
        device_of_rect[right["name"]] = right["name"]
    for single in singles:
        x, y, rotation = single_state[single["name"]]
        w, h = _dims(single, rotation)
        rects.append(Rect(single["name"], x - w / 2, y - h / 2,
                          x + w / 2, y + h / 2))
        device_of_rect[single["name"]] = single["name"]
    return rects, device_of_rect


def _device_centers(rects, device_of_rect):
    groups: dict[str, list[tuple[float, float]]] = {}
    for rect in rects:
        device = device_of_rect[rect.name]
        groups.setdefault(device, []).append(
            ((rect.x0 + rect.x1) / 2.0, (rect.y0 + rect.y1) / 2.0)
        )
    return {
        device: (
            sum(point[0] for point in points) / len(points),
            sum(point[1] for point in points) / len(points),
        )
        for device, points in groups.items()
    }


def _evaluate(
    unit_state,
    macro_state,
    single_state,
    unit_pairs,
    macro_pairs,
    singles,
    gap,
    hpwl_info,
    lambda_hpwl,
    lambda_compact,
    max_width,
):
    rects, device_of_rect = _build_rects(
        unit_state, macro_state, single_state, unit_pairs, macro_pairs, singles
    )
    centers = _device_centers(rects, device_of_rect)

    min_x = min(r.x0 for r in rects)
    max_x = max(r.x1 for r in rects)
    min_y = min(r.y0 for r in rects)
    max_y = max(r.y1 for r in rects)
    width = max_x - min_x
    height = max_y - min_y

    overlap = 0.0
    for i in range(len(rects)):
        for j in range(i + 1, len(rects)):
            overlap += _overlap_area(rects[i], rects[j], gap)

    hpwl = 0.0
    for name, nets in hpwl_info.items():
        center = centers.get(name)
        if center is None:
            continue
        for others in nets.values():
            points = [center]
            for other in others:
                if other in centers:
                    points.append(centers[other])
            if len(points) > 1:
                xs = [point[0] for point in points]
                ys = [point[1] for point in points]
                hpwl += max(xs) - min(xs) + max(ys) - min(ys)

    cx0 = (min_x + max_x) / 2.0
    cy0 = (min_y + max_y) / 2.0
    compactness = sum(
        (center[0] - cx0) ** 2 + (center[1] - cy0) ** 2
        for center in centers.values()
    )

    aspect = width / height if height else 1.0
    aspect_penalty = 1.0 + 0.8 * abs(math.log(max(aspect, 1e-9) / 1.0))
    width_penalty = 500.0 * max(0.0, width - max_width) ** 2
    score = (
        width * height * aspect_penalty
        + 1000.0 * overlap
        + lambda_hpwl * hpwl
        + lambda_compact * compactness
        + width_penalty
    )
    return {
        "score": score,
        "area": width * height,
        "overlap": overlap,
        "hpwl": hpwl,
        "compactness": compactness,
        "width": width,
        "height": height,
        "rects": rects,
        "device_of_rect": device_of_rect,
        "centers": centers,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--init", choices=("checkerboard", "random", "central"), default="checkerboard"
    )
    parser.add_argument("--iterations", type=int, default=40000)
    parser.add_argument(
        "--out-dir", default=str(ROOT / "build" / "ota3" / "unit_level")
    )
    parser.add_argument("--name", default="ota3_unit_level")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in out_dir.glob(f"{args.name}.*"):
        stale.unlink()
    rng = random.Random(args.seed)

    devices, _ = load_case(
        ROOT / "examples" / "ota3.sp",
        ROOT / "examples" / "ota3.sym",
        subckt="OTA2_my",
    )

    gap = 0.6
    unit_pairs = []
    macro_pairs = []

    for spec in GROUPS:
        devs = [devices[name] for name in spec["devices"]]
        if len(devs) != 2:
            continue
        left_dev, right_dev = devs

        if left_dev.multi > 1 or right_dev.multi > 1:
            # Unit-level pair: split each physical device into units of
            # half the original finger count (10 -> 5 fingers), so M0/M7
            # become 8 units each instead of 4 long 10-finger units.
            if left_dev.multi != right_dev.multi or left_dev.nf != right_dev.nf:
                raise NotImplementedError("unit-level pair must have equal multi/nf")
            fingers_per_unit = left_dev.nf // 2 if left_dev.nf % 2 == 0 else left_dev.nf
            unit_count = left_dev.multi * (left_dev.nf // fingers_per_unit)
            unit_cell, unit_w, unit_h = device_cell(
                left_dev, nf_override=fingers_per_unit
            )
            unit_rot90 = rotate_cell_90(unit_cell, f"{left_dev.name}_unit_rot90")
            unit_mirror = mirror_cell_horizontally(
                unit_cell, f"{left_dev.name}_unit_mirror"
            )
            unit_mirror_rot90 = rotate_cell_90(
                unit_mirror, f"{left_dev.name}_unit_mirror_rot90"
            )
            unit_pairs.append(
                {
                    "name": spec["name"],
                    "left": {"name": left_dev.name},
                    "right": {"name": right_dev.name},
                    "unit_count": unit_count,
                    "unit_w": unit_w,
                    "unit_h": unit_h,
                    "unit_cell": unit_cell,
                    "unit_rot90": unit_rot90,
                    "unit_mirror": unit_mirror,
                    "unit_mirror_rot90": unit_mirror_rot90,
                }
            )
            print(
                f"unit-level pair {spec['name']}: {left_dev.name}/{right_dev.name} "
                f"{unit_count} units each ({fingers_per_unit} fingers, "
                f"{unit_w:.2f}x{unit_h:.2f})"
            )
        else:
            left_cell, left_w, left_h = device_macro(left_dev)
            right_cell, right_w, right_h = device_macro(right_dev)
            right_cell = mirror_cell_horizontally(
                right_cell, f"{right_dev.name}_mirror"
            )
            macro_pairs.append(
                {
                    "name": spec["name"],
                    "left": _device_entry(left_dev.name, left_cell, left_w, left_h),
                    "right": _device_entry(right_dev.name, right_cell, right_w, right_h),
                }
            )

    singles = []
    for name in ("M9", "M12"):
        cell, width, height = device_macro(devices[name])
        singles.append(_device_entry(name, cell, width, height))

    matched_names = {name for spec in GROUPS for name in spec["devices"]}
    net_to_devices: dict[str, list[str]] = {}
    for name, dev in devices.items():
        for net in dev.nodes:
            net_to_devices.setdefault(net, []).append(name)

    power_nets = {"gnd", "vdd", "vss", "0"}
    hpwl_info: dict[str, dict[str, list[str]]] = {}
    for name in ("M9", "M12", "C1", "C2"):
        targets: dict[str, list[str]] = {}
        for net in devices[name].nodes:
            if net.lower() in power_nets:
                continue
            others = [
                other
                for other in net_to_devices.get(net, [])
                if other != name and other in matched_names
            ]
            if others:
                targets[net] = others
        hpwl_info[name] = targets

    # Initial unit-level placement: two units above and two below the device
    # centroid, leaving a central band for other devices.
    unit_state = {}
    placed: list[Rect] = []
    for pair in unit_pairs:
        unit_w, unit_h = pair["unit_w"], pair["unit_h"]

        if args.init == "checkerboard":
            left_offsets, right_offsets = _checkerboard_offsets(
                unit_w, unit_h, pair["unit_count"], gap
            )
        elif args.init == "central":
            left_offsets = _band_offsets(
                unit_w, unit_h, pair["unit_count"], gap
            )
            shifted = [
                (ox, oy + 0.5 * (unit_h + gap)) for ox, oy in left_offsets
            ]
            mx = sum(point[0] for point in shifted) / len(shifted)
            my = sum(point[1] for point in shifted) / len(shifted)
            right_offsets = [(x - mx, y - my) for x, y in shifted]
        else:
            left_offsets, right_offsets = _checkerboard_offsets(
                unit_w, unit_h, pair["unit_count"], gap
            )

            def jitter(points):
                pts = [
                    (
                        x + rng.uniform(-0.5 * unit_w, 0.5 * unit_w),
                        y + rng.uniform(-0.5 * unit_h, 0.5 * unit_h),
                    )
                    for x, y in points
                ]
                mx = sum(point[0] for point in pts) / len(pts)
                my = sum(point[1] for point in pts) / len(pts)
                return [(x - mx, y - my) for x, y in pts]

            left_offsets = jitter(left_offsets)
            right_offsets = jitter(right_offsets)

        state = (0.0, 0, left_offsets, right_offsets)
        unit_state[pair["name"]] = state
        pair_rects, _ = _unit_rects(pair, state)
        placed.extend(pair_rects)

    # Greedy placement for the remaining macro pairs and singles.
    macro_state = {}
    single_state = {}

    def bbox_score(rects):
        width = max(r.x1 for r in rects) - min(r.x0 for r in rects)
        height = max(r.y1 for r in rects) - min(r.y0 for r in rects)
        return width * height

    def candidates_y(height):
        values = {0.0, height / 2.0}
        for r in placed:
            values.update({r.y1 + gap + height / 2, r.y0 - gap - height / 2, (r.y0 + r.y1) / 2})
        return values

    def candidates_x(width):
        values = {0.0, width / 2, -width / 2}
        for r in placed:
            values.update({r.x1 + gap + width / 2, r.x0 - gap - width / 2, (r.x0 + r.x1) / 2})
        return values

    items = sorted(
        [("pair", p) for p in macro_pairs] + [("single", s) for s in singles],
        key=lambda item: (
            item[1]["left"]["w"] * item[1]["left"]["h"]
            + item[1]["right"]["w"] * item[1]["right"]["h"]
            if item[0] == "pair"
            else item[1]["w"] * item[1]["h"]
        ),
        reverse=True,
    )

    max_half_width = 45.0
    for kind, item in items:
        if kind == "pair":
            left, right = item["left"], item["right"]
            best = None
            for rotation in (0, 90):
                lw, lh = _dims(left, rotation)
                rw, rh = _dims(right, rotation)
                max_w = max(lw, rw)
                max_h = max(lh, rh)
                d_values = {max_w / 2 + gap / 2}
                for r in placed:
                    d_values.add(abs(r.x1) + gap + max_w / 2)
                    d_values.add(abs(r.x0) + gap + max_w / 2)
                for y in candidates_y(max_h):
                    for d in d_values:
                        if d > max_half_width - max_w / 2:
                            continue
                        lr = Rect(left["name"], -d - lw / 2, y - lh / 2,
                                  -d + lw / 2, y + lh / 2)
                        rr = Rect(right["name"], d - rw / 2, y - rh / 2,
                                  d + rw / 2, y + rh / 2)
                        if any(_overlap_area(lr, r, gap) > 0 or _overlap_area(rr, r, gap) > 0
                               for r in placed):
                            continue
                        if _overlap_area(lr, rr, gap) > 0:
                            continue
                        score = bbox_score(placed + [lr, rr])
                        if best is None or score < best[0]:
                            best = (score, rotation, d, y, lr, rr)
            if best is None:
                lw, lh = left["w"], left["h"]
                rw, rh = right["w"], right["h"]
                y = max((r.y1 for r in placed), default=0.0) + gap + max(lh, rh) / 2
                d = max(lw, rw) / 2 + gap / 2
                lr = Rect(left["name"], -d - lw / 2, y - lh / 2,
                          -d + lw / 2, y + lh / 2)
                rr = Rect(right["name"], d - rw / 2, y - rh / 2,
                          d + rw / 2, y + rh / 2)
                best = (0, 0, d, y, lr, rr)
            _, rotation, d, y, lr, rr = best
            placed.extend([lr, rr])
            macro_state[item["name"]] = (d, y, rotation)
        else:
            best = None
            for rotation in (0, 90):
                w, h = _dims(item, rotation)
                for y in candidates_y(h):
                    for x in candidates_x(w):
                        if abs(x) > max_half_width - w / 2:
                            continue
                        r = Rect(item["name"], x - w / 2, y - h / 2,
                                 x + w / 2, y + h / 2)
                        if any(_overlap_area(r, other, gap) > 0 for other in placed):
                            continue
                        score = bbox_score(placed + [r])
                        if best is None or score < best[0]:
                            best = (score, rotation, x, y, r)
            if best is None:
                w, h = item["w"], item["h"]
                y = max((r.y1 for r in placed), default=0.0) + gap + h / 2
                best = (0, 0, 0.0, y, Rect(item["name"], -w / 2, y - h / 2, w / 2, y + h / 2))
            _, rotation, x, y, r = best
            placed.append(r)
            single_state[item["name"]] = (x, y, rotation)

    lambda_hpwl = 0.5
    lambda_compact = 0.2
    max_width = 2.0 * max_half_width
    initial = _evaluate(
        unit_state, macro_state, single_state, unit_pairs, macro_pairs, singles,
        gap, hpwl_info, lambda_hpwl, lambda_compact, max_width,
    )

    current_unit = dict(unit_state)
    current_macro = dict(macro_state)
    current_single = dict(single_state)
    current_score = initial["score"]
    best = initial
    best_unit = dict(unit_state)
    best_macro = dict(macro_state)
    best_single = dict(single_state)

    iterations = args.iterations
    t_start = max(initial["score"] * 0.02, 1e-6)
    t_end = max(initial["score"] * 1e-5, 1e-9)
    temperature = t_start
    cooling = (t_end / t_start) ** (1.0 / max(iterations - 1, 1))

    for step in range(iterations):
        u_state = dict(current_unit)
        m_state = dict(current_macro)
        s_state = dict(current_single)
        scale = max(0.05, 8.0 * (1.0 - step / iterations))

        roll = rng.random()
        if unit_pairs and roll < 0.45:
            pair = rng.choice(unit_pairs)
            y, rotation, left_offsets, right_offsets = u_state[pair["name"]]
            left_offsets = list(left_offsets)
            right_offsets = list(right_offsets)

            if rng.random() < 0.15:
                # Rotate the whole pair; all units keep the same orientation.
                rotation = 90 if rotation == 0 else 0
            elif rng.random() < 0.35:
                y += rng.uniform(-scale, scale)
            else:
                # Centroid-preserving move inside one device: move two of
                # its unit offsets in opposite directions.  The other
                # device is untouched, so both device centroids stay at
                # their current points.
                offsets = left_offsets if rng.random() < 0.5 else right_offsets
                i, j = rng.sample(range(len(offsets)), 2)
                dx = rng.uniform(-scale, scale)
                dy = rng.uniform(-scale, scale)
                offsets[i] = (offsets[i][0] + dx, offsets[i][1] + dy)
                offsets[j] = (offsets[j][0] - dx, offsets[j][1] - dy)

            u_state[pair["name"]] = (y, rotation, left_offsets, right_offsets)
        elif macro_pairs and roll < 0.75:
            pair = rng.choice(macro_pairs)
            name = pair["name"]
            d, y, rotation = m_state[name]
            if rng.random() < 0.12:
                rotation = 90 if rotation == 0 else 0
            elif rng.random() < 0.5:
                lw, lh = _dims(pair["left"], rotation)
                rw, rh = _dims(pair["right"], rotation)
                d = max(max(lw, rw) / 2 + gap / 2, d + rng.uniform(-scale, scale))
            else:
                y += rng.uniform(-scale, scale)
            lw, lh = _dims(pair["left"], rotation)
            rw, rh = _dims(pair["right"], rotation)
            d = max(d, max(lw, rw) / 2 + gap / 2)
            m_state[name] = (d, y, rotation)
        elif singles:
            single = rng.choice(singles)
            name = single["name"]
            x, y, rotation = s_state[name]
            if rng.random() < 0.12:
                rotation = 90 if rotation == 0 else 0
            elif rng.random() < 0.5:
                x += rng.uniform(-scale, scale)
            else:
                y += rng.uniform(-scale, scale)
            s_state[name] = (x, y, rotation)

        candidate = _evaluate(
            u_state, m_state, s_state, unit_pairs, macro_pairs, singles,
            gap, hpwl_info, lambda_hpwl, lambda_compact, max_width,
        )
        delta = candidate["score"] - current_score
        if delta <= 0 or rng.random() < math.exp(-delta / max(temperature, 1e-12)):
            current_unit = u_state
            current_macro = m_state
            current_single = s_state
            current_score = candidate["score"]
            if candidate["overlap"] <= 1e-9 and candidate["score"] < best["score"]:
                best = candidate
                best_unit = dict(u_state)
                best_macro = dict(m_state)
                best_single = dict(s_state)
        temperature *= cooling

    rects = best["rects"]
    device_of_rect = best["device_of_rect"]
    centers = best["centers"]

    min_x = min(r.x0 for r in rects)
    min_y = min(r.y0 for r in rects)
    max_x = max(r.x1 for r in rects)
    max_y = max(r.y1 for r in rects)
    margin = 5.0
    top_width = max_x - min_x + 2 * margin
    top_height = max_y - min_y + 2 * margin
    shift_x = -min_x + margin
    shift_y = -min_y + margin

    top_lib = gdspy.GdsLibrary()
    top = gdspy.Cell(args.name)
    top_lib.add(top)
    top.add(gdspy.Rectangle((0, 0), (top_width, top_height), layer=99, datatype=0))

    cell_by_name = {}
    for pair in unit_pairs:
        _, rotation, _, _ = best_unit[pair["name"]]
        left_cell = pair["unit_rot90"] if rotation else pair["unit_cell"]
        right_cell = (
            pair["unit_mirror_rot90"] if rotation else pair["unit_mirror"]
        )
        for index in range(pair["unit_count"]):
            left_name = f"{pair['left']['name']}#{index}"
            right_name = f"{pair['right']['name']}#{index}"
            cell_by_name[left_name] = left_cell
            cell_by_name[right_name] = right_cell
    for pair in macro_pairs:
        _, _, rotation = best_macro[pair["name"]]
        cell_by_name[pair["left"]["name"]] = (
            pair["left"]["cell_rot90"] if rotation else pair["left"]["cell"]
        )
        cell_by_name[pair["right"]["name"]] = (
            pair["right"]["cell_rot90"] if rotation else pair["right"]["cell"]
        )
    for single in singles:
        _, _, rotation = best_single[single["name"]]
        cell_by_name[single["name"]] = (
            single["cell_rot90"] if rotation else single["cell"]
        )

    for rect in rects:
        cell = cell_by_name[rect.name]
        top.add(gdspy.CellReference(cell, origin=(rect.x0 + shift_x, rect.y0 + shift_y)))
        top.add(gdspy.Label(rect.name, (rect.x0 + shift_x + 0.4, rect.y0 + shift_y + 0.4), layer=100))
        top_lib.add(cell)

    gds_path = out_dir / f"{args.name}.gds"
    top_lib.write_gds(str(gds_path))
    info = render_gds(gds_path, out_dir / f"{args.name}.png", dpi=180)

    print("=== unit-level global symmetric layout ===")
    print(f"initial area : {initial['area']:.1f} um^2 "
          f"(overlap {initial['overlap']:.4f}, hpwl {initial['hpwl']:.2f})")
    print(f"final area   : {best['area']:.1f} um^2 "
          f"(overlap {best['overlap']:.4f}, hpwl {best['hpwl']:.2f})")
    print(f"top size     : {top_width:.2f} x {top_height:.2f} um "
          f"(area {top_width * top_height:.1f} um^2)")
    print(f"GDS          : {gds_path}")
    print(f"PNG          : {info['png']}")

    for pair in unit_pairs:
        left_center = centers[pair["left"]["name"]]
        right_center = centers[pair["right"]["name"]]
        print(
            f"  {pair['name']:28s} "
            f"left={left_center}, right={right_center}, "
            f"delta=({left_center[0] - right_center[0]:+.3f},"
            f"{left_center[1] - right_center[1]:+.3f})"
        )
    for pair in macro_pairs:
        left_center = centers[pair["left"]["name"]]
        right_center = centers[pair["right"]["name"]]
        print(
            f"  {pair['name']:28s} "
            f"centroid-sum=({left_center[0] + right_center[0]:+.3f},"
            f"{left_center[1] - right_center[1]:+.3f})"
        )

    # Group-colored schematic.
    pair_colors = plt.get_cmap("tab10").colors
    color_by_device = {}
    group_names = []
    for index, pair in enumerate(unit_pairs + macro_pairs):
        color = pair_colors[index % len(pair_colors)]
        color_by_device[pair["left"]["name"]] = color
        color_by_device[pair["right"]["name"]] = color
        group_names.append(pair["name"])
    for single in singles:
        color_by_device[single["name"]] = "#cccccc"

    fig, ax = plt.subplots(figsize=(max(6.0, top_width / 6.0), max(6.0, top_height / 6.0)))
    for rect in rects:
        device = device_of_rect[rect.name]
        ax.add_patch(
            MplRectangle(
                (rect.x0 + shift_x, rect.y0 + shift_y),
                rect.w, rect.h,
                facecolor=color_by_device[device],
                edgecolor="black", linewidth=0.8, alpha=0.75,
            )
        )
        ax.text(
            rect.x0 + shift_x + rect.w / 2.0,
            rect.y0 + shift_y + rect.h / 2.0,
            rect.name, ha="center", va="center", fontsize=6,
        )
    ax.set_xlim(0, top_width)
    ax.set_ylim(0, top_height)
    ax.set_aspect("equal")
    ax.axis("off")
    legend = [
        MplRectangle((0, 0), 1, 1, facecolor=pair_colors[i % len(pair_colors)],
                     edgecolor="black", alpha=0.75)
        for i in range(len(group_names))
    ] + [MplRectangle((0, 0), 1, 1, facecolor="#cccccc", edgecolor="black", alpha=0.75)]
    ax.legend(legend, group_names + ["unmatched"], loc="upper right", fontsize=6)
    schematic_path = out_dir / f"{args.name}_groups.png"
    fig.savefig(schematic_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"group schematic: {schematic_path}")

    summary = {
        "name": args.name,
        "init": args.init,
        "seed": args.seed,
        "iterations": args.iterations,
        "area": best["area"],
        "width": best["width"],
        "height": best["height"],
        "aspect": best["width"] / best["height"] if best["height"] else 0.0,
        "hpwl": best["hpwl"],
        "overlap": best["overlap"],
        "top_width": top_width,
        "top_height": top_height,
        "top_area": top_width * top_height,
        "gds": str(gds_path),
        "png": str(info["png"]),
        "groups_png": str(schematic_path),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"summary: {out_dir / 'summary.json'}")


if __name__ == "__main__":
    main()
