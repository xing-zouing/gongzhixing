#!/usr/bin/env python3
"""Step 9: simulated-annealing refinement of the symmetric layout.

State:
* each matching pair: ``(d, y)`` -> devices at ``(-d, y)`` and ``(+d, y)``;
* each unmatched device: free ``(x, y)``.

Objective:
* bounding-box area;
* overlap penalty (devices closer than the requested spacing);
* mild aspect-ratio penalty.

The greedy layout from Step 8 is used as the initial state; SA then moves
pairs/singles continuously to reduce area.  The best zero-overlap state is
written as GDS/PNG.

Run from the project root:

    .venv/bin/python scripts/ota3_step9_sa_refine.py
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
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
from cc.gds.device_cells import device_macro, mirror_cell_horizontally  # noqa: E402
from cc.gds.render import render_gds  # noqa: E402
from cc.io.netlist_parser import load_case  # noqa: E402
from cc.plan.ota3_groups import GROUPS  # noqa: E402


@dataclass
class Rect:
    name: str
    x0: float
    y0: float
    x1: float
    y1: float

    @property
    def w(self) -> float:
        return self.x1 - self.x0

    @property
    def h(self) -> float:
        return self.y1 - self.y0


def _overlap_area(a: Rect, b: Rect, gap: float) -> float:
    # Expand both rectangles by half the desired spacing.
    ax0, ay0, ax1, ay1 = a.x0 - gap / 2, a.y0 - gap / 2, a.x1 + gap / 2, a.y1 + gap / 2
    bx0, by0, bx1, by1 = b.x0 - gap / 2, b.y0 - gap / 2, b.x1 + gap / 2, b.y1 + gap / 2
    w = min(ax1, bx1) - max(ax0, bx0)
    h = min(ay1, by1) - max(ay0, by0)
    return max(0.0, w) * max(0.0, h)


def _device_entry(name: str, cell, width: float, height: float) -> dict:
    return {
        "name": name,
        "cell": cell,
        "cell_rot90": rotate_cell_90(cell, f"{name}_rot90"),
        "w": width,
        "h": height,
    }


def _dims(entry: dict, rotation: int) -> tuple[float, float]:
    if rotation:
        return entry["h"], entry["w"]
    return entry["w"], entry["h"]


def _build_rects(pair_state, single_state, pairs, singles) -> list[Rect]:
    rects = []
    for pair in pairs:
        d, y, rotation = pair_state[pair["name"]]
        left, right = pair["left"], pair["right"]
        lw, lh = _dims(left, rotation)
        rw, rh = _dims(right, rotation)
        rects.append(Rect(left["name"], -d - lw / 2, y - lh / 2,
                          -d + lw / 2, y + lh / 2))
        rects.append(Rect(right["name"], d - rw / 2, y - rh / 2,
                          d + rw / 2, y + rh / 2))
    for single in singles:
        x, y, rotation = single_state[single["name"]]
        w, h = _dims(single, rotation)
        rects.append(Rect(single["name"], x - w / 2, y - h / 2, x + w / 2, y + h / 2))
    return rects


def _evaluate(
    pair_state,
    single_state,
    pairs,
    singles,
    gap: float,
    hpwl_info: dict | None = None,
    lambda_hpwl: float = 0.0,
    max_width: float | None = None,
    lambda_compact: float = 0.0,
):
    rects = _build_rects(pair_state, single_state, pairs, singles)
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
    if hpwl_info:
        rect_by_name = {r.name: r for r in rects}
        for name, nets in hpwl_info.items():
            rect = rect_by_name.get(name)
            if rect is None:
                continue
            cx = (rect.x0 + rect.x1) / 2.0
            cy = (rect.y0 + rect.y1) / 2.0
            for others in nets.values():
                points = [(cx, cy)]
                for other in others:
                    other_rect = rect_by_name.get(other)
                    if other_rect is not None:
                        points.append(
                            (
                                (other_rect.x0 + other_rect.x1) / 2.0,
                                (other_rect.y0 + other_rect.y1) / 2.0,
                            )
                        )
                if len(points) > 1:
                    xs = [point[0] for point in points]
                    ys = [point[1] for point in points]
                    hpwl += max(xs) - min(xs) + max(ys) - min(ys)

    center_x = (min_x + max_x) / 2.0
    center_y = (min_y + max_y) / 2.0
    compactness = 0.0
    for rect in rects:
        rx = (rect.x0 + rect.x1) / 2.0
        ry = (rect.y0 + rect.y1) / 2.0
        compactness += (rx - center_x) ** 2 + (ry - center_y) ** 2

    aspect = width / height if height else 1.0
    # A 1:1 target is often impossible here because the widest mirror pair
    # sets a lower bound on the layout width.  Target a practical 2:1 shape
    # and penalise deviations more strongly than before.
    target_aspect = 1.0
    aspect_penalty = 1.0 + 0.8 * abs(math.log(max(aspect, 1e-9) / target_aspect))
    width_penalty = 0.0
    if max_width is not None and width > max_width:
        width_penalty = 500.0 * (width - max_width) ** 2
    score = (
        width * height * aspect_penalty
        + 1000.0 * overlap
        + lambda_hpwl * hpwl
        + width_penalty
        + lambda_compact * compactness
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
    }


def _greedy_initial(
    pairs,
    singles,
    gap: float,
    max_half_width: float | None = None,
    allow_rotation: bool = True,
):
    """Greedy largest-first initial layout (same idea as Step 8)."""
    rotation_options = (0, 90) if allow_rotation else (0,)
    placed: list[Rect] = []
    pair_state = {}
    single_state = {}

    def bbox_score(rects):
        w = max(r.x1 for r in rects) - min(r.x0 for r in rects)
        h = max(r.y1 for r in rects) - min(r.y0 for r in rects)
        return w * h

    def candidates_y(h):
        values = {0.0, h / 2.0}
        for r in placed:
            values.update({r.y1 + gap + h / 2, r.y0 - gap - h / 2, (r.y0 + r.y1) / 2})
        return values

    def candidates_x(w):
        values = {0.0, w / 2, -w / 2}
        for r in placed:
            values.update({r.x1 + gap + w / 2, r.x0 - gap - w / 2, (r.x0 + r.x1) / 2})
        return values

    items = sorted(
        [("pair", p) for p in pairs] + [("single", s) for s in singles],
        key=lambda item: (
            item[1]["left"]["w"] * item[1]["left"]["h"]
            + item[1]["right"]["w"] * item[1]["right"]["h"]
            if item[0] == "pair"
            else item[1]["w"] * item[1]["h"]
        ),
        reverse=True,
    )

    for kind, item in items:
        if kind == "pair":
            left, right = item["left"], item["right"]
            best = None
            for rotation in rotation_options:
                lw, lh = _dims(left, rotation)
                rw, rh = _dims(right, rotation)
                max_w = max(lw, rw)
                max_h = max(lh, rh)

                d_values = {max_w / 2 + gap / 2}
                for r in placed:
                    d_values.add(abs(r.x1) + gap + max_w / 2)
                    d_values.add(abs(r.x0) + gap + max_w / 2)
                max_d = (
                    max_half_width - max_w / 2.0
                    if max_half_width is not None
                    else None
                )
                for y in candidates_y(max_h):
                    for d in d_values:
                        if max_d is not None and d > max_d + 1e-9:
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
                max_w = max(lw, rw)
                max_h = max(lh, rh)
                y = max((r.y1 for r in placed), default=0.0) + gap + max_h / 2
                d = max_w / 2 + gap / 2
                lr = Rect(left["name"], -d - lw / 2, y - lh / 2,
                          -d + lw / 2, y + lh / 2)
                rr = Rect(right["name"], d - rw / 2, y - rh / 2,
                          d + rw / 2, y + rh / 2)
                best = (0, 0, d, y, lr, rr)
            _, rotation, d, y, lr, rr = best
            placed.extend([lr, rr])
            pair_state[item["name"]] = (d, y, rotation)
        else:
            best = None
            for rotation in rotation_options:
                w, h = _dims(item, rotation)
                for y in candidates_y(h):
                    for x in candidates_x(w):
                        if max_half_width is not None and abs(x) > max_half_width - w / 2 + 1e-9:
                            continue
                        r = Rect(item["name"], x - w / 2, y - h / 2, x + w / 2, y + h / 2)
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

    return pair_state, single_state


def main() -> None:
    out_dir = ROOT / "build" / "ota3" / "final"
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in out_dir.glob("ota3_global_final.*"):
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

        left_dev, right_dev = devs
        left_cell, left_w, left_h = device_macro(left_dev)
        right_cell, right_w, right_h = device_macro(right_dev)
        right_cell = mirror_cell_horizontally(
            right_cell, f"{right_dev.name}_mirror"
        )
        pairs.append(
            {
                "name": spec["name"],
                "left": _device_entry(left_dev.name, left_cell, left_w, left_h),
                "right": _device_entry(right_dev.name, right_cell, right_w, right_h),
            }
        )
        print(
            f"{spec['name']}: {left_dev.name} multi={left_dev.multi} "
            f"({left_w:.2f}x{left_h:.2f}), {right_dev.name} multi={right_dev.multi} "
            f"({right_w:.2f}x{right_h:.2f})"
        )

    singles = []
    for name in ("M9", "M12"):
        cell, width, height = device_macro(devices[name])
        singles.append(_device_entry(name, cell, width, height))

    # Connectivity-aware cost for unmatched devices: connect M9 to net19
    # (M0/M7) and M12 to net17 (M10/M11).  Power/ground nets are ignored.
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
        print(f"HPWL targets {name}: {targets}")

    lambda_hpwl = 2.0
    lambda_compact = 1.0
    spacing = 0.6

    # The widest pair sets a lower bound on the layout width; use it as a
    # shelf width so the initial layout stacks rows instead of becoming a
    # very long strip.
    max_half_width = max(
        max(
            pair["left"]["w"],
            pair["left"]["h"],
            pair["right"]["w"],
            pair["right"]["h"],
        )
        for pair in pairs
    ) + spacing
    max_width = 2.0 * max_half_width
    print(f"max layout width (shelf): {max_width:.2f} um")

    rot_state = _greedy_initial(
        pairs, singles, spacing, max_half_width=max_half_width, allow_rotation=True
    )
    rot_initial = _evaluate(
        *rot_state, pairs, singles, spacing,
        hpwl_info=hpwl_info, lambda_hpwl=lambda_hpwl, max_width=max_width,
        lambda_compact=lambda_compact,
    )
    flat_state = _greedy_initial(
        pairs, singles, spacing, max_half_width=max_half_width, allow_rotation=False
    )
    flat_initial = _evaluate(
        *flat_state, pairs, singles, spacing,
        hpwl_info=hpwl_info, lambda_hpwl=lambda_hpwl, max_width=max_width,
        lambda_compact=lambda_compact,
    )
    if flat_initial["score"] <= rot_initial["score"]:
        pair_state, single_state = flat_state
        initial = flat_initial
        print("initial: rotation disabled (better greedy start)")
    else:
        pair_state, single_state = rot_state
        initial = rot_initial
        print("initial: rotation enabled (better greedy start)")

    current_pair = dict(pair_state)
    current_single = dict(single_state)
    current_score = initial["score"]

    best = initial
    best_pair = dict(pair_state)
    best_single = dict(single_state)

    rng = random.Random(0)
    iterations = 40000
    t_start = max(initial["score"] * 0.02, 1e-6)
    t_end = max(initial["score"] * 1e-5, 1e-9)
    temperature = t_start
    cooling = (t_end / t_start) ** (1.0 / max(iterations - 1, 1))

    for step in range(iterations):
        p_state = dict(current_pair)
        s_state = dict(current_single)
        scale = max(0.05, 8.0 * (1.0 - step / iterations))

        if rng.random() < 0.7 and pairs:
            pair = rng.choice(pairs)
            name = pair["name"]
            d, y, rotation = p_state[name]

            if rng.random() < 0.12:
                # Rotate both devices of the pair together; matching devices
                # must keep the same orientation.
                rotation = 90 if rotation == 0 else 0
            elif rng.random() < 0.5:
                lw, lh = _dims(pair["left"], rotation)
                rw, rh = _dims(pair["right"], rotation)
                d = max(max(lw, rw) / 2 + spacing / 2,
                        d + rng.uniform(-scale, scale))
            else:
                y += rng.uniform(-scale, scale)

            lw, lh = _dims(pair["left"], rotation)
            rw, rh = _dims(pair["right"], rotation)
            d = max(d, max(lw, rw) / 2 + spacing / 2)
            p_state[name] = (d, y, rotation)
        elif singles:
            single = rng.choice(singles)
            name = single["name"]
            x, y, rotation = s_state[name]

            if rng.random() < 0.12:
                rotation = 90 if rotation == 0 else 0
            else:
                # With some probability, jump next to a device connected by
                # a relevant net; otherwise perturb locally.
                targets = [
                    other
                    for others in hpwl_info.get(name, {}).values()
                    for other in others
                ]
                if targets and rng.random() < 0.35:
                    current_rects = _build_rects(current_pair, current_single, pairs, singles)
                    rect_by_name = {rect.name: rect for rect in current_rects}
                    valid_targets = [other for other in targets if other in rect_by_name]
                    if valid_targets:
                        target = rect_by_name[rng.choice(valid_targets)]
                        x = (target.x0 + target.x1) / 2.0 + rng.uniform(-6.0, 6.0)
                        y = (target.y0 + target.y1) / 2.0 + rng.uniform(-6.0, 6.0)
                    else:
                        x += rng.uniform(-scale, scale)
                        y += rng.uniform(-scale, scale)
                else:
                    if rng.random() < 0.5:
                        x += rng.uniform(-scale, scale)
                    else:
                        y += rng.uniform(-scale, scale)
            s_state[name] = (x, y, rotation)

        candidate = _evaluate(
            p_state, s_state, pairs, singles, spacing,
            hpwl_info=hpwl_info, lambda_hpwl=lambda_hpwl, max_width=max_width,
            lambda_compact=lambda_compact,
        )
        delta = candidate["score"] - current_score
        if delta <= 0 or rng.random() < math.exp(-delta / max(temperature, 1e-12)):
            current_pair = p_state
            current_single = s_state
            current_score = candidate["score"]
            if candidate["overlap"] <= 1e-9 and candidate["score"] < best["score"]:
                best = candidate
                best_pair = dict(p_state)
                best_single = dict(s_state)
        temperature *= cooling

    # The greedy start is valid; if SA never found a better valid state,
    # it remains the best.
    if best["overlap"] > 1e-9:
        best = initial
        best_pair = dict(pair_state)
        best_single = dict(single_state)

    rects = _build_rects(best_pair, best_single, pairs, singles)
    min_x = min(r.x0 for r in rects)
    min_y = min(r.y0 for r in rects)
    max_x = max(r.x1 for r in rects)
    max_y = max(r.y1 for r in rects)
    margin = 5.0
    top_width = max_x - min_x + 2 * margin
    top_height = max_y - min_y + 2 * margin
    shift_x = -min_x + margin
    shift_y = -min_y + margin

    cell_by_name = {}
    for pair in pairs:
        _, _, rotation = best_pair[pair["name"]]
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

    top_lib = gdspy.GdsLibrary()
    top = gdspy.Cell("ota3_global_final")
    top_lib.add(top)
    top.add(gdspy.Rectangle((0, 0), (top_width, top_height), layer=99, datatype=0))
    for rect in rects:
        cell = cell_by_name[rect.name]
        top.add(gdspy.CellReference(cell, origin=(rect.x0 + shift_x, rect.y0 + shift_y)))
        top.add(gdspy.Label(rect.name, (rect.x0 + shift_x + 0.4, rect.y0 + shift_y + 0.4), layer=100))
        top_lib.add(cell)

    gds_path = out_dir / "ota3_global_final.gds"
    top_lib.write_gds(str(gds_path))
    info = render_gds(gds_path, out_dir / "ota3_global_final.png", dpi=180)

    print("=== final global symmetric layout (multi expanded + HPWL) ===")
    print(f"initial area : {initial['area']:.1f} um^2 "
          f"(overlap {initial['overlap']:.4f}, hpwl {initial['hpwl']:.2f}, "
          f"compact {initial['compactness']:.1f})")
    print(f"final area   : {best['area']:.1f} um^2 "
          f"(overlap {best['overlap']:.4f}, hpwl {best['hpwl']:.2f}, "
          f"compact {best['compactness']:.1f})")
    print(f"top size     : {top_width:.2f} x {top_height:.2f} um "
          f"(area {top_width * top_height:.1f} um^2)")
    print(f"GDS          : {gds_path}")
    print(f"PNG          : {info['png']}")

    # Group-colored schematic: both devices of a matching pair share a color;
    # unmatched devices are neutral gray.
    pair_colors = plt.get_cmap("tab10").colors
    color_by_device = {}
    for index, pair in enumerate(pairs):
        color = pair_colors[index % len(pair_colors)]
        color_by_device[pair["left"]["name"]] = color
        color_by_device[pair["right"]["name"]] = color
    for single in singles:
        color_by_device[single["name"]] = "#cccccc"

    fig_w = max(6.0, top_width / 6.0)
    fig_h = max(6.0, top_height / 6.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    for rect in rects:
        ax.add_patch(
            MplRectangle(
                (rect.x0 + shift_x, rect.y0 + shift_y),
                rect.w,
                rect.h,
                facecolor=color_by_device[rect.name],
                edgecolor="black",
                linewidth=0.8,
                alpha=0.75,
            )
        )
        ax.text(
            rect.x0 + shift_x + rect.w / 2.0,
            rect.y0 + shift_y + rect.h / 2.0,
            rect.name,
            ha="center",
            va="center",
            fontsize=7,
        )
    ax.set_xlim(0, top_width)
    ax.set_ylim(0, top_height)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("matching groups: same color per pair; unmatched gray")
    legend = [
        MplRectangle((0, 0), 1, 1, facecolor=pair_colors[i % len(pair_colors)],
                     edgecolor="black", alpha=0.75)
        for i in range(len(pairs))
    ] + [MplRectangle((0, 0), 1, 1, facecolor="#cccccc", edgecolor="black", alpha=0.75)]
    legend_labels = [pair["name"] for pair in pairs] + ["unmatched"]
    ax.legend(legend, legend_labels, loc="upper right", fontsize=6)
    schematic_path = out_dir / "ota3_global_final_groups.png"
    fig.savefig(schematic_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"group schematic: {schematic_path}")

    for pair in pairs:
        left = next(r for r in rects if r.name == pair["left"]["name"])
        right = next(r for r in rects if r.name == pair["right"]["name"])
        print(
            f"  {pair['name']:28s} "
            f"center-x-sum={(left.x0 + left.x1 + right.x0 + right.x1) / 2:+.3f} "
            f"y-diff={(left.y0 + left.y1 - right.y0 - right.y1) / 2:+.3f}"
        )

    print("  rotations:")
    for pair in pairs:
        print(f"    {pair['name']:28s} {best_pair[pair['name']][2]:3d} deg")
    for single in singles:
        print(f"    {single['name']:28s} {best_single[single['name']][2]:3d} deg")

    rect_by_name = {r.name: r for r in rects}
    for name, nets in hpwl_info.items():
        rect = rect_by_name[name]
        cx = (rect.x0 + rect.x1) / 2.0
        cy = (rect.y0 + rect.y1) / 2.0
        print(f"  device {name}: center=({cx:.2f}, {cy:.2f})")
        for net, others in nets.items():
            for other in others:
                other_rect = rect_by_name[other]
                ox = (other_rect.x0 + other_rect.x1) / 2.0
                oy = (other_rect.y0 + other_rect.y1) / 2.0
                print(f"    {net} -> {other}: manhattan={abs(cx - ox) + abs(cy - oy):.2f}")


if __name__ == "__main__":
    main()
