#!/usr/bin/env python3
"""Step 8: largest-first, symmetry-aware greedy packing.

Each matching pair is still a left/right mirror about the global vertical
axis at x = 0, with both devices at the same y.  Unlike Step 7, a pair does
not own a row: the packer tries candidate y positions and symmetric
distances in the existing layout, and picks the placement with the smallest
resulting bounding-box area.  Large items are placed first; small unmatched
devices fill the remaining gaps.

Run from the project root:

    .venv/bin/python scripts/ota3_step8_greedy_symmetric.py
"""

from __future__ import annotations

from dataclasses import dataclass
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


@dataclass
class Rect:
    name: str
    x0: float
    y0: float
    x1: float
    y1: float

    @property
    def width(self) -> float:
        return self.x1 - self.x0

    @property
    def height(self) -> float:
        return self.y1 - self.y0


def _overlap(a: Rect, b: Rect, gap: float = 0.0) -> bool:
    return not (
        a.x1 + gap <= b.x0
        or b.x1 + gap <= a.x0
        or a.y1 + gap <= b.y0
        or b.y1 + gap <= a.y0
    )


def _bbox(rects: list[Rect]) -> tuple[float, float, float, float]:
    return (
        min(r.x0 for r in rects),
        min(r.y0 for r in rects),
        max(r.x1 for r in rects),
        max(r.y1 for r in rects),
    )


def _score(rects: list[Rect]) -> float:
    x0, y0, x1, y1 = _bbox(rects)
    width = x1 - x0
    height = y1 - y0
    aspect = width / height if height else 1.0
    return width * height * (1.0 + 0.15 * abs(aspect - 1.0))


def main() -> None:
    out_dir = ROOT / "build" / "ota3" / "greedy_symmetric"
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in out_dir.glob("ota3_greedy_symmetric.*"):
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
        pairs.append(
            {
                "kind": "pair",
                "name": spec["name"],
                "left": cells[0],
                "right": cells[1],
                "area": cells[0]["w"] * cells[0]["h"] + cells[1]["w"] * cells[1]["h"],
            }
        )

    singles = []
    for name in ("M9", "M12"):
        cell, width, height = device_cell(devices[name])
        singles.append(
            {
                "kind": "single",
                "name": name,
                "cell": cell,
                "w": width,
                "h": height,
                "area": width * height,
            }
        )

    items = sorted(pairs + singles, key=lambda item: item["area"], reverse=True)
    gap = 2.0
    placed: list[Rect] = []
    pair_positions: dict[str, tuple[float, float]] = {}

    def candidate_ys(height: float) -> set[float]:
        values = {0.0, height / 2.0}
        for rect in placed:
            values.add(rect.y1 + gap + height / 2.0)
            values.add(rect.y0 - gap - height / 2.0)
            values.add((rect.y0 + rect.y1) / 2.0)
        return values

    def candidate_xs(width: float) -> set[float]:
        values = {0.0, width / 2.0, -width / 2.0}
        for rect in placed:
            values.add(rect.x1 + gap + width / 2.0)
            values.add(rect.x0 - gap - width / 2.0)
            values.add((rect.x0 + rect.x1) / 2.0)
        return values

    for item in items:
        if item["kind"] == "pair":
            left = item["left"]
            right = item["right"]
            max_w = max(left["w"], right["w"])
            max_h = max(left["h"], right["h"])

            # Candidate symmetric distances: outside existing rectangles and
            # a coarse grid to fill larger gaps.
            d_values = {max_w / 2.0 + gap / 2.0}
            for rect in placed:
                d_values.add(abs(rect.x1) + gap + max_w / 2.0)
                d_values.add(abs(rect.x0) + gap + max_w / 2.0)
            max_extent = max((abs(r.x0) for r in placed), default=0.0)
            step = max(max_w / 2.0 + gap, 1.0)
            x = max_w / 2.0 + gap / 2.0
            while x <= max_extent + max_w:
                d_values.add(x)
                x += step

            best = None
            for cy in candidate_ys(max_h):
                for d in sorted(v for v in d_values if v > 0):
                    left_rect = Rect(
                        left["name"],
                        -d - left["w"] / 2.0,
                        cy - left["h"] / 2.0,
                        -d + left["w"] / 2.0,
                        cy + left["h"] / 2.0,
                    )
                    right_rect = Rect(
                        right["name"],
                        d - right["w"] / 2.0,
                        cy - right["h"] / 2.0,
                        d + right["w"] / 2.0,
                        cy + right["h"] / 2.0,
                    )
                    if _overlap(left_rect, right_rect, gap):
                        continue
                    if any(_overlap(left_rect, r, gap) or _overlap(right_rect, r, gap) for r in placed):
                        continue
                    score = _score(placed + [left_rect, right_rect])
                    if best is None or score < best[0]:
                        best = (score, d, cy, left_rect, right_rect)

            if best is None:
                # Fallback: put the pair above everything.
                y1 = max((r.y1 for r in placed), default=0.0)
                cy = y1 + gap + max_h / 2.0
                d = max_w / 2.0 + gap / 2.0
                left_rect = Rect(
                    left["name"], -d - left["w"] / 2.0, cy - left["h"] / 2.0,
                    -d + left["w"] / 2.0, cy + left["h"] / 2.0,
                )
                right_rect = Rect(
                    right["name"], d - right["w"] / 2.0, cy - right["h"] / 2.0,
                    d + right["w"] / 2.0, cy + right["h"] / 2.0,
                )
                best = (0.0, d, cy, left_rect, right_rect)

            _, d, cy, left_rect, right_rect = best
            placed.extend([left_rect, right_rect])
            pair_positions[item["name"]] = (d, cy)
        else:
            width, height = item["w"], item["h"]
            best = None
            for cy in candidate_ys(height):
                for cx in candidate_xs(width):
                    rect = Rect(
                        item["name"],
                        cx - width / 2.0,
                        cy - height / 2.0,
                        cx + width / 2.0,
                        cy + height / 2.0,
                    )
                    if any(_overlap(rect, r, gap) for r in placed):
                        continue
                    score = _score(placed + [rect])
                    if best is None or score < best[0]:
                        best = (score, rect)
            if best is None:
                y1 = max((r.y1 for r in placed), default=0.0)
                rect = Rect(item["name"], -width / 2.0, y1 + gap, width / 2.0, y1 + gap + height)
                best = (0.0, rect)
            placed.append(best[1])

    min_x, min_y, max_x, max_y = _bbox(placed)
    margin = 5.0
    top_width = max_x - min_x + 2 * margin
    top_height = max_y - min_y + 2 * margin
    shift_x = -min_x + margin
    shift_y = -min_y + margin

    device_cells = {}
    for pair in pairs:
        device_cells[pair["left"]["name"]] = pair["left"]
        device_cells[pair["right"]["name"]] = pair["right"]
    for single in singles:
        device_cells[single["name"]] = single

    top_lib = gdspy.GdsLibrary()
    top = gdspy.Cell("ota3_greedy_symmetric")
    top_lib.add(top)
    top.add(gdspy.Rectangle((0, 0), (top_width, top_height), layer=99, datatype=0))

    for rect in placed:
        cell_info = device_cells[rect.name]
        top.add(
            gdspy.CellReference(
                cell_info["cell"],
                origin=(rect.x0 + shift_x, rect.y0 + shift_y),
            )
        )
        top.add(gdspy.Label(rect.name, (rect.x0 + shift_x + 0.4, rect.y0 + shift_y + 0.4), layer=100))
        top_lib.add(cell_info["cell"])

    gds_path = out_dir / "ota3_greedy_symmetric.gds"
    top_lib.write_gds(str(gds_path))
    info = render_gds(gds_path, out_dir / "ota3_greedy_symmetric.png", dpi=180)

    print("=== largest-first symmetric greedy placement ===")
    print(f"top size: {top_width:.2f} x {top_height:.2f} um "
          f"(area {top_width * top_height:.1f} um^2)")
    print(f"GDS: {gds_path}")
    print(f"PNG: {info['png']}")

    for pair in pairs:
        x0, x1 = pair["left"]["name"], pair["right"]["name"]
        left_rect = next(r for r in placed if r.name == x0)
        right_rect = next(r for r in placed if r.name == x1)
        print(
            f"  {pair['name']:28s} "
            f"center-x-sum={(left_rect.x0 + left_rect.x1 + right_rect.x0 + right_rect.x1) / 2:+.3f} "
            f"y-diff={(left_rect.y0 + left_rect.y1 - right_rect.y0 - right_rect.y1) / 2:+.3f}"
        )


if __name__ == "__main__":
    main()
