"""Heterogeneous, symmetry-preserving array writer with real dummy fill.

Compared with ``generate_mixed_group_gds`` (uniform max-cell grid), this
writer:

* sizes each column by the widest cell in that column and each row by the
  tallest cell in that row;
* symmetrizes column widths and row heights about the array center, so
  mirrored logical positions remain physically mirrored;
* fills every ``0`` slot with the unit cell of the nearest real device
  instead of leaving an empty grid position.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import gdspy
import numpy as np

from cc.gds.cell_ops import rotate_cell_90
from cc.gds.compaction import compact_placements, pair_move_compact
from cc.gds.mixed_writer import _make_cell


def _device_type(dev: dict) -> tuple:
    """Electrical/environment type used for dummy selection."""
    if dev["kind"] == "cap":
        return ("cap",)
    if dev["kind"] == "mos":
        return ("mos", bool(dev.get("nch")), tuple(sorted(dev.get("attr", []))))
    return (dev["kind"],)


def _neighborhood_positions(arr: np.ndarray, row: int, col: int, radius: int):
    positions = np.argwhere(arr != 0)
    return [
        tuple(pos)
        for pos in positions
        if max(abs(pos[0] - row), abs(pos[1] - col)) <= radius
    ]


def _choose_dummy_label(
    arr: np.ndarray,
    devices: list[dict],
    row: int,
    col: int,
    radius: int = 1,
) -> int:
    """Choose a dummy by neighborhood majority *device type*, then nearest."""
    positions = []
    for current_radius in range(radius, max(arr.shape) + 1):
        positions = _neighborhood_positions(arr, row, col, current_radius)
        if positions:
            break
    if not positions:
        raise ValueError("matrix has no real devices")

    type_counts = Counter(
        _device_type(devices[int(arr[r, c]) - 1]) for r, c in positions
    )
    best_type = max(
        type_counts,
        key=lambda candidate: (type_counts[candidate], tuple(map(str, candidate))),
    )
    candidates = [
        (r, c)
        for r, c in positions
        if _device_type(devices[int(arr[r, c]) - 1]) == best_type
    ]
    candidates.sort(key=lambda pos: (abs(pos[0] - row) + abs(pos[1] - col), pos))
    r, c = candidates[0]
    return int(arr[r, c])


def generate_mixed_table_gds(
    devices: list[dict],
    matrix,
    output_path: str | Path,
    array_name: str = "mixed_table",
    gap: float = 0.3,
    fill_dummies: bool = True,
    groups=None,
    compact: bool = False,
    compact_min_gap: float = 0.0,
    compact_mode: str = "scale",
    compact_iterations: int = 200_000,
    compact_seed: int = 0,
    square_reservation: bool = False,
    rotation: str = "horizontal",
) -> dict:
    arr = np.asarray(matrix, dtype=int)
    if arr.ndim == 1:
        side = int(np.ceil(np.sqrt(arr.size)))
        arr = np.pad(arr, (0, side * side - arr.size)).reshape(side, side)
    rows, cols = arr.shape

    top_lib = gdspy.GdsLibrary()
    cell_info = {}
    for group_id, dev in enumerate(devices, start=1):
        cell, width, height = _make_cell(dev)
        rotated = rotate_cell_90(cell, f"{cell.name}_rot90")
        cell_info[group_id] = {
            "cell": cell,
            "cell_rot90": rotated,
            "w": width,
            "h": height,
            "side": max(width, height),
        }
        top_lib.add(cell)
        top_lib.add(rotated)

    # Per-slot cell choice.  Dummy slots use the majority device *type* in
    # the local 3x3 neighborhood, then the nearest device of that type.
    slot_label = np.zeros((rows, cols), dtype=int)
    for row in range(rows):
        for col in range(cols):
            label = int(arr[row, col])
            if label == 0 and fill_dummies:
                label = _choose_dummy_label(arr, devices, row, col, radius=1)
            slot_label[row, col] = label

    # Mirror-symmetric dummy types: if two mirrored 0 slots chose different
    # types, pick the majority type in the union of their neighborhoods and
    # use one representative for both.
    if fill_dummies:
        for row in range(rows):
            for col in range(cols):
                if arr[row, col] != 0:
                    continue
                mirror_row, mirror_col = rows - 1 - row, cols - 1 - col
                if arr[mirror_row, mirror_col] != 0:
                    continue
                if (row, col) > (mirror_row, mirror_col):
                    continue

                first = int(slot_label[row, col])
                second = int(slot_label[mirror_row, mirror_col])
                if _device_type(devices[first - 1]) == _device_type(devices[second - 1]):
                    continue

                positions = set(
                    _neighborhood_positions(arr, row, col, radius=2)
                ) | set(
                    _neighborhood_positions(arr, mirror_row, mirror_col, radius=2)
                )
                if not positions:
                    continue

                type_counts = Counter(
                    _device_type(devices[int(arr[r, c]) - 1])
                    for r, c in positions
                )
                best_type = max(
                    type_counts,
                    key=lambda candidate: (
                        type_counts[candidate],
                        tuple(map(str, candidate)),
                    ),
                )
                candidates = [
                    (r, c)
                    for r, c in positions
                    if _device_type(devices[int(arr[r, c]) - 1]) == best_type
                ]
                if not candidates:
                    continue
                mid_row = (row + mirror_row) / 2.0
                mid_col = (col + mirror_col) / 2.0
                candidates.sort(
                    key=lambda pos: (
                        abs(pos[0] - mid_row) + abs(pos[1] - mid_col),
                        pos,
                    )
                )
                representative = int(arr[candidates[0]])
                slot_label[row, col] = representative
                slot_label[mirror_row, mirror_col] = representative

    if rotation not in ("horizontal", "vertical", "auto"):
        raise ValueError(f"unsupported rotation mode: {rotation}")
    if rotation == "auto" and not square_reservation:
        raise ValueError("rotation='auto' requires square_reservation=True")

    def _reservation(label):
        info = cell_info[label]
        if square_reservation:
            return info["side"], info["side"]
        if rotation == "vertical":
            return info["h"], info["w"]
        return info["w"], info["h"]

    col_widths = [
        max(_reservation(int(slot_label[r, c]))[0] for r in range(rows))
        for c in range(cols)
    ]
    row_heights = [
        max(_reservation(int(slot_label[r, c]))[1] for c in range(cols))
        for r in range(rows)
    ]

    # Symmetrize the bands about the array center.
    for i in range(cols // 2):
        j = cols - 1 - i
        width = max(col_widths[i], col_widths[j])
        col_widths[i] = col_widths[j] = width
    for i in range(rows // 2):
        j = rows - 1 - i
        height = max(row_heights[i], row_heights[j])
        row_heights[i] = row_heights[j] = height

    total_w = sum(col_widths) + gap * (cols - 1)
    total_h = sum(row_heights) + gap * (rows - 1)
    start_x = -total_w / 2.0
    start_y = -total_h / 2.0

    col_centers = []
    cursor = start_x
    for width in col_widths:
        col_centers.append(cursor + width / 2.0)
        cursor += width + gap
    row_centers = []
    cursor = start_y
    for height in row_heights:
        row_centers.append(cursor + height / 2.0)
        cursor += height + gap

    placements = []
    for row in range(rows):
        for col in range(cols):
            label = int(slot_label[row, col])
            info = cell_info[label]
            device_name = devices[label - 1]["name"]
            is_dummy = bool(arr[row, col] == 0)
            label_text = f"D({device_name})" if is_dummy else device_name
            res_w, res_h = _reservation(label)
            rotation_value = 90 if rotation == "vertical" else 0
            if rotation_value == 90:
                collision_w, collision_h = info["h"], info["w"]
            else:
                collision_w, collision_h = info["w"], info["h"]
            placements.append(
                {
                    "row": row,
                    "col": col,
                    "label": label,
                    "name": device_name,
                    "w": res_w,
                    "h": res_h,
                    "collision_w": collision_w,
                    "collision_h": collision_h,
                    "actual_w": info["w"],
                    "actual_h": info["h"],
                    "rotation": rotation_value,
                    "dummy": is_dummy,
                    "cx": col_centers[col],
                    "cy": row_centers[rows - 1 - row],
                    "text": label_text,
                }
            )

    compact_stats = None
    if compact:
        if compact_mode == "pairs":
            placements, compact_stats = pair_move_compact(
                placements,
                groups=groups,
                min_gap=compact_min_gap,
                iterations=compact_iterations,
                seed=compact_seed,
            )
        elif compact_mode == "scale":
            placements, compact_stats = compact_placements(
                placements,
                groups=groups,
                min_gap=compact_min_gap,
            )
        else:
            raise ValueError(f"unsupported compact_mode: {compact_mode}")

    def _placed_dims(placement):
        if placement.get("rotation", 0) == 90:
            return placement["actual_h"], placement["actual_w"]
        return placement["actual_w"], placement["actual_h"]

    def _actual_bbox(pls):
        left = min(p["cx"] - _placed_dims(p)[0] / 2.0 for p in pls)
        right = max(p["cx"] + _placed_dims(p)[0] / 2.0 for p in pls)
        bottom = min(p["cy"] - _placed_dims(p)[1] / 2.0 for p in pls)
        top = max(p["cy"] + _placed_dims(p)[1] / 2.0 for p in pls)
        return left, bottom, right, top

    def _actual_overlap(pls, gap):
        for first in range(len(pls)):
            first_w, first_h = _placed_dims(pls[first])
            for second in range(first + 1, len(pls)):
                second_w, second_h = _placed_dims(pls[second])
                if (
                    abs(pls[first]["cx"] - pls[second]["cx"])
                    < (first_w + second_w) / 2.0 + gap - 1e-9
                    and abs(pls[first]["cy"] - pls[second]["cy"])
                    < (first_h + second_h) / 2.0 + gap - 1e-9
                ):
                    return True
        return False

    if rotation == "auto":
        for _ in range(4):
            current = _actual_bbox(placements)
            current_area = (current[2] - current[0]) * (current[3] - current[1])
            improved = False
            for placement in placements:
                if abs(placement["actual_w"] - placement["actual_h"]) < 1e-9:
                    continue
                old_rotation = placement.get("rotation", 0)
                placement["rotation"] = 90 if old_rotation == 0 else 0
                candidate = _actual_bbox(placements)
                candidate_area = (candidate[2] - candidate[0]) * (
                    candidate[3] - candidate[1]
                )
                if (
                    candidate_area < current_area - 1e-9
                    and not _actual_overlap(placements, compact_min_gap)
                ):
                    current_area = candidate_area
                    improved = True
                else:
                    placement["rotation"] = old_rotation
            if not improved:
                break

    top = gdspy.Cell(array_name)
    top_lib.add(top)
    for placement in placements:
        info = cell_info[placement["label"]]
        cx = placement["cx"]
        cy = placement["cy"]
        placed_w, placed_h = _placed_dims(placement)
        cell = (
            info["cell_rot90"]
            if placement.get("rotation", 0) == 90
            else info["cell"]
        )
        top.add(
            gdspy.CellReference(
                cell,
                origin=(cx - placed_w / 2.0, cy - placed_h / 2.0),
            )
        )
        top.add(gdspy.Label(placement["text"], (cx, cy), layer=100))

    actual_left, actual_bottom, actual_right, actual_top = _actual_bbox(placements)
    width = actual_right - actual_left
    height = actual_top - actual_bottom

    res_left = min(p["cx"] - p["w"] / 2.0 for p in placements)
    res_right = max(p["cx"] + p["w"] / 2.0 for p in placements)
    res_bottom = min(p["cy"] - p["h"] / 2.0 for p in placements)
    res_top = max(p["cy"] + p["h"] / 2.0 for p in placements)
    reserved_width = res_right - res_left
    reserved_height = res_top - res_bottom

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    top_lib.write_gds(str(output_path))
    return {
        "path": str(output_path),
        "rows": rows,
        "cols": cols,
        "width": width,
        "height": height,
        "reserved_width": reserved_width,
        "reserved_height": reserved_height,
        "band_width": total_w,
        "band_height": total_h,
        "square_reservation": bool(square_reservation),
        "rotation": rotation,
        "placements": placements,
        "dummies": int(np.count_nonzero(arr == 0)),
        "compact": compact_stats,
    }
