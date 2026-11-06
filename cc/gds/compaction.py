"""Centroid-preserving placement compaction.

The table writer lays unit cells out on symmetric row/column bands.  Because
unit cells have different sizes, those bands leave whitespace.  This module
moves cell *centers* -- never cell sizes -- toward the common matching center.

The key property is that a matching constraint is preserved when a set of
cell centers is scaled about its own centroid:

    c = mean(p_i)
    p_i' = c + s * (p_i - c)
    mean(p_i') = c

The same holds for an axial pair when ``c`` is the pair midpoint.  All matched
groups in the table layout share one physical center, so the common case is a
single ``(sx, sy)`` scale about that center.  If the groups do not share a
center, each constraint set is scaled about its own center instead.

Only positions are compressed; cell width/height and orientation are
unchanged.  The search rejects any scale that creates an overlap (with an
optional minimum clearance).
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Sequence

import numpy as np


def _normalize_groups(groups):
    """Accept ``[(name, (a, b)), ...]`` or ``[(a, b), ...]``."""
    if groups is None:
        return None
    normalized = []
    for item in groups:
        if (
            isinstance(item, (tuple, list))
            and len(item) == 2
            and isinstance(item[1], (tuple, list))
        ):
            names = tuple(item[1])
        elif isinstance(item, (tuple, list)) and all(
            isinstance(value, str) for value in item
        ):
            names = tuple(item)
        else:
            raise ValueError(f"unsupported group specification: {item!r}")
        if names:
            normalized.append(names)
    return normalized


def _real_cells_by_name(placements):
    by_name: dict[str, list[int]] = {}
    for index, placement in enumerate(placements):
        if placement.get("dummy"):
            continue
        by_name.setdefault(placement["name"], []).append(index)
    return by_name


def _constraint_sets(placements, groups):
    """Return cell-index sets whose centroid/midpoint must be preserved.

    * A device with two or more cells is one common-centroid set.
    * A group whose devices each have exactly one cell is one axial-pair set.
    * ``groups=None`` falls back to per-device sets.
    """
    groups = _normalize_groups(groups)
    by_name = _real_cells_by_name(placements)
    sets: list[list[int]] = []

    if groups is None:
        for indices in by_name.values():
            if len(indices) >= 2:
                sets.append(list(indices))
    else:
        for names in groups:
            if any(name not in by_name for name in names):
                continue
            if all(len(by_name[name]) == 1 for name in names) and len(names) >= 2:
                sets.append([index for name in names for index in by_name[name]])
            else:
                for name in names:
                    indices = by_name[name]
                    if indices:
                        # A single-cell device is kept as a one-element set:
                        # it cannot move, but its centroid still participates
                        # in the verification.
                        sets.append(list(indices))

    unique_sets = []
    seen = set()
    for indices in sets:
        key = frozenset(indices)
        if key not in seen:
            seen.add(key)
            unique_sets.append(indices)
    return unique_sets


def _center(placements, indices):
    points = np.array(
        [(placements[index]["cx"], placements[index]["cy"]) for index in indices],
        dtype=float,
    )
    return points.mean(axis=0)


def _common_center(placements, sets, tolerance=1e-6):
    if not sets:
        return None
    centers = np.array([_center(placements, indices) for indices in sets])
    center = centers.mean(axis=0)
    if np.all(np.linalg.norm(centers - center, axis=1) <= tolerance):
        return center
    return None


def _dummy_center(placements):
    real = [
        (placement["cx"], placement["cy"])
        for placement in placements
        if not placement.get("dummy")
    ]
    if not real:
        return np.zeros(2, dtype=float)
    return np.array(real, dtype=float).mean(axis=0)


def _scaled_placements(
    placements,
    sets,
    common_center,
    dummy_center,
    sx,
    sy,
):
    center_of: dict[int, np.ndarray] = {}
    for indices in sets:
        center = (
            common_center
            if common_center is not None
            else _center(placements, indices)
        )
        for index in indices:
            center_of[index] = center

    scaled = []
    for index, placement in enumerate(placements):
        updated = dict(placement)
        center = center_of.get(index, dummy_center)
        updated["cx"] = float(center[0] + sx * (placement["cx"] - center[0]))
        updated["cy"] = float(center[1] + sy * (placement["cy"] - center[1]))
        scaled.append(updated)
    return scaled


def _has_overlap(placements, min_gap=0.0):
    count = len(placements)
    for first in range(count):
        a = placements[first]
        for second in range(first + 1, count):
            b = placements[second]
            if (
                abs(a["cx"] - b["cx"])
                < (a["w"] + b["w"]) / 2.0 + min_gap - 1e-9
                and abs(a["cy"] - b["cy"])
                < (a["h"] + b["h"]) / 2.0 + min_gap - 1e-9
            ):
                return True
    return False


def _bbox(placements):
    left = min(p["cx"] - p["w"] / 2.0 for p in placements)
    right = max(p["cx"] + p["w"] / 2.0 for p in placements)
    bottom = min(p["cy"] - p["h"] / 2.0 for p in placements)
    top = max(p["cy"] + p["h"] / 2.0 for p in placements)
    return left, bottom, right, top


def _area(placements):
    left, bottom, right, top = _bbox(placements)
    return (right - left) * (top - bottom), right - left, top - bottom


def _centroid_error(placements, sets, initial_centers):
    if not sets:
        return 0.0
    errors = []
    for indices, initial in zip(sets, initial_centers):
        errors.append(np.linalg.norm(_center(placements, indices) - initial))
    return float(max(errors))


@dataclass
class CompactionStats:
    applied: bool
    sx: float
    sy: float
    initial_area: float
    final_area: float
    initial_width: float
    initial_height: float
    final_width: float
    final_height: float
    centroid_error: float
    constraint_count: int
    common_center: tuple[float, float] | None

    def as_dict(self):
        return {
            "applied": self.applied,
            "sx": self.sx,
            "sy": self.sy,
            "initial_area": self.initial_area,
            "final_area": self.final_area,
            "initial_width": self.initial_width,
            "initial_height": self.initial_height,
            "final_width": self.final_width,
            "final_height": self.final_height,
            "centroid_error": self.centroid_error,
            "constraint_count": self.constraint_count,
            "common_center": self.common_center,
        }


def compact_placements(
    placements: Sequence[dict],
    groups=None,
    *,
    min_gap: float = 0.0,
    sy_values: Sequence[float] | None = None,
    max_scale: float = 1.0,
    min_scale: float = 0.05,
    center_tolerance: float = 1e-6,
    binary_iterations: int = 40,
) -> tuple[list[dict], dict]:
    """Compress placements toward their matching centers.

    For each candidate ``sy`` the smallest valid ``sx`` is found by binary
    search.  Validity means "no overlap" and, by construction, all matching
    centroids/midpoints remain fixed.  ``max_scale=1.0`` means compaction
    only; values above 1 would stretch an axis.
    """
    placements = [dict(placement) for placement in placements]
    sets = _constraint_sets(placements, groups)
    if not sets:
        area, width, height = _area(placements)
        return placements, CompactionStats(
            applied=False,
            sx=1.0,
            sy=1.0,
            initial_area=area,
            final_area=area,
            initial_width=width,
            initial_height=height,
            final_width=width,
            final_height=height,
            centroid_error=0.0,
            constraint_count=0,
            common_center=None,
        ).as_dict()

    common = _common_center(placements, sets, center_tolerance)
    dummy_center = _dummy_center(placements)
    initial_centers = [_center(placements, indices) for indices in sets]
    initial_area, initial_width, initial_height = _area(placements)

    def build(sx, sy):
        return _scaled_placements(
            placements, sets, common, dummy_center, sx, sy
        )

    if sy_values is None:
        sy_values = np.linspace(1.0, 0.7, 31)

    best = (initial_area, 1.0, 1.0, placements)
    for sy in sy_values:
        sy = float(sy)
        if max_scale < sy or sy < min_scale:
            continue
        if _has_overlap(build(max_scale, sy), min_gap):
            continue

        low, high = min_scale, max_scale
        for _ in range(binary_iterations):
            middle = (low + high) / 2.0
            if _has_overlap(build(middle, sy), min_gap):
                low = middle
            else:
                high = middle
        sx = high

        candidate = build(sx, sy)
        area, _, _ = _area(candidate)
        if area < best[0] - 1e-9:
            best = (area, sx, sy, candidate)

    area, sx, sy, final = best
    final_width = _area(final)[1]
    final_height = _area(final)[2]
    centroid_error = _centroid_error(final, sets, initial_centers)

    stats = CompactionStats(
        applied=True,
        sx=float(sx),
        sy=float(sy),
        initial_area=float(initial_area),
        final_area=float(area),
        initial_width=float(initial_width),
        initial_height=float(initial_height),
        final_width=float(final_width),
        final_height=float(final_height),
        centroid_error=centroid_error,
        constraint_count=len(sets),
        common_center=None if common is None else tuple(map(float, common)),
    )
    return final, stats.as_dict()


def pair_move_compact(
    placements: Sequence[dict],
    groups=None,
    *,
    min_gap: float = 0.3,
    iterations: int = 200_000,
    seed: int = 0,
    step_start: float = 2.0,
    step_end: float = 0.15,
    temp_start: float = 20.0,
    temp_end: float = 0.05,
) -> tuple[list[dict], dict]:
    """Compact with cell-pair moves that keep every centroid fixed.

    A move picks two cells from one constraint set (one device for a
    common-centroid group, or the two cells of an axial pair) and applies
    ``+d`` to one cell and ``-d`` to the other.  The set centroid/midpoint is
    therefore unchanged exactly.

    Unlike :func:`compact_placements`, this is a local search: it can move a
    cell across row/column bands and fill a local gap while its partner moves
    the other way.  Moves that violate ``min_gap`` are rejected, so the
    original minimum clearance is preserved.  Dummy cells are unconstrained
    and may be moved individually to fill remaining holes.

    The objective is the bounding-box area; simulated annealing allows
    temporary area increases to escape a local minimum.
    """
    placements = [dict(placement) for placement in placements]
    sets = _constraint_sets(placements, groups)
    dummy_indices = [
        index for index, placement in enumerate(placements) if placement["dummy"]
    ]
    if not sets and not dummy_indices:
        area, width, height = _area(placements)
        return placements, {
            "applied": False,
            "mode": "pair_moves",
            "initial_area": area,
            "final_area": area,
            "initial_width": width,
            "initial_height": height,
            "final_width": width,
            "final_height": height,
            "centroid_error": 0.0,
            "constraint_count": 0,
            "min_gap": min_gap,
            "iterations": iterations,
            "seed": seed,
            "max_displacement": 0.0,
            "mean_displacement": 0.0,
        }

    def _collision_dims(placement):
        return (
            float(placement.get("collision_w", placement["w"])),
            float(placement.get("collision_h", placement["h"])),
        )

    def _collision_area(pls):
        left = min(p["cx"] - _collision_dims(p)[0] / 2.0 for p in pls)
        right = max(p["cx"] + _collision_dims(p)[0] / 2.0 for p in pls)
        bottom = min(p["cy"] - _collision_dims(p)[1] / 2.0 for p in pls)
        top = max(p["cy"] + _collision_dims(p)[1] / 2.0 for p in pls)
        return (right - left) * (top - bottom), right - left, top - bottom

    cells = [
        {
            "x": float(placement["cx"]),
            "y": float(placement["cy"]),
            "w": _collision_dims(placement)[0],
            "h": _collision_dims(placement)[1],
        }
        for placement in placements
    ]
    pair_sets = [indices for indices in sets if len(indices) >= 2]
    initial_centers = [_center(placements, indices) for indices in sets]
    initial_area, initial_width, initial_height = _collision_area(placements)
    best_positions = [(cell["x"], cell["y"]) for cell in cells]
    best_area = initial_area

    def current_area():
        left = min(cell["x"] - cell["w"] / 2.0 for cell in cells)
        right = max(cell["x"] + cell["w"] / 2.0 for cell in cells)
        bottom = min(cell["y"] - cell["h"] / 2.0 for cell in cells)
        top = max(cell["y"] + cell["h"] / 2.0 for cell in cells)
        return (right - left) * (top - bottom)

    def overlaps(first, second):
        a = cells[first]
        b = cells[second]
        return (
            abs(a["x"] - b["x"])
            < (a["w"] + b["w"]) / 2.0 + min_gap - 1e-9
            and abs(a["y"] - b["y"])
            < (a["h"] + b["h"]) / 2.0 + min_gap - 1e-9
        )

    def valid_moved(indices):
        moved = set(indices)
        for first in indices:
            for second in range(len(cells)):
                if second in moved:
                    continue
                if overlaps(first, second):
                    return False
        if len(indices) == 2 and overlaps(indices[0], indices[1]):
            return False
        return True

    rng = random.Random(seed)
    for iteration in range(iterations):
        fraction = iteration / max(iterations, 1)
        temperature = temp_start * (temp_end / temp_start) ** fraction
        step = step_start * (step_end / step_start) ** fraction

        use_pair = bool(pair_sets) and rng.random() < 0.9
        if use_pair:
            indices = rng.choice(pair_sets)
            first, second = rng.sample(indices, 2)
            dx = rng.uniform(-step, step)
            dy = rng.uniform(-step, step)
            if abs(dx) < 1e-12 and abs(dy) < 1e-12:
                continue
            old_area = current_area()
            old_first = (cells[first]["x"], cells[first]["y"])
            old_second = (cells[second]["x"], cells[second]["y"])
            cells[first]["x"] += dx
            cells[first]["y"] += dy
            cells[second]["x"] -= dx
            cells[second]["y"] -= dy
            if valid_moved([first, second]):
                new_area = current_area()
                delta = new_area - old_area
                if delta < 0.0 or rng.random() < math.exp(-delta / temperature):
                    if new_area < best_area:
                        best_area = new_area
                        best_positions = [(cell["x"], cell["y"]) for cell in cells]
                    continue
            cells[first]["x"], cells[first]["y"] = old_first
            cells[second]["x"], cells[second]["y"] = old_second
        elif dummy_indices:
            index = rng.choice(dummy_indices)
            dx = rng.uniform(-step, step)
            dy = rng.uniform(-step, step)
            if abs(dx) < 1e-12 and abs(dy) < 1e-12:
                continue
            old_area = current_area()
            old = (cells[index]["x"], cells[index]["y"])
            cells[index]["x"] += dx
            cells[index]["y"] += dy
            if valid_moved([index]):
                new_area = current_area()
                delta = new_area - old_area
                if delta < 0.0 or rng.random() < math.exp(-delta / temperature):
                    if new_area < best_area:
                        best_area = new_area
                        best_positions = [(cell["x"], cell["y"]) for cell in cells]
                    continue
            cells[index]["x"], cells[index]["y"] = old

    for cell, (x, y) in zip(cells, best_positions):
        cell["x"], cell["y"] = x, y

    final = []
    for placement, cell in zip(placements, cells):
        updated = dict(placement)
        updated["cx"] = cell["x"]
        updated["cy"] = cell["y"]
        final.append(updated)

    final_area, final_width, final_height = _collision_area(final)
    centroid_error = max(
        (np.linalg.norm(_center(final, indices) - initial) for indices, initial in zip(sets, initial_centers)),
        default=0.0,
    )
    displacements = [
        math.hypot(
            cell["x"] - placement["cx"], cell["y"] - placement["cy"]
        )
        for cell, placement in zip(cells, placements)
    ]
    stats = {
        "applied": final_area < initial_area - 1e-9,
        "mode": "pair_moves",
        "initial_area": float(initial_area),
        "final_area": float(final_area),
        "initial_width": float(initial_width),
        "initial_height": float(initial_height),
        "final_width": float(final_width),
        "final_height": float(final_height),
        "centroid_error": float(centroid_error),
        "constraint_count": len(sets),
        "min_gap": float(min_gap),
        "iterations": int(iterations),
        "seed": int(seed),
        "max_displacement": float(max(displacements, default=0.0)),
        "mean_displacement": float(
            sum(displacements) / len(displacements) if displacements else 0.0
        ),
    }
    return final, stats
