"""Block-level floorplanning for matched-group GDS blocks.

Two packers are provided:

``shelf_pack`` / ``best_shelf_pack``
    Simple deterministic row packing.  Kept for comparison and debugging.

``sa_skyline_pack``
    Skyline bottom-left packing with optional 90-degree block rotation,
    driven by simulated annealing.  This is the default for the ota3 flow.

The ``Block`` interface is independent of GDS; the caller provides the
original cell and, optionally, a pre-rotated cell.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import random
from typing import Iterable


@dataclass
class Block:
    name: str
    width: float
    height: float
    cell: object | None = None
    cell_rot90: object | None = None
    x: float = 0.0
    y: float = 0.0
    rotation: int = 0

    @property
    def area(self) -> float:
        return self.width * self.height

    @property
    def eff_width(self) -> float:
        return self.height if self.rotation else self.width

    @property
    def eff_height(self) -> float:
        return self.width if self.rotation else self.height


@dataclass
class Floorplan:
    blocks: list[Block]
    width: float
    height: float
    margin: float = 0.0

    @property
    def area(self) -> float:
        return self.width * self.height


def shelf_pack(
    blocks: Iterable[Block],
    max_width: float | None = None,
    spacing: float = 2.0,
    target_aspect: float = 1.4,
) -> Floorplan:
    """Pack fixed-size blocks into rows (deterministic shelf packing)."""
    blocks = list(blocks)
    if not blocks:
        return Floorplan([], 0.0, 0.0)

    total_area = sum(block.area for block in blocks)
    widest = max(block.width for block in blocks)
    if max_width is None:
        max_width = max(widest, math.sqrt(total_area) * target_aspect)

    ordered = sorted(blocks, key=lambda b: (b.height, b.width), reverse=True)

    cursor_x = 0.0
    row_bottom = 0.0
    row_height = 0.0

    for block in ordered:
        if cursor_x > 0.0 and cursor_x + block.width > max_width:
            cursor_x = 0.0
            # Rows grow upward; the next row starts above the previous one.
            row_bottom += row_height + spacing
            row_height = 0.0

        block.x = cursor_x
        block.y = row_bottom
        cursor_x += block.width + spacing
        row_height = max(row_height, block.height)

    max_x = max(block.x + block.width for block in ordered)
    max_y = max(block.y + block.height for block in ordered)

    return Floorplan(ordered, max_x, max_y)


def best_shelf_pack(
    blocks: Iterable[Block],
    spacing: float = 2.0,
    target_aspect: float = 1.2,
    samples: int = 30,
) -> Floorplan:
    """Try several row widths and keep the most compact shelf plan."""
    blocks = list(blocks)
    if not blocks:
        return Floorplan([], 0.0, 0.0)

    total_area = sum(block.area for block in blocks)
    widest = max(block.width for block in blocks)
    lo = widest
    hi = max(widest, math.sqrt(total_area) * 3.0)

    best_width = lo
    best_score = float("inf")
    for index in range(samples):
        max_width = lo + (hi - lo) * index / max(samples - 1, 1)
        plan = shelf_pack(blocks, max_width=max_width, spacing=spacing)
        aspect = plan.width / plan.height if plan.height else 1.0
        penalty = 1.0 + 0.25 * abs(math.log(max(aspect, 1e-9) / target_aspect))
        score = plan.area * penalty
        if score < best_score:
            best_score = score
            best_width = max_width

    return shelf_pack(blocks, max_width=best_width, spacing=spacing)


def _skyline_pack(
    items: list[tuple[float, float]],
    spacing: float = 0.0,
) -> tuple[list[tuple[float, float]], float, float]:
    """Bottom-left skyline packing for ``(width, height)`` items.

    Returns ``(positions, bbox_width, bbox_height)``.  The skyline reserves an
    extra ``spacing`` to the right and above each block; the returned bbox is
    based on the real block extents.
    """
    if not items:
        return [], 0.0, 0.0

    # Each segment is [x0, x1, top_y].
    skyline: list[list[float]] = [[0.0, float("inf"), 0.0]]
    positions: list[tuple[float, float]] = []

    for width, height in items:
        candidates = {0.0}
        for x0, x1, _ in skyline:
            if x0 < 1e11:
                candidates.add(x0)

        best_key = None
        best_x = 0.0
        best_y = 0.0
        for x in sorted(candidates):
            y = 0.0
            for x0, x1, top in skyline:
                if x1 > x and x0 < x + width:
                    y = max(y, top)
            key = (y, x)
            if best_key is None or key < best_key:
                best_key = key
                best_x = x
                best_y = y

        positions.append((best_x, best_y))

        # Reserve [x, x + width + spacing] at top y + height + spacing.
        left = best_x
        right = best_x + width + spacing
        new_top = best_y + height + spacing

        updated: list[list[float]] = []
        for x0, x1, top in skyline:
            if x1 <= left or x0 >= right:
                updated.append([x0, x1, top])
            else:
                if x0 < left:
                    updated.append([x0, left, top])
                if x1 > right:
                    updated.append([right, x1, top])
        updated.append([left, right, new_top])
        updated.sort(key=lambda segment: segment[0])

        merged: list[list[float]] = []
        for segment in updated:
            if (
                merged
                and abs(merged[-1][2] - segment[2]) < 1e-9
                and abs(merged[-1][1] - segment[0]) < 1e-9
            ):
                merged[-1][1] = segment[1]
            else:
                merged.append(segment)
        skyline = merged

    bbox_width = max(x + w for (x, _), (w, _) in zip(positions, items))
    bbox_height = max(y + h for (_, y), (_, h) in zip(positions, items))
    return positions, bbox_width, bbox_height


def _pack_state(
    blocks: list[Block],
    order: list[int],
    rotation: list[int],
    spacing: float,
) -> tuple[list[tuple[float, float]], float, float]:
    """Pack one (order, rotation) state and return consistent geometry.

    The bbox is recomputed from the returned positions rather than trusted
    from an earlier call.  This keeps the SA cost and the final placement in
    sync.
    """
    items: list[tuple[float, float]] = []
    # ``rotation`` is aligned with the *placement positions*, not with the
    # block indices, so use enumerate(order) and rotation[position].
    for position, idx in enumerate(order):
        block = blocks[idx]
        rotated = rotation[position]
        items.append(
            (block.height, block.width) if rotated else (block.width, block.height)
        )

    positions, _, _ = _skyline_pack(items, spacing=spacing)
    width = max(x + w for (x, _), (w, _) in zip(positions, items))
    height = max(y + h for (_, y), (_, h) in zip(positions, items))
    return positions, width, height


def axis_aligned_column(
    blocks: Iterable[Block],
    spacing: float = 3.0,
) -> Floorplan:
    """Stack blocks so every block center lies on one vertical axis.

    This cannot put all block centroids at the same 2-D point (the blocks
    would overlap), but it gives every matching block the same vertical
    symmetry axis.  It is the practical "global alignment" layout for
    heterogeneous groups that cannot share one unit-cell array.
    """
    blocks = list(blocks)
    if not blocks:
        return Floorplan([], 0.0, 0.0)

    width = max(block.width for block in blocks)
    cursor_y = 0.0

    for block in blocks:
        # The aligned layout uses the original orientation; a rotation chosen
        # by the compact floorplanner must not leak into this layout.
        block.rotation = 0
        block.x = (width - block.width) / 2.0
        block.y = cursor_y
        cursor_y += block.height + spacing

    height = cursor_y - spacing
    return Floorplan(blocks, width, height)


def _shelf_pack_state(
    blocks: list[Block],
    order: list[int],
    rotation: list[int],
    max_width: float,
    spacing: float,
) -> tuple[list[tuple[float, float]], float, float]:
    """Pack one (order, rotation, max_width) shelf state."""
    positions: list[tuple[float, float]] = []
    dimensions: list[tuple[float, float]] = []

    cursor_x = 0.0
    row_bottom = 0.0
    row_height = 0.0

    for position, idx in enumerate(order):
        block = blocks[idx]
        rotated = rotation[position]
        width, height = (block.height, block.width) if rotated else (block.width, block.height)
        dimensions.append((width, height))

        if cursor_x > 0.0 and cursor_x + width > max_width:
            cursor_x = 0.0
            # Rows grow upward; the next row starts above the previous one.
            row_bottom += row_height + spacing
            row_height = 0.0

        positions.append((cursor_x, row_bottom))
        cursor_x += width + spacing
        row_height = max(row_height, height)

    bbox_width = max(x + w for (x, _), (w, _) in zip(positions, dimensions))
    bbox_height = max(y + h for (_, y), (_, h) in zip(positions, dimensions))
    return positions, bbox_width, bbox_height


def sa_shelf_pack(
    blocks: Iterable[Block],
    spacing: float = 2.0,
    target_aspect: float = 1.3,
    iterations: int = 4000,
    seed: int = 0,
    width_samples: int = 20,
) -> Floorplan:
    """Shelf packing with 90-degree rotation and order search by SA.

    Unlike the unconstrained skyline packer, this keeps a finite row width,
    so it produces a balanced 2-D layout instead of one very long row.
    """
    blocks = list(blocks)
    if not blocks:
        return Floorplan([], 0.0, 0.0)

    count = len(blocks)
    rng = random.Random(seed)

    total_area = sum(block.area for block in blocks)
    widest = max(max(block.width, block.height) for block in blocks)
    width_hi = max(widest, math.sqrt(total_area) * 2.5)
    width_candidates = [
        widest + (width_hi - widest) * index / max(width_samples - 1, 1)
        for index in range(width_samples)
    ]

    def evaluate(order: list[int], rotation: list[int]):
        best_score = float("inf")
        best_positions: list[tuple[float, float]] = []
        best_width = 0.0
        best_height = 0.0

        for max_width in width_candidates:
            positions, width, height = _shelf_pack_state(
                blocks, order, rotation, max_width, spacing
            )
            aspect = width / height if height else 1.0
            penalty = 1.0 + 0.35 * abs(math.log(max(aspect, 1e-9) / target_aspect))
            score = width * height * penalty
            if score < best_score:
                best_score = score
                best_positions = positions
                best_width = width
                best_height = height

        return best_score, best_width, best_height, best_positions

    order = sorted(range(count), key=lambda idx: blocks[idx].area, reverse=True)
    rotation = [0] * count
    current_cost, _, _, _ = evaluate(order, rotation)
    best = (current_cost, list(order), list(rotation))

    t_start = max(current_cost * 0.05, 1e-9)
    t_end = max(current_cost * 1e-4, 1e-12)
    temperature = t_start
    cooling = (t_end / t_start) ** (1.0 / max(iterations - 1, 1))

    for _ in range(iterations):
        new_order = list(order)
        new_rotation = list(rotation)
        roll = rng.random()

        if roll < 0.4 and count >= 2:
            i, j = rng.sample(range(count), 2)
            new_order[i], new_order[j] = new_order[j], new_order[i]
        elif roll < 0.7 and count >= 2:
            i, j = rng.sample(range(count), 2)
            new_order.insert(j, new_order.pop(i))
        else:
            i = rng.randrange(count)
            new_rotation[i] ^= 1

        cost, _, _, _ = evaluate(new_order, new_rotation)
        delta = cost - current_cost
        if delta <= 0.0 or rng.random() < math.exp(-delta / max(temperature, 1e-12)):
            order = new_order
            rotation = new_rotation
            current_cost = cost
            if cost < best[0]:
                best = (cost, list(order), list(rotation))
        temperature *= cooling

    _, best_order, best_rotation = best
    _, width, height, positions = evaluate(best_order, best_rotation)

    placed: list[Block] = []
    for idx, (x, y), rot in zip(best_order, positions, best_rotation):
        block = blocks[idx]
        block.x = x
        block.y = y
        block.rotation = rot
        placed.append(block)

    return Floorplan(placed, width, height)


def sa_skyline_pack(
    blocks: Iterable[Block],
    spacing: float = 2.0,
    target_aspect: float = 1.2,
    iterations: int = 6000,
    seed: int = 0,
) -> Floorplan:
    """Skyline packing with 90-degree rotation, optimized by simulated annealing.

    The SA state is only ``(order, rotation)``.  The final positions and
    bounding box are recomputed from the winning state, so the returned
    ``Floorplan`` is always self-consistent.
    """
    blocks = list(blocks)
    if not blocks:
        return Floorplan([], 0.0, 0.0)

    count = len(blocks)
    rng = random.Random(seed)

    def evaluate(order: list[int], rotation: list[int]) -> tuple[float, float, float, list]:
        positions, width, height = _pack_state(blocks, order, rotation, spacing)
        aspect = width / height if height else 1.0
        penalty = 1.0 + 0.2 * abs(math.log(max(aspect, 1e-9) / target_aspect))
        return width * height * penalty, width, height, positions

    order = sorted(range(count), key=lambda idx: blocks[idx].area, reverse=True)
    rotation = [0] * count
    current_cost, _, _, _ = evaluate(order, rotation)
    best = (current_cost, list(order), list(rotation))

    t_start = max(current_cost * 0.05, 1e-9)
    t_end = max(current_cost * 1e-4, 1e-12)
    temperature = t_start
    cooling = (t_end / t_start) ** (1.0 / max(iterations - 1, 1))

    for _ in range(iterations):
        new_order = list(order)
        new_rotation = list(rotation)
        roll = rng.random()

        if roll < 0.4 and count >= 2:
            i, j = rng.sample(range(count), 2)
            new_order[i], new_order[j] = new_order[j], new_order[i]
        elif roll < 0.7 and count >= 2:
            i, j = rng.sample(range(count), 2)
            new_order.insert(j, new_order.pop(i))
        else:
            i = rng.randrange(count)
            new_rotation[i] ^= 1

        cost, _, _, _ = evaluate(new_order, new_rotation)
        delta = cost - current_cost
        if delta <= 0.0 or rng.random() < math.exp(-delta / max(temperature, 1e-12)):
            order = new_order
            rotation = new_rotation
            current_cost = cost
            if cost < best[0]:
                best = (cost, list(order), list(rotation))
        temperature *= cooling

    _, best_order, best_rotation = best
    positions, width, height = _pack_state(blocks, best_order, best_rotation, spacing)

    placed: list[Block] = []
    for idx, (x, y), rot in zip(best_order, positions, best_rotation):
        block = blocks[idx]
        block.x = x
        block.y = y
        block.rotation = rot
        placed.append(block)

    return Floorplan(placed, width, height)
