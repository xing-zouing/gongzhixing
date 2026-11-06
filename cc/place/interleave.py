"""Interleaved multi-group matching patterns.

The block-level flow places each matching group in its own GDS block.  A
shared array instead interleaves several groups so that all of them share
the same symmetry center, e.g. two differential pairs:

    a a b b        (separate blocks / concatenation)
    a b b a        (interleaved, both pair midpoints at the array center)
    b a a b

For two pairs with one unit per device this module provides the classic
1-D ABBA/BAAB and 2-D checkerboard patterns.  The same idea generalizes to
more groups and multiple units via symmetry-orbit allocation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


@dataclass(frozen=True)
class PairGroup:
    """A two-device matching group (e.g. a differential pair)."""

    name: str
    devices: tuple[str, str]


@dataclass(frozen=True)
class Group:
    """A matching group with per-device unit counts.

    ``counts`` allows ratioed devices (e.g. a 1:2 current mirror).  The
    current allocator supports two-device groups.
    """

    name: str
    devices: tuple[str, ...]
    counts: tuple[int, ...]

    def __post_init__(self) -> None:
        if len(self.devices) != len(self.counts):
            raise ValueError(f"{self.name}: devices and counts length mismatch")
        if any(count < 0 for count in self.counts):
            raise ValueError(f"{self.name}: counts must be non-negative")
        if sum(self.counts) == 0:
            raise ValueError(f"{self.name}: group has no cells")


def _orbits_1d(length: int) -> tuple[list[tuple[tuple[int, int], tuple[int, int]]], tuple[int, int] | None]:
    """Reflection orbits of a 1-D array, plus the center cell if odd."""
    if length < 1:
        raise ValueError("length must be >= 1")
    orbits = [((0, i), (0, length - 1 - i)) for i in range(length // 2)]
    center = (0, length // 2) if length % 2 else None
    return orbits, center


def _orbits_2d(
    rows: int, cols: int
) -> tuple[list[tuple[tuple[int, int], tuple[int, int]]], tuple[int, int] | None]:
    """180-degree rotation orbits of a 2-D array, plus the center cell."""
    if rows < 1 or cols < 1:
        raise ValueError("rows/cols must be >= 1")
    seen: set[tuple[int, int]] = set()
    orbits: list[tuple[tuple[int, int], tuple[int, int]]] = []
    center_cell = None
    for row in range(rows):
        for col in range(cols):
            cell = (row, col)
            if cell in seen:
                continue
            partner = (rows - 1 - row, cols - 1 - col)
            if partner == cell:
                center_cell = cell
                seen.add(cell)
            else:
                orbits.append((cell, partner))
                seen.add(cell)
                seen.add(partner)
    return orbits, center_cell


def _allocate_groups(
    groups: list[Group],
    orbits: list[tuple[tuple[int, int], tuple[int, int]]],
    center: tuple[int, int] | None,
) -> dict[tuple[int, int], str]:
    """Allocate symmetry orbits to two-device groups.

    Each group's cells are a union of whole orbits (plus the center cell when
    its total count is odd), so the group's centroid is the array center.
    Within a group, whole orbits give device-level common centroid; a split
    orbit is used when both device counts are odd.
    """
    needs_center = 0
    total_orbits = 0
    for group in groups:
        total = sum(group.counts)
        total_orbits += total // 2
        needs_center += total % 2
        if len(group.devices) != 2:
            raise NotImplementedError("only two-device groups are supported")
    if total_orbits > len(orbits):
        raise ValueError(
            f"not enough orbits: need {total_orbits}, have {len(orbits)}"
        )
    if needs_center > (1 if center is not None else 0):
        raise ValueError("odd total count requires an available center cell")

    # Round-robin orbit allocation spreads each group across the array.
    allocation: dict[str, list[tuple[tuple[int, int], tuple[int, int]]]] = {
        group.name: [] for group in groups
    }
    remaining = {group.name: sum(group.counts) // 2 for group in groups}
    orbit_index = 0
    while any(value > 0 for value in remaining.values()):
        progressed = False
        for group in groups:
            if remaining[group.name] > 0:
                allocation[group.name].append(orbits[orbit_index])
                orbit_index += 1
                remaining[group.name] -= 1
                progressed = True
        if not progressed:
            break

    center_group = next(
        (group for group in groups if sum(group.counts) % 2 == 1), None
    )
    layout: dict[tuple[int, int], str] = {}

    for group in groups:
        counts = list(group.counts)
        assigned = allocation[group.name]
        orbit_cursor = 0
        device_cells: list[list[tuple[int, int]]] = [[], []]
        odd = [i for i, count in enumerate(counts) if count % 2 == 1]

        if not odd:
            for device_index, count in enumerate(counts):
                for _ in range(count // 2):
                    orbit = assigned[orbit_cursor]
                    orbit_cursor += 1
                    device_cells[device_index].extend(orbit)
        elif len(odd) == 1:
            if center is None or center_group is not group:
                raise ValueError(
                    f"{group.name}: odd device count needs a center cell"
                )
            odd_index = odd[0]
            device_cells[odd_index].append(center)
            for _ in range((counts[odd_index] - 1) // 2):
                orbit = assigned[orbit_cursor]
                orbit_cursor += 1
                device_cells[odd_index].extend(orbit)
            even_index = 1 - odd_index
            for _ in range(counts[even_index] // 2):
                orbit = assigned[orbit_cursor]
                orbit_cursor += 1
                device_cells[even_index].extend(orbit)
        else:
            # Both counts odd: whole orbits plus one split orbit.
            for device_index, count in enumerate(counts):
                for _ in range((count - 1) // 2):
                    orbit = assigned[orbit_cursor]
                    orbit_cursor += 1
                    device_cells[device_index].extend(orbit)
            orbit = assigned[orbit_cursor]
            orbit_cursor += 1
            device_cells[0].append(orbit[0])
            device_cells[1].append(orbit[1])

        for device_index, count in enumerate(counts):
            if len(device_cells[device_index]) != count:
                raise RuntimeError(
                    f"{group.name}: device {group.devices[device_index]} "
                    f"got {len(device_cells[device_index])} cells, expected {count}"
                )
            for cell in device_cells[device_index]:
                if cell in layout:
                    raise RuntimeError(f"cell {cell} assigned twice")
                layout[cell] = group.devices[device_index]

    return layout


def allocate_orbits_1d(groups: list[Group], length: int) -> list[str | None]:
    """Allocate groups to a 1-D shared array (reflection symmetry)."""
    orbits, center = _orbits_1d(length)
    layout_map = _allocate_groups(groups, orbits, center)
    return [layout_map.get((0, index)) for index in range(length)]


def allocate_orbits_2d(groups: list[Group], rows: int, cols: int) -> list[list[str | None]]:
    """Allocate groups to a 2-D shared array (180-degree rotation)."""
    orbits, center = _orbits_2d(rows, cols)
    layout_map = _allocate_groups(groups, orbits, center)
    return [
        [layout_map.get((row, col)) for col in range(cols)]
        for row in range(rows)
    ]


def group_metrics(matrix: list[list[str | None]], groups: list[Group]) -> list[dict]:
    """Centroid and symmetry metrics for arbitrary two-device groups."""
    rows = len(matrix)
    cols = len(matrix[0]) if rows else 0
    center = ((rows - 1) / 2.0, (cols - 1) / 2.0)
    results = []

    for group in groups:
        centroids = []
        for device in group.devices:
            positions = _device_positions(matrix, device)
            if not positions:
                raise ValueError(f"{group.name}: device {device} not found")
            centroids.append(_centroid(positions))

        total = sum(group.counts)
        weighted = (
            sum(count * centroid[0] for count, centroid in zip(group.counts, centroids)) / total,
            sum(count * centroid[1] for count, centroid in zip(group.counts, centroids)) / total,
        )
        midpoint = (
            (centroids[0][0] + centroids[1][0]) / 2.0,
            (centroids[0][1] + centroids[1][1]) / 2.0,
        )
        results.append(
            {
                "group": group.name,
                "devices": group.devices,
                "counts": group.counts,
                "centroids": tuple(centroids),
                "weighted_centroid": weighted,
                "midpoint": midpoint,
                "array_center": center,
                "weighted_error": (
                    abs(weighted[0] - center[0]),
                    abs(weighted[1] - center[1]),
                ),
                "midpoint_error": (
                    abs(midpoint[0] - center[0]),
                    abs(midpoint[1] - center[1]),
                ),
                "device_center_errors": tuple(
                    (abs(c[0] - center[0]), abs(c[1] - center[1]))
                    for c in centroids
                ),
            }
        )
    return results


def print_group_metrics(title: str, metrics: list[dict]) -> None:
    print(f"\n=== {title} ===")
    for item in metrics:
        print(
            f"  {item['group']} devices={item['devices']} counts={item['counts']}"
        )
        for device, centroid, error in zip(
            item["devices"], item["centroids"], item["device_center_errors"]
        ):
            print(f"    {device}: centroid={centroid}, center_error={error}")
        print(f"    weighted centroid={item['weighted_centroid']}, error={item['weighted_error']}")
        print(f"    pair midpoint={item['midpoint']}, error={item['midpoint_error']}")


def abba_1d(pair_a: PairGroup, pair_b: PairGroup) -> list[str]:
    """``a0 b0 b1 a1`` — both pair midpoints are the array center."""
    a0, a1 = pair_a.devices
    b0, b1 = pair_b.devices
    return [a0, b0, b1, a1]


def baab_1d(pair_a: PairGroup, pair_b: PairGroup) -> list[str]:
    """``b0 a0 a1 b1`` — the mirror variant of ABBA."""
    a0, a1 = pair_a.devices
    b0, b1 = pair_b.devices
    return [b0, a0, a1, b1]


def checkerboard_2d(pair_a: PairGroup, pair_b: PairGroup) -> list[list[str]]:
    """2-D interleaving::

        a0 b0
        b1 a1

    Both pairs have their centroid at the array center (1/2, 1/2), and the
    pattern is invariant under 180-degree rotation.
    """
    a0, a1 = pair_a.devices
    b0, b1 = pair_b.devices
    return [[a0, b0], [b1, a1]]


def interleave_pairs_multiunit_1d(
    pair_a: PairGroup,
    pair_b: PairGroup,
    units_per_device: int,
) -> list[str]:
    """Interleave two pairs when every device has ``units_per_device`` units.

    The array length is ``4 * units_per_device``.  Reflection orbits are
    ``{i, L-1-i}``.  Orbits are assigned round-robin to the two pairs; within
    a pair, whole orbits are assigned alternately to the two devices.

    * even ``units_per_device``: every device gets whole orbits, so each
      device's centroid is exactly the array center (true device-level
      common centroid).
    * odd ``units_per_device``: one orbit is split between the two devices,
      so the pair midpoint is the center, but an individual device centroid
      is not necessarily the center (pair-level mirror symmetry).
    """
    k = int(units_per_device)
    if k < 1:
        raise ValueError("units_per_device must be >= 1")

    length = 4 * k
    layout: list[str | None] = [None] * length
    pair_counts = {pair_a.name: 0, pair_b.name: 0}

    for orbit_index in range(length // 2):
        left = orbit_index
        right = length - 1 - orbit_index
        pair = pair_a if orbit_index % 2 == 0 else pair_b
        pair_orbit = pair_counts[pair.name]
        pair_counts[pair.name] += 1

        if k % 2 == 0 or pair_orbit < k - 1:
            # Whole orbit to one device.
            device = pair.devices[pair_orbit % 2]
            layout[left] = device
            layout[right] = device
        else:
            # Odd k: split the last orbit between the two devices.
            layout[left] = pair.devices[0]
            layout[right] = pair.devices[1]

    if any(name is None for name in layout):
        raise RuntimeError("internal error: unfilled interleave cell")
    return [name for name in layout if name is not None]


def _device_positions(matrix: list[list[str]], device: str) -> list[tuple[float, float]]:
    return [
        (float(row), float(col))
        for row, line in enumerate(matrix)
        for col, name in enumerate(line)
        if name == device
    ]


def _centroid(positions: list[tuple[float, float]]) -> tuple[float, float]:
    if not positions:
        raise ValueError("device not found in matrix")
    return (
        sum(p[0] for p in positions) / len(positions),
        sum(p[1] for p in positions) / len(positions),
    )


def pair_metrics(
    matrix: list[list[str]],
    pairs: Iterable[PairGroup],
) -> list[dict]:
    """Report centroid/midpoint/symmetry-error for each pair.

    For a pair of devices, the ideal shared-array layout has the two device
    centroids mirror-symmetric about the array center, i.e. their midpoint is
    the center and their offsets cancel.
    """
    rows = len(matrix)
    cols = len(matrix[0]) if rows else 0
    center = ((rows - 1) / 2.0, (cols - 1) / 2.0)

    results = []
    for pair in pairs:
        c0 = _centroid(_device_positions(matrix, pair.devices[0]))
        c1 = _centroid(_device_positions(matrix, pair.devices[1]))
        midpoint = ((c0[0] + c1[0]) / 2.0, (c0[1] + c1[1]) / 2.0)
        center_error = (
            abs(midpoint[0] - center[0]),
            abs(midpoint[1] - center[1]),
        )
        # (c0 - center) + (c1 - center); zero means mirror symmetry.
        symmetry_error = (
            abs(c0[0] + c1[0] - 2.0 * center[0]),
            abs(c0[1] + c1[1] - 2.0 * center[1]),
        )
        results.append(
            {
                "pair": pair.name,
                "devices": pair.devices,
                "centroids": (c0, c1),
                "midpoint": midpoint,
                "array_center": center,
                "center_error": center_error,
                "symmetry_error": symmetry_error,
            }
        )
    return results


def print_metrics(title: str, metrics: list[dict]) -> None:
    print(f"\n=== {title} ===")
    for item in metrics:
        c0, c1 = item["centroids"]
        print(f"  {item['pair']} ({item['devices'][0]}, {item['devices'][1]})")
        print(f"    centroids     : {c0} , {c1}")
        print(f"    midpoint      : {item['midpoint']}")
        print(f"    array center  : {item['array_center']}")
        print(f"    center error  : {item['center_error']}")
        print(f"    symmetry error: {item['symmetry_error']}")
