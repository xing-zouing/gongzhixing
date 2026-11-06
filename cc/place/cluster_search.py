"""Search for the best arrangement of a compatible matching cluster.

For small clusters we can afford to enumerate:

* group-level arrangements produced by the project's common-centroid
  algorithm (``cc_main_flow``) and the constructive algorithm;
* every partition of each group's cells among its devices.

The valid arrangements are ranked by:

1. pair midpoint error,
2. weighted-centroid error,
3. dispersion (sum of squared device-centroid distances from the array
   center).

For the ota3 cluster (two pairs, one unit per device) this chooses a 2x2
checkerboard over the 1x4 ABBA strip because both are perfectly symmetric
but the checkerboard has lower dispersion.

This module is used by ``cc.gds.cluster_writer`` so Step 3 automatically
writes the best cluster block.
"""

from __future__ import annotations

from itertools import combinations, product
import contextlib
import io

import numpy as np

from cc.ca_sym import construction_algorithm_symmetry
from cc.main import cc_main_flow
from cc.place.compatibility import build_cluster_layout, group_unit_counts
from cc.place.interleave import Group, group_metrics


def _matrix_candidates(value) -> list[np.ndarray]:
    arr = np.asarray(value, dtype=int)
    candidates: list[np.ndarray] = []
    if arr.ndim == 1:
        candidates.append(arr.reshape(1, -1))
        side = int(round(arr.size ** 0.5))
        if side * side == arr.size:
            candidates.append(arr.reshape(side, side))
    elif arr.ndim == 2:
        candidates.append(arr)
        candidates.append(arr.flatten().reshape(1, -1))
    return candidates


def _valid_counts(matrix: np.ndarray, totals: list[int]) -> bool:
    return all(
        int(np.count_nonzero(matrix == group_id)) == total
        for group_id, total in enumerate(totals, start=1)
    )


def _assign_group_cells(cells, devices, counts):
    if len(devices) != 2:
        raise NotImplementedError("only two-device groups are supported")
    for chosen in combinations(cells, counts[0]):
        chosen_set = set(chosen)
        rest = [cell for cell in cells if cell not in chosen_set]
        if len(rest) == counts[1]:
            yield {devices[0]: list(chosen), devices[1]: rest}


def _build_layouts(matrix: np.ndarray, groups: list[Group]):
    rows, cols = matrix.shape
    group_cells = {
        index + 1: [tuple(map(int, pos)) for pos in np.argwhere(matrix == index + 1)]
        for index in range(len(groups))
    }

    options = []
    for index, group in enumerate(groups, start=1):
        cells = group_cells[index]
        if len(cells) != sum(group.counts):
            return []
        options.append(list(_assign_group_cells(cells, group.devices, group.counts)))

    layouts = []
    for combination in product(*options):
        layout: list[list[str | None]] = [[None] * cols for _ in range(rows)]
        for assignment in combination:
            for device, cells in assignment.items():
                for row, col in cells:
                    layout[row][col] = device
        layouts.append(layout)
    return layouts


def _evaluate(layout, groups: list[Group]) -> dict:
    metrics = group_metrics(layout, groups)
    center = metrics[0]["array_center"]
    midpoint_error = max(max(item["midpoint_error"]) for item in metrics)
    weighted_error = max(max(item["weighted_error"]) for item in metrics)
    dispersion = sum(
        (centroid[0] - center[0]) ** 2 + (centroid[1] - center[1]) ** 2
        for item in metrics
        for centroid in item["centroids"]
    )
    return {
        "layout": layout,
        "metrics": metrics,
        "midpoint_error": midpoint_error,
        "weighted_error": weighted_error,
        "dispersion": dispersion,
        "valid": midpoint_error < 1e-9 and weighted_error < 1e-9,
        "shape": (len(layout), len(layout[0]) if layout else 0),
    }


def _cluster_devices(devices: dict, cluster: dict) -> tuple[list[Group], list[dict]]:
    groups: list[Group] = []
    gds_devices: list[dict] = []
    seen: set[str] = set()

    for spec in cluster["groups"]:
        counts = group_unit_counts(devices, spec)
        groups.append(Group(spec["name"], tuple(spec["devices"]), counts))
        seg_nf = spec.get("seg_nf")
        for name in spec["devices"]:
            if name in seen:
                continue
            seen.add(name)
            dev = devices[name]
            gds_devices.append(
                {
                    "name": dev.name,
                    "nch": dev.nch,
                    "model": dev.model,
                    "w": dev.unit_width,
                    "l": dev.unit_length,
                    "nf": int(seg_nf) if seg_nf is not None else dev.nf,
                    "multi": dev.multi,
                    "attr": dev.attr,
                }
            )
    return groups, gds_devices


def search_arrangements(devices: dict, cluster: dict) -> list[dict]:
    """Return all unique valid/invalid device-level arrangements, ranked."""
    groups, _ = _cluster_devices(devices, cluster)
    totals = [sum(group.counts) for group in groups]

    candidates: list[np.ndarray] = []
    calls = [
        dict(square_array=True, row_numbers=0),
        dict(square_array=True, row_numbers=2),
        dict(square_array=False, orientation="ver", row_numbers=0),
        dict(square_array=False, orientation="hor", row_numbers=0),
    ]
    for kwargs in calls:
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                result = cc_main_flow(
                    totals, random_state=0, max_iter=200000, **kwargs
                )
        except Exception:
            continue
        candidates.extend(_matrix_candidates(result))

    try:
        constructive = construction_algorithm_symmetry(totals, square_array=True)
        candidates.extend(_matrix_candidates(constructive))
    except Exception:
        pass

    unique_matrices: dict[tuple, np.ndarray] = {}
    for matrix in candidates:
        if _valid_counts(matrix, totals):
            unique_matrices[tuple(map(tuple, matrix.tolist()))] = matrix

    arrangements: dict[tuple, dict] = {}
    for matrix in unique_matrices.values():
        for layout in _build_layouts(matrix, groups):
            key = tuple(
                tuple("" if cell is None else cell for cell in row)
                for row in layout
            )
            if key not in arrangements:
                arrangements[key] = _evaluate(layout, groups)

    return sorted(
        arrangements.values(),
        key=lambda item: (
            not item["valid"],
            item["midpoint_error"],
            item["weighted_error"],
            item["dispersion"],
        ),
    )


def best_arrangement(devices: dict, cluster: dict) -> dict | None:
    for arrangement in search_arrangements(devices, cluster):
        if arrangement["valid"]:
            return arrangement
    return None


def build_cluster_layout_best(devices: dict, cluster: dict) -> dict | None:
    """Best valid arrangement, with the orbit allocator as a fallback."""
    best = best_arrangement(devices, cluster)
    if best is None:
        return build_cluster_layout(devices, cluster)

    groups, gds_devices = _cluster_devices(devices, cluster)
    name_to_id = {device["name"]: index + 1 for index, device in enumerate(gds_devices)}
    matrix = [
        [0 if cell is None else name_to_id[cell] for cell in row]
        for row in best["layout"]
    ]
    return {
        "gds_devices": gds_devices,
        "matrix": matrix,
        "layout": best["layout"],
        "groups": groups,
        "length": len(matrix) * (len(matrix[0]) if matrix else 0),
        "arrangement": best,
    }
