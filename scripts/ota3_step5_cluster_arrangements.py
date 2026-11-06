#!/usr/bin/env python3
"""Step 5: enumerate all valid arrangements of the compatible cluster.

For the ota3 cluster M2/M5 + M4/M3:

* call the project's common-centroid algorithm at the *group* level
  (counts [2,2]), producing ABBA/BAAB group patterns;
* also include the constructive 2-D result and row-major reshapes;
* assign each group's two unit cells to its two devices in every way;
* evaluate pair midpoint / weighted centroid / dispersion;
* write a contact sheet of all unique device-level arrangements and GDS
  files for the valid ones.

This is an inspection script: it does not change the top-level block used by
Step 3.

Run from the project root:

    .venv/bin/python scripts/ota3_step5_cluster_arrangements.py
"""

from __future__ import annotations

from itertools import combinations, product
from math import ceil
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from cc.ca_sym import construction_algorithm_symmetry  # noqa: E402
from cc.eva import EvA  # noqa: E402
from cc.gds.array_writer import generate_group_gds  # noqa: E402
from cc.gds.render import render_gds  # noqa: E402
from cc.io.netlist_parser import load_case  # noqa: E402
from cc.main import cc_main_flow  # noqa: E402
from cc.place.compatibility import cluster_groups, group_unit_counts  # noqa: E402
from cc.place.interleave import Group, group_metrics  # noqa: E402
from cc.plan.ota3_groups import GROUPS  # noqa: E402


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
    for group_id, total in enumerate(totals, start=1):
        if int(np.count_nonzero(matrix == group_id)) != total:
            return False
    return True


def _assign_group_cells(cells, devices, counts):
    """Yield every partition of ``cells`` matching the two device counts."""
    if len(devices) != 2:
        raise NotImplementedError("only two-device groups are supported")
    for chosen in combinations(cells, counts[0]):
        chosen_set = set(chosen)
        rest = [cell for cell in cells if cell not in chosen_set]
        if len(rest) != counts[1]:
            continue
        yield {devices[0]: list(chosen), devices[1]: rest}


def _build_layout(matrix: np.ndarray, groups: list[Group]) -> list[list[str | None]]:
    rows, cols = matrix.shape
    group_cells = {
        index + 1: [tuple(map(int, pos)) for pos in np.argwhere(matrix == index + 1)]
        for index in range(len(groups))
    }

    per_group_options = []
    for index, group in enumerate(groups, start=1):
        cells = group_cells[index]
        if len(cells) != sum(group.counts):
            return []
        per_group_options.append(
            list(_assign_group_cells(cells, group.devices, group.counts))
        )

    layouts = []
    for combination in product(*per_group_options):
        layout: list[list[str | None]] = [[None] * cols for _ in range(rows)]
        for assignment in combination:
            for device, cells in assignment.items():
                for row, col in cells:
                    layout[row][col] = device
        layouts.append(layout)
    return layouts


def _layout_key(layout) -> tuple:
    return tuple(tuple("" if cell is None else cell for cell in row) for row in layout)


def _evaluate(layout, groups, gds_devices):
    metrics = group_metrics(layout, groups)
    center = metrics[0]["array_center"]
    midpoint_error = max(max(item["midpoint_error"]) for item in metrics)
    weighted_error = max(max(item["weighted_error"]) for item in metrics)
    dispersion = sum(
        (centroid[0] - center[0]) ** 2 + (centroid[1] - center[1]) ** 2
        for item in metrics
        for centroid in item["centroids"]
    )
    valid = midpoint_error < 1e-9 and weighted_error < 1e-9
    return {
        "metrics": metrics,
        "midpoint_error": midpoint_error,
        "weighted_error": weighted_error,
        "dispersion": dispersion,
        "valid": valid,
        "shape": (len(layout), len(layout[0]) if layout else 0),
    }


def _layout_text(layout) -> str:
    return " ".join(
        f"[{' '.join('.' if cell is None else cell for cell in row)}]"
        for row in layout
    )


def _contact_sheet(arrangements, output: Path) -> None:
    count = len(arrangements)
    columns = 4
    rows = ceil(count / columns)
    fig, axes = plt.subplots(rows, columns, figsize=(columns * 3.2, rows * 3.0))
    axes = np.atleast_1d(axes).reshape(rows, columns)

    for index, item in enumerate(arrangements):
        ax = axes[index // columns][index % columns]
        layout = item["layout"]
        n_rows, n_cols = item["shape"]
        for row in range(n_rows):
            for col in range(n_cols):
                device = layout[row][col]
                if device is None:
                    face = "#eeeeee"
                elif device.startswith("M2") or device.startswith("M5"):
                    face = "#f4a261"
                else:
                    face = "#8ecae6"
                ax.add_patch(
                    Rectangle((col, n_rows - 1 - row), 1, 1, facecolor=face, edgecolor="black")
                )
                if device is not None:
                    ax.text(
                        col + 0.5, n_rows - 1 - row + 0.5, device,
                        ha="center", va="center", fontsize=7,
                    )
        ax.set_xlim(0, n_cols)
        ax.set_ylim(0, n_rows)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(
            f"#{index:02d} {n_rows}x{n_cols} {'valid' if item['valid'] else 'INVALID'}\n"
            f"disp={item['dispersion']:.2f}  mid_err={item['midpoint_error']:.2f}",
            fontsize=7,
        )

    for index in range(count, rows * columns):
        axes[index // columns][index % columns].axis("off")

    fig.tight_layout()
    fig.savefig(output, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    netlist = ROOT / "examples" / "ota3.sp"
    sym = ROOT / "examples" / "ota3.sym"
    out_dir = ROOT / "build" / "ota3" / "arrangements"
    out_dir.mkdir(parents=True, exist_ok=True)

    for stale in out_dir.glob("arrangement_*"):
        stale.unlink()
    for stale in out_dir.glob("cluster_best.*"):
        stale.unlink()

    devices, _ = load_case(netlist, sym, subckt="OTA2_my")
    clusters = cluster_groups(devices, GROUPS)
    cluster = next(c for c in clusters if c["interleavable"])

    group_objects: list[Group] = []
    gds_devices: list[dict] = []
    for spec in cluster["groups"]:
        counts = group_unit_counts(devices, spec)
        group_objects.append(
            Group(spec["name"], tuple(spec["devices"]), counts)
        )
        for name in spec["devices"]:
            dev = devices[name]
            gds_devices.append(
                {
                    "name": dev.name,
                    "nch": dev.nch,
                    "model": dev.model,
                    "w": dev.unit_width,
                    "l": dev.unit_length,
                    "nf": 1,
                    "multi": dev.multi,
                    "attr": dev.attr,
                }
            )

    totals = [sum(group.counts) for group in group_objects]
    print("cluster groups:", [(g.name, g.devices, g.counts) for g in group_objects])
    print("group-level totals:", totals)

    # Candidate group-level arrangements from the project algorithm.
    candidates: list[tuple[str, np.ndarray]] = []
    calls = [
        ("cc_main_flow square/auto", dict(square_array=True, row_numbers=0)),
        ("cc_main_flow square/rows2", dict(square_array=True, row_numbers=2)),
        ("cc_main_flow nonsquare/ver", dict(square_array=False, orientation="ver", row_numbers=0)),
        ("cc_main_flow nonsquare/hor", dict(square_array=False, orientation="hor", row_numbers=0)),
    ]
    for label, kwargs in calls:
        try:
            result = cc_main_flow(
                totals, random_state=0, max_iter=200000, **kwargs
            )
        except Exception as exc:  # pragma: no cover - defensive
            print(f"  {label}: skipped ({type(exc).__name__}: {exc})")
            continue
        for matrix in _matrix_candidates(result):
            if _valid_counts(matrix, totals):
                candidates.append((label, matrix))

    constructive = construction_algorithm_symmetry(totals, square_array=True)
    for matrix in _matrix_candidates(constructive):
        if _valid_counts(matrix, totals):
            candidates.append(("constructive ca_sym", matrix))

    # Deduplicate group-level matrices, then enumerate device assignments.
    unique_group_matrices: dict[tuple, tuple[set[str], np.ndarray]] = {}
    for label, matrix in candidates:
        key = tuple(map(tuple, matrix.tolist()))
        if key not in unique_group_matrices:
            unique_group_matrices[key] = (set(), matrix)
        unique_group_matrices[key][0].add(label)

    arrangements: dict[tuple, dict] = {}
    for key, (sources, matrix) in unique_group_matrices.items():
        print(f"\ngroup matrix {matrix.shape} from {sorted(sources)}:\n{matrix}")
        print(f"  group-level TMOC (with 0=dummy): {EvA(matrix.flatten(), matrix.shape)}")
        for layout in _build_layout(matrix, group_objects):
            layout_key = _layout_key(layout)
            if layout_key not in arrangements:
                arrangements[layout_key] = {
                    "layout": layout,
                    "sources": set(),
                    **_evaluate(layout, group_objects, gds_devices),
                }
            arrangements[layout_key]["sources"].update(sources)

    ordered = sorted(
        arrangements.values(),
        key=lambda item: (
            not item["valid"],
            item["midpoint_error"],
            item["weighted_error"],
            item["dispersion"],
        ),
    )

    print("\n=== all unique device-level arrangements ===")
    for index, item in enumerate(ordered):
        print(
            f"#{index:02d} {item['shape'][0]}x{item['shape'][1]} "
            f"{'valid' if item['valid'] else 'INVALID'} "
            f"mid_err={item['midpoint_error']:.3f} "
            f"weighted_err={item['weighted_error']:.3f} "
            f"disp={item['dispersion']:.3f} "
            f"src={sorted(item['sources'])}"
        )
        print(f"    {_layout_text(item['layout'])}")

    contact = out_dir / "cluster_arrangements.png"
    _contact_sheet(ordered, contact)
    print(f"\ncontact sheet: {contact}")

    # Write GDS/PNG for every valid arrangement.
    name_to_id = {device["name"]: index + 1 for index, device in enumerate(gds_devices)}
    valid_written = 0
    for index, item in enumerate(ordered):
        if not item["valid"]:
            continue
        matrix = [
            [0 if cell is None else name_to_id[cell] for cell in row]
            for row in item["layout"]
        ]
        name = f"arrangement_{index:02d}_{item['shape'][0]}x{item['shape'][1]}"
        output = out_dir / f"{name}.gds"
        generate_group_gds(
            devices=gds_devices,
            cc_matrix=matrix,
            output_path=output,
            array_name=name,
        )
        render_gds(output, out_dir / f"{name}.png", dpi=200)
        valid_written += 1
    print(f"valid GDS written: {valid_written}")

    # The first valid entry is the best by midpoint/weighted error/dispersion.
    best = next(item for item in ordered if item["valid"])
    best_matrix = [
        [0 if cell is None else name_to_id[cell] for cell in row]
        for row in best["layout"]
    ]
    best_gds = out_dir / "cluster_best.gds"
    generate_group_gds(
        devices=gds_devices,
        cc_matrix=best_matrix,
        output_path=best_gds,
        array_name="cluster_best",
    )
    render_gds(best_gds, out_dir / "cluster_best.png", dpi=200)
    print(f"\nbest valid arrangement: {_layout_text(best['layout'])}")
    print(f"best GDS: {best_gds}")


if __name__ == "__main__":
    main()
