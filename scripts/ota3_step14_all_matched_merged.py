#!/usr/bin/env python3
"""Step 14: all matched groups in ONE algorithm array (merged axial pairs).

Input to the project algorithm:

    labels:
      1 M0      8 units (5 fingers/unit)
      2 M7      8 units
      3 M11     2 units (3 fingers/unit)
      4 M10     2 units
      5 M8+M6   merged axial pair, 2 cells
      6 M2+M5   merged axial pair, 2 cells
      7 M4+M3   merged axial pair, 2 cells
      8 C2+C1   merged capacitor pair, 2 cells

    [8, 8, 2, 2, 2, 2, 2, 2] -> 6x6, TMOC = 0

After the algorithm places the merged labels, each merged label's two cells
are split between the two physical devices.  The logical group centroid is
the array center, so the physical pair is an axial pair (pair midpoint at
the center).

M9/M12 (unmatched) are excluded.

Run from the project root:

    .venv/bin/python scripts/ota3_step14_all_matched_merged.py
"""

from __future__ import annotations

import argparse
import numpy as np

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle as MplRectangle  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cc.eva import EvA  # noqa: E402
from cc.gds.render import render_gds  # noqa: E402
from cc.gds.table_writer import generate_mixed_table_gds  # noqa: E402
from cc.io.netlist_parser import load_case  # noqa: E402
from cc.main import cc_main_flow  # noqa: E402


LABEL_DEVICES = {
    1: ["M0"],
    2: ["M7"],
    3: ["M11"],
    4: ["M10"],
    5: ["M8", "M6"],
    6: ["M2", "M5"],
    7: ["M4", "M3"],
    8: ["C2", "C1"],
}

GROUPS = [
    ("input_pair_M0_M7", ("M0", "M7")),
    ("mirror_load_M11_M10", ("M11", "M10")),
    ("pair_M8_M6", ("M8", "M6")),
    ("pair_M2_M5", ("M2", "M5")),
    ("pair_M4_M3", ("M4", "M3")),
    ("cap_pair_C2_C1", ("C2", "C1")),
]


def _as_matrix(result) -> np.ndarray:
    matrix = np.asarray(result, dtype=int)
    if matrix.ndim == 1:
        side = int(np.ceil(np.sqrt(matrix.size)))
        matrix = np.pad(matrix, (0, side * side - matrix.size)).reshape(side, side)
    return matrix


def _mos(dev, nf):
    return {
        "kind": "mos",
        "name": dev.name,
        "nch": dev.nch,
        "w": dev.unit_width,
        "l": dev.unit_length,
        "nf": nf,
        "attr": dev.attr,
    }


def _cap(dev):
    return {
        "kind": "cap",
        "name": dev.name,
        "w": dev.unit_width,
        "l": dev.unit_length,
        "sp": float(dev.params.get("s", 0.07)),
        "nf": dev.nr,
        "m_bot": int(dev.params.get("stm", 3)),
        "m_top": int(dev.params.get("spm", 5)),
        "attr": ["2t"] if "2t" in dev.model.lower() else [],
        "f_tip": float(dev.params.get("ftip", 0.14)),
    }


def _centroid(matrix, device_id):
    positions = np.argwhere(matrix == device_id)
    return positions.mean(axis=0) if len(positions) else None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m0-units", type=int, choices=(4, 8), default=8)
    parser.add_argument("--name", default="ota3_all_matched")
    parser.add_argument(
        "--no-compact",
        action="store_true",
        help="keep the original row/column band positions",
    )
    parser.add_argument(
        "--min-gap",
        type=float,
        default=0.3,
        help="minimum cell clearance during compaction",
    )
    parser.add_argument(
        "--compact-mode",
        choices=("pairs", "scale", "none"),
        default="pairs",
        help="pairs = centroid-preserving cell-pair SA; scale = global scale",
    )
    parser.add_argument("--compact-iterations", type=int, default=200_000)
    parser.add_argument("--compact-seed", type=int, default=0)
    parser.add_argument(
        "--square-reservation",
        action="store_true",
        help="reserve a max(w,h) x max(w,h) square for every unit cell",
    )
    parser.add_argument(
        "--rotation",
        choices=("horizontal", "vertical", "auto"),
        default="horizontal",
        help="unit-cell orientation; auto toggles cells to shrink the actual bbox",
    )
    args = parser.parse_args()

    out_dir = ROOT / "build" / "ota3" / "all_matched"
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in out_dir.glob(f"{args.name}.*"):
        stale.unlink()

    devices, _ = load_case(
        ROOT / "examples" / "ota3.sp",
        ROOT / "examples" / "ota3.sym",
        subckt="OTA2_my",
    )

    m0_total_fingers = devices["M0"].nf * devices["M0"].multi
    if m0_total_fingers % args.m0_units != 0:
        raise ValueError("M0 unit count must divide the total finger count")
    m0_nf = m0_total_fingers // args.m0_units
    counts = [args.m0_units, args.m0_units, 2, 2, 2, 2, 2, 2]
    print(f"M0/M7: {args.m0_units} units of {m0_nf} fingers")
    result = cc_main_flow(
        counts,
        square_array=True,
        row_numbers=0,
        random_state=0,
        max_iter=200000,
    )
    label_matrix = _as_matrix(result)
    print("input:", counts)
    print("label matrix:")
    print(label_matrix)
    print("TMOC:", EvA(label_matrix, label_matrix.shape))

    # physical device id order
    device_order = [
        "M0", "M7", "M11", "M10",
        "M8", "M6", "M2", "M5", "M4", "M3",
        "C2", "C1",
    ]
    id_of = {name: index + 1 for index, name in enumerate(device_order)}
    device_matrix = np.zeros_like(label_matrix)

    for label, names in LABEL_DEVICES.items():
        cells = sorted(tuple(pos) for pos in np.argwhere(label_matrix == label))
        if len(names) == 1:
            for row, col in cells:
                device_matrix[row, col] = id_of[names[0]]
        else:
            if len(cells) != len(names):
                raise ValueError(
                    f"label {label} has {len(cells)} cells for {len(names)} devices"
                )
            for (row, col), name in zip(cells, names):
                device_matrix[row, col] = id_of[name]

    print("\nphysical device matrix:")
    print(device_matrix)

    # Build mixed MOS/cap device dicts.
    gds_devices = []
    nf_by_device = {
        "M0": m0_nf, "M7": m0_nf,
        "M11": 3, "M10": 3,
        "M8": 1, "M6": 1,
        "M2": 1, "M5": 1,
        "M4": 1, "M3": 1,
    }
    for name in device_order:
        dev = devices[name]
        if dev.kind == "cap":
            gds_devices.append(_cap(dev))
        else:
            gds_devices.append(_mos(dev, nf_by_device[name]))

    output = out_dir / f"{args.name}.gds"
    compact_mode = "none" if args.no_compact else args.compact_mode
    info = generate_mixed_table_gds(
        devices=gds_devices,
        matrix=device_matrix,
        output_path=output,
        array_name=args.name,
        gap=0.3,
        fill_dummies=True,
        groups=GROUPS,
        compact=compact_mode != "none",
        compact_min_gap=args.min_gap,
        compact_mode=compact_mode if compact_mode != "none" else "scale",
        compact_iterations=args.compact_iterations,
        compact_seed=args.compact_seed,
        square_reservation=args.square_reservation,
        rotation=args.rotation,
    )
    print("\nGDS:", {k: v for k, v in info.items() if k != "placements"})
    if info.get("compact"):
        stats = info["compact"]
        actual_area = info["width"] * info["height"]
        if stats.get("mode") == "pair_moves":
            print(
                "pair-move compaction (actual cell footprint): "
                f"area {stats['initial_area']:.1f} -> {stats['final_area']:.1f} um^2 "
                f"({(1.0 - stats['final_area'] / stats['initial_area']) * 100:.1f}% smaller), "
                f"min_gap={stats['min_gap']:.2f}, "
                f"centroid error {stats['centroid_error']:.3e} um, "
                f"max/mean cell displacement "
                f"{stats['max_displacement']:.2f}/{stats['mean_displacement']:.2f} um"
            )
        else:
            print(
                "scale compaction: "
                f"sx={stats['sx']:.4f} sy={stats['sy']:.4f} "
                f"area {stats['initial_area']:.1f} -> {stats['final_area']:.1f} um^2 "
                f"({(1.0 - stats['final_area'] / stats['initial_area']) * 100:.1f}% smaller), "
                f"centroid error {stats['centroid_error']:.3e} um"
            )
        print(
            f"actual GDS bbox: {info['width']:.2f} x {info['height']:.2f} "
            f"= {actual_area:.1f} um^2"
        )
        if info.get("square_reservation"):
            reserved_area = info["reserved_width"] * info["reserved_height"]
            print(
                f"reserved square bbox: {info['reserved_width']:.2f} x "
                f"{info['reserved_height']:.2f} = {reserved_area:.1f} um^2"
            )
    render_gds(output, out_dir / f"{args.name}.png", dpi=180)

    # Annotated schematic: one color per matching group; dummies gray and
    # explicitly labeled as D(<device>).
    group_of = {}
    group_names = []
    for group_name, names in GROUPS:
        group_names.append(group_name)
        for name in names:
            group_of[name] = group_name
    palette = plt.get_cmap("tab10").colors
    color_of_group = {
        group_name: palette[index % len(palette)]
        for index, group_name in enumerate(group_names)
    }

    def _placed_dims(placement):
        if placement.get("rotation", 0) == 90:
            return placement["actual_h"], placement["actual_w"]
        return placement["actual_w"], placement["actual_h"]

    placements_for_plot = info["placements"]
    xs = [p["cx"] for p in placements_for_plot]
    ys = [p["cy"] for p in placements_for_plot]
    plot_w = [
        max(_placed_dims(p)[0], p["w"]) for p in placements_for_plot
    ]
    plot_h = [
        max(_placed_dims(p)[1], p["h"]) for p in placements_for_plot
    ]
    left = min(p["cx"] - w / 2.0 for p, w in zip(placements_for_plot, plot_w))
    right = max(p["cx"] + w / 2.0 for p, w in zip(placements_for_plot, plot_w))
    bottom = min(p["cy"] - h / 2.0 for p, h in zip(placements_for_plot, plot_h))
    top = max(p["cy"] + h / 2.0 for p, h in zip(placements_for_plot, plot_h))
    width = right - left
    height = top - bottom

    fig, ax = plt.subplots(figsize=(max(7.0, width / 5.0), max(7.0, height / 5.0)))
    for p in placements_for_plot:
        if p["dummy"]:
            face = "#d9d9d9"
        else:
            face = color_of_group[group_of[p["name"]]]
        if info.get("square_reservation"):
            ax.add_patch(
                MplRectangle(
                    (p["cx"] - p["w"] / 2, p["cy"] - p["h"] / 2),
                    p["w"],
                    p["h"],
                    facecolor="none",
                    edgecolor="#777777",
                    linewidth=0.5,
                    linestyle="--",
                    alpha=0.7,
                )
            )
        placed_w, placed_h = _placed_dims(p)
        ax.add_patch(
            MplRectangle(
                (p["cx"] - placed_w / 2, p["cy"] - placed_h / 2),
                placed_w,
                placed_h,
                facecolor=face,
                edgecolor="black",
                linewidth=0.6,
                alpha=0.8,
            )
        )
        ax.text(
            p["cx"], p["cy"], p["text"],
            ha="center", va="center",
            fontsize=5.5,
            color="black",
        )
    ax.set_xlim(left - width * 0.05, right + width * 0.05)
    ax.set_ylim(bottom - height * 0.05, top + height * 0.05)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(
        "one color per matching group; D(...) = dummy; "
        f"rotation={info.get('rotation')}, "
        f"square_reservation={info.get('square_reservation')}"
    )

    legend_handles = [
        MplRectangle((0, 0), 1, 1, facecolor=color_of_group[name],
                     edgecolor="black", alpha=0.8)
        for name in group_names
    ] + [MplRectangle((0, 0), 1, 1, facecolor="#d9d9d9",
                      edgecolor="black", alpha=0.8)]
    ax.legend(
        legend_handles,
        group_names + ["dummy"],
        loc="upper right",
        fontsize=6,
    )
    annotated = out_dir / f"{args.name}_groups.png"
    fig.savefig(annotated, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("annotated schematic:", annotated)

    # Physical centroids from the actual placed cell centers.
    positions = {}
    for placement in info["placements"]:
        if placement["dummy"]:
            continue
        label = placement["label"]
        positions.setdefault(label, []).append((placement["cx"], placement["cy"]))

    def physical_centroid(device_name):
        label = id_of[device_name]
        points = positions.get(label, [])
        return np.mean(points, axis=0) if points else None

    all_points = [
        (placement["cx"], placement["cy"])
        for placement in info["placements"]
    ]
    center = np.mean(np.array(all_points), axis=0)
    print(f"\ntable center (mean of all slots): {center}")
    for group_name, names in GROUPS:
        c0 = physical_centroid(names[0])
        c1 = physical_centroid(names[1])
        midpoint = (c0 + c1) / 2.0
        print(
            f"  {group_name:28s} "
            f"{names[0]}={c0}, {names[1]}={c1}, "
            f"midpoint={midpoint}, delta={c0 - c1}"
        )


if __name__ == "__main__":
    main()
