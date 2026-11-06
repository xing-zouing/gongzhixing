#!/usr/bin/env python3
"""Generic matched-device layout flow.

The cases are configured in ``CASES`` below.  For each matching pair in the
``.sym`` file:

* if both devices decompose into one unit, they are merged into one axial-pair
  label (count 2);
* otherwise each device gets its own label with its unit count.

Unmatched single-device groups are ignored.  The labels are fed to the
project common-centroid algorithm, the resulting matrix is written as a
heterogeneous GDS table, and centroid-preserving pair moves compact the
placement.  Outputs are:

    <case>_matched.gds          layout
    <case>_matched.png          GDS render
    <case>_matched_groups.png   matching-group / dummy annotation

Run:

    .venv/bin/python scripts/run_matched_layout.py --case ota3
    .venv/bin/python scripts/run_matched_layout.py --case all
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Rectangle as MplRectangle  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cc.eva import EvA  # noqa: E402
from cc.gds.device_cells import device_cell  # noqa: E402
from cc.gds.render import render_gds  # noqa: E402
from cc.gds.table_writer import generate_mixed_table_gds  # noqa: E402
from cc.io.netlist_parser import load_case  # noqa: E402
from cc.main import cc_main_flow  # noqa: E402


# ---------------------------------------------------------------------------
# Case configuration: netlist, matching constraints, and optional unit splits.
# ---------------------------------------------------------------------------
CASES = {
    "ota3": {
        "netlist": ROOT / "examples" / "ota3.sp",
        "sym": ROOT / "examples" / "ota3.sym",
        "subckt": "OTA2_my",
        # ota3 uses 8 units of 5 fingers for M0/M7 and 2 units of 3 fingers
        # for M11/M10.  The other matched devices are single-unit.
        "units": {"M0": 8, "M7": 8, "M11": 2, "M10": 2},
    },
    "ota1": {
        "netlist": ROOT / "mycase" / "ota1.sp",
        "sym": ROOT / "mycase" / "ota1.sym",
        "subckt": "OPA_my",
        "units": {},
    },
    "ota2": {
        "netlist": ROOT / "mycase" / "ota2.sp",
        "sym": ROOT / "mycase" / "ota2.sym",
        "subckt": "ota2",
        "units": {},
    },
}


def _as_matrix(result) -> np.ndarray:
    matrix = np.asarray(result, dtype=int)
    if matrix.ndim == 1:
        side = int(np.ceil(np.sqrt(matrix.size)))
        matrix = np.pad(matrix, (0, side * side - matrix.size)).reshape(side, side)
    return matrix


def _mos_spec(dev, nf: int) -> dict:
    return {
        "kind": "mos",
        "name": dev.name,
        "nch": dev.nch,
        "w": dev.unit_width,
        "l": dev.unit_length,
        "nf": nf,
        "attr": dev.attr,
    }


def _cap_spec(dev) -> dict:
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


def _total_fingers(dev) -> int:
    if dev.kind == "mos":
        return dev.multi * dev.nf
    return 1


def _unit_aspect(dev, nf: int) -> float:
    _, width, height = device_cell(dev, nf_override=nf)
    return max(width, height) / min(width, height)


def _choose_units(
    dev,
    name: str,
    overrides: dict,
    max_units: int = 8,
    aspect_threshold: float = 3.5,
) -> int:
    """Pick a unit count, splitting further when the unit cell is too elongated.

    The initial split is the explicit override, or the largest divisor of the
    total finger count up to ``max_units``.  If its unit-cell aspect exceeds
    ``aspect_threshold``, try splits with more units (fewer fingers per unit),
    in increasing unit count, and take the first one that passes.  If none
    passes, keep the initial split.
    """
    if dev.kind != "mos":
        return 1
    total = _total_fingers(dev)
    if total <= 1:
        return 1

    override = overrides.get(name)
    if override is not None:
        initial = int(override)
        if initial < 1 or total % initial != 0:
            raise ValueError(
                f"{name}: unit count {initial} does not divide {total} fingers"
            )
    else:
        divisors = [value for value in range(1, total + 1) if total % value == 0 and value <= max_units]
        initial = max(divisors) if divisors else 1

    nf_initial = total // initial
    aspect_initial = _unit_aspect(dev, nf_initial)
    if aspect_threshold is None or aspect_threshold <= 0:
        return initial
    if aspect_initial <= aspect_threshold:
        return initial

    candidates = sorted(
        value
        for value in range(initial + 1, total + 1)
        if total % value == 0 and value <= max_units
    )
    for units in candidates:
        nf = total // units
        aspect = _unit_aspect(dev, nf)
        if aspect <= aspect_threshold:
            print(
                f"  auto-split {name}: {initial} units of nf={nf_initial} "
                f"(aspect {aspect_initial:.2f}) -> {units} units of nf={nf} "
                f"(aspect {aspect:.2f})"
            )
            return units

    print(
        f"  keep {name}: {initial} units of nf={nf_initial} "
        f"(aspect {aspect_initial:.2f}; no split passes)"
    )
    return initial


def _collect_records(devices, groups, overrides, max_units=8, aspect_threshold=3.5):
    records = []
    for group in groups:
        names = [name for name in group if name in devices]
        if len(names) != len(group):
            print(f"  skip {group}: missing device(s)")
            continue
        if len(names) == 1:
            # Unmatched single-device group; deliberately ignored.
            continue
        if len(names) != 2:
            print(f"  skip {group}: only two-device groups are supported")
            continue
        if any(devices[name].kind not in ("mos", "cap") for name in names):
            print(f"  skip {group}: unsupported device kind")
            continue
        units = [
            _choose_units(
                devices[name],
                name,
                overrides,
                max_units=max_units,
                aspect_threshold=aspect_threshold,
            )
            for name in names
        ]
        records.append(
            {
                "names": names,
                "units": units,
                "group_name": f"pair_{names[0]}_{names[1]}",
            }
        )
    return records


def _build_labels(records):
    """Separate-label groups first, then merged axial pairs."""
    separate = [
        record
        for record in records
        if not (record["units"][0] == 1 and record["units"][1] == 1)
    ]
    merged = [
        record
        for record in records
        if record["units"][0] == 1 and record["units"][1] == 1
    ]

    counts = []
    label_devices = []
    groups_for_compaction = []
    device_order = []
    for record in separate + merged:
        names = record["names"]
        units = record["units"]
        groups_for_compaction.append((record["group_name"], tuple(names)))
        for name in names:
            if name not in device_order:
                device_order.append(name)
        if record in merged:
            counts.append(2)
            label_devices.append(list(names))
        else:
            for name, unit_count in zip(names, units):
                counts.append(unit_count)
                label_devices.append([name])
    return counts, label_devices, groups_for_compaction, device_order


def _run_case(case_name: str, config: dict, out_root: Path, args) -> dict:
    out_dir = out_root / case_name
    out_dir.mkdir(parents=True, exist_ok=True)

    devices, groups = load_case(
        config["netlist"], config["sym"], subckt=config.get("subckt")
    )
    records = _collect_records(
        devices,
        groups,
        config.get("units", {}),
        max_units=args.max_units,
        aspect_threshold=args.unit_aspect_threshold,
    )
    if not records:
        raise ValueError(f"{case_name}: no supported matching pairs")

    counts, label_devices, groups_for_compaction, device_order = _build_labels(records)
    print(f"\n=== {case_name} ===")
    print("matching groups:", [name for name, _ in groups_for_compaction])
    print("label counts:", counts)
    print("device order:", device_order)

    best = None
    for seed in range(3):
        result = cc_main_flow(
            counts,
            square_array=True,
            row_numbers=0,
            random_state=seed,
            max_iter=args.max_iter,
        )
        matrix = _as_matrix(result)
        tmoc = EvA(matrix, matrix.shape)
        if best is None or tmoc < best[0]:
            best = (tmoc, matrix, seed)
        if tmoc == 0.0:
            break
    tmoc, label_matrix, seed = best
    rows, cols = label_matrix.shape
    print("label matrix:")
    print(label_matrix)
    print(f"TMOC: {tmoc} (GA seed {seed})")
    if tmoc != 0.0:
        print("WARNING: algorithm did not reach TMOC=0 for these counts")

    id_of = {name: index + 1 for index, name in enumerate(device_order)}
    device_matrix = np.zeros_like(label_matrix)
    for label, names in enumerate(label_devices, start=1):
        cells = sorted(tuple(pos) for pos in np.argwhere(label_matrix == label))
        expected = counts[label - 1]
        if len(cells) != expected:
            raise ValueError(
                f"label {label}: matrix has {len(cells)} cells, expected {expected}"
            )
        if len(names) == 1:
            for row, col in cells:
                device_matrix[row, col] = id_of[names[0]]
        else:
            if len(cells) != len(names):
                raise ValueError(
                    f"label {label}: {len(cells)} cells for {len(names)} devices"
                )
            for (row, col), name in zip(cells, names):
                device_matrix[row, col] = id_of[name]

    gds_devices = []
    for name in device_order:
        dev = devices[name]
        units = next(
            record["units"][record["names"].index(name)]
            for record in records
            if name in record["names"]
        )
        if dev.kind == "mos":
            gds_devices.append(_mos_spec(dev, _total_fingers(dev) // units))
        else:
            gds_devices.append(_cap_spec(dev))

    gds_path = out_dir / f"{case_name}_matched.gds"
    info = generate_mixed_table_gds(
        devices=gds_devices,
        matrix=device_matrix,
        output_path=gds_path,
        array_name=f"{case_name}_matched",
        gap=0.3,
        fill_dummies=True,
        groups=groups_for_compaction,
        compact=True,
        compact_min_gap=args.min_gap,
        compact_mode="pairs",
        compact_iterations=args.iterations,
        compact_seed=args.seed,
        square_reservation=False,
        rotation="horizontal",
    )

    render_png = out_dir / f"{case_name}_matched.png"
    preview = render_gds(gds_path, render_png, dpi=180)
    annotated_png = out_dir / f"{case_name}_matched_groups.png"
    _draw_annotated(info, groups_for_compaction, case_name, annotated_png)

    positions = {}
    for placement in info["placements"]:
        if placement["dummy"]:
            continue
        positions.setdefault(placement["name"], []).append(
            (placement["cx"], placement["cy"])
        )
    effective_centers = []
    for record in records:
        device_centers = [
            np.array(positions[name], dtype=float).mean(axis=0)
            for name in record["names"]
        ]
        if all(unit == 1 for unit in record["units"]):
            effective_centers.append(np.mean(device_centers, axis=0))
        else:
            effective_centers.extend(device_centers)
    effective_centers = np.array(effective_centers)
    group_center = effective_centers.mean(axis=0)
    group_center_error = float(
        np.max(np.linalg.norm(effective_centers - group_center, axis=1))
    )

    cell_area = sum(
        placement["actual_w"] * placement["actual_h"]
        for placement in info["placements"]
    )
    actual_area = info["width"] * info["height"]
    stats = info.get("compact") or {}
    result = {
        "case": case_name,
        "subckt": config.get("subckt"),
        "groups": len(groups_for_compaction),
        "counts": counts,
        "matrix": f"{rows}x{cols}",
        "tmoc": float(tmoc),
        "width": info["width"],
        "height": info["height"],
        "area": actual_area,
        "utilization": cell_area / actual_area if actual_area else 0.0,
        "centroid_error": stats.get("centroid_error", 0.0),
        "group_center_error": group_center_error,
        "min_gap": stats.get("min_gap", args.min_gap),
        "gds": str(gds_path),
        "png": str(preview.get("png", render_png)),
        "groups_png": str(annotated_png),
    }
    print(
        f"{case_name}: bbox {result['width']:.2f} x {result['height']:.2f} "
        f"= {result['area']:.1f} um^2, utilization {result['utilization'] * 100:.1f}%, "
        f"set centroid error {result['centroid_error']:.2e} um, "
        f"group-center error {result['group_center_error']:.2e} um"
    )
    print("GDS:", result["gds"])
    print("preview:", result["png"])
    print("groups preview:", result["groups_png"])
    return result


def _draw_annotated(info: dict, groups, case_name: str, output_path: Path) -> None:
    group_of = {}
    group_names = []
    for group_name, names in groups:
        group_names.append(group_name)
        for name in names:
            group_of[name] = group_name
    palette = plt.get_cmap("tab10").colors
    color_of_group = {
        group_name: palette[index % len(palette)]
        for index, group_name in enumerate(group_names)
    }

    placements = info["placements"]
    left = min(p["cx"] - p["actual_w"] / 2.0 for p in placements)
    right = max(p["cx"] + p["actual_w"] / 2.0 for p in placements)
    bottom = min(p["cy"] - p["actual_h"] / 2.0 for p in placements)
    top = max(p["cy"] + p["actual_h"] / 2.0 for p in placements)
    width = right - left
    height = top - bottom

    fig, ax = plt.subplots(figsize=(max(7.0, width / 5.0), max(7.0, height / 5.0)))
    for placement in placements:
        face = (
            "#d9d9d9"
            if placement["dummy"]
            else color_of_group[group_of[placement["name"]]]
        )
        ax.add_patch(
            MplRectangle(
                (
                    placement["cx"] - placement["actual_w"] / 2.0,
                    placement["cy"] - placement["actual_h"] / 2.0,
                ),
                placement["actual_w"],
                placement["actual_h"],
                facecolor=face,
                edgecolor="black",
                linewidth=0.6,
                alpha=0.8,
            )
        )
        ax.text(
            placement["cx"],
            placement["cy"],
            placement["text"],
            ha="center",
            va="center",
            fontsize=5.5,
            color="black",
        )
    ax.set_xlim(left - width * 0.05, right + width * 0.05)
    ax.set_ylim(bottom - height * 0.05, top + height * 0.05)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"{case_name}: one color per matching group; D(...) = dummy")

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
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _write_summary(results: list[dict], out_root: Path, args) -> None:
    lines = [
        "# Matched layout summary",
        "",
        f"> horizontal layout; min gap {args.min_gap:g} um; "
        f"unit aspect threshold {args.unit_aspect_threshold:g}; "
        f"max units {args.max_units}",
        "",
        "| case | groups | counts | matrix | TMOC | bbox (um) | area (um^2) | util | set centroid err (um) | group-center err (um) | min gap (um) |",
        "|---|---:|---|---:|---:|---|---:|---:|---:|---:|---:|",
    ]
    for result in results:
        lines.append(
            f"| {result['case']} | {result['groups']} | {result['counts']} | "
            f"{result['matrix']} | {result['tmoc']:.3g} | "
            f"{result['width']:.2f} x {result['height']:.2f} | "
            f"{result['area']:.1f} | {result['utilization'] * 100:.1f}% | "
            f"{result['centroid_error']:.2e} | "
            f"{result['group_center_error']:.2e} | {result['min_gap']:.2f} |"
        )
    lines.append("")
    for result in results:
        lines.append(f"## {result['case']}")
        lines.append(f"- GDS: `{result['gds']}`")
        lines.append(f"- render: `{result['png']}`")
        lines.append(f"- groups: `{result['groups_png']}`")
        lines.append("")
    summary_path = out_root / "summary.md"
    summary_path.write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary:", summary_path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        choices=sorted(CASES) + ["all"],
        default="all",
    )
    parser.add_argument("--out-root", default=str(ROOT / "build" / "generic"))
    parser.add_argument("--min-gap", type=float, default=0.3)
    parser.add_argument("--iterations", type=int, default=200_000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-iter", type=int, default=200_000)
    parser.add_argument(
        "--unit-aspect-threshold",
        type=float,
        default=3.5,
        help="split a unit cell further when max(w,h)/min(w,h) exceeds this; <=0 disables",
    )
    parser.add_argument(
        "--max-units",
        type=int,
        default=8,
        help="maximum unit count per device considered by the auto-split policy",
    )
    args = parser.parse_args()

    case_names = sorted(CASES) if args.case == "all" else [args.case]
    out_root = Path(args.out_root)
    results = []
    for case_name in case_names:
        results.append(_run_case(case_name, CASES[case_name], out_root, args))
    _write_summary(results, out_root, args)


if __name__ == "__main__":
    main()
