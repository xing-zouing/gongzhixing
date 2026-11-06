#!/usr/bin/env python3
"""Step 4: interleaving comparison and the 2-D checkerboard demo.

The 1-D ABBA block is now generated automatically by the compatibility
cluster flow in Step 3, so this script only writes the 2-D checkerboard
GDS.  It still prints the aabb/abba/baab metrics for comparison.

    a a b b   -> separate / concatenated, pair midpoints differ
    a b b a   -> ABBA, both pair midpoints at the array center
    b a a b   -> BAAB
    2-D:
        a b
        b a   -> both pairs share the 2-D center

Run from the project root:

    .venv/bin/python scripts/ota3_step4_interleave.py
"""

from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cc.io.netlist_parser import load_case  # noqa: E402
from cc.place.interleave import (  # noqa: E402
    PairGroup,
    abba_1d,
    baab_1d,
    checkerboard_2d,
    pair_metrics,
    print_metrics,
)
from cc.gds.array_writer import generate_group_gds  # noqa: E402
from cc.gds.render import render_gds  # noqa: E402


def _gds_device(dev) -> dict:
    return {
        "name": dev.name,
        "nch": dev.nch,
        "model": dev.model,
        "w": dev.unit_width,
        "l": dev.unit_length,
        "nf": dev.nf,
        "multi": dev.multi,
        "attr": dev.attr,
    }


def _clean_interleave_outputs(out_dir: Path) -> None:
    """Remove stale outputs of *this* demo scripts so the directory does not grow.

    Other demos (e.g. multiunit_pair_k2) are intentionally left alone.
    """
    if not out_dir.exists():
        return
    names = ("interleave_concatenated_aabb", "interleave_abba", "interleave_baab",
             "interleave_checkerboard_2d")
    for name in names:
        for suffix in (".gds", ".png"):
            path = out_dir / f"{name}{suffix}"
            if path.exists():
                path.unlink()


def _to_ids(matrix, name_to_id):
    if matrix and isinstance(matrix[0], list):
        return [[name_to_id[name] for name in row] for row in matrix]
    # Keep an explicit 1 x N strip; do not let the GDS writer reshape it
    # into a square.
    return [[name_to_id[name] for name in matrix]]


def main() -> None:
    netlist = ROOT / "examples" / "ota3.sp"
    sym = ROOT / "examples" / "ota3.sym"
    out_dir = ROOT / "build" / "ota3"
    devices, _ = load_case(netlist, sym, subckt="OTA2_my")

    pair_a = PairGroup("A_M2_M5", ("M2", "M5"))
    pair_b = PairGroup("B_M4_M3", ("M4", "M3"))
    names = [pair_a.devices[0], pair_a.devices[1], pair_b.devices[0], pair_b.devices[1]]
    devs = [devices[name] for name in names]

    # All four unit cells must be compatible for a shared grid.
    ref = devs[0]
    for dev in devs[1:]:
        if (dev.unit_width, dev.unit_length, dev.nf, dev.nch) != (
            ref.unit_width,
            ref.unit_length,
            ref.nf,
            ref.nch,
        ):
            raise ValueError(f"unit-cell mismatch: {dev.name} vs {ref.name}")
    print("unit cells:", [(d.name, d.unit_width, d.unit_length, d.nf) for d in devs])

    name_to_id = {dev.name: i + 1 for i, dev in enumerate(devs)}
    gds_devices = [_gds_device(dev) for dev in devs]

    # Evaluate all patterns, but only write the two canonical GDS demos.
    patterns = [
        ("concatenated_aabb", ["M2", "M5", "M4", "M3"]),
        ("abba", abba_1d(pair_a, pair_b)),
        ("baab", baab_1d(pair_a, pair_b)),
        ("checkerboard_2d", checkerboard_2d(pair_a, pair_b)),
    ]
    # The 1-D ABBA GDS is produced automatically by the Step 3 cluster flow.
    write_patterns = {"checkerboard_2d"}

    _clean_interleave_outputs(out_dir)

    all_metrics = {}
    for name, matrix in patterns:
        metric_matrix = [matrix] if matrix and isinstance(matrix[0], str) else matrix
        metrics = pair_metrics(metric_matrix, [pair_a, pair_b])
        all_metrics[name] = metrics
        print_metrics(name, metrics)

        if name not in write_patterns:
            print("  (metrics only; GDS not written)")
            continue

        matrix_ids = _to_ids(matrix, name_to_id)
        output = out_dir / f"interleave_{name}.gds"
        info = generate_group_gds(
            devices=gds_devices,
            cc_matrix=matrix_ids,
            output_path=output,
            array_name=f"interleave_{name}",
        )
        print(f"  GDS: {info['path']}")
        render_gds(output, out_dir / f"interleave_{name}.png", dpi=200)

    # A compact comparison of pair-center errors.
    print("\n=== pair center-error summary ===")
    for name, metrics in all_metrics.items():
        errors = [f"{m['pair']}={m['center_error']}" for m in metrics]
        print(f"  {name:20s} " + "  ".join(errors))


if __name__ == "__main__":
    main()
