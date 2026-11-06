#!/usr/bin/env python3
"""Step 4b: multi-unit pair interleaving and device-level common centroid.

For two pairs (A: M2/M5, B: M4/M3) with identical unit cells:

* k=1 unit/device:  ABBA / BAAB (pair-level mirror symmetry)
* k=2 units/device: every device gets whole mirror orbits, so every device
  centroid is the array center (device-level common centroid)
* k=3 units/device: one orbit is split, so only the pair midpoint is exactly
  the center (pair-level)

This script prints the layouts and metrics for k=1,2,3 and writes the k=2
GDS/PNG.

Run from the project root:

    .venv/bin/python scripts/ota3_step4b_multiunit.py
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
    interleave_pairs_multiunit_1d,
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


def main() -> None:
    netlist = ROOT / "examples" / "ota3.sp"
    sym = ROOT / "examples" / "ota3.sym"
    out_dir = ROOT / "build" / "ota3"
    devices, _ = load_case(netlist, sym, subckt="OTA2_my")

    pair_a = PairGroup("A_M2_M5", ("M2", "M5"))
    pair_b = PairGroup("B_M4_M3", ("M4", "M3"))
    names = [pair_a.devices[0], pair_a.devices[1], pair_b.devices[0], pair_b.devices[1]]
    devs = [devices[name] for name in names]

    ref = devs[0]
    for dev in devs[1:]:
        if (dev.unit_width, dev.unit_length, dev.nf, dev.nch) != (
            ref.unit_width,
            ref.unit_length,
            ref.nf,
            ref.nch,
        ):
            raise ValueError(f"unit-cell mismatch: {dev.name} vs {ref.name}")

    name_to_id = {dev.name: i + 1 for i, dev in enumerate(devs)}
    gds_devices = [_gds_device(dev) for dev in devs]

    print(f"unit cells: w={ref.unit_width}, l={ref.unit_length}, nf={ref.nf}, nch={ref.nch}")

    for k in (1, 2, 3):
        layout = interleave_pairs_multiunit_1d(pair_a, pair_b, k)
        print(f"\n=== k={k} unit(s)/device, 1x{len(layout)} array ===")
        print("  layout:", " ".join(layout))
        metrics = pair_metrics([layout], [pair_a, pair_b])
        print_metrics(f"k={k}", metrics)
        if k % 2 == 0:
            print("  -> device-level common centroid (each device centroid = array center)")
        else:
            print("  -> pair-level mirror symmetry (pair midpoint = array center)")

    # Canonical GDS: k=2, where all four devices share the center.
    k = 2
    layout = interleave_pairs_multiunit_1d(pair_a, pair_b, k)
    matrix_ids = [[name_to_id[name] for name in layout]]

    for suffix in (".gds", ".png"):
        stale = out_dir / f"multiunit_pair_k2{suffix}"
        if stale.exists():
            stale.unlink()

    output = out_dir / "multiunit_pair_k2.gds"
    info = generate_group_gds(
        devices=gds_devices,
        cc_matrix=matrix_ids,
        output_path=output,
        array_name="multiunit_pair_k2",
    )
    print(f"\nGDS: {info['path']}")
    render_gds(output, out_dir / "multiunit_pair_k2.png", dpi=200)
    print(f"PNG: {out_dir / 'multiunit_pair_k2.png'}")


if __name__ == "__main__":
    main()
