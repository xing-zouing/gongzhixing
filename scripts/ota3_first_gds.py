#!/usr/bin/env python3
"""First end-to-end slice: ota3 netlist + .sym -> common-centroid -> GDS.

This uses the multi-level unit convention:

    unit cell = one ``multi`` instance, containing ``nf`` fingers
    unit count per device = multi

For the ota3 differential pair M0/M7 (w=2u, l=1.5u, nf=10, multi=4) this
gives ``input_list = [4, 4]`` and an 8-cell matching group plus dummies.

Run from the project root:

    .venv/bin/python scripts/ota3_first_gds.py
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cc.io.netlist_parser import load_case  # noqa: E402
from cc.main import cc_main_flow  # noqa: E402
from cc.eva import EvA  # noqa: E402
from cc.gds.array_writer import generate_group_gds  # noqa: E402


def main() -> None:
    netlist = ROOT / "examples" / "ota3.sp"
    sym = ROOT / "examples" / "ota3.sym"
    # Legacy tutorial output; kept separate so it does not mix with the
    # canonical per-group block outputs.
    output = ROOT / "build" / "ota3" / "legacy" / "M0_M7_cc_array.gds"

    devices, groups = load_case(netlist, sym, subckt="OTA2_my")

    print("=== parsed devices ===")
    for name, dev in devices.items():
        print(
            f"  {name:4s} {dev.kind:3s} {dev.model:8s} "
            f"w={dev.unit_width:g} l={dev.unit_length:g} "
            f"multi={dev.multi} nf={dev.nf} nr={dev.nr}"
        )

    print("\n=== matching groups ===")
    for idx, group in enumerate(groups, start=1):
        print(f"  G{idx}: {group}")

    # Pick the ota3 differential pair.
    group = next(g for g in groups if "M0" in g and "M7" in g)
    group_devices = [devices[name] for name in group]
    print(f"\nselected group: {group}")

    # Multi-level unit convention: one multi instance per unit cell.
    input_list = [dev.multi for dev in group_devices]
    print(f"input_list (unit counts): {input_list}")

    cc_result = cc_main_flow(input_list, square_array=True, row_numbers=0)
    cc_array = np.asarray(cc_result, dtype=int)
    if cc_array.ndim == 1:
        side = int(np.ceil(np.sqrt(cc_array.size)))
        cc_array = np.pad(cc_array, (0, side * side - cc_array.size)).reshape(side, side)

    print("\n=== common-centroid matrix ===")
    print(cc_array)
    print(f"TMOC = {EvA(cc_array, cc_array.shape)}")
    for idx, dev in enumerate(group_devices, start=1):
        print(f"  id {idx} -> {dev.name}")

    gds_devices = [
        {
            "name": dev.name,
            "nch": dev.nch,
            "model": dev.model,
            "w": dev.unit_width,
            "l": dev.unit_length,
            "nf": dev.nf,
            "multi": dev.multi,
            "attr": dev.attr,
        }
        for dev in group_devices
    ]

    info = generate_group_gds(
        devices=gds_devices,
        cc_matrix=cc_array,
        output_path=output,
        array_name="M0_M7_cc_array",
        space_x=0.3,
        space_y=0.3,
    )

    print("\n=== GDS written ===")
    for key, value in info.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
