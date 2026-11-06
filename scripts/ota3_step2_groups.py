#!/usr/bin/env python3
"""Step 2: classify ota3 matching groups and generate per-group GDS blocks.

Types
-----
common_centroid
    Every device is split into ``unit_count = multi * nf / seg_nf`` identical
    unit cells; all devices share one centroid.  Example: M0/M7 (seg_nf=10,
    4 units each) and M11/M10 (seg_nf=1, 6 units each).

axial
    Each device is one unit cell; the second cell is mirrored about a
    vertical symmetry axis.  Used for the single-unit pairs M8/M6, M2/M5,
    M4/M3.

Groups that share a compatible unit cell are merged by the automatic
cluster interleaver in Step 3; Step 2 generates only the standalone blocks.

Run from the project root:

    .venv/bin/python scripts/ota3_step2_groups.py
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cc.io.netlist_parser import Device, load_case  # noqa: E402
from cc.main import cc_main_flow  # noqa: E402
from cc.eva import EvA  # noqa: E402
from cc.gds.array_writer import (  # noqa: E402
    generate_axial_cap_pair_gds,
    generate_axial_pair_gds,
    generate_group_gds,
)
from cc.place.compatibility import cluster_groups  # noqa: E402
from cc.plan.ota3_groups import GROUPS  # noqa: E402


def _gds_device(dev: Device, seg_nf: int) -> dict:
    return {
        "name": dev.name,
        "nch": dev.nch,
        "model": dev.model,
        "w": dev.unit_width,
        "l": dev.unit_length,
        "nf": seg_nf,
        "multi": dev.multi,
        "attr": dev.attr,
    }


def _cap_device(dev: Device) -> dict:
    return {
        "name": dev.name,
        "w": dev.unit_width,
        "l": dev.unit_length,
        "sp": float(dev.params.get("s", 0.07)),
        "nf": dev.nr,
        "m_bot": int(dev.params.get("stm", 3)),
        "m_top": int(dev.params.get("spm", 5)),
        "attr": ["2t"] if "2t" in dev.model.lower() else [],
        "f_tip": float(dev.params.get("ftip", 0.14)),
        "multi": dev.multi,
    }


def _check_cap_compatible(devs: list[Device]) -> None:
    for dev in devs:
        if dev.kind != "cap":
            raise ValueError(f"{dev.name} is not a capacitor")
    ref = devs[0]
    for dev in devs[1:]:
        if (dev.unit_width, dev.unit_length, dev.nr, dev.multi) != (
            ref.unit_width,
            ref.unit_length,
            ref.nr,
            ref.multi,
        ):
            raise ValueError(
                f"capacitor geometry mismatch: {dev.name} "
                f"({dev.unit_width},{dev.unit_length},nr={dev.nr},multi={dev.multi}) "
                f"vs {ref.name} ({ref.unit_width},{ref.unit_length},nr={ref.nr},multi={ref.multi})"
            )


def _check_compatible(devs: list[Device], seg_nf: int) -> None:
    for dev in devs:
        if dev.kind != "mos":
            raise ValueError(f"{dev.name} is not a MOS device")
        if dev.nf % seg_nf != 0:
            raise ValueError(
                f"{dev.name}: nf={dev.nf} is not divisible by seg_nf={seg_nf}"
            )
    ref = devs[0]
    for dev in devs[1:]:
        if (dev.unit_width, dev.unit_length) != (ref.unit_width, ref.unit_length):
            raise ValueError(
                f"group geometry mismatch: {dev.name} "
                f"({dev.unit_width},{dev.unit_length}) vs "
                f"{ref.name} ({ref.unit_width},{ref.unit_length})"
            )


def _reshape_matrix(cc_result) -> np.ndarray:
    arr = np.asarray(cc_result, dtype=int)
    if arr.ndim == 1:
        side = int(np.ceil(np.sqrt(arr.size)))
        arr = np.pad(arr, (0, side * side - arr.size)).reshape(side, side)
    return arr


def main() -> None:
    netlist = ROOT / "examples" / "ota3.sp"
    sym = ROOT / "examples" / "ota3.sym"
    out_dir = ROOT / "build" / "ota3"

    devices, _ = load_case(netlist, sym, subckt="OTA2_my")

    clusters = cluster_groups(devices, GROUPS)
    interleaved_names = {
        spec["name"]
        for cluster in clusters
        if cluster["interleavable"]
        for spec in cluster["groups"]
    }

    print("\n=== automatic compatibility clusters ===")
    for cluster in clusters:
        names = [spec["name"] for spec in cluster["groups"]]
        mode = "INTERLEAVE" if cluster["interleavable"] else "standalone"
        print(f"  {mode:10s} {names}")

    for group in GROUPS:
        name = group["name"]
        if name in interleaved_names:
            print(
                f"\n=== {name}: skipped "
                f"(merged into an automatic interleaved cluster) ==="
            )
            continue
        devs = [devices[dev_name] for dev_name in group["devices"]]
        print(f"\n=== {name} ({group['type']}) ===")
        print(f"  devices      : {[d.name for d in devs]}")

        if group["type"] == "axial_cap":
            _check_cap_compatible(devs)
            unit_counts = [dev.multi for dev in devs]
            print(f"  unit         : multi-level, counts={unit_counts}")
            if any(count != 1 for count in unit_counts):
                print(f"  SKIP: axial_cap currently supports one unit per device, got {unit_counts}")
                continue
            info = generate_axial_cap_pair_gds(
                left_dev=_cap_device(devs[0]),
                right_dev=_cap_device(devs[1]),
                output_path=out_dir / f"{name}_axial.gds",
                array_name=name,
                space=0.6,
                mirror_right=True,
            )
            print(f"  GDS          : {info['path']}")
            for key in ("step", "mirrored", "nf"):
                if key in info:
                    print(f"    {key}: {info[key]}")
            continue

        seg_nf = int(group["seg_nf"])
        _check_compatible(devs, seg_nf)
        unit_counts = [dev.multi * dev.nf // seg_nf for dev in devs]
        print(f"  seg_nf       : {seg_nf}")
        print(f"  unit counts  : {unit_counts}")

        if group["type"] == "common_centroid":
            cc_result = cc_main_flow(
                unit_counts,
                square_array=True,
                row_numbers=0,
                random_state=0,
                max_iter=200000,
            )
            arr = _reshape_matrix(cc_result)
            print(f"  matrix       :\n{arr}")
            print(f"  TMOC         : {EvA(arr, arr.shape)}")
            info = generate_group_gds(
                devices=[_gds_device(d, seg_nf) for d in devs],
                cc_matrix=arr,
                output_path=out_dir / f"{name}_cc.gds",
                array_name=name,
            )
        elif group["type"] == "axial":
            if any(count != 1 for count in unit_counts):
                print(f"  SKIP: axial currently supports one unit per device, got {unit_counts}")
                continue
            info = generate_axial_pair_gds(
                left_dev=_gds_device(devs[0], seg_nf),
                right_dev=_gds_device(devs[1], seg_nf),
                output_path=out_dir / f"{name}_axial.gds",
                array_name=name,
                space=0.6,
                mirror_right=True,
            )
        else:
            raise ValueError(f"unknown group type: {group['type']}")

        print(f"  GDS          : {info['path']}")
        for key in ("rows", "cols", "unit_width", "unit_height", "placements", "dummies", "step", "mirrored"):
            if key in info:
                print(f"    {key}: {info[key]}")


if __name__ == "__main__":
    main()
