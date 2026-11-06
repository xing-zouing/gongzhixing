#!/usr/bin/env python3
"""Step 12: per-group matched blocks from the project algorithm.

For every two-device matching group:

* if each device can be split into multiple unit cells, enumerate the
  possible "fingers per unit" values and call ``cc_main_flow([k, k])``;
* keep the candidates whose output is a true common-centroid pattern
  (TMOC == 0);
* choose the candidate whose physical block aspect is closest to a target;
* if each device has only one unit (e.g. nf=1), fall back to an axial pair.

The group arrangement is exactly the algorithm output matrix; the outer loop
only controls the finger split to adjust the physical aspect ratio.

Run from the project root:

    .venv/bin/python scripts/ota3_step12_algo_blocks.py
"""

from __future__ import annotations

import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cc.eva import EvA  # noqa: E402
from cc.gds.array_writer import (  # noqa: E402
    generate_axial_cap_pair_gds,
    generate_axial_pair_gds,
    generate_group_gds,
)
from cc.gds.device_cells import device_cell  # noqa: E402
from cc.gds.render import render_gds  # noqa: E402
from cc.io.netlist_parser import Device, load_case  # noqa: E402
from cc.main import cc_main_flow  # noqa: E402
from cc.plan.ota3_groups import GROUPS  # noqa: E402


def _divisors(value: int) -> list[int]:
    return [d for d in range(1, value + 1) if value % d == 0]


def _as_matrix(result) -> np.ndarray:
    matrix = np.asarray(result, dtype=int)
    if matrix.ndim == 1:
        side = int(math.ceil(math.sqrt(matrix.size)))
        matrix = np.pad(matrix, (0, side * side - matrix.size)).reshape(side, side)
    return matrix


def _gds_device(dev: Device, nf: int) -> dict:
    return {
        "name": dev.name,
        "nch": dev.nch,
        "model": dev.model,
        "w": dev.unit_width,
        "l": dev.unit_length,
        "nf": nf,
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


def main() -> None:
    out_dir = ROOT / "build" / "ota3" / "algo_blocks"
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in out_dir.glob("*"):
        if stale.is_file():
            stale.unlink()

    devices, _ = load_case(
        ROOT / "examples" / "ota3.sp",
        ROOT / "examples" / "ota3.sym",
        subckt="OTA2_my",
    )

    target_aspect = 1.0
    rows = []

    for spec in GROUPS:
        devs = [devices[name] for name in spec["devices"]]
        if len(devs) != 2:
            continue
        left, right = devs

        if left.kind != "mos" or right.kind != "mos":
            # Capacitor pair: axial for now.
            if spec["type"] == "axial_cap":
                path = out_dir / f"{spec['name']}_axial.gds"
                generate_axial_cap_pair_gds(
                    _cap_device(left), _cap_device(right), path,
                    array_name=spec["name"], space=0.6,
                )
                render_gds(path, path.with_suffix(".png"), dpi=200)
                print(f"{spec['name']}: axial cap pair")
            continue

        total_fingers = left.multi * left.nf
        if total_fingers != right.multi * right.nf:
            raise ValueError(f"{spec['name']}: unequal finger totals")

        if total_fingers <= 1:
            # Two single-unit devices cannot share a centroid: axial pair.
            path = out_dir / f"{spec['name']}_axial.gds"
            generate_axial_pair_gds(
                _gds_device(left, left.nf),
                _gds_device(right, right.nf),
                path,
                array_name=spec["name"],
                space=0.6,
            )
            render_gds(path, path.with_suffix(".png"), dpi=200)
            rows.append(
                {
                    "group": spec["name"],
                    "mode": "axial",
                    "input": (1, 1),
                    "fingers_per_unit": left.nf,
                    "matrix": "1x2",
                    "tmoc": "-",
                    "block": "-",
                    "aspect": "-",
                    "score": "-",
                    "path": path.name,
                }
            )
            continue

        candidates = []
        for fingers_per_unit in _divisors(left.nf):
            units = total_fingers // fingers_per_unit
            result = cc_main_flow(
                [units, units],
                square_array=True,
                row_numbers=0,
                random_state=0,
                max_iter=200000,
            )
            matrix = _as_matrix(result)
            tmoc = EvA(matrix, matrix.shape)
            if tmoc != 0.0:
                continue

            unit_cell, unit_w, unit_h = device_cell(
                left, nf_override=fingers_per_unit
            )
            rows_count, cols_count = matrix.shape
            gap = 0.6
            block_w = cols_count * (unit_w + gap) - gap
            block_h = rows_count * (unit_h + gap) - gap
            aspect = block_w / block_h if block_h else 0.0
            score = abs(math.log(max(aspect, 1e-9) / target_aspect))

            path = out_dir / f"{spec['name']}_k{units}_f{fingers_per_unit}.gds"
            generate_group_gds(
                devices=[
                    _gds_device(left, fingers_per_unit),
                    _gds_device(right, fingers_per_unit),
                ],
                cc_matrix=matrix,
                output_path=path,
                array_name=f"{spec['name']}_k{units}_f{fingers_per_unit}",
            )
            render_gds(path, path.with_suffix(".png"), dpi=200)

            candidates.append(
                {
                    "group": spec["name"],
                    "mode": "common_centroid",
                    "input": (units, units),
                    "fingers_per_unit": fingers_per_unit,
                    "matrix": f"{rows_count}x{cols_count}",
                    "tmoc": tmoc,
                    "block": (block_w, block_h),
                    "aspect": aspect,
                    "score": score,
                    "path": path.name,
                    "matrix_values": matrix,
                }
            )

        if not candidates:
            # No common-centroid candidate; axial fallback.
            path = out_dir / f"{spec['name']}_axial.gds"
            generate_axial_pair_gds(
                _gds_device(left, left.nf),
                _gds_device(right, right.nf),
                path,
                array_name=spec["name"],
                space=0.6,
            )
            render_gds(path, path.with_suffix(".png"), dpi=200)
            rows.append(
                {
                    "group": spec["name"],
                    "mode": "axial",
                    "input": "-",
                    "fingers_per_unit": left.nf,
                    "matrix": "1x2",
                    "tmoc": "-",
                    "block": "-",
                    "aspect": "-",
                    "score": "-",
                    "path": path.name,
                }
            )
            continue

        candidates.sort(key=lambda item: item["score"])
        best = candidates[0]
        best["chosen"] = True
        for candidate in candidates[1:]:
            candidate["chosen"] = False
        rows.extend(candidates)

        print(f"\n{spec['name']}:")
        for candidate in candidates:
            mark = "  <== chosen" if candidate["chosen"] else ""
            print(
                f"  input={candidate['input']} fpr={candidate['fingers_per_unit']} "
                f"matrix={candidate['matrix']} TMOC={candidate['tmoc']} "
                f"block={candidate['block'][0]:.2f}x{candidate['block'][1]:.2f} "
                f"aspect={candidate['aspect']:.2f} score={candidate['score']:.3f}{mark}"
            )
            print(f"    {candidate['matrix_values'].tolist()}")

    print("\n=== summary ===")
    print(f"{'group':28s} {'mode':16s} {'input':>10s} {'fpr':>4s} {'matrix':>7s} {'aspect':>7s} {'chosen':>6s}")
    for row in rows:
        print(
            f"{row['group']:28s} {row['mode']:16s} {str(row['input']):>10s} "
            f"{str(row['fingers_per_unit']):>4s} {row['matrix']:>7s} "
            f"{str(row['aspect']):>7s} {str(row.get('chosen', False)):>6s}"
        )

    print(f"\nblocks: {out_dir}")


if __name__ == "__main__":
    main()
