#!/usr/bin/env python3
"""Step 13: matched-device-only layout.

All common-centroid-capable MOS groups go into ONE algorithm array:

    M0/M7  (8 units each, 5 fingers/unit)
    M11/M10 (2 units each, 3 fingers/unit)
    input [8, 8, 2, 2] -> 5x5, TMOC = 0

Single-finger pairs (M8/M6, M2/M5, M4/M3) and the capacitor pair (C1/C2)
cannot be common-centroid in the same array; they are built as axial pair
blocks.  Unmatched devices M9/M12 are deliberately excluded.

The blocks are then packed as super-cells (no global axis requirement).

Run from the project root:

    .venv/bin/python scripts/ota3_step13_all_matched_layout.py
"""

from __future__ import annotations

from pathlib import Path
import sys

import gdspy

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cc.eva import EvA  # noqa: E402
from cc.gds.array_writer import (  # noqa: E402
    generate_axial_cap_pair_gds,
    generate_axial_pair_gds,
    generate_group_gds,
)
from cc.gds.cell_ops import rotate_cell_90  # noqa: E402
from cc.gds.device_cells import device_cell  # noqa: E402
from cc.gds.render import render_gds  # noqa: E402
from cc.io.netlist_parser import load_case  # noqa: E402
from cc.main import cc_main_flow  # noqa: E402
from cc.place.floorplan import Block, sa_shelf_pack  # noqa: E402
from cc.plan.ota3_groups import GROUPS  # noqa: E402


def _mos_device(dev, nf: int) -> dict:
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


def _cap_device(dev) -> dict:
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


def _load_block(path: Path, name: str) -> Block:
    library = gdspy.GdsLibrary()
    library.read_gds(str(path))
    top = max(library.top_level(), key=lambda cell: (len(cell.references), len(cell.polygons)))
    top.flatten()
    bbox = top.get_bounding_box()
    ox, oy = bbox[0]
    width = float(bbox[1][0] - bbox[0][0])
    height = float(bbox[1][1] - bbox[0][1])
    cell = gdspy.Cell(name, True)
    for (layer, datatype), polygons in top.get_polygons(by_spec=True).items():
        for polygon in polygons:
            shifted = [(x - ox, y - oy) for x, y in polygon]
            cell.add(gdspy.Polygon(shifted, layer=layer, datatype=datatype))
    rotated = rotate_cell_90(cell, f"{name}_rot90")
    return Block(name=name, width=width, height=height, cell=cell, cell_rot90=rotated)


def main() -> None:
    out_dir = ROOT / "build" / "ota3" / "matched_only"
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in out_dir.glob("*"):
        if stale.is_file():
            stale.unlink()

    devices, _ = load_case(
        ROOT / "examples" / "ota3.sp",
        ROOT / "examples" / "ota3.sym",
        subckt="OTA2_my",
    )

    # ---- one common-centroid array for M0/M7 + M11/M10 ----
    core_counts = [8, 8, 2, 2]  # M0, M7, M11, M10
    matrix = cc_main_flow(
        core_counts,
        square_array=True,
        row_numbers=0,
        random_state=0,
        max_iter=200000,
    )
    import numpy as np

    matrix = np.asarray(matrix, dtype=int)
    if matrix.ndim == 1:
        side = int(np.ceil(np.sqrt(matrix.size)))
        matrix = np.pad(matrix, (0, side * side - matrix.size)).reshape(side, side)
    tmoc = EvA(matrix, matrix.shape)
    print("core input:", core_counts, "matrix:", matrix.shape, "TMOC:", tmoc)
    print(matrix)

    core_path = out_dir / "cc_core_M0_M7_M11_M10.gds"
    generate_group_gds(
        devices=[
            _mos_device(devices["M0"], 5),
            _mos_device(devices["M7"], 5),
            _mos_device(devices["M11"], 3),
            _mos_device(devices["M10"], 3),
        ],
        cc_matrix=matrix,
        output_path=core_path,
        array_name="cc_core_M0_M7_M11_M10",
    )
    render_gds(core_path, core_path.with_suffix(".png"), dpi=200)

    blocks = [_load_block(core_path, "cc_core")]

    # ---- axial pairs ----
    axial_specs = [
        ("pair_M8_M6", "M8", "M6", "mos"),
        ("pair_M2_M5", "M2", "M5", "mos"),
        ("pair_M4_M3", "M4", "M3", "mos"),
        ("cap_pair_C2_C1", "C2", "C1", "cap"),
    ]
    for group_name, left_name, right_name, kind in axial_specs:
        path = out_dir / f"{group_name}_axial.gds"
        if kind == "cap":
            generate_axial_cap_pair_gds(
                _cap_device(devices[left_name]),
                _cap_device(devices[right_name]),
                path,
                array_name=group_name,
                space=0.6,
            )
        else:
            generate_axial_pair_gds(
                _mos_device(devices[left_name], devices[left_name].nf),
                _mos_device(devices[right_name], devices[right_name].nf),
                path,
                array_name=group_name,
                space=0.6,
            )
        render_gds(path, path.with_suffix(".png"), dpi=200)
        blocks.append(_load_block(path, group_name))

    # ---- pack the matched blocks; no global axis ----
    floorplan = sa_shelf_pack(
        blocks,
        spacing=2.0,
        target_aspect=1.2,
        iterations=3000,
        seed=0,
    )
    margin = 5.0
    top_width = floorplan.width + 2 * margin
    top_height = floorplan.height + 2 * margin

    top_lib = gdspy.GdsLibrary()
    top = gdspy.Cell("ota3_matched_only")
    top_lib.add(top)
    top.add(gdspy.Rectangle((0, 0), (top_width, top_height), layer=99, datatype=0))
    for block in floorplan.blocks:
        cell = block.cell_rot90 if block.rotation else block.cell
        x = block.x + margin
        y = block.y + margin
        top.add(gdspy.CellReference(cell, origin=(x, y)))
        top.add(gdspy.Label(block.name, (x + 0.4, y + 0.4), layer=100))
        top_lib.add(cell)

    top_path = out_dir / "ota3_matched_only.gds"
    top_lib.write_gds(str(top_path))
    info = render_gds(top_path, out_dir / "ota3_matched_only.png", dpi=180)

    print("\n=== matched blocks ===")
    for block in floorplan.blocks:
        print(
            f"  {block.name:28s} {block.eff_width:8.2f} x {block.eff_height:8.2f} "
            f"rot={block.rotation} at ({block.x + margin:.2f}, {block.y + margin:.2f})"
        )
    print(f"\ntop size: {top_width:.2f} x {top_height:.2f} um "
          f"(area {top_width * top_height:.1f} um^2)")
    print(f"top GDS: {top_path}")
    print(f"top PNG: {info['png']}")


if __name__ == "__main__":
    main()
