"""Write a mixed MOS/capacitor unit-cell array to GDS.

``generate_group_gds`` only builds MOSFET cells.  This writer accepts device
dicts with ``kind`` = ``"mos"`` or ``"cap"`` and places the corresponding
unit cell at every matrix position carrying that device id.  Matrix value 0
is a dummy slot and is left empty (no cell is placed).
"""

from __future__ import annotations

from pathlib import Path

import gdspy

from cc.device_generation.Capacitor import Capacitor
from cc.device_generation.Mosfet import Mosfet
from cc.gds.device_cells import normalized_cell


def _make_cell(dev: dict):
    if dev["kind"] == "mos":
        cell = Mosfet(
            nch=dev["nch"],
            name=dev["name"],
            w=dev["w"],
            l=dev["l"],
            nf=dev["nf"],
            attr=list(dev.get("attr", [])),
            spectre=False,
            pinConType=None,
            bulkCon=[0],
        ).cell
    elif dev["kind"] == "cap":
        cell = Capacitor(
            name=dev["name"],
            w=dev["w"],
            sp=dev.get("sp", 0.07),
            nf=dev["nf"],
            l=dev["l"],
            m_bot=dev.get("m_bot", 3),
            m_top=dev.get("m_top", 5),
            attr=list(dev.get("attr", [])),
            f_tip=dev.get("f_tip", 0.14),
            flatten=True,
        ).cell
    else:
        raise ValueError(f"unsupported kind: {dev['kind']}")
    return normalized_cell(cell, dev["name"])


def generate_mixed_group_gds(
    devices: list[dict],
    matrix,
    output_path: str | Path,
    array_name: str = "mixed_array",
    gap: float = 0.3,
) -> dict:
    import numpy as np

    arr = np.asarray(matrix, dtype=int)
    if arr.ndim == 1:
        side = int(np.ceil(np.sqrt(arr.size)))
        arr = np.pad(arr, (0, side * side - arr.size)).reshape(side, side)
    rows, cols = arr.shape

    top_lib = gdspy.GdsLibrary()
    unit_info = {}
    max_w = 0.0
    max_h = 0.0

    for group_id, dev in enumerate(devices, start=1):
        cell, width, height = _make_cell(dev)
        unit_info[group_id] = {
            "cell": cell,
            "w": width,
            "h": height,
            "ox": 0.0,
            "oy": 0.0,
        }
        max_w = max(max_w, width)
        max_h = max(max_h, height)
        top_lib.add(cell)

    step_x = max_w + gap
    step_y = max_h + gap
    total_w = cols * step_x - gap
    total_h = rows * step_y - gap
    start_x = -total_w / 2.0
    start_y = -total_h / 2.0

    top = gdspy.Cell(array_name)
    top_lib.add(top)
    for row in range(rows):
        for col in range(cols):
            group_id = int(arr[row, col])
            if group_id == 0 or group_id not in unit_info:
                continue
            info = unit_info[group_id]
            target_x = start_x + col * step_x
            target_y = start_y + (rows - 1 - row) * step_y
            origin = (
                target_x - info["ox"],
                target_y - info["oy"],
            )
            top.add(gdspy.CellReference(info["cell"], origin=origin))

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    top_lib.write_gds(str(output_path))
    return {
        "path": str(output_path),
        "rows": rows,
        "cols": cols,
        "unit_width": max_w,
        "unit_height": max_h,
        "placements": int(np.count_nonzero(arr)),
        "dummies": int(np.count_nonzero(arr == 0)),
    }
