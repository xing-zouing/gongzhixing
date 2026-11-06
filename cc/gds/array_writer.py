"""Write a common-centroid unit-cell array to GDS.

The first version uses one ``multi`` instance as the unit cell: a unit cell
contains ``nf`` fingers of width ``w`` (Spectre convention, per-finger width).
When a device has ``multi=4, nf=10``, it therefore contributes 4 unit cells,
each with 10 fingers.

The GDS layer stack comes from ``cc.device_generation.glovar`` (TSMC40-style
mock PDK).  This is a *layout prototype*: replace the layer stack / PCell with
the real PDK before tape-out.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import gdspy
import numpy as np

from cc.device_generation.Capacitor import Capacitor
from cc.device_generation.Mosfet import Mosfet


def _as_matrix(cc_matrix, fill: int = 0) -> np.ndarray:
    arr = np.asarray(cc_matrix, dtype=int)
    if arr.ndim == 1:
        side = int(np.ceil(np.sqrt(arr.size)))
        arr = np.pad(arr, (0, side * side - arr.size), constant_values=fill)
        arr = arr.reshape(side, side)
    if arr.ndim != 2:
        raise ValueError(f"cc_matrix must be 1-D or 2-D, got shape {arr.shape}")
    return arr


def _make_mosfet(dev: dict, name: str) -> Mosfet:
    """Create a MOSFET unit cell.

    ``spectre=False`` means ``w`` is already the per-finger width.  This
    matches the Spectre netlists in ``mycase`` (e.g. ``w=2u nf=10`` means
    10 fingers of 2 um each).
    """
    return Mosfet(
        nch=bool(dev["nch"]),
        name=name,
        w=float(dev["w"]),
        l=float(dev["l"]),
        nf=int(dev["nf"]),
        attr=list(dev.get("attr", [])),
        spectre=False,
        pinConType=None,
        bulkCon=[0],
    )


def generate_group_gds(
    devices: Iterable[dict],
    cc_matrix,
    output_path: str | Path,
    array_name: str = "cc_array",
    space_x: float = 0.3,
    space_y: float = 0.3,
    dummy_name: str = "dummy_unit",
) -> dict:
    """Place the unit cells of one matching group into a GDS.

    Parameters
    ----------
    devices:
        One dict per non-zero group id, in the order used by ``cc_matrix``.
        Required keys: ``name``, ``nch``, ``w``, ``l``, ``nf``; optional:
        ``multi``, ``attr``.
    cc_matrix:
        2-D (or 1-D) matrix of group ids.  ``0`` is a dummy cell.
    output_path:
        Destination ``.gds`` file.

    Returns
    -------
    dict with the matrix shape, unit size, and output path.
    """
    devices = list(devices)
    if not devices:
        raise ValueError("devices must not be empty")

    arr = _as_matrix(cc_matrix)
    rows, cols = arr.shape

    top_lib = gdspy.GdsLibrary()
    unit_info: dict[int, dict] = {}

    # Real unit cells: group id 1..N.
    for group_id, dev in enumerate(devices, start=1):
        name = f"{dev['name']}_nf{int(dev['nf'])}"
        mos = _make_mosfet(dev, name)
        bbox = mos.cell.get_bounding_box()
        unit_info[group_id] = {
            "cell": mos.cell,
            "width": float(bbox[1][0] - bbox[0][0]),
            "height": float(bbox[1][1] - bbox[0][1]),
            "offset_x": float(bbox[0][0]),
            "offset_y": float(bbox[0][1]),
        }
        top_lib.add(mos.cell)

    # Add a dummy cell only when the matrix actually uses id 0.
    if np.any(arr == 0):
        ref = devices[0]
        dummy = _make_mosfet(ref, dummy_name)
        bbox = dummy.cell.get_bounding_box()
        unit_info[0] = {
            "cell": dummy.cell,
            "width": float(bbox[1][0] - bbox[0][0]),
            "height": float(bbox[1][1] - bbox[0][1]),
            "offset_x": float(bbox[0][0]),
            "offset_y": float(bbox[0][1]),
        }
        top_lib.add(dummy.cell)

    # Use the largest cell footprint so different unit geometries do not overlap.
    unit_w = max(info["width"] for info in unit_info.values())
    unit_h = max(info["height"] for info in unit_info.values())
    step_x = unit_w + space_x
    step_y = unit_h + space_y

    total_w = cols * step_x - space_x
    total_h = rows * step_y - space_y
    start_x = -total_w / 2.0
    start_y = -total_h / 2.0

    array_cell = gdspy.Cell(array_name)
    top_lib.add(array_cell)

    for row in range(rows):
        for col in range(cols):
            group_id = int(arr[row, col])
            if group_id not in unit_info:
                continue
            info = unit_info[group_id]
            target_x = start_x + col * step_x
            # Keep the matrix row order visually top-to-bottom.
            target_y = start_y + (rows - 1 - row) * step_y
            origin = (
                target_x - info["offset_x"],
                target_y - info["offset_y"],
            )
            array_cell.add(gdspy.CellReference(info["cell"], origin=origin))

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    top_lib.write_gds(str(output_path))

    return {
        "path": str(output_path),
        "rows": rows,
        "cols": cols,
        "unit_width": unit_w,
        "unit_height": unit_h,
        "placements": int(np.count_nonzero(arr)),
        "dummies": int(np.count_nonzero(arr == 0)),
    }


def _make_capacitor(dev: dict, name: str) -> Capacitor:
    """Create a MOM-capacitor unit cell from a parsed netlist device.

    ``nf`` here is the number of fingers (the netlist ``nr``).  The
    ``multi`` multiplier is represented by instantiating the unit cell
    multiple times at the array level.
    """
    return Capacitor(
        name=name,
        w=float(dev["w"]),
        sp=float(dev.get("sp", 0.07)),
        nf=int(dev["nf"]),
        l=float(dev["l"]),
        m_bot=int(dev.get("m_bot", 3)),
        m_top=int(dev.get("m_top", 5)),
        attr=list(dev.get("attr", [])),
        f_tip=float(dev.get("f_tip", 0.14)),
        flatten=True,
    )


def _mirror_cell_horizontally(cell: gdspy.Cell, new_name: str) -> gdspy.Cell:
    """Return a new cell mirrored about its vertical center line."""
    bbox = cell.get_bounding_box()
    axis = (bbox[0][0] + bbox[1][0]) / 2.0
    mirrored = gdspy.Cell(new_name, True)
    for (layer, datatype), polygons in cell.get_polygons(by_spec=True).items():
        for polygon in polygons:
            points = [(2.0 * axis - x, y) for x, y in polygon]
            mirrored.add(gdspy.Polygon(points, layer=layer, datatype=datatype))
    return mirrored


def generate_axial_cap_pair_gds(
    left_dev: dict,
    right_dev: dict,
    output_path: str | Path,
    array_name: str = "axial_cap_pair",
    space: float = 0.6,
    mirror_right: bool = True,
) -> dict:
    """Place two capacitor cells as an axial mirror pair."""
    top_lib = gdspy.GdsLibrary()

    left = _make_capacitor(left_dev, f"{left_dev['name']}_nr{int(left_dev['nf'])}")
    right = _make_capacitor(right_dev, f"{right_dev['name']}_nr{int(right_dev['nf'])}")
    right_cell = (
        _mirror_cell_horizontally(right.cell, f"{right.name}_mirror")
        if mirror_right
        else right.cell
    )

    def _info(cell: gdspy.Cell) -> dict:
        bbox = cell.get_bounding_box()
        return {
            "cell": cell,
            "width": float(bbox[1][0] - bbox[0][0]),
            "height": float(bbox[1][1] - bbox[0][1]),
            "cx": float((bbox[0][0] + bbox[1][0]) / 2.0),
            "cy": float((bbox[0][1] + bbox[1][1]) / 2.0),
        }

    left_info = _info(left.cell)
    right_info = _info(right_cell)
    step = max(left_info["width"], right_info["width"]) + space

    pair_cell = gdspy.Cell(array_name)
    for info, target_x in ((left_info, -step / 2.0), (right_info, +step / 2.0)):
        origin = (target_x - info["cx"], 0.0 - info["cy"])
        pair_cell.add(gdspy.CellReference(info["cell"], origin=origin))
        top_lib.add(info["cell"])

    top_lib.add(pair_cell)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    top_lib.write_gds(str(output_path))

    return {
        "path": str(output_path),
        "left": left_dev["name"],
        "right": right_dev["name"],
        "mirrored": bool(mirror_right),
        "axis_x": 0.0,
        "step": step,
        "left_width": left_info["width"],
        "right_width": right_info["width"],
        "nf": int(left_dev["nf"]),
    }


def generate_axial_pair_gds(
    left_dev: dict,
    right_dev: dict,
    output_path: str | Path,
    array_name: str = "axial_pair",
    space: float = 0.6,
    mirror_right: bool = True,
) -> dict:
    """Place two unit cells as a mirror pair about the vertical axis x=0.

    ``left_dev`` and ``right_dev`` must describe compatible unit cells
    (same ``w``, ``l``, ``nf``).  The right cell is flipped horizontally by
    ``Mosfet.flip_vert()`` so the pair is truly axial-symmetric.
    """
    top_lib = gdspy.GdsLibrary()

    left = _make_mosfet(left_dev, f"{left_dev['name']}_nf{int(left_dev['nf'])}")
    right = _make_mosfet(right_dev, f"{right_dev['name']}_nf{int(right_dev['nf'])}")
    if mirror_right:
        right.flip_vert()

    def _info(mos: Mosfet) -> dict:
        bbox = mos.cell.get_bounding_box()
        return {
            "cell": mos.cell,
            "width": float(bbox[1][0] - bbox[0][0]),
            "height": float(bbox[1][1] - bbox[0][1]),
            "cx": float((bbox[0][0] + bbox[1][0]) / 2.0),
            "cy": float((bbox[0][1] + bbox[1][1]) / 2.0),
        }

    left_info = _info(left)
    right_info = _info(right)
    step = max(left_info["width"], right_info["width"]) + space

    pair_cell = gdspy.Cell(array_name)
    for info, target_x in ((left_info, -step / 2.0), (right_info, +step / 2.0)):
        origin = (target_x - info["cx"], 0.0 - info["cy"])
        pair_cell.add(gdspy.CellReference(info["cell"], origin=origin))
        top_lib.add(info["cell"])

    top_lib.add(pair_cell)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    top_lib.write_gds(str(output_path))

    return {
        "path": str(output_path),
        "left": left_dev["name"],
        "right": right_dev["name"],
        "mirrored": bool(mirror_right),
        "axis_x": 0.0,
        "step": step,
        "left_width": left_info["width"],
        "right_width": right_info["width"],
    }
