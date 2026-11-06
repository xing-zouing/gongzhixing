"""Build one GDS cell per netlist device (prototype granularity)."""

from __future__ import annotations

import math

import gdspy

from cc.device_generation.Capacitor import Capacitor
from cc.device_generation.Mosfet import Mosfet
from cc.io.netlist_parser import Device


def normalized_cell(cell: gdspy.Cell, name: str) -> tuple[gdspy.Cell, float, float]:
    """Flatten a cell and return a copy whose lower-left is (0, 0)."""
    cell.flatten()
    bbox = cell.get_bounding_box()
    if bbox is None:
        raise ValueError(f"cell {name} has no geometry")
    ox, oy = bbox[0]
    width = float(bbox[1][0] - bbox[0][0])
    height = float(bbox[1][1] - bbox[0][1])
    normalized = gdspy.Cell(name, True)
    for (layer, datatype), polygons in cell.get_polygons(by_spec=True).items():
        for polygon in polygons:
            shifted = [(x - ox, y - oy) for x, y in polygon]
            normalized.add(gdspy.Polygon(shifted, layer=layer, datatype=datatype))
    return normalized, width, height


def device_cell(
    dev: Device,
    nf_override: int | None = None,
) -> tuple[gdspy.Cell, float, float]:
    """Create a normalized cell for one netlist device.

    ``nf_override`` is used by the macro folder to build sub-cells with a
    smaller finger count (e.g. 5 fingers instead of 10).
    """
    if dev.kind == "mos":
        mos = Mosfet(
            nch=dev.nch,
            name=dev.name,
            w=dev.unit_width,
            l=dev.unit_length,
            nf=dev.nf if nf_override is None else int(nf_override),
            attr=dev.attr,
            spectre=False,
            pinConType=None,
            bulkCon=[0],
        )
        return normalized_cell(mos.cell, dev.name)
    if dev.kind == "cap":
        cap = Capacitor(
            name=dev.name,
            w=dev.unit_width,
            sp=float(dev.params.get("s", 0.07)),
            nf=dev.nr,
            l=dev.unit_length,
            m_bot=int(dev.params.get("stm", 3)),
            m_top=int(dev.params.get("spm", 5)),
            attr=["2t"] if "2t" in dev.model.lower() else [],
            f_tip=float(dev.params.get("ftip", 0.14)),
            flatten=True,
        )
        return normalized_cell(cap.cell, dev.name)
    raise NotImplementedError(f"unsupported device kind: {dev.kind}")


def mirror_cell_horizontally(cell: gdspy.Cell, new_name: str) -> gdspy.Cell:
    """Mirror a cell about its vertical center line."""
    bbox = cell.get_bounding_box()
    if bbox is None:
        raise ValueError(f"cell {cell.name} has no geometry")
    axis = (bbox[0][0] + bbox[1][0]) / 2.0
    mirrored = gdspy.Cell(new_name, True)
    for (layer, datatype), polygons in cell.get_polygons(by_spec=True).items():
        for polygon in polygons:
            points = [(2.0 * axis - x, y) for x, y in polygon]
            mirrored.add(gdspy.Polygon(points, layer=layer, datatype=datatype))
    return mirrored


def device_macro(
    dev: Device,
    spacing: float = 1.0,
    target_aspect: float = 1.0,
) -> tuple[gdspy.Cell, float, float]:
    """Create a macro for one physical device, folding fingers if useful.

    The total finger count is ``dev.nf * dev.multi``.  For every divisor
    ``fpr`` of ``dev.nf`` we build sub-cells with ``fpr`` fingers and arrange
    ``dev.nf * dev.multi / fpr`` of them in a grid.  The grid whose aspect is
    closest to ``target_aspect`` is selected.

    Example: M0/M7 have ``nf=10, multi=4`` (40 fingers).  Instead of four
    10-finger cells in a 2x2 block (39 x 8.3 um), the folder can choose four
    10-finger cells stacked 1x4 (19 x 17.6 um), or eight 5-finger cells in a
    2x4 grid, whichever is closer to square.
    """
    if dev.kind != "mos":
        return device_cell(dev)
    if dev.multi == 1 and dev.nf == 1:
        return device_cell(dev)

    candidates = []
    # Only the full finger count or half of it are allowed.  This avoids
    # exploding a 10-finger device into twenty 2-finger cells; with nf=10
    # the choices are 10-finger or 5-finger sub-cells.
    finger_options = {dev.nf}
    if dev.nf >= 2 and dev.nf % 2 == 0:
        finger_options.add(dev.nf // 2)

    for fingers_per_cell in sorted(finger_options, reverse=True):
        if dev.nf % fingers_per_cell != 0:
            continue
        units = dev.multi * (dev.nf // fingers_per_cell)
        if units < 1 or units > 256:
            continue

        unit_cell, unit_w, unit_h = device_cell(dev, nf_override=fingers_per_cell)

        for columns in range(1, units + 1):
            rows = math.ceil(units / columns)
            occupied_columns = min(columns, units)
            width = occupied_columns * (unit_w + spacing) - spacing
            height = rows * (unit_h + spacing) - spacing
            if height <= 0:
                continue
            aspect = width / height
            score = abs(math.log(max(aspect, 1e-9) / target_aspect))
            # Prefer fewer sub-cells when the aspect is comparable.
            score += 1e-4 * units
            candidates.append(
                (score, fingers_per_cell, columns, rows, unit_cell, unit_w, unit_h)
            )

    if not candidates:
        return device_cell(dev)

    _, fingers_per_cell, columns, rows, unit_cell, unit_w, unit_h = min(
        candidates, key=lambda item: item[0]
    )
    units = dev.multi * (dev.nf // fingers_per_cell)

    macro = gdspy.Cell(f"{dev.name}_macro", True)
    for index in range(units):
        row, col = divmod(index, columns)
        macro.add(
            gdspy.CellReference(
                unit_cell,
                origin=(col * (unit_w + spacing), row * (unit_h + spacing)),
            )
        )
    return normalized_cell(macro, f"{dev.name}_macro")
