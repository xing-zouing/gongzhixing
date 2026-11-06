"""Geometry helpers for GDS cells (rotation, normalization)."""

from __future__ import annotations

import gdspy


def rotate_cell_90(cell: gdspy.Cell, new_name: str) -> gdspy.Cell:
    """Return a 90-degree counter-clockwise rotated copy of ``cell``.

    The result is normalized so its bounding box starts at (0, 0).  A 90-degree
    rotation swaps the bounding-box width and height, which is what the
    floorplanner uses for block orientation search.
    """
    rotated: list[tuple[int, int, list[tuple[float, float]]]] = []
    min_x = float("inf")
    min_y = float("inf")

    for (layer, datatype), polygons in cell.get_polygons(by_spec=True).items():
        for polygon in polygons:
            points = [(-y, x) for x, y in polygon]
            rotated.append((layer, datatype, points))
            for x, y in points:
                min_x = min(min_x, x)
                min_y = min(min_y, y)

    if not rotated:
        raise ValueError(f"cell {cell.name} has no polygons")

    new_cell = gdspy.Cell(new_name, True)
    for layer, datatype, points in rotated:
        shifted = [(x - min_x, y - min_y) for x, y in points]
        new_cell.add(gdspy.Polygon(shifted, layer=layer, datatype=datatype))
    return new_cell
