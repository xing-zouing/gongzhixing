"""Render a GDS file to PNG for quick inspection.

This is a lightweight verification aid, not a replacement for KLayout.
It reads the polygons of the top cell and draws them layer by layer.
"""

from __future__ import annotations

from pathlib import Path

import gdspy
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
import numpy as np  # noqa: E402


_PALETTE = plt.get_cmap("tab20").colors


def _color_for(layer: int, datatype: int) -> tuple[float, float, float]:
    return _PALETTE[(int(layer) * 7 + int(datatype)) % len(_PALETTE)]


def render_gds(
    gds_path: str | Path,
    png_path: str | Path | None = None,
    dpi: int = 160,
    max_fig: float = 14.0,
    min_fig: float = 4.0,
    alpha: float = 0.85,
) -> dict:
    """Render the top cell of ``gds_path`` to ``png_path``.

    Returns a small metadata dict (bbox, polygon count, PNG path).
    """
    gds_path = Path(gds_path)
    if png_path is None:
        png_path = gds_path.with_suffix(".png")
    png_path = Path(png_path)
    png_path.parent.mkdir(parents=True, exist_ok=True)

    library = gdspy.GdsLibrary()
    library.read_gds(str(gds_path))
    tops = library.top_level()
    if not tops:
        raise ValueError(f"{gds_path} has no top-level cell")
    # Prefer the cell that instantiates the array; an unused helper cell
    # (e.g. dummy_unit) must not win just because it was inserted first.
    top = max(tops, key=lambda cell: (len(cell.references), len(cell.polygons)))

    bbox = top.get_bounding_box()
    if bbox is None:
        raise ValueError(f"{gds_path} top cell has no geometry")
    (x0, y0), (x1, y1) = bbox
    width = max(float(x1 - x0), 1e-9)
    height = max(float(y1 - y0), 1e-9)

    # Keep the aspect ratio while keeping the figure a reasonable size.
    scale = max_fig / max(width, height)
    scale = max(scale, min_fig / max(width, height))
    fig_w = max(min_fig, width * scale)
    fig_h = max(min_fig, height * scale)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    polygons_by_spec = top.get_polygons(by_spec=True)

    total = 0
    # Layers 99/100 are used for the top boundary and block labels; they are
    # not layout geometry and would otherwise fill the whole image.
    ignored_layers = {99, 100}
    for (layer, datatype), polygons in sorted(polygons_by_spec.items()):
        if not polygons or layer in ignored_layers:
            continue
        points = [np.asarray(poly, dtype=float) for poly in polygons]
        collection = PolyCollection(
            points,
            facecolors=_color_for(layer, datatype),
            edgecolors="none",
            alpha=alpha,
        )
        ax.add_collection(collection)
        total += len(points)

    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"{top.name}  ({width:.2f} x {height:.2f} um, {total} polygons)", fontsize=10)
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    return {
        "gds": str(gds_path),
        "png": str(png_path),
        "cell": top.name,
        "width": width,
        "height": height,
        "polygons": total,
        "layers": len(polygons_by_spec),
    }
