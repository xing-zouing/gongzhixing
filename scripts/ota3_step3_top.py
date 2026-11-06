#!/usr/bin/env python3
"""Step 3: block-level top layout for ota3.

* each matching group GDS becomes a fixed-size block;
* the unmatched devices M9/M12 become individual blocks;
* a shelf packer places all blocks;
* a top-level GDS + PNG are written.

Run from the project root:

    .venv/bin/python scripts/ota3_step3_top.py

The block abstraction (``Block``) is intentionally independent of the
packer, so the shelf packer can later be replaced by simulated annealing,
sequence-pair, or the continuous optimizer in
``CEDA/analog_placement/LAYOUT``.
"""

from __future__ import annotations

from pathlib import Path
import sys

import gdspy

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cc.device_generation.Mosfet import Mosfet  # noqa: E402
from cc.gds.cell_ops import rotate_cell_90  # noqa: E402
from cc.gds.cluster_writer import write_cluster_gds  # noqa: E402
from cc.gds.render import render_gds  # noqa: E402
from cc.io.netlist_parser import load_case  # noqa: E402
from cc.place.compatibility import cluster_groups  # noqa: E402
from cc.place.floorplan import Block, axis_aligned_column, sa_shelf_pack  # noqa: E402
from cc.plan.ota3_groups import GROUPS  # noqa: E402


UNMATCHED = [
    {"name": "M9", "nch": True, "w": 3.0, "l": 4.0, "nf": 1, "attr": ["lvt"]},
    {"name": "M12", "nch": True, "w": 3.0, "l": 4.0, "nf": 1, "attr": ["lvt"]},
]


def _normalized_cell(cell: gdspy.Cell, name: str) -> tuple[gdspy.Cell, float, float]:
    """Flatten ``cell`` and return a copy whose lower-left is (0, 0)."""
    cell.flatten()
    bbox = cell.get_bounding_box()
    if bbox is None:
        raise ValueError(f"cell {cell.name} has no geometry")
    ox, oy = bbox[0]
    width = float(bbox[1][0] - bbox[0][0])
    height = float(bbox[1][1] - bbox[0][1])

    normalized = gdspy.Cell(name, True)
    for (layer, datatype), polygons in cell.get_polygons(by_spec=True).items():
        for polygon in polygons:
            shifted = [(x - ox, y - oy) for x, y in polygon]
            normalized.add(gdspy.Polygon(shifted, layer=layer, datatype=datatype))
    return normalized, width, height


def _make_block(cell: gdspy.Cell, name: str) -> Block:
    normalized, width, height = _normalized_cell(cell, name)
    rotated = rotate_cell_90(normalized, f"{name}_rot90")
    return Block(
        name=name,
        width=width,
        height=height,
        cell=normalized,
        cell_rot90=rotated,
    )


def load_block(path: Path, name: str) -> Block:
    library = gdspy.GdsLibrary()
    library.read_gds(str(path))
    tops = library.top_level()
    if not tops:
        raise ValueError(f"{path} has no top-level cell")
    # Pick the array cell, not an unused dummy/helper cell.
    top = max(tops, key=lambda cell: (len(cell.references), len(cell.polygons)))
    return _make_block(top, name)


def unmatched_block(params: dict) -> Block:
    mos = Mosfet(
        nch=params["nch"],
        name=params["name"],
        w=params["w"],
        l=params["l"],
        nf=params["nf"],
        attr=list(params.get("attr", [])),
        spectre=False,
        pinConType=None,
        bulkCon=[0],
    )
    return _make_block(mos.cell, params["name"])


def write_top(floorplan, output_gds: Path, output_png: Path, top_name: str,
              title: str, margin: float = 5.0) -> dict:
    top_width = floorplan.width + 2 * margin
    top_height = floorplan.height + 2 * margin

    top_lib = gdspy.GdsLibrary()
    top = gdspy.Cell(top_name)
    top_lib.add(top)

    # Boundary on layer 99/0, labels on layer 100.
    top.add(gdspy.Rectangle((0, 0), (top_width, top_height), layer=99, datatype=0))

    for block in floorplan.blocks:
        x = block.x + margin
        y = block.y + margin
        cell = block.cell_rot90 if block.rotation else block.cell
        top.add(gdspy.CellReference(cell, origin=(x, y)))
        top.add(gdspy.Label(block.name, (x + 0.5, y + 0.5), layer=100))
        top_lib.add(cell)

    top_lib.write_gds(str(output_gds))

    print(f"\n=== {title} ===")
    print(f"{'name':28s} {'rot':>4s} {'w':>8s} {'h':>8s} {'x':>8s} {'y':>8s}")
    for block in floorplan.blocks:
        print(
            f"{block.name:28s} {block.rotation:4d} {block.eff_width:8.2f} {block.eff_height:8.2f} "
            f"{block.x + margin:8.2f} {block.y + margin:8.2f}"
        )
    print(f"top size: {top_width:.2f} x {top_height:.2f} um "
          f"(area {top_width * top_height:.1f} um^2)")
    print(f"top GDS: {output_gds}")

    info = render_gds(output_gds, output_png, dpi=160)
    print(f"top PNG: {info['png']}  ({info['width']:.2f}x{info['height']:.2f} um, "
          f"{info['polygons']} polygons, {info['layers']} layers)")
    return info


def main() -> None:
    block_dir = ROOT / "build" / "ota3"
    netlist = ROOT / "examples" / "ota3.sp"
    sym = ROOT / "examples" / "ota3.sym"
    devices, _ = load_case(netlist, sym, subckt="OTA2_my")

    blocks: list[Block] = []
    clusters = cluster_groups(devices, GROUPS)
    print("=== automatic compatibility clusters ===")
    for cluster in clusters:
        names = [spec["name"] for spec in cluster["groups"]]
        mode = "INTERLEAVE" if cluster["interleavable"] else "standalone"
        print(f"  {mode:10s} {names}")

    for cluster in clusters:
        names = [spec["name"] for spec in cluster["groups"]]
        if cluster["interleavable"]:
            path = write_cluster_gds(devices, cluster, block_dir)
            if path is not None:
                blocks.append(load_block(path, path.stem))
                continue
            print(f"  fallback: could not interleave {names}; using standalone blocks")

        for spec in cluster["groups"]:
            block_file = spec.get("block_file")
            if block_file:
                path = block_dir / block_file
                if not path.exists():
                    raise FileNotFoundError(
                        f"missing block GDS: {path}; run step 2 first"
                    )
                blocks.append(load_block(path, Path(block_file).stem))
                continue

            if spec["type"] == "axial_cap":
                raise FileNotFoundError(
                    f"capacitor group {spec['name']} needs a standalone block"
                )

            # Fallback for a MOS group that has no standalone file.
            single = {"groups": [spec], "interleavable": False}
            path = write_cluster_gds(devices, single, block_dir)
            if path is None:
                raise RuntimeError(f"cannot build a block for {spec['name']}")
            blocks.append(load_block(path, path.stem))

    for params in UNMATCHED:
        blocks.append(unmatched_block(params))

    # 1) Compact placement (area-oriented).
    compact = sa_shelf_pack(
        blocks,
        spacing=3.0,
        target_aspect=1.3,
        iterations=4000,
        seed=0,
    )
    write_top(
        compact,
        block_dir / "ota3_top.gds",
        block_dir / "ota3_top.png",
        "ota3_top",
        "compact placement",
    )

    # 2) Global axis-aligned placement: every block center on one vertical axis.
    #    This is the practical "global matching" layout for heterogeneous groups.
    aligned = axis_aligned_column(blocks, spacing=3.0)
    write_top(
        aligned,
        block_dir / "ota3_top_aligned.gds",
        block_dir / "ota3_top_aligned.png",
        "ota3_top_aligned",
        "axis-aligned placement (all block centers on one vertical axis)",
    )

    # Also render a couple of representative blocks at higher detail.
    for name in ("input_pair_M0_M7_cc", "mirror_load_M11_M10_cc", "cap_pair_C2_C1_axial"):
        gds_path = block_dir / f"{name}.gds"
        if gds_path.exists():
            render_gds(gds_path, block_dir / f"{name}.png", dpi=200)


if __name__ == "__main__":
    main()
