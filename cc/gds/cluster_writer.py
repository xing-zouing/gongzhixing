"""Write an automatically clustered shared-array GDS block."""

from __future__ import annotations

from pathlib import Path

from cc.gds.array_writer import generate_group_gds
from cc.io.netlist_parser import Device
from cc.place.cluster_search import build_cluster_layout_best
from cc.place.compatibility import cluster_slug


def write_cluster_gds(
    devices: dict[str, Device],
    cluster: dict,
    output_dir: str | Path,
) -> Path | None:
    """Build and write the shared-array GDS for one compatible cluster.

    Returns the GDS path, or ``None`` if the cluster's counts cannot be
    represented by the orbit allocator (caller should fall back).
    """
    result = build_cluster_layout_best(devices, cluster)
    if result is None:
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    slug = cluster_slug(cluster)
    output_path = output_dir / f"{slug}.gds"

    generate_group_gds(
        devices=result["gds_devices"],
        cc_matrix=result["matrix"],
        output_path=output_path,
        array_name=slug,
    )
    return output_path
