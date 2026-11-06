"""Automatic compatibility clustering for matching groups.

Groups can share one interleaved unit-cell array only when their unit cells
are physically identical: same device kind, same polarity/model attributes,
same unit width/length, same finger count, etc.  Groups that do not match
stay as standalone blocks.

The six ota3 groups therefore split into:

* one cluster: M2/M5 + M4/M3 (both pch, 2.5u/5u, nf=1);
* four standalone singletons (M0/M7, M8/M6, M11/M10, C2/C1).

This is a physical constraint, not a limitation of the orbit allocator.
"""

from __future__ import annotations

from typing import Iterable

from cc.io.netlist_parser import Device
from cc.place.interleave import Group, allocate_orbits_1d


def unit_key(device: Device) -> tuple:
    """A comparable key describing the physical unit cell."""
    attrs = tuple(sorted(device.attr))
    if device.kind == "mos":
        return (
            "mos",
            bool(device.nch),
            round(device.unit_width, 6),
            round(device.unit_length, 6),
            int(device.nf),
            attrs,
        )
    if device.kind == "cap":
        return (
            "cap",
            round(device.unit_width, 6),
            round(device.unit_length, 6),
            int(device.nr),
            int(device.params.get("stm", 3)),
            int(device.params.get("spm", 5)),
            attrs,
        )
    if device.kind == "res":
        return (
            "res",
            round(device.unit_width, 6),
            round(device.unit_length, 6),
            int(device.series),
            attrs,
        )
    return (device.kind,)


def group_unit_counts(devices: dict[str, Device], spec: dict) -> tuple[int, ...]:
    """Unit-cell count for each device in a group."""
    devs = [devices[name] for name in spec["devices"]]
    if spec["type"] == "axial_cap":
        return tuple(dev.multi for dev in devs)
    seg_nf = int(spec["seg_nf"])
    counts = []
    for dev in devs:
        total_fingers = dev.multi * dev.nf
        if total_fingers % seg_nf != 0:
            raise ValueError(
                f"{spec['name']}: {dev.name} total units {total_fingers} "
                f"not divisible by seg_nf={seg_nf}"
            )
        counts.append(total_fingers // seg_nf)
    return tuple(counts)


def cluster_groups(
    devices: dict[str, Device],
    group_specs: Iterable[dict],
) -> list[dict]:
    """Cluster group specs by physical unit-cell compatibility.

    Returns a list of ``{"key", "groups", "interleavable"}`` dicts, in the
    order the first member group appears in ``group_specs``.
    """
    clusters: list[dict] = []
    key_to_index: dict[tuple, int] = {}

    for spec in group_specs:
        devs = [devices[name] for name in spec["devices"]]

        eligible = (
            len(devs) == 2
            and all(dev.kind == "mos" for dev in devs)
            and spec.get("type") in ("axial", "common_centroid")
        )
        key = None
        if eligible:
            keys = [unit_key(dev) for dev in devs]
            if all(candidate == keys[0] for candidate in keys):
                key = keys[0]

        if key is None:
            clusters.append({"key": None, "groups": [spec], "interleavable": False})
            continue

        if key in key_to_index:
            clusters[key_to_index[key]]["groups"].append(spec)
        else:
            key_to_index[key] = len(clusters)
            clusters.append({"key": key, "groups": [spec], "interleavable": False})

    for cluster in clusters:
        cluster["interleavable"] = (
            len(cluster["groups"]) >= 2
            and all(len(spec["devices"]) == 2 for spec in cluster["groups"])
        )
    return clusters


def cluster_slug(cluster: dict) -> str:
    return "cluster_" + "__".join(spec["name"] for spec in cluster["groups"])


def build_cluster_layout(devices: dict[str, Device], cluster: dict) -> dict | None:
    """Build a 1-D shared-array layout for a compatible cluster.

    Returns ``None`` if the counts cannot be represented by reflection
    orbits (e.g. two groups with odd total counts competing for one center
    cell).  In that case the caller should fall back to standalone blocks.
    """
    group_objects: list[Group] = []
    gds_devices: list[dict] = []
    seen: set[str] = set()

    for spec in cluster["groups"]:
        counts = group_unit_counts(devices, spec)
        group_objects.append(
            Group(name=spec["name"], devices=tuple(spec["devices"]), counts=counts)
        )
        seg_nf = spec.get("seg_nf")
        for name in spec["devices"]:
            if name in seen:
                continue
            seen.add(name)
            dev = devices[name]
            gds_devices.append(
                {
                    "name": dev.name,
                    "nch": dev.nch,
                    "model": dev.model,
                    "w": dev.unit_width,
                    "l": dev.unit_length,
                    "nf": int(seg_nf) if seg_nf is not None else dev.nf,
                    "multi": dev.multi,
                    "attr": dev.attr,
                }
            )

    total_cells = sum(sum(group.counts) for group in group_objects)
    if total_cells < 1:
        return None

    # Odd total counts need a center cell.  Try the natural length first,
    # then add one dummy/center cell if needed.
    lengths = [total_cells]
    if total_cells % 2 == 0:
        lengths.append(total_cells + 1)

    layout = None
    used_length = None
    for length in lengths:
        try:
            layout = allocate_orbits_1d(group_objects, length)
            used_length = length
            break
        except ValueError:
            continue
    if layout is None or used_length is None:
        return None

    name_to_id = {device["name"]: index + 1 for index, device in enumerate(gds_devices)}
    matrix = [
        [0 if device is None else name_to_id[device] for device in layout]
    ]
    return {
        "gds_devices": gds_devices,
        "matrix": matrix,
        "layout": layout,
        "groups": group_objects,
        "length": used_length,
    }
