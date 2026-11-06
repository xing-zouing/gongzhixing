# Matched layout final results

This document records the current final flow: netlist + matching constraints
-> unit-cell decomposition -> common-centroid algorithm -> horizontal
centroid-preserving pair-move compaction -> GDS and two preview images.

The older macro/global-axis experiment is superseded by this flow.

## Flow

1. Parse the Spectre netlist and the `.sym` matching groups.
2. For every two-device matching group:
   - if both devices decompose into one unit, merge them into one axial-pair
     label (count 2);
   - otherwise give each device its own label with its unit count.
3. Ignore unmatched single-device groups for now.
4. Feed all labels to `cc_main_flow`; require `TMOC = 0`.
5. Split merged labels back into the two physical devices.
6. Write a heterogeneous table GDS.
7. Compact with centroid-preserving cell-pair moves:
   - every move is `+d` on one cell and `-d` on another cell of the same
     matching set;
   - the device centroid / axial-pair midpoint is therefore unchanged;
   - `min_gap = 0.3 um` is a hard constraint.
8. Emit:
   - `<case>_matched.gds`
   - `<case>_matched.png` (GDS render)
   - `<case>_matched_groups.png` (one color per matching group, dummies gray)

Horizontal orientation is retained.  Vertical rotation was tested but not
kept: it reduced area only because it changed the aspect ratio and the
pair-move iterations, not because rotation itself helped matching.  Square
reservation (`max(w,h)^2` per cell) was also tested and rejected; it
over-reserves wide cells and gives a worse final area.

## Automatic unit split

For each MOS device the flow starts with the explicit `units` override, or
with the largest divisor of the total finger count up to `--max-units`
(default 8).  It then checks the unit-cell aspect

```text
aspect = max(cell_w, cell_h) / min(cell_w, cell_h)
```

If `aspect > --unit-aspect-threshold` (default 3.5), it tries splits with
more units, in increasing unit count, and takes the first one that passes.
If no split passes, it keeps the initial split.

Example for ota3 M0/M7:

```text
4 units of 10 fingers : 19.05 x 3.65, aspect 5.22  -> fails
8 units of 5 fingers  : 11.85 x 3.65, aspect 3.25  -> accepted
```

The threshold is deliberately 3.5: 3.0 would force the 5-finger unit
(aspect 3.25) to split further into many 1-finger units, while 4.0 is too
permissive for elongated cells.

## Results

All three cases reach `TMOC = 0`, preserve the physical centroids to floating
point precision, and keep a 0.3 um minimum clearance.

| case | subckt | groups | counts | matrix | TMOC | bbox (um) | area (um^2) | utilization | set centroid err (um) | group-center err (um) |
|---|---|---:|---|---:|---:|---|---:|---:|---:|---:|
| ota1 | OPA_my | 6 | [4, 4, 2, 2, 2, 2, 2] | 5x5 | 0 | 33.22 x 21.00 | 697.8 | 66.3% | 6.27e-15 | 5.17e-15 |
| ota2 | ota2 | 3 | [4, 4, 4, 4, 2] | 5x5 | 0 | 28.95 x 17.95 | 519.7 | 74.3% | 0.00e+00 | 2.39e-16 |
| ota3 | OTA2_my | 6 | [8, 8, 2, 2, 2, 2, 2, 2] | 6x6 | 0 | 70.22 x 27.73 | 1947.0 | 78.1% | 2.29e-15 | 2.26e-15 |

Notes:

- ota1 groups: M13/M14 use 4 units each; all other matched pairs are
  single-unit axial pairs.  C0 and R0 are unmatched and ignored.
- ota2 groups: M8/M9 and M1/M0 use 4 units each; M10/M11 is an axial pair.
  C0 and M2 are unmatched and ignored.
- ota3 groups: M0/M7 use 8 units of 5 fingers; M11/M10 use 2 units of 3
  fingers; the remaining four pairs are axial pairs.  M9/M12 are unmatched
  and ignored.

## Outputs

```text
build/generic/ota1/ota1_matched.gds
build/generic/ota1/ota1_matched.png
build/generic/ota1/ota1_matched_groups.png

build/generic/ota2/ota2_matched.gds
build/generic/ota2/ota2_matched.png
build/generic/ota2/ota2_matched_groups.png

build/generic/ota3/ota3_matched.gds
build/generic/ota3/ota3_matched.png
build/generic/ota3/ota3_matched_groups.png

build/generic/summary.md
```

## Configure and run

The case table is at the top of `scripts/run_matched_layout.py`:

```python
CASES = {
    "ota3": {
        "netlist": ROOT / "examples" / "ota3.sp",
        "sym": ROOT / "examples" / "ota3.sym",
        "subckt": "OTA2_my",
        "units": {"M0": 8, "M7": 8, "M11": 2, "M10": 2},
    },
    "ota1": {
        "netlist": ROOT / "mycase" / "ota1.sp",
        "sym": ROOT / "mycase" / "ota1.sym",
        "subckt": "OPA_my",
        "units": {},
    },
    "ota2": {
        "netlist": ROOT / "mycase" / "ota2.sp",
        "sym": ROOT / "mycase" / "ota2.sym",
        "subckt": "ota2",
        "units": {},
    },
}
```

Add or edit a case, then run:

```bash
# all configured cases
.venv/bin/python scripts/run_matched_layout.py --case all

# one case
.venv/bin/python scripts/run_matched_layout.py --case ota1

# larger pair-move search
.venv/bin/python scripts/run_matched_layout.py --case ota3 \
  --iterations 500000 --seed 0

# change or disable the automatic unit split
.venv/bin/python scripts/run_matched_layout.py --case ota3 \
  --unit-aspect-threshold 3.5 --max-units 8
# --unit-aspect-threshold 0 disables the split check
```

## Caveats

- The GDS uses the mock TSMC40-style layer stack from
  `cc/device_generation`; it is a placement prototype, not a PDK/DRC/LVS
  clean tape-out layout.
- Routing is not generated.
- Unmatched devices are intentionally excluded in this phase.
- The pair-move compaction preserves the centroid and minimum spacing but
  may change the original unit-cell pattern; the amount of cell displacement
  is reported by the flow.
