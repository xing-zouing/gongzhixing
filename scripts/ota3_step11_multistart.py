#!/usr/bin/env python3
"""Step 11: multi-start + Pareto search for the unit-level layout.

Runs Step 10 for several initialization modes and random seeds, collects
``(top_area, |aspect-1|, hpwl)``, and keeps the non-dominated solutions.
The Pareto GDS/PNG files are copied to ``build/ota3/pareto/``.

Run from the project root:

    .venv/bin/python scripts/ota3_step11_multistart.py
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
STEP10 = ROOT / "scripts" / "ota3_step10_unit_level.py"
BASE = ROOT / "build" / "ota3" / "multistart"
PARETO_DIR = ROOT / "build" / "ota3" / "pareto"


def run_one(mode: str, seed: int, iterations: int) -> dict:
    run_dir = BASE / f"{mode}_s{seed}"
    run_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(STEP10),
        "--seed", str(seed),
        "--init", mode,
        "--iterations", str(iterations),
        "--out-dir", str(run_dir),
        "--name", "unit",
    ]
    subprocess.run(
        cmd,
        cwd=ROOT,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
    )
    return json.loads((run_dir / "summary.json").read_text())


def dominates(a: dict, b: dict) -> bool:
    keys = ("top_area", "aspect_error", "hpwl")
    return all(a[key] <= b[key] + 1e-12 for key in keys) and any(
        a[key] < b[key] - 1e-12 for key in keys
    )


def main() -> None:
    modes = ("checkerboard", "random", "central")
    seeds = (0, 1)
    iterations = 20000

    rows = []
    for mode in modes:
        for seed in seeds:
            print(f"running {mode} seed={seed} ...", flush=True)
            summary = run_one(mode, seed, iterations)
            summary["aspect_error"] = abs(summary["aspect"] - 1.0)
            rows.append(summary)

    pareto = []
    for row in rows:
        if not any(dominates(other, row) for other in rows if other is not row):
            pareto.append(row)

    pareto.sort(key=lambda item: item["top_area"])
    PARETO_DIR.mkdir(parents=True, exist_ok=True)
    for old in PARETO_DIR.glob("pareto_*"):
        old.unlink()

    for index, row in enumerate(pareto):
        stem = f"pareto_{index:02d}_{row['init']}_s{row['seed']}"
        for key, suffix in (("gds", ".gds"), ("png", ".png"), ("groups_png", "_groups.png")):
            source = Path(row[key])
            if source.exists():
                shutil.copy2(source, PARETO_DIR / f"{stem}{suffix}")
        row["pareto_stem"] = stem

    print("\n=== all runs ===")
    print(f"{'init':12s} {'seed':>4s} {'top_area':>10s} {'WxH':>16s} {'aspect':>7s} {'hpwl':>8s} {'pareto':>6s}")
    for row in sorted(rows, key=lambda item: item["top_area"]):
        is_pareto = row in pareto
        print(
            f"{row['init']:12s} {row['seed']:4d} {row['top_area']:10.1f} "
            f"{row['top_width']:.2f}x{row['top_height']:.2f} "
            f"{row['aspect']:7.3f} {row['hpwl']:8.2f} {str(is_pareto):>6s}"
        )

    csv_path = PARETO_DIR / "pareto.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "init", "seed", "top_area", "top_width", "top_height",
                "aspect", "hpwl", "overlap", "area", "width", "height",
                "pareto_stem",
            ],
        )
        writer.writeheader()
        for row in pareto:
            writer.writerow({key: row.get(key) for key in writer.fieldnames})

    fig, ax = plt.subplots(figsize=(7, 5))
    for row in rows:
        is_pareto = row in pareto
        ax.scatter(
            row["top_area"], row["aspect_error"],
            s=60 if is_pareto else 30,
            c="#d62728" if is_pareto else "#999999",
            marker="o",
            zorder=3 if is_pareto else 2,
        )
    for row in pareto:
        ax.annotate(
            f"{row['init']}-s{row['seed']}",
            (row["top_area"], row["aspect_error"]),
            fontsize=7,
            xytext=(4, 4),
            textcoords="offset points",
        )
    ax.set_xlabel("top area (um^2)")
    ax.set_ylabel("|aspect - 1|")
    ax.set_title("multi-start results (red = Pareto)")
    fig.tight_layout()
    plot_path = PARETO_DIR / "pareto.png"
    fig.savefig(plot_path, dpi=160)
    plt.close(fig)

    print(f"\nPareto solutions: {len(pareto)}")
    print(f"CSV: {csv_path}")
    print(f"plot: {plot_path}")
    for row in pareto:
        print(
            f"  {row['pareto_stem']}: {row['top_area']:.1f} um^2, "
            f"{row['top_width']:.2f}x{row['top_height']:.2f}, "
            f"aspect={row['aspect']:.3f}, hpwl={row['hpwl']:.2f}"
        )


if __name__ == "__main__":
    main()
