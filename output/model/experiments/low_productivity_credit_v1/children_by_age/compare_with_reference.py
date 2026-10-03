"""Compare age-specific children-ever-born shares in saved reference and low-z runs."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
REFERENCE = ROOT / "output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/fable_analysis/credit_mechanism/credit_relaxation/children_by_age/children_by_age.csv"
LOW = HERE / "figures/children_by_age.csv"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_rows(path: Path, regime: str) -> list[dict]:
    with path.open(newline="") as stream:
        raw = list(csv.DictReader(stream))
    out = []
    for row in raw:
        item = {"regime": regime, "phi": float(row["arm"].split("_")[1]) / 100,
                "age_start": float(row["age_start"]), "age_end_exclusive": float(row["age_end_exclusive"]),
                "household_mass": float(row["household_mass"])}
        for count in ("0", "1", "2", "3plus"):
            item[f"share_n_{count}"] = float(row[f"share_n_{count}"])
        item["share_n_2plus"] = item["share_n_2"] + item["share_n_3plus"]
        out.append(item)
    return out


def main() -> None:
    rows = read_rows(REFERENCE, "heterogeneous_productivity") + read_rows(LOW, "permanent_lowest_productivity")
    by_key = {(r["regime"], r["phi"], r["age_start"]): r for r in rows}
    ages = sorted({r["age_start"] for r in rows})
    if len(ages) != 17 or len(by_key) != 4 * len(ages):
        raise ValueError("Expected four complete arms on the same 17 lifecycle ages")
    for r in rows:
        if r["age_end_exclusive"] != r["age_start"] + 4:
            raise ValueError("Mismatched four-year age cells")
        if abs(sum(r[f"share_n_{k}"] for k in ("0", "1", "2", "3plus")) - 1) > 2e-14:
            raise ValueError("Count shares do not sum to one")

    path = HERE / "comparison_by_age.csv"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7), sharex=True, layout="constrained")
    for regime, color, name in (("heterogeneous_productivity", "#456990", "Reference productivity"),
                                ("permanent_lowest_productivity", "#c24b3a", "Everyone at lowest productivity")):
        for phi, style in ((0.8, "-"), (1.0, "--")):
            subset = [by_key[(regime, phi, age)] for age in ages]
            for ax, col in ((axes[0], "share_n_0"), (axes[1], "share_n_2plus")):
                ax.plot(ages, [r[col] for r in subset], style, color=color, linewidth=2.3,
                        label=rf"{name}, $\phi={phi:.1f}$")
    axes[0].set_title("No children ever born")
    axes[1].set_title("Two or more children ever born")
    for ax in axes:
        ax.set_xlabel("Age cell start (years)")
        ax.set_ylabel("Share of all households at this age")
        ax.set_ylim(-0.02, 1.02)
        ax.grid(alpha=0.2)
    axes[1].legend(frameon=False, loc="center left", fontsize=9)
    fig.suptitle("Children ever born by household age: reference versus permanently low productivity\n"
                 "Fixed price; four-year age cells; the two low-productivity lines overlap", fontsize=12)
    fig.savefig(HERE / "reference_vs_low_productivity_by_age.png", dpi=180)
    plt.close(fig)

    low_parent_max = max(1 - by_key[("permanent_lowest_productivity", phi, age)]["share_n_0"]
                         for phi in (0.8, 1.0) for age in ages)
    receipt = {"reference_csv": str(REFERENCE), "reference_csv_sha256": sha256(REFERENCE),
               "low_csv": str(LOW), "low_csv_sha256": sha256(LOW),
               "max_low_productivity_share_with_any_children": low_parent_max,
               "source": "Saved post-birth g distributions; no model solve or policy interpolation",
               "measure": "Within each arm and four-year age cell, all households including childless"}
    (HERE / "comparison_provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()
