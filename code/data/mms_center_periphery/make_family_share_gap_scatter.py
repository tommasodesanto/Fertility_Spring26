"""Cross-metro scatter: periphery-center fertility gap vs. periphery-center
share-of-3+-bedroom-units gap.

The clean availability measure. Says directly: "the periphery has X pp more
3+ bedroom units than the center."

Reads two ACS-derived CSVs and produces one two-panel figure.
"""

from __future__ import annotations
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
SUPPLY = ROOT / "code/data/mms_center_periphery/output_couillard_bedroom_supply/acs_bedroom_premia_by_metro.csv"
FERT   = ROOT / "code/data/mms_center_periphery/output_income_fertility_cross_section/acs_fertility_gaps_by_metro.csv"
OUT    = ROOT / "docs/model/family_space_empirical_packet_2026-05-24/figures/family_share_gap_vs_fertility_gap.png"


def main():
    supply = {}
    with open(SUPPLY) as f:
        for row in csv.DictReader(f):
            try:
                m = row["met2013"]
                cs = float(row["family_stock_share_center"])
                ps = float(row["family_stock_share_periphery"])
                # gap in percentage points
                supply[m] = 100 * (ps - cs)
            except (ValueError, KeyError):
                continue

    rows = []
    with open(FERT) as f:
        for row in csv.DictReader(f):
            try:
                m = row["met2013"]
                if m not in supply: continue
                rows.append({
                    "m": m,
                    "x": supply[m],
                    "yb": float(row["periphery_center_recent_birth_gap"]),
                    "yc": float(row["periphery_center_mean_nchild_gap"]),
                    "w": float(row["analysis_weight"]),
                })
            except (ValueError, KeyError):
                continue

    print(f"matched metros: {len(rows)}")
    x  = np.array([r["x"]  for r in rows])
    yb = np.array([r["yb"] for r in rows])
    yc = np.array([r["yc"] for r in rows])
    w  = np.array([r["w"]  for r in rows])

    print(f"x family-share gap (pp): mean={x.mean():.2f} sd={x.std():.2f} range=[{x.min():.2f}, {x.max():.2f}]")

    def wls(x, y, w):
        xb = np.average(x, weights=w); yb_ = np.average(y, weights=w)
        sxx = np.sum(w * (x - xb)**2)
        b = np.sum(w * (x - xb) * (y - yb_)) / sxx
        a = yb_ - b * xb
        resid = y - (a + b * x)
        var_b = np.sum((w * (x - xb) * resid)**2) / sxx**2
        return a, b, np.sqrt(var_b)

    ab, bb, seb = wls(x, yb, w)
    ac, bc, sec = wls(x, yc, w)
    print(f"birth gap slope: {bb:+.5f} (SE {seb:.5f})")
    print(f"child gap slope: {bc:+.5f} (SE {sec:.5f})")

    sizes = 350 * np.sqrt(w / w.max())

    plt.rcParams.update({"font.size": 13, "axes.titlesize": 13, "axes.labelsize": 12})
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4), constrained_layout=True)
    for ax, y, slope, intercept, se, ylab, ttl in [
        (axes[0], yb, bb, ab, seb,
         "Periphery $-$ Center recent-birth rate",
         "Recent-Birth Gap"),
        (axes[1], yc, bc, ac, sec,
         "Periphery $-$ Center mean own children",
         "Mean Own Children Gap"),
    ]:
        ax.scatter(x, y, s=sizes, alpha=0.55, edgecolor="white", linewidth=0.7,
                   color="#5276A6")
        xx = np.linspace(x.min(), x.max(), 50)
        ax.plot(xx, intercept + slope * xx, color="#B23A48", linewidth=2)
        ax.axhline(0, color="gray", linewidth=0.5, linestyle="--", alpha=0.6)
        ax.axvline(0, color="gray", linewidth=0.5, linestyle="--", alpha=0.6)
        ax.set_xlabel("Share of 3+ bedroom units: periphery $-$ center (pp)")
        ax.set_ylabel(ylab)
        ax.set_title(f"{ttl}  ($\\beta={slope:+.4f}$, robust SE $={se:.4f}$)")
        ax.grid(True, alpha=0.3)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=180, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
