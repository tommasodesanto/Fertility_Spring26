"""Cross-metro scatter: periphery-center fertility gap vs periphery-center mean-rooms gap.

The analog of the existing bedroom-scarcity plot using actual room counts (continuous,
interpretable in 'rooms' units) instead of the small/big bedroom-share log ratio.

Reads acs_fertility_gaps_by_metro.csv from output_income_fertility_cross_section
(original MMS sample, 42 large metros) and writes two-panel PNG to the May 24
family-space packet figures folder.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
SRC = ROOT / "code/data/mms_center_periphery/output_income_fertility_cross_section/acs_fertility_gaps_by_metro.csv"
OUTDIR = ROOT / "docs/model/family_space_empirical_packet_2026-05-24/figures"
OUTDIR.mkdir(parents=True, exist_ok=True)


def main():
    rooms_gap, birth_gap, child_gap, weights = [], [], [], []
    with open(SRC) as f:
        for row in csv.DictReader(f):
            try:
                rg = float(row["periphery_center_rooms_gap"])
                bg = float(row["periphery_center_recent_birth_gap"])
                cg = float(row["periphery_center_mean_nchild_gap"])
                w  = float(row["analysis_weight"])
            except (ValueError, KeyError):
                continue
            if any(not np.isfinite(x) for x in [rg, bg, cg, w]):
                continue
            rooms_gap.append(rg); birth_gap.append(bg); child_gap.append(cg); weights.append(w)

    x = np.array(rooms_gap)
    yb = np.array(birth_gap); yc = np.array(child_gap); w = np.array(weights)
    print(f"N metros = {len(x)}")
    print(f"x rooms gap: mean={x.mean():.3f} sd={x.std():.3f} range=[{x.min():.3f}, {x.max():.3f}]")
    print(f"y birth gap: mean={yb.mean():.4f} sd={yb.std():.4f}")
    print(f"y child gap: mean={yc.mean():.4f} sd={yc.std():.4f}")

    # Weighted OLS slopes with proper White SE
    def wls_slope(x, y, w):
        xb = np.average(x, weights=w); yb_ = np.average(y, weights=w)
        sxx = np.sum(w * (x - xb) ** 2)
        b = np.sum(w * (x - xb) * (y - yb_)) / sxx
        a = yb_ - b * xb
        resid = y - (a + b * x)
        # Sandwich SE: meat = sum(w^2 * resid^2 * (x-xb)^2), bread = sxx
        meat = np.sum((w * (x - xb) * resid) ** 2)
        var_b = meat / sxx**2
        return a, b, np.sqrt(var_b)

    a_b, slope_b, se_b = wls_slope(x, yb, w)
    a_c, slope_c, se_c = wls_slope(x, yc, w)
    print(f"birth gap slope on rooms gap: {slope_b:+.5f} (SE {se_b:.5f}) intercept {a_b:+.5f}")
    print(f"child gap slope on rooms gap: {slope_c:+.5f} (SE {se_c:.5f}) intercept {a_c:+.5f}")

    # Marker size by sqrt(weight) for visibility
    smax = 350
    sizes = smax * np.sqrt(w / w.max())

    plt.rcParams.update({"font.size": 12, "axes.titlesize": 12, "axes.labelsize": 11})
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4), constrained_layout=True)

    for ax, y, slope, intercept, se, ylab, ttl in [
        (axes[0], yb, slope_b, a_b, se_b,
         "Periphery $-$ Center recent-birth rate",
         "Recent-Birth Gap"),
        (axes[1], yc, slope_c, a_c, se_c,
         "Periphery $-$ Center mean own children",
         "Mean Own Children Gap"),
    ]:
        ax.scatter(x, y, s=sizes, alpha=0.55, edgecolor="white", linewidth=0.7,
                   color="#5276A6")
        xx = np.linspace(x.min(), x.max(), 50)
        ax.plot(xx, intercept + slope * xx, color="#B23A48", linewidth=2)
        ax.axhline(0, color="gray", linewidth=0.5, linestyle="--", alpha=0.6)
        ax.set_xlabel("Periphery $-$ Center mean rooms (count)")
        ax.set_ylabel(ylab)
        ax.set_title(f"{ttl}  ($\\beta={slope:+.3f}$, robust SE $={se:.3f}$)", fontsize=12)
        ax.grid(True, alpha=0.3)

    out = OUTDIR / "rooms_gap_vs_fertility_gap.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
