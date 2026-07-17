"""Cross-metro scatter: log rent gap (center/periphery) vs.
non-parent minus new-parent center destination share among movers.

Matches the style of rooms_gap_vs_fertility_gap.png (scatter, WLS line, beta + robust SE).
"""

from __future__ import annotations
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
SRC = ROOT / "output/model/acs_income_crossmetro_v1/acs_cross_metro_dest_gap_v1.csv"
OUT = ROOT / "output/model/acs_supporting_facts_v1/acs_cross_metro_dest_gap_scatter.png"


def main():
    x_log, y, w = [], [], []
    with open(SRC) as f:
        for row in csv.DictReader(f):
            try:
                lr = float(row["log_rent_gap"])
                gp = float(row["gap_pp"])
                mw = float(row["mover_obs"])
            except (ValueError, KeyError):
                continue
            if not all(np.isfinite([lr, gp, mw])) or mw <= 0:
                continue
            x_log.append(lr); y.append(gp); w.append(mw)

    x = np.array(x_log); y = np.array(y); w = np.array(w)
    print(f"N metros = {len(x)}")
    print(f"x log-rent gap: mean={x.mean():.3f} sd={x.std():.3f}")
    print(f"y gap_pp:        mean={y.mean():.3f} sd={y.std():.3f}")

    def wls(x, y, w):
        xb = np.average(x, weights=w); yb = np.average(y, weights=w)
        sxx = np.sum(w * (x - xb) ** 2)
        b = np.sum(w * (x - xb) * (y - yb)) / sxx
        a = yb - b * xb
        resid = y - (a + b * x)
        var_b = np.sum((w * (x - xb) * resid) ** 2) / sxx ** 2
        return a, b, np.sqrt(var_b)

    a, slope, se = wls(x, y, w)
    print(f"slope = {slope:+.3f}  robust SE = {se:.3f}  intercept = {a:+.3f}")

    sizes = 350 * np.sqrt(w / w.max())
    plt.rcParams.update({"font.size": 12, "axes.titlesize": 12, "axes.labelsize": 11})
    fig, ax = plt.subplots(figsize=(6.7, 5.4), constrained_layout=True)
    ax.scatter(x, y, s=sizes, alpha=0.55, edgecolor="white", linewidth=0.7, color="#5276A6")
    xx = np.linspace(x.min(), x.max(), 50)
    ax.plot(xx, a + slope * xx, color="#B23A48", linewidth=2)
    ax.axhline(0, color="gray", linewidth=0.5, linestyle="--", alpha=0.6)
    ax.set_xlabel(r"$\log(r_C / r_P)$ (metro rent gap)")
    ax.set_ylabel(r"Non-parent $-$ new-parent center destination (pp)")
    ax.set_title(rf"Center Sorting Gap  ($\beta={slope:+.2f}$, robust SE $={se:.2f}$)", fontsize=12)
    ax.grid(True, alpha=0.3)
    fig.savefig(OUT, dpi=180, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
