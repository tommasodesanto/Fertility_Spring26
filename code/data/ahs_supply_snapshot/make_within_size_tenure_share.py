"""AHS bedrooms-by-tenure, normalized WITHIN size bin.

The existing plot normalizes within tenure category. This one normalizes within
bedroom bin (columns sum to 100), so we see directly: of all 4+ bedroom units,
what fraction are owner-occupied?
"""

from __future__ import annotations
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
SRC = ROOT / "code/data/ahs_supply_snapshot/output_ahs_family_unit_menu_national/ahs_stock_by_bedroom_tenure.csv"
OUTDIR = ROOT / "code/data/ahs_supply_snapshot/output_ahs_family_unit_menu_national"


def main():
    # Aggregate to three groups: owner, renter, vacant.
    # Roll occupied_no_rent and unknown into "other" but ignore for cleanliness, since
    # owner+renter+vacant is the standard split. Keep all categories for full 100%.
    units = {}  # (bin, group) -> units
    keep_groups = {"owner": "owner", "renter": "renter", "vacant": "vacant"}
    with open(SRC) as f:
        for row in csv.DictReader(f):
            t = row["tenure"]
            b = row["bedroom_bin"]
            u = float(row["units"])
            g = keep_groups.get(t)
            if g is None:
                continue
            units[(b, g)] = units.get((b, g), 0) + u

    bins = ["0-1 bedrooms", "2 bedrooms", "3 bedrooms", "4+ bedrooms"]
    groups = ["owner", "renter", "vacant"]
    shares = {g: [] for g in groups}
    for b in bins:
        tot = sum(units.get((b, g), 0) for g in groups)
        for g in groups:
            shares[g].append(100 * units.get((b, g), 0) / tot if tot > 0 else 0)

    plt.rcParams.update({"font.size": 13, "axes.titlesize": 14, "axes.labelsize": 12})
    fig, ax = plt.subplots(figsize=(7.5, 5.0), constrained_layout=True)
    x = np.arange(len(bins))
    width = 0.62
    colors = {"owner": "#E1812C", "renter": "#5276A6", "vacant": "#7B9C76"}
    labels = {"owner": "Owner-occupied", "renter": "Renter-occupied", "vacant": "Vacant"}

    bottom = np.zeros(len(bins))
    for g in groups:
        vals = np.array(shares[g])
        ax.bar(x, vals, width, bottom=bottom, color=colors[g], label=labels[g])
        for i, v in enumerate(vals):
            if v >= 4:
                ax.text(i, bottom[i] + v / 2, f"{v:.0f}%", ha="center", va="center",
                        color="white", fontsize=12, fontweight="bold")
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(bins)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Share of units in size bin")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=3, frameon=False)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", alpha=0.3)

    out = OUTDIR / "ahs_tenure_share_within_size.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    print(f"wrote {out}")
    for g in groups:
        print(f"{g} share by bin:", dict(zip(bins, [f'{v:.1f}%' for v in shares[g]])))


if __name__ == "__main__":
    main()
