"""AHS bedrooms-by-kids, normalized WITHIN size bin.

The existing plot normalizes within household type (rows sum to 100). This one
normalizes within bedroom bin (columns sum to 100), so we can see directly: of
all 4+ bedroom units, what fraction are occupied by households with kids?
"""

from __future__ import annotations
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
SRC = ROOT / "code/data/ahs_supply_snapshot/output_ahs_family_unit_menu_national/ahs_occupied_by_kids_and_bedrooms.csv"
OUTDIR = ROOT / "code/data/ahs_supply_snapshot/output_ahs_family_unit_menu_national"


def main():
    units = {}  # (bin, any_kids) -> units
    with open(SRC) as f:
        for row in csv.DictReader(f):
            b = row["bedroom_bin"]
            k = row["any_kids"].strip().upper() == "TRUE"
            u = float(row["units"])
            units[(b, k)] = u

    bins = ["0-1 bedrooms", "2 bedrooms", "3 bedrooms", "4+ bedrooms"]
    kids_share = []
    nokids_share = []
    for b in bins:
        kids = units.get((b, True), 0)
        nokids = units.get((b, False), 0)
        tot = kids + nokids
        kids_share.append(100 * kids / tot if tot > 0 else 0)
        nokids_share.append(100 * nokids / tot if tot > 0 else 0)

    plt.rcParams.update({"font.size": 13, "axes.titlesize": 14, "axes.labelsize": 12})
    fig, ax = plt.subplots(figsize=(7.5, 5.0), constrained_layout=True)
    x = np.arange(len(bins))
    width = 0.62
    p1 = ax.bar(x, kids_share, width, color="#D9534F", label="Household has children")
    p2 = ax.bar(x, nokids_share, width, bottom=kids_share, color="#5276A6",
                label="No children in household")

    # Annotate kids share
    for i, v in enumerate(kids_share):
        ax.text(i, v / 2, f"{v:.0f}%", ha="center", va="center", color="white",
                fontsize=13, fontweight="bold")
    for i, v in enumerate(nokids_share):
        bot = kids_share[i]
        ax.text(i, bot + v / 2, f"{v:.0f}%", ha="center", va="center", color="white",
                fontsize=13, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(bins)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Share of occupied units in size bin")
    ax.set_title("Who occupies each size? (within-bedroom-bin share)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=2, frameon=False)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", alpha=0.3)

    out = OUTDIR / "ahs_bedroom_share_within_size.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    print(f"wrote {out}")
    print("kids share by bin:", dict(zip(bins, [f"{v:.1f}%" for v in kids_share])))


if __name__ == "__main__":
    main()
