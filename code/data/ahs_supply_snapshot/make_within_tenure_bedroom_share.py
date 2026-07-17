"""AHS bedrooms-by-tenure (within-tenure normalization), no title.

Mirrors the existing R ahs_bedroom_menu_by_tenure.png but built fresh in
Python and without a title so it can be used on a stripped slide.
"""

from __future__ import annotations
import csv
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
SRC = ROOT / "code/data/ahs_supply_snapshot/output_ahs_family_unit_menu_national/ahs_stock_by_bedroom_tenure.csv"
OUT = ROOT / "code/data/ahs_supply_snapshot/output_ahs_family_unit_menu_national/ahs_bedroom_menu_by_tenure_notitle.png"


def main():
    cells = {}
    with open(SRC) as f:
        for row in csv.DictReader(f):
            t = row["tenure"]
            b = row["bedroom_bin"]
            try:
                share = float(row["share_within_tenure"]) * 100
            except ValueError:
                continue
            cells[(t, b)] = share

    bins = ["0-1 bedrooms", "2 bedrooms", "3 bedrooms", "4+ bedrooms"]
    keep = [("renter", "Renter", "#5276A6"),
            ("owner",  "Owner",  "#E1812C"),
            ("vacant", "Vacant", "#7B9C76")]

    plt.rcParams.update({"font.size": 13, "axes.labelsize": 12})
    fig, ax = plt.subplots(figsize=(7.5, 5.0), constrained_layout=True)
    x = np.arange(len(bins))
    width = 0.27
    for i, (t, label, color) in enumerate(keep):
        vals = [cells.get((t, b), 0) for b in bins]
        ax.bar(x + (i - 1) * width, vals, width, color=color, label=label)
    ax.set_xticks(x)
    ax.set_xticklabels(bins)
    ax.set_ylabel("Share within tenure")
    ax.set_yticks([0, 10, 20, 30, 40, 50])
    ax.set_yticklabels(["0%", "10%", "20%", "30%", "40%", "50%"])
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=3, frameon=False)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", alpha=0.3)
    fig.savefig(OUT, dpi=120, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
