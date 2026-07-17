"""Plot the Redfin-style horizontal bar chart from the ACS replica CSV."""

from __future__ import annotations
import csv
from pathlib import Path
import matplotlib.pyplot as plt

REPO = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
CSV = REPO / "code/data/mms_center_periphery/output_redfin_big_homes_replica/redfin_big_homes_replica.csv"
OUT = REPO / "code/data/mms_center_periphery/output_redfin_big_homes_replica/redfin_big_homes_replica.png"


def main():
    rows = []
    with open(CSV) as f:
        for r in csv.DictReader(f):
            rows.append(dict(
                generation=r["generation"],
                hh_type=r["hh_type"],
                share=float(r["share"]),
            ))
    rows.sort(key=lambda x: -x["share"])
    labels = [f"{r['generation']}, {r['hh_type']}" for r in rows]
    shares = [r["share"] * 100 for r in rows]

    # Color by hh_type
    color_map = {
        "1-2 adults": "#5276A6",
        "with own children": "#E1812C",
        "3+ people, no children": "#7B9C76",
    }
    colors = [color_map.get(r["hh_type"], "#888888") for r in rows]

    plt.rcParams.update({"font.size": 12, "axes.titlesize": 14, "axes.labelsize": 12})
    fig, ax = plt.subplots(figsize=(11, 7.5), constrained_layout=True)
    y = list(range(len(labels)))
    bars = ax.barh(y, shares, color=colors, edgecolor="white", linewidth=0.6)
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.invert_yaxis()  # largest on top
    ax.set_xlabel("Share of 3+ bedroom owner-occupied homes (%)")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.grid(True, axis="x", alpha=0.3)
    for i, v in enumerate(shares):
        ax.text(v + 0.4, i, f"{v:.1f}%", va="center", fontsize=11)
    ax.set_xlim(0, max(shares) * 1.18)

    # Legend
    handles = [plt.Rectangle((0,0), 1, 1, color=c) for c in color_map.values()]
    ax.legend(handles, list(color_map.keys()), loc="lower right", frameon=False)

    fig.savefig(OUT, dpi=120, bbox_inches="tight")
    print(f"wrote {OUT}")
    # Also print summary aggregates
    by_type = {}
    for r in rows:
        by_type[r["hh_type"]] = by_type.get(r["hh_type"], 0) + r["share"]
    print("\nAggregated by household type:")
    for k, v in sorted(by_type.items(), key=lambda kv: -kv[1]):
        print(f"  {k}: {v*100:.1f}%")
    print(f"\nBoomer + Silent (60+) '1-2 adults' = {(by_type.get('1-2 adults', 0) - sum(r['share'] for r in rows if r['hh_type']=='1-2 adults' and r['generation'] in ('Gen X','Millennial','Gen Z')))*100:.1f}%")
    empty_nest_old = sum(r["share"] for r in rows if r["hh_type"] == "1-2 adults" and r["generation"] in ("Boomer", "Silent"))
    print(f"Empty-nester (Boomer + Silent, 1-2 adults) share = {empty_nest_old*100:.1f}%")


if __name__ == "__main__":
    main()
