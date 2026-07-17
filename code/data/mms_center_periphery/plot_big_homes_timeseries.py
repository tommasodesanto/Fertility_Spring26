"""Plot the time series of intergenerational allocation of large homes.

Three plots:
  (1) Share of 3+ bed owner-occupied homes by household type, 2012-2023
      (1-2 adults vs with-minor-kids vs 3+ no kids)
  (2) Within 'with minor kids', shrinking? rising?
  (3) Within '1-2 adults', age composition over time (is the
      empty-nester share being driven by 60+ heads, or by 40-59 heads?)
  (4) Fixed-cohort generation shares over time
"""

from __future__ import annotations
import csv
from collections import defaultdict
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

REPO = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
OUTDIR = REPO / "code/data/mms_center_periphery/output_big_homes_timeseries"

HHTYPE_COLORS = {
    "1-2 adults": "#5276A6",
    "with own children": "#E1812C",
    "3+ people, no children": "#7B9C76",
}
AGE_COLORS = {"22-39": "#5276A6", "40-59": "#E1812C", "60+": "#C73E3A"}
GEN_COLORS = {
    "Silent": "#5C5C5C",
    "Boomer": "#C73E3A",
    "Gen X": "#E1812C",
    "Millennial": "#5276A6",
    "Gen Z": "#7B9C76",
}


def load_csv(path):
    rows = []
    with open(path) as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


def plot_hhtype():
    rows = load_csv(OUTDIR / "big_homes_hhtype_by_year.csv")
    by_yt = defaultdict(dict)
    for r in rows:
        by_yt[int(r["year"])][r["hh_type"]] = float(r["share"])
    years = sorted(by_yt.keys())
    fig, ax = plt.subplots(figsize=(10, 5.2), constrained_layout=True)
    for t in ["1-2 adults", "with own children", "3+ people, no children"]:
        vals = [by_yt[y].get(t, np.nan) * 100 for y in years]
        ax.plot(years, vals, marker="o", linewidth=2.5, markersize=7,
                color=HHTYPE_COLORS.get(t, "gray"), label=t)
    ax.set_xlabel("Year"); ax.set_ylabel("Share of 3+ bed owner-occupied homes (%)")
    ax.grid(True, alpha=0.3); ax.legend(loc="best", frameon=False)
    out = OUTDIR / "ts_hhtype.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print(f"wrote {out}")
    # Print summary
    print("\nHousehold-type composition over time (% of 3+ bed owner-occ):")
    print(f"{'year':>6} {'1-2 adults':>15} {'with kids':>15} {'3+ no kids':>15}")
    for y in years:
        print(f"{y:>6} {by_yt[y].get('1-2 adults',0)*100:>14.1f}% {by_yt[y].get('with own children',0)*100:>14.1f}% {by_yt[y].get('3+ people, no children',0)*100:>14.1f}%")


def plot_emptynester_age():
    rows = load_csv(OUTDIR / "big_homes_emptynester_by_age_year.csv")
    by_ya = defaultdict(dict)
    for r in rows:
        by_ya[int(r["year"])][r["age_bin"]] = float(r["share_within_1_2adults"])
    years = sorted(by_ya.keys())
    fig, ax = plt.subplots(figsize=(10, 5.2), constrained_layout=True)
    for age in ["22-39", "40-59", "60+"]:
        vals = [by_ya[y].get(age, np.nan) * 100 for y in years]
        ax.plot(years, vals, marker="o", linewidth=2.5, markersize=7,
                color=AGE_COLORS.get(age, "gray"), label=age)
    ax.set_xlabel("Year"); ax.set_ylabel("Share of '1-2 adults' big-home owners (%)")
    ax.grid(True, alpha=0.3); ax.legend(title="Head age", loc="best", frameon=False)
    out = OUTDIR / "ts_emptynester_age.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print(f"wrote {out}")
    print("\nWithin '1-2 adults', age composition over time:")
    print(f"{'year':>6} {'22-39':>10} {'40-59':>10} {'60+':>10}")
    for y in years:
        print(f"{y:>6} {by_ya[y].get('22-39',0)*100:>9.1f}% {by_ya[y].get('40-59',0)*100:>9.1f}% {by_ya[y].get('60+',0)*100:>9.1f}%")


def plot_generation():
    rows = load_csv(OUTDIR / "big_homes_generation_by_year.csv")
    # Aggregate over hh_type within (year, generation) -> total share by generation
    by_yg = defaultdict(lambda: defaultdict(float))
    for r in rows:
        by_yg[int(r["year"])][r["generation"]] += float(r["share"])
    years = sorted(by_yg.keys())
    fig, ax = plt.subplots(figsize=(10, 5.2), constrained_layout=True)
    for gen in ["Silent", "Boomer", "Gen X", "Millennial", "Gen Z"]:
        vals = [by_yg[y].get(gen, np.nan) * 100 for y in years]
        ax.plot(years, vals, marker="o", linewidth=2.5, markersize=7,
                color=GEN_COLORS.get(gen, "gray"), label=gen)
    ax.set_xlabel("Year"); ax.set_ylabel("Share of 3+ bed owner-occupied homes (%)")
    ax.grid(True, alpha=0.3); ax.legend(title="Generation", loc="best", frameon=False)
    out = OUTDIR / "ts_generation_share.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print(f"wrote {out}")
    print("\nFixed-cohort generation share of 3+ bed owner-occ homes over time:")
    print(f"{'year':>6} " + " ".join(f"{g:>11}" for g in ["Silent","Boomer","Gen X","Millennial","Gen Z"]))
    for y in years:
        print(f"{y:>6} " + " ".join(f"{by_yg[y].get(g,0)*100:>10.1f}%" for g in ["Silent","Boomer","Gen X","Millennial","Gen Z"]))


def plot_boomer_emptynester():
    """Track the key Redfin cell: Boomer + 1-2 adults, share over time."""
    rows = load_csv(OUTDIR / "big_homes_generation_by_year.csv")
    cells = {}
    for r in rows:
        key = (int(r["year"]), r["generation"], r["hh_type"])
        cells[key] = float(r["share"])
    years = sorted(set(int(r["year"]) for r in rows))
    fig, ax = plt.subplots(figsize=(10, 5.2), constrained_layout=True)
    headlines = [
        ("Boomer + Silent, 1-2 adults", ["Boomer","Silent"], "1-2 adults", "#C73E3A"),
        ("Boomer + Silent, with minor kids", ["Boomer","Silent"], "with own children", "#7B9C76"),
        ("Millennial + Gen Z, 1-2 adults", ["Millennial","Gen Z"], "1-2 adults", "#1f5fa6"),
        ("Millennial + Gen Z, with minor kids", ["Millennial","Gen Z"], "with own children", "#E1812C"),
    ]
    for label, gens, typ, color in headlines:
        ys = []
        for y in years:
            total = sum(cells.get((y, g, typ), 0) for g in gens)
            ys.append(total * 100)
        ax.plot(years, ys, marker="o", linewidth=2.5, markersize=7, color=color, label=label)
    ax.set_xlabel("Year"); ax.set_ylabel("Share of 3+ bed owner-occupied homes (%)")
    ax.grid(True, alpha=0.3); ax.legend(loc="best", fontsize=10, frameon=False)
    out = OUTDIR / "ts_headline_cells.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print(f"wrote {out}")
    print("\nHeadline cells over time:")
    for label, gens, typ, color in headlines:
        print(f"  {label}:")
        for y in years:
            total = sum(cells.get((y, g, typ), 0) for g in gens)
            print(f"    {y}: {total*100:.1f}%")


def main():
    plot_hhtype()
    plot_emptynester_age()
    plot_generation()
    plot_boomer_emptynester()


if __name__ == "__main__":
    main()
