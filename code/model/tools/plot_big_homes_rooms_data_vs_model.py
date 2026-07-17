"""Side-by-side bars: rooms-based 'who owns the big homes' for data and model.

Apples-to-apples with the model: both restrict to owner units with at least
6 rooms (the model's H_own grid is in rooms, so H >= 6 lines up with the data
filter rooms >= 6).
"""

from __future__ import annotations
import csv
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

REPO = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
DATA_CSV = REPO / "code/data/mms_center_periphery/output_big_homes_rooms/big_homes_rooms_cells.csv"
MODEL_CSV = REPO / "output/model/house_size_by_age_v1/redfin_big_homes_model_analog.csv"
OUT = REPO / "output/model/house_size_by_age_v1/big_homes_rooms_data_vs_model.png"


def main():
    gen_to_bin = {"Boomer": "Old (60+)", "Silent": "Old (60+)",
                  "Gen X": "Middle (40-59)",
                  "Millennial": "Young (22-39)", "Gen Z": "Young (22-39)"}
    typ_to_bin = {"1-2 adults": "no children <18",
                  "3+ people, no children": "no children <18",
                  "with own children": "children <18"}

    data_cells = {}
    with open(DATA_CSV) as f:
        for r in csv.DictReader(f):
            g = gen_to_bin.get(r["generation"], r["generation"])
            t = typ_to_bin.get(r["hh_type"], r["hh_type"])
            data_cells[(g, t)] = data_cells.get((g, t), 0) + float(r["share"])

    model_cells = {}
    with open(MODEL_CSV) as f:
        for r in csv.DictReader(f):
            g = r["generation_bin"]
            t = "no children <18" if r["hh_type"] == "1-2 adults" else "children <18"
            model_cells[(g, t)] = model_cells.get((g, t), 0) + float(r["share"])

    GENS = ["Young (22-39)", "Middle (40-59)", "Old (60+)"]
    TYPES = ["no children <18", "children <18"]
    bar_labels = [f"{g}\n{t}" for g in GENS for t in TYPES]

    data_vals = [data_cells.get((g, t), 0) * 100 for g in GENS for t in TYPES]
    model_vals = [model_cells.get((g, t), 0) * 100 for g in GENS for t in TYPES]

    print("Cell      | Data | Model")
    for lbl, dv, mv in zip(bar_labels, data_vals, model_vals):
        print(f"  {lbl.replace(chr(10),'/'):<40s} {dv:5.1f}%   {mv:5.1f}%")

    plt.rcParams.update({"font.size": 12, "axes.titlesize": 13, "axes.labelsize": 12})
    fig, ax = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
    x = np.arange(len(bar_labels))
    width = 0.38
    bars_d = ax.bar(x - width/2, data_vals, width, color="#1f5fa6", label="ACS 2022-2023, 6+ rooms, owner-occ")
    bars_m = ax.bar(x + width/2, model_vals, width, color="#C73E3A", label=r"Model, $H \geq 6$ owner units")
    for b, v in zip(bars_d, data_vals):
        ax.text(b.get_x() + b.get_width()/2, v + 0.5, f"{v:.0f}%", ha="center", fontsize=10, color="#1f5fa6")
    for b, v in zip(bars_m, model_vals):
        ax.text(b.get_x() + b.get_width()/2, v + 0.5, f"{v:.0f}%", ha="center", fontsize=10, color="#C73E3A")
    ax.set_xticks(x); ax.set_xticklabels(bar_labels, fontsize=10)
    ax.set_ylabel("Share of large homes (%)")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper right", frameon=False)
    ax.set_ylim(0, max(max(data_vals), max(model_vals)) * 1.18)
    fig.savefig(OUT, dpi=120, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
