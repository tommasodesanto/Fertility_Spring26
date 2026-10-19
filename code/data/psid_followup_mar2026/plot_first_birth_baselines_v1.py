"""Levels of the first-birth housing outcomes by event window, A2h sample.

Reads output/sa_first_birth_baselines_v1/baseline_means.csv (written by
sa_first_birth_baselines_v1.do) and writes one PDF for the mock appendix.
"""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
SRC = HERE / "output" / "sa_first_birth_baselines_v1" / "baseline_means.csv"
OUT = HERE.parents[2] / "latex" / "JMP_DS_mock" / "figures" / "first_birth_baselines_A2h.pdf"

ORDER = ["Wleft", "Wm2", "Wm1", "Wref", "Wp1", "Wp2", "Wp3", "Wp4", "Wp5", "Wp6", "Wright"]
LABELS = ["$\\leq-8$", "$-7/-6$", "$-5/-4$", "$-3/-2$", "$-1/0$", "$+1/+2$", "$+3/+4$",
          "$+5/+6$", "$+7/+8$", "$+9/+10$", "$\\geq+11$"]

means = {}
for r in csv.DictReader(open(SRC)):
    means.setdefault(r["outcome"], {})[r["window"]] = float(r["mean"])


def series(outcome, scale=1.0):
    return [means[outcome][w] * scale for w in ORDER]


x = list(range(len(ORDER)))
ref = ORDER.index("Wref")
panels = [
    ("(a) Rooms", [("rooms", 1.0, "Rooms", "#1f4e79")], "Rooms"),
    ("(b) Owns the dwelling", [("own", 100.0, "Owns", "#1f4e79")], "Percent"),
    ("(c) Moved since last interview", [("moved", 100.0, "Moved", "#1f4e79")], "Percent"),
    ("(d) Moved, by first reason",
     [("moved_space", 100.0, "For more space", "#1f4e79"), ("moved_nbhd", 100.0, "For neighborhood", "#b5651d")],
     "Percent"),
]

fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.4))
for ax, (title, lines, ylab) in zip(axes.flat, panels):
    for outcome, scale, lab, color in lines:
        y = series(outcome, scale)
        ax.plot(x, y, color=color, marker="o", markersize=4, linewidth=1.3, label=lab)
        ax.plot([ref], [y[ref]], marker="o", markersize=7, markerfacecolor="white",
                markeredgecolor=color, linestyle="none")
        ax.annotate(f"{y[ref]:.2f}" if outcome == "rooms" else f"{y[ref]:.1f}",
                    (ref, y[ref]), textcoords="offset points", xytext=(0, 8),
                    ha="center", fontsize=8, color=color)
    ax.axvspan(ref - 0.4, ref + 0.4, color="0.9", zorder=0)
    ax.set_title(title, fontsize=10, loc="left")
    ax.set_ylabel(ylab, fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(LABELS, fontsize=7, rotation=45)
    ax.tick_params(axis="y", labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    if len(lines) > 1:
        ax.legend(frameon=False, fontsize=8)
for ax in axes[1]:
    ax.set_xlabel("Years relative to first birth", fontsize=9)
fig.tight_layout()
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT)
print(OUT)
