"""Plot phi sweep results: joint(own AND parent), TFR, ownership, by scenario."""

from __future__ import annotations
import csv
from collections import defaultdict
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

REPO = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
CSV = REPO / "output/model/phi_sweep_v1/phi_sweep_results.csv"
OUTDIR = REPO / "output/model/phi_sweep_v1"


def load():
    data = defaultdict(lambda: defaultdict(list))  # scenario -> col -> values
    with open(CSV) as f:
        for row in csv.DictReader(f):
            sc = row.get("scenario")
            if not sc: continue
            for k, v in row.items():
                if k in ("scenario",): continue
                try:
                    data[sc][k].append(float(v) if v not in ("", "nan", "NaN") else float("nan"))
                except (ValueError, TypeError):
                    data[sc][k].append(float("nan"))
    return data


def main():
    data = load()
    if not data:
        print("no data yet"); return

    plt.rcParams.update({"font.size": 13, "axes.titlesize": 13, "axes.labelsize": 12})

    SC_LABELS = {"universal": "Universal $\\phi$", "parent_targeted": "Parent-targeted $\\phi$"}
    SC_COLORS = {"universal": "#5276A6", "parent_targeted": "#E1812C"}

    # Headline: joint(own AND parent) vs phi
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0), constrained_layout=True)
    for sc, d in data.items():
        x = np.array(d["phi"])
        idx = np.argsort(x)
        x = x[idx]
        y = np.array(d["joint_own_par"])[idx]
        axes[0].plot(x, y * 100, marker="o", color=SC_COLORS[sc], label=SC_LABELS[sc], linewidth=2, markersize=8)
        tfr = np.array(d["tfr"])[idx]
        axes[1].plot(x, tfr, marker="o", color=SC_COLORS[sc], label=SC_LABELS[sc], linewidth=2, markersize=8)

    axes[0].set_xlabel("$\\phi$ (LTV cap; $1-\\phi$ = DP requirement)")
    axes[0].set_ylabel("Joint $P(\\mathrm{own} \\cap \\mathrm{parent})$, %")
    axes[0].set_title("Joint Margin vs. $\\phi$")
    axes[0].grid(True, alpha=0.3); axes[0].legend()

    axes[1].set_xlabel("$\\phi$ (LTV cap)")
    axes[1].set_ylabel("TFR")
    axes[1].set_title("TFR vs. $\\phi$")
    axes[1].grid(True, alpha=0.3); axes[1].legend()

    out = OUTDIR / "phi_sweep_joint_and_tfr.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    print(f"wrote {out}")

    # Second figure: ownership and conditional fertility
    fig2, axes2 = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
    for sc, d in data.items():
        x = np.array(d["phi"])
        idx = np.argsort(x)
        x = x[idx]
        axes2[0].plot(x, np.array(d["own_rate"])[idx] * 100, marker="o", color=SC_COLORS[sc], label=SC_LABELS[sc], linewidth=2, markersize=7)
        axes2[1].plot(x, np.array(d["P_par_own"])[idx] * 100, marker="o", color=SC_COLORS[sc], label=SC_LABELS[sc], linewidth=2, markersize=7)
        axes2[2].plot(x, np.array(d["P_par_rent"])[idx] * 100, marker="o", color=SC_COLORS[sc], label=SC_LABELS[sc], linewidth=2, markersize=7)
    for ax, ttl, ylab in [
        (axes2[0], "Ownership Rate", "% owners"),
        (axes2[1], "$P(\\mathrm{parent} \\mid \\mathrm{owner})$", "%"),
        (axes2[2], "$P(\\mathrm{parent} \\mid \\mathrm{renter})$", "%"),
    ]:
        ax.set_xlabel("$\\phi$")
        ax.set_ylabel(ylab)
        ax.set_title(ttl)
        ax.grid(True, alpha=0.3); ax.legend()
    out2 = OUTDIR / "phi_sweep_ownership_and_conditional.png"
    fig2.savefig(out2, dpi=180, bbox_inches="tight")
    print(f"wrote {out2}")


if __name__ == "__main__":
    main()
