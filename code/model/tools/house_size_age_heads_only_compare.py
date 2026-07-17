"""Heads-only model-vs-ACS plots:
  A. Ownership rate and mean rooms by age, model vs heads-only ACS.
  B. Inverted view: rooms (size) on x-axis, age distribution stacked,
     side-by-side panels for ACS and model.
"""

from __future__ import annotations
import csv
from collections import defaultdict
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

REPO = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
MODEL_CSV = REPO / "output/model/house_size_by_age_v1/house_size_by_age_bins.csv"
ACS_AGE = REPO / "code/data/mms_center_periphery/output_heads_only_age_rooms/acs_heads_age_profile.csv"
ACS_CELLS = REPO / "code/data/mms_center_periphery/output_heads_only_age_rooms/acs_heads_rooms_x_age_by_tenure.csv"
OUTDIR = REPO / "output/model/house_size_by_age_v1"


def load_model_bins():
    rows = []
    with open(MODEL_CSV) as f:
        for row in csv.DictReader(f):
            lo, hi = row["age_bin"].split("-")
            lo, hi = int(lo), int(hi)
            r = float(row["renter_mass"])
            owner_dict = {}
            for k, v in row.items():
                if k.startswith("owner_H") and "mass_H=" in k:
                    H = float(k.split("H=")[1])
                    owner_dict[H] = float(v)
            rows.append((lo, hi, r, owner_dict))
    return rows


def load_acs_age():
    out = {}
    with open(ACS_AGE) as f:
        for row in csv.DictReader(f):
            age = int(row["age"])
            out[age] = dict(
                owner_rate=float(row["owner_rate"]),
                mean_rooms=float(row["mean_rooms"]),
                pop_weight=float(row["pop_weight"]),
            )
    return out


def load_acs_cells():
    cells = []
    with open(ACS_CELLS) as f:
        for row in csv.DictReader(f):
            try:
                cells.append(dict(
                    rooms_bin=row["rooms_label"],
                    age_bin=row["age_label"],
                    tenure=row["tenure_label"],
                    weight=float(row["cell_weight"]),
                ))
            except (ValueError, KeyError):
                continue
    return cells


def figA_heads_only_lifecycle():
    rows = load_model_bins()
    acs = load_acs_age()

    centers, mod_own, mod_mr = [], [], []
    acs_own, acs_mr = [], []
    RENTER_PROXY = 4.8
    for lo, hi, r, owner_dict in rows:
        c = (lo + hi) / 2
        owner_mass = sum(owner_dict.values())
        total = r + owner_mass
        centers.append(c)
        mod_own.append(owner_mass / total if total > 0 else np.nan)
        if total > 0:
            mod_mr.append((r * RENTER_PROXY + sum(m * H for H, m in owner_dict.items())) / total)
        else:
            mod_mr.append(np.nan)
        # ACS: average across integer ages in bin, weighted by pop_weight
        ws, sums_own, sums_mr = 0.0, 0.0, 0.0
        for age in range(lo, hi + 1):
            if age in acs:
                w = acs[age]["pop_weight"]
                ws += w
                sums_own += acs[age]["owner_rate"] * w
                sums_mr += acs[age]["mean_rooms"] * w
        acs_own.append(sums_own / ws if ws > 0 else np.nan)
        acs_mr.append(sums_mr / ws if ws > 0 else np.nan)

    centers = np.array(centers)
    plt.rcParams.update({"font.size": 13, "axes.titlesize": 13, "axes.labelsize": 12})
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0), constrained_layout=True)
    axes[0].plot(centers, np.array(mod_own) * 100, marker="o", color="#C73E3A", linewidth=2.5, label="Model", markersize=8)
    axes[0].plot(centers, np.array(acs_own) * 100, marker="s", color="#1f5fa6", linewidth=2.5, label="ACS (heads only)", markersize=8, linestyle="--")
    axes[0].set_xlabel("Age (bin center)"); axes[0].set_ylabel("Ownership rate (%)")
    axes[0].set_title("Ownership by Age: Model vs ACS (heads-only)")
    axes[0].grid(True, alpha=0.3); axes[0].legend(); axes[0].set_ylim(0, 100)

    axes[1].plot(centers, mod_mr, marker="o", color="#C73E3A", linewidth=2.5, label="Model", markersize=8)
    axes[1].plot(centers, acs_mr, marker="s", color="#1f5fa6", linewidth=2.5, label="ACS (heads only)", markersize=8, linestyle="--")
    axes[1].set_xlabel("Age (bin center)"); axes[1].set_ylabel("Mean rooms (head's unit)")
    axes[1].set_title("Mean Rooms by Age: Model vs ACS (heads-only)")
    axes[1].grid(True, alpha=0.3); axes[1].legend()
    out = OUTDIR / "model_vs_acs_lifecycle_heads_only.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print(f"wrote {out}")

    # Print headline numbers
    print("\nOwnership rate at key ages:")
    for c, mo, ao in zip(centers, mod_own, acs_own):
        print(f"  age {c:.0f}: model {mo*100:5.1f}%   ACS {ao*100:5.1f}%")
    print("\nMean rooms at key ages:")
    for c, mr_m, mr_a in zip(centers, mod_mr, acs_mr):
        print(f"  age {c:.0f}: model {mr_m:.2f}   ACS {mr_a:.2f}")


def figB_inverted_size_to_age():
    """For each size (or rooms bin), show stacked age distribution.
    Side-by-side: ACS (data) and Model.
    """
    AGES = ["22-39", "40-59", "60+"]
    # Map fine ACS bins to coarse
    ACS_FINE_TO_COARSE = {
        "22-29": "22-39", "30-39": "22-39",
        "40-49": "40-59", "50-59": "40-59",
        "60-69": "60+",   "70+": "60+",
    }
    AGE_COLORS = ["#5276A6", "#E1812C", "#C73E3A"]  # blue young, orange mid, red old

    # --- ACS side ---
    acs_cells = load_acs_cells()
    # All-tenure shares: sum owner + renter
    # Group by rooms_bin -> coarse age_bin -> total weight
    acs_table = defaultdict(lambda: defaultdict(float))
    for c in acs_cells:
        coarse = ACS_FINE_TO_COARSE.get(c["age_bin"], c["age_bin"])
        acs_table[c["rooms_bin"]][coarse] += c["weight"]
    ACS_BINS = ["<=4", "5", "6", "7-8", "9-10", "11+"]
    acs_shares = np.zeros((len(ACS_BINS), len(AGES)))
    for i, rb in enumerate(ACS_BINS):
        tot = sum(acs_table[rb].values())
        for j, ab in enumerate(AGES):
            acs_shares[i, j] = (acs_table[rb][ab] / tot * 100) if tot > 0 else 0

    # --- Model side ---
    rows = load_model_bins()
    # Map model age bins -> ACS age bins (5-yr -> 10-yr)
    # 20-24 + 25-29 -> 22-29 (approx)
    # 30-34 + 35-39 -> 30-39
    # 40-44 + 45-49 -> 40-49
    # 50-54 + 55-59 -> 50-59
    # 60-64 + 65-69 -> 60-69
    # 70-74 + 75-79 -> 70+
    bin_map = {
        "20-24": "22-39", "25-29": "22-39",
        "30-34": "22-39", "35-39": "22-39",
        "40-44": "40-59", "45-49": "40-59",
        "50-54": "40-59", "55-59": "40-59",
        "60-64": "60+",   "65-69": "60+",
        "70-74": "60+",   "75-79": "60+",
    }
    # Model size labels: renter, H=6.0, H=8.0, H=9.5 (active slots)
    H_total = defaultdict(float)
    for lo, hi, r, owner_dict in rows:
        for H, m in owner_dict.items():
            H_total[H] += m
    active_H = sorted([H for H, t in H_total.items() if t > 1e-5])
    MODEL_LABELS = ["Renter"] + [f"H={H:g}" for H in active_H]
    n_model_cols = len(MODEL_LABELS)

    # Per-model-bin, per ACS-age-bin, sum mass
    model_mass = defaultdict(lambda: defaultdict(float))  # model_col -> age_bin -> mass
    for lo, hi, r, owner_dict in rows:
        ab = bin_map[f"{lo}-{hi}"]
        model_mass["Renter"][ab] += r
        for H in active_H:
            model_mass[f"H={H:g}"][ab] += owner_dict[H]

    model_shares = np.zeros((n_model_cols, len(AGES)))
    for i, col in enumerate(MODEL_LABELS):
        tot = sum(model_mass[col].values())
        for j, ab in enumerate(AGES):
            model_shares[i, j] = (model_mass[col][ab] / tot * 100) if tot > 0 else 0

    # --- Plot ---
    plt.rcParams.update({"font.size": 12, "axes.titlesize": 13, "axes.labelsize": 11})
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5), constrained_layout=True)

    for ax, labels, shares, title in [
        (axes[0], ACS_BINS, acs_shares, "ACS (heads only): age within each room bin"),
        (axes[1], MODEL_LABELS, model_shares, "Model: age within each housing type"),
    ]:
        x = np.arange(len(labels))
        bottoms = np.zeros(len(labels))
        for j, ab in enumerate(AGES):
            ax.bar(x, shares[:, j], bottom=bottoms, color=AGE_COLORS[j], label=ab,
                   edgecolor="white", linewidth=0.4, width=0.72)
            # Annotate share if > 7%
            for i, v in enumerate(shares[:, j]):
                if v >= 7:
                    ax.text(x[i], bottoms[i] + v/2, f"{v:.0f}%", ha="center", va="center",
                            fontsize=9, color="white", fontweight="bold")
            bottoms += shares[:, j]
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_xlabel("Housing type / rooms bin")
        ax.set_ylabel("Share within column (%)")
        ax.set_title(title)
        ax.set_ylim(0, 100)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    axes[1].legend(title="Age", loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=10, frameon=False)
    out = OUTDIR / "rooms_on_x_age_stacked_model_vs_acs.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print(f"wrote {out}")

    # Print: 60+ share by housing type
    print("\nACS: 60+ share within each rooms bin:")
    for i, rb in enumerate(ACS_BINS):
        print(f"  {rb}: {acs_shares[i, AGES.index('60+')]:.1f}%")
    print("\nModel: 60+ share within each housing type:")
    for i, col in enumerate(MODEL_LABELS):
        print(f"  {col}: {model_shares[i, AGES.index('60+')]:.1f}%")


def main():
    figA_heads_only_lifecycle()
    figB_inverted_size_to_age()


if __name__ == "__main__":
    main()
