"""Model-vs-ACS comparison of lifecycle housing and the inverted 'who owns
the big houses' view from the model.

Two figures:
  A. Ownership rate by age: model vs ACS (clean apples-to-apples)
  B. Within-size age distribution from the model: share of each owner H slot
     by age. Says directly which ages occupy which size of unit.
"""

from __future__ import annotations
import csv
from collections import defaultdict
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

REPO = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
MODEL_CSV = REPO / "output/model/house_size_by_age_v1/house_size_by_age_bins.csv"
ACS_CSV = REPO / "code/data/mms_center_periphery/output_middle_center/mms_age_profiles_full.csv"
OUTDIR = REPO / "output/model/house_size_by_age_v1"


def load_model_bins():
    """Returns list of (age_bin_start, age_bin_end, masses dict with H value: mass)"""
    rows = []
    with open(MODEL_CSV) as f:
        for row in csv.DictReader(f):
            bin_label = row["age_bin"]
            lo, hi = bin_label.split("-")
            lo, hi = int(lo), int(hi)
            r = float(row["renter_mass"])
            owner_keys = sorted(k for k in row if k.startswith("owner_H"))
            owner_dict = {}
            for k in owner_keys:
                # parse H value out of column name 'owner_H3_mass_H=6.0'
                H = float(k.split("H=")[1])
                owner_dict[H] = float(row[k])
            rows.append((lo, hi, r, owner_dict))
    return rows


def load_acs_by_age():
    """Returns dict age -> (owner_rate, mean_rooms, weight) aggregated over location."""
    by_age = defaultdict(lambda: dict(num_own=0.0, num_rooms=0.0, weight=0.0))
    with open(ACS_CSV) as f:
        for row in csv.DictReader(f):
            try:
                age = int(row["age"])
                w = float(row["pop_weight"])
                own = float(row["owner_rate"])
                mr = float(row["mean_rooms"])
            except (ValueError, KeyError):
                continue
            if not all(np.isfinite([w, own, mr])): continue
            d = by_age[age]
            d["num_own"] += own * w
            d["num_rooms"] += mr * w
            d["weight"] += w
    out = {}
    for age, d in by_age.items():
        if d["weight"] > 0:
            out[age] = dict(owner_rate=d["num_own"]/d["weight"], mean_rooms=d["num_rooms"]/d["weight"])
    return out


def model_ownership_by_bin(rows):
    """Returns lists of (bin_center, model_own_rate)"""
    centers, owns = [], []
    for lo, hi, r, owner_dict in rows:
        owner_mass = sum(owner_dict.values())
        total = r + owner_mass
        centers.append((lo + hi) / 2)
        owns.append(owner_mass / total if total > 0 else 0)
    return np.array(centers), np.array(owns)


def acs_ownership_by_bin(acs, rows):
    """Aggregate ACS owner_rate to the same bins, weighted by ACS weights."""
    centers, owns = [], []
    for lo, hi, _r, _owner_dict in rows:
        total_w = 0.0
        sum_own_w = 0.0
        for age in range(lo, hi + 1):
            if age in acs:
                # weight proportional to pop_weight; use 1 as a fallback for equality
                # we don't have per-age weights here, so use unit weights
                # (each age contributes equally within the bin)
                sum_own_w += acs[age]["owner_rate"]
                total_w += 1
        centers.append((lo + hi) / 2)
        owns.append(sum_own_w / total_w if total_w > 0 else np.nan)
    return np.array(centers), np.array(owns)


def acs_meanrooms_by_bin(acs, rows):
    centers, mr = [], []
    for lo, hi, _r, _owner_dict in rows:
        total_w = 0
        sum_mr = 0.0
        for age in range(lo, hi + 1):
            if age in acs:
                sum_mr += acs[age]["mean_rooms"]
                total_w += 1
        centers.append((lo + hi) / 2)
        mr.append(sum_mr / total_w if total_w > 0 else np.nan)
    return np.array(centers), np.array(mr)


def model_meanrooms_by_bin(rows):
    """All-household mean rooms by bin (renter ~ 4.8 from baseline diagnostic; owner sizes from H)."""
    # Use the renter mean from baseline diagnostic. The earlier housing_size_fit.png shows
    # the model's prime-age renter median rooms = 4.8 (single bin spike). Use 4.8 as proxy.
    RENTER_PROXY = 4.8
    centers, mr = [], []
    for lo, hi, r, owner_dict in rows:
        owner_mass = sum(owner_dict.values())
        total = r + owner_mass
        if total <= 0:
            centers.append((lo+hi)/2); mr.append(np.nan); continue
        renter_contrib = r * RENTER_PROXY
        owner_contrib = sum(m * H for H, m in owner_dict.items())
        centers.append((lo+hi)/2)
        mr.append((renter_contrib + owner_contrib) / total)
    return np.array(centers), np.array(mr)


def figA_model_vs_acs():
    rows = load_model_bins()
    acs = load_acs_by_age()
    centers, mod_own = model_ownership_by_bin(rows)
    _, acs_own = acs_ownership_by_bin(acs, rows)
    _, mod_mr = model_meanrooms_by_bin(rows)
    _, acs_mr = acs_meanrooms_by_bin(acs, rows)

    plt.rcParams.update({"font.size": 13, "axes.titlesize": 13, "axes.labelsize": 12})
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0), constrained_layout=True)
    axes[0].plot(centers, mod_own * 100, marker="o", color="#C73E3A", linewidth=2.5, label="Model", markersize=8)
    axes[0].plot(centers, acs_own * 100, marker="s", color="#1f5fa6", linewidth=2.5, label="ACS 2012-2023", markersize=8, linestyle="--")
    axes[0].set_xlabel("Age (bin center)")
    axes[0].set_ylabel("Ownership rate (%)")
    axes[0].set_title("Ownership by Age: Model vs ACS")
    axes[0].grid(True, alpha=0.3); axes[0].legend()
    axes[0].set_ylim(0, 105)

    axes[1].plot(centers, mod_mr, marker="o", color="#C73E3A", linewidth=2.5, label="Model", markersize=8)
    axes[1].plot(centers, acs_mr, marker="s", color="#1f5fa6", linewidth=2.5, label="ACS 2012-2023", markersize=8, linestyle="--")
    axes[1].set_xlabel("Age (bin center)")
    axes[1].set_ylabel("Mean rooms (all-household)")
    axes[1].set_title("Mean Rooms by Age: Model vs ACS")
    axes[1].grid(True, alpha=0.3); axes[1].legend()

    out = OUTDIR / "model_vs_acs_lifecycle.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print(f"wrote {out}")


def figB_inverted_age_within_size():
    """For each owner H slot, show the age distribution of its occupants.
    This is 'who owns each house size' - the inverted cut."""
    rows = load_model_bins()
    bin_labels = [f"{lo}-{hi}" for lo, hi, _, _ in rows]
    H_keys = sorted(rows[0][3].keys())
    # Drop empty slots
    H_total = {H: sum(row[3][H] for row in rows) for H in H_keys}
    active_H = [H for H in H_keys if H_total[H] > 1e-5]
    print(f"Active H slots: {active_H}")
    # Also include renter
    renter_total = sum(r for _, _, r, _ in rows)

    plt.rcParams.update({"font.size": 13, "axes.titlesize": 13, "axes.labelsize": 12})
    fig, ax = plt.subplots(figsize=(11.5, 5.5), constrained_layout=True)
    x = np.arange(len(rows))
    width = 0.8 / (len(active_H) + 1)
    offsets = np.linspace(-0.4 + width/2, 0.4 - width/2, len(active_H) + 1)

    # Bars: renter
    renter_shares = np.array([r/renter_total for _, _, r, _ in rows]) * 100
    ax.bar(x + offsets[0], renter_shares, width, color="#5276A6", label=f"Renter (total mass {renter_total:.3f})")
    # Bars: each active H
    colors = plt.cm.YlOrRd(np.linspace(0.3, 0.95, len(active_H)))
    for i, H in enumerate(active_H):
        shares = np.array([row[3][H] / H_total[H] for row in rows]) * 100
        ax.bar(x + offsets[i+1], shares, width, color=colors[i],
               label=f"Owner H={H} (total mass {H_total[H]:.3f})")

    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels, rotation=0, fontsize=10)
    ax.set_xlabel("Age bin")
    ax.set_ylabel("Share of this housing type occupied by age bin (%)")
    ax.set_title("Who Occupies Each Housing Type? (within-size age distribution)")
    ax.legend(loc="upper right", fontsize=10, frameon=True)
    ax.grid(True, axis="y", alpha=0.3)

    out = OUTDIR / "age_within_house_size.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print(f"wrote {out}")

    # Print summary: for the largest active H, which age bin holds most?
    print("\nWho owns the biggest unit (H=8.0)?")
    H = 8.0
    if H in active_H:
        for (lo, hi, _, owner_dict), share_print in zip(rows, [owner_dict[H] / H_total[H] for _, _, _, owner_dict in rows]):
            print(f"  age {lo}-{hi}: {share_print*100:.1f}%")


def main():
    figA_model_vs_acs()
    figB_inverted_age_within_size()


if __name__ == "__main__":
    main()
