"""Lifecycle distribution of house sizes from the model.

For each 5-year age bin, compute the share of population in
{renter, owner H_1, owner H_2, ..., owner H_K}. Plot stacked area
and supplementary summaries (mean rooms by age, share family-sized).
"""

from __future__ import annotations
import json, sys, time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

REPO = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
sys.path.insert(0, str(REPO / "code/model"))

from dt_cp_model.direct_calibration import build_direct_calibration_setup, evaluate_direct_theta
from dt_cp_model.solver import run_model_cp_dt
from dt_cp_model.parameters import asdict
from dt_cp_model.theta import apply_theta
from dt_cp_model.direct_calibration import DIRECT_GEOMETRY_NAMES
import copy
from types import SimpleNamespace

BEST_JSON = REPO / "code/cluster/results_python_direct_geometry_py_direct_renewal_calibrated_global_12h_20260506/direct_geometry_best.json"
OUTDIR = REPO / "output/model/house_size_by_age_v1"
OUTDIR.mkdir(parents=True, exist_ok=True)


def load_theta():
    data = json.load(open(BEST_JSON))
    best = min(data, key=lambda x: float(x["loss"]))
    return np.asarray(best["theta"], dtype=float), best


def solve_baseline():
    setup = build_direct_calibration_setup(
        setup_mode="benchmark",
        population_closure="renewal_valve_calibrated",
        geo_weight=100.0,
    )
    theta, _ = load_theta()
    print(f"Solving baseline at May 6 best theta (phi={setup.P_base.phi[0]})", flush=True)
    t0 = time.perf_counter()
    res = evaluate_direct_theta(theta, setup, verbose=False)
    print(f"Solved in {res.elapsed_sec:.1f}s. Loss={res.loss:.3f}, converged={res.converged}", flush=True)
    # Re-solve to get sol object (evaluate_direct_theta only returns moments).
    # Apply theta then run_model_cp_dt directly.
    P_base = SimpleNamespace(**copy.deepcopy(asdict(setup.P_base)))
    structural_names = [
        n for n in setup.names
        if n not in DIRECT_GEOMETRY_NAMES and n not in ("outside_value", "outside_entry_flow")
    ]
    theta_dict = {n: float(v) for n, v in zip(setup.names, theta)}
    structural_theta = [theta_dict[n] for n in structural_names]
    P = apply_theta(P_base, structural_theta, structural_names)
    P.E_loc = np.array([float(P.E_loc[0]), theta_dict["E_C"]])
    P.r_bar = np.array([float(P.r_bar[0]), theta_dict["r_bar_C"]])
    print("Re-solving for full sol object...", flush=True)
    t0 = time.perf_counter()
    sol, P_out, p_eq = run_model_cp_dt(P, verbose=False)
    print(f"Got sol in {time.perf_counter()-t0:.1f}s", flush=True)
    return sol, P_out, p_eq


def extract_age_size(sol, P):
    """For each integer age a in [0, J-1], compute mass by (tenure, h_index).
    h_index: 0 = renter, 1..n_house = owner with H_own[h_index-1].
    Returns: ages (real years), mass_array shape (J, n_house+1), mean_rooms_by_age (J,)
    """
    g = sol.g  # forward distribution
    # Shape inspection. The forward distribution iterates over states.
    # We'll integrate over (b, i, z, n, s) for each (a, h, ten).
    # First find which axes.
    print(f"g.shape = {g.shape}", flush=True)
    # In this model the conventional shape is (Nb, Nh+1 or similar, I, Nz, Nn, Ns, J).
    # Tenure is encoded by h slot: index 0 = renter, 1..Nh = owner with H_own[i-1].
    # Check from solver compute_statistics for the axis layout.
    # Use a heuristic: last dim is J (age).
    J = P.J
    n_house = P.n_house
    age_start = P.age_start

    n_slots = n_house + 1
    # Find J axis (== P.J) and h-slot axis (== n_house+1)
    j_candidates = [ax for ax, dim in enumerate(g.shape) if dim == J]
    h_candidates = [ax for ax, dim in enumerate(g.shape) if dim == n_slots]
    if not j_candidates or not h_candidates:
        raise RuntimeError(f"Could not find J={J} or h={n_slots} axis in g.shape={g.shape}")
    j_axis = j_candidates[0]
    h_axis = h_candidates[0]
    print(f"j_axis = {j_axis}, h_axis = {h_axis}", flush=True)

    # Sum over all dims except (j_axis, h_axis)
    sum_axes = tuple(ax for ax in range(g.ndim) if ax not in (j_axis, h_axis))
    mass = g.sum(axis=sum_axes)
    # Order so result is (J, n_slots)
    if j_axis > h_axis:
        mass = mass.T
    if mass.shape != (J, n_slots):
        raise RuntimeError(f"Unexpected mass shape {mass.shape}, wanted ({J},{n_slots})")

    # Compute mean rooms per age = sum (size * mass) / sum(mass)
    sizes = np.zeros(n_slots)
    # renter: use renter housing policy average (approximate via average renter h^R)
    # Try to get hR-realized policy from sol
    h_renter_mean = None
    if hasattr(sol, "h_realized_renter_by_age"):
        h_renter_mean = sol.h_realized_renter_by_age
    # Fallback: use h_bar_0 (renter floor) — known underestimate
    sizes[0] = float(getattr(P, "hR_max", 8.0)) / 2  # rough placeholder; we'll override with mean rooms if available
    sizes[1:] = np.asarray(P.H_own)

    # mean rooms per age, with renter row using the model's renter rooms if computable
    # For renters we approximate by P.hR_max/2 if no policy known.
    # We'll separately compute mean owner size per age (mass-weighted within owners).
    age_mass = mass.sum(axis=1)
    with np.errstate(invalid="ignore"):
        mean_rooms_per_age = np.where(age_mass > 0, (mass * sizes[None, :]).sum(axis=1) / age_mass, np.nan)
        owner_mass = mass[:, 1:].sum(axis=1)
        mean_owner_size_per_age = np.where(owner_mass > 0,
                                           (mass[:, 1:] * sizes[1:][None, :]).sum(axis=1) / owner_mass,
                                           np.nan)

    real_ages = np.arange(J) + age_start
    return real_ages, mass, sizes, mean_rooms_per_age, mean_owner_size_per_age


def aggregate_to_bins(real_ages, mass, sizes, bin_starts=None, bin_width=5):
    if bin_starts is None:
        bin_starts = list(range(20, 80, bin_width))
    bin_labels = [f"{b}-{b+bin_width-1}" for b in bin_starts]
    n_bins = len(bin_starts)
    n_slots = mass.shape[1]
    binned = np.zeros((n_bins, n_slots))
    for i, b in enumerate(bin_starts):
        mask = (real_ages >= b) & (real_ages < b + bin_width)
        binned[i, :] = mass[mask, :].sum(axis=0)
    return bin_labels, binned


def plot(bin_labels, binned, sizes, P, real_ages, mean_rooms_per_age, mean_owner_size_per_age):
    n_bins, n_slots = binned.shape
    n_house = n_slots - 1
    # Normalize each bin to share
    bin_totals = binned.sum(axis=1, keepdims=True)
    shares = np.where(bin_totals > 0, binned / bin_totals, 0)

    # Two-panel figure
    plt.rcParams.update({"font.size": 13, "axes.titlesize": 13, "axes.labelsize": 12})
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5), constrained_layout=True)

    # Panel 1: stacked bars for tenure x size by age bin
    ax = axes[0]
    H_own = np.asarray(P.H_own)
    labels = ["Renter"] + [f"Owner H={H_own[k]:.1f}" for k in range(n_house)]
    # Colors: renter in blue gradient, owners in warm gradient
    colors = ["#5276A6"] + list(plt.cm.YlOrRd(np.linspace(0.25, 0.95, n_house)))
    bottoms = np.zeros(n_bins)
    x = np.arange(n_bins)
    for s in range(n_slots):
        ax.bar(x, shares[:, s] * 100, bottom=bottoms * 100,
               color=colors[s], label=labels[s], width=0.78,
               edgecolor="white", linewidth=0.4)
        bottoms += shares[:, s]
    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels, rotation=0, fontsize=10)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Share within age bin (%)")
    ax.set_xlabel("Age")
    ax.set_title("Housing tenure and size by age")
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=9, frameon=False)

    # Panel 2: mean rooms / mean owner-size by age
    ax2 = axes[1]
    real_age_mask = (real_ages >= 20) & (real_ages <= 80)
    ax2.plot(real_ages[real_age_mask], mean_owner_size_per_age[real_age_mask],
             color="#C73E3A", linewidth=2.5, label="Owners only: mean house size", marker="o", markersize=4)
    ax2.set_xlabel("Age")
    ax2.set_ylabel("Mean owner house size (model H units)")
    ax2.set_title("Mean owner house size over lifecycle")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="best")
    # Add horizontal reference lines for H grid
    for h in H_own:
        ax2.axhline(h, color="gray", alpha=0.2, linewidth=0.5, linestyle=":")
    out = OUTDIR / "house_size_by_age.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    print(f"wrote {out}", flush=True)

    # Print numerical summary
    print("\nShare within bin (rows sum to 1):", flush=True)
    print(f"{'bin':>10}  " + "  ".join(f"{lbl:>10}" for lbl in labels), flush=True)
    for i, lbl in enumerate(bin_labels):
        print(f"{lbl:>10}  " + "  ".join(f"{shares[i,s]*100:>9.1f}%" for s in range(n_slots)), flush=True)


def main():
    sol, P, p_eq = solve_baseline()
    real_ages, mass, sizes, mean_rooms_per_age, mean_owner_size_per_age = extract_age_size(sol, P)
    bin_labels, binned = aggregate_to_bins(real_ages, mass, sizes)
    plot(bin_labels, binned, sizes, P, real_ages, mean_rooms_per_age, mean_owner_size_per_age)
    # Save CSV too
    out_csv = OUTDIR / "house_size_by_age_bins.csv"
    with open(out_csv, "w") as f:
        n_slots = mass.shape[1]
        H_own = np.asarray(P.H_own)
        cols = ["age_bin", "renter_mass"] + [f"owner_H{i+1}_mass_H={H_own[i]:.1f}" for i in range(P.n_house)]
        f.write(",".join(cols) + "\n")
        for i, lbl in enumerate(bin_labels):
            row = [lbl] + [f"{binned[i,s]:.6f}" for s in range(n_slots)]
            f.write(",".join(row) + "\n")
    print(f"wrote {out_csv}", flush=True)


if __name__ == "__main__":
    main()
