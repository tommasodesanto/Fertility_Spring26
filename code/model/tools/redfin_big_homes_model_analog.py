"""Build the model analog of the Redfin big-homes-by-generation chart.

Redfin chart: for owner-occupied 3+ bedroom homes, share by (generation,
household type). Model analog: for owner units with H >= 6, share by
(age-bin, current-children).

Mapping:
  Generation       <- Age bin (proxy, since model has no calendar year)
    22-39  ~  Millennial
    40-59  ~  Gen X
    60+    ~  Boomer/Silent
  Household type   <- current_children flag
    has kids at home  ~  "with own children"
    no kids at home   ~  "1-2 adults" (model has no spouses)

Output: two bars: (data, model) side by side for headline cells, or
horizontal bar chart matched to Redfin style.
"""

from __future__ import annotations
import argparse
import json, sys, time, csv, copy
from pathlib import Path
from types import SimpleNamespace
import numpy as np

REPO = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
sys.path.insert(0, str(REPO / "code/model"))

from dt_cp_model.direct_calibration import build_direct_calibration_setup, DIRECT_GEOMETRY_NAMES
from dt_cp_model.direct_calibration import OUTSIDE_VALUE_NAME, RENEWAL_FLOW_NAME
from dt_cp_model.solver import run_model_cp_dt
from dt_cp_model.parameters import asdict
from dt_cp_model.theta import apply_theta

BEST_JSON = REPO / "code/cluster/results_python_direct_geometry_py_direct_renewal_calibrated_global_12h_20260506/direct_geometry_best.json"
OUTDIR = REPO / "output/model/house_size_by_age_v1"


def load_theta():
    data = json.load(open(BEST_JSON))
    best = min(data, key=lambda x: float(x["loss"]))
    return np.asarray(best["theta"], dtype=float)


def solve_baseline():
    setup = build_direct_calibration_setup(
        setup_mode="benchmark",
        population_closure="renewal_valve_calibrated",
        geo_weight=100.0,
    )
    theta = load_theta()
    P_base = SimpleNamespace(**copy.deepcopy(asdict(setup.P_base)))
    structural_names = [n for n in setup.names if n not in DIRECT_GEOMETRY_NAMES and n not in ("outside_value", "outside_entry_flow")]
    theta_dict = {n: float(v) for n, v in zip(setup.names, theta)}
    P = apply_theta(P_base, [theta_dict[n] for n in structural_names], structural_names)
    P.E_loc = np.array([float(P.E_loc[0]), theta_dict["E_C"]])
    P.r_bar = np.array([float(P.r_bar[0]), theta_dict["r_bar_C"]])
    print("Solving baseline at May 6 best theta...", flush=True)
    t0 = time.perf_counter()
    sol, P_out, p_eq = run_model_cp_dt(P, verbose=False)
    print(f"Solved in {time.perf_counter()-t0:.1f}s", flush=True)
    return sol, P_out


def solve_record(record_path: Path):
    record = json.loads(record_path.read_text())
    alpha_cons_bounds = record.get("alpha_cons_bounds")
    if alpha_cons_bounds is not None:
        alpha_cons_bounds = tuple(float(x) for x in alpha_cons_bounds)
    setup = build_direct_calibration_setup(
        setup_mode="benchmark",
        population_closure="outside_option_benchmark_normalized",
        geo_weight=100.0,
        scale_target=1.0,
        scale_weight=100.0,
        hR_max=record.get("hR_max", 8.0),
        alpha_cons=record.get("alpha_cons"),
        owner_h_bar_scale=record.get("owner_h_bar_scale"),
        owner_size_cost=record.get("owner_size_cost"),
        owner_size_cost_ref=record.get("owner_size_cost_ref"),
        owner_size_cost_power=record.get("owner_size_cost_power"),
        tenure_choice_kappa=record.get("tenure_choice_kappa"),
        weight_overrides=record.get("weight_overrides") or None,
        extra_targets=record.get("extra_targets") or None,
        calibrate_alpha_cons=bool(record.get("calibrate_alpha_cons", False)),
        alpha_cons_bounds=alpha_cons_bounds,
        H_own=record.get("H_own"),
    )
    theta = np.asarray(record["theta"], dtype=float)
    P_base = SimpleNamespace(**copy.deepcopy(asdict(setup.P_base)))
    theta_dict = {n: float(v) for n, v in zip(setup.names, theta)}
    structural_names = [
        n for n in setup.names
        if n not in DIRECT_GEOMETRY_NAMES and n not in (OUTSIDE_VALUE_NAME, RENEWAL_FLOW_NAME)
    ]
    P = apply_theta(P_base, [theta_dict[n] for n in structural_names], structural_names)
    P.E_loc = np.array([float(P.E_loc[0]), theta_dict["E_C"]])
    P.r_bar = np.array([float(P.r_bar[0]), theta_dict["r_bar_C"]])
    if OUTSIDE_VALUE_NAME in theta_dict:
        P.outside_value = theta_dict[OUTSIDE_VALUE_NAME]
        P.outside_value_is_calibrated = True
        P.allow_uncalibrated_outside_value = False
    if RENEWAL_FLOW_NAME in theta_dict:
        P.outside_entry_flow = theta_dict[RENEWAL_FLOW_NAME]
    print(f"Solving record {record_path}...", flush=True)
    t0 = time.perf_counter()
    sol, P_out, p_eq = run_model_cp_dt(P, verbose=False)
    print(f"Solved in {time.perf_counter()-t0:.1f}s", flush=True)
    return sol, P_out


def current_children(nn: int, cs: int, P) -> int:
    if nn <= 0 or cs == 0 or cs > P.n_child_stages:
        return 0
    return int(nn)


def extract_cells(sol, P):
    """Return masses by (age, h_slot, has_kids).
    g.shape = (Nb, h_slot, I, J, n, s)
    """
    g = sol.g
    Nb, n_h_slots, I, J, n_parity, n_cs = g.shape
    age_start = P.age_start
    # Pre-compute has_kids for each (n, s)
    has_kids = np.zeros((n_parity, n_cs), dtype=bool)
    for nn in range(n_parity):
        for cs in range(n_cs):
            has_kids[nn, cs] = current_children(nn, cs, P) > 0
    # Sum g over (Nb, I) -> (n_h_slots, J, n_parity, n_cs)
    g_red = g.sum(axis=0).sum(axis=1)  # axes 0 (Nb) and 2 (I); after first sum axis index of I becomes 1
    # Actually careful: after first sum over axis 0, axes become (n_h_slots, I, J, n, s).
    # Sum axis 1 next gives (n_h_slots, J, n, s).
    # h_slot axis: 0 (renter), 1..n_house (owner with H_own[i-1])
    # For each (h_slot, age) compute mass with kids vs without
    mass_with = np.zeros((n_h_slots, J))
    mass_without = np.zeros((n_h_slots, J))
    for nn in range(n_parity):
        for cs in range(n_cs):
            if has_kids[nn, cs]:
                mass_with += g_red[:, :, nn, cs]
            else:
                mass_without += g_red[:, :, nn, cs]
    # Compute age (real year)
    ages = np.arange(J) + age_start
    return ages, mass_with, mass_without


def main():
    parser = argparse.ArgumentParser(description="Build model analog of Redfin big-home ownership chart.")
    parser.add_argument("--record", type=Path, help="Calibration best.json record to solve instead of the legacy default.")
    parser.add_argument("--outdir", type=Path, default=OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.record is None:
        sol, P = solve_baseline()
    else:
        sol, P = solve_record(args.record)
    ages, mass_with, mass_without = extract_cells(sol, P)
    H_own = np.asarray(P.H_own)
    print(f"H_own = {H_own.tolist()}")
    # Owner slots index 1..n_house correspond to H_own[0..n_house-1]
    # Identify "big" owner slots: H >= 6
    big_slots = [k for k, H in enumerate(H_own, start=1) if H >= 6]
    print(f"Big slots (h_slot index, H value): {[(k, float(H_own[k-1])) for k in big_slots]}")

    # For each big slot and age (in real years), sum mass_with and mass_without
    # Aggregate to age generation bins
    def gen_bin(a):
        if a < 40: return "Young (22-39)"
        if a < 60: return "Middle (40-59)"
        return "Old (60+)"

    cells = {}  # (gen, hh_type) -> total mass
    for k in big_slots:
        for j, a in enumerate(ages):
            if a < 22: continue
            g = gen_bin(int(a))
            cells[(g, "with own children")] = cells.get((g, "with own children"), 0) + mass_with[k, j]
            cells[(g, "1-2 adults")] = cells.get((g, "1-2 adults"), 0) + mass_without[k, j]

    total = sum(cells.values())
    print(f"\nTotal big-owner mass: {total:.4f}")
    print("\nModel analog: share of big-owner units (H>=6) by (age-bin, household-type):")
    items = sorted(cells.items(), key=lambda kv: -kv[1])
    for (g, t), m in items:
        print(f"  {g:>20s}, {t:>20s}: {m/total*100:5.1f}%")

    # Save CSV
    csv_path = args.outdir / "redfin_big_homes_model_analog.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["generation_bin", "hh_type", "share"])
        for (g, t), m in items:
            w.writerow([g, t, m/total])
    print(f"\nWrote {csv_path}")

    # Also: compute the analog of "Boomer 1-2 adults" share
    boom_empty = cells.get(("Old (60+)", "1-2 adults"), 0) / total
    young_w_kids = cells.get(("Young (22-39)", "with own children"), 0) / total
    print(f"\nModel 'empty nester old (60+)' share of big-owner units: {boom_empty*100:.1f}%")
    print(f"Model 'young families' share of big-owner units: {young_w_kids*100:.1f}%")
    print(f"\nFor comparison: Redfin headline 'Boomer 1-2 adults' = 27.8%; my ACS replica = 22.4%")


if __name__ == "__main__":
    main()
