"""Sweep phi across many levels, both UNIVERSAL relaxation (phi[*] = phi) and
PARENT-TARGETED relaxation (phi[0] = 0.80 fixed for non-parents, phi[1:] = phi
for parities >=1). Checkpoint after each solve so partial results survive
hangs.

Outcomes per solve: ownership, TFR, childless, joint(own AND parent), prices.
"""

from __future__ import annotations
import csv, json, sys, time
from pathlib import Path
import numpy as np

REPO = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
sys.path.insert(0, str(REPO / "code/model"))

from dt_cp_model.direct_calibration import build_direct_calibration_setup, evaluate_direct_theta

BEST_JSON = REPO / "code/cluster/results_python_direct_geometry_py_direct_renewal_calibrated_global_12h_20260506/direct_geometry_best.json"
OUTDIR = REPO / "output/model/phi_sweep_v1"
OUTDIR.mkdir(parents=True, exist_ok=True)
CSV = OUTDIR / "phi_sweep_results.csv"


def load_theta():
    data = json.load(open(BEST_JSON))
    best = min(data, key=lambda x: float(x["loss"]))
    return np.asarray(best["theta"], dtype=float)


def tenure_stratified(m):
    own = m.get("own_rate", float("nan"))
    cless = m.get("childless_rate", float("nan"))
    gap = m.get("own_family_gap", float("nan"))
    if not all(np.isfinite([own, cless, gap])):
        return float("nan"), float("nan"), float("nan"), float("nan")
    p_par = 1 - cless
    P_own_npar = own - gap * p_par
    P_own_par = P_own_npar + gap
    joint = P_own_par * p_par
    P_par_own = P_own_par * p_par / own if own > 0 else float("nan")
    P_par_rent = (1 - P_own_par) * p_par / (1 - own) if own < 1 else float("nan")
    return joint, P_par_own, P_par_rent, P_own_par


def solve_one(theta, phi_value, scenario):
    """scenario in {'universal', 'parent_targeted'}"""
    setup = build_direct_calibration_setup(
        setup_mode="benchmark",
        population_closure="renewal_valve_calibrated",
        geo_weight=100.0,
    )
    n_parity = len(setup.P_base.phi)
    if scenario == "universal":
        setup.P_base.phi = float(phi_value) * np.ones(n_parity)
    elif scenario == "parent_targeted":
        new_phi = np.array(setup.P_base.phi, dtype=float)
        new_phi[0] = 0.80
        new_phi[1:] = float(phi_value)
        setup.P_base.phi = new_phi
    else:
        raise ValueError(scenario)
    t0 = time.perf_counter()
    res = evaluate_direct_theta(theta, setup, verbose=False)
    elapsed = time.perf_counter() - t0
    return res, elapsed


def append_row(row):
    write_header = not CSV.exists()
    with open(CSV, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(row.keys()))
        if write_header: w.writeheader()
        w.writerow(row)


def main():
    if CSV.exists():
        CSV.unlink()  # fresh start
    theta = load_theta()
    phi_grid = [0.80, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.93]
    t0 = time.perf_counter()

    for phi in phi_grid:
        for scenario in ["universal", "parent_targeted"]:
            print(f"\n>>> phi={phi:.3f}  scenario={scenario}", flush=True)
            try:
                res, elapsed = solve_one(theta, phi, scenario)
            except Exception as e:
                print(f"   FAILED: {e}", flush=True)
                append_row(dict(phi=phi, scenario=scenario, error=repr(e), elapsed=float("nan")))
                continue
            joint, p_po, p_pr, p_op = tenure_stratified(res.moments)
            m = res.moments
            row = dict(
                phi=phi, scenario=scenario,
                converged=int(bool(res.converged)),
                elapsed=round(elapsed, 1),
                tfr=float(m.get("tfr", float("nan"))),
                own_rate=float(m.get("own_rate", float("nan"))),
                childless_rate=float(m.get("childless_rate", float("nan"))),
                own_family_gap=float(m.get("own_family_gap", float("nan"))),
                joint_own_par=joint,
                P_par_own=p_po,
                P_par_rent=p_pr,
                p_P=float(res.p_eq[0]),
                p_C=float(res.p_eq[1]),
                inv_rent_ratio=float(m.get("inv_rent_ratio_C_over_P", float("nan"))),
                housing_increment_0to1=float(m.get("housing_increment_0to1", float("nan"))),
            )
            append_row(row)
            print(f"   elapsed {elapsed:.1f}s  TFR={row['tfr']:.4f}  own={row['own_rate']:.4f}  joint={joint:.4f}", flush=True)

    print(f"\nTotal elapsed: {time.perf_counter()-t0:.1f}s", flush=True)
    print(f"CSV saved to {CSV}", flush=True)


if __name__ == "__main__":
    main()
