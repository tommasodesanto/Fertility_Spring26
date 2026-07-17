"""Two-in-one model exercise:
 (1) Hacamo quantitative analog: vary phi at fixed theta; trace ownership,
     joint P(own AND parent), TFR, prices.
 (2) Ablation: chi-bump under phi=0.80 vs phi=0.99 (near-free borrowing). If
     ownership response to chi shock collapses under phi=0.99, the constraint
     was doing the work.

Baseline theta = May 6 best (loss 8.13); preliminary calibration caveat.
"""

from __future__ import annotations
import json, sys, time, copy
from pathlib import Path
from dataclasses import asdict
import numpy as np

REPO = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
sys.path.insert(0, str(REPO / "code/model"))

from dt_cp_model.direct_calibration import build_direct_calibration_setup, evaluate_direct_theta

BEST_JSON = REPO / "code/cluster/results_python_direct_geometry_py_direct_renewal_calibrated_global_12h_20260506/direct_geometry_best.json"


def load_theta():
    data = json.load(open(BEST_JSON))
    best = min(data, key=lambda x: float(x["loss"]))
    return np.asarray(best["theta"], dtype=float)


def tenure_stratified(moments):
    own = moments.get("own_rate", float("nan"))
    cless = moments.get("childless_rate", float("nan"))
    gap = moments.get("own_family_gap", float("nan"))
    if not all(np.isfinite([own, cless, gap])):
        return None
    p_par = 1 - cless
    P_own_npar = own - gap * p_par
    P_own_par = P_own_npar + gap
    P_par_own = P_own_par * p_par / own if own > 0 else float("nan")
    P_par_rent = (1 - P_own_par) * p_par / (1 - own) if own < 1 else float("nan")
    joint_own_par = P_own_par * p_par  # P(own AND parent) = P(own|parent)*P(parent)
    return dict(
        P_own_par=P_own_par, P_own_npar=P_own_npar,
        P_par_own=P_par_own, P_par_rent=P_par_rent,
        joint_own_par=joint_own_par,
    )


def solve_one(theta, phi_value, chi_mult=1.0, label=""):
    CHI_IX = 7
    setup = build_direct_calibration_setup(
        setup_mode="benchmark",
        population_closure="renewal_valve_calibrated",
        geo_weight=100.0,
    )
    setup.P_base.phi = float(phi_value) * np.ones_like(setup.P_base.phi)
    setup.ub = setup.ub.copy()
    setup.ub[CHI_IX] = 2.0
    theta_use = theta.copy()
    theta_use[CHI_IX] = theta[CHI_IX] * chi_mult
    t0 = time.perf_counter()
    res = evaluate_direct_theta(theta_use, setup, verbose=False)
    print(f"[{label}] elapsed {res.elapsed_sec:.1f}s converged={res.converged} loss={res.loss:.3f}")
    return res


def fmt(d):
    keys = ["tfr", "childless_rate", "own_rate", "own_family_gap",
            "housing_increment_0to1", "housing_increment_1to2",
            "inv_pop_share_C", "inv_rent_ratio_C_over_P"]
    out = {k: float(d.get(k)) for k in keys if k in d and np.isfinite(float(d[k]))}
    return out


def report(label, result):
    t = tenure_stratified(result.moments)
    print(f"\n== {label} ==")
    print(f"  prices (P, C): {result.p_eq}")
    print(f"  TFR={result.moments.get('tfr', float('nan')):.4f}  own={result.moments.get('own_rate', float('nan')):.4f}  childless={result.moments.get('childless_rate', float('nan')):.4f}")
    if t:
        print(f"  P(parent|owner)={t['P_par_own']:.4f}  P(parent|renter)={t['P_par_rent']:.4f}")
        print(f"  joint P(own AND parent)={t['joint_own_par']:.4f}")
    return dict(label=label, p_eq=[float(x) for x in result.p_eq],
                moments=fmt(result.moments), tenure=t,
                converged=result.converged, elapsed=result.elapsed_sec)


def main():
    theta = load_theta()
    t0 = time.perf_counter()

    out = {}

    # (1) Hacamo quantitative analog: vary phi at baseline chi
    out["baseline_phi080"] = report("baseline phi=0.80", solve_one(theta, 0.80, label="phi080"))
    out["phi085"] = report("phi=0.85 (Hacamo-style relaxation)", solve_one(theta, 0.85, label="phi085"))
    out["phi095"] = report("phi=0.95 (FHA-like)", solve_one(theta, 0.95, label="phi095"))

    # (2) Ablation: chi shock under tight vs near-no constraint
    out["chi105_phi080"] = report("chi*1.05 under phi=0.80", solve_one(theta, 0.80, chi_mult=1.05, label="chi*1.05@phi080"))
    out["chi105_phi099"] = report("chi*1.05 under phi=0.99 (ablation)", solve_one(theta, 0.99, chi_mult=1.05, label="chi*1.05@phi099"))
    out["baseline_phi099"] = report("baseline phi=0.99 (ablation reference)", solve_one(theta, 0.99, label="phi099"))

    print("\n" + "="*70)
    print("HACAMO QUANTITATIVE ANALOG: shift in joint(own AND parent) as phi rises")
    print("="*70)
    j080 = out["baseline_phi080"]["tenure"]["joint_own_par"]
    j085 = out["phi085"]["tenure"]["joint_own_par"]
    j095 = out["phi095"]["tenure"]["joint_own_par"]
    print(f"  phi=0.80 -> 0.85: delta joint = {(j085-j080)*100:+.2f}pp")
    print(f"  phi=0.80 -> 0.95: delta joint = {(j095-j080)*100:+.2f}pp")
    print(f"  Hacamo reduced-form headline: +6 pp")

    print("\n" + "="*70)
    print("ABLATION: chi*1.05 response under phi=0.80 vs phi=0.99")
    print("="*70)
    own080_base = out["baseline_phi080"]["moments"]["own_rate"]
    own080_chi  = out["chi105_phi080"]["moments"]["own_rate"]
    own099_base = out["baseline_phi099"]["moments"]["own_rate"]
    own099_chi  = out["chi105_phi099"]["moments"]["own_rate"]
    print(f"  phi=0.80: own_rate baseline {own080_base:.4f} -> chi*1.05 {own080_chi:.4f}  (delta {(own080_chi-own080_base)*100:+.2f}pp)")
    print(f"  phi=0.99: own_rate baseline {own099_base:.4f} -> chi*1.05 {own099_chi:.4f}  (delta {(own099_chi-own099_base)*100:+.2f}pp)")
    j080_chi = out["chi105_phi080"]["tenure"]["joint_own_par"]
    j099_chi = out["chi105_phi099"]["tenure"]["joint_own_par"]
    j099_base = out["baseline_phi099"]["tenure"]["joint_own_par"]
    print(f"  phi=0.80: joint(own,parent) baseline {j080:.4f} -> chi*1.05 {j080_chi:.4f}  (delta {(j080_chi-j080)*100:+.2f}pp)")
    print(f"  phi=0.99: joint(own,parent) baseline {j099_base:.4f} -> chi*1.05 {j099_chi:.4f}  (delta {(j099_chi-j099_base)*100:+.2f}pp)")

    out_path = REPO / "output/model/dk_lm_asymmetry_test/phi_shift_and_ablation_v1.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(out_path, "w"), indent=2)
    print(f"\nWrote {out_path}\nTotal elapsed: {time.perf_counter()-t0:.1f}s")


if __name__ == "__main__":
    main()
