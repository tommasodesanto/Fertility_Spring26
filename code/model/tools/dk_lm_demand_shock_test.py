"""DK / LM demand-shock variant test.

Baseline + two demand-side counterfactuals:
  (A) chi bumped 5%   -- uniform preference shift toward ownership; raises
      both p_C and p_P; closest single-knob analog to a Saiz-style demand shock.
  (B) E_C bumped +0.10 -- positive amenity shift in the center; raises p_C only,
      generating a within-location asymmetry similar to a single-MSA shock.

The DK / LM prediction is that the demand shock raises owner fertility and
lowers renter fertility. Earlier supply-contraction test failed in strict form
because supply contraction also reduces total housing services; this test
isolates the demand-side capital-gain channel.

Note: the May 6 calibration is not the live benchmark closure. Read results
as suggestive comparative statics, not final benchmark counterfactuals.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
sys.path.insert(0, str(REPO / "code/model"))

from dt_cp_model.direct_calibration import build_direct_calibration_setup, evaluate_direct_theta

BEST_JSON = REPO / "code/cluster/results_python_direct_geometry_py_direct_renewal_calibrated_global_12h_20260506/direct_geometry_best.json"


def load_best_theta():
    data = json.load(open(BEST_JSON))
    best = min(data, key=lambda x: float(x.get("loss", float("inf"))))
    return np.asarray(best["theta"], dtype=float), best


def tenure_stratified_parent_rates(moments: dict):
    own_rate = float(moments.get("own_rate", float("nan")))
    childless = float(moments.get("childless_rate", float("nan")))
    gap = float(moments.get("own_family_gap", float("nan")))
    if not all(np.isfinite([own_rate, childless, gap])):
        return None, None, None, None
    p_parent = 1.0 - childless
    P_own_given_nonparent = own_rate - gap * p_parent
    P_own_given_parent = P_own_given_nonparent + gap
    P_rent = 1.0 - own_rate
    if P_rent <= 0 or own_rate <= 0:
        return None, None, None, None
    P_parent_given_own = P_own_given_parent * p_parent / own_rate
    P_parent_given_rent = (1.0 - P_own_given_parent) * p_parent / P_rent
    return P_parent_given_own, P_parent_given_rent, P_own_given_parent, P_own_given_nonparent


def summarize(label: str, result):
    print(f"\n========== {label} ==========")
    print(f"  Elapsed: {result.elapsed_sec:.1f}s   ok={result.solve_ok}   converged={result.converged}")
    print(f"  p_eq (P, C): {result.p_eq}")
    P_par_own, P_par_rent, P_own_par, P_own_nonpar = tenure_stratified_parent_rates(result.moments)
    if P_par_own is not None:
        print(f"  P(parent | owner)  = {P_par_own:.4f}")
        print(f"  P(parent | renter) = {P_par_rent:.4f}")
        print(f"  Owner-renter wedge = {P_par_own - P_par_rent:+.4f}")
    keys = ["tfr", "childless_rate", "own_rate", "own_family_gap",
            "housing_increment_0to1", "housing_increment_1to2",
            "inv_pop_share_C", "inv_rent_ratio_C_over_P"]
    for k in keys:
        v = result.moments.get(k, None)
        if v is not None and np.isfinite(float(v)):
            print(f"  {k}: {float(v):.4f}")
    return {
        "label": label,
        "p_eq": [float(x) for x in result.p_eq],
        "P_parent_given_owner": P_par_own,
        "P_parent_given_renter": P_par_rent,
        "P_own_given_parent": P_own_par,
        "P_own_given_nonparent": P_own_nonpar,
        "moments": {k: float(result.moments[k]) for k in keys if k in result.moments and np.isfinite(float(result.moments[k]))},
        "elapsed_sec": result.elapsed_sec,
        "converged": result.converged,
    }


def main():
    t0 = time.perf_counter()
    theta, best_entry = load_best_theta()
    print(f"Loaded May 6 best theta: job_{best_entry.get('job_id')} eval {best_entry.get('eval_id')} loss {best_entry.get('loss'):.4f}")

    # parameter indices (renewal-valve direct-geometry vector, 15 entries)
    # 0 beta, 1 b_entry_fixed, 2 psi_child, 3 h_bar_jump, 4 h_bar_n, 5 c_bar_n,
    # 6 kappa_fert, 7 chi, 8 kappa_loc, 9 mu_move, 10 theta0, 11 theta_n,
    # 12 h_bar_0, 13 E_C, 14 r_bar_C
    CHI_IX = 7
    E_C_IX = 13

    # BASELINE
    setup_base = build_direct_calibration_setup(
        setup_mode="benchmark",
        population_closure="renewal_valve_calibrated",
        geo_weight=100.0,
    )
    print(f"\nBaseline chi = {theta[CHI_IX]:.4f}; E_C = {theta[E_C_IX]:.4f}")
    base = summarize("BASELINE", evaluate_direct_theta(theta, setup_base, verbose=False))

    # COUNTERFACTUAL A: chi bumped 5%
    setup_a = build_direct_calibration_setup(
        setup_mode="benchmark",
        population_closure="renewal_valve_calibrated",
        geo_weight=100.0,
    )
    # Expand calibration-search bounds for the counterfactual (we are not searching;
    # we are doing comparative statics, so the bounds are not binding economically).
    setup_a.ub = setup_a.ub.copy()
    setup_a.ub[CHI_IX] = 2.0
    theta_a = theta.copy()
    theta_a[CHI_IX] = theta[CHI_IX] * 1.05  # chi up 5%
    print(f"\nCounterfactual A: chi bumped 5% to {theta_a[CHI_IX]:.4f}")
    cfa = summarize("COUNTERFACTUAL A: chi x 1.05", evaluate_direct_theta(theta_a, setup_a, verbose=False))

    # COUNTERFACTUAL B: E_C bumped +0.10
    setup_b = build_direct_calibration_setup(
        setup_mode="benchmark",
        population_closure="renewal_valve_calibrated",
        geo_weight=100.0,
    )
    theta_b = theta.copy()
    theta_b[E_C_IX] = theta[E_C_IX] + 0.10  # E_C +0.10 (center amenity up)
    print(f"\nCounterfactual B: E_C bumped +0.10 to {theta_b[E_C_IX]:.4f}")
    cfb = summarize("COUNTERFACTUAL B: E_C + 0.10", evaluate_direct_theta(theta_b, setup_b, verbose=False))

    # ---- VERDICT ----
    print("\n============================================================")
    print("DK/LM DEMAND-SHOCK ASYMMETRY CHECK")
    print("============================================================")
    for label, cf in [("A: chi x 1.05", cfa), ("B: E_C + 0.10", cfb)]:
        if cf["P_parent_given_owner"] is None or base["P_parent_given_owner"] is None: continue
        d_own = cf["P_parent_given_owner"] - base["P_parent_given_owner"]
        d_rent = cf["P_parent_given_renter"] - base["P_parent_given_renter"]
        dp_P = cf["p_eq"][0] / base["p_eq"][0] - 1.0
        dp_C = cf["p_eq"][1] / base["p_eq"][1] - 1.0
        print(f"\n  {label}")
        print(f"    %Delta p_P = {dp_P*100:+.2f}%   %Delta p_C = {dp_C*100:+.2f}%")
        print(f"    Delta P(parent|owner)  : {d_own*100:+.2f}pp  (DK predict POSITIVE)")
        print(f"    Delta P(parent|renter) : {d_rent*100:+.2f}pp  (DK predict NEGATIVE)")
        print(f"    Owner-renter wedge     : {(d_own - d_rent)*100:+.2f}pp")
        verdict = "PASS" if (d_own > 0 and d_rent < 0) else ("PARTIAL" if d_own * d_rent < 0 else "FAIL")
        print(f"    Verdict: {verdict}")

    out = {"baseline": base, "counterfactual_chi_x1.05": cfa, "counterfactual_E_C_plus_0.10": cfb,
           "theta": theta.tolist(), "best_entry_meta": {k: best_entry.get(k) for k in ["job_id","eval_id","loss","run_tag"]}}
    out_path = REPO / "output/model/dk_lm_asymmetry_test"
    json.dump(out, open(out_path / "dk_lm_demand_shock_test_v1.json", "w"), indent=2)
    print(f"\nTotal elapsed: {time.perf_counter()-t0:.1f}s")


if __name__ == "__main__":
    main()
