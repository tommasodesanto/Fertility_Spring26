"""DK / LM owner-renter sign reversal test.

Comparative-statics exercise on the May 6 best calibration (renewal_valve_calibrated closure).

Solves two stationary equilibria:
  (1) BASELINE -- May 6 best theta, default H0.
  (2) COUNTERFACTUAL -- same theta, H0 reduced by 10% in both locations.

The supply contraction raises equilibrium owner-house prices and rents.
DK / LM prediction: in the counterfactual, owners' parent rate should rise
(housing-wealth channel) and renters' parent rate should fall (rent-cost
channel).

The model exposes own_family_gap and own_rate, from which we back out
parent rate conditional on tenure. We also dump prices, ownership, housing
increments, and other relevant moments.
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
    """Back out P(parent|own) and P(parent|rent) from own_rate, childless_rate, own_family_gap."""
    own_rate = float(moments.get("own_rate", float("nan")))
    childless = float(moments.get("childless_rate", float("nan")))
    gap = float(moments.get("own_family_gap", float("nan")))
    if not all(np.isfinite([own_rate, childless, gap])):
        return None, None, None, None
    p_parent = 1.0 - childless
    p_nonparent = childless
    # gap = P(own|parent) - P(own|nonparent)
    # own_rate = P(own|parent)*p_parent + P(own|nonparent)*p_nonparent
    # Let x = P(own|nonparent). Then P(own|parent) = x + gap.
    # own_rate = (x + gap) * p_parent + x * p_nonparent = x + gap * p_parent
    # So x = own_rate - gap * p_parent
    P_own_given_nonparent = own_rate - gap * p_parent
    P_own_given_parent = P_own_given_nonparent + gap
    # By Bayes:
    P_rent = 1.0 - own_rate
    if P_rent <= 0 or own_rate <= 0:
        return None, None, None, None
    P_parent_given_own = P_own_given_parent * p_parent / own_rate
    P_parent_given_rent = (1.0 - P_own_given_parent) * p_parent / P_rent
    # Sanity: alternative via second tenure cell
    P_parent_given_rent_alt = (1.0 - P_own_given_nonparent) * p_nonparent  # this would be wrong; just sanity
    return P_parent_given_own, P_parent_given_rent, P_own_given_parent, P_own_given_nonparent


def report(label: str, result, P_base_H0):
    print(f"\n========== {label} ==========")
    print(f"  Solve elapsed: {result.elapsed_sec:.1f}s   solve_ok={result.solve_ok}   converged={result.converged}")
    print(f"  Loss: {result.loss:.4f}")
    print(f"  Prices p_eq (P, C): {result.p_eq}")
    P_par_own, P_par_rent, P_own_par, P_own_nonpar = tenure_stratified_parent_rates(result.moments)
    print(f"  Tenure-stratified parent rates:")
    print(f"    P(parent | owner)   = {P_par_own:.4f}" if P_par_own is not None else "    P(parent|owner) N/A")
    print(f"    P(parent | renter)  = {P_par_rent:.4f}" if P_par_rent is not None else "    P(parent|renter) N/A")
    print(f"    P(own | parent)     = {P_own_par:.4f}" if P_own_par is not None else "")
    print(f"    P(own | non-parent) = {P_own_nonpar:.4f}" if P_own_nonpar is not None else "")
    keys = [
        "tfr", "childless_rate", "mean_age_first_birth", "tfr_gradient",
        "own_rate", "own_gradient", "own_family_gap",
        "prime_childless_renter_median_rooms", "prime_childless_owner_median_rooms",
        "housing_increment_0to1", "housing_increment_1to2",
        "young_liquid_wealth_to_income",
        "center_share_nonparents", "center_share_newparents",
        "migration_rate", "old_age_own_rate", "old_age_parent_childless_gap",
        "inv_pop_share_C", "inv_rent_ratio_C_over_P",
    ]
    print(f"  Headline moments:")
    for k in keys:
        v = result.moments.get(k, None)
        if v is None or not np.isfinite(float(v)): continue
        print(f"    {k}: {float(v):.4f}")
    print(f"  H0 used: {P_base_H0.tolist()}")
    return {
        "label": label,
        "p_eq": [float(x) for x in result.p_eq],
        "P_parent_given_owner": P_par_own,
        "P_parent_given_renter": P_par_rent,
        "P_own_given_parent": P_own_par,
        "P_own_given_nonparent": P_own_nonpar,
        "moments": {k: float(result.moments[k]) for k in keys if k in result.moments and np.isfinite(float(result.moments[k]))},
        "H0": [float(x) for x in P_base_H0],
        "solve_elapsed_sec": result.elapsed_sec,
        "converged": result.converged,
    }


def main():
    t0 = time.perf_counter()
    theta, best_entry = load_best_theta()
    print(f"Loaded May 6 best theta from job_{best_entry.get('job_id')}, eval {best_entry.get('eval_id')}, loss {best_entry.get('loss'):.4f}")
    print(f"theta = {theta.tolist()}")

    # BASELINE
    setup_base = build_direct_calibration_setup(
        setup_mode="benchmark",
        population_closure="renewal_valve_calibrated",
        geo_weight=100.0,
    )
    H0_base = setup_base.P_base.H0.copy()
    print(f"\nBaseline H0: {H0_base.tolist()}")
    result_base = evaluate_direct_theta(theta, setup_base, verbose=False)
    base_summary = report("BASELINE", result_base, H0_base)

    # COUNTERFACTUAL: H0 * 0.9
    setup_cf = build_direct_calibration_setup(
        setup_mode="benchmark",
        population_closure="renewal_valve_calibrated",
        geo_weight=100.0,
    )
    setup_cf.P_base.H0 = setup_cf.P_base.H0 * 0.9
    H0_cf = setup_cf.P_base.H0.copy()
    print(f"\nCounterfactual H0 (x0.9): {H0_cf.tolist()}")
    result_cf = evaluate_direct_theta(theta, setup_cf, verbose=False)
    cf_summary = report("COUNTERFACTUAL  H0 x 0.9", result_cf, H0_cf)

    # Compare
    print("\n============================================================")
    print("DK/LM ASYMMETRY CHECK")
    print("============================================================")
    if base_summary["P_parent_given_owner"] is not None and cf_summary["P_parent_given_owner"] is not None:
        d_own = cf_summary["P_parent_given_owner"] - base_summary["P_parent_given_owner"]
        d_rent = cf_summary["P_parent_given_renter"] - base_summary["P_parent_given_renter"]
        print(f"  Delta P(parent|owner)   : {d_own:+.4f}   (DK/LM predict POSITIVE)")
        print(f"  Delta P(parent|renter)  : {d_rent:+.4f}   (DK/LM predict NEGATIVE)")
        print(f"  Owner-renter wedge      : {d_own - d_rent:+.4f}")
        verdict = "PASS" if (d_own > 0 and d_rent < 0) else ("PARTIAL" if d_own * d_rent < 0 else "FAIL")
        print(f"  Verdict: {verdict}")

    out = {"baseline": base_summary, "counterfactual_H0_x0.9": cf_summary,
           "theta": theta.tolist(), "best_entry_meta": {k: best_entry.get(k) for k in ["job_id","eval_id","loss","run_tag"]}}
    out_path = REPO / "output/model/dk_lm_asymmetry_test"
    out_path.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(out_path / "dk_lm_asymmetry_test_v1.json", "w"), indent=2)
    print(f"\nWrote {out_path}/dk_lm_asymmetry_test_v1.json")
    print(f"\nTotal elapsed: {time.perf_counter()-t0:.1f}s")


if __name__ == "__main__":
    main()
