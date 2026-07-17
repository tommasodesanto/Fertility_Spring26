"""Hacamo quantitative analog only: vary phi at fixed theta and report joint(own, parent).

Trimmed from phi_shift_and_ablation.py to skip the ablation and the extreme
phi=0.99 case (which previously hung). Solves at phi in {0.80, 0.85, 0.90}.
Output flushed eagerly so we can see progress.
"""

from __future__ import annotations
import json, sys, time
from pathlib import Path
import numpy as np

REPO = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
sys.path.insert(0, str(REPO / "code/model"))

from dt_cp_model.direct_calibration import build_direct_calibration_setup, evaluate_direct_theta

BEST_JSON = REPO / "code/cluster/results_python_direct_geometry_py_direct_renewal_calibrated_global_12h_20260506/direct_geometry_best.json"


def load_theta():
    data = json.load(open(BEST_JSON))
    best = min(data, key=lambda x: float(x["loss"]))
    return np.asarray(best["theta"], dtype=float)


def tenure_stratified(m):
    own = m.get("own_rate", float("nan"))
    cless = m.get("childless_rate", float("nan"))
    gap = m.get("own_family_gap", float("nan"))
    if not all(np.isfinite([own, cless, gap])):
        return None
    p_par = 1 - cless
    P_own_npar = own - gap * p_par
    P_own_par = P_own_npar + gap
    return dict(
        P_own_par=P_own_par,
        P_par_own=P_own_par * p_par / own if own > 0 else float("nan"),
        P_par_rent=(1 - P_own_par) * p_par / (1 - own) if own < 1 else float("nan"),
        joint_own_par=P_own_par * p_par,
    )


def solve_one(theta, phi_value, label=""):
    setup = build_direct_calibration_setup(
        setup_mode="benchmark",
        population_closure="renewal_valve_calibrated",
        geo_weight=100.0,
    )
    setup.P_base.phi = float(phi_value) * np.ones_like(setup.P_base.phi)
    t0 = time.perf_counter()
    print(f"[{label}] starting phi={phi_value}...", flush=True)
    res = evaluate_direct_theta(theta, setup, verbose=False)
    print(f"[{label}] done: {time.perf_counter()-t0:.1f}s   converged={res.converged}   loss={res.loss:.3f}", flush=True)
    return res


def report(label, result):
    t = tenure_stratified(result.moments)
    print(f"\n== {label} ==", flush=True)
    print(f"  prices (P, C): {[round(x,4) for x in result.p_eq]}", flush=True)
    print(f"  TFR={result.moments.get('tfr', float('nan')):.4f}  own={result.moments.get('own_rate', float('nan')):.4f}  childless={result.moments.get('childless_rate', float('nan')):.4f}", flush=True)
    if t:
        print(f"  P(parent|owner)={t['P_par_own']:.4f}  P(parent|renter)={t['P_par_rent']:.4f}", flush=True)
        print(f"  joint P(own AND parent)={t['joint_own_par']:.4f}", flush=True)
    return dict(label=label, p_eq=[float(x) for x in result.p_eq],
                moments={k: float(v) for k, v in result.moments.items() if np.isfinite(float(v) if isinstance(v, (int, float)) else 0)},
                tenure=t, converged=result.converged, elapsed=result.elapsed_sec)


def main():
    theta = load_theta()
    t0 = time.perf_counter()
    out = {}
    for phi in [0.80, 0.85, 0.90]:
        label = f"phi={phi}"
        out[f"phi{int(phi*100):03d}"] = report(label, solve_one(theta, phi, label))
    print("\n" + "="*70, flush=True)
    print("HACAMO QUANTITATIVE ANALOG", flush=True)
    print("="*70, flush=True)
    j080 = out["phi080"]["tenure"]["joint_own_par"]
    j085 = out["phi085"]["tenure"]["joint_own_par"]
    j090 = out["phi090"]["tenure"]["joint_own_par"]
    print(f"  joint(own AND parent):", flush=True)
    print(f"    phi=0.80: {j080:.4f}", flush=True)
    print(f"    phi=0.85: {j085:.4f}  (delta vs baseline: {(j085-j080)*100:+.2f}pp)", flush=True)
    print(f"    phi=0.90: {j090:.4f}  (delta vs baseline: {(j090-j080)*100:+.2f}pp)", flush=True)
    print(f"  Hacamo reduced-form headline (credit shock): +6 pp", flush=True)
    out_path = REPO / "output/model/dk_lm_asymmetry_test/phi_hacamo_analog_v1.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(out_path, "w"), indent=2)
    print(f"\nTotal elapsed: {time.perf_counter()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
