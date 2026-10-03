"""Direct mini-sweep over parameter specifications.

No SMM objective and no geography inversion. Each case is one GE solve at a
chosen theta, followed by the scalar population-accounting diagnostic.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dt_cp_model.parameters import build_calibration_setup
from dt_cp_model.solver import accounting_population_scale, run_model_cp_dt
from dt_cp_model.theta import apply_theta
from dt_cp_model.utils import make_grid


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a direct mini-sweep of GE solves at hand-picked theta specs.")
    parser.add_argument("--setup", choices=["fast", "benchmark"], default="fast")
    parser.add_argument("--max-iter-eq", type=int, default=120)
    parser.add_argument("--n", type=int, default=10)
    parser.add_argument("--json", type=Path, default=Path("benchmarks/mini_parameter_solve_2026_05_05.json"))
    args = parser.parse_args()

    setup = build_calibration_setup(args.setup)
    cases = build_cases(setup.x0, setup.names)[: args.n]
    if not cases:
        raise SystemExit("no cases requested")

    out: dict[str, object] = {
        "setup": args.setup,
        "max_iter_eq": args.max_iter_eq,
        "theta_names": list(setup.names),
        "cases": [],
    }

    outside_value = None
    t_all = time.perf_counter()
    for idx, case in enumerate(cases, start=1):
        label = str(case["label"])
        theta = np.asarray(case["theta"], dtype=float)
        P = build_calibration_setup(args.setup).P_base
        P = apply_theta(P, theta, setup.names)
        P.max_iter_eq = args.max_iter_eq

        print(f"[{idx:02d}/{len(cases):02d}] {label}: solving", flush=True)
        t0 = time.perf_counter()
        try:
            sol, P_out, p_eq = run_model_cp_dt(P, verbose=False)
            elapsed = time.perf_counter() - t0
            scale = accounting_population_scale(
                sol,
                P_out,
                make_grid(P_out),
                outside_value=outside_value,
                calibrate_outside_value=outside_value is None,
            )
            if outside_value is None:
                outside_value = float(scale.outside_value)
                out["baseline_outside_value"] = outside_value
            rec = summarize_case(label, theta, sol, P_out, p_eq, scale, elapsed, error=None)
            print(
                f"    ok {elapsed:6.1f}s  err={rec['best_eq_error']:.4g} "
                f"reason={rec['convergence_reason']} TFR={rec['tfr']:.3f} "
                f"own={rec['own_rate']:.3f} scale={rec['scale_factor']:.3f}",
                flush=True,
            )
        except Exception as exc:  # pragma: no cover - diagnostic script
            elapsed = time.perf_counter() - t0
            rec = {
                "label": label,
                "theta": theta.tolist(),
                "elapsed_sec": elapsed,
                "error": repr(exc),
            }
            print(f"    FAILED {elapsed:6.1f}s  {exc!r}", flush=True)
        out["cases"].append(rec)
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(out, indent=2, sort_keys=True))

    out["elapsed_sec"] = time.perf_counter() - t_all
    args.json.write_text(json.dumps(out, indent=2, sort_keys=True))
    print(f"wrote {args.json} in {out['elapsed_sec']:.1f}s", flush=True)


def build_cases(x0: np.ndarray, names: list[str]) -> list[dict[str, object]]:
    base = np.asarray(x0, dtype=float).copy()
    idx = {name: i for i, name in enumerate(names)}

    def case(label: str, **updates: float) -> dict[str, object]:
        theta = base.copy()
        for name, value in updates.items():
            theta[idx[name]] = float(value)
        return {"label": label, "theta": theta}

    return [
        case("baseline_x0"),
        case("lower_fertility_value", psi_child=0.040),
        case("higher_fertility_value", psi_child=0.085),
        case("patient_low_child_value", beta=0.970, psi_child=0.045),
        case("impatient_high_child_value", beta=0.940, psi_child=0.085),
        case("bigger_first_child_housing", h_bar_jump=1.70, h_bar_n=0.70),
        case("flatter_child_housing", h_bar_jump=0.80, h_bar_n=0.35),
        case("lower_owner_premium", chi=1.045),
        case("mobile_high_loc_noise", mu_move=2.0, kappa_loc=2.5),
        case("sticky_low_loc_noise", mu_move=8.0, kappa_loc=0.9),
        case("stronger_bequests", theta0=0.90, theta_n=0.45),
        case("weaker_bequests", theta0=0.25, theta_n=0.10),
        case("lower_entry_wealth", b_entry_fixed=0.05),
        case("higher_entry_wealth", b_entry_fixed=0.35),
        case("low_fertility_dispersion", kappa_fert=1.60),
        case("high_fertility_dispersion", kappa_fert=5.50),
        case("low_baseline_housing_floor", h_bar_0=3.0),
        case("high_baseline_housing_floor", h_bar_0=4.4),
        case("low_child_costs", c_bar_n=0.07, h_bar_n=0.45),
        case("high_child_costs", c_bar_n=0.20, h_bar_n=0.90),
    ]


def summarize_case(
    label: str,
    theta: np.ndarray,
    sol: SimpleNamespace,
    P: SimpleNamespace,
    p_eq: np.ndarray,
    scale: SimpleNamespace,
    elapsed: float,
    error: str | None,
) -> dict[str, object]:
    timings = getattr(sol, "timings", {})
    return {
        "label": label,
        "theta": theta.tolist(),
        "elapsed_sec": elapsed,
        "error": error,
        "p_eq": [float(x) for x in p_eq],
        "best_eq_error": float(timings.get("best_eq_error", np.nan)),
        "final_eq_error": float(timings.get("final_eq_error", np.nan)),
        "iterations_completed": int(timings.get("iterations_completed", 0)),
        "accepted": bool(timings.get("accepted", False)),
        "strict_converged": bool(timings.get("strict_converged", False)),
        "convergence_reason": timings.get("convergence_reason", ""),
        "tfr": float(2 * sol.mean_parity),
        "childless_rate": float(getattr(sol, "parity_dist", np.array([np.nan]))[0]),
        "own_rate": float(getattr(sol, "own_rate", np.nan)),
        "own_rate_3055": float(getattr(sol, "own_rate_3055", np.nan)),
        "pop_share": [float(x) for x in sol.pop_share],
        "mean_age_first_birth": float(getattr(sol, "mean_age_first_birth", np.nan)),
        "migration_rate_2245": float(getattr(sol, "migration_rate_2245", np.nan)),
        "housing_increment_0to1": float(getattr(sol, "housing_increment_0to1_eventstudy_t3", np.nan)),
        "housing_increment_1to2": float(getattr(sol, "housing_increment_1to2_proxy_t3", np.nan)),
        "prime_childless_renter_median_rooms": float(getattr(sol, "prime_childless_renter_median_rooms", np.nan)),
        "prime_childless_owner_median_rooms": float(getattr(sol, "prime_childless_owner_median_rooms", np.nan)),
        "scale_factor": float(scale.scale_factor),
        "scale_denominator": float(scale.denominator),
        "scale_residual": float(scale.stationary_scale_residual),
        "entry_residual": float(scale.stationary_entry_residual),
        "entry_relative_residual": float(scale.stationary_entry_relative_residual),
        "outside_value": float(scale.outside_value),
        "outside_entry_prob": float(scale.outside_entry_prob),
        "city_entry_prob_total": float(scale.city_entry_prob_total),
        "total_mass": float(getattr(sol, "total_mass", np.nan)),
        "E_loc": [float(x) for x in np.asarray(P.E_loc).reshape(-1)],
        "r_bar": [float(x) for x in np.asarray(P.r_bar).reshape(-1)],
    }


if __name__ == "__main__":
    main()
