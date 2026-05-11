"""Dump curated solver slices for the interactive tour.

Run once at the calibrated theta to produce tour_data.json, which the
assemble step inlines into tour.html.

Usage:
    python tools/dump_explorer_data.py
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from dt_cp_model.parameters import build_calibration_setup
from dt_cp_model.solver import run_model_cp_dt
from dt_cp_model.theta import apply_theta


# Use the project's x0 starting point. We could overlay the four memorized
# calibrated values from `calibration_dt_april` (beta=0.947, chi=1.08,
# kappa_loc=1.92, h_bar_jump=1.43), but the other 9 entries paired with those
# during the actual calibration aren't memorized — and substituting just 4
# into the x0 defaults yields nonsense moments. The honest choice is x0.
CALIBRATED_OVERRIDES = {}


TARGET_TO_SOL_FIELD = {
    "tfr": ("mean_parity", lambda sol: 2 * float(sol.mean_parity)),
    "childless_rate": ("parity_dist", lambda sol: float(sol.parity_dist[0])),
    "mean_age_first_birth": ("mean_age_first_birth", lambda sol: float(sol.mean_age_first_birth)),
    "tfr_gradient": ("mean_parity_by_loc", lambda sol: 2 * float(sol.mean_parity_by_loc[0] - sol.mean_parity_by_loc[1])),
    "own_rate": ("own_rate_3055", lambda sol: float(sol.own_rate_3055)),
    "own_gradient": ("own_gradient_3055", lambda sol: float(sol.own_gradient_3055)),
    "own_family_gap": ("own_gap_newparent_nonparent_3055", lambda sol: float(sol.own_gap_newparent_nonparent_3055)),
    "prime_childless_renter_median_rooms": ("prime_childless_renter_median_rooms", lambda sol: float(sol.prime_childless_renter_median_rooms)),
    "prime_childless_owner_median_rooms": ("prime_childless_owner_median_rooms", lambda sol: float(sol.prime_childless_owner_median_rooms)),
    "housing_increment_0to1": ("housing_increment_0to1_eventstudy_t3", lambda sol: float(getattr(sol, "housing_increment_0to1_eventstudy_t3", np.nan))),
    "housing_increment_1to2": ("housing_increment_1to2_proxy_t3", lambda sol: float(getattr(sol, "housing_increment_1to2_proxy_t3", np.nan))),
    "young_liquid_wealth_to_income": ("young_liquid_wealth_to_income", lambda sol: float(sol.young_liquid_wealth_to_income)),
    "center_share_nonparents": ("center_share_nonparents_2245", lambda sol: float(sol.center_share_nonparents_2245)),
    "center_share_newparents": ("center_share_newparents_2245", lambda sol: float(sol.center_share_newparents_2245)),
    "migration_rate": ("migration_rate_2245", lambda sol: float(sol.migration_rate_2245)),
    "old_age_own_rate": ("old_age_own_rate_6575", lambda sol: float(sol.old_age_own_rate_6575)),
    "old_age_parent_childless_gap": ("old_age_parent_childless_gap_6575", lambda sol: float(sol.old_age_parent_childless_gap_6575)),
}


def build_theta_from_x0(setup) -> np.ndarray:
    """Start from setup.x0; override the four memory-recorded calibrated entries."""
    names = list(setup.names)
    theta = np.array(setup.x0, dtype=float).copy()
    for k, v in CALIBRATED_OVERRIDES.items():
        if k in names:
            theta[names.index(k)] = v
    return theta


def slice_policies(sol, P) -> dict:
    """Pick a representative slice of the 6-D arrays for the explorer.

    Array shapes (verified by probe):
      bp_pol, c_pol, tenure_choice, hR_pol : (Nb, nt, I_origin, J, npar, ncs)
      loc_probs                            : (Nb, nt, I_origin, I_dest, J, npar, ncs)
      fert_probs                           : (Nb, nt, I_origin, J, npar)  (no cs axis)

    Slicing axes:
      - j (age): every 5th age, plus j=0 and j=J-1.
      - n (parity): all.
      - cs (child stage): {0, 1, ncs // 2, ncs - 1}.
      - ten_origin: all.
      - i_origin: all.
    """
    bp = sol.bp_pol
    Nb, nt, I, J, npar, ncs = bp.shape
    age_idx = sorted(set(list(range(0, J, 5)) + [J - 1]))
    cs_idx = sorted(set([0, 1, ncs // 2, ncs - 1]))

    out = {}
    for j in age_idx:
        for n in range(npar):
            for cs in cs_idx:
                for to in range(nt):
                    for i in range(I):
                        key = f"{j}|{n}|{cs}|{to}|{i}"
                        slice_data = {
                            "bp": sol.bp_pol[:, to, i, j, n, cs].astype(float).round(5).tolist(),
                            "c":  sol.c_pol[:, to, i, j, n, cs].astype(float).round(5).tolist(),
                            "tn": sol.tenure_choice[:, to, i, j, n, cs].astype(int).tolist(),
                        }
                        # loc_probs[b, ten, i_orig, i_dest, j, n, cs] -> P(i' | i)
                        if hasattr(sol, "loc_probs") and sol.loc_probs is not None:
                            slice_data["loc_probs"] = sol.loc_probs[:, to, i, :, j, n, cs].astype(float).round(4).tolist()
                        # fert_probs[b, ten, i, j, n'] only meaningful at (n=0, cs=0)
                        if hasattr(sol, "fert_probs") and sol.fert_probs is not None and n == 0 and cs == 0:
                            try:
                                slice_data["fert_probs"] = sol.fert_probs[:, to, i, j, :].astype(float).round(4).tolist()
                            except Exception:
                                pass
                        out[key] = slice_data
    return {
        "slices": out,
        "axes": {
            "j": age_idx,
            "cs": cs_idx,
            "npar": int(npar),
            "nt": int(nt),
            "I": int(I),
            "Nb": int(Nb),
            "J": int(J),
            "ncs": int(ncs),
        },
    }


def extract_moment_dict(sol):
    out = {}
    for target_name, (sol_attr, getter) in TARGET_TO_SOL_FIELD.items():
        try:
            v = getter(sol)
            if np.isfinite(v):
                out[target_name] = v
        except Exception:
            pass
    return out


def _to_serializable(obj):
    if obj is None:
        return None
    if isinstance(obj, dict):
        return {k: _to_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_serializable(x) for x in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, (float, int, str, bool)):
        return obj
    return str(obj)


def dump_tour_data(setup_mode: str = "fast", out_path: Path | None = None) -> Path:
    if out_path is None:
        out_path = Path(__file__).resolve().parents[1] / "docs" / "tour_data.json"
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    setup = build_calibration_setup(setup_mode)
    names = list(setup.names)
    theta = build_theta_from_x0(setup)
    P = apply_theta(setup.P_base, theta, names)
    P.collect_ge_trace = True

    print(f"Solving at theta:")
    for n, v in zip(names, theta):
        print(f"  {n:18s} = {v:.4f}")

    t0 = time.perf_counter()
    sol, P_out, p_eq = run_model_cp_dt(P, verbose=True)
    elapsed = time.perf_counter() - t0
    print(f"Solved in {elapsed:.1f}s. p_eq = {p_eq.tolist()}")

    print("Slicing policies...")
    policies = slice_policies(sol, P_out)
    print(f"  {len(policies['slices'])} slices.")

    moments_target = setup.targets
    moments_weight = setup.weights
    moments_model = extract_moment_dict(sol)

    # Build summary block (always there even if some moments fail)
    summary = {
        "TFR": float(2 * sol.mean_parity),
        "own_rate": float(sol.own_rate),
        "p_eq": p_eq.tolist(),
        "pop_share": sol.pop_share.tolist() if hasattr(sol.pop_share, "tolist") else list(sol.pop_share),
        "mean_parity_by_loc": sol.mean_parity_by_loc.tolist() if hasattr(sol.mean_parity_by_loc, "tolist") else list(sol.mean_parity_by_loc),
    }

    # b_grid lives on P_out
    b_grid = None
    if hasattr(P_out, "b_grid"):
        b_grid = np.asarray(P_out.b_grid).tolist()
    else:
        # rebuild it
        try:
            from dt_cp_model.utils import make_grid
            b_grid = make_grid(P_out).tolist()
        except Exception:
            pass

    H_own = np.asarray(P_out.H_own).tolist() if hasattr(P_out, "H_own") else None

    data = {
        "meta": {
            "setup_mode": setup_mode,
            "theta_names": names,
            "theta": theta.tolist(),
            "calibrated_overrides": CALIBRATED_OVERRIDES,
            "p_eq": p_eq.tolist(),
            "elapsed_sec": float(elapsed),
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "n_ge_iterations": len(getattr(sol, "ge_trace", []) or []),
        },
        "grids": {
            "b_grid": b_grid,
            "H_own": H_own,
            "J": int(P_out.J),
            "I": int(P_out.I),
            "n_parity": int(P_out.n_parity),
            "ncs": int(getattr(P_out, "n_child_states", policies["axes"]["ncs"])),
        },
        "policies": policies,
        "ge_trace": _to_serializable(getattr(sol, "ge_trace", [])),
        "summary": _to_serializable(summary),
        "moments": {
            "targets": _to_serializable(moments_target),
            "weights": _to_serializable(moments_weight),
            "model": _to_serializable(moments_model),
        },
    }

    with open(out_path, "w") as f:
        json.dump(data, f, separators=(",", ":"))
    size_kb = out_path.stat().st_size / 1024
    print(f"Wrote {out_path} ({size_kb:.1f} KB)")
    return out_path


if __name__ == "__main__":
    import sys
    mode = "fast"
    if len(sys.argv) > 1 and sys.argv[1] in ("fast", "benchmark"):
        mode = sys.argv[1]
    dump_tour_data(setup_mode=mode)
