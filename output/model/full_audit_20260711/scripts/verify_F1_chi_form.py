#!/usr/bin/env python3
"""Verifier check for finding F1 (owner surplus chi*(H-hbar) vs paper chi*H-hbar).

Exploits the existing owner_h_bar_scale override: setting it to 1/chi turns the
code form chi*(H - scale*hbar) into exactly the paper form chi*H - hbar.
Runs the current_bound_best theta at tiny Nb under both forms and reports the
per-moment differences. Read-only w.r.t. production code.
"""
import json
import time
from pathlib import Path

import numpy as np

import intergen_housing_fertility.local_panel as lp
from intergen_housing_fertility.production_profile import (
    PRODUCTION_J,
    PRODUCTION_TARGET_SET,
    production_profile_overrides,
)

MATCHED_ANNUAL_RHO = 0.9601845894041878
MATCHED_ANNUAL_INNOVATION_SD = 0.06453733259357768
ROOMS_TARGET = 5.779970481941968
PERIOD_YEARS = 4.0
NB = 40
MAX_ITER_EQ = 10

THETA = {
    "H0": 9.999672059746388,
    "alpha_cons": 0.6731064483929385,
    "beta": 0.8411423842686492,
    "c_bar_0": 1.28,
    "c_bar_n": 0.4559787014242588,
    "chi": 1.15,
    "h_bar_0": 1.0,
    "h_bar_jump": 1.4758939442666474,
    "h_bar_n": 0.984106767291167,
    "kappa_fert": 1.8906304000000005,
    "psi_child": 0.26907999999999993,
    "tenure_choice_kappa": 0.0,
    "theta0": 0.1318301350511569,
    "theta_n": 0.7679202774006243,
}

OUTDIR = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/full_audit_20260711/verify_F1")


def build_env(owner_h_bar_scale=None):
    targets, weights = lp.get_target_set(PRODUCTION_TARGET_SET)
    targets["aggregate_mean_occupied_rooms_18_85"] = ROOMS_TARGET
    weights["aggregate_mean_occupied_rooms_18_85"] = 6.0
    income = lp.income_process_overrides(
        5, "rouwenhorst", MATCHED_ANNUAL_INNOVATION_SD, MATCHED_ANNUAL_RHO
    )
    extra = production_profile_overrides()
    extra.update(
        {
            "q": (1.0 + 0.02) ** PERIOD_YEARS - 1.0,
            "delta": 1.0 - (1.0 - 0.011) ** PERIOD_YEARS,
            "eta_supply": np.array([1.75]),
            "normalize_bequest_utility": True,
            "max_iter_eq": MAX_ITER_EQ,
        }
    )
    if owner_h_bar_scale is not None:
        extra["owner_h_bar_scale"] = float(owner_h_bar_scale)
    return targets, weights, income, extra


def run(label, owner_h_bar_scale=None):
    targets, weights, income, extra = build_env(owner_h_bar_scale)
    t0 = time.perf_counter()
    rec = lp.run_local_panel_case(
        0, {"label": label, "theta": dict(THETA)}, PRODUCTION_J, NB, 5,
        MAX_ITER_EQ, income, targets, weights, extra,
    )
    rec["wall"] = time.perf_counter() - t0
    return rec, targets, weights


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rec_code, targets, weights = run("code_form")           # scale=1 default: chi*(H-hbar)
    rec_paper, _, _ = run("paper_form", owner_h_bar_scale=1.0 / THETA["chi"])  # chi*H-hbar

    rows = []
    m_c = rec_code.get("moments") or {}
    m_p = rec_paper.get("moments") or {}
    for name in sorted(targets):
        mc = m_c.get(name)
        mp = m_p.get(name)
        rows.append({
            "moment": name,
            "target": float(targets[name]),
            "weight": float(weights[name]),
            "code_form": None if mc is None else float(mc),
            "paper_form": None if mp is None else float(mp),
            "diff_paper_minus_code": None if (mc is None or mp is None) else float(mp) - float(mc),
        })
    out = {
        "nb": NB,
        "note": "Nb=40 tiny grid; levels not comparable to Nb=120, differences are the object",
        "loss_code_form": float(rec_code.get("rank_loss", float("nan"))),
        "loss_paper_form": float(rec_paper.get("rank_loss", float("nan"))),
        "residual_code": float(rec_code.get("market_residual", float("nan"))),
        "residual_paper": float(rec_paper.get("market_residual", float("nan"))),
        "strict_code": rec_code.get("strict_converged"),
        "strict_paper": rec_paper.get("strict_converged"),
        "wall_code": rec_code.get("wall"),
        "wall_paper": rec_paper.get("wall"),
        "fit": rows,
    }
    (OUTDIR / "F1_chi_form_moment_diffs.json").write_text(json.dumps(out, indent=2))
    for r in rows:
        d = r["diff_paper_minus_code"]
        print(f"{r['moment']:.<52} code={r['code_form']!s:>12} paper={r['paper_form']!s:>12} diff={d!s}")
    print(f"loss code={out['loss_code_form']:.4f} paper={out['loss_paper_form']:.4f}")


if __name__ == "__main__":
    main()
