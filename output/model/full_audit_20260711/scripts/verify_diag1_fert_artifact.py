#!/usr/bin/env python3
"""VERIFIER for DIAG-1: independent counterfactual at Nb=40.

Solve PE at the candidate theta, then re-run the production KFE
(forward_distribution_markov_income) twice:
  (a) baseline fert_probs from the solver;
  (b) fert_probs forced to [1,0,0] (stay childless) at nodes where the
      post-fertility-choice childless value V[...,0,0] <= -1e9 (all options
      carry the -1e10 infeasibility sentinel, so the softmax is ~uniform).
Compare mean_parity / TFR / childless_rate / total_births_kfe. Audit-only:
no production code touched; the mask lives in this script's copy of fp.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from audit_diag_packet import build_overrides, VALID_V  # noqa: E402

from intergen_housing_fertility.solver import (  # noqa: E402
    run_model_cp_dt,
    forward_distribution_markov_income,
    precompute_shared,
)

BASE = Path(
    "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/full_audit_20260711"
)
NB = 40


def main() -> None:
    theta = dict(json.loads((BASE / "repro/current_repro_exact.json").read_text())["theta"])
    overrides, _, _ = build_overrides(NB, 10, theta)
    overrides["solve_mode"] = "pe"
    overrides["p_fixed"] = np.array([0.6897624769313723])
    sol, P, _ = run_model_cp_dt(overrides, verbose=False)

    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    SD = precompute_shared(P, b_grid)
    p = np.asarray(sol.p_eq, dtype=float).reshape(-1)
    r = P.user_cost_rate * p

    V = np.asarray(sol.V)
    fp0 = np.asarray(sol.fert_probs).copy()

    # Baseline KFE rerun with the solver's own policies (should reproduce sol).
    g_b, st_b = forward_distribution_markov_income(
        sol.bp_pol, sol.hR_pol, sol.tenure_choice, sol.loc_probs, fp0,
        r, p, P, b_grid, SD, fast_stats=False, tenure_probs=sol.tenure_probs,
    )

    # Masked fert probs: stay childless where the childless post-choice value
    # is at the sentinel scale (all fertility options infeasible).
    fp1 = fp0.copy()
    npar = fp1.shape[-1]
    invalid_childless = V[:, :, :, :, :, 0, 0] <= VALID_V  # (Nb, nt, I, J, Nz)
    stay = np.zeros(npar)
    stay[0] = 1.0
    fp1[invalid_childless, :] = stay
    n_masked = int(np.sum(invalid_childless))

    g_m, st_m = forward_distribution_markov_income(
        sol.bp_pol, sol.hR_pol, sol.tenure_choice, sol.loc_probs, fp1,
        r, p, P, b_grid, SD, fast_stats=False, tenure_probs=sol.tenure_probs,
    )

    # Direct artifact-mass measurement at j=0 (pre-reshuffle entrants on
    # invalid childless nodes): reconstruct from masked run where the mass
    # stays at n=0.
    inv_j0_mass_baseline = float(
        np.sum((g_b[:, :, :, 0] * (V[:, :, :, 0] <= VALID_V))[..., 0, 0])
    )

    out = {
        "Nb": NB,
        "n_masked_node_rows": n_masked,
        "solver_sol": {
            "mean_parity": float(sol.mean_parity),
            "tfr": 2.0 * float(sol.mean_parity),
            "childless_rate": float(sol.parity_dist[0]),
            "total_births_kfe": float(getattr(sol, "total_births_kfe", np.nan)),
        },
        "baseline_rerun": {
            "mean_parity": float(st_b.mean_parity),
            "tfr": 2.0 * float(st_b.mean_parity),
            "childless_rate": float(st_b.parity_dist[0]),
            "total_births_kfe": float(getattr(st_b, "total_births_kfe", np.nan)),
        },
        "masked_infeasible_fert": {
            "mean_parity": float(st_m.mean_parity),
            "tfr": 2.0 * float(st_m.mean_parity),
            "childless_rate": float(st_m.parity_dist[0]),
            "total_births_kfe": float(getattr(st_m, "total_births_kfe", np.nan)),
        },
        "invalid_childless_j0_mass_baseline_postreshuffle_n0": inv_j0_mass_baseline,
        "artifact_birth_share_of_total": (
            1.0
            - float(getattr(st_m, "total_births_kfe", np.nan))
            / float(getattr(st_b, "total_births_kfe", np.nan))
        ),
    }
    (BASE / "verify_F1" / "verify_diag1_fert_artifact_nb40.json").write_text(
        json.dumps(out, indent=2)
    )
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
