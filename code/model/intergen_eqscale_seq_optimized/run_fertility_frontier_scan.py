"""Fertility-frontier reachability probe under literal parity (L4).

Grid over (psi_child, kappa_fert) x gamma_e with all other coordinates fixed
at the certified E3 winner, L4 conventions on. Asks whether any region
delivers completed fertility ~1.918 AND childlessness ~0.188 jointly — i.e.
whether the two-parameter fertility block plus a real child cost can
generate the required bimodality (19% childless, parents averaging ~2.4
births) under honest income risk, or whether a structural change (per-parity
psi, permanent taste heterogeneity) is needed. One tight GE solve per cell.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from intergen_eqscale_seq_optimized.calibration import extract_moments
from intergen_eqscale_seq_optimized.e1_profile import e1_overrides
from intergen_eqscale_seq_optimized.solver import run_model_cp_dt

E3_THETA = {
    "H0": 7.026, "alpha_cons": 0.6979, "beta": 0.955, "chi": 1.0676,
    "delta_alpha": 0.0236, "delta_alpha_jump": 0.0713, "gamma_e": 0.0085,
    "kappa_fert": 15.9109, "psi_child": 2.7233, "tenure_choice_kappa": 0.0173,
    "theta0": 0.1678, "theta1": 0.0959, "theta_n": 0.0,
}
E3_PRICE_GUESS = 0.85

L4 = {
    "n_parity": 4,
    "fertility_units": "literal_topcode",
    "tfr_top_bin_weight": 3.4,
    "entrant_conversion_factor": 0.5,
    "child_bin_high_cutoff": 3,
}

PSI_GRID = [-0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.7233]
KAPPA_GRID = [0.5, 1.0, 1.8, 3.0, 6.0, 10.0, 15.9109]
GAMMA_GRID = [0.0085, 0.2, 0.5]

KEEP = [
    "tfr", "mean_completed_fertility", "childless_rate",
    "parity_share_0", "parity_share_1", "parity_share_2plus", "parity_share_3plus",
    "parity_progression_1to2_flow", "parity_progression_2to3_flow",
    "mean_age_first_birth", "share_first_births_age30plus",
    "childless_chosen_45", "childless_clock_45",
    "own_rate", "own_family_gap", "market_residual",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell", type=int, required=True, help="1..147")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    ncell = len(PSI_GRID) * len(KAPPA_GRID) * len(GAMMA_GRID)
    if not 1 <= args.cell <= ncell:
        raise SystemExit(f"cell must be in 1..{ncell}")
    idx = args.cell - 1
    ig, rem = divmod(idx, len(PSI_GRID) * len(KAPPA_GRID))
    ip, ik = divmod(rem, len(KAPPA_GRID))
    psi, kap, gam = PSI_GRID[ip], KAPPA_GRID[ik], GAMMA_GRID[ig]

    overrides = {
        **e1_overrides(tight=True, optimized=True),
        **E3_THETA,
        **L4,
        "psi_child": psi,
        "kappa_fert": kap,
        "gamma_e": gam,
        "p_init_override": np.array([E3_PRICE_GUESS]),
    }
    sol, P, _ = run_model_cp_dt(overrides, verbose=False)
    moments = extract_moments(sol, P)
    record = {
        "cell": args.cell,
        "psi_child": psi,
        "kappa_fert": kap,
        "gamma_e": gam,
        "moments": {k: float(moments.get(k, float("nan"))) for k in KEEP},
        "joint_distance": float(
            ((float(moments.get("tfr", float("nan"))) - 1.918) / 1.918) ** 2
            + ((float(moments.get("childless_rate", float("nan"))) - 0.188) / 0.188) ** 2
        ),
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / f"cell_{args.cell:03d}.json").write_text(json.dumps(record, indent=1))
    print(json.dumps(record))


if __name__ == "__main__":
    main()
