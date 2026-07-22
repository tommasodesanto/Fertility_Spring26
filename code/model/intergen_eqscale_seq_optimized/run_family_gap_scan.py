"""T4 mechanism probe: family-ownership-gap reachability over (d_jump, d_a).

Fixed-theta GE evaluations at the collected E2 winner with only the two
share-tilt parameters varied on a 5x5 grid. Asks whether the eqscale
preference architecture can produce a substantial parent/nonparent ownership
gap anywhere in the admissible box (data target 0.168; E2 winner delivers
0.042). Runs under the E2 architecture (n_parity=3 bins), matching the
documented failure point. One tight solve per cell; no optimizer; no
production import changed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from intergen_eqscale_seq_optimized.calibration import extract_moments, get_target_system
from intergen_eqscale_seq_optimized.e1_profile import e1_overrides
from intergen_eqscale_seq_optimized.solver import run_model_cp_dt

E2_THETA = {
    "H0": 7.155359847249007,
    "alpha_cons": 0.6594894135201201,
    "beta": 0.9845586356367257,
    "chi": 1.0759980336241914,
    "delta_alpha": 0.02164113100132686,
    "delta_alpha_jump": 0.03884456373416531,
    "gamma_e": 0.20249058184087665,
    "kappa_fert": 1.7977544685025562,
    "psi_child": -0.3000089859271815,
    "tenure_choice_kappa": 0.023120596155365748,
    "theta0": 0.14801064516847956,
    "theta1": 0.28514054670326805,
    "theta_n": 0.0,
}
E2_PRICE = 0.8515432012274374

# First entries reproduce the E2 winner exactly (continuity anchor).
DELTA_JUMP_GRID = [0.03884456373416531, 0.10, 0.16, 0.21, 0.25]
DELTA_ALPHA_GRID = [0.02164113100132686, 0.08, 0.14, 0.20, 0.25]

KEEP_MOMENTS = [
    "own_family_gap",
    "own_rate",
    "own_rate_nonparents_3055",
    "own_rate_newparents_3055",
    "tfr",
    "childless_rate",
    "housing_increment_0to1",
    "prime30_55_parent_3plus_minus_1to2_mean_rooms",
    "prime30_55_childless_renter_mean_rooms",
    "aggregate_mean_occupied_rooms_18_85",
    "market_residual",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell", type=int, required=True, help="1..25")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.cell <= 25:
        raise SystemExit("cell must be in 1..25")
    dj = DELTA_JUMP_GRID[(args.cell - 1) // 5]
    da = DELTA_ALPHA_GRID[(args.cell - 1) % 5]

    overrides = {
        **e1_overrides(tight=True, optimized=True),
        **E2_THETA,
        "delta_alpha_jump": dj,
        "delta_alpha": da,
        "p_init_override": np.array([E2_PRICE]),
    }
    sol, P, _ = run_model_cp_dt(overrides, verbose=False)
    moments = extract_moments(sol, P)
    system = get_target_system("candidate_replacement_income_disciplined_v1")
    targets, weights = system.targets_dict(), system.weights_dict()
    loss = sum(
        float(weights[k]) * (float(moments.get(k, float("nan"))) - float(v)) ** 2
        for k, v in targets.items()
    )
    record = {
        "cell": args.cell,
        "delta_alpha_jump": dj,
        "delta_alpha": da,
        "loss_15_moment_fixed_theta": loss,
        "moments": {k: float(moments.get(k, float("nan"))) for k in KEEP_MOMENTS},
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / f"cell_{args.cell:02d}.json"
    out.write_text(json.dumps(record, indent=1))
    print(json.dumps(record))


if __name__ == "__main__":
    main()
