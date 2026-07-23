"""Fertility-frontier reachability probe under literal parity (L4).

Grid over (psi_child, kappa_fert) x gamma_e with all other coordinates fixed
at the certified E3 winner, L4 conventions on. Asks whether any region
delivers completed fertility ~1.918 AND childlessness ~0.188 jointly — i.e.
whether the two-parameter fertility block plus a real child cost can
generate the required bimodality (19% childless, parents averaging ~2.4
births) under honest income risk, or whether a structural change (per-parity
psi, permanent taste heterogeneity) is needed. One tight GE solve per cell.

The v2 design replaces the scanned linear scale with imposed power, square-root,
and linear-gamma-0.5 forms under the Floden-Linde x HSV after-tax income process.
It pre-registers whether a concave imposed scale can jointly reach completed
fertility 1.918 and childlessness 0.188, with linear_g05 isolating the income change.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from intergen_eqscale_seq_optimized.calibration import extract_moments
from intergen_eqscale_seq_optimized.e1_profile import e1_overrides
from intergen_eqscale_seq_optimized.externals import flhsv_income_overrides
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
FORM_GRID = ("power", "sqrt", "linear_g05")
PSI_GRID_V3 = [0.5, 1.0, 1.5, 2.0, 3.0, 4.5, 6.0]
DJUMP_GRID_V3 = [0.087, 0.15, 0.22, 0.30, 0.38, 0.46, 0.55]
KCFG_GRID_V3 = [
    ("shared", 15.9109, None),
    ("split_lo", 3.0, 0.3),
    ("split_hi", 10.0, 0.3),
]

KEEP = [
    "tfr", "mean_completed_fertility", "childless_rate",
    "parity_share_0", "parity_share_1", "parity_share_2plus", "parity_share_3plus",
    "parity_progression_1to2_flow", "parity_progression_2to3_flow",
    "mean_age_first_birth", "share_first_births_age30plus",
    "childless_chosen_45", "childless_clock_45",
    "own_rate", "own_family_gap", "market_residual",
]


def cell_spec(cell: int, design: str) -> tuple[float, float, float | str | tuple[str, float, float | None]]:
    """Return the form-major grid coordinates for one one-indexed frontier cell."""
    third_grid: tuple[float, ...] | tuple[str, ...]
    if design == "v3_entry_noise":
        ncell = len(KCFG_GRID_V3) * len(PSI_GRID_V3) * len(DJUMP_GRID_V3)
        if not 1 <= cell <= ncell:
            raise ValueError(f"cell must be in 1..{ncell}")
        kcfg_idx, rem = divmod(cell - 1, len(PSI_GRID_V3) * len(DJUMP_GRID_V3))
        psi_idx, djump_idx = divmod(rem, len(DJUMP_GRID_V3))
        return PSI_GRID_V3[psi_idx], DJUMP_GRID_V3[djump_idx], KCFG_GRID_V3[kcfg_idx]
    if design == "v1":
        third_grid = tuple(GAMMA_GRID)
    elif design == "v2_imposed_scale":
        third_grid = FORM_GRID
    else:
        raise ValueError(f"unknown frontier design: {design}")
    ncell = len(PSI_GRID) * len(KAPPA_GRID) * len(third_grid)
    if not 1 <= cell <= ncell:
        raise ValueError(f"cell must be in 1..{ncell}")
    third_idx, rem = divmod(cell - 1, len(PSI_GRID) * len(KAPPA_GRID))
    psi_idx, kappa_idx = divmod(rem, len(KAPPA_GRID))
    return PSI_GRID[psi_idx], KAPPA_GRID[kappa_idx], third_grid[third_idx]


def cell_overrides(cell: int, design: str = "v1") -> dict:
    """Build one fixed-coordinate frontier cell without solving the model."""
    psi, kap, third = cell_spec(cell, design)
    if design == "v3_entry_noise":
        kcfg, kappa_fert, kappa_fert_continuation = third
        overrides = {
            **e1_overrides(tight=True, optimized=True),
            **E3_THETA,
            **L4,
            **flhsv_income_overrides(),
            "eqscale_form": "power",
            "gamma_e": 0.0,
            "psi_child": psi,
            "delta_alpha_jump": kap,
            "kappa_fert": kappa_fert,
            "p_init_override": np.array([E3_PRICE_GUESS]),
        }
        if kappa_fert_continuation is not None:
            overrides["kappa_fert_continuation"] = kappa_fert_continuation
        return overrides
    if design == "v1":
        return {
            **e1_overrides(tight=True, optimized=True),
            **E3_THETA,
            **L4,
            "psi_child": psi,
            "kappa_fert": kap,
            "gamma_e": third,
            "p_init_override": np.array([E3_PRICE_GUESS]),
        }
    form = third
    form_overrides = (
        {"eqscale_form": form, "gamma_e": 0.0}
        if form in {"power", "sqrt"}
        else {"eqscale_form": "linear", "gamma_e": 0.5}
    )
    return {
        **e1_overrides(tight=True, optimized=True),
        **E3_THETA,
        **L4,
        **flhsv_income_overrides(),
        "psi_child": psi,
        "kappa_fert": kap,
        **form_overrides,
        "p_init_override": np.array([E3_PRICE_GUESS]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell", type=int, required=True, help="1..147")
    parser.add_argument("--design", choices=("v1", "v2_imposed_scale", "v3_entry_noise"), default="v1")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    try:
        psi, kap, third = cell_spec(args.cell, args.design)
        overrides = cell_overrides(args.cell, args.design)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    sol, P, _ = run_model_cp_dt(overrides, verbose=False)
    moments = extract_moments(sol, P)
    record = {
        "cell": args.cell,
        "psi_child": psi,
        "kappa_fert": kap,
        "gamma_e": overrides["gamma_e"],
        "moments": {k: float(moments.get(k, float("nan"))) for k in KEEP},
        "joint_distance": float(
            ((float(moments.get("tfr", float("nan"))) - 1.918) / 1.918) ** 2
            + ((float(moments.get("childless_rate", float("nan"))) - 0.188) / 0.188) ** 2
        ),
    }
    if args.design == "v2_imposed_scale":
        record.update(
            {
                "design": args.design,
                "form": third,
                "income": "flhsv",
                "eqscale_form": overrides["eqscale_form"],
                "gamma_e": overrides["gamma_e"],
            }
        )
    elif args.design == "v3_entry_noise":
        kcfg, kappa_fert, kappa_fert_continuation = third
        record.update(
            {
                "design": args.design,
                "form": "power",
                "income": "flhsv",
                "eqscale_form": overrides["eqscale_form"],
                "gamma_e": overrides["gamma_e"],
                "kcfg": kcfg,
                "delta_alpha_jump": kap,
                "kappa_fert": kappa_fert,
                "kappa_fert_continuation": kappa_fert_continuation,
            }
        )
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / f"cell_{args.cell:03d}.json").write_text(json.dumps(record, indent=1))
    print(json.dumps(record))


if __name__ == "__main__":
    main()
