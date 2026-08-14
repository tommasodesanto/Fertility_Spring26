#!/usr/bin/env python3
"""Build the small quantitative incidence audit for the population closure.

The script reads the maintained E5b reproductive-closure audit.  It does not
resolve or recalibrate the household model.  At the audited benchmark it asks:

1. What fixed outside inflow is required for a stationary open population?
2. How do local supply and reproduction shifts affect price and population
   once stationary scale is endogenous?

The local calculation combines

    S = M / (E - rho B),
    d log(B/E) / d log(p) = -eta_R,
    d log(H^s / D) / d log(p) = eta_F.

For an outside-origin share s=M/E, baseline stationarity implies
rho B/E=1-s and lambda=(1-s)/s.  The formulas deliberately hold the required
entrant flow E fixed, matching the maintained stationary solver.  A direct
reproduction shift also holds per-household housing demand fixed; this makes
that row a transparent incidence diagnostic, not a policy counterfactual.
"""
from __future__ import annotations

import argparse
import csv
import importlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "output/model/closed_reproductive_closure_audit_20260811/summary.json"
WINNER_SOURCE = (
    ROOT
    / "output/model/closed_reproductive_closure_audit_20260811/source/"
    "incumbent_current_contract_preflight_summary.json"
)
OUTDIR = ROOT / "output/model/population_closure_update"
POLICY_CASES: tuple[tuple[str, str, dict[str, Any]], ...] = (
    ("baseline", "Baseline", {}),
    (
        "birth_grant_A0p4_Hge6",
        "Birth grant, family-size homes",
        {
            "birth_entry_grant": True,
            "birth_entry_grant_amount": 0.4,
            "birth_entry_grant_owner_rungs": np.array([3, 4, 5]),
            "birth_entry_grant_locations": np.array([], dtype=int),
        },
    ),
    (
        "tax2_grant_A0p4_Hge6",
        "Property tax plus birth grant",
        {
            "birth_entry_grant": True,
            "birth_entry_grant_amount": 0.4,
            "birth_entry_grant_owner_rungs": np.array([3, 4, 5]),
            "birth_entry_grant_locations": np.array([], dtype=int),
            "tau_H": 0.08,
        },
    ),
)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def solve_legacy_policy_cases() -> list[dict[str, Any]]:
    """Resolve the retained household block under normalized and renewal scale."""

    os.environ["E3_L4"] = "1"
    os.environ["E5"] = "1"
    os.environ["E3_TFR_TOP_BIN_WEIGHT"] = "3.602359422009"
    for name in ("E5_MATURATION_REPAIR", "E5F", "E6A", "E6B", "E6C"):
        os.environ.pop(name, None)
    model_root = ROOT / "code/model"
    if str(model_root) not in sys.path:
        sys.path.insert(0, str(model_root))
    from intergen_eqscale_seq_optimized import run_e1_chain

    chain = importlib.reload(run_e1_chain)
    chain.load_runtime()
    source = json.loads(WINNER_SOURCE.read_text(encoding="utf-8"))
    winner = source["best_tight"]
    theta = {str(key): float(value) for key, value in winner["theta"].items()}
    runtime = argparse.Namespace(J=17, Nb=120, max_iter_eq=40, tol_eq=2.5e-5)
    base = chain.common_overrides(runtime)
    base.update(theta)

    normalized: dict[str, tuple[Any, Any, np.ndarray, float]] = {}
    for case, _, policy in POLICY_CASES:
        started = time.perf_counter()
        solution, parameters, price = chain.run_model_cp_dt({**base, **policy}, verbose=False)
        normalized[case] = (
            solution,
            parameters,
            np.asarray(price, dtype=float).reshape(-1),
            time.perf_counter() - started,
        )

    baseline = normalized["baseline"][0]
    outside_flow = float(baseline.entry_rate) - float(baseline.entrants_mature_total)
    if outside_flow <= 0.0:
        raise RuntimeError(f"Maintained full-retention closure requires positive outside inflow: {outside_flow}")

    rows: list[dict[str, Any]] = []
    for case, label, policy in POLICY_CASES:
        fixed, _, fixed_price, fixed_seconds = normalized[case]
        renewal_overrides = {
            **base,
            **policy,
            "population_closure": "renewal_valve",
            "renewal_retention": 1.0,
            "outside_entry_flow": outside_flow,
            "renewal_calibrate_outside_flow": False,
        }
        started = time.perf_counter()
        renewal, parameters, renewal_price = chain.run_model_cp_dt(renewal_overrides, verbose=False)
        renewal_seconds = time.perf_counter() - started
        scale = renewal.accounting_scale
        moments = chain.extract_moments(renewal, parameters)
        rows.append(
            {
                "case": case,
                "label": label,
                "fixed_population_price": float(fixed_price[0]),
                "renewal_price": float(np.asarray(renewal_price).reshape(-1)[0]),
                "renewal_population_scale": float(scale.scale_factor),
                "tfr_per_household": float(moments["tfr"]),
                "current_births_per_unit_scale": float(renewal.total_births_kfe),
                "total_current_births": float(scale.scale_factor * renewal.total_births_kfe),
                "mature_local_entrants_per_unit_scale": float(renewal.entrants_mature_total),
                "effective_reproduction_B_over_E": float(renewal.entrants_mature_total / renewal.entry_rate),
                "total_mature_local_entrants": float(scale.scale_factor * renewal.entrants_mature_total),
                "outside_entry_flow": float(scale.outside_entry_flow),
                "renewal_identity_residual": float(scale.stationary_scale_residual),
                "renewal_market_residual": float(renewal.best_max_abs_rel_excess),
                "fixed_population_solve_seconds": float(fixed_seconds),
                "renewal_solve_seconds": float(renewal_seconds),
            }
        )

    base_row = rows[0]
    if abs(base_row["renewal_population_scale"] - 1.0) > 1e-6:
        raise RuntimeError(f"Renewal baseline failed S=1 nesting: {base_row}")
    if abs(base_row["renewal_price"] - base_row["fixed_population_price"]) > 2e-6:
        raise RuntimeError(f"Renewal baseline failed price nesting: {base_row}")
    for row in rows:
        row["tfr_change_percent"] = 100.0 * (row["tfr_per_household"] / base_row["tfr_per_household"] - 1.0)
        row["fixed_population_price_change_percent"] = 100.0 * (
            row["fixed_population_price"] / base_row["fixed_population_price"] - 1.0
        )
        row["population_change_percent"] = 100.0 * (row["renewal_population_scale"] - 1.0)
        row["renewal_price_change_percent"] = 100.0 * (row["renewal_price"] / base_row["renewal_price"] - 1.0)
        row["total_current_births_change_percent"] = 100.0 * (
            row["total_current_births"] / base_row["total_current_births"] - 1.0
        )
        row["total_mature_entrants_change_percent"] = 100.0 * (
            row["total_mature_local_entrants"] / base_row["total_mature_local_entrants"] - 1.0
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--solve-policies",
        action="store_true",
        help="Also resolve the retained legacy policy cases under the renewal closure.",
    )
    args = parser.parse_args()
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    baseline = source["baseline"]
    grid = source["maintained_grid"]

    entry = float(baseline["entry_households_per_adult"])
    mature = float(baseline["mature_entrant_households_per_adult"])
    reproduction = float(baseline["reproduction_ratio_B_over_E"])
    eta_reproduction = -float(grid["local_elasticity_B_over_E_to_user_cost"])
    eta_capacity = float(grid["local_elasticity_scale_map_to_user_cost"])
    minimum_outside_share = 1.0 - reproduction
    minimum_outside_flow = entry - mature

    shares = [minimum_outside_share, 0.20, 0.25, 0.50]
    rows: list[dict[str, Any]] = []
    for outside_share in shares:
        if outside_share + 1e-12 < minimum_outside_share:
            continue
        retained_native_fraction = (1.0 - outside_share) / reproduction
        renewal_multiplier = (1.0 - outside_share) / outside_share
        demographic_price_slope = renewal_multiplier * eta_reproduction
        denominator = eta_capacity + demographic_price_slope
        rows.append(
            {
                "outside_share_of_required_entrants": outside_share,
                "implied_native_retention": retained_native_fraction,
                "renewal_multiplier": renewal_multiplier,
                "price_elasticity_to_supply_intercept": -1.0 / denominator,
                "population_elasticity_to_supply_intercept": demographic_price_slope / denominator,
                "price_elasticity_to_direct_mature_birth_flow": renewal_multiplier / denominator,
                "population_elasticity_to_direct_mature_birth_flow": eta_capacity * renewal_multiplier / denominator,
            }
        )

    summary = {
        "status": "local incidence diagnostic; not a transition and not a policy run",
        "source": str(SOURCE),
        "period_years": 4,
        "entry_households_per_adult": entry,
        "mature_local_entrants_per_adult": mature,
        "effective_reproduction_B_over_E": reproduction,
        "minimum_outside_share_for_rho_at_most_one": minimum_outside_share,
        "minimum_outside_flow_per_adult_per_four_years": minimum_outside_flow,
        "simple_annualized_minimum_outside_flow_per_adult": minimum_outside_flow / 4.0,
        "local_reproduction_elasticity_to_housing_cost": -eta_reproduction,
        "local_housing_capacity_elasticity_to_housing_cost": eta_capacity,
        "closed_positive_root_found": bool(grid["raw_reproduction_root_found"]),
        "interpretation": {
            "open_stationary": (
                "The maintained normalized point is algebraically an open stationary equilibrium "
                "if fixed outside inflow and native retention satisfy the renewal identity."
            ),
            "transition": (
                "This calculation does not establish that the normalized KFE is a point on a "
                "calendar-time transition; that requires an observed initial cohort distribution, "
                "a child pipeline, migration by age, and a path of rents and asset prices."
            ),
        },
    }

    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_csv(OUTDIR / "local_incidence.csv", rows)
    policy_rows = solve_legacy_policy_cases() if args.solve_policies else []
    if policy_rows:
        write_csv(OUTDIR / "legacy_policy_comparison.csv", policy_rows)

    first = rows[0]
    readme = f"""# Population-closure update

This is a small algebraic audit built from the maintained E5b fixed-price
schedule. It is not a recalibration, policy run, or calendar-time transition.

## Benchmark verdict

- Effective reproduction is `{reproduction:.6f}`.
- With full native retention, maintaining scale requires outside inflow equal
  to `{minimum_outside_share:.3%}` of required entrants, or
  `{minimum_outside_flow:.6f}` per adult every four years.
- The maintained fertility response is too flat to generate a positive closed
  reproductive steady state over the audited price range.
- Consequently the current calibration can be represented as an **open
  stationary point**, but it has not yet been shown to be a point on a
  calendar-time transition.

## Local incidence at the minimum-inflow closure

For a one-percent proportional supply expansion, price changes by
`{first['price_elasticity_to_supply_intercept']:.3f}` percent and population
by `{first['population_elasticity_to_supply_intercept']:.3f}` percent. For a
one-percent direct increase in mature locally born entrants, abstracting from
its direct housing-demand effect, price changes by
`{first['price_elasticity_to_direct_mature_birth_flow']:.3f}` percent and
population by
`{first['population_elasticity_to_direct_mature_birth_flow']:.3f}` percent.
These are elasticities, so the values are responses to a one-percent—not a
one-unit—change.

The contrast is the economic result: with weak housing-to-fertility feedback,
extra supply mainly lowers price, whereas a direct fertility response is
amplified into stationary population scale. Both magnitudes depend strongly on
the outside-inflow anchor; `local_incidence.csv` reports the sensitivity.

## Retained policy cases

`legacy_policy_comparison.csv` is written when the driver is run with
`--solve-policies`. It resolves the retained shared-clock household model under
both normalized demand and the fixed-inflow renewal closure. These legacy
policy cases are diagnostics, not promoted counterfactuals: the empirical
rooms target remains under hold, and the tax-and-grant exercises do not settle
the new transition question.

## Regeneration

From the repository root:

```bash
code/model/.venv/bin/python code/model/tools/build_population_closure_update.py --solve-policies
```
"""
    (OUTDIR / "README.md").write_text(readme, encoding="utf-8")


if __name__ == "__main__":
    main()
