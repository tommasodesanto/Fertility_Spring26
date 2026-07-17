#!/usr/bin/env python3
"""Audit the treated margin for parent-targeted credit relief.

The policy comparison showed only a tiny fertility response to raising the
financed share for new parents. This script separates three objects:

1. feasibility exposure of fertile childless renters under the baseline
   distribution;
2. policy-on-baseline-state changes in fertility and birth-plus-owner-entry
   probabilities;
3. full equilibrium moments after re-solving each policy case.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError:  # Figures are optional; CSV/JSON diagnostics are primary.
    plt = None


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from intergen_housing_fertility.calibration import (  # noqa: E402
    PERIOD_YEARS,
    base_overrides,
    extract_moments,
    jsonable,
)
from intergen_housing_fertility.local_panel import income_process_overrides  # noqa: E402
from intergen_housing_fertility.solver import get_phi_choice_tensor, income_at_state, run_model_cp_dt  # noqa: E402


DEFAULT_SUMMARY = (
    ROOT
    / "output/model/intergen_loss6_current_20260630/best.json"
)
DEFAULT_OUTDIR = ROOT / "output/model/intergen_parent_credit_margin_audit_20260630"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--J", type=int, default=17)
    parser.add_argument("--Nb", type=int, default=60)
    parser.add_argument("--n-house", type=int, default=5)
    parser.add_argument("--income-states", type=int, default=5)
    parser.add_argument("--max-iter-eq", type=int, default=10)
    parser.add_argument(
        "--fixed-price-from-record",
        action="store_true",
        help="Solve all cases in partial equilibrium at the record price.",
    )
    parser.add_argument(
        "--fixed-price",
        type=float,
        default=None,
        help="Solve all cases in partial equilibrium at this price.",
    )
    parser.add_argument("--skip-ltv100-all-child", action="store_true")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    source = load_source_record(args.summary)
    theta = normalize_theta(dict(source["theta"]))
    fixed_price = resolve_fixed_price(source, args)
    cases = default_cases(include_ltv100_all_child=not bool(args.skip_ltv100_all_child))

    solved: dict[str, dict[str, Any]] = {}
    for case in cases:
        t0 = time.perf_counter()
        sol, P, p_eq = solve_case(
            case,
            theta=theta,
            J=args.J,
            Nb=args.Nb,
            n_house=args.n_house,
            income_states=args.income_states,
            max_iter_eq=args.max_iter_eq,
            fixed_price=fixed_price,
        )
        moments = extract_moments(sol, P)
        solved[case["case"]] = {
            "case": case,
            "sol": sol,
            "P": P,
            "p_eq": p_eq,
            "moments": moments,
            "elapsed_sec": time.perf_counter() - t0,
        }
        print(
            f"{case['case']}: TFR={moments['tfr']:.4f}, "
            f"childless={moments['childless_rate']:.4f}, "
            f"own={moments['own_rate']:.4f}, "
            f"resid={getattr(sol, 'best_max_abs_rel_excess', np.nan):.2e}, "
            f"elapsed={solved[case['case']]['elapsed_sec']:.1f}s",
            flush=True,
        )

    baseline = solved["baseline"]
    state_weights = baseline_state_weights(baseline["sol"], baseline["P"])
    feasibility_rows = []
    response_rows = []
    response_by_age_rows = []
    for case_name, bundle in solved.items():
        feasibility_rows.extend(feasibility_audit(case_name, bundle["P"], bundle["p_eq"], state_weights))
        response = policy_response_on_baseline_states(
            case_name,
            bundle["sol"],
            bundle["P"],
            state_weights,
        )
        response_rows.append(response["summary"])
        response_by_age_rows.extend(response["by_age"])

    decomposition_payload = mechanism_decomposition(solved, state_weights)
    equilibrium_rows = [equilibrium_row(name, bundle) for name, bundle in solved.items()]
    write_csv(args.outdir / "equilibrium_moments.csv", equilibrium_rows)
    write_csv(args.outdir / "feasibility_exposure.csv", feasibility_rows)
    write_csv(args.outdir / "policy_on_baseline_state_response.csv", response_rows)
    write_csv(args.outdir / "policy_on_baseline_state_response_by_age.csv", response_by_age_rows)
    write_csv(args.outdir / "mechanism_decomposition.csv", decomposition_payload["summary"])
    write_csv(args.outdir / "susceptible_mass_bins.csv", decomposition_payload["bins"])
    write_csv(args.outdir / "value_path_unit_test.csv", decomposition_payload["unit_test"])
    state_dir = args.outdir / "state_decomposition"
    for case_name, rows in decomposition_payload["state_rows"].items():
        write_csv(state_dir / f"{case_name}.csv", rows)

    payload = {
        "source_record": source,
        "fixed_price": None if fixed_price is None else float(fixed_price),
        "cases": [case["case"] for case in cases],
        "equilibrium_moments": equilibrium_rows,
        "feasibility_exposure": feasibility_rows,
        "policy_on_baseline_state_response": response_rows,
        "policy_on_baseline_state_response_by_age": response_by_age_rows,
        "mechanism_decomposition": decomposition_payload["summary"],
        "susceptible_mass_bins": decomposition_payload["bins"],
        "value_path_unit_test": decomposition_payload["unit_test"],
    }
    write_json(args.outdir / "summary.json", payload)
    write_figures(args.outdir / "figures", feasibility_rows, response_rows, response_by_age_rows)
    write_readme(args.outdir, payload)


def load_source_record(path: Path) -> dict[str, Any]:
    raw = json.loads(path.read_text())
    if isinstance(raw, dict) and isinstance(raw.get("best"), dict):
        return raw["best"]
    if isinstance(raw, dict) and "theta" in raw:
        return raw
    raise ValueError(f"{path} must be a best-record JSON or a summary JSON with a 'best' record.")


def normalize_theta(theta: dict[str, Any]) -> dict[str, Any]:
    out = dict(theta)
    if "beta" not in out and "beta_annual" in out:
        out["beta"] = float(out["beta_annual"]) ** PERIOD_YEARS
    out.pop("beta_annual", None)
    return out


def resolve_fixed_price(source: dict[str, Any], args: argparse.Namespace) -> float | None:
    if args.fixed_price is not None:
        return float(args.fixed_price)
    if not bool(args.fixed_price_from_record):
        return None
    vals = source.get("p_eq", source.get("price", None))
    if vals is None:
        raise ValueError("--fixed-price-from-record requested but source record has no p_eq/price.")
    return float(np.asarray(vals, dtype=float).reshape(-1)[0])


def default_cases(*, include_ltv100_all_child: bool) -> list[dict[str, Any]]:
    cases = [
        {
            "case": "baseline",
            "label": "Baseline",
            "overrides": {},
        },
        {
            "case": "parent_ltv95_birth",
            "label": "Parent LTV 95, birth state",
            "overrides": {
                "parent_dp_waiver": True,
                "parent_dp_waiver_phi": 0.95,
                "parent_dp_waiver_birth_state_only": True,
            },
        },
        {
            "case": "parent_ltv100_birth",
            "label": "Parent LTV 100, birth state",
            "overrides": {
                "parent_dp_waiver": True,
                "parent_dp_waiver_phi": 1.00,
                "parent_dp_waiver_birth_state_only": True,
            },
        },
        {
            "case": "parent_ltv100_family_rungs_birth",
            "label": "Parent LTV 100, birth state, H>=6 only",
            "overrides": {
                "parent_dp_waiver": True,
                "parent_dp_waiver_phi": 1.00,
                "parent_dp_waiver_birth_state_only": True,
                "parent_dp_waiver_owner_rungs": np.array([3, 4, 5]),
            },
        },
    ]
    if include_ltv100_all_child:
        cases.append(
            {
                "case": "parent_ltv100_all_child",
                "label": "Parent LTV 100, all child-at-home states",
                "overrides": {
                    "parent_dp_waiver": True,
                    "parent_dp_waiver_phi": 1.00,
                    "parent_dp_waiver_birth_state_only": False,
                },
            }
        )
    return cases


def solve_case(
    case: dict[str, Any],
    *,
    theta: dict[str, float],
    J: int,
    Nb: int,
    n_house: int,
    income_states: int,
    max_iter_eq: int,
    fixed_price: float | None,
) -> tuple[Any, Any, Any]:
    overrides = {
        **base_overrides(J=J, Nb=Nb, n_house=n_house, max_iter_eq=max_iter_eq),
        **income_process_overrides(income_states),
        **theta,
        **dict(case["overrides"]),
    }
    if fixed_price is not None:
        overrides.update({"solve_mode": "pe", "p_fixed": np.array([float(fixed_price)])})
    return run_model_cp_dt(overrides, verbose=False)


def baseline_state_weights(sol: Any, P: Any) -> list[dict[str, float]]:
    g = np.asarray(sol.g, dtype=float)
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    z_grid = np.asarray(getattr(sol, "type_values", getattr(P, "z_grid", [1.0])), dtype=float).reshape(-1)
    rows: list[dict[str, float]] = []
    for j in fertile_age_indices(P):
        age = float(P.age_start + j * P.da)
        for i in range(int(P.I)):
            for zz, z_value in enumerate(z_grid):
                mass_b = g[:, 0, i, j, zz, 0, 0]
                nz = np.flatnonzero(mass_b > 1e-14)
                y = income_at_state(P, i, j, float(z_value))
                for ib in nz:
                    rows.append(
                        {
                            "b_index": int(ib),
                            "b": float(b_grid[ib]),
                            "j": int(j),
                            "age": age,
                            "i": int(i),
                            "z_index": int(zz),
                            "z": float(z_value),
                            "income": float(y),
                            "mass": float(mass_b[ib]),
                        }
                    )
    return rows


def fertile_age_indices(P: Any) -> list[int]:
    return [
        j
        for j in range(int(P.J))
        if (j + 1 >= int(P.A_f_start)) and (j + 1 <= int(P.A_f_end))
    ]


def family_need_rungs(P: Any, parity: int) -> np.ndarray:
    parity = int(parity)
    if str(getattr(P, "child_housing_spec", "jump_plus_linear")).lower() == "linear_only":
        need = float(P.h_bar_0 + P.h_bar_n * parity)
    else:
        need = float(P.h_bar_0 + P.h_bar_jump + P.h_bar_n * parity)
    h = np.asarray(P.H_own, dtype=float).reshape(-1)
    idx = np.flatnonzero(h >= need) + 1
    if idx.size == 0:
        idx = np.array([len(h)], dtype=int)
    return idx


def feasibility_audit(case_name: str, P: Any, p_eq: Any, states: list[dict[str, float]]) -> list[dict[str, Any]]:
    p = np.asarray(p_eq, dtype=float).reshape(-1)
    h = np.asarray(P.H_own, dtype=float).reshape(-1)
    phi = get_phi_choice_tensor(P)
    rows: list[dict[str, Any]] = []
    for parity in range(1, int(P.n_parity)):
        rung_sets = {
            "any_owner": np.arange(1, int(P.n_house) + 1),
            "family_need": family_need_rungs(P, parity),
        }
        for rung_scope, tenures in rung_sets.items():
            rows.append(feasibility_row(case_name, P, p, h, phi, states, parity, rung_scope, tenures))
    return rows


def feasibility_row(
    case_name: str,
    P: Any,
    p: np.ndarray,
    h: np.ndarray,
    phi: np.ndarray,
    states: list[dict[str, float]],
    parity: int,
    rung_scope: str,
    tenures: np.ndarray,
) -> dict[str, Any]:
    total_mass = sum(s["mass"] for s in states)
    dp_feasible = 0.0
    pti_feasible = 0.0
    both_feasible = 0.0
    dp_fail = 0.0
    pti_fail = 0.0
    dp_pass_pti_fail = 0.0
    dp_fail_pti_pass = 0.0
    near_5 = 0.0
    near_10 = 0.0
    near_25 = 0.0
    mean_min_dp_shortfall_income = 0.0
    mean_min_pti_ratio = 0.0

    for s in states:
        i = int(s["i"])
        y = max(float(s["income"]), 1e-12)
        b = float(s["b"])
        mass = float(s["mass"])
        dp_ok_any = False
        pti_ok_any = False
        both_ok_any = False
        min_shortfall = np.inf
        min_pti_ratio = np.inf
        for ten in tenures:
            ph = float(phi[i, int(ten), parity, 1])
            hcost = float(p[i] * h[int(ten) - 1])
            dp_threshold = (1.0 - ph) * hcost
            payment = (float(P.q) * ph + float(P.tau_H)) * hcost
            dp_ok = b >= dp_threshold
            pti_ok = payment <= float(P.pti_limit) * y
            dp_ok_any = dp_ok_any or dp_ok
            pti_ok_any = pti_ok_any or pti_ok
            both_ok_any = both_ok_any or (dp_ok and pti_ok)
            min_shortfall = min(min_shortfall, max(dp_threshold - b, 0.0))
            min_pti_ratio = min(min_pti_ratio, payment / max(float(P.pti_limit) * y, 1e-12))
        dp_feasible += mass * float(dp_ok_any)
        pti_feasible += mass * float(pti_ok_any)
        both_feasible += mass * float(both_ok_any)
        dp_fail += mass * float(not dp_ok_any)
        pti_fail += mass * float(not pti_ok_any)
        dp_pass_pti_fail += mass * float(dp_ok_any and not pti_ok_any)
        dp_fail_pti_pass += mass * float((not dp_ok_any) and pti_ok_any)
        short_y = min_shortfall / y
        near_5 += mass * float(0.0 < short_y <= 0.05)
        near_10 += mass * float(0.0 < short_y <= 0.10)
        near_25 += mass * float(0.0 < short_y <= 0.25)
        mean_min_dp_shortfall_income += mass * short_y
        mean_min_pti_ratio += mass * min_pti_ratio

    den = max(total_mass, 1e-12)
    return {
        "case": case_name,
        "parity": int(parity),
        "rung_scope": rung_scope,
        "rungs": json.dumps([int(x) for x in tenures]),
        "rung_services": json.dumps([float(h[int(x) - 1]) for x in tenures]),
        "mass": total_mass,
        "dp_feasible_share": dp_feasible / den,
        "pti_feasible_share": pti_feasible / den,
        "both_feasible_share": both_feasible / den,
        "dp_fail_share": dp_fail / den,
        "pti_fail_share": pti_fail / den,
        "dp_pass_pti_fail_share": dp_pass_pti_fail / den,
        "dp_fail_pti_pass_share": dp_fail_pti_pass / den,
        "near_dp_shortfall_le_5pct_income_share": near_5 / den,
        "near_dp_shortfall_le_10pct_income_share": near_10 / den,
        "near_dp_shortfall_le_25pct_income_share": near_25 / den,
        "mean_min_dp_shortfall_over_income": mean_min_dp_shortfall_income / den,
        "mean_min_pti_ratio_to_limit": mean_min_pti_ratio / den,
    }


def policy_response_on_baseline_states(case_name: str, sol: Any, P: Any, states: list[dict[str, float]]) -> dict[str, Any]:
    fp = np.asarray(sol.fert_probs, dtype=float)
    tp = getattr(sol, "tenure_probs", None)
    tc = getattr(sol, "tenure_choice", None)
    rows_by_age: dict[float, dict[str, float]] = {}
    totals = init_response_accumulator(case_name, age=np.nan)
    for s in states:
        age = float(s["age"])
        if age not in rows_by_age:
            rows_by_age[age] = init_response_accumulator(case_name, age=age)
        ib = int(s["b_index"])
        j = int(s["j"])
        i = int(s["i"])
        zz = int(s["z_index"])
        mass = float(s["mass"])
        for acc in (totals, rows_by_age[age]):
            acc["mass"] += mass
        fert = fp[ib, 0, i, j, zz, :]
        expected_n = float(np.sum(np.arange(int(P.n_parity)) * fert))
        birth_prob = float(np.sum(fert[1:]))
        for acc in (totals, rows_by_age[age]):
            acc["expected_completed_children_sum"] += mass * expected_n
            acc["birth_prob_sum"] += mass * birth_prob
            acc["childless_prob_sum"] += mass * float(fert[0])
            if int(P.n_parity) > 2:
                acc["two_plus_prob_sum"] += mass * float(np.sum(fert[2:]))

        for parity in range(1, int(P.n_parity)):
            birth_mass = mass * float(fert[parity])
            any_owner = owner_entry_prob(sol, P, tp, tc, ib, i, j, zz, parity, "any_owner")
            family_owner = owner_entry_prob(sol, P, tp, tc, ib, i, j, zz, parity, "family_need")
            for acc in (totals, rows_by_age[age]):
                acc["birth_and_owner_sum"] += birth_mass * any_owner
                acc["birth_and_family_owner_sum"] += birth_mass * family_owner
                acc["owner_if_birth_sum"] += birth_mass * any_owner
                acc["family_owner_if_birth_sum"] += birth_mass * family_owner
                acc["birth_mass_sum"] += birth_mass

    return {
        "summary": finalize_response_accumulator(totals),
        "by_age": [finalize_response_accumulator(rows_by_age[age]) for age in sorted(rows_by_age)],
    }


def init_response_accumulator(case_name: str, *, age: float) -> dict[str, float | str]:
    return {
        "case": case_name,
        "age": age,
        "mass": 0.0,
        "expected_completed_children_sum": 0.0,
        "birth_prob_sum": 0.0,
        "childless_prob_sum": 0.0,
        "two_plus_prob_sum": 0.0,
        "birth_and_owner_sum": 0.0,
        "birth_and_family_owner_sum": 0.0,
        "owner_if_birth_sum": 0.0,
        "family_owner_if_birth_sum": 0.0,
        "birth_mass_sum": 0.0,
    }


def finalize_response_accumulator(acc: dict[str, float | str]) -> dict[str, Any]:
    mass = max(float(acc["mass"]), 1e-12)
    birth_mass = max(float(acc["birth_mass_sum"]), 1e-12)
    return {
        "case": acc["case"],
        "age": float(acc["age"]),
        "mass": float(acc["mass"]),
        "expected_completed_children": float(acc["expected_completed_children_sum"]) / mass,
        "tfr_equivalent_on_states": 2.0 * float(acc["expected_completed_children_sum"]) / mass,
        "birth_probability": float(acc["birth_prob_sum"]) / mass,
        "childless_probability": float(acc["childless_prob_sum"]) / mass,
        "two_plus_probability": float(acc["two_plus_prob_sum"]) / mass,
        "birth_and_owner_entry_probability": float(acc["birth_and_owner_sum"]) / mass,
        "birth_and_family_owner_entry_probability": float(acc["birth_and_family_owner_sum"]) / mass,
        "owner_entry_probability_conditional_on_birth": float(acc["owner_if_birth_sum"]) / birth_mass,
        "family_owner_entry_probability_conditional_on_birth": float(acc["family_owner_if_birth_sum"]) / birth_mass,
    }


def owner_entry_prob(
    sol: Any,
    P: Any,
    tenure_probs: Any,
    tenure_choice: Any,
    ib: int,
    i: int,
    j: int,
    zz: int,
    parity: int,
    rung_scope: str,
) -> float:
    if rung_scope == "family_need":
        tenures = family_need_rungs(P, parity)
    else:
        tenures = np.arange(1, int(P.n_house) + 1)
    if tenure_probs is not None:
        pr = np.asarray(tenure_probs)[ib, 0, i, j, zz, parity, 1, :]
        return float(np.sum(pr[tenures]))
    choice = int(np.asarray(tenure_choice)[ib, 0, i, j, zz, parity, 1])
    return float(choice in set(int(t) for t in tenures))


def equilibrium_row(case_name: str, bundle: dict[str, Any]) -> dict[str, Any]:
    sol = bundle["sol"]
    moments = bundle["moments"]
    return {
        "case": case_name,
        "elapsed_sec": float(bundle["elapsed_sec"]),
        "market_residual": float(getattr(sol, "best_max_abs_rel_excess", np.nan)),
        "owner_asset_price": float(np.asarray(sol.owner_asset_price).reshape(-1)[0]),
        "owner_user_cost": float(np.asarray(sol.owner_user_cost).reshape(-1)[0]),
        "tfr": float(moments.get("tfr", np.nan)),
        "childless_rate": float(moments.get("childless_rate", np.nan)),
        "own_rate": float(moments.get("own_rate", np.nan)),
        "own_rate_2534": float(moments.get("own_rate_2534", np.nan)),
        "own_rate_3544": float(moments.get("own_rate_3544", np.nan)),
        "own_family_gap": float(moments.get("own_family_gap", np.nan)),
        "housing_increment_0to1": float(moments.get("housing_increment_0to1", np.nan)),
        "housing_increment_1to2": float(moments.get("housing_increment_1to2", np.nan)),
        "mean_age_first_birth": float(moments.get("mean_age_first_birth", np.nan)),
    }


def mechanism_decomposition(solved: dict[str, dict[str, Any]], states: list[dict[str, float]]) -> dict[str, Any]:
    if "baseline" not in solved:
        raise ValueError("mechanism_decomposition requires a 'baseline' case.")
    baseline = solved["baseline"]
    state_rows_by_case: dict[str, list[dict[str, Any]]] = {}
    summary_rows: list[dict[str, Any]] = []
    bin_rows: list[dict[str, Any]] = []
    unit_rows: list[dict[str, Any]] = []
    selected_unit_states = select_unit_test_states(baseline["P"], baseline["p_eq"], states)

    for case_name, bundle in solved.items():
        if case_name == "baseline":
            continue
        rows = []
        for state in states:
            row = state_mechanism_row(case_name, baseline, bundle, state)
            rows.append(row)
        state_rows_by_case[case_name] = rows
        for group_name, predicate in [
            ("all_fertile_childless_renters", lambda r: True),
            ("young25_35_childless_renters", lambda r: 25.0 <= float(r["age"]) <= 35.0),
        ]:
            group_rows = [r for r in rows if predicate(r)]
            summary_rows.append(summarize_mechanism_rows(case_name, group_name, group_rows))
            bin_rows.extend(susceptible_bin_rows(case_name, group_name, group_rows))
        for state in selected_unit_states:
            unit_rows.append(unit_test_row(case_name, baseline, bundle, state))

    return {
        "summary": summary_rows,
        "bins": bin_rows,
        "unit_test": unit_rows,
        "state_rows": state_rows_by_case,
    }


def state_mechanism_row(
    case_name: str,
    baseline: dict[str, Any],
    policy: dict[str, Any],
    state: dict[str, float],
) -> dict[str, Any]:
    bsol, bP = baseline["sol"], baseline["P"]
    psol, pP = policy["sol"], policy["P"]
    ib = int(state["b_index"])
    j = int(state["j"])
    i = int(state["i"])
    zz = int(state["z_index"])
    mass = float(state["mass"])

    bvals = fertility_values_for_state(bsol, bP, ib, i, j, zz)
    pvals = fertility_values_for_state(psol, pP, ib, i, j, zz)
    b_owner_any = owner_probability_for_birth(bsol, bP, ib, i, j, zz, 1, "any_owner")
    p_owner_any = owner_probability_for_birth(psol, pP, ib, i, j, zz, 1, "any_owner")
    b_owner_family = owner_probability_for_birth(bsol, bP, ib, i, j, zz, 1, "family_need")
    p_owner_family = owner_probability_for_birth(psol, pP, ib, i, j, zz, 1, "family_need")
    b_starter_owner = max(b_owner_any - b_owner_family, 0.0)
    p_starter_owner = max(p_owner_any - p_owner_family, 0.0)
    b_h_birth = expected_housing_services_for_birth(bsol, bP, ib, i, j, zz, 1)
    p_h_birth = expected_housing_services_for_birth(psol, pP, ib, i, j, zz, 1)
    b_h_ge6 = probability_services_at_least_for_birth(bsol, bP, ib, i, j, zz, 1, 6.0)
    p_h_ge6 = probability_services_at_least_for_birth(psol, pP, ib, i, j, zz, 1, 6.0)
    shortfall = min_family_dp_shortfall_over_income(bP, baseline["p_eq"], state, parity=1)

    d_birth_value = pvals["birth_logsum_value"] - bvals["birth_logsum_value"]
    d_no_value = pvals["no_child_value"] - bvals["no_child_value"]
    d_rel_value = d_birth_value - d_no_value
    d_birth_prob = pvals["birth_probability"] - bvals["birth_probability"]
    d_expected_n = pvals["expected_children"] - bvals["expected_children"]
    d_owner_any = p_owner_any - b_owner_any
    d_owner_family = p_owner_family - b_owner_family
    d_starter_owner = p_starter_owner - b_starter_owner
    d_h_birth = p_h_birth - b_h_birth
    d_h_ge6 = p_h_ge6 - b_h_ge6
    actual_derivative = bvals["birth_probability"] * (1.0 - bvals["birth_probability"]) / max(float(bP.kappa_fert), 1e-12)
    max_derivative = 0.25 / max(float(bP.kappa_fert), 1e-12)
    linearized_actual_dp = actual_derivative * d_rel_value
    upper_abs_dp = max_derivative * abs(d_rel_value)

    return {
        "case": case_name,
        "age": float(state["age"]),
        "j": j,
        "b_index": ib,
        "b": float(state["b"]),
        "income": float(state["income"]),
        "z_index": zz,
        "z": float(state["z"]),
        "mass": mass,
        "dp_shortfall_family_parity1_over_income": shortfall,
        "baseline_birth_probability": bvals["birth_probability"],
        "policy_birth_probability": pvals["birth_probability"],
        "delta_birth_probability": d_birth_prob,
        "baseline_expected_children": bvals["expected_children"],
        "policy_expected_children": pvals["expected_children"],
        "delta_expected_children": d_expected_n,
        "baseline_logit_derivative_pb_1mpb": bvals["birth_probability"] * (1.0 - bvals["birth_probability"]),
        "baseline_logit_derivative_per_value": actual_derivative,
        "max_logit_derivative_per_value": max_derivative,
        "linearized_actual_derivative_delta_birth_probability": linearized_actual_dp,
        "upper_bound_abs_delta_birth_probability_max_derivative": upper_abs_dp,
        "upper_bound_abs_delta_birth_probability_max_derivative_capped": min(upper_abs_dp, 1.0),
        "baseline_no_child_value": bvals["no_child_value"],
        "policy_no_child_value": pvals["no_child_value"],
        "delta_no_child_value": d_no_value,
        "baseline_birth_logsum_value": bvals["birth_logsum_value"],
        "policy_birth_logsum_value": pvals["birth_logsum_value"],
        "delta_birth_logsum_value": d_birth_value,
        "baseline_birth_minus_no_value": bvals["birth_logsum_value"] - bvals["no_child_value"],
        "policy_birth_minus_no_value": pvals["birth_logsum_value"] - pvals["no_child_value"],
        "delta_birth_minus_no_value": d_rel_value,
        "baseline_owner_prob_if_birth_any": b_owner_any,
        "policy_owner_prob_if_birth_any": p_owner_any,
        "delta_owner_prob_if_birth_any": d_owner_any,
        "baseline_owner_prob_if_birth_family": b_owner_family,
        "policy_owner_prob_if_birth_family": p_owner_family,
        "delta_owner_prob_if_birth_family": d_owner_family,
        "baseline_starter_owner_prob_if_birth": b_starter_owner,
        "policy_starter_owner_prob_if_birth": p_starter_owner,
        "delta_starter_owner_prob_if_birth": d_starter_owner,
        "baseline_housing_services_if_birth": b_h_birth,
        "policy_housing_services_if_birth": p_h_birth,
        "delta_housing_services_if_birth": d_h_birth,
        "baseline_prob_h_ge6_if_birth": b_h_ge6,
        "policy_prob_h_ge6_if_birth": p_h_ge6,
        "delta_prob_h_ge6_if_birth": d_h_ge6,
        "ownership_mover": int(abs(d_owner_any) > 1e-7 or abs(d_owner_family) > 1e-7),
        "fertility_mover": int(abs(d_birth_prob) > 1e-6 or abs(d_expected_n) > 1e-6),
        "both_mover": int((abs(d_owner_any) > 1e-7 or abs(d_owner_family) > 1e-7) and (abs(d_birth_prob) > 1e-6 or abs(d_expected_n) > 1e-6)),
    }


def fertility_values_for_state(sol: Any, P: Any, ib: int, i: int, j: int, zz: int) -> dict[str, float]:
    fp = np.asarray(sol.fert_probs, dtype=float)
    fv = np.asarray(sol.fert_value, dtype=float)
    if fp.ndim == 6:
        probs = fp[ib, 0, i, j, zz, :]
        logsum_value = float(fv[ib, 0, i, j, zz])
    else:
        probs = fp[ib, 0, i, j, :]
        logsum_value = float(fv[ib, 0, i, j])
    probs = np.asarray(probs, dtype=float).reshape(-1)
    eps = 1e-300
    kappa = float(P.kappa_fert)
    no_prob = float(probs[0]) if probs.size > 0 else 1.0
    birth_prob = float(np.sum(probs[1:])) if probs.size > 1 else 0.0
    expected_n = float(np.sum(np.arange(probs.size) * probs))
    no_value = logsum_value + kappa * math.log(max(no_prob, eps))
    birth_logsum_value = logsum_value + kappa * math.log(max(birth_prob, eps))
    return {
        "no_child_probability": no_prob,
        "birth_probability": birth_prob,
        "expected_children": expected_n,
        "no_child_value": no_value,
        "birth_logsum_value": birth_logsum_value,
    }


def tenure_probability_vector(sol: Any, P: Any, ib: int, i: int, j: int, zz: int, parity: int, cs: int) -> np.ndarray:
    nt = 1 + int(P.n_house)
    tp = getattr(sol, "tenure_probs", None)
    if tp is not None:
        arr = np.asarray(tp, dtype=float)
        if arr.ndim == 8:
            pr = arr[ib, 0, i, j, zz, parity, cs, :]
        elif arr.ndim == 7:
            pr = arr[ib, 0, i, j, parity, cs, :]
        else:
            pr = np.zeros(nt)
        pr = np.asarray(pr, dtype=float).reshape(-1)
        if pr.size == nt and np.isfinite(pr).all() and float(np.sum(pr)) > 0.5:
            return pr / max(float(np.sum(pr)), 1e-12)
    tc = np.asarray(getattr(sol, "tenure_choice"), dtype=int)
    if tc.ndim == 7:
        choice = int(tc[ib, 0, i, j, zz, parity, cs])
    else:
        choice = int(tc[ib, 0, i, j, parity, cs])
    out = np.zeros(nt)
    out[min(max(choice, 0), nt - 1)] = 1.0
    return out


def owner_probability_for_birth(sol: Any, P: Any, ib: int, i: int, j: int, zz: int, parity: int, rung_scope: str) -> float:
    pr = tenure_probability_vector(sol, P, ib, i, j, zz, parity, 1)
    if rung_scope == "family_need":
        tenures = family_need_rungs(P, parity)
    else:
        tenures = np.arange(1, int(P.n_house) + 1)
    return float(np.sum(pr[[int(t) for t in tenures]]))


def expected_housing_services_for_birth(sol: Any, P: Any, ib: int, i: int, j: int, zz: int, parity: int) -> float:
    pr = tenure_probability_vector(sol, P, ib, i, j, zz, parity, 1)
    services = np.zeros(1 + int(P.n_house))
    hR = np.asarray(sol.hR_pol, dtype=float)
    if hR.ndim == 7:
        services[0] = float(hR[ib, 0, i, j, zz, parity, 1])
    else:
        services[0] = float(hR[ib, 0, i, j, parity, 1])
    services[1:] = np.asarray(P.H_own, dtype=float).reshape(-1)
    return float(np.sum(pr * services))


def probability_services_at_least_for_birth(
    sol: Any,
    P: Any,
    ib: int,
    i: int,
    j: int,
    zz: int,
    parity: int,
    threshold: float,
) -> float:
    pr = tenure_probability_vector(sol, P, ib, i, j, zz, parity, 1)
    services = np.zeros(1 + int(P.n_house))
    hR = np.asarray(sol.hR_pol, dtype=float)
    if hR.ndim == 7:
        services[0] = float(hR[ib, 0, i, j, zz, parity, 1])
    else:
        services[0] = float(hR[ib, 0, i, j, parity, 1])
    services[1:] = np.asarray(P.H_own, dtype=float).reshape(-1)
    return float(np.sum(pr[services >= float(threshold)]))


def min_family_dp_shortfall_over_income(P: Any, p_eq: Any, state: dict[str, float], *, parity: int) -> float:
    p = np.asarray(p_eq, dtype=float).reshape(-1)
    phi = get_phi_choice_tensor(P)
    h = np.asarray(P.H_own, dtype=float).reshape(-1)
    i = int(state["i"])
    b = float(state["b"])
    y = max(float(state["income"]), 1e-12)
    shortfall = np.inf
    for ten in family_need_rungs(P, parity):
        ten = int(ten)
        hcost = float(p[i] * h[ten - 1])
        dp = (1.0 - float(phi[i, ten, parity, 1])) * hcost
        shortfall = min(shortfall, dp - b)
    return float(shortfall / y)


def summarize_mechanism_rows(case_name: str, group_name: str, rows: list[dict[str, Any]]) -> dict[str, Any]:
    mass = sum(float(r["mass"]) for r in rows)
    den = max(mass, 1e-12)

    def wavg(key: str) -> float:
        return float(sum(float(r["mass"]) * float(r[key]) for r in rows) / den)

    return {
        "case": case_name,
        "group": group_name,
        "mass": mass,
        "baseline_birth_probability": wavg("baseline_birth_probability"),
        "policy_birth_probability": wavg("policy_birth_probability"),
        "delta_birth_probability": wavg("delta_birth_probability"),
        "delta_expected_children": wavg("delta_expected_children"),
        "delta_tfr_equivalent_on_states": 2.0 * wavg("delta_expected_children"),
        "baseline_logit_derivative_pb_1mpb": wavg("baseline_logit_derivative_pb_1mpb"),
        "baseline_logit_derivative_per_value": wavg("baseline_logit_derivative_per_value"),
        "linearized_actual_derivative_delta_tfr_equivalent": 2.0 * wavg("linearized_actual_derivative_delta_birth_probability"),
        "upper_bound_tfr_equivalent_max_derivative": 2.0 * wavg("upper_bound_abs_delta_birth_probability_max_derivative"),
        "upper_bound_tfr_equivalent_max_derivative_capped": 2.0 * wavg("upper_bound_abs_delta_birth_probability_max_derivative_capped"),
        "delta_birth_minus_no_value": wavg("delta_birth_minus_no_value"),
        "delta_no_child_value": wavg("delta_no_child_value"),
        "delta_birth_logsum_value": wavg("delta_birth_logsum_value"),
        "delta_owner_prob_if_birth_any": wavg("delta_owner_prob_if_birth_any"),
        "delta_owner_prob_if_birth_family": wavg("delta_owner_prob_if_birth_family"),
        "delta_starter_owner_prob_if_birth": wavg("delta_starter_owner_prob_if_birth"),
        "delta_housing_services_if_birth": wavg("delta_housing_services_if_birth"),
        "delta_prob_h_ge6_if_birth": wavg("delta_prob_h_ge6_if_birth"),
        "ownership_mover_mass_share": sum(float(r["mass"]) * int(r["ownership_mover"]) for r in rows) / den,
        "fertility_mover_mass_share": sum(float(r["mass"]) * int(r["fertility_mover"]) for r in rows) / den,
        "both_mover_mass_share": sum(float(r["mass"]) * int(r["both_mover"]) for r in rows) / den,
    }


def susceptible_bin_rows(case_name: str, group_name: str, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    bins = [
        ("feasible", -np.inf, 0.0),
        ("shortfall_0_5pct_income", 0.0, 0.05),
        ("shortfall_5_10pct_income", 0.05, 0.10),
        ("shortfall_10_25pct_income", 0.10, 0.25),
        ("shortfall_gt_25pct_income", 0.25, np.inf),
    ]
    out = []
    total_mass = sum(float(r["mass"]) for r in rows)
    for label, lo, hi in bins:
        sub = [
            r
            for r in rows
            if float(r["dp_shortfall_family_parity1_over_income"]) > lo
            and float(r["dp_shortfall_family_parity1_over_income"]) <= hi
        ]
        row = summarize_mechanism_rows(case_name, group_name, sub) if sub else empty_bin_summary(case_name, group_name)
        row["dp_shortfall_bin"] = label
        row["mass_share_in_group"] = row["mass"] / max(total_mass, 1e-12)
        out.append(row)
    return out


def empty_bin_summary(case_name: str, group_name: str) -> dict[str, Any]:
    keys = summarize_mechanism_rows(case_name, group_name, [])
    keys["mass"] = 0.0
    return keys


def select_unit_test_states(P: Any, p_eq: Any, states: list[dict[str, float]], n: int = 5) -> list[dict[str, float]]:
    young = [s for s in states if 25.0 <= float(s["age"]) <= 35.0 and float(s["mass"]) > 1e-14]
    if not young:
        young = [s for s in states if float(s["mass"]) > 1e-14]
    scored = []
    for state in young:
        shortfall = min_family_dp_shortfall_over_income(P, p_eq, state, parity=1)
        scored.append((abs(shortfall), -float(state["mass"]), state))
    scored.sort(key=lambda x: (x[0], x[1]))
    return [s for _, _, s in scored[:n]]


def unit_test_row(case_name: str, baseline: dict[str, Any], policy: dict[str, Any], state: dict[str, float]) -> dict[str, Any]:
    row = state_mechanism_row(case_name, baseline, policy, state)
    bP = baseline["P"]
    pP = policy["P"]
    i = int(state["i"])
    p_base = np.asarray(baseline["p_eq"], dtype=float).reshape(-1)
    p_pol = np.asarray(policy["p_eq"], dtype=float).reshape(-1)
    b_phi = get_phi_choice_tensor(bP)
    p_phi = get_phi_choice_tensor(pP)
    ten = int(family_need_rungs(bP, 1)[0])
    h_base = float(np.asarray(bP.H_own, dtype=float)[ten - 1])
    h_pol = float(np.asarray(pP.H_own, dtype=float)[ten - 1])
    row.update(
        {
            "tested_owner_rung": ten,
            "tested_owner_rooms_baseline": h_base,
            "tested_owner_rooms_policy": h_pol,
            "baseline_phi_birth": float(b_phi[i, ten, 1, 1]),
            "policy_phi_birth": float(p_phi[i, ten, 1, 1]),
            "baseline_dp_birth": float((1.0 - b_phi[i, ten, 1, 1]) * p_base[i] * h_base),
            "policy_dp_birth": float((1.0 - p_phi[i, ten, 1, 1]) * p_pol[i] * h_pol),
            "baseline_borrowing_floor_birth": float(-b_phi[i, ten, 1, 1] * p_base[i] * h_base),
            "policy_borrowing_floor_birth": float(-p_phi[i, ten, 1, 1] * p_pol[i] * h_pol),
        }
    )
    return row


def write_figures(
    outdir: Path,
    feasibility_rows: list[dict[str, Any]],
    response_rows: list[dict[str, Any]],
    response_by_age_rows: list[dict[str, Any]],
) -> None:
    if plt is None:
        return
    outdir.mkdir(parents=True, exist_ok=True)
    fam1 = [r for r in feasibility_rows if r["parity"] == 1 and r["rung_scope"] == "family_need"]
    cases = [r["case"] for r in fam1]
    x = np.arange(len(cases))
    fig, ax = plt.subplots(figsize=(9.5, 4.5))
    width = 0.22
    for offset, key, label in [
        (-width, "dp_feasible_share", "DP feasible"),
        (0.0, "pti_feasible_share", "PTI feasible"),
        (width, "both_feasible_share", "Both feasible"),
    ]:
        ax.bar(x + offset, [r[key] for r in fam1], width=width, label=label)
    ax.set_title("Baseline fertile childless renters: one-child family-rung feasibility")
    ax.set_ylabel("share of baseline-state mass")
    ax.set_xticks(x, cases, rotation=25, ha="right")
    ax.set_ylim(0.0, 1.05)
    ax.grid(axis="y", alpha=0.2)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "01_feasibility_family_rung_one_child.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.5, 4.5))
    ax.bar(x, [r["dp_pass_pti_fail_share"] for r in fam1])
    ax.set_title("DP-pass but PTI-fail mass after policy")
    ax.set_ylabel("share of baseline-state mass")
    ax.set_xticks(x, cases, rotation=25, ha="right")
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir / "02_dp_pass_pti_fail_family_rung_one_child.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    cases = [r["case"] for r in response_rows]
    x = np.arange(len(cases))
    axes[0].bar(x, [r["tfr_equivalent_on_states"] for r in response_rows])
    axes[0].set_title("Fertility on baseline states")
    axes[0].set_ylabel("TFR-equivalent")
    axes[1].bar(x, [r["birth_and_family_owner_entry_probability"] for r in response_rows])
    axes[1].set_title("Birth and family-rung owner entry")
    axes[1].set_ylabel("probability")
    for ax in axes:
        ax.set_xticks(x, cases, rotation=25, ha="right")
        ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir / "03_policy_response_on_baseline_states.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    for case in cases:
        rows = [r for r in response_by_age_rows if r["case"] == case]
        rows = sorted(rows, key=lambda r: r["age"])
        axes[0].plot([r["age"] for r in rows], [r["tfr_equivalent_on_states"] for r in rows], marker="o", label=case)
        axes[1].plot(
            [r["age"] for r in rows],
            [r["birth_and_family_owner_entry_probability"] for r in rows],
            marker="o",
            label=case,
        )
    axes[0].set_title("Fertility on baseline states by age")
    axes[0].set_ylabel("TFR-equivalent")
    axes[1].set_title("Birth and family-rung owner entry by age")
    axes[1].set_ylabel("probability")
    for ax in axes:
        ax.set_xlabel("age")
        ax.grid(alpha=0.2)
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "04_policy_response_by_age.png", dpi=180)
    plt.close(fig)


def write_readme(outdir: Path, payload: dict[str, Any]) -> None:
    rows = payload["equilibrium_moments"]
    response = payload["policy_on_baseline_state_response"]
    lines = [
        "# Parent Credit Treated-Margin Audit",
        "",
        "This diagnostic uses the final global-DE toy best and asks why parent-targeted credit relief barely moves fertility.",
        "",
        "The audit reports feasibility exposure for fertile childless renters under the baseline distribution, then compares policy functions on the same baseline states.",
        "",
        "## Equilibrium Moments",
        "",
        "| Case | TFR | Childless | Prime own | Own 25-34 | Own 35-44 | Family own gap | Residual |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        lines.append(
            f"| `{r['case']}` | {r['tfr']:.3f} | {r['childless_rate']:.3f} | "
            f"{r['own_rate']:.3f} | {r['own_rate_2534']:.3f} | {r['own_rate_3544']:.3f} | "
            f"{r['own_family_gap']:.3f} | {r['market_residual']:.2e} |"
        )
    lines.extend(
        [
            "",
            "## Policy On Baseline Fertile Childless Renter States",
            "",
            "| Case | TFR-equivalent | Birth prob | Birth and owner entry | Birth and family-rung owner entry | Owner entry if birth |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for r in response:
        lines.append(
            f"| `{r['case']}` | {r['tfr_equivalent_on_states']:.3f} | {r['birth_probability']:.3f} | "
            f"{r['birth_and_owner_entry_probability']:.3f} | "
            f"{r['birth_and_family_owner_entry_probability']:.3f} | "
            f"{r['owner_entry_probability_conditional_on_birth']:.3f} |"
        )
    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `equilibrium_moments.csv`",
            "- `feasibility_exposure.csv`",
            "- `policy_on_baseline_state_response.csv`",
            "- `policy_on_baseline_state_response_by_age.csv`",
            "- `figures/`",
        ]
    )
    (outdir / "README.md").write_text("\n".join(lines) + "\n")


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(jsonable(obj), indent=2, sort_keys=True))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
