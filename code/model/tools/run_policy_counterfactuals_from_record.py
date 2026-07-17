#!/usr/bin/env python3
"""Run policy counterfactuals from a saved direct-geometry record.

This is deliberately separate from the diagnostic-note builders. It rebuilds
the model from one JSON record, checks that the baseline re-solve is close to
the saved moments, then solves policy variants with the outside option active.

Counterfactual closure:
    1. solve the record under ``outside_option_benchmark_normalized``;
    2. recover the benchmark outside value and outside-origin entrant pool;
    3. solve baseline and policies under ``accounting_scale_prices`` holding
       those outside objects fixed, so GE prices and scale move.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Callable

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from dt_cp_model.direct_calibration import (  # noqa: E402
    DIRECT_GEOMETRY_NAMES,
    OUTSIDE_VALUE_NAME,
    RENEWAL_FLOW_NAME,
    build_direct_calibration_setup,
)
from dt_cp_model.objective import extract_moments  # noqa: E402
from dt_cp_model.parameters import apply_overrides, asdict, setup_parameters  # noqa: E402
from dt_cp_model.solver import run_model_cp_dt  # noqa: E402
from dt_cp_model.theta import apply_theta  # noqa: E402


DEFAULT_RECORD = (
    ROOT
    / "output/model/reduced_target_overnight_20260527/records/hR8_default_best.json"
)
DEFAULT_OUTDIR = ROOT / "output/model/policy_counterfactuals_20260528"

DEFAULT_CASES = [
    "parent_phi_relief",
    "center_xi_relaxed",
]

CHECK_KEYS = [
    "tfr",
    "childless_rate",
    "own_rate",
    "own_gradient",
    "own_family_gap",
    "housing_increment_0to1",
    "housing_increment_1to2",
    "center_share_nonparents",
    "center_share_newparents",
    "old_age_own_rate",
    "old_age_parent_childless_gap",
    "inv_pop_C",
    "inv_rent_ratio",
]

SUMMARY_KEYS = [
    "tfr",
    "childless_rate",
    "own_rate",
    "own_gradient",
    "own_family_gap",
    "housing_increment_0to1",
    "housing_increment_1to2",
    "young_liquid_wealth_to_income",
    "center_share_nonparents",
    "center_share_newparents",
    "migration_rate",
    "old_age_own_rate",
    "old_age_parent_childless_gap",
]


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    record = json.loads(args.record.read_text())
    theta = np.asarray(record["theta"], dtype=float)

    setup = build_setup(record)
    print(f"Loaded record {args.record}", flush=True)
    print("Solving saved record under benchmark-normalized outside-option closure...", flush=True)
    normalized = solve_record_case(
        "baseline_normalized_resolve",
        "Baseline re-solve under benchmark-normalized outside option",
        record,
        setup,
        theta,
        modifier=no_change,
        p_init=record.get("p_eq"),
        max_iter_eq=args.normalized_baseline_max_iter_eq,
    )
    write_json(args.outdir / "baseline_normalized_resolve.json", normalized)

    check_rows = baseline_check_rows(record.get("moments", {}), normalized["moments"])
    write_csv(args.outdir / "baseline_resolve_check.csv", check_rows)
    max_gap = max((abs(float(row["gap"])) for row in check_rows), default=0.0)
    if max_gap > args.baseline_tol:
        message = (
            f"baseline re-solve differs from saved record; max abs gap={max_gap:.6g}. "
            f"See {args.outdir / 'baseline_resolve_check.csv'}"
        )
        if args.baseline_mismatch == "fail":
            raise RuntimeError(message)
        if args.baseline_mismatch == "warn":
            print(f"WARNING: {message}", file=sys.stderr, flush=True)

    scale = SimpleNamespace(**normalized["scale"])
    outside_value = float(scale.outside_value)
    outside_entry_flow = float(scale.outside_entry_flow)
    print(
        "Fixed outside objects for CFs: "
        f"outside_value={outside_value:.6f}, outside_entry_flow={outside_entry_flow:.6f}",
        flush=True,
    )

    fixed_modifier = fixed_outside_modifier(outside_value, outside_entry_flow, float(scale.kappa_entry))
    cases = build_cases(args, fixed_modifier)

    results = []
    for case_id, description, modifier in cases:
        print(f"Solving {case_id}: {description}", flush=True)
        result = solve_record_case(
            case_id,
            description,
            record,
            setup,
            theta,
            modifier=modifier,
            p_init=record.get("p_eq"),
            max_iter_eq=args.policy_max_iter_eq if case_id != "baseline_fixed_outside" else args.fixed_baseline_max_iter_eq,
        )
        results.append(result)
        write_json(args.outdir / "partial_results.json", results)
        write_csv(args.outdir / "partial_summary.csv", summary_rows(results))
        m = result["moments"]
        eq = result["equilibrium"]
        sc = result["scale"]
        print(
            f"  {case_id}: TFR={m.get('tfr', math.nan):.3f}, "
            f"own={m.get('own_rate', math.nan):.3f}, "
            f"popC={eq.get('pop_share_C', math.nan):.3f}, "
            f"rr={eq.get('rent_ratio_C_over_P', math.nan):.3f}, "
            f"scale={sc.get('scale_factor', math.nan):.3f}",
            flush=True,
        )

    write_json(args.outdir / "counterfactual_results.json", results)
    write_json(args.outdir / "benchmark_outside_objects.json", normalized["scale"])
    rows = summary_rows(results)
    write_csv(args.outdir / "summary.csv", rows)
    write_csv(args.outdir / "equilibrium.csv", equilibrium_rows(results))
    write_markdown_summary(args.outdir / "summary.md", rows, args.record, normalized["scale"])
    write_latex_table(args.outdir / "summary_table.tex", rows)
    write_summary_plot(args.outdir / "summary_effects.png", rows)
    write_json(args.outdir / "policy_metadata.json", policy_metadata(args))
    write_validation_note(args.outdir / "validation_note.md", results, args.record, normalized["scale"], policy_metadata(args))
    print(f"Wrote {args.outdir / 'summary.csv'}", flush=True)
    print(f"Wrote {args.outdir / 'equilibrium.csv'}", flush=True)
    print(f"Wrote {args.outdir / 'summary_table.tex'}", flush=True)
    print(f"Wrote {args.outdir / 'summary_effects.png'}", flush=True)
    print(f"Wrote {args.outdir / 'validation_note.md'}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run policy counterfactuals from a saved calibration record.")
    parser.add_argument("--record", type=Path, default=DEFAULT_RECORD)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--baseline-tol", type=float, default=1e-3)
    parser.add_argument("--normalized-baseline-max-iter-eq", type=int, default=None)
    parser.add_argument("--fixed-baseline-max-iter-eq", type=int, default=60)
    parser.add_argument(
        "--policy-max-iter-eq",
        type=int,
        default=80,
        help="Maximum GE iterations for fixed-outside baseline and policy cases.",
    )
    parser.add_argument(
        "--baseline-mismatch",
        choices=("warn", "fail", "ignore"),
        default="warn",
        help="What to do if the local baseline re-solve differs from saved record moments.",
    )
    parser.add_argument(
        "--allow-baseline-mismatch",
        action="store_true",
        help="Deprecated alias for --baseline-mismatch warn.",
    )
    parser.add_argument(
        "--case",
        action="append",
        choices=sorted(available_case_ids()),
        help="Policy case to run. Repeat for multiple cases. Defaults to the presentation policy cases.",
    )
    parser.add_argument("--supply-plus", type=float, default=0.10)
    parser.add_argument("--center-supply-plus", type=float, default=0.10)
    parser.add_argument(
        "--center-supply-elasticity",
        type=float,
        default=2.0,
        help="Counterfactual center housing supply elasticity xi_C for de-zoning.",
    )
    parser.add_argument(
        "--center-amenity-plus",
        type=float,
        default=0.05,
        help="Additive increase in the center location shifter E_C for demand-pressure counterfactuals.",
    )
    parser.add_argument("--universal-phi", type=float, default=0.95)
    parser.add_argument("--parent-phi", type=float, default=0.95)
    parser.add_argument("--retiree-tax-rate", type=float, default=0.10)
    parser.add_argument(
        "--child-cost-subsidy",
        type=float,
        default=0.03,
        help="Reduction in c_bar_n used as a quick proxy for transfers to young parents.",
    )
    args = parser.parse_args()
    if args.allow_baseline_mismatch:
        args.baseline_mismatch = "warn"
    return args


def build_setup(record: dict):
    alpha_cons_bounds = record.get("alpha_cons_bounds")
    if alpha_cons_bounds is not None:
        alpha_cons_bounds = tuple(float(x) for x in alpha_cons_bounds)
    return build_direct_calibration_setup(
        "benchmark",
        geo_weight=100.0,
        population_closure="outside_option_benchmark_normalized",
        scale_target=1.0,
        scale_weight=100.0,
        hR_max=record.get("hR_max"),
        alpha_cons=record.get("alpha_cons"),
        owner_h_bar_scale=record.get("owner_h_bar_scale"),
        owner_size_cost=record.get("owner_size_cost"),
        owner_size_cost_ref=record.get("owner_size_cost_ref"),
        owner_size_cost_power=record.get("owner_size_cost_power"),
        tenure_choice_kappa=record.get("tenure_choice_kappa"),
        weight_overrides=record.get("weight_overrides") or None,
        extra_targets=record.get("extra_targets") or None,
        parent_dp_waiver=record.get("parent_dp_waiver"),
        parent_dp_waiver_phi=record.get("parent_dp_waiver_phi"),
        calibrate_alpha_cons=bool(record.get("calibrate_alpha_cons", False)),
        alpha_cons_bounds=alpha_cons_bounds,
        H_own=record.get("H_own") or None,
    )


def available_case_ids() -> set[str]:
    return {
        "all_H0_plus10",
        "center_H0_plus10",
        "center_H0_relaxed",
        "center_xi_relaxed",
        "center_amenity_plus",
        "center_amenity_plus_eta_relaxed",
        "universal_phi_relief",
        "parent_phi_relief",
        "parent_no_downpayment",
        "property_tax_plus_1pp",
        "retiree_tax_young_parent_transfer",
    }


def build_cases(
    args: argparse.Namespace,
    fixed_modifier: Callable[[SimpleNamespace], None],
) -> list[tuple[str, str, Callable[[SimpleNamespace], None]]]:
    selected = args.case or DEFAULT_CASES
    registry: dict[str, tuple[str, Callable[[SimpleNamespace], None]]] = {
        "all_H0_plus10": (
            f"Increase baseline housing supply H0 in both locations by {100 * args.supply_plus:.1f}%",
            supply_H0_plus(args.supply_plus, center_only=False),
        ),
        "center_H0_plus10": (
            f"Increase baseline housing supply H0 in the center by {100 * args.center_supply_plus:.1f}%",
            supply_H0_plus(args.center_supply_plus, center_only=True),
        ),
        "center_H0_relaxed": (
            f"Center de-zoning: increase baseline center housing supply H0_C by {100 * args.center_supply_plus:.1f}%",
            supply_H0_plus(args.center_supply_plus, center_only=True),
        ),
        "center_xi_relaxed": (
            "De-zoning: hold center supply normalization fixed and "
            f"raise elasticity xi_C to {args.center_supply_elasticity:.3f}",
            center_supply_elasticity(args.center_supply_elasticity),
        ),
        "center_amenity_plus": (
            f"Center demand pressure: raise E_C by {args.center_amenity_plus:.3f}",
            center_amenity_plus(args.center_amenity_plus),
        ),
        "center_amenity_plus_eta_relaxed": (
            f"Center demand pressure plus relaxed eta_C={args.center_supply_elasticity:.3f}, "
            "supply re-anchored at benchmark",
            compose_modifiers(
                center_supply_elasticity_anchored(args.center_supply_elasticity),
                center_amenity_plus(args.center_amenity_plus),
            ),
        ),
        "universal_phi_relief": (
            f"Universal down-payment relief: financed share phi={args.universal_phi:.3f}",
            universal_phi(args.universal_phi),
        ),
        "parent_phi_relief": (
            f"Parent-targeted down-payment relief: financed share phi={args.parent_phi:.3f}",
            parent_phi(args.parent_phi),
        ),
        "parent_no_downpayment": (
            "Remove down-payment constraint for parents",
            parent_no_downpayment,
        ),
        "property_tax_plus_1pp": (
            "Increase property tax by 1 percentage point",
            property_tax_plus_1pp,
        ),
        "retiree_tax_young_parent_transfer": (
            f"{100 * args.retiree_tax_rate:.1f}% pension cut plus "
            f"c_bar_n subsidy proxy of {args.child_cost_subsidy:.4f}",
            retiree_tax_young_parent_transfer(args.retiree_tax_rate, args.child_cost_subsidy),
        ),
    }
    cases: list[tuple[str, str, Callable[[SimpleNamespace], None]]] = [
        ("baseline_fixed_outside", "Baseline with fixed benchmark outside option", fixed_modifier)
    ]
    for case_id in selected:
        description, modifier = registry[case_id]
        cases.append((case_id, description, compose_modifiers(fixed_modifier, modifier)))
    return cases


def solve_record_case(
    case_id: str,
    description: str,
    record: dict,
    setup,
    theta: np.ndarray,
    *,
    modifier: Callable[[SimpleNamespace], None],
    p_init: list[float] | None,
    max_iter_eq: int | None = None,
) -> dict:
    P = record_to_parameters(record, setup, theta)
    if p_init is not None:
        P.p_init_override = np.asarray(p_init, dtype=float)
    modifier(P)
    if max_iter_eq is not None:
        P.max_iter_eq = int(max_iter_eq)
    t0 = time.perf_counter()
    sol, P_out, p_eq = run_model_cp_dt(P, verbose=False)
    elapsed = time.perf_counter() - t0
    moments = extract_moments(sol, P_out, p_eq, 0.0, 0.0, 0.0, True)
    scale = getattr(sol, "accounting_scale", None)
    rents = P_out.user_cost_rate * np.asarray(p_eq, dtype=float)
    return {
        "case": case_id,
        "description": description,
        "elapsed_sec": float(elapsed),
        "moments": namespace_scalars(moments),
        "equilibrium": {
            "p_P": float(p_eq[0]),
            "p_C": float(p_eq[1]),
            "rent_P": float(rents[0]),
            "rent_C": float(rents[1]),
            "rent_ratio_C_over_P": float(rents[1] / rents[0]),
            "pop_share_P": float(sol.pop_share[0]),
            "pop_share_C": float(sol.pop_share[1]),
            "entry_share_P": float(P_out.entry_shares[0]),
            "entry_share_C": float(P_out.entry_shares[1]),
            "entry_total": float(P_out.E_total),
            "tau_H": float(P_out.tau_H),
            "user_cost_rate": float(P_out.user_cost_rate),
            "xi_supply_P": float(np.asarray(P_out.xi_supply, dtype=float).reshape(-1)[0]),
            "xi_supply_C": float(np.asarray(P_out.xi_supply, dtype=float).reshape(-1)[1]),
            "H0_P": float(np.asarray(P_out.H0, dtype=float).reshape(-1)[0]),
            "H0_C": float(np.asarray(P_out.H0, dtype=float).reshape(-1)[1]),
            "E_loc_P": float(np.asarray(P_out.E_loc, dtype=float).reshape(-1)[0]),
            "E_loc_C": float(np.asarray(P_out.E_loc, dtype=float).reshape(-1)[1]),
            "phi_nonparent": float(np.asarray(P_out.phi).reshape(-1)[0]),
            "phi_parent": float(np.asarray(P_out.phi).reshape(-1)[1]) if len(np.asarray(P_out.phi).reshape(-1)) > 1 else math.nan,
            "parent_dp_waiver": int(bool(getattr(P_out, "parent_dp_waiver", False))),
            "parent_dp_waiver_phi": float(getattr(P_out, "parent_dp_waiver_phi", math.nan)),
            "pension": float(getattr(P_out, "pension", math.nan)),
            "c_bar_n": float(getattr(P_out, "c_bar_n", math.nan)),
            "population_closure": str(getattr(P_out, "population_closure", "")),
        },
        "scale": namespace_scalars(scale) if scale is not None else {},
        "timings": jsonable_mapping(getattr(sol, "timings", {})),
    }


def record_to_parameters(record: dict, setup, theta: np.ndarray) -> SimpleNamespace:
    theta_dict = {name: float(value) for name, value in zip(setup.names, theta)}
    P = SimpleNamespace(**copy.deepcopy(asdict(setup.P_base)))
    structural_names = [
        name
        for name in setup.names
        if name not in DIRECT_GEOMETRY_NAMES and name not in (OUTSIDE_VALUE_NAME, RENEWAL_FLOW_NAME)
    ]
    P = apply_theta(P, [theta_dict[name] for name in structural_names], structural_names)
    P.E_loc = np.array([float(P.E_loc[0]), theta_dict["E_C"]])
    P.r_bar = np.array([float(P.r_bar[0]), theta_dict["r_bar_C"]])
    if record.get("phi") is not None:
        P.phi = np.asarray(record["phi"], dtype=float).reshape(-1)
    if OUTSIDE_VALUE_NAME in theta_dict:
        P.outside_value = theta_dict[OUTSIDE_VALUE_NAME]
        P.outside_value_is_calibrated = True
    if RENEWAL_FLOW_NAME in theta_dict:
        P.outside_entry_flow = theta_dict[RENEWAL_FLOW_NAME]
    return P


def no_change(_P: SimpleNamespace) -> None:
    return None


def fixed_outside_modifier(outside_value: float, outside_entry_flow: float, kappa_entry: float):
    def _modifier(P: SimpleNamespace) -> None:
        P.population_closure = "accounting_scale_prices"
        P.outside_value = float(outside_value)
        P.outside_value_is_calibrated = True
        P.allow_uncalibrated_outside_value = False
        P.calibrate_outside_value_to_entry_prob = False
        P.outside_entry_flow = float(outside_entry_flow)
        P.kappa_entry = float(kappa_entry)

    return _modifier


def parent_no_downpayment(P: SimpleNamespace) -> None:
    P.parent_dp_waiver = True
    P.parent_dp_waiver_phi = 1.0
    P.parent_dp_waiver_locations = np.array([], dtype=int)
    P.parent_dp_waiver_owner_rungs = np.array([], dtype=int)
    P.parent_dp_waiver_birth_state_only = False


def universal_phi(phi_value: float):
    def _modifier(P: SimpleNamespace) -> None:
        P.phi = float(phi_value) * np.ones_like(np.asarray(P.phi, dtype=float).reshape(-1))

    return _modifier


def parent_phi(phi_value: float):
    def _modifier(P: SimpleNamespace) -> None:
        phi = np.asarray(P.phi, dtype=float).reshape(-1).copy()
        if phi.size <= 1:
            phi = np.repeat(phi, int(getattr(P, "n_parity", 1)))
        if phi.size > 1:
            phi[1:] = np.maximum(phi[1:], float(phi_value))
        P.phi = phi

    return _modifier


def supply_H0_plus(percent: float, *, center_only: bool):
    def _modifier(P: SimpleNamespace) -> None:
        H0 = np.asarray(P.H0, dtype=float).reshape(-1).copy()
        if center_only:
            H0[1] *= 1.0 + float(percent)
        else:
            H0 *= 1.0 + float(percent)
        P.H0 = H0

    return _modifier


def center_supply_elasticity(value: float):
    def _modifier(P: SimpleNamespace) -> None:
        base = getattr(P, "eta_supply", None)
        if base is None:
            base = getattr(P, "xi_supply")
        eta = np.asarray(base, dtype=float).reshape(-1).copy()
        if eta.size < 2:
            eta = np.resize(eta, 2)
        eta[1] = float(value)
        P.eta_supply = eta
        P.xi_supply = eta.copy()

    return _modifier


def center_supply_elasticity_anchored(value: float):
    def _modifier(P: SimpleNamespace) -> None:
        base = getattr(P, "eta_supply", None)
        if base is None:
            base = getattr(P, "xi_supply")
        eta = np.asarray(base, dtype=float).reshape(-1).copy()
        if eta.size < 2:
            eta = np.resize(eta, 2)
        H0 = np.asarray(P.H0, dtype=float).reshape(-1).copy()
        r_bar = np.asarray(P.r_bar, dtype=float).reshape(-1)
        p_ref = np.asarray(getattr(P, "p_init_override", np.full(2, np.nan)), dtype=float).reshape(-1)
        if H0.size >= 2 and r_bar.size >= 2 and p_ref.size >= 2:
            user_cost_rate = getattr(P, "user_cost_rate", None)
            if user_cost_rate is None or not math.isfinite(float(user_cost_rate)):
                user_cost_rate = float(getattr(P, "q", 0.0)) + float(getattr(P, "delta", 0.0)) + float(getattr(P, "tau_H", 0.0))
            rent_ref_C = float(user_cost_rate) * float(p_ref[1])
            ratio = rent_ref_C / max(float(r_bar[1]), 1e-12)
            if math.isfinite(ratio) and ratio > 0.0:
                H0[1] *= ratio ** (float(eta[1]) - float(value))
                P.H0 = H0
        eta[1] = float(value)
        P.eta_supply = eta
        P.xi_supply = eta.copy()

    return _modifier


def center_amenity_plus(delta: float):
    def _modifier(P: SimpleNamespace) -> None:
        E = np.asarray(P.E_loc, dtype=float).reshape(-1).copy()
        if E.size < 2:
            E = np.resize(E, 2)
        E[1] += float(delta)
        P.E_loc = E

    return _modifier


def property_tax_plus_1pp(P: SimpleNamespace) -> None:
    P.tau_H = float(P.tau_H) + 0.01
    P.user_cost_rate = float(P.q) + float(P.delta) + float(P.tau_H)


def retiree_tax_young_parent_transfer(tax_rate: float, child_cost_subsidy: float):
    def _modifier(P: SimpleNamespace) -> None:
        base_pension = float(getattr(P, "pension", math.nan))
        if not math.isfinite(base_pension):
            base_pension = float(apply_overrides(setup_parameters(), P).pension)
        P.pension_mode = "manual"
        P.pension = (1.0 - float(tax_rate)) * base_pension
        P.c_bar_n = max(0.0, float(P.c_bar_n) - float(child_cost_subsidy))

    return _modifier


def compose_modifiers(*modifiers):
    def _modifier(P: SimpleNamespace) -> None:
        for modifier in modifiers:
            modifier(P)

    return _modifier


def baseline_check_rows(saved: dict, resolved: dict) -> list[dict]:
    rows = []
    for key in CHECK_KEYS:
        if key not in saved or key not in resolved:
            continue
        saved_val = float(saved[key])
        resolved_val = float(resolved[key])
        rows.append(
            {
                "moment": key,
                "saved": saved_val,
                "resolved": resolved_val,
                "gap": resolved_val - saved_val,
                "abs_gap": abs(resolved_val - saved_val),
            }
        )
    return rows


def summary_rows(results: list[dict]) -> list[dict]:
    baseline = results[0]["moments"] if results else {}
    rows = []
    for result in results:
        row = {"case": result["case"], "description": result["description"]}
        for key in SUMMARY_KEYS:
            val = result["moments"].get(key, math.nan)
            base_val = baseline.get(key, math.nan)
            row[key] = val
            row[f"delta_{key}"] = val - base_val if math.isfinite(val) and math.isfinite(base_val) else math.nan
        tfr = float(row.get("tfr", math.nan))
        delta_tfr = float(row.get("delta_tfr", math.nan))
        base_tfr = float(baseline.get("tfr", math.nan))
        row["births_per_100_women"] = 100.0 * tfr if math.isfinite(tfr) else math.nan
        row["delta_births_per_100_women"] = 100.0 * delta_tfr if math.isfinite(delta_tfr) else math.nan
        row["tfr_pct_change"] = delta_tfr / base_tfr if math.isfinite(delta_tfr) and math.isfinite(base_tfr) and abs(base_tfr) > 1e-12 else math.nan
        row.update(
            {
                "pop_share_C": result["equilibrium"].get("pop_share_C", math.nan),
                "rent_ratio_C_over_P": result["equilibrium"].get("rent_ratio_C_over_P", math.nan),
                "scale_factor": result["scale"].get("scale_factor", math.nan),
                "implied_total_population": result["scale"].get("implied_total_population", math.nan),
                "city_entry_prob_total": result["scale"].get("city_entry_prob_total", math.nan),
                "outside_entry_prob": result["scale"].get("outside_entry_prob", math.nan),
                "elapsed_sec": result.get("elapsed_sec", math.nan),
                "best_eq_error": result["timings"].get("best_eq_error", math.nan),
                "convergence_reason": result["timings"].get("convergence_reason", ""),
            }
        )
        rows.append(row)
    return rows


def equilibrium_rows(results: list[dict]) -> list[dict]:
    rows = []
    scale_keys = [
        "scale_factor",
        "implied_total_population",
        "implied_entry_total",
        "implied_potential_entrant_mass",
        "implied_outside_entry_mass",
        "outside_entry_flow",
        "city_entry_prob_total",
        "outside_entry_prob",
        "denominator",
        "stationary_entry_relative_residual",
        "reference_housing_demand_P",
        "reference_housing_demand_C",
        "implied_housing_demand_P",
        "implied_housing_demand_C",
    ]
    for result in results:
        row = {"case": result["case"], "description": result["description"]}
        row.update(result["equilibrium"])
        for key in scale_keys:
            if key in result["scale"]:
                row[key] = result["scale"][key]
        rows.append(row)
    return rows


def namespace_scalars(obj) -> dict[str, float]:
    if obj is None:
        return {}
    out = {}
    for key, value in vars(obj).items():
        arr = np.asarray(value)
        if arr.size == 1:
            try:
                val = float(arr.reshape(-1)[0])
            except (TypeError, ValueError):
                continue
            if math.isfinite(val):
                out[key] = val
    return out


def jsonable_mapping(mapping: dict) -> dict:
    out = {}
    for key, value in mapping.items():
        if isinstance(value, np.ndarray):
            out[key] = value.tolist()
        elif isinstance(value, np.generic):
            out[key] = value.item()
        else:
            out[key] = value
    return out


def policy_metadata(args: argparse.Namespace) -> dict:
    return {
        "record": str(args.record),
        "selected_cases": args.case or DEFAULT_CASES,
        "baseline_tolerance": float(args.baseline_tol),
        "baseline_mismatch": str(args.baseline_mismatch),
        "supply_plus": float(args.supply_plus),
        "center_supply_plus": float(args.center_supply_plus),
        "center_supply_elasticity": float(args.center_supply_elasticity),
        "center_amenity_plus": float(args.center_amenity_plus),
        "universal_phi": float(args.universal_phi),
        "parent_phi": float(args.parent_phi),
        "retiree_tax_rate": float(args.retiree_tax_rate),
        "child_cost_subsidy": float(args.child_cost_subsidy),
        "retiree_transfer_policy_object": (
            "presentation diagnostic: reduce manual pension by retiree_tax_rate "
            "and lower c_bar_n by child_cost_subsidy as a child-resource subsidy proxy"
        ),
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, sort_keys=True))


def write_markdown_summary(path: Path, rows: list[dict], record_path: Path, scale: dict) -> None:
    fields = [
        ("tfr", "TFR"),
        ("childless_rate", "Childless"),
        ("own_rate", "Own 30-55"),
        ("center_share_newparents", "Center new parents"),
        ("housing_increment_0to1", "H01"),
        ("housing_increment_1to2", "H12"),
        ("scale_factor", "Scale"),
        ("rent_ratio_C_over_P", "Rent ratio"),
    ]
    lines = [
        "# Policy Counterfactuals",
        "",
        f"- Record: `{record_path}`",
        "- Closure: baseline outside option objects fixed; counterfactuals solved with `accounting_scale_prices`.",
        f"- Fixed outside value: `{float(scale.get('outside_value', math.nan)):.6g}`",
        f"- Fixed kappa entry: `{float(scale.get('kappa_entry', math.nan)):.6g}`",
        f"- Fixed outside entry flow: `{float(scale.get('outside_entry_flow', math.nan)):.6g}`",
        "",
        "| Case | " + " | ".join(label for _key, label in fields) + " |",
        "|---" + "|---:" * len(fields) + "|",
    ]
    for row in rows:
        vals = []
        for key, _label in fields:
            val = row.get(key, math.nan)
            vals.append(format_cell(val))
        lines.append(f"| {row.get('case', '')} | " + " | ".join(vals) + " |")
    path.write_text("\n".join(lines) + "\n")


def write_validation_note(path: Path, results: list[dict], record_path: Path, fixed_scale: dict, metadata: dict) -> None:
    fixed_keys = ["outside_value", "kappa_entry", "outside_entry_flow"]
    lines = [
        "# Counterfactual Validation Note",
        "",
        f"- Record: `{record_path}`",
        "- Closure: benchmark outside-option objects are recovered once, then all cases are solved with `accounting_scale_prices`.",
        f"- De-zoning implementation: set center housing supply elasticity `xi_C` to `{float(metadata.get('center_supply_elasticity', math.nan)):.6g}`.",
        f"- Parent down-payment implementation: set parent-state financed share `phi` to `{float(metadata.get('parent_phi', math.nan)):.6g}` when selected.",
        "",
        "| Case | Outside value | Gap | Kappa entry | Gap | Outside entry flow | Gap |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    fixed = {key: float(fixed_scale.get(key, math.nan)) for key in fixed_keys}
    for result in results:
        scale = result.get("scale", {})
        vals = []
        for key in fixed_keys:
            val = float(scale.get(key, math.nan))
            gap = val - fixed[key] if math.isfinite(val) and math.isfinite(fixed[key]) else math.nan
            vals.extend([format_cell_precise(val), format_cell_precise(gap)])
        lines.append(f"| {result.get('case', '')} | " + " | ".join(vals) + " |")
    path.write_text("\n".join(lines) + "\n")


def write_latex_table(path: Path, rows: list[dict]) -> None:
    display = [
        ("tfr", "TFR", ".3f"),
        ("delta_tfr", "$\\Delta$ TFR", ".3f"),
        ("delta_births_per_100_women", "$\\Delta$ births/100", ".2f"),
        ("own_rate", "Own", ".3f"),
        ("center_share_newparents", "Center new parents", ".3f"),
        ("housing_increment_0to1", "$H_{01}$", ".3f"),
        ("housing_increment_1to2", "$H_{12}$", ".3f"),
        ("scale_factor", "Scale", ".3f"),
    ]
    lines = [
        "\\begin{tabular}{l" + "r" * len(display) + "}",
        "\\toprule",
        "Scenario & " + " & ".join(label for _key, label, _fmt in display) + " \\\\",
        "\\midrule",
    ]
    for row in rows:
        label = latex_escape(row.get("case", ""))
        vals = []
        for key, _label, fmt in display:
            val = row.get(key, math.nan)
            vals.append(format_latex_number(val, fmt))
        lines.append(label + " & " + " & ".join(vals) + " \\\\")
    lines.extend(["\\bottomrule", "\\end{tabular}", ""])
    path.write_text("\n".join(lines))


def write_summary_plot(path: Path, rows: list[dict]) -> None:
    policy_rows = [row for row in rows if row.get("case") != "baseline_fixed_outside"]
    if not policy_rows:
        return
    labels = [pretty_case_label(str(row.get("case", ""))) for row in policy_rows]
    metrics = [
        ("delta_tfr", "$\\Delta$ TFR"),
        ("delta_own_rate", "$\\Delta$ ownership"),
        ("scale_factor", "City scale"),
        ("rent_ratio_C_over_P", "Rent ratio C/P"),
    ]
    fig, axes = plt.subplots(1, len(metrics), figsize=(13.5, 3.6), constrained_layout=True)
    colors = ["#1f5fa6", "#c73e3a", "#3f6f43"]
    for ax, (key, title) in zip(axes, metrics):
        vals = [float(row.get(key, math.nan)) for row in policy_rows]
        if key in ("scale_factor", "rent_ratio_C_over_P"):
            base = float(rows[0].get(key, math.nan))
            vals = [val - base if math.isfinite(val) and math.isfinite(base) else math.nan for val in vals]
            title = "$\\Delta$ " + title
        ax.axhline(0.0, color="#444444", linewidth=0.8)
        ax.bar(np.arange(len(vals)), vals, color=colors[: len(vals)], width=0.68)
        ax.set_title(title)
        ax.set_xticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=25, ha="right")
        ax.grid(True, axis="y", alpha=0.3)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def pretty_case_label(case_id: str) -> str:
    labels = {
        "baseline_fixed_outside": "Baseline",
        "parent_no_downpayment": "No parent down payment",
        "property_tax_plus_1pp": "Property tax +1pp",
        "retiree_tax_young_parent_transfer": "Old-to-young transfer",
        "all_H0_plus10": "All supply +10%",
        "center_H0_plus10": "Center supply +10%",
        "center_H0_relaxed": "Center de-zoning",
        "center_xi_relaxed": "Center de-zoning",
        "center_amenity_plus": "Center demand",
        "center_amenity_plus_eta_relaxed": "Demand + relaxed eta",
        "universal_phi_relief": "Universal DP relief",
        "parent_phi_relief": "Parent DP relief",
    }
    return labels.get(case_id, case_id.replace("_", " "))


def format_cell(value) -> str:
    try:
        val = float(value)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(val):
        return ""
    return f"{val:.3f}"


def format_cell_precise(value) -> str:
    try:
        val = float(value)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(val):
        return ""
    return f"{val:.8g}"


def format_latex_number(value, fmt: str) -> str:
    try:
        val = float(value)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(val):
        return ""
    return format(val, fmt)


def latex_escape(value: str) -> str:
    return str(value).replace("_", "\\_")


if __name__ == "__main__":
    main()
