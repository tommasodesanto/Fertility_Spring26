#!/usr/bin/env python3
"""Write the E2 winner's standard diagnostics, timing supplement, and policy packet.

This is a deliberately single-process inspection script: it solves the winner
once, reuses that solution as the policy baseline, and then solves each policy
arm serially in the same Python process.  It does not alter calibration code.
"""

from __future__ import annotations

import argparse
import csv
import importlib
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from intergen_eqscale_seq_optimized.calibration import (  # noqa: E402
    diagnostic_loss,
    extract_moments,
    get_target_set,
    jsonable,
)
from intergen_eqscale_seq_optimized.diagnostics import write_diagnostics  # noqa: E402
from intergen_eqscale_seq_optimized.e1_profile import E1_TARGET_SET, e1_overrides  # noqa: E402
from intergen_eqscale_seq_optimized.externals import flhsv_income_overrides  # noqa: E402
from intergen_eqscale_seq_optimized.parameters import get_fecundity_by_age  # noqa: E402
from intergen_eqscale_seq_optimized.solver import run_model_cp_dt  # noqa: E402

SOURCE = ROOT / "output/model/eqscale_seq_optimized_recalibration_20260719/report/results.json"
DEFAULT_OUTDIR = ROOT / "output/model/eqscale_seq_policy_packet_20260720"
E4_SOURCE = ROOT / "output/model/eqscale_seq_e4_split_recalibration_20260723/report/results.json"
E4_DEFAULT_OUTDIR = ROOT / "output/model/eqscale_seq_e4_policy_packet_20260724"
E5_SOURCE = ROOT / "output/model/eqscale_seq_e5b_recalibration_20260725/report/results.json"
E5_DEFAULT_OUTDIR = ROOT / "output/model/eqscale_seq_e5b_policy_packet_20260725"
ROOM_TARGET = "aggregate_mean_occupied_rooms_18_85"
ROOM_TARGET_VALUE = 5.779970481941968
ROOM_WEIGHT = 6.0
E4_EXTERNALS = {
    "n_parity": 4,
    "fertility_units": "literal_topcode",
    "tfr_top_bin_weight": 3.4,
    "entrant_conversion_factor": 0.5,
    "child_bin_high_cutoff": 3,
    "eqscale_form": "power",
    "gamma_e": 0.0,
}

POLICY_CASES: tuple[dict[str, Any], ...] = (
    {"case": "baseline", "label": "Baseline", "overrides": {}},
    {
        "case": "birth_grant_A0p4_Hge6",
        "label": "Birth grant, family-size homes",
        "overrides": {
            "birth_entry_grant": True,
            "birth_entry_grant_amount": 0.4,
            "birth_entry_grant_owner_rungs": np.array([3, 4, 5]),
            "birth_entry_grant_locations": np.array([], dtype=int),
        },
    },
    {
        "case": "tax2_grant_A0p4_Hge6",
        "label": "Property tax plus birth grant",
        "overrides": {
            "birth_entry_grant": True,
            "birth_entry_grant_amount": 0.4,
            "birth_entry_grant_owner_rungs": np.array([3, 4, 5]),
            "birth_entry_grant_locations": np.array([], dtype=int),
            "tau_H": 0.08,
        },
    },
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--arm", choices=("e2", "e4", "e5"), default="e2")
    args = parser.parse_args()
    if args.arm == "e4":
        if args.source == SOURCE:
            args.source = E4_SOURCE
        if args.outdir == DEFAULT_OUTDIR:
            args.outdir = E4_DEFAULT_OUTDIR
    if args.arm == "e5":
        if args.source == SOURCE:
            args.source = E5_SOURCE
        if args.outdir == DEFAULT_OUTDIR:
            args.outdir = E5_DEFAULT_OUTDIR
    return args


def load_winner_theta(source: Path, arm: str) -> dict[str, float]:
    payload = json.loads(source.read_text())
    theta = ((payload.get("winners") or {}).get("E1") or {}).get("theta")
    if not isinstance(theta, dict):
        raise ValueError(f"{source} does not contain winners.E1.theta")
    required = {"H0", "alpha_cons", "beta", "chi", "delta_alpha", "delta_alpha_jump",
                "gamma_e", "kappa_fert", "psi_child", "tenure_choice_kappa", "theta0",
                "theta1", "theta_n"}
    if arm == "e4":
        required.remove("gamma_e")
        required.add("kappa_fert_continuation")
    if set(theta) != required:
        raise ValueError(f"unexpected {arm.upper()} winner keys: {sorted(theta)}")
    return {name: float(value) for name, value in theta.items()}


def arm_externals(arm: str) -> dict[str, Any]:
    if arm == "e2":
        return {}
    if arm == "e4":
        return {**E4_EXTERNALS, **flhsv_income_overrides()}
    raise ValueError(f"unknown packet arm: {arm}")


def target_system() -> tuple[dict[str, float], dict[str, float]]:
    targets, weights = get_target_set(E1_TARGET_SET)
    targets[ROOM_TARGET] = ROOM_TARGET_VALUE
    weights[ROOM_TARGET] = ROOM_WEIGHT
    if len(targets) != 15 or set(targets) != set(weights):
        raise RuntimeError("E2 packet requires the unchanged 15-target E1 system")
    return targets, weights


def age_grid(P: Any, length: int) -> np.ndarray:
    return float(P.age_start) + float(P.da) * np.arange(length, dtype=float)


def as_array(sol: Any, name: str) -> np.ndarray:
    return np.asarray(getattr(sol, name, np.array([], dtype=float)), dtype=float).reshape(-1)


def solve_case(theta: dict[str, float], policy_overrides: dict[str, Any], targets: dict[str, float], weights: dict[str, float], arm: str) -> dict[str, Any]:
    # Policy settings are intentionally last, matching the requested legacy arms.
    overrides = {**e1_overrides(tight=True, optimized=True), **arm_externals(arm), **theta, **policy_overrides}
    started = time.perf_counter()
    sol, P, p_eq = run_model_cp_dt(overrides, verbose=False)
    moments = extract_moments(sol, P)
    return {
        "sol": sol,
        "P": P,
        "p_eq": np.asarray(p_eq, dtype=float).reshape(-1),
        "moments": moments,
        "loss": float(diagnostic_loss(moments, targets=targets, weights=weights)),
        "market_residual": float(getattr(sol, "best_max_abs_rel_excess", math.nan)),
        "elapsed_sec": time.perf_counter() - started,
        "overrides": overrides,
    }


def e5_chain_overrides(source: Path) -> tuple[dict[str, Any], dict[str, float], dict[str, float], dict[str, float]]:
    """Build the E5 runtime exactly through its calibration-chain entry point."""
    os.environ["E3_L4"] = "1"
    os.environ["E5"] = "1"
    os.environ["E3_TFR_TOP_BIN_WEIGHT"] = "3.602359422009"
    from intergen_eqscale_seq_optimized import run_e1_chain

    chain = importlib.reload(run_e1_chain)
    chain.load_runtime()
    chain_args = argparse.Namespace(J=17, Nb=120, max_iter_eq=40, tol_eq=2.5e-5)
    overrides = chain.common_overrides(chain_args)
    payload = json.loads(source.read_text())
    theta = ((payload.get("winners") or {}).get("E1") or {}).get("theta")
    if not isinstance(theta, dict):
        raise ValueError(f"{source} does not contain winners.E1.theta")
    theta = {name: float(value) for name, value in theta.items()}
    overrides.update(theta)
    # The E5 contract is not in the calibration TARGET_SETS registry; the
    # chain itself routes E5 through the profile's target system.
    from intergen_eqscale_seq_optimized.e5_profile import e5_target_system

    system = e5_target_system()
    targets, weights = system.targets_dict(), system.weights_dict()
    return overrides, theta, targets, weights


def solve_e5_case(overrides: dict[str, Any], policy_overrides: dict[str, Any], targets: dict[str, float], weights: dict[str, float]) -> dict[str, Any]:
    """Solve an E5 policy case from the chain-built baseline overrides."""
    case_overrides = {**overrides, **policy_overrides}
    started = time.perf_counter()
    sol, P, p_eq = run_model_cp_dt(case_overrides, verbose=False)
    moments = extract_moments(sol, P)
    return {
        "sol": sol,
        "P": P,
        "p_eq": np.asarray(p_eq, dtype=float).reshape(-1),
        "moments": moments,
        "loss": float(diagnostic_loss(moments, targets=targets, weights=weights)),
        "market_residual": float(getattr(sol, "best_max_abs_rel_excess", math.nan)),
        "elapsed_sec": time.perf_counter() - started,
        "overrides": case_overrides,
    }


def verify_e5_winner(moments: dict[str, Any], source: Path) -> None:
    """Require exact reproduction of every certified E5 target-fit row."""
    target_fit_path = source.parent / "target_fit_full.csv"
    with target_fit_path.open(newline="", encoding="utf-8") as handle:
        certified_rows = [(row["moment"], float(row["model"])) for row in csv.DictReader(handle)]
    if len(certified_rows) != 12 or len({name for name, _ in certified_rows}) != 12:
        raise ValueError(f"Expected 12 distinct certified E5 rows in {target_fit_path}, found {len(certified_rows)}")
    rows = []
    for name, certified_value in certified_rows:
        if name not in moments:
            raise ValueError(f"E5 solve does not report certified moment {name!r}")
        model_value = float(moments[name])
        rows.append((name, model_value, certified_value, abs(model_value - certified_value)))
    mismatches = [row for row in rows if row[3] > 1e-6]
    if mismatches:
        table = "\n".join(
            ["moment,solved,certified,abs_diff"]
            + [f"{name},{model:.15g},{expected:.15g},{diff:.3e}" for name, model, expected, diff in mismatches]
        )
        raise RuntimeError(f"Certified E5 moment verification failed; no packet output written.\n{table}")


def write_policy_births_summary(path: Path, records: list[dict[str, Any]]) -> None:
    """Write the compact fertility and ownership policy readout for E5."""
    base = records[0]
    columns = (
        ("tfr", "completed_fertility", "tfr"),
        ("childless_rate", "completed_fertility", "childless_rate"),
        ("mean_age_first_birth", "present_fertility", "mean_age_first_birth"),
        ("share_first_births_age30plus", "present_fertility", "share_first_births_age30plus"),
        ("own_rate", "equilibrium", "own_rate"),
    )
    fieldnames = ["case", "label"] + [name for name, _, _ in columns] + [f"{name}_change_vs_baseline" for name, _, _ in columns]
    rows = []
    for record in records:
        row: dict[str, Any] = {"case": record["case"], "label": record["label"]}
        for name, group, key in columns:
            value = float(record[group][key])
            baseline = float(base[group][key])
            row[name] = value
            row[f"{name}_change_vs_baseline"] = value - baseline
        rows.append(row)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def target_fit(moments: dict[str, Any], targets: dict[str, float], weights: dict[str, float]) -> list[dict[str, float | str]]:
    rows = []
    for name, target in targets.items():
        model = float(moments.get(name, math.nan))
        gap = model - float(target)
        rows.append({"moment": name, "target": float(target), "model": model, "gap": gap,
                     "weight": float(weights[name]), "loss_contribution": float(weights[name]) * gap ** 2})
    return rows


def fertility_record(result: dict[str, Any]) -> dict[str, Any]:
    sol, P, moments = result["sol"], result["P"], result["moments"]
    shares = as_array(sol, "parity_dist")
    return {
        "present_fertility": {
            "total_births_kfe": float(getattr(sol, "total_births_kfe", math.nan)),
            "attempt_hazard_by_age": as_array(sol, "attempt_hazard_by_age"),
            "first_birth_hazard_by_age": as_array(sol, "first_birth_hazard_by_age"),
            "second_attempt_hazard_by_age": as_array(sol, "second_attempt_hazard_by_age"),
            "second_birth_hazard_by_age": as_array(sol, "second_birth_hazard_by_age"),
            "mean_age_first_birth": float(getattr(sol, "mean_age_first_birth", math.nan)),
            "share_first_births_age30plus": float(getattr(sol, "share_first_births_age30plus", math.nan)),
        },
        "completed_fertility": {
            "tfr": float(moments.get("tfr", math.nan)),
            "childless_rate": float(moments.get("childless_rate", math.nan)),
            "childless_chosen_45": float(getattr(sol, "childless_chosen_45", math.nan)),
            "childless_clock_45": float(getattr(sol, "childless_clock_45", math.nan)),
            "family_size_shares": shares,
            "progression_1to2": float(getattr(sol, "parity_progression_1to2_flow", math.nan)),
        },
        "equilibrium": {
            "p_eq": result["p_eq"],
            "own_rate": float(moments.get("own_rate", math.nan)),
            "old_age_own_rate": float(moments.get("old_age_own_rate", math.nan)),
            "market_residual": float(result["market_residual"]),
            "loss_15_target": float(result["loss"]),
        },
        "fecundity_by_age": get_fecundity_by_age(P),
    }


def plot_line(path: Path, ages: np.ndarray, series: list[tuple[str, np.ndarray, str]], title: str, ylabel: str = "Hazard") -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    for label, values, style in series:
        n = min(len(ages), len(values))
        ax.plot(ages[:n], values[:n], style, linewidth=2, label=label)
    ax.set(title=title, xlabel="Age", ylabel=ylabel)
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_timing_supplement(outdir: Path, result: dict[str, Any]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    sol, P = result["sol"], result["P"]
    ages = age_grid(P, int(P.J))
    attempts, first = as_array(sol, "attempt_hazard_by_age"), as_array(sol, "first_birth_hazard_by_age")
    second_attempts, second = as_array(sol, "second_attempt_hazard_by_age"), as_array(sol, "second_birth_hazard_by_age")
    plot_line(outdir / "first_birth_attempt_realized_fecundity.png", ages, [("Attempt", attempts, "-"), ("Realized", first, "-"), ("Fecundity $\\pi_j$", get_fecundity_by_age(P), "--")], "First-birth attempt and realized hazards")
    plot_line(outdir / "second_birth_attempt_realized.png", ages, [("Attempt", second_attempts, "-"), ("Realized", second, "-")], "Second-birth attempt and realized hazards")
    dist = as_array(sol, "first_birth_age_distribution")
    n = min(len(ages), len(dist))
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.bar(ages[:n], dist[:n], width=max(float(P.da) * 0.8, 0.5), color="#4c78a8")
    mean = float(getattr(sol, "mean_age_first_birth", math.nan))
    if np.isfinite(mean):
        ax.axvline(mean, color="#e45756", linewidth=2, label=f"Mean = {mean:.1f}")
        ax.legend(frameon=False)
    ax.set(title="First-birth age distribution", xlabel="Age", ylabel="Share of first births")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout(); fig.savefig(outdir / "first_birth_age_distribution.png", dpi=180); plt.close(fig)
    chosen, clock = float(getattr(sol, "childless_chosen_45", math.nan)), float(getattr(sol, "childless_clock_45", math.nan))
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    ax.bar(["Childlessness at 45"], [chosen], label="Chosen", color="#4c78a8")
    ax.bar(["Childlessness at 45"], [clock], bottom=[chosen], label="Clock", color="#f58518")
    ax.set(ylabel="Share", title="Chosen versus clock childlessness")
    ax.legend(frameon=False); ax.grid(axis="y", alpha=0.25)
    fig.tight_layout(); fig.savefig(outdir / "childlessness_chosen_clock.png", dpi=180); plt.close(fig)
    shares = as_array(sol, "parity_dist")
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    labels = ["0 children", "1 child", "2 children"]
    ax.bar(labels[:len(shares)], shares[:len(labels)], color=["#4c78a8", "#72b7b2", "#54a24b"][:len(shares)])
    ax.set(ylabel="Share", title="Completed family-size shares")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout(); fig.savefig(outdir / "completed_family_size_shares.png", dpi=180); plt.close(fig)


def write_hazard_delta(path: Path, baseline: dict[str, Any], case: dict[str, Any]) -> None:
    base_hazard = as_array(baseline["sol"], "first_birth_hazard_by_age")
    case_hazard = as_array(case["sol"], "first_birth_hazard_by_age")
    n = min(len(base_hazard), len(case_hazard), int(case["P"].J))
    ages = age_grid(case["P"], n)
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.plot(ages, case_hazard[:n] - base_hazard[:n], linewidth=2, color="#e45756")
    ax.set(title="First-birth realized-hazard change relative to baseline", xlabel="Age", ylabel="Hazard difference")
    ax.grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig)


def flatten_record(record: dict[str, Any]) -> dict[str, Any]:
    present, completed, equilibrium = record["present_fertility"], record["completed_fertility"], record["equilibrium"]
    return {"case": record["case"], "label": record["label"],
            **{f"present_{k}": json.dumps(jsonable(v)) if isinstance(v, np.ndarray) else v for k, v in present.items()},
            **{f"completed_{k}": json.dumps(jsonable(v)) if isinstance(v, np.ndarray) else v for k, v in completed.items()},
            **equilibrium}


def pct_delta(value: float, baseline: float) -> str:
    if not (np.isfinite(value) and np.isfinite(baseline)) or baseline == 0:
        return "NA"
    return f"{100.0 * (value / baseline - 1.0):+.2f}%"


def write_policy_table(path: Path, records: list[dict[str, Any]]) -> None:
    base = records[0]
    rows = ["| Case | Present: total births | Present: mean first-birth age | Present: first births age 30+ | Completed: TFR | Completed: childlessness | Completed: chosen childlessness | Completed: clock childlessness | Completed: family size 0 | Completed: family size 1 | Completed: family size 2 | Completed: progression 1to2 | Own rate | Old-age own rate | Loss |", "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    keys = [("present_fertility", "total_births_kfe"), ("present_fertility", "mean_age_first_birth"), ("present_fertility", "share_first_births_age30plus"), ("completed_fertility", "tfr"), ("completed_fertility", "childless_rate"), ("completed_fertility", "childless_chosen_45"), ("completed_fertility", "childless_clock_45"), ("completed_fertility", "progression_1to2"), ("equilibrium", "own_rate"), ("equilibrium", "old_age_own_rate"), ("equilibrium", "loss_15_target")]
    for record in records:
        values = [pct_delta(float(record[group][key]), float(base[group][key])) for group, key in keys[:7]]
        for index in range(3):
            value = np.asarray(record["completed_fertility"]["family_size_shares"], dtype=float)
            baseline_value = np.asarray(base["completed_fertility"]["family_size_shares"], dtype=float)
            values.append(pct_delta(float(value[index]) if index < len(value) else math.nan,
                                    float(baseline_value[index]) if index < len(baseline_value) else math.nan))
        values.extend(pct_delta(float(record[group][key]), float(base[group][key])) for group, key in keys[7:])
        rows.append("| " + str(record["case"]) + " | " + " | ".join(values) + " |")
    path.write_text("# Policy cases: percent change relative to baseline\n\nPresent and completed fertility are reported separately. Values are percent changes; baseline is therefore `+0.00%`.\n\n" + "\n".join(rows) + "\n")


def main() -> None:
    args = parse_args()
    outdir = args.outdir
    if args.arm != "e5":
        outdir.mkdir(parents=True, exist_ok=True)
    results: dict[str, dict[str, Any]] = {}
    records: list[dict[str, Any]] = []

    if args.arm == "e5":
        e5_overrides, theta, targets, weights = e5_chain_overrides(args.source)
        baseline_policy = POLICY_CASES[0]
        baseline = solve_e5_case(e5_overrides, dict(baseline_policy["overrides"]), targets, weights)
        verify_e5_winner(baseline["moments"], args.source)
        results["baseline"] = baseline
        baseline_record = {"case": baseline_policy["case"], "label": baseline_policy["label"],
                           "policy_overrides": baseline_policy["overrides"], **fertility_record(baseline),
                           "target_fit": target_fit(baseline["moments"], targets, weights),
                           "elapsed_sec": baseline["elapsed_sec"]}
        records.append(baseline_record)
        remaining_policies = POLICY_CASES[1:]
    else:
        theta, (targets, weights) = load_winner_theta(args.source, args.arm), target_system()
        remaining_policies = POLICY_CASES

    for policy in remaining_policies:
        case = str(policy["case"])
        result = (solve_e5_case(e5_overrides, dict(policy["overrides"]), targets, weights)
                  if args.arm == "e5"
                  else solve_case(theta, dict(policy["overrides"]), targets, weights, args.arm))
        results[case] = result
        record = {"case": case, "label": policy["label"], "policy_overrides": policy["overrides"],
                  **fertility_record(result), "target_fit": target_fit(result["moments"], targets, weights),
                  "elapsed_sec": result["elapsed_sec"]}
        records.append(record)
    if args.arm == "e5":
        outdir.mkdir(parents=True, exist_ok=True)
    write_diagnostics(results["baseline"]["sol"], results["baseline"]["P"], outdir / "standard")
    write_timing_supplement(outdir / "supplemental_timing", results["baseline"])
    for record in records[1:]:
        case_dir = outdir / "supplemental_timing" / str(record["case"])
        case_dir.mkdir(parents=True, exist_ok=True)
        write_hazard_delta(case_dir / "first_birth_realized_hazard_delta_vs_baseline.png", results["baseline"], results[str(record["case"])])
    (outdir / "policy_cases.json").write_text(json.dumps(jsonable({"source": str(args.source), "theta": theta, "target_set": E1_TARGET_SET, "targets": targets, "weights": weights, "cases": records}), indent=2, sort_keys=True) + "\n")
    flat = [flatten_record(record) for record in records]
    with (outdir / "policy_cases.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flat[0]))
        writer.writeheader(); writer.writerows(flat)
    write_policy_table(outdir / "policy_table.md", records)
    if args.arm == "e5":
        write_policy_births_summary(outdir / "policy_births_summary.csv", records)
    print(f"Wrote E2 policy packet to {outdir}", flush=True)


if __name__ == "__main__":
    main()
