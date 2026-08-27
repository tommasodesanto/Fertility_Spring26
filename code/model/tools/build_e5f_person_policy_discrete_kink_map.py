#!/usr/bin/env python3
"""Audit and plot the H128 deterministic-tenure path-blend diagnostic."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


MARKET_TOLERANCE = 2.0e-4
FISCAL_TOLERANCE = 2.5e-5
EXPECTED_CASE = "rebated-tax2-reform"
OUTPUT_PREFIX = "e5f_pf_person_policy_rebated-tax2-reform_20260827a_kink_"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coarse-manifest", type=Path, required=True)
    parser.add_argument("--bracket-manifest", type=Path, required=True)
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-driver-sha256", required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_manifest(path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    cases = manifest.get("cases")
    if not isinstance(cases, list) or not cases:
        raise RuntimeError(f"Manifest has no cases: {path}")
    return manifest, cases


def find_output_for_seed(results_root: Path, seed_sha256: str) -> Path:
    matches: list[Path] = []
    for candidate in sorted(results_root.glob(f"{OUTPUT_PREFIX}*_smoke")):
        contract_path = candidate / "initial_path_seed_contract.json"
        if not contract_path.exists():
            continue
        contract = json.loads(contract_path.read_text(encoding="utf-8"))
        if contract.get("source_sha256") == seed_sha256:
            matches.append(candidate)
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one output for seed {seed_sha256}, found {len(matches)}"
        )
    return matches[0]


def maximum_row(
    rows: list[dict[str, str]], column: str
) -> tuple[float, dict[str, str]]:
    row = max(rows, key=lambda item: abs(float(item[column])))
    return abs(float(row[column])), row


def audit_case(
    case: dict[str, Any], results_root: Path, expected_driver_sha256: str
) -> dict[str, Any]:
    lambda_value = float(case["lambda"])
    expected_seed_sha256 = str(case["sha256"])
    output_dir = find_output_for_seed(results_root, expected_seed_sha256)
    summary_path = output_dir / "summary.json"
    transition_path = output_dir / "transition_path.csv"
    seed_contract_path = output_dir / "initial_path_seed_contract.json"
    launch_contract_path = output_dir / "launch_contract.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    seed_contract = json.loads(seed_contract_path.read_text(encoding="utf-8"))
    launch_contract = json.loads(launch_contract_path.read_text(encoding="utf-8"))
    rows = read_csv(transition_path)

    horizon = int(summary["horizon"])
    periods = [int(row["period"]) for row in rows]
    maximum_market, market_row = maximum_row(rows, "relative_market_residual")
    maximum_fiscal, fiscal_row = maximum_row(rows, "government_budget_residual")
    maximum_feasibility = max(
        abs(float(row["feasibility_frontier_projection_mass"])) for row in rows
    )
    year_rows = {int(row["calendar_year"]): row for row in rows}
    required_years = (2351, 2363, 2367)
    if any(year not in year_rows for year in required_years):
        raise RuntimeError(f"Missing diagnostic year in {transition_path}")

    accounting_gates = summary.get("accounting_gates") or {}
    terminal_convergence = summary.get("terminal_convergence") or {}
    checks = {
        "complete_status": summary.get("status")
        == "complete_unpromoted_person_demography_policy_path",
        "correct_case": summary.get("case") == EXPECTED_CASE,
        "correct_driver": (summary.get("source_hashes") or {}).get("driver")
        == expected_driver_sha256,
        "launch_driver_matches": (launch_contract.get("hashes") or {}).get(
            "E5F_PF_PERSON_POLICY_EXPECTED_DRIVER_SHA256"
        )
        == expected_driver_sha256,
        "seed_hash_matches": seed_contract.get("source_sha256")
        == expected_seed_sha256,
        "contiguous_path": len(rows) == horizon and periods == list(range(horizon)),
        "market_residual_matches_summary": math.isclose(
            maximum_market,
            float(summary["maximum_market_residual"]),
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ),
        "fiscal_residual_matches_summary": math.isclose(
            maximum_fiscal,
            float(summary["maximum_fiscal_residual"]),
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ),
        "zero_feasibility_adjustment": maximum_feasibility == 0.0,
        "accounting_gates_pass": bool(accounting_gates)
        and all(bool(value) for value in accounting_gates.values()),
        "terminal_gates_pass": bool(terminal_convergence.get("all_checks_pass")),
        "history_reproduction_pass": summary.get("history_reproduction_status")
        == "passed",
    }
    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise RuntimeError(f"Audit failed for lambda={lambda_value}: {failed}")

    result: dict[str, Any] = {
        "lambda": lambda_value,
        "output_dir": str(output_dir.resolve()),
        "input_path_sha256": expected_seed_sha256,
        "summary_sha256": sha256(summary_path),
        "transition_path_sha256": sha256(transition_path),
        "maximum_market_residual": maximum_market,
        "maximum_market_residual_year": int(market_row["calendar_year"]),
        "maximum_fiscal_residual": maximum_fiscal,
        "maximum_fiscal_residual_year": int(fiscal_row["calendar_year"]),
        "market_gate_pass": maximum_market <= MARKET_TOLERANCE,
        "fiscal_gate_pass": maximum_fiscal <= FISCAL_TOLERANCE,
        "path_gate_pass": maximum_market <= MARKET_TOLERANCE
        and maximum_fiscal <= FISCAL_TOLERANCE,
        "maximum_feasibility_adjustment": maximum_feasibility,
        "all_source_accounting_terminal_checks_pass": True,
    }
    for year in required_years:
        row = year_rows[year]
        result[f"asset_price_{year}"] = float(row["asset_price"])
        result[f"equal_transfer_{year}"] = float(
            row["equal_transfer_period_units"]
        )
        result[f"owner_rate_{year}"] = float(row["owner_rate"])
        result[f"market_residual_{year}"] = float(
            row["relative_market_residual"]
        )
        result[f"fiscal_residual_{year}"] = float(
            row["government_budget_residual"]
        )
    return result


def branch_switch(results: list[dict[str, Any]]) -> dict[str, Any]:
    sorted_results = sorted(results, key=lambda row: float(row["lambda"]))
    differences = [
        abs(
            float(right["owner_rate_2351"])
            - float(left["owner_rate_2351"])
        )
        for left, right in zip(sorted_results[:-1], sorted_results[1:])
    ]
    index = max(range(len(differences)), key=differences.__getitem__)
    left = sorted_results[index]
    right = sorted_results[index + 1]
    return {
        "lower_lambda": float(left["lambda"]),
        "upper_lambda": float(right["lambda"]),
        "owner_rate_jump_percentage_points": 100.0
        * (
            float(right["owner_rate_2351"])
            - float(left["owner_rate_2351"])
        ),
        "asset_price_change": float(right["asset_price_2351"])
        - float(left["asset_price_2351"]),
        "equal_transfer_change": float(right["equal_transfer_2351"])
        - float(left["equal_transfer_2351"]),
        "pre_switch_market_residual": float(left["maximum_market_residual"]),
        "pre_switch_fiscal_residual": float(left["maximum_fiscal_residual"]),
        "post_switch_market_residual": float(right["maximum_market_residual"]),
        "post_switch_fiscal_residual": float(right["maximum_fiscal_residual"]),
        "pre_switch_path_gate_pass": bool(left["path_gate_pass"]),
        "post_switch_path_gate_pass": bool(right["path_gate_pass"]),
    }


def write_results(path: Path, results: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)


def build_figure(
    output_dir: Path,
    results: list[dict[str, Any]],
    switch: dict[str, Any],
) -> None:
    lambdas = [float(row["lambda"]) for row in results]
    lower = float(switch["lower_lambda"])
    upper = float(switch["upper_lambda"])
    figure, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), constrained_layout=True)

    axes[0, 0].plot(
        lambdas,
        [float(row["maximum_market_residual"]) for row in results],
        marker="o",
    )
    axes[0, 0].axhline(MARKET_TOLERANCE, color="tab:red", linestyle="--")
    axes[0, 0].set(title="Maximum housing-market residual", ylabel="Absolute residual")

    axes[0, 1].plot(
        lambdas,
        [float(row["maximum_fiscal_residual"]) for row in results],
        marker="o",
    )
    axes[0, 1].axhline(FISCAL_TOLERANCE, color="tab:red", linestyle="--")
    axes[0, 1].set(title="Maximum fiscal residual", ylabel="Absolute residual")

    axes[1, 0].plot(
        lambdas,
        [100.0 * float(row["owner_rate_2351"]) for row in results],
        marker="o",
    )
    axes[1, 0].set(title="Ownership at the switching date (2351)", ylabel="Percent")

    prices = [float(row["asset_price_2351"]) for row in results]
    price_changes = [1.0e6 * (price - prices[0]) for price in prices]
    owners = [100.0 * float(row["owner_rate_2351"]) for row in results]
    axes[1, 1].plot(price_changes, owners, marker="o")
    for lambda_value, price_change, owner in zip(lambdas, price_changes, owners):
        if lower <= lambda_value <= upper:
            axes[1, 1].annotate(
                f"$\\lambda={lambda_value:.3f}$",
                (price_change, owner),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8,
            )
    axes[1, 1].set(
        title="Discrete ownership response to a tiny price change",
        xlabel="Asset-price change in 2351 (millionths from $\\lambda=0.05$)",
        ylabel="Ownership rate in 2351 (percent)",
    )

    for axis in axes.flat:
        axis.axvspan(lower, upper, color="0.85", alpha=0.45) if axis is not axes[1, 1] else None
        axis.grid(alpha=0.2)
        if axis is not axes[1, 1]:
            axis.set_xlabel("Blend weight $\\lambda$")
    figure.suptitle("H128 deterministic-tenure discontinuity map")
    figure.savefig(output_dir / "kink_map.png", dpi=220)
    figure.savefig(output_dir / "kink_map.pdf")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    coarse_manifest, coarse_cases = load_manifest(args.coarse_manifest.resolve())
    bracket_manifest, bracket_cases = load_manifest(args.bracket_manifest.resolve())
    all_cases = coarse_cases + bracket_cases
    input_hashes = [str(case["sha256"]) for case in all_cases]
    if len(input_hashes) != len(set(input_hashes)):
        raise RuntimeError("Duplicate input-path hash across manifests")

    results = sorted(
        [
            audit_case(case, args.results_root.resolve(), args.expected_driver_sha256)
            for case in all_cases
        ],
        key=lambda row: float(row["lambda"]),
    )
    switch = branch_switch(results)
    if switch["pre_switch_path_gate_pass"] or switch["post_switch_path_gate_pass"]:
        raise RuntimeError("A switch-bracketing path unexpectedly passes both gates")
    write_results(output_dir / "results.csv", results)
    build_figure(output_dir, results, switch)

    summary = {
        "status": "complete_unaccepted_deterministic_tenure_kink_map",
        "accepted_policy_transition": False,
        "all_source_accounting_terminal_checks_pass": all(
            bool(row["all_source_accounting_terminal_checks_pass"])
            for row in results
        ),
        "all_evaluated_paths_pass_market_and_fiscal_gates": all(
            bool(row["path_gate_pass"]) for row in results
        ),
        "any_evaluated_path_passes_market_and_fiscal_gates": any(
            bool(row["path_gate_pass"]) for row in results
        ),
        "evaluated_case_count": len(results),
        "market_tolerance": MARKET_TOLERANCE,
        "fiscal_tolerance": FISCAL_TOLERANCE,
        "expected_driver_sha256": args.expected_driver_sha256,
        "coarse_manifest": str(args.coarse_manifest.resolve()),
        "coarse_manifest_sha256": sha256(args.coarse_manifest.resolve()),
        "bracket_manifest": str(args.bracket_manifest.resolve()),
        "bracket_manifest_sha256": sha256(args.bracket_manifest.resolve()),
        "branch_switch": switch,
        "conclusion": (
            "Along the audited one-dimensional blend between the immutable best "
            "path and the adverse low-damping branch, the pure deterministic-tenure "
            "response switches discontinuously and neither bracketing side clears "
            "the declared market and fiscal gates. This diagnoses the failed local "
            "fixed-point update; it is not a global proof that no pure equilibrium "
            "exists elsewhere."
        ),
        "next_options": [
            "Mixed or convexified tenure choice at the switching mass.",
            "Finer asset-grid verification with both policy cases rerun.",
            "Positive tenure-choice dispersion as an identified calibration change.",
        ],
        "results_csv_sha256": sha256(output_dir / "results.csv"),
        "kink_map_png_sha256": sha256(output_dir / "kink_map.png"),
        "kink_map_pdf_sha256": sha256(output_dir / "kink_map.pdf"),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
