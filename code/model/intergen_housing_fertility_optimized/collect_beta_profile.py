#!/usr/bin/env python3
"""Collect a conditional-beta diagnostic for the repaired one-shot objective."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

from .calibration_collect import repeat_is_exact
from .new_moment_profile import new_moment_target_system


EXPECTED_BETAS = (0.98, 0.99, 0.995, 0.998, 0.9995, 0.9997, 0.9999)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--expected-betas",
        type=float,
        nargs="+",
        default=EXPECTED_BETAS,
        help="Annual beta cells expected in the run (two chains per cell).",
    )
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    expected_betas = tuple(float(value) for value in args.expected_betas)
    if not expected_betas or len(set(expected_betas)) != len(expected_betas):
        raise ValueError("--expected-betas must contain distinct values")
    system = new_moment_target_system()
    grouped: dict[float, list[tuple[Path, dict[str, Any], dict[str, Any]]]] = {}
    task_rows: list[dict[str, Any]] = []
    for summary_path in sorted(args.run_root.glob("task_*/summary.json")):
        summary = json.loads(summary_path.read_text())
        metadata = dict(summary.get("metadata") or {})
        fingerprint = metadata.get("target_fingerprint")
        if fingerprint != system.fingerprint:
            raise RuntimeError(f"target fingerprint mismatch in {summary_path.parent.name}")
        beta = float(metadata["fixed_beta_annual"])
        repeats = list(summary.get("tight_repeats") or [])
        exact = repeat_is_exact(repeats)
        tight = repeats[0] if exact else None
        eligible = bool(exact and summary.get("eligible"))
        row = {
            "task": summary_path.parent.name,
            "beta_annual": beta,
            "method": metadata.get("method"),
            "cases_completed": summary.get("n_cases_completed"),
            "status": summary.get("status"),
            "eligible": eligible,
            "strict_loss": tight.get("rank_loss") if tight else None,
            "market_residual": tight.get("market_residual") if tight else None,
        }
        task_rows.append(row)
        if eligible:
            grouped.setdefault(beta, []).append((summary_path.parent, summary, tight))

    if len(task_rows) != 2 * len(expected_betas):
        write_csv(args.output / "all_tasks.csv", task_rows)
        raise RuntimeError(
            f"expected {2 * len(expected_betas)} completed tasks, found {len(task_rows)}"
        )
    profile_rows: list[dict[str, Any]] = []
    selected_target_fit_rows: list[dict[str, Any]] = []
    selected_parameter_rows: list[dict[str, Any]] = []
    selected_records: dict[str, Any] = {}
    for expected_beta in expected_betas:
        matching_beta = next(
            (beta for beta in grouped if math.isclose(beta, expected_beta, abs_tol=1e-12)),
            None,
        )
        if matching_beta is None or len(grouped[matching_beta]) != 2:
            raise RuntimeError(f"beta={expected_beta} lacks two eligible strict chains")
        selected_dir, selected_summary, selected = min(
            grouped[matching_beta],
            key=lambda item: float(item[2]["rank_loss"]),
        )
        parameters = list(selected.get("parameters") or [])
        if sum(row.get("role") == "estimated" for row in parameters) != 13:
            raise RuntimeError(f"beta={expected_beta} does not report 13 nuisance estimates")
        if sum(row.get("role") == "profile_fixed" for row in parameters) != 1:
            raise RuntimeError(f"beta={expected_beta} does not mark beta profile-fixed")
        fit = list(selected.get("target_fit") or [])
        if len(fit) != system.count:
            raise RuntimeError(f"beta={expected_beta} lacks the complete target-fit table")
        strict_loss = float(selected["rank_loss"])
        if not math.isclose(
            sum(float(row["loss_contribution"]) for row in fit),
            strict_loss,
            rel_tol=1e-12,
            abs_tol=1e-12,
        ):
            raise RuntimeError(f"beta={expected_beta} target-fit contributions do not sum")
        theta = dict(selected["theta"])
        moments = dict(selected["moments"])
        components = dict(selected.get("wealth_components") or {})
        selected_target_fit_rows.extend(
            {"beta_annual": expected_beta, **fit_row} for fit_row in fit
        )
        selected_parameter_rows.extend(
            {"beta_annual": expected_beta, **parameter_row}
            for parameter_row in parameters
        )
        row = {
            "beta_annual": expected_beta,
            "beta_period": float(theta["beta"]),
            "selected_task": selected_dir.name,
            "strict_loss": strict_loss,
            "market_residual": float(selected["market_residual"]),
            "wealth_to_gross_labor_earnings": float(
                moments["aggregate_wealth_to_annual_gross_labor_earnings"]
            ),
            "aggregate_wealth": float(components.get("aggregate_wealth", math.nan)),
            "annual_gross_labor_earnings": float(
                components.get("aggregate_annual_gross_labor_earnings", math.nan)
            ),
            "wealth_to_gross_labor_earnings_26_35": float(
                components.get(
                    "aggregate_wealth_to_annual_gross_labor_earnings_26_35",
                    math.nan,
                )
            ),
            "wealth_to_gross_labor_earnings_36_45": float(
                components.get(
                    "aggregate_wealth_to_annual_gross_labor_earnings_36_45",
                    math.nan,
                )
            ),
            "wealth_to_gross_labor_earnings_46_55": float(
                components.get(
                    "aggregate_wealth_to_annual_gross_labor_earnings_46_55",
                    math.nan,
                )
            ),
            "wealth_to_gross_labor_earnings_56_65": float(
                components.get(
                    "aggregate_wealth_to_annual_gross_labor_earnings_56_65",
                    math.nan,
                )
            ),
            "annual_bequest_flow": float(
                components.get("annual_bequest_flow", math.nan)
            ),
            "bequest_flow_to_wealth": float(
                moments["annual_bequest_flow_to_aggregate_wealth"]
            ),
            "tfr": float(moments["tfr"]),
            "childless_rate": float(moments["childless_rate"]),
            "theta0": float(theta["theta0"]),
            "theta1": float(theta["theta1"]),
            "psi_child": float(theta["psi_child"]),
            "kappa_fert": float(theta["kappa_fert"]),
        }
        profile_rows.append(row)
        selected_records[f"{expected_beta:.4f}"] = {
            "row": row,
            "selected": selected,
            "selected_summary": selected_summary,
        }

    write_csv(args.output / "all_tasks.csv", task_rows)
    write_csv(args.output / "beta_profile.csv", profile_rows)
    write_csv(args.output / "selected_target_fits.csv", selected_target_fit_rows)
    write_csv(args.output / "selected_parameters.csv", selected_parameter_rows)
    best = min(profile_rows, key=lambda row: float(row["strict_loss"]))
    payload = {
        "status": "conditional_beta_profile_not_a_calibration",
        "target_system": system.name,
        "target_fingerprint": system.fingerprint,
        "profile_fixed_parameter": "beta_annual",
        "profile_values": list(expected_betas),
        "nuisance_free_parameter_count": 13,
        "hard_moment_count": system.count,
        "best_profile_row": best,
        "profile_rows": profile_rows,
        "selected_records": selected_records,
    }
    (args.output / "results.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    lines = [
        "# Conditional annual-beta profile",
        "",
        "This is a diagnostic profile of the timing-repaired one-shot objective,",
        "not a reweighted or promoted calibration.",
        "Each beta cell reoptimizes the other 13 parameters with two independent chains",
        "and requires two exact strict repeats per chain.",
        "",
        f"- Best profiled beta: `{float(best['beta_annual']):.6g}`",
        f"- Best canonical strict loss: `{float(best['strict_loss']):.12g}`",
        f"- Wealth / gross labor earnings there: `{float(best['wealth_to_gross_labor_earnings']):.6g}`",
        f"- TFR there: `{float(best['tfr']):.6g}`",
        "",
        "The complete curve is in `beta_profile.csv`; complete target-fit and",
        "parameter tables are in `selected_target_fits.csv` and",
        "`selected_parameters.csv`; task status is in `all_tasks.csv`.",
    ]
    (args.output / "REPORT.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({key: value for key, value in payload.items() if key != "selected_records"}, indent=2))


if __name__ == "__main__":
    main()
