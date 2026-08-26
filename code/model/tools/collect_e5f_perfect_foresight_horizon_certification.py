#!/usr/bin/env python3
"""Certify independently solved E5F perfect-foresight horizons.

Each input directory must be a complete output from
``run_e5f_perfect_foresight_transition.py`` containing exactly one horizon.
The collector is read-only with respect to those inputs.  It checks their
common scientific and code provenance, applies the numerical equilibrium and
terminal-convergence gates, and measures early-path stability against the
longest horizon.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_MARKET_TOLERANCE = 2e-4
DEFAULT_STABILITY_TOLERANCE = 1e-3
PROVENANCE_KEYS = (
    "selected_report_sha256",
    "selected_transition_sha256",
    "stationary_source_sha256",
    "selected_code_bundle_sha256",
    "perfect_foresight_solver_sha256",
    "calendar_bellman_sha256",
    "closed_endpoint_solver_sha256",
)
TERMINAL_KEYS = (
    "psi_child",
    "asset_price",
    "renter_price",
    "adult_population",
    "entry_flow_E",
    "B_over_E",
)
PATH_VARIABLES = (
    "asset_price",
    "renter_price",
    "birth_children_topcode_adjusted",
    "adult_population",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dirs", type=Path, nargs="+", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--comparison-periods", type=int, default=5)
    parser.add_argument(
        "--stability-tolerance", type=float, default=DEFAULT_STABILITY_TOLERANCE
    )
    return parser.parse_args()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, payload: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty CSV: {path}")
    fields = list(dict.fromkeys(key for row in rows for key in row))
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def load_run(run_dir: Path) -> dict[str, Any]:
    directory = run_dir.resolve()
    summary_path = directory / "summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("status") != "complete_isolated_perfect_foresight_diagnostic":
        raise RuntimeError(f"Incomplete transition packet: {summary_path}")
    horizons = dict(summary.get("horizons") or {})
    if len(horizons) != 1:
        raise RuntimeError(
            f"Each certification input must contain exactly one horizon: {directory}"
        )
    horizon = int(next(iter(horizons)))
    path_file = directory / f"horizon_{horizon:02d}" / "transition_path.csv"
    if not path_file.is_file():
        raise FileNotFoundError(path_file)
    with path_file.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != horizon:
        raise RuntimeError(
            f"Horizon {horizon} has {len(rows)} transition rows in {path_file}."
        )
    return {
        "directory": directory,
        "summary_path": summary_path,
        "summary_sha256": file_sha256(summary_path),
        "summary": summary,
        "horizon": horizon,
        "diagnostics": horizons[str(horizon)],
        "path_file": path_file,
        "path_sha256": file_sha256(path_file),
        "rows": rows,
    }


def common_contract_gate(runs: Sequence[dict[str, Any]]) -> dict[str, Any]:
    reference = runs[-1]["summary"]
    reference_provenance = dict(reference.get("provenance") or {})
    reference_terminal = dict(reference.get("terminal_steady_state") or {})
    checks: dict[str, bool] = {}
    for run in runs:
        summary = run["summary"]
        provenance = dict(summary.get("provenance") or {})
        terminal = dict(summary.get("terminal_steady_state") or {})
        missing_terminal = [key for key in TERMINAL_KEYS if key not in terminal]
        if missing_terminal:
            raise RuntimeError(
                f"Horizon {run['horizon']} predates the scaled-terminal schema; "
                f"missing terminal fields: {missing_terminal}."
            )
        for key in PROVENANCE_KEYS:
            checks[f"h{run['horizon']}_provenance_{key}"] = (
                provenance.get(key) == reference_provenance.get(key)
                and provenance.get(key) is not None
            )
        for key in TERMINAL_KEYS:
            checks[f"h{run['horizon']}_terminal_{key}"] = math.isclose(
                float(terminal[key]),
                float(reference_terminal[key]),
                rel_tol=0.0,
                abs_tol=2e-12,
            )
        checks[f"h{run['horizon']}_psi_persistence"] = math.isclose(
            float(summary["psi_persistence_per_four_year_period"]),
            float(reference["psi_persistence_per_four_year_period"]),
            rel_tol=0.0,
            abs_tol=0.0,
        )
    return {
        "status": "passed" if all(checks.values()) else "failed",
        "checks": checks,
        "provenance": {
            key: reference_provenance.get(key) for key in PROVENANCE_KEYS
        },
        "terminal_steady_state": {
            key: reference_terminal.get(key) for key in TERMINAL_KEYS
        },
        "psi_persistence_per_four_year_period": reference[
            "psi_persistence_per_four_year_period"
        ],
    }


def horizon_gate(run: dict[str, Any]) -> dict[str, Any]:
    summary = run["summary"]
    diagnostics = run["diagnostics"]
    terminal = dict(diagnostics.get("terminal_convergence") or {})
    numerical_checks = {
        "zero_shock_nesting": (summary.get("zero_shock_test") or {}).get("status")
        == "passed",
        "terminal_fixed_point": (
            (summary.get("terminal_steady_state") or {}).get("fixed_point_gates")
            or {}
        ).get("status")
        == "passed",
        "history_reproduction": (
            summary.get("history_reproduction_gate") or {}
        ).get("status")
        == "passed",
        "housing_market": bool(diagnostics.get("converged", False))
        and float(diagnostics["maximum_market_residual"])
        <= DEFAULT_MARKET_TOLERANCE,
        "policy_reproduction": float(
            diagnostics["maximum_policy_reproduction_error"]
        )
        <= 2e-8,
        "mass_accounting": float(diagnostics["maximum_mass_accounting_error"])
        <= 2e-8,
        "feasibility_projection": float(
            diagnostics["maximum_feasibility_projection_mass"]
        )
        <= 1e-6,
    }
    return {
        "horizon": run["horizon"],
        "numerical_status": (
            "passed" if all(numerical_checks.values()) else "failed"
        ),
        "numerical_checks": numerical_checks,
        "terminal_status": terminal.get("status"),
        "terminal_all_checks_pass": bool(terminal.get("all_checks_pass", False)),
        "iterations": int(diagnostics["iterations"]),
        "maximum_market_residual": float(
            diagnostics["maximum_market_residual"]
        ),
        "maximum_log_price_residual": float(
            diagnostics["maximum_log_price_residual"]
        ),
        "maximum_policy_reproduction_error": float(
            diagnostics["maximum_policy_reproduction_error"]
        ),
        "maximum_mass_accounting_error": float(
            diagnostics["maximum_mass_accounting_error"]
        ),
        "maximum_feasibility_projection_mass": float(
            diagnostics["maximum_feasibility_projection_mass"]
        ),
        "terminal_metrics": terminal.get("metrics"),
        "elapsed_seconds": float(diagnostics["elapsed_seconds"]),
        "summary_path": str(run["summary_path"]),
        "summary_sha256": run["summary_sha256"],
        "transition_path": str(run["path_file"]),
        "transition_path_sha256": run["path_sha256"],
    }


def stability_gate(
    runs: Sequence[dict[str, Any]], comparison_periods: int, tolerance: float
) -> dict[str, Any]:
    if len(runs) < 2:
        raise ValueError("Horizon certification requires at least two runs.")
    reference = runs[-1]
    records: list[dict[str, Any]] = []
    for candidate in runs[:-1]:
        common = min(
            int(comparison_periods), candidate["horizon"], reference["horizon"]
        )
        candidate_years = [
            int(row["calendar_year"]) for row in candidate["rows"][:common]
        ]
        reference_years = [
            int(row["calendar_year"]) for row in reference["rows"][:common]
        ]
        if candidate_years != reference_years:
            raise RuntimeError("Horizon paths do not share the same calendar dates.")
        for variable in PATH_VARIABLES:
            left = np.array(
                [float(row[variable]) for row in candidate["rows"][:common]]
            )
            right = np.array(
                [float(row[variable]) for row in reference["rows"][:common]]
            )
            scale = np.maximum(np.abs(right), 1e-12)
            records.append(
                {
                    "short_horizon": candidate["horizon"],
                    "reference_horizon": reference["horizon"],
                    "comparison_periods": common,
                    "variable": variable,
                    "maximum_absolute_gap": float(np.max(np.abs(left - right))),
                    "maximum_relative_gap": float(
                        np.max(np.abs(left - right) / scale)
                    ),
                }
            )
    maximum = max(float(row["maximum_relative_gap"]) for row in records)
    return {
        "status": "passed" if maximum <= float(tolerance) else "failed",
        "reference_horizon": reference["horizon"],
        "comparison_periods": int(comparison_periods),
        "relative_tolerance": float(tolerance),
        "maximum_relative_gap": maximum,
        "rows": records,
    }


def make_comparison_figure(
    runs: Sequence[dict[str, Any]], terminal: dict[str, Any], output_dir: Path
) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(10.5, 6.8), constrained_layout=True)
    specifications = (
        ("asset_price", "Asset price", axes[0, 0]),
        ("renter_price", "Renter price", axes[0, 1]),
        ("birth_children_topcode_adjusted", "Adjusted births", axes[1, 0]),
        ("adult_population", "Adult population", axes[1, 1]),
    )
    for run in runs:
        years = np.array([int(row["calendar_year"]) for row in run["rows"]])
        for variable, _, axis in specifications:
            values = np.array([float(row[variable]) for row in run["rows"]])
            axis.plot(years, values, marker="o", label=f"H={run['horizon']}")
    axes[0, 0].axhline(
        float(terminal["asset_price"]), color="0.35", linestyle="--", linewidth=1
    )
    axes[0, 1].axhline(
        float(terminal["renter_price"]), color="0.35", linestyle="--", linewidth=1
    )
    axes[1, 1].axhline(
        float(terminal["adult_population"]),
        color="0.35",
        linestyle="--",
        linewidth=1,
    )
    for _, title, axis in specifications:
        axis.set_title(title)
        axis.set_xlabel("Year")
        axis.grid(alpha=0.2)
        axis.legend(frameon=False)
    for suffix in ("png", "pdf"):
        figure.savefig(output_dir / f"horizon_comparison.{suffix}", dpi=180)
    plt.close(figure)


def main() -> None:
    args = parse_args()
    if int(args.comparison_periods) < 1:
        raise ValueError("--comparison-periods must be positive.")
    if not 0.0 < float(args.stability_tolerance) < 1.0:
        raise ValueError("--stability-tolerance must lie in (0,1).")
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise RuntimeError(f"Refusing to overwrite nonempty output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    runs = sorted(
        (load_run(path) for path in args.run_dirs), key=lambda run: run["horizon"]
    )
    if len({run["horizon"] for run in runs}) != len(runs):
        raise RuntimeError("Duplicate horizons were supplied to the collector.")
    contract = common_contract_gate(runs)
    horizon_results = [horizon_gate(run) for run in runs]
    stability = stability_gate(
        runs, int(args.comparison_periods), float(args.stability_tolerance)
    )
    longest_horizon_gate = horizon_results[-1]
    certified = (
        contract["status"] == "passed"
        and stability["status"] == "passed"
        and all(row["numerical_status"] == "passed" for row in horizon_results)
        and bool(longest_horizon_gate["terminal_all_checks_pass"])
    )
    payload = {
        "status": (
            "certified_horizon_stable_perfect_foresight_transition"
            if certified
            else "failed_perfect_foresight_certification"
        ),
        "promotion_status": "not_promoted",
        "interpretation": (
            "Conditional on the declared exogenous fertility-preference path, "
            "the closed national perfect-foresight equilibrium is certified only "
            "when every horizon clears housing markets and passes numerical-health "
            "gates, the longest horizon converges to the jointly solved terminal "
            "fixed point, and the horizons share a stable early path."
        ),
        "common_contract_gate": contract,
        "horizon_gates": horizon_results,
        "horizon_stability": {
            key: value for key, value in stability.items() if key != "rows"
        },
    }
    write_csv(output_dir / "horizon_gates.csv", horizon_results)
    write_csv(output_dir / "horizon_stability.csv", stability["rows"])
    write_json(output_dir / "certification.json", payload)
    make_comparison_figure(runs, contract["terminal_steady_state"], output_dir)
    print(
        f"PF_CERTIFICATION status={payload['status']} output={output_dir}",
        flush=True,
    )


if __name__ == "__main__":
    main()
