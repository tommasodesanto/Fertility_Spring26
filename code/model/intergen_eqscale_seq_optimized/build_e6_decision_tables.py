#!/usr/bin/env python3
"""Build comparable E5b/E6 calibration tables from certified collector reports."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from intergen_eqscale_seq_optimized.e5_profile import E5_DOMAIN  # noqa: E402


DEFAULT_BASELINE = (
    ROOT / "output/model/eqscale_seq_e5b_recalibration_20260725/report/results.json"
)
DEFAULT_OUTDIR = ROOT / "output/model/eqscale_seq_e6_decision_tables_20260727"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument(
        "--variant",
        action="append",
        default=[],
        metavar="LABEL=RESULTS_JSON",
        help="certified collector result; repeat for each E6 variant",
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def parse_variant(value: str) -> tuple[str, Path]:
    label, separator, path = value.partition("=")
    if not separator or not label.strip() or not path.strip():
        raise ValueError(f"variant must be LABEL=RESULTS_JSON, received {value!r}")
    return label.strip(), Path(path.strip())


def read_report(label: str, path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    winner = ((payload.get("winners") or {}).get("E1"))
    if not isinstance(winner, dict):
        raise ValueError(f"{path} does not contain winners.E1")
    rows = winner.get("target_fit")
    if not isinstance(rows, list) or len(rows) != 12:
        raise ValueError(f"{path} must contain the complete twelve-row target fit")
    names = [str(row["moment"]) for row in rows]
    if len(set(names)) != 12:
        raise ValueError(f"{path} target-fit moments are not unique")
    theta = winner.get("theta")
    if not isinstance(theta, dict):
        raise ValueError(f"{path} winner omits theta")
    metadata = payload.get("winner_metadata") or {}
    return {
        "label": label,
        "path": path,
        "payload": payload,
        "winner": winner,
        "rows": rows,
        "theta": theta,
        "metadata": metadata,
    }


def validate_common_contract(records: list[dict[str, Any]]) -> None:
    baseline_rows = {str(row["moment"]): row for row in records[0]["rows"]}
    for record in records[1:]:
        rows = {str(row["moment"]): row for row in record["rows"]}
        if set(rows) != set(baseline_rows):
            raise ValueError(f"{record['label']} uses a different target set")
        for name, baseline in baseline_rows.items():
            for key in ("target", "weight"):
                if not math.isclose(
                    float(rows[name][key]),
                    float(baseline[key]),
                    rel_tol=0.0,
                    abs_tol=1e-12,
                ):
                    raise ValueError(
                        f"{record['label']} changes {name} {key}: "
                        f"{rows[name][key]} versus {baseline[key]}"
                    )


def derived_parameter_rows(record: dict[str, Any]) -> list[dict[str, Any]]:
    """Fallback for the historical E5b report, which predates parameter CSVs."""
    theta = record["theta"]
    rows: list[dict[str, Any]] = []
    for name, lower, upper, transform in E5_DOMAIN:
        theta_name = "beta" if name == "beta_annual" else name
        estimate = float(theta[theta_name])
        if name == "beta_annual":
            estimate **= 0.25
        span = float(upper) - float(lower)
        distance = min(estimate - float(lower), float(upper) - estimate)
        near = distance <= 0.02 * span
        rows.append(
            {
                "parameter": name,
                "estimate": estimate,
                "lower_bound": float(lower),
                "upper_bound": float(upper),
                "transform": transform,
                "distance_to_nearest_bound": distance,
                "distance_as_share_of_range": distance / span,
                "near_bound_2pct": near,
                "near_bound_side": (
                    "lower" if near and estimate - float(lower) <= float(upper) - estimate
                    else "upper" if near
                    else ""
                ),
                "external_restriction": ">= 0.94" if name == "beta_annual" else "",
                "external_restriction_satisfied": (
                    estimate >= 0.94 if name == "beta_annual" else ""
                ),
            }
        )
    return rows


def read_parameter_rows(record: dict[str, Any]) -> list[dict[str, Any]]:
    path = record["path"].parent / "parameter_table_full.csv"
    if path.exists():
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
    else:
        rows = derived_parameter_rows(record)
    expected = int(record["metadata"].get("free_parameter_count", 10))
    if len(rows) != expected or len({str(row["parameter"]) for row in rows}) != expected:
        raise ValueError(
            f"{record['label']} must report exactly {expected} free parameters"
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def build_summary(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    baseline_loss = float(records[0]["winner"]["rank_loss"])
    keys = (
        "tfr",
        "childless_rate",
        "mean_age_first_birth",
        "share_first_births_age30plus",
        "own_family_gap",
        "own_rate",
        "aggregate_wealth_to_annual_gross_labor_earnings",
        "old_total_wealth_to_annual_income_p90_p50_7684",
    )
    rows = []
    for record in records:
        winner = record["winner"]
        moments = winner.get("moments") or {}
        theta = winner["theta"]
        row = {
            "variant": record["label"],
            "arm": record["metadata"].get("arm", winner.get("arm")),
            "loss": float(winner["rank_loss"]),
            "loss_change_vs_E5b": float(winner["rank_loss"]) - baseline_loss,
            "market_residual": float(winner["market_residual"]),
            "strict_converged": bool(winner["strict_converged"]),
            "tight_solve_seconds": float(winner.get("elapsed_sec", math.nan)),
            "beta_annual": float(theta["beta"]) ** 0.25,
            "eligible_chains": record["payload"].get("eligible_E1_chains"),
            "production_chains": record["payload"].get("production_chain_count"),
            "source": str(record["path"]),
        }
        row.update({key: float(moments.get(key, math.nan)) for key in keys})
        rows.append(row)
    return sorted(rows, key=lambda row: float(row["loss"]))


def build_target_comparison(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    baseline = {str(row["moment"]): row for row in records[0]["rows"]}
    rows = []
    for record in records:
        for target_row in record["rows"]:
            name = str(target_row["moment"])
            base = baseline[name]
            rows.append(
                {
                    "variant": record["label"],
                    "moment": name,
                    "target": float(target_row["target"]),
                    "model": float(target_row["model"]),
                    "model_change_vs_E5b": float(target_row["model"]) - float(base["model"]),
                    "gap": float(target_row["gap"]),
                    "absolute_gap": abs(float(target_row["gap"])),
                    "absolute_gap_change_vs_E5b": (
                        abs(float(target_row["gap"])) - abs(float(base["gap"]))
                    ),
                    "weight": float(target_row["weight"]),
                    "loss_contribution": float(target_row["loss_contribution"]),
                    "loss_contribution_change_vs_E5b": (
                        float(target_row["loss_contribution"])
                        - float(base["loss_contribution"])
                    ),
                }
            )
    return rows


def build_parameter_comparison(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_record = {
        record["label"]: {
            str(row["parameter"]): row for row in read_parameter_rows(record)
        }
        for record in records
    }
    baseline = by_record[records[0]["label"]]
    rows = []
    for record in records:
        for name, parameter in by_record[record["label"]].items():
            estimate = float(parameter["estimate"])
            baseline_estimate = (
                float(baseline[name]["estimate"]) if name in baseline else math.nan
            )
            rows.append(
                {
                    "variant": record["label"],
                    "parameter": name,
                    "estimate": estimate,
                    "estimate_change_vs_E5b": estimate - baseline_estimate,
                    "lower_bound": parameter["lower_bound"],
                    "upper_bound": parameter["upper_bound"],
                    "near_bound_2pct": parameter["near_bound_2pct"],
                    "near_bound_side": parameter["near_bound_side"],
                    "external_restriction": parameter.get("external_restriction", ""),
                    "external_restriction_satisfied": parameter.get(
                        "external_restriction_satisfied", ""
                    ),
                }
            )
    return rows


def main() -> None:
    args = parse_args()
    variants = [parse_variant(value) for value in args.variant]
    records = [read_report("E5b", args.baseline)]
    records.extend(read_report(label, path) for label, path in variants)
    validate_common_contract(records)
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "variant_summary_ranked.csv", build_summary(records))
    write_csv(args.outdir / "target_fit_comparison_full.csv", build_target_comparison(records))
    write_csv(args.outdir / "parameter_comparison_full.csv", build_parameter_comparison(records))
    (args.outdir / "README.md").write_text(
        "\n".join(
            [
                "# E6 certified comparison tables",
                "",
                "These tables compare only strict collector winners on the common signed",
                "twelve-row target system. The baseline is E5b. Model deltas, absolute-gap",
                "deltas, loss-contribution deltas, all ten free parameters, bounds, and",
                "external-restriction checks are reported without running the model.",
                "",
                "The historical E5b report predates `parameter_table_full.csv`; its table",
                "is reconstructed from the signed E5 domain and winner vector.",
            ]
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
