#!/usr/bin/env python3
"""Collect and report the best strict-converged intergenerational calibration."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable

try:
    from intergen_housing_fertility.calibration import get_target_set as _get_target_set
    from intergen_housing_fertility.production_profile import (
        PRODUCTION_SEARCH_BOUNDS as _PRODUCTION_SEARCH_BOUNDS,
        PRODUCTION_TARGET_SET as _PRODUCTION_TARGET_SET,
    )

    CONSTANT_SOURCE = "active Python modules"
except ImportError:
    # Keep collection usable on a lightweight login machine that lacks the
    # numerical model dependencies. These literals mirror the named active
    # source objects and are checked against them whenever imports are available.
    _get_target_set = None
    _PRODUCTION_TARGET_SET = "candidate_replacement_post_audit_v1"
    _PRODUCTION_SEARCH_BOUNDS = (
        ("beta_annual", 0.940, 0.995),
        ("alpha_cons", 0.400, 0.950),
        ("c_bar_0", 0.080, 1.280),
        ("c_bar_n", 0.050, 1.500),
        ("h_bar_0", 1.000, 6.000),
        ("h_bar_jump", 0.050, 2.500),
        ("h_bar_n", 0.020, 2.000),
        ("psi_child", 0.000, 0.350),
        ("kappa_fert", 1.000, 12.000),
        ("tenure_choice_kappa", 0.000, 0.120),
        ("chi", 0.400, 1.150),
        ("theta0", 0.000, 2.000),
        ("theta_n", 0.000, 1.500),
    )
    CONSTANT_SOURCE = "standard-library fallback mirroring active source"


PRODUCTION_TARGET_SET = _PRODUCTION_TARGET_SET
PRODUCTION_SEARCH_BOUNDS = _PRODUCTION_SEARCH_BOUNDS

FALLBACK_TARGETS = {
    "tfr": 1.918,
    "childless_rate": 0.188,
    "own_rate": 0.57547241,
    "own_family_gap": 0.16766167,
    "housing_increment_0to1": 0.66443467,
    "old_parent_childless_nonhousing_wealth_to_income_gap_6575": 1.00744952,
    "prime30_55_childless_renter_mean_rooms": 3.80528810,
    "prime30_55_childless_owner_share_rooms_ge6": 0.59613112,
    "old_nonhousing_wealth_to_income_median_6575": 2.23046078,
    "young_childless_renter_liquid_wealth_to_annual_gross_income_2535": 0.17922556,
    "prime30_55_childless_owner_minus_renter_mean_rooms": 2.41876173,
    "old_age_own_rate": 0.76426097,
    "own_rate_2534": 0.34116609,
    "prime30_55_parent_3plus_minus_1to2_mean_rooms": 0.36769955881,
}
FALLBACK_WEIGHTS = {
    "tfr": 20.0,
    "childless_rate": 20.0,
    "own_rate": 100.0,
    "own_family_gap": 45.0,
    "housing_increment_0to1": 14.0,
    "old_parent_childless_nonhousing_wealth_to_income_gap_6575": 2.0,
    "prime30_55_childless_renter_mean_rooms": 6.0,
    "prime30_55_childless_owner_share_rooms_ge6": 25.0,
    "old_nonhousing_wealth_to_income_median_6575": 0.8,
    "young_childless_renter_liquid_wealth_to_annual_gross_income_2535": 12.0,
    "prime30_55_childless_owner_minus_renter_mean_rooms": 12.0,
    "old_age_own_rate": 160.0,
    "own_rate_2534": 80.0,
    "prime30_55_parent_3plus_minus_1to2_mean_rooms": 8.0,
}


ROOMS_MOMENT = "aggregate_mean_occupied_rooms_18_85"
ROOMS_TARGET = 5.779970481941968
ROOMS_WEIGHT = 6.0
H0_BOUND = ("H0", 1.0, 20.0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-root",
        type=Path,
        required=True,
        help="Collected task tree containing cases.jsonl and/or best.json files.",
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--incumbent",
        type=Path,
        help="Optional incumbent JSON (a candidate record or an object with a 'best' record).",
    )
    return parser.parse_args()


def finite_number(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def eligible(record: Any) -> bool:
    return bool(
        isinstance(record, dict)
        and record.get("status") == "ok"
        and record.get("strict_converged") is True
        and finite_number(record.get("rank_loss"))
    )


def candidate_fingerprint(record: dict[str, Any]) -> str:
    payload = {
        "rank_loss": record.get("rank_loss"),
        "theta": record.get("theta"),
        "moments": record.get("moments"),
        "status": record.get("status"),
        "strict_converged": record.get("strict_converged"),
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=True)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def iter_jsonl(path: Path) -> Iterable[tuple[dict[str, Any], int]]:
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            payload = json.loads(line)
            if isinstance(payload, dict):
                yield payload, line_number


def scan_results(root: Path) -> tuple[list[dict[str, Any]], dict[str, int]]:
    if not root.is_dir():
        raise FileNotFoundError(f"results root is not a directory: {root}")

    found: list[dict[str, Any]] = []
    counts = {
        "files_scanned": 0,
        "records_parsed": 0,
        "parse_errors": 0,
        "eligible_records": 0,
        "unique_eligible_candidates": 0,
    }
    seen: set[str] = set()

    paths = sorted({*root.rglob("cases.jsonl"), *root.rglob("best.json")})
    for path in paths:
        counts["files_scanned"] += 1
        try:
            if path.name == "cases.jsonl":
                records = iter_jsonl(path)
            else:
                payload = json.loads(path.read_text())
                if not isinstance(payload, dict):
                    raise ValueError("top-level JSON value is not an object")
                record = payload.get("best", payload)
                if not isinstance(record, dict):
                    raise ValueError("'best' is not an object")
                records = ((record, 0),)

            for record, line_number in records:
                counts["records_parsed"] += 1
                if not eligible(record):
                    continue
                counts["eligible_records"] += 1
                fingerprint = candidate_fingerprint(record)
                if fingerprint in seen:
                    continue
                seen.add(fingerprint)
                found.append(
                    {
                        "record": record,
                        "source": str(path.resolve()),
                        "line": line_number or None,
                        "pool": "refinement",
                        "fingerprint": fingerprint,
                    }
                )
        except (json.JSONDecodeError, OSError, ValueError):
            counts["parse_errors"] += 1

    counts["unique_eligible_candidates"] = len(found)
    return found, counts


def load_incumbent(path: Path | None) -> dict[str, Any] | None:
    if path is None:
        return None
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError("incumbent JSON must contain an object")
    record = payload.get("best", payload)
    if not isinstance(record, dict):
        raise ValueError("incumbent 'best' value must be an object")
    if not eligible(record):
        raise ValueError("incumbent is not a finite, status-ok, strict-converged candidate")
    return {
        "record": record,
        "source": str(path.resolve()),
        "line": None,
        "pool": "incumbent",
        "fingerprint": candidate_fingerprint(record),
    }


def active_targets() -> tuple[dict[str, float], dict[str, float]]:
    if _get_target_set is None:
        targets, weights = dict(FALLBACK_TARGETS), dict(FALLBACK_WEIGHTS)
    else:
        targets, weights = _get_target_set(PRODUCTION_TARGET_SET)
        targets = dict(targets)
        weights = dict(weights)
        if targets != FALLBACK_TARGETS or weights != FALLBACK_WEIGHTS:
            raise RuntimeError("reporting constants have drifted from the active target set")
    targets[ROOMS_MOMENT] = ROOMS_TARGET
    weights[ROOMS_MOMENT] = ROOMS_WEIGHT
    if len(targets) != 15 or set(targets) != set(weights):
        raise RuntimeError("active combined target system is not the expected 15-moment system")
    return targets, weights


def target_rows(
    record: dict[str, Any], targets: dict[str, float], weights: dict[str, float]
) -> list[dict[str, Any]]:
    moments = record.get("moments")
    if not isinstance(moments, dict):
        raise ValueError("selected candidate has no moments object")
    missing = [name for name in targets if not finite_number(moments.get(name))]
    if missing:
        raise ValueError(f"selected candidate lacks finite active moments: {missing}")

    rows: list[dict[str, Any]] = []
    for name, target_value in targets.items():
        model_value = float(moments[name])
        gap = model_value - float(target_value)
        weight = float(weights[name])
        rows.append(
            {
                "moment": name,
                "target": float(target_value),
                "model": model_value,
                "gap_model_minus_target": gap,
                "weight": weight,
                "loss_contribution": weight * gap * gap,
            }
        )
    return rows


def parameter_rows(record: dict[str, Any]) -> list[dict[str, Any]]:
    theta = record.get("theta")
    if not isinstance(theta, dict):
        raise ValueError("selected candidate has no theta object")

    rows: list[dict[str, Any]] = []
    for name, lower, upper in (*PRODUCTION_SEARCH_BOUNDS, H0_BOUND):
        model_name = "beta" if name == "beta_annual" else name
        if not finite_number(theta.get(model_name)):
            raise ValueError(f"selected candidate lacks finite parameter {model_name!r}")
        model_value = float(theta[model_name])
        estimate = model_value ** 0.25 if name == "beta_annual" else model_value
        span = float(upper - lower)
        lower_distance = estimate - float(lower)
        upper_distance = float(upper) - estimate
        normalized_position = lower_distance / span
        near_lower = lower_distance <= 0.02 * span
        near_upper = upper_distance <= 0.02 * span
        if near_lower and near_upper:
            bound_side = "both"
        elif near_lower:
            bound_side = "lower"
        elif near_upper:
            bound_side = "upper"
        else:
            bound_side = ""
        rows.append(
            {
                "parameter": name,
                "model_key": model_name,
                "estimate": estimate,
                "lower_bound": float(lower),
                "upper_bound": float(upper),
                "normalized_position": normalized_position,
                "distance_to_lower": lower_distance,
                "distance_to_upper": upper_distance,
                "near_bound_2pct": near_lower or near_upper,
                "near_bound_side": bound_side,
            }
        )
    if len(rows) != 14:
        raise RuntimeError("active search vector is not the expected 14-parameter vector")
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def format_cell(value: Any) -> str:
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.10g}"
    return str(value).replace("|", "\\|")


def write_markdown(path: Path, title: str, rows: list[dict[str, Any]]) -> None:
    columns = list(rows[0])
    lines = [f"# {title}", "", "| " + " | ".join(columns) + " |"]
    lines.append("| " + " | ".join("---" for _ in columns) + " |")
    lines.extend(
        "| " + " | ".join(format_cell(row[column]) for column in columns) + " |"
        for row in rows
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    candidates, counts = scan_results(args.results_root.resolve())
    incumbent = load_incumbent(args.incumbent.resolve() if args.incumbent else None)
    if incumbent is not None and all(
        item["fingerprint"] != incumbent["fingerprint"] for item in candidates
    ):
        candidates.append(incumbent)
    if not candidates:
        raise RuntimeError("no finite, status-ok, strict-converged candidate found")

    selected = min(candidates, key=lambda item: float(item["record"]["rank_loss"]))
    best = selected["record"]
    targets, weights = active_targets()
    fit = target_rows(best, targets, weights)
    parameters = parameter_rows(best)
    recomputed_loss = sum(float(row["loss_contribution"]) for row in fit)
    reported_loss = float(best["rank_loss"])
    incumbent_loss = float(incumbent["record"]["rank_loss"]) if incumbent else None
    improvement = incumbent_loss - reported_loss if incumbent_loss is not None else None

    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "best.json").write_text(json.dumps(best, indent=2, sort_keys=True) + "\n")
    write_csv(args.outdir / "calibration_target_fit.csv", fit)
    write_markdown(args.outdir / "calibration_target_fit.md", "Calibration Target Fit", fit)
    write_csv(args.outdir / "calibration_parameters.csv", parameters)
    write_markdown(args.outdir / "calibration_parameters.md", "Calibration Parameters", parameters)

    summary = {
        "results_root": str(args.results_root.resolve()),
        "target_set": PRODUCTION_TARGET_SET,
        "constant_source": CONSTANT_SOURCE,
        "target_count": len(fit),
        "parameter_count": len(parameters),
        "selected_pool": selected["pool"],
        "selected_source": selected["source"],
        "selected_source_line": selected["line"],
        "selected_fingerprint_sha256": selected["fingerprint"],
        "reported_rank_loss": reported_loss,
        "recomputed_rank_loss": recomputed_loss,
        "reported_minus_recomputed_loss": reported_loss - recomputed_loss,
        "market_residual": best.get("market_residual"),
        "status": best.get("status"),
        "strict_converged": best.get("strict_converged"),
        "incumbent_source": incumbent["source"] if incumbent else None,
        "incumbent_rank_loss": incumbent_loss,
        "improvement_vs_incumbent": improvement,
        "improvement_percent_vs_incumbent": (
            100.0 * improvement / incumbent_loss
            if improvement is not None and incumbent_loss != 0.0
            else None
        ),
        "scan_counts": counts,
    }
    (args.outdir / "calibration_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )

    improvement_text = (
        "not computed (no incumbent supplied)"
        if improvement is None
        else f"{improvement:.10g} loss points ({summary['improvement_percent_vs_incumbent']:.4g}%)"
    )
    readme = f"""# Intergenerational Calibration Refinement Report

The selected candidate has reported loss **{reported_loss:.10g}**, recomputed
15-moment loss **{recomputed_loss:.10g}**, market residual
**{format_cell(best.get('market_residual'))}**, status `{best.get('status')}`, and
strict convergence `{best.get('strict_converged')}`.

- Selected pool: `{selected['pool']}`
- Selected source: `{selected['source']}`{f":{selected['line']}" if selected['line'] else ""}
- Improvement versus incumbent: {improvement_text}
- Refinement files scanned: {counts['files_scanned']}
- Parsed refinement records: {counts['records_parsed']}
- Eligible refinement records: {counts['eligible_records']}
- Unique eligible refinement candidates: {counts['unique_eligible_candidates']}
- Parse errors: {counts['parse_errors']}

## Artifacts

- `best.json`: the complete selected candidate record.
- `calibration_target_fit.csv` and `.md`: every active target, model moment, gap,
  weight, and loss contribution.
- `calibration_parameters.csv` and `.md`: all 14 searched parameters, active
  bounds, normalized bound position, and the two-percent near-bound flag.
- `calibration_summary.json`: machine-readable provenance, counts, losses, and
  incumbent comparison.

`beta_annual` is reported as the annual discount factor, obtained from the
model-period value as $\\beta^{1/4}$. The fourteenth search coordinate is
$H_0 \\in [1,20]$; the other thirteen bounds come from
`PRODUCTION_SEARCH_BOUNDS`.
"""
    (args.outdir / "README.md").write_text(readme)


if __name__ == "__main__":
    main()
