#!/usr/bin/env python3
"""Targeted mechanism grid audit for the one-market intergen model.

This is a diagnostic, not an SMM calibration. It takes saved frontier points
and tests whether local combinations of access, finance, rental-cap, and
owner/renter service-wedge levers can jointly improve young ownership, old-age
exit, and owner-renter room separation.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from intergen_housing_fertility.calibration import (  # noqa: E402
    base_overrides,
    diagnostic_loss,
    extract_moments,
    get_target_set,
    jsonable,
)
from intergen_housing_fertility.local_panel import (  # noqa: E402
    GLOBAL_DE_BOUNDS,
    income_process_overrides,
)
from intergen_housing_fertility.solver import run_model_cp_dt  # noqa: E402


DEFAULT_POINT_DIR = ROOT / "code/cluster/audit_points_20260621_frontier"
DEFAULT_OUTDIR = ROOT / "output/model/intergen_mechanism_grid_20260622"

INTERNAL_PARAMETERS = [
    "beta",
    "alpha_cons",
    "b_entry_fixed",
    "c_bar_0",
    "c_bar_n",
    "h_bar_0",
    "h_bar_jump",
    "h_bar_n",
    "psi_child",
    "kappa_fert",
    "chi",
    "theta0",
    "theta_n",
]

KEY_MOMENTS = [
    "tfr",
    "childless_rate",
    "own_rate",
    "own_rate_2534",
    "old_age_own_rate",
    "prime30_55_childless_renter_mean_rooms",
    "prime30_55_childless_owner_mean_rooms",
    "prime30_55_childless_owner_minus_renter_mean_rooms",
    "old_nonhousing_wealth_to_income_median_6575",
    "old_parent_childless_nonhousing_wealth_to_income_gap_6575",
    "young_liquid_wealth_to_income",
    "housing_increment_0to1",
    "housing_increment_1to2",
    "market_residual",
]

STRICT_SCREEN = {
    "young_own_min": 0.25,
    "old_own_max": 0.85,
    "room_gap_min": 1.50,
    "tfr_min": 1.55,
    "tfr_max": 1.90,
    "childless_max": 0.30,
}

SOFT_SCREEN = {
    "young_own_min": 0.20,
    "old_own_max": 0.88,
    "room_gap_min": 1.00,
    "tfr_min": 1.55,
    "tfr_max": 1.90,
    "childless_max": 0.30,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--point-json", action="append", default=[], help="Saved point JSON, optionally LABEL=PATH")
    parser.add_argument("--point-label", action="append", default=[], help="Point label under --point-dir")
    parser.add_argument("--point-dir", type=Path, default=DEFAULT_POINT_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--target-set", default="", help="Override target set; default uses point source_target_set")
    parser.add_argument("--J", type=int, default=16)
    parser.add_argument("--Nb", type=int, default=60)
    parser.add_argument("--n-house", type=int, default=5)
    parser.add_argument("--max-iter-eq", type=int, default=3)
    parser.add_argument("--income-states", type=int, default=5)
    parser.add_argument("--max-cases", type=int, default=0, help="Smoke-test cap; 0 runs all generated cases")
    args = parser.parse_args()

    os.environ.setdefault("NUMBA_NUM_THREADS", "1")
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

    points = load_points(args)
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_json(
        args.outdir / "audit_config.json",
        {
            "status": "mechanism_grid_diagnostic_not_calibration",
            "points": [p["point_label"] for p in points],
            "target_set_override": args.target_set,
            "J": args.J,
            "Nb": args.Nb,
            "n_house": args.n_house,
            "max_iter_eq": args.max_iter_eq,
            "income_states": args.income_states,
            "strict_screen": STRICT_SCREEN,
            "soft_screen": SOFT_SCREEN,
        },
    )

    summaries = []
    for point in points:
        summaries.append(run_point(point, args))
    write_json(args.outdir / "audit_summary.json", {"points": summaries})


def load_points(args: argparse.Namespace) -> list[dict[str, Any]]:
    specs: list[str] = []
    specs.extend(args.point_json)
    specs.extend(f"{label}={args.point_dir / (label + '.json')}" for label in args.point_label)
    if not specs:
        specs = [
            f"{path.stem}={path}"
            for path in sorted(args.point_dir.glob("*.json"))
        ]
    return [load_point_json(spec) for spec in specs]


def load_point_json(spec: str) -> dict[str, Any]:
    if "=" in spec:
        label, path_text = spec.split("=", 1)
        point_label = label.strip()
        path = Path(path_text).expanduser()
    else:
        path = Path(spec).expanduser()
        point_label = path.stem
    if not path.exists():
        raise FileNotFoundError(path)
    obj = json.loads(path.read_text())
    rec = dict(obj.get("best") or obj)
    theta = dict(rec.get("theta") or {})
    if not theta:
        raise ValueError(f"No theta in {path}")
    return {
        "point_label": point_label,
        "source_target_set": rec.get("source_target_set", ""),
        "source_rank_loss": float(rec.get("rank_loss", math.nan)),
        "source_path": str(path),
        "theta": {k: float(v) for k, v in theta.items() if k in INTERNAL_PARAMETERS},
        "source_moments": dict(rec.get("moments", {})),
    }


def run_point(point: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    label = str(point["point_label"])
    outdir = args.outdir / label
    outdir.mkdir(parents=True, exist_ok=True)
    cases_path = outdir / "cases.jsonl"
    cases_path.write_text("")

    target_set = str(args.target_set or point.get("source_target_set") or "candidate_replacement_old_retention_v1")
    targets, weights = get_target_set(target_set)
    income = income_process_overrides(args.income_states)
    base = {
        **base_overrides(J=args.J, Nb=args.Nb, n_house=args.n_house, max_iter_eq=args.max_iter_eq),
        **income,
    }

    candidates = generate_candidates(point["theta"])
    if args.max_cases > 0:
        candidates = candidates[: int(args.max_cases)]

    metadata = {
        "point_label": label,
        "source_target_set": point.get("source_target_set"),
        "run_target_set": target_set,
        "source_rank_loss": point.get("source_rank_loss"),
        "source_path": point.get("source_path"),
        "n_cases": len(candidates),
        "J": args.J,
        "Nb": args.Nb,
        "n_house": args.n_house,
        "max_iter_eq": args.max_iter_eq,
        "income_states": args.income_states,
        "strict_screen": STRICT_SCREEN,
        "soft_screen": SOFT_SCREEN,
        "targets": targets,
        "weights": weights,
    }
    write_json(outdir / "metadata.json", metadata)

    best_rank: dict[str, Any] | None = None
    best_joint: dict[str, Any] | None = None
    best_soft: dict[str, Any] | None = None
    best_strict: dict[str, Any] | None = None
    rows: list[dict[str, Any]] = []
    start = time.perf_counter()

    for idx, candidate in enumerate(candidates):
        record = evaluate_case(idx, candidate, base, targets, weights)
        record["point_label"] = label
        append_jsonl(cases_path, record)
        write_json(outdir / "latest_case.json", record)

        if is_better(record, best_rank, "rank_loss"):
            best_rank = record
            write_json(outdir / "best_by_rank_loss.json", record)
        if is_better(record, best_joint, "joint_distance"):
            best_joint = record
            write_json(outdir / "best_by_joint_distance.json", record)
        if record.get("soft_pass") and is_better(record, best_soft, "joint_distance"):
            best_soft = record
            write_json(outdir / "best_soft_pass.json", record)
        if record.get("strict_pass") and is_better(record, best_strict, "rank_loss"):
            best_strict = record
            write_json(outdir / "best_strict_pass.json", record)

        rows.append(case_csv_row(record))
        write_csv(outdir / "cases_summary.csv", rows)
        print(
            f"{label} {idx + 1}/{len(candidates)} {record['case_label']}: "
            f"rank={record['rank_loss']:.3g} joint={record['joint_distance']:.3g} "
            f"young={record['moments_key'].get('own_rate_2534', math.nan):.3f} "
            f"old={record['moments_key'].get('old_age_own_rate', math.nan):.3f} "
            f"gap={record['moments_key'].get('prime30_55_childless_owner_minus_renter_mean_rooms', math.nan):.3f}",
            flush=True,
        )

    summary = {
        "point_label": label,
        "target_set": target_set,
        "elapsed_sec": time.perf_counter() - start,
        "n_cases": len(candidates),
        "strict_count": sum(1 for r in rows if truthy(r["strict_pass"])),
        "soft_count": sum(1 for r in rows if truthy(r["soft_pass"])),
        "best_rank": compact_record(best_rank),
        "best_joint": compact_record(best_joint),
        "best_soft": compact_record(best_soft),
        "best_strict": compact_record(best_strict),
    }
    write_json(outdir / "point_summary.json", summary)
    return summary


def generate_candidates(base_theta: dict[str, float]) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
    seen: set[str] = set()

    def add(block: str, label: str, overrides: dict[str, Any]) -> None:
        theta = bounded_theta({**base_theta, **overrides})
        key = json.dumps(jsonable(theta), sort_keys=True)
        if key in seen:
            return
        seen.add(key)
        candidates.append({"block": block, "label": label, "theta": theta})

    add("baseline", "baseline", {})

    for mult in [0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30]:
        add("chi_only", f"chi_x{mult:.2f}", {"chi": base_theta["chi"] * mult})

    for delta in [1.0, 2.0, 4.0, 6.0]:
        add("entry_wealth_only", f"b_entry_plus{delta:g}", {"b_entry_fixed": base_theta["b_entry_fixed"] + delta})
    for phi in [0.85, 0.90, 0.95, 0.98]:
        add("finance_phi_only", f"phi_{phi:.2f}", {"phi": phi})
    for pti in [0.45, 0.60]:
        add("finance_pti_only", f"pti_{pti:.2f}", {"pti_limit": pti})

    for chi_mult in [0.75, 0.85, 0.95, 1.00]:
        for delta in [0.0, 2.0, 4.0, 6.0]:
            for phi in [0.80, 0.90, 0.95]:
                for pti in [0.30, 0.45]:
                    add(
                        "chi_access_compensation",
                        f"chi{chi_mult:.2f}_b{delta:g}_phi{phi:.2f}_pti{pti:.2f}",
                        {
                            "chi": base_theta["chi"] * chi_mult,
                            "b_entry_fixed": base_theta["b_entry_fixed"] + delta,
                            "phi": phi,
                            "pti_limit": pti,
                        },
                    )

    for hR in [4.0, 4.5, 5.0, 5.5, 6.0]:
        for delta in [0.0, 2.0, 4.0]:
            for phi in [0.80, 0.90, 0.95]:
                add(
                    "rental_cap_access",
                    f"hR{hR:.1f}_b{delta:g}_phi{phi:.2f}",
                    {
                        "hR_max": hR,
                        "b_entry_fixed": base_theta["b_entry_fixed"] + delta,
                        "phi": phi,
                    },
                )

    for hR in [4.0, 4.5, 5.0]:
        for chi_mult in [0.85, 1.00]:
            for delta in [0.0, 2.0, 4.0]:
                add(
                    "rental_cap_chi_access",
                    f"hR{hR:.1f}_chi{chi_mult:.2f}_b{delta:g}",
                    {
                        "hR_max": hR,
                        "chi": base_theta["chi"] * chi_mult,
                        "b_entry_fixed": base_theta["b_entry_fixed"] + delta,
                    },
                )

    return candidates


def bounded_theta(theta: dict[str, Any]) -> dict[str, Any]:
    out = dict(theta)
    bounds = parameter_bounds()
    for key, (lo, hi) in bounds.items():
        if key in out:
            out[key] = min(max(float(out[key]), lo), hi)
    if "phi" in out:
        out["phi"] = min(max(float(out["phi"]), 0.50), 0.995)
    if "pti_limit" in out:
        out["pti_limit"] = min(max(float(out["pti_limit"]), 0.05), 1.00)
    if "hR_max" in out:
        out["hR_max"] = min(max(float(out["hR_max"]), 2.0), 10.0)
    return out


def parameter_bounds() -> dict[str, tuple[float, float]]:
    bounds: dict[str, tuple[float, float]] = {}
    for name, lo, hi in GLOBAL_DE_BOUNDS:
        if name == "beta_annual":
            bounds["beta"] = (float(lo) ** 4.0, float(hi) ** 4.0)
        else:
            bounds[name] = (float(lo), float(hi))
    return bounds


def evaluate_case(
    idx: int,
    candidate: dict[str, Any],
    base: dict[str, Any],
    targets: dict[str, float],
    weights: dict[str, float],
) -> dict[str, Any]:
    theta = dict(candidate["theta"])
    overrides = {**base, **theta}
    t0 = time.perf_counter()
    try:
        sol, P, p_eq = run_model_cp_dt(overrides, verbose=False)
        moments = extract_moments(sol, P)
        loss = diagnostic_loss(moments, targets=targets, weights=weights)
        status = "ok"
        market_residual = float(getattr(sol, "best_max_abs_rel_excess", np.nan))
        timings = getattr(sol, "timings", {})
    except Exception as exc:  # noqa: BLE001 - diagnostic should preserve failures.
        moments = {}
        loss = math.inf
        status = f"failed: {type(exc).__name__}: {exc}"
        p_eq = np.array([np.nan])
        market_residual = math.inf
        timings = {}
    moments_key = selected_moments(moments)
    moments_key["market_residual"] = market_residual
    joint_distance, strict_pass, soft_pass = screen_metrics(moments_key)
    return {
        "case": int(idx),
        "case_label": str(candidate["label"]),
        "block": str(candidate["block"]),
        "status": status,
        "rank_loss": float(loss),
        "joint_distance": float(joint_distance),
        "strict_pass": bool(strict_pass),
        "soft_pass": bool(soft_pass),
        "theta": jsonable(theta),
        "moments": jsonable(moments),
        "moments_key": jsonable(moments_key),
        "p_eq": jsonable(p_eq),
        "market_residual": float(market_residual),
        "elapsed_sec": float(time.perf_counter() - t0),
        "timings": jsonable(timings),
    }


def selected_moments(moments: dict[str, Any]) -> dict[str, float]:
    out = {key: finite_float(moments.get(key)) for key in KEY_MOMENTS}
    if not math.isfinite(out["prime30_55_childless_owner_minus_renter_mean_rooms"]):
        owner = out["prime30_55_childless_owner_mean_rooms"]
        renter = out["prime30_55_childless_renter_mean_rooms"]
        out["prime30_55_childless_owner_minus_renter_mean_rooms"] = owner - renter
    return out


def screen_metrics(m: dict[str, float]) -> tuple[float, bool, bool]:
    strict = screen_pass(m, STRICT_SCREEN)
    soft = screen_pass(m, SOFT_SCREEN)
    v = [
        max(0.0, STRICT_SCREEN["young_own_min"] - m["own_rate_2534"]) / STRICT_SCREEN["young_own_min"],
        max(0.0, m["old_age_own_rate"] - STRICT_SCREEN["old_own_max"]) / 0.15,
        max(0.0, STRICT_SCREEN["room_gap_min"] - m["prime30_55_childless_owner_minus_renter_mean_rooms"])
        / STRICT_SCREEN["room_gap_min"],
        max(0.0, STRICT_SCREEN["tfr_min"] - m["tfr"]) / 0.35,
        max(0.0, m["tfr"] - STRICT_SCREEN["tfr_max"]) / 0.35,
        max(0.0, m["childless_rate"] - STRICT_SCREEN["childless_max"]) / STRICT_SCREEN["childless_max"],
    ]
    if any(not math.isfinite(x) for x in v):
        return math.inf, strict, soft
    return float(sum(x * x for x in v)), strict, soft


def screen_pass(m: dict[str, float], screen: dict[str, float]) -> bool:
    return (
        m["own_rate_2534"] >= screen["young_own_min"]
        and m["old_age_own_rate"] <= screen["old_own_max"]
        and m["prime30_55_childless_owner_minus_renter_mean_rooms"] >= screen["room_gap_min"]
        and screen["tfr_min"] <= m["tfr"] <= screen["tfr_max"]
        and m["childless_rate"] <= screen["childless_max"]
    )


def is_better(record: dict[str, Any], best: dict[str, Any] | None, key: str) -> bool:
    val = float(record.get(key, math.inf))
    if not math.isfinite(val):
        return False
    if best is None:
        return True
    return val < float(best.get(key, math.inf))


def case_csv_row(record: dict[str, Any]) -> dict[str, Any]:
    m = record["moments_key"]
    row = {
        "case": record["case"],
        "case_label": record["case_label"],
        "block": record["block"],
        "status": record["status"],
        "rank_loss": record["rank_loss"],
        "joint_distance": record["joint_distance"],
        "strict_pass": record["strict_pass"],
        "soft_pass": record["soft_pass"],
        "market_residual": record["market_residual"],
        "elapsed_sec": record["elapsed_sec"],
    }
    row.update({key: m.get(key, math.nan) for key in KEY_MOMENTS})
    theta = record["theta"]
    for key in ["chi", "b_entry_fixed", "phi", "pti_limit", "hR_max", "alpha_cons", "c_bar_n", "h_bar_0"]:
        row[f"theta_{key}"] = theta.get(key, math.nan)
    return row


def compact_record(record: dict[str, Any] | None) -> dict[str, Any] | None:
    if record is None:
        return None
    return {
        "case": record["case"],
        "case_label": record["case_label"],
        "block": record["block"],
        "rank_loss": record["rank_loss"],
        "joint_distance": record["joint_distance"],
        "strict_pass": record["strict_pass"],
        "soft_pass": record["soft_pass"],
        "theta": record["theta"],
        "moments_key": record["moments_key"],
    }


def finite_float(value: Any) -> float:
    try:
        out = float(value)
    except Exception:
        return math.nan
    return out if math.isfinite(out) else math.nan


def truthy(value: Any) -> bool:
    if isinstance(value, str):
        return value.lower() == "true"
    return bool(value)


def append_jsonl(path: Path, record: dict[str, Any]) -> None:
    with path.open("a") as fh:
        fh.write(json.dumps(jsonable(record), sort_keys=True) + "\n")


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(jsonable(obj), indent=2, sort_keys=True) + "\n")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
