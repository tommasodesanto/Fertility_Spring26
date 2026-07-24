#!/usr/bin/env python3
"""Reproduce the two established draft figures at the current M estimate."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import pickle
import time
from pathlib import Path
from typing import Any

import numpy as np

from intergen_housing_fertility_optimized.calibration import diagnostic_loss, extract_moments
from intergen_housing_fertility_optimized.new_moment_profile import (
    NEW_MOMENT_PROFILE_NAME,
    new_moment_overrides,
    new_moment_target_system,
)
from intergen_housing_fertility_optimized.solver import run_model_cp_dt
from plot_intergen_draft_figures_from_cache import (
    acs_age_profiles,
    decision_rule_rows,
    model_age_profiles,
    plot_decision_rules,
    plot_lifecycle,
    write_rows,
)


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SOURCE = (
    ROOT
    / "output/model/intergen_new_moment_unrestricted_overnight_20260723_w2/report/results.json"
)
DEFAULT_OUTDIR = ROOT / "output/model/intergen_current_m_figures_policy_20260724/equilibrium"
ACS_AGE = ROOT / "code/data/mms_center_periphery/output/mms_age_profiles_full.csv"
ACS_OWNERSHIP = (
    ROOT
    / "code/data/mms_center_periphery/output_ownership_audit/acs_ownership_age_profiles.csv"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_selected_source(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("profile") != NEW_MOMENT_PROFILE_NAME:
        raise ValueError(
            f"Expected profile {NEW_MOMENT_PROFILE_NAME!r}, found {payload.get('profile')!r}."
        )
    selected = payload.get("selected")
    if not isinstance(selected, dict) or selected.get("status") != "ok":
        raise ValueError(f"Source does not contain a successful selected result: {path}")
    theta = selected.get("theta")
    if not isinstance(theta, dict) or not theta:
        raise ValueError(f"Selected result does not contain theta: {path}")
    return payload, selected


def main() -> None:
    args = parse_args()
    source_path = args.source.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    payload, selected = load_selected_source(source_path)
    theta = {str(key): float(value) for key, value in selected["theta"].items()}
    overrides = new_moment_overrides(tight=True, optimized=True)
    overrides.update(theta)

    started = time.perf_counter()
    solution, parameters, price = run_model_cp_dt(overrides, verbose=False)
    elapsed = time.perf_counter() - started
    moments = extract_moments(solution, parameters)
    target_system = new_moment_target_system()
    loss = diagnostic_loss(
        moments,
        targets=target_system.targets_dict(),
        weights=target_system.weights_dict(),
    )

    target_rows: list[dict[str, Any]] = []
    for name, target, weight in zip(
        target_system.moment_names,
        target_system.target_values,
        target_system.weights,
    ):
        model = float(moments[name])
        gap = model - float(target)
        target_rows.append(
            {
                "moment": name,
                "target": float(target),
                "model": model,
                "gap": gap,
                "weight": float(weight),
                "loss_contribution": float(weight) * gap * gap,
            }
        )
    write_csv(outdir / "target_fit.csv", target_rows)
    write_csv(outdir / "parameters.csv", list(selected.get("parameters", [])))

    classic_dir = outdir / "classic_draft"
    lifecycle_png = classic_dir / "quant_lifecycle_equilibrium_repaired_nb120.png"
    decision_png = classic_dir / "quant_decision_rules_repaired_nb120.png"
    model_profiles = model_age_profiles(solution, parameters)
    data_profiles = acs_age_profiles(ACS_AGE, ACS_OWNERSHIP)
    decision_rows = decision_rule_rows(solution, parameters, requested_age=42.0)
    write_rows(
        lifecycle_png.with_suffix(".csv"),
        [{"source": "model", **row} for row in model_profiles]
        + [{"source": "ACS_2012_2023", **row} for row in data_profiles],
    )
    write_rows(decision_png.with_suffix(".csv"), decision_rows)
    plot_lifecycle(model_profiles, data_profiles, lifecycle_png)
    plot_decision_rules(decision_rows, decision_png, wealth_min=-2.2, wealth_max=7.2)

    cache_payload = {
        "status": "local_diagnostic_pickle_only_load_if_trusted",
        "created_at": dt.datetime.now().isoformat(timespec="seconds"),
        "source": str(source_path),
        "profile": NEW_MOMENT_PROFILE_NAME,
        "theta": theta,
        "baseline": {
            "label": "current_M_calibrated_baseline",
            "sol": solution,
            "P": parameters,
            "p_eq": np.asarray(price, dtype=float).reshape(-1),
            "moments": moments,
            "rank_loss": float(loss),
            "market_residual": float(solution.best_max_abs_rel_excess),
            "elapsed_sec": float(elapsed),
            "overrides": overrides,
        },
    }
    with (outdir / "solution_cache.pkl").open("wb") as handle:
        pickle.dump(cache_payload, handle, protocol=pickle.HIGHEST_PROTOCOL)

    source_price = float(selected["price"])
    source_loss = float(selected["rank_loss"])
    metadata = {
        "status": "complete",
        "purpose": "current-M reproduction of the two established draft figures",
        "source": str(source_path),
        "profile": NEW_MOMENT_PROFILE_NAME,
        "theta": theta,
        "elapsed_seconds": elapsed,
        "price": float(np.asarray(price).reshape(-1)[0]),
        "source_price": source_price,
        "price_difference": float(np.asarray(price).reshape(-1)[0]) - source_price,
        "loss": float(loss),
        "source_loss": source_loss,
        "loss_difference": float(loss) - source_loss,
        "market_residual": float(solution.best_max_abs_rel_excess),
        "strict_converged": bool(
            getattr(solution, "timings", {}).get("strict_converged", False)
        ),
        "calibration_baseline_fiscal_convention": (
            "the estimated current-M baseline, with the model's 1 percent property tax "
            "and no lump-sum rebate"
        ),
        "figures": [str(lifecycle_png), str(decision_png)],
        "selected_result_status": payload.get("status"),
    }
    (outdir / "metadata.json").write_text(
        json.dumps(jsonable(metadata), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        f"reproduced current M: loss={loss:.9f}, price={metadata['price']:.9f}, "
        f"market residual={metadata['market_residual']:.3e}, elapsed={elapsed:.1f}s",
        flush=True,
    )
    print(f"wrote {lifecycle_png}", flush=True)
    print(f"wrote {decision_png}", flush=True)


if __name__ == "__main__":
    main()
