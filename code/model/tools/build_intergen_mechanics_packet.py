#!/usr/bin/env python3
"""Build a diagnostic mechanics packet for the one-market intergen model.

This script is intentionally a non-production inspection tool. It reads a
saved candidate theta, re-solves the active one-market intergenerational model,
and writes plots/tables for policy-function and target-fit inspection. It does
not change model logic or calibration targets.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import json
import math
import pickle
import shutil
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
from intergen_housing_fertility.diagnostics import write_diagnostics  # noqa: E402
from intergen_housing_fertility.local_panel import income_process_overrides  # noqa: E402
from intergen_housing_fertility.solver import run_model_cp_dt  # noqa: E402
from intergen_housing_fertility.solver import precompute_shared  # noqa: E402


DEFAULT_SOURCE = ROOT / "output/model/intergen_room_distribution_current_best_20260623/summary.json"
DEFAULT_TARGET_SET = "candidate_replacement_young_old_roomgap_v1"
SOLUTION_CACHE_NAME = "solution_cache.pkl"

ROOM_BIN_DATA_SOURCE = (
    "ACS/IPUMS extract27 household heads, 2012-2023, matched MMS metro sample; "
    "prime age 30-55 childless heads; see "
    "output/model/intergen_room_distribution_current_best_20260623/"
    "room_bin_shares_prime30_55_childless.csv"
)
ROOM_BIN_DATA = {
    ("renter", "S_<=4"): 0.719075271037314,
    ("renter", "M_4to6"): 0.22056262451180755,
    ("renter", "L_>6"): 0.0603621044508785,
    ("owner", "S_<=4"): 0.20064321188507736,
    ("owner", "M_4to6"): 0.41890838244888606,
    ("owner", "L_>6"): 0.3804484056660366,
}
ROOM_BIN_ORDER = ["S_<=4", "M_4to6", "L_>6"]
TENURE_ORDER = ["renter", "owner"]


def main() -> None:
    args = parse_args()
    outdir = args.outdir or default_outdir()
    if args.clean_outdir and outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cache_path = solution_cache_path(args, outdir)

    targets, weights = get_target_set(args.target_set)
    policy_records: list[dict[str, Any]] = []

    if args.refresh_plots_from_cache:
        cache_payload, baseline, source, grid = load_solution_cache(cache_path)
        if cache_payload.get("target_set") == args.target_set:
            targets = dict(cache_payload.get("targets", targets))
            weights = dict(cache_payload.get("weights", weights))
        write_packet_metadata(
            outdir,
            source=source,
            target_set=args.target_set,
            targets=targets,
            weights=weights,
            grid=grid,
            args=args,
            cache_path=cache_path,
            refreshed_from_cache=True,
        )
        if not args.skip_standard_diagnostics:
            write_diagnostics(baseline["sol"], baseline["P"], outdir / "diagnostics")
        write_core_outputs(
            outdir,
            baseline,
            source,
            targets,
            weights,
            include_contact_sheet=not args.skip_contact_sheet,
            quick_first_look_only=bool(args.quick_first_look_only),
            write_csv_sidecars=not args.no_csv,
        )
        write_readme(outdir, source, baseline, args.target_set, targets, policy_records)
        print(
            f"Refreshed intergen mechanics packet from {cache_path} "
            f"(loss={baseline['rank_loss']:.6g}, residual={baseline['market_residual']:.2e})",
            flush=True,
        )
        return

    source = load_source_record(args.source, row=int(args.source_row))
    grid = resolve_grid(args, source)
    write_packet_metadata(
        outdir,
        source=source,
        target_set=args.target_set,
        targets=targets,
        weights=weights,
        grid=grid,
        args=args,
        cache_path=cache_path,
        refreshed_from_cache=False,
    )

    baseline = solve_candidate(
        theta=source["theta"],
        grid=grid,
        extra_overrides={},
        targets=targets,
        weights=weights,
        label="baseline",
    )
    if not args.no_save_solution_cache:
        write_solution_cache(
            cache_path,
            baseline=baseline,
            source=source,
            target_set=args.target_set,
            targets=targets,
            weights=weights,
            grid=grid,
        )
    baseline_dir = outdir / "baseline"
    baseline_dir.mkdir(parents=True, exist_ok=True)

    if not args.skip_standard_diagnostics:
        write_diagnostics(baseline["sol"], baseline["P"], outdir / "diagnostics")

    write_core_outputs(
        outdir,
        baseline,
        source,
        targets,
        weights,
        include_contact_sheet=not args.skip_contact_sheet,
        quick_first_look_only=bool(args.quick_first_look_only),
        write_csv_sidecars=not args.no_csv,
    )
    if args.run_policy_cases:
        policy_records = write_policy_cases(
            outdir / "policy_cases",
            source=source,
            grid=grid,
            baseline=baseline,
            targets=targets,
            weights=weights,
            write_case_diagnostics=bool(args.policy_diagnostics),
        )

    write_readme(outdir, source, baseline, args.target_set, targets, policy_records)
    print(
        f"Wrote intergen mechanics packet to {outdir} "
        f"(loss={baseline['rank_loss']:.6g}, residual={baseline['market_residual']:.2e})",
        flush=True,
    )


def write_packet_metadata(
    outdir: Path,
    *,
    source: dict[str, Any],
    target_set: str,
    targets: dict[str, float],
    weights: dict[str, float],
    grid: dict[str, Any],
    args: argparse.Namespace,
    cache_path: Path,
    refreshed_from_cache: bool,
) -> None:
    write_json(
        outdir / "metadata.json",
        {
            "status": "diagnostic_only_not_production_policy_or_calibration",
            "source": source["source_meta"],
            "target_set": target_set,
            "targets": targets,
            "weights": weights,
            "grid": grid,
            "policy_cases_requested": bool(args.run_policy_cases),
            "solution_cache": str(cache_path),
            "refreshed_from_solution_cache": bool(refreshed_from_cache),
            "room_bin_data_source": ROOM_BIN_DATA_SOURCE,
        },
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--source-row", type=int, default=1, help="1-based row for CSV sources.")
    parser.add_argument("--outdir", type=Path, default=None)
    parser.add_argument("--target-set", default=DEFAULT_TARGET_SET)
    parser.add_argument("--J", type=int, default=None)
    parser.add_argument("--Nb", type=int, default=None)
    parser.add_argument("--income-states", type=int, default=None)
    parser.add_argument("--n-house", type=int, default=None)
    parser.add_argument("--hR-max", type=float, default=None)
    parser.add_argument("--H-own", default=None, help="Comma-separated owner ladder, e.g. 2,4,6,8,10.")
    parser.add_argument("--max-iter-eq", type=int, default=10)
    # Grid geometry: dense-core + sparse-buffer controls (None => model defaults).
    parser.add_argument("--b-min", type=float, default=None, help="Liquid-wealth grid lower bound (default -35).")
    parser.add_argument("--b-max", type=float, default=None, help="Liquid-wealth grid upper bound (default 100).")
    parser.add_argument("--b-core-lo", type=float, default=None, help="Dense-core lower edge (default -3).")
    parser.add_argument("--b-core-hi", type=float, default=None, help="Dense-core upper edge (default 6).")
    parser.add_argument("--b-mid-hi", type=float, default=None, help="Mid-band upper edge (default 20).")
    parser.add_argument("--b-frac-low", type=float, default=None, help="Node fraction in the lower tail (default 0.15).")
    parser.add_argument("--b-frac-core", type=float, default=None, help="Node fraction in the dense core (default 0.45).")
    parser.add_argument("--b-frac-mid", type=float, default=None, help="Node fraction in the mid band (default 0.15).")
    parser.add_argument("--interp-method", default=None, choices=["linear", "monotone_cubic"],
                        help="Value-function interpolation for the savings continuation (default linear).")
    parser.add_argument(
        "--entry-wealth-spread-nodes",
        type=int,
        default=None,
        help="Diagnostic entrant-wealth distribution width. Default leaves the model's point-mass entry closure.",
    )
    parser.add_argument("--skip-standard-diagnostics", action="store_true")
    parser.add_argument("--skip-contact-sheet", action="store_true")
    parser.add_argument("--quick-first-look-only", action="store_true")
    parser.add_argument(
        "--no-csv",
        action="store_true",
        help="Skip CSV sidecars for plot source data; keep JSON/Markdown summaries and figures.",
    )
    parser.add_argument(
        "--clean-outdir",
        action="store_true",
        help="Delete OUTDIR before writing, useful for the stable current quick-review folder.",
    )
    parser.add_argument(
        "--refresh-plots-from-cache",
        action="store_true",
        help=f"Skip the model solve and rebuild plots/tables from {SOLUTION_CACHE_NAME}.",
    )
    parser.add_argument(
        "--solution-cache",
        type=Path,
        default=None,
        help=f"Path for the local solved-object cache. Defaults to OUTDIR/{SOLUTION_CACHE_NAME}.",
    )
    parser.add_argument(
        "--no-save-solution-cache",
        action="store_true",
        help="Do not write the local solved-object pickle after solving.",
    )
    parser.add_argument("--run-policy-cases", action="store_true")
    parser.add_argument("--policy-diagnostics", action="store_true")
    return parser.parse_args()


def default_outdir() -> Path:
    return ROOT / "output/model/intergen_mechanics_packet_current"


def solution_cache_path(args: argparse.Namespace, outdir: Path) -> Path:
    return args.solution_cache if args.solution_cache is not None else outdir / SOLUTION_CACHE_NAME


def load_source_record(path: Path, *, row: int) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"source record not found: {path}")
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith(".json"):
        payload = json.loads(path.read_text())
        theta = extract_theta_from_json(payload)
        grid = dict(payload.get("grid", {})) if isinstance(payload, dict) else {}
        return {
            "theta": theta,
            "grid": grid,
            "raw": payload,
            "source_meta": {
                "path": str(path),
                "format": "json",
                "source_record": payload.get("source_record") if isinstance(payload, dict) else None,
            },
        }
    if suffixes.endswith(".csv") or suffixes.endswith(".csv.gz"):
        rows = read_csv_rows(path)
        if row < 1 or row > len(rows):
            raise IndexError(f"CSV source row {row} outside 1..{len(rows)} for {path}")
        record = rows[row - 1]
        if "theta" not in record:
            raise KeyError(f"CSV source {path} has no theta column")
        theta = json.loads(record["theta"])
        return {
            "theta": coerce_theta(theta),
            "grid": {},
            "raw": record,
            "source_meta": {
                "path": str(path),
                "format": "csv",
                "row": int(row),
                "label": record.get("label"),
                "rank_loss": maybe_float(record.get("rank_loss")),
            },
        }
    raise ValueError(f"unsupported source file type: {path}")


def extract_theta_from_json(payload: Any) -> dict[str, float]:
    if not isinstance(payload, dict):
        raise TypeError("JSON source must be an object")
    if isinstance(payload.get("theta"), dict):
        return coerce_theta(payload["theta"])
    if isinstance(payload.get("best"), dict) and isinstance(payload["best"].get("theta"), dict):
        return coerce_theta(payload["best"]["theta"])
    if isinstance(payload.get("source_record"), dict) and isinstance(payload["source_record"].get("theta"), dict):
        return coerce_theta(payload["source_record"]["theta"])
    raise KeyError("could not find theta in JSON source")


def coerce_theta(theta: dict[str, Any]) -> dict[str, float]:
    return {str(k): float(v) for k, v in theta.items()}


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if "".join(path.suffixes).lower().endswith(".csv.gz") else open
    with opener(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh))


def resolve_grid(args: argparse.Namespace, source: dict[str, Any]) -> dict[str, Any]:
    source_grid = dict(source.get("grid", {}))
    h_own = parse_float_list(args.H_own) if args.H_own else source_grid.get("H_own")
    if h_own is not None:
        h_own = [float(x) for x in h_own]
    n_house = args.n_house
    if n_house is None:
        n_house = len(h_own) if h_own is not None else int(source_grid.get("n_house", 5))
    return {
        "J": int(args.J if args.J is not None else source_grid.get("J", 16)),
        "Nb": int(args.Nb if args.Nb is not None else source_grid.get("Nb", 60)),
        "income_states": int(args.income_states if args.income_states is not None else source_grid.get("Nz", 5)),
        "n_house": int(n_house),
        "H_own": h_own,
        "hR_max": float(args.hR_max if args.hR_max is not None else source_grid.get("hR_max", 8.0)),
        "max_iter_eq": int(args.max_iter_eq),
        # Optional grid-geometry / interpolation overrides (None => use defaults).
        "b_min": args.b_min,
        "b_max": args.b_max,
        "b_core_lo": args.b_core_lo,
        "b_core_hi": args.b_core_hi,
        "b_mid_hi": args.b_mid_hi,
        "b_frac_low": args.b_frac_low,
        "b_frac_core": args.b_frac_core,
        "b_frac_mid": args.b_frac_mid,
        "interp_method": args.interp_method,
        "entry_wealth_spread_nodes": args.entry_wealth_spread_nodes,
    }


def parse_float_list(text: str) -> list[float]:
    vals = [float(x.strip()) for x in text.split(",") if x.strip()]
    if not vals:
        raise ValueError("--H-own must contain at least one number")
    return vals


def solve_candidate(
    *,
    theta: dict[str, float],
    grid: dict[str, Any],
    extra_overrides: dict[str, Any],
    targets: dict[str, float],
    weights: dict[str, float],
    label: str,
) -> dict[str, Any]:
    overrides = {
        **base_overrides(
            J=int(grid["J"]),
            Nb=int(grid["Nb"]),
            n_house=int(grid["n_house"]),
            max_iter_eq=int(grid["max_iter_eq"]),
        ),
        **income_process_overrides(int(grid["income_states"])),
        **theta,
        "hR_max": float(grid["hR_max"]),
        "diagnostic_policy_ages": np.array([30.0, 42.0]),
        **extra_overrides,
    }
    if grid.get("H_own") is not None:
        overrides["H_own"] = np.asarray(grid["H_own"], dtype=float)
        overrides["n_house"] = len(overrides["H_own"])
    for _gk in ("b_min", "b_max", "b_core_lo", "b_core_hi", "b_mid_hi",
                "b_frac_low", "b_frac_core", "b_frac_mid", "interp_method",
                "entry_wealth_spread_nodes"):
        if grid.get(_gk) is not None:
            overrides[_gk] = grid[_gk]
    t0 = time.perf_counter()
    sol, P, p_eq = run_model_cp_dt(overrides, verbose=False)
    elapsed = time.perf_counter() - t0
    moments = extract_moments(sol, P)
    rank_loss = diagnostic_loss(moments, targets=targets, weights=weights)
    return {
        "label": label,
        "sol": sol,
        "P": P,
        "p_eq": np.asarray(p_eq, dtype=float).reshape(-1),
        "moments": moments,
        "rank_loss": float(rank_loss),
        "market_residual": float(getattr(sol, "best_max_abs_rel_excess", np.nan)),
        "elapsed_sec": float(elapsed),
        "overrides": overrides,
    }


def write_solution_cache(
    path: Path,
    *,
    baseline: dict[str, Any],
    source: dict[str, Any],
    target_set: str,
    targets: dict[str, float],
    weights: dict[str, float],
    grid: dict[str, Any],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "status": "local_diagnostic_pickle_only_load_if_trusted",
        "created_at": dt.datetime.now().isoformat(timespec="seconds"),
        "source_meta": source["source_meta"],
        "theta": source["theta"],
        "target_set": target_set,
        "targets": dict(targets),
        "weights": dict(weights),
        "grid": dict(grid),
        "baseline": {
            "label": baseline["label"],
            "sol": baseline["sol"],
            "P": baseline["P"],
            "p_eq": np.asarray(baseline["p_eq"], dtype=float).reshape(-1),
            "moments": dict(baseline["moments"]),
            "rank_loss": float(baseline["rank_loss"]),
            "market_residual": float(baseline["market_residual"]),
            "elapsed_sec": float(baseline["elapsed_sec"]),
            "overrides": dict(baseline["overrides"]),
        },
    }
    with path.open("wb") as fh:
        pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)


def load_solution_cache(path: Path) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    if not path.exists():
        raise FileNotFoundError(f"solution cache not found: {path}")
    with path.open("rb") as fh:
        payload = pickle.load(fh)
    if not isinstance(payload, dict) or not isinstance(payload.get("baseline"), dict):
        raise ValueError(f"solution cache has unexpected format: {path}")
    baseline = dict(payload["baseline"])
    baseline["p_eq"] = np.asarray(baseline["p_eq"], dtype=float).reshape(-1)
    source = {
        "theta": dict(payload.get("theta", {})),
        "source_meta": {
            **dict(payload.get("source_meta", {})),
            "solution_cache": str(path),
            "format": "solution_cache_pickle",
        },
    }
    grid = dict(payload.get("grid", {}))
    return payload, baseline, source, grid


def write_core_outputs(
    outdir: Path,
    baseline: dict[str, Any],
    source: dict[str, Any],
    targets: dict[str, float],
    weights: dict[str, float],
    *,
    include_contact_sheet: bool = True,
    quick_first_look_only: bool = False,
    write_csv_sidecars: bool = True,
) -> None:
    sol = baseline["sol"]
    P = baseline["P"]
    moments = baseline["moments"]
    write_json(
        outdir / "solution_summary.json",
        {
            "status": "diagnostic_only_not_production_calibration",
            "source": source["source_meta"],
            "theta": source["theta"],
            "p_eq": baseline["p_eq"],
            "rank_loss": baseline["rank_loss"],
            "market_residual": baseline["market_residual"],
            "elapsed_sec": baseline["elapsed_sec"],
            "solver_timings": jsonable(getattr(sol, "timings", {})),
            "H_own": np.asarray(P.H_own, dtype=float),
            "hR_max": float(P.hR_max),
            "entry_wealth_distribution": {
                "spread_nodes": int(getattr(P, "entry_wealth_spread_nodes", 1)),
                "grid_values": jsonable(getattr(sol, "entry_wealth_grid_values", np.array([]))),
                "weights": jsonable(getattr(sol, "entry_wealth_grid_weights", np.array([]))),
                "mean": maybe_float(getattr(sol, "entry_wealth_grid_mean", np.nan)),
                "target_b_entry_fixed": float(getattr(P, "b_entry_fixed", np.nan)),
            },
            "moments": moments,
        },
    )
    write_json(outdir / "moments.json", moments)

    target_rows = target_fit_rows(moments, targets, weights)
    if write_csv_sidecars:
        write_csv(outdir / "target_fit.csv", target_rows)
    write_markdown_table(outdir / "target_fit.md", target_rows, "Target Fit")

    if not quick_first_look_only:
        room_rows, room_fit_rows = room_bin_outputs(sol, P)
        if write_csv_sidecars:
            write_csv(outdir / "room_bin_shares_prime30_55_childless.csv", room_rows)
            write_csv(outdir / "room_bin_fit_prime30_55_childless.csv", room_fit_rows)
        plot_room_bins(room_fit_rows, outdir / "room_bin_shares_prime30_55_childless.png")

        rung_rows = owner_rung_rows(sol, P)
        if write_csv_sidecars:
            write_csv(outdir / "owner_rung_shares_prime30_55_childless.csv", rung_rows)
        plot_owner_rungs(
            rung_rows,
            outdir / "owner_rung_shares_prime30_55_childless.png",
            ylabel="share of prime-age childless owners",
        )
        all_owner_rung_rows = owner_rung_rows_all(sol, P)
        if write_csv_sidecars:
            write_csv(outdir / "owner_rung_shares_all_owners.csv", all_owner_rung_rows)
        plot_owner_rungs(all_owner_rung_rows, outdir / "owner_rung_shares_all_owners.png", ylabel="share of all owners")

        age_rows = age_profile_rows(sol, P)
        if write_csv_sidecars:
            write_csv(outdir / "age_profiles.csv", age_rows)
        plot_age_profiles(age_rows, outdir / "age_profiles.png")
        tenure_age_rows = tenure_by_age_rows(sol, P)
        if write_csv_sidecars:
            write_csv(outdir / "tenure_by_age.csv", tenure_age_rows)
        plot_tenure_by_age(tenure_age_rows, outdir / "tenure_by_age.png")

    first_look_policy_rows, first_look_market_rows = first_look_rows(sol, P)
    first_look_density_rows = wealth_density_rows(sol, P, wealth_measure="liquid")
    first_look_total_wealth_density_rows = wealth_density_rows(sol, P, wealth_measure="total")
    wealth_grid_rows = wealth_grid_coverage_rows(first_look_density_rows)
    total_wealth_grid_rows = wealth_grid_coverage_rows(first_look_total_wealth_density_rows)
    policy_xlim = density_xlim_from_rows(first_look_density_rows)
    total_wealth_policy_xlim = policy_xlim_from_rows(
        first_look_policy_rows,
        wealth_key="expected_total_wealth_after_tenure",
        base_xlim=policy_xlim,
    )
    if write_csv_sidecars:
        write_csv(outdir / "first_look_policy_lines.csv", first_look_policy_rows)
        write_csv(outdir / "first_look_market_summary.csv", first_look_market_rows)
        write_csv(outdir / "first_look_wealth_density.csv", first_look_density_rows)
        write_csv(outdir / "first_look_total_wealth_density.csv", first_look_total_wealth_density_rows)
        write_csv(outdir / "wealth_grid_coverage.csv", wealth_grid_rows)
        write_csv(outdir / "total_wealth_grid_coverage.csv", total_wealth_grid_rows)
    plot_first_look(
        first_look_policy_rows,
        first_look_market_rows,
        outdir / "first_look_policies_markets.png",
        mode="simple",
        xlim=policy_xlim,
        fast=quick_first_look_only,
    )
    if not quick_first_look_only:
        plot_first_look(
            first_look_policy_rows,
            first_look_market_rows,
            outdir / "first_look_policies_markets_on_path.png",
            mode="simple",
            xlim=policy_xlim,
            mass_filter_min=1e-10,
        )
        plot_first_look(
            first_look_policy_rows,
            first_look_market_rows,
            outdir / "first_look_policies_markets_full.png",
            mode="full",
            xlim=policy_xlim,
        )
    plot_first_look(
        first_look_policy_rows,
        first_look_market_rows,
        outdir / "first_look_policies_markets_total_wealth.png",
        mode="simple",
        xlim=total_wealth_policy_xlim,
        wealth_key="expected_total_wealth_after_tenure",
        wealth_label="expected total wealth after tenure choice",
        fast=quick_first_look_only,
    )
    if not quick_first_look_only:
        plot_first_look(
            first_look_policy_rows,
            first_look_market_rows,
            outdir / "first_look_policies_markets_total_wealth_on_path.png",
            mode="simple",
            xlim=total_wealth_policy_xlim,
            wealth_key="expected_total_wealth_after_tenure",
            wealth_label="expected total wealth after tenure choice",
            mass_filter_min=1e-10,
        )
        plot_first_look(
            first_look_policy_rows,
            first_look_market_rows,
            outdir / "first_look_policies_markets_total_wealth_full.png",
            mode="full",
            xlim=total_wealth_policy_xlim,
            wealth_key="expected_total_wealth_after_tenure",
            wealth_label="expected total wealth after tenure choice",
        )
    plot_first_look_wealth_density(
        first_look_density_rows,
        outdir / "first_look_wealth_density.png",
        wealth_label="liquid wealth b",
        title="Ergodic mass over the occupied liquid-wealth grid",
        fast=quick_first_look_only,
    )
    plot_first_look_wealth_density(
        first_look_total_wealth_density_rows,
        outdir / "first_look_total_wealth_density.png",
        wealth_label="total wealth = b + liquidated housing value",
        title="Ergodic mass over occupied total wealth",
        fast=quick_first_look_only,
    )
    if quick_first_look_only:
        return

    det_policy_rows = ergodic_deterministic_policy_rows(sol, P)
    if write_csv_sidecars:
        write_csv(outdir / "ergodic_deterministic_policy_states.csv", det_policy_rows)
    plot_ergodic_deterministic_policy_slices(
        det_policy_rows,
        outdir / "ergodic_deterministic_policy_slices.png",
    )
    plot_ergodic_deterministic_policy_top_slice(
        det_policy_rows,
        outdir / "ergodic_deterministic_policy_top_slice.png",
    )
    plot_ergodic_deterministic_policy_bins(
        det_policy_rows,
        outdir / "ergodic_deterministic_policy_bins.png",
    )

    threshold_rows = owner_entry_threshold_rows(sol, P)
    if write_csv_sidecars:
        write_csv(outdir / "owner_entry_thresholds.csv", threshold_rows)
    plot_owner_entry_thresholds(threshold_rows, outdir / "owner_entry_thresholds.png")

    policy_rows = owner_entry_policy_rows(sol, P, ages=(30.0, 42.0))
    if write_csv_sidecars:
        write_csv(outdir / "owner_entry_policy_childless_renter_age30_42.csv", policy_rows)

    if include_contact_sheet:
        write_contact_sheet(outdir)


def target_fit_rows(
    moments: dict[str, float],
    targets: dict[str, float],
    weights: dict[str, float],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for name, target in targets.items():
        model = maybe_float(moments.get(name))
        gap = model - float(target) if math.isfinite(model) else math.nan
        weight = float(weights.get(name, 1.0))
        contribution = weight * gap * gap if math.isfinite(gap) else math.inf
        rows.append(
            {
                "moment": name,
                "target": float(target),
                "model": model,
                "gap_model_minus_target": gap,
                "weight": weight,
                "loss_contribution": contribution,
            }
        )
    rows.sort(key=lambda r: float(r["loss_contribution"]), reverse=True)
    return rows


def room_bin_outputs(sol: Any, P: Any) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    model = compute_room_bin_shares(sol, P)
    rows: list[dict[str, Any]] = []
    fit_rows: list[dict[str, Any]] = []
    for tenure in TENURE_ORDER:
        for bin_name in ROOM_BIN_ORDER:
            data_share = ROOM_BIN_DATA[(tenure, bin_name)]
            model_share = float(model.get((tenure, bin_name), math.nan))
            rows.append(
                {
                    "source": "data",
                    "tenure": tenure,
                    "bin": bin_name,
                    "share": data_share,
                    "data_source": ROOM_BIN_DATA_SOURCE,
                }
            )
            rows.append(
                {
                    "source": "model",
                    "tenure": tenure,
                    "bin": bin_name,
                    "share": model_share,
                    "data_source": "current model solve",
                }
            )
            fit_rows.append(
                {
                    "tenure": tenure,
                    "bin": bin_name,
                    "data_share": data_share,
                    "model_share": model_share,
                    "gap_model_minus_data": model_share - data_share,
                }
            )
    return rows, fit_rows


def compute_room_bin_shares(sol: Any, P: Any) -> dict[tuple[str, str], float]:
    bins = {(tenure, bin_name): 0.0 for tenure in TENURE_ORDER for bin_name in ROOM_BIN_ORDER}
    totals = {tenure: 0.0 for tenure in TENURE_ORDER}
    g = np.asarray(sol.g, dtype=float)
    hR = np.asarray(sol.hR_pol, dtype=float)
    ages = np.asarray(P.age_start + np.arange(P.J) * P.da, dtype=float)
    age_idx = np.where((ages >= 30.0) & (ages <= 55.0))[0]
    if g.ndim != 7 or hR.shape != g.shape:
        return {key: math.nan for key in bins}
    for j in age_idx:
        for zz in range(g.shape[4]):
            renter_mass = g[:, 0, 0, j, zz, 0, 0]
            renter_h = hR[:, 0, 0, j, zz, 0, 0]
            for bin_name in ROOM_BIN_ORDER:
                mass = float(np.sum(renter_mass[room_bin_mask(renter_h, bin_name)]))
                bins[("renter", bin_name)] += mass
                totals["renter"] += 0.0
            totals["renter"] += float(np.sum(renter_mass))
            for ten in range(1, 1 + int(P.n_house)):
                owner_mass = float(np.sum(g[:, ten, 0, j, zz, 0, 0]))
                bin_name = room_bin_name(float(P.H_own[ten - 1]))
                bins[("owner", bin_name)] += owner_mass
                totals["owner"] += owner_mass
    return {
        key: float(value / totals[key[0]]) if totals[key[0]] > 1e-14 else math.nan
        for key, value in bins.items()
    }


def room_bin_mask(values: np.ndarray, bin_name: str) -> np.ndarray:
    vals = np.asarray(values, dtype=float)
    if bin_name == "S_<=4":
        return vals <= 4.0
    if bin_name == "M_4to6":
        return (vals > 4.0) & (vals <= 6.0)
    if bin_name == "L_>6":
        return vals > 6.0
    raise ValueError(bin_name)


def room_bin_name(value: float) -> str:
    if value <= 4.0:
        return "S_<=4"
    if value <= 6.0:
        return "M_4to6"
    return "L_>6"


def owner_rung_rows(sol: Any, P: Any) -> list[dict[str, Any]]:
    g = np.asarray(sol.g, dtype=float)
    ages = np.asarray(P.age_start + np.arange(P.J) * P.da, dtype=float)
    age_idx = np.where((ages >= 30.0) & (ages <= 55.0))[0]
    masses = np.zeros(int(P.n_house), dtype=float)
    if g.ndim == 7:
        for j in age_idx:
            for zz in range(g.shape[4]):
                for ten in range(1, 1 + int(P.n_house)):
                    masses[ten - 1] += float(np.sum(g[:, ten, 0, j, zz, 0, 0]))
    total = float(np.sum(masses))
    rows = []
    for idx, rooms in enumerate(np.asarray(P.H_own, dtype=float)):
        rows.append(
            {
                "rung": idx + 1,
                "rooms": float(rooms),
                "mass": float(masses[idx]),
                "share": float(masses[idx] / total) if total > 1e-14 else math.nan,
            }
        )
    return rows


def owner_rung_rows_all(sol: Any, P: Any) -> list[dict[str, Any]]:
    g = np.asarray(sol.g, dtype=float)
    masses = np.zeros(int(P.n_house), dtype=float)
    if g.ndim == 7:
        for ten in range(1, 1 + int(P.n_house)):
            masses[ten - 1] = float(np.sum(g[:, ten, :, :, :, :, :]))
    owner_total = float(np.sum(masses))
    population_total = float(np.sum(g)) if g.ndim == 7 else math.nan
    rows = []
    for idx, rooms in enumerate(np.asarray(P.H_own, dtype=float)):
        mass = float(masses[idx])
        rows.append(
            {
                "rung": idx + 1,
                "rooms": float(rooms),
                "mass": mass,
                "share": float(mass / owner_total) if owner_total > 1e-14 else math.nan,
                "population_share": float(mass / population_total) if population_total > 1e-14 else math.nan,
            }
        )
    return rows


def age_profile_rows(sol: Any, P: Any) -> list[dict[str, Any]]:
    g = np.asarray(sol.g, dtype=float)
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    ages = np.asarray(P.age_start + np.arange(P.J) * P.da, dtype=float)
    rows: list[dict[str, Any]] = []
    for j, age in enumerate(ages):
        gj = g[:, :, :, j, :, :, :]
        mass = float(np.sum(gj))
        mass_b = np.sum(gj, axis=tuple(range(1, gj.ndim))) if g.ndim == 7 else np.zeros_like(b_grid)
        mean_wealth = float(np.sum(mass_b * b_grid) / max(mass, 1e-14))
        rows.append(
            {
                "age": float(age),
                "ownership_rate": float(np.asarray(sol.own_by_age, dtype=float)[j]),
                "fertility": float(np.asarray(getattr(sol, "fert_by_age", np.zeros(P.J)), dtype=float)[j]),
                "mean_housing_services": mean_housing_at_age(sol, P, j),
                "mean_liquid_wealth": mean_wealth,
                "mass": mass,
            }
        )
    return rows


def tenure_by_age_rows(sol: Any, P: Any) -> list[dict[str, Any]]:
    g = np.asarray(sol.g, dtype=float)
    ages = np.asarray(P.age_start + np.arange(P.J) * P.da, dtype=float)
    own_by_age = np.asarray(getattr(sol, "own_by_age", np.zeros(P.J)), dtype=float).reshape(-1)
    rows: list[dict[str, Any]] = []
    for j, age in enumerate(ages):
        gj = g[:, :, :, j, :, :, :] if g.ndim == 7 else np.asarray([])
        mass = float(np.sum(gj)) if g.ndim == 7 else math.nan
        owner_mass = float(np.sum(g[:, 1:, :, j, :, :, :])) if g.ndim == 7 else math.nan
        renter_mass = float(np.sum(g[:, 0, :, j, :, :, :])) if g.ndim == 7 else math.nan
        owner_share = float(own_by_age[j]) if j < own_by_age.size else owner_mass / max(mass, 1e-14)
        renter_share = float(renter_mass / max(mass, 1e-14)) if np.isfinite(mass) else 1.0 - owner_share
        rows.append(
            {
                "age": float(age),
                "renter_share": renter_share,
                "owner_share": owner_share,
                "renter_mass": renter_mass,
                "owner_mass": owner_mass,
                "mass": mass,
            }
        )
    return rows


def mean_housing_at_age(sol: Any, P: Any, j: int) -> float:
    g = np.asarray(sol.g, dtype=float)
    hR = np.asarray(sol.hR_pol, dtype=float)
    if g.ndim != 7 or hR.shape != g.shape:
        return math.nan
    total_h = 0.0
    total_mass = 0.0
    for zz in range(g.shape[4]):
        renter_mass = g[:, 0, 0, j, zz, :, :]
        total_h += float(np.sum(renter_mass * hR[:, 0, 0, j, zz, :, :]))
        total_mass += float(np.sum(renter_mass))
        for ten in range(1, 1 + int(P.n_house)):
            mass = float(np.sum(g[:, ten, 0, j, zz, :, :]))
            total_h += mass * float(P.H_own[ten - 1])
            total_mass += mass
    return float(total_h / max(total_mass, 1e-14))


def owner_entry_threshold_rows(sol: Any, P: Any) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    z_grid = np.asarray(getattr(sol, "type_values", getattr(P, "z_grid", [1.0])), dtype=float).reshape(-1)
    ages = np.asarray(P.age_start + np.arange(P.J) * P.da, dtype=float)
    for j, age in enumerate(ages):
        for zz, z_value in enumerate(z_grid):
            line = owner_entry_line(sol, P, j, zz)
            probs = line[:, 1]
            rows.append(
                {
                    "age": float(age),
                    "z_index": int(zz),
                    "z": float(z_value),
                    "threshold_prob_0p10": first_threshold(line, 0.10),
                    "threshold_prob_0p50": first_threshold(line, 0.50),
                    "threshold_prob_0p90": first_threshold(line, 0.90),
                    "min_owner_probability": safe_nanmin(probs),
                    "max_owner_probability": safe_nanmax(probs),
                    "owner_probability_lowest_wealth": maybe_float(probs[0]) if probs.size else math.nan,
                    "owner_probability_highest_wealth": maybe_float(probs[-1]) if probs.size else math.nan,
                }
            )
    return rows


def owner_entry_policy_rows(sol: Any, P: Any, *, ages: tuple[float, ...]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    z_grid = np.asarray(getattr(sol, "type_values", getattr(P, "z_grid", [1.0])), dtype=float).reshape(-1)
    for age in ages:
        j = int(np.clip(round((age - float(P.age_start)) / max(float(P.da), 1e-12)), 0, int(P.J) - 1))
        actual_age = float(P.age_start + j * P.da)
        for zz, z_value in enumerate(z_grid):
            line = owner_entry_line(sol, P, j, zz)
            for wealth, prob in line:
                rows.append(
                    {
                        "requested_age": float(age),
                        "age": actual_age,
                        "z_index": int(zz),
                        "z": float(z_value),
                        "liquid_wealth": float(wealth),
                        "owner_entry_probability": float(prob),
                    }
                )
    return rows


def owner_entry_line(sol: Any, P: Any, j: int, zz: int) -> np.ndarray:
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    tp = getattr(sol, "tenure_probs", None)
    if tp is None:
        probs = (np.asarray(sol.tenure_choice)[:, 0, 0, j, zz, 0, 0] > 0).astype(float)
    else:
        probs = np.sum(np.asarray(tp)[:, 0, 0, j, zz, 0, 0, 1:], axis=1)
    V = np.asarray(getattr(sol, "V", np.empty_like(probs)))
    if V.shape[:7] == np.asarray(sol.hR_pol).shape:
        valid = V[:, 0, 0, j, zz, 0, 0] > -1e9
        probs = np.where(valid, probs, np.nan)
    return np.column_stack([b_grid, probs])


def first_look_rows(sol: Any, P: Any) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    ages = selected_policy_ages(P)
    z_indices = selected_income_indices(sol, P)
    z_grid = np.asarray(getattr(sol, "type_values", getattr(P, "z_grid", [1.0])), dtype=float).reshape(-1)
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    c_pol = np.asarray(sol.c_pol, dtype=float)
    hR_pol = np.asarray(sol.hR_pol, dtype=float)
    bp_pol = np.asarray(sol.bp_pol, dtype=float)
    tenure_choice = np.asarray(sol.tenure_choice)
    tenure_probs = getattr(sol, "tenure_probs", None)
    tenure_probs = None if tenure_probs is None else np.asarray(tenure_probs, dtype=float)
    fert_probs = np.asarray(sol.fert_probs, dtype=float)
    valid_source = np.asarray(getattr(sol, "V", np.empty_like(c_pol)), dtype=float)
    g = np.asarray(getattr(sol, "g", np.zeros(0)), dtype=float)
    nvec = np.arange(int(P.n_parity), dtype=float)
    asset_price = owner_asset_price_vector(sol, P)
    liquidation_discount = 1.0 - float(getattr(P, "psi", 0.0))
    policy_rows: list[dict[str, Any]] = []
    for age in ages:
        j = age_to_index_for_packet(P, float(age))
        actual_age = float(P.age_start + j * P.da)
        for zz in z_indices:
            z_value = float(z_grid[zz]) if zz < len(z_grid) else float(zz)
            mass_line = childless_renter_mass_line(g, len(b_grid), j, zz)
            slice_mass = float(np.nansum(mass_line))
            for bb, wealth in enumerate(b_grid):
                if not is_valid_policy_point(valid_source, bb, j, zz):
                    continue
                tchoice = int(tenure_choice[bb, 0, 0, j, zz, 0, 0])
                branch_wealth = modal_branch_wealth_from_childless_renter(float(wealth), tchoice, P, asset_price)
                branch_wealth_clipped = float(np.clip(branch_wealth, b_grid[0], b_grid[-1]))
                modal_consumption = interp_policy_scalar(b_grid, c_pol[:, tchoice, 0, j, zz, 0, 0], branch_wealth)
                modal_next_liquid_wealth = interp_policy_scalar(
                    b_grid, bp_pol[:, tchoice, 0, j, zz, 0, 0], branch_wealth
                )
                renter_housing_if_renter = interp_policy_scalar(b_grid, hR_pol[:, 0, 0, j, zz, 0, 0], float(wealth))
                chosen_renter_housing = (
                    interp_policy_scalar(b_grid, hR_pol[:, 0, 0, j, zz, 0, 0], branch_wealth)
                    if tchoice <= 0
                    else math.nan
                )
                housing = (
                    chosen_renter_housing
                    if tchoice <= 0
                    else float(np.asarray(P.H_own, dtype=float)[max(0, min(tchoice - 1, int(P.n_house) - 1))])
                )
                owner_prob = owner_probability_from_tenure_probs(tenure_probs, bb, j, zz, tchoice)
                ergodic_mass = maybe_float(mass_line[bb]) if bb < len(mass_line) else math.nan
                liquidated_value_after_tenure = liquidated_housing_value(P, asset_price, tchoice, 0)
                total_wealth_after_tenure = (
                    float(branch_wealth + liquidated_value_after_tenure)
                    if math.isfinite(liquidated_value_after_tenure)
                    else math.nan
                )
                expected = expected_childless_renter_tenure_policy(
                    b_grid=b_grid,
                    c_pol=c_pol,
                    hR_pol=hR_pol,
                    bp_pol=bp_pol,
                    tenure_probs=tenure_probs,
                    bb=bb,
                    j=j,
                    zz=zz,
                    modal_tenure=tchoice,
                    liquid_wealth=float(wealth),
                    P=P,
                    asset_price=asset_price,
                )
                policy_rows.append(
                    {
                        "initial_state": "childless_renter",
                        "age": actual_age,
                        "z_index": int(zz),
                        "z": z_value,
                        "liquid_wealth": float(wealth),
                        "total_wealth_before_choice": float(wealth),
                        "total_wealth_after_tenure": total_wealth_after_tenure,
                        "branch_liquid_wealth_after_transaction": branch_wealth,
                        "branch_liquid_wealth_clipped_to_grid": branch_wealth_clipped,
                        "expected_branch_liquid_wealth_after_transaction": expected["branch_liquid_wealth"],
                        "expected_branch_liquid_wealth_clipped_to_grid": expected["branch_liquid_wealth_clipped"],
                        "liquidated_housing_value_after_tenure": liquidated_value_after_tenure,
                        "expected_liquidated_housing_value_after_tenure": expected["liquidated_housing_value"],
                        "liquidated_housing_discount_1_minus_psi": liquidation_discount,
                        "consumption": expected["consumption"],
                        "next_liquid_wealth": expected["next_liquid_wealth"],
                        "modal_consumption": modal_consumption,
                        "modal_next_liquid_wealth": modal_next_liquid_wealth,
                        "renter_housing_policy": renter_housing_if_renter,
                        "chosen_renter_housing": chosen_renter_housing,
                        "chosen_tenure_index": int(tchoice),
                        "chosen_tenure": "renter" if tchoice <= 0 else "owner",
                        "owner_entry_probability": owner_prob,
                        "housing_services_after_tenure": expected["housing_services"],
                        "modal_housing_services_after_tenure": housing,
                        "expected_total_wealth_after_tenure": expected["total_wealth"],
                        "expected_children": fertility_expectation(fert_probs, bb, j, zz, nvec),
                        "ergodic_mass": ergodic_mass,
                        "mass_within_age_z_childless_renter": (
                            ergodic_mass / slice_mass
                            if math.isfinite(ergodic_mass) and slice_mass > 1e-14
                            else math.nan
                        ),
                    }
                )

    markets = np.arange(int(P.I))
    owner_user_cost = np.asarray(getattr(sol, "owner_user_cost", np.full(P.I, np.nan)), dtype=float).reshape(-1)
    rental_demand = np.asarray(getattr(sol, "rental_demand_by_market", np.zeros(P.I)), dtype=float).reshape(-1)
    owner_demand = np.asarray(getattr(sol, "owner_demand_by_market", np.zeros(P.I)), dtype=float).reshape(-1)
    supply = np.asarray(getattr(sol, "housing_supply", np.full(P.I, np.nan)), dtype=float).reshape(-1)
    user_cost_rate = float(getattr(P, "user_cost_rate", math.nan))
    q_rate = float(getattr(P, "q", math.nan))
    depreciation = float(getattr(P, "delta", math.nan))
    property_tax = float(getattr(P, "tau_H", math.nan))
    rent_normalizer = np.asarray(getattr(P, "r_bar", np.full(P.I, np.nan)), dtype=float).reshape(-1)
    market_rows: list[dict[str, Any]] = []
    for i in markets:
        rental = maybe_vector_value(rental_demand, i)
        owner = maybe_vector_value(owner_demand, i)
        total = rental + owner if math.isfinite(rental) and math.isfinite(owner) else math.nan
        market_rows.append(
            {
                "market": int(i + 1),
                "rental_demand": rental,
                "owner_demand": owner,
                "total_demand": total,
                "supply": maybe_vector_value(supply, i),
                "house_price": maybe_vector_value(asset_price, i),
                "asset_price": maybe_vector_value(asset_price, i),
                "flow_rent_or_user_cost": maybe_vector_value(owner_user_cost, i),
                "user_cost_rate": user_cost_rate,
                "interest_rate_q": q_rate,
                "depreciation_delta": depreciation,
                "property_tax_tau_H": property_tax,
                "rent_normalizer_r_bar": maybe_vector_value(rent_normalizer, i),
            }
        )
    return policy_rows, market_rows


def selected_policy_ages(P: Any) -> list[float]:
    ages = np.asarray(P.age_start + np.arange(P.J) * P.da, dtype=float)
    fertile_start = float(P.age_start + (int(P.A_f_start) - 1) * P.da)
    fertile_end = float(P.age_start + (int(P.A_f_end) - 1) * P.da)
    requested = np.array(
        [
            fertile_start,
            0.5 * (fertile_start + fertile_end),
            fertile_end,
        ],
        dtype=float,
    )
    out: list[float] = []
    for age in requested:
        actual = float(ages[np.argmin(np.abs(ages - age))])
        if actual not in out:
            out.append(actual)
    return out


def selected_income_indices(sol: Any, P: Any) -> list[int]:
    z_grid = np.asarray(getattr(sol, "type_values", getattr(P, "z_grid", [1.0])), dtype=float).reshape(-1)
    if len(z_grid) <= 3:
        return list(range(len(z_grid)))
    return sorted({0, len(z_grid) // 2, len(z_grid) - 1})


def age_to_index_for_packet(P: Any, age: float) -> int:
    idx = int(round((float(age) - float(P.age_start)) / max(float(P.da), 1e-12)))
    return int(np.clip(idx, 0, int(P.J) - 1))


def is_valid_policy_point(values: np.ndarray, bb: int, j: int, zz: int) -> bool:
    if values.ndim < 7:
        return True
    try:
        return bool(values[bb, 0, 0, j, zz, 0, 0] > -1e9)
    except IndexError:
        return True


def fertility_expectation(fert_probs: np.ndarray, bb: int, j: int, zz: int, nvec: np.ndarray) -> float:
    if fert_probs.ndim == 8:
        probs = fert_probs[bb, 0, 0, j, zz, 0, 0, :]
    elif fert_probs.ndim == 6:
        probs = fert_probs[bb, 0, 0, j, zz, :]
    elif fert_probs.ndim == 5:
        probs = fert_probs[bb, 0, 0, j, :]
    else:
        return math.nan
    probs = np.asarray(probs, dtype=float).reshape(-1)
    n = min(len(probs), len(nvec))
    return float(np.dot(probs[:n], nvec[:n])) if n else math.nan


def childless_renter_mass_line(g: np.ndarray, nb: int, j: int, zz: int) -> np.ndarray:
    if g.ndim != 7:
        return np.full(nb, math.nan)
    try:
        return np.asarray(g[:, 0, 0, j, zz, 0, 0], dtype=float).reshape(-1)
    except IndexError:
        return np.full(nb, math.nan)


def owner_probability_from_tenure_probs(
    tenure_probs: np.ndarray | None,
    bb: int,
    j: int,
    zz: int,
    tenure_choice: int,
) -> float:
    if tenure_probs is None:
        return float(tenure_choice > 0)
    if tenure_probs.ndim < 8:
        return math.nan
    try:
        return float(np.sum(tenure_probs[bb, 0, 0, j, zz, 0, 0, 1:]))
    except IndexError:
        return math.nan


def tenure_probability_vector(
    tenure_probs: np.ndarray | None,
    bb: int,
    j: int,
    zz: int,
    modal_tenure: int,
    nt: int,
) -> np.ndarray:
    if tenure_probs is None or tenure_probs.ndim < 8:
        probs = np.zeros(nt, dtype=float)
        probs[int(np.clip(modal_tenure, 0, nt - 1))] = 1.0
        return probs
    try:
        probs = np.asarray(tenure_probs[bb, 0, 0, j, zz, 0, 0, :nt], dtype=float).reshape(-1)
    except IndexError:
        probs = np.zeros(nt, dtype=float)
        probs[int(np.clip(modal_tenure, 0, nt - 1))] = 1.0
        return probs
    total = float(np.nansum(probs))
    if not math.isfinite(total) or total <= 0.0:
        out = np.zeros(nt, dtype=float)
        out[int(np.clip(modal_tenure, 0, nt - 1))] = 1.0
        return out
    return probs / total


def expected_childless_renter_tenure_policy(
    *,
    b_grid: np.ndarray,
    c_pol: np.ndarray,
    hR_pol: np.ndarray,
    bp_pol: np.ndarray,
    tenure_probs: np.ndarray | None,
    bb: int,
    j: int,
    zz: int,
    modal_tenure: int,
    liquid_wealth: float,
    P: Any,
    asset_price: np.ndarray,
) -> dict[str, float]:
    nt = 1 + int(getattr(P, "n_house", 0))
    probs = tenure_probability_vector(tenure_probs, bb, j, zz, modal_tenure, nt)
    h_own = np.asarray(getattr(P, "H_own", []), dtype=float).reshape(-1)
    totals = {
        "consumption": 0.0,
        "next_liquid_wealth": 0.0,
        "housing_services": 0.0,
        "branch_liquid_wealth": 0.0,
        "branch_liquid_wealth_clipped": 0.0,
        "liquidated_housing_value": 0.0,
        "total_wealth": 0.0,
    }
    weight = 0.0
    for tn, prob in enumerate(probs):
        p = float(prob)
        if p <= 1e-14 or not math.isfinite(p):
            continue
        branch_wealth = modal_branch_wealth_from_childless_renter(liquid_wealth, tn, P, asset_price)
        if not math.isfinite(branch_wealth):
            continue
        clipped = float(np.clip(branch_wealth, b_grid[0], b_grid[-1]))
        c_val = interp_policy_scalar(b_grid, c_pol[:, tn, 0, j, zz, 0, 0], branch_wealth)
        bp_val = interp_policy_scalar(b_grid, bp_pol[:, tn, 0, j, zz, 0, 0], branch_wealth)
        if tn <= 0:
            h_val = interp_policy_scalar(b_grid, hR_pol[:, 0, 0, j, zz, 0, 0], branch_wealth)
        elif tn - 1 < h_own.size:
            h_val = float(h_own[tn - 1])
        else:
            h_val = math.nan
        liq_val = liquidated_housing_value(P, asset_price, tn, 0)
        if not all(math.isfinite(x) for x in [c_val, bp_val, h_val, liq_val]):
            continue
        totals["consumption"] += p * c_val
        totals["next_liquid_wealth"] += p * bp_val
        totals["housing_services"] += p * h_val
        totals["branch_liquid_wealth"] += p * branch_wealth
        totals["branch_liquid_wealth_clipped"] += p * clipped
        totals["liquidated_housing_value"] += p * liq_val
        totals["total_wealth"] += p * (branch_wealth + liq_val)
        weight += p
    if weight <= 0.0:
        return {key: math.nan for key in totals}
    return {key: value / weight for key, value in totals.items()}


def maybe_vector_value(values: np.ndarray, index: int) -> float:
    arr = np.asarray(values, dtype=float).reshape(-1)
    if index >= len(arr):
        return math.nan
    return maybe_float(arr[index])


def interp_policy_scalar(grid: np.ndarray, values: np.ndarray, query: float) -> float:
    x = np.asarray(grid, dtype=float).reshape(-1)
    y = np.asarray(values, dtype=float).reshape(-1)
    if x.size == 0 or y.size != x.size or not math.isfinite(query):
        return math.nan
    q = float(np.clip(query, x[0], x[-1]))
    return float(np.interp(q, x, y))


def owner_asset_price_vector(sol: Any, P: Any) -> np.ndarray:
    fallback = np.full(int(getattr(P, "I", 1)), math.nan)
    values = getattr(sol, "owner_asset_price", getattr(sol, "p_eq", fallback))
    arr = np.asarray(values, dtype=float).reshape(-1)
    return arr if arr.size else fallback


def liquidated_housing_value(P: Any, asset_price: np.ndarray, tenure_index: int, market_index: int = 0) -> float:
    if int(tenure_index) <= 0:
        return 0.0
    h_own = np.asarray(getattr(P, "H_own", []), dtype=float).reshape(-1)
    if h_own.size == 0:
        return math.nan
    house_idx = int(np.clip(int(tenure_index) - 1, 0, h_own.size - 1))
    price = maybe_vector_value(asset_price, market_index)
    if not math.isfinite(price) and asset_price.size:
        price = maybe_vector_value(asset_price, 0)
    if not math.isfinite(price):
        return math.nan
    liquidation_discount = 1.0 - float(getattr(P, "psi", 0.0))
    return float(liquidation_discount * price * h_own[house_idx])


def modal_branch_wealth_from_childless_renter(
    liquid_wealth: float,
    tenure_index: int,
    P: Any,
    asset_price: np.ndarray,
) -> float:
    if int(tenure_index) <= 0:
        return float(liquid_wealth)
    h_own = np.asarray(getattr(P, "H_own", []), dtype=float).reshape(-1)
    if h_own.size == 0:
        return math.nan
    house_idx = int(np.clip(int(tenure_index) - 1, 0, h_own.size - 1))
    price = maybe_vector_value(asset_price, 0)
    if not math.isfinite(price):
        return math.nan
    return float(liquid_wealth - price * h_own[house_idx])


def first_threshold(line: np.ndarray, cutoff: float) -> float:
    if line.size == 0:
        return math.nan
    valid = np.isfinite(line[:, 1])
    hit = valid & (line[:, 1] >= float(cutoff))
    if not np.any(hit):
        return math.nan
    return float(line[np.argmax(hit), 0])


def ergodic_deterministic_policy_rows(sol: Any, P: Any, *, mass_min: float = 1e-12) -> list[dict[str, Any]]:
    g = np.asarray(getattr(sol, "g", np.zeros(0)), dtype=float)
    if g.ndim != 7:
        return []
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    c_pol = np.asarray(sol.c_pol, dtype=float)
    hR_pol = np.asarray(sol.hR_pol, dtype=float)
    bp_pol = np.asarray(sol.bp_pol, dtype=float)
    tenure_choice = np.asarray(sol.tenure_choice)
    fert_probs = np.asarray(sol.fert_probs, dtype=float)
    price = owner_asset_price_vector(sol, P)
    h_own = np.asarray(P.H_own, dtype=float).reshape(-1)
    z_grid = np.asarray(getattr(sol, "type_values", getattr(P, "z_grid", [1.0])), dtype=float).reshape(-1)
    age_grid = np.asarray(float(P.age_start) + np.arange(int(P.J)) * float(P.da), dtype=float)
    shared = precompute_shared(P, b_grid)
    total_mass = float(np.sum(g))
    rows: list[dict[str, Any]] = []
    nz = np.argwhere(g > float(mass_min))
    for bb, to, i, j, zz, nn, cs in nz:
        mass = float(g[bb, to, i, j, zz, nn, cs])
        b = float(b_grid[bb])
        tn = int(tenure_choice[bb, to, i, j, zz, nn, cs])
        bw = deterministic_branch_wealth(
            b=b,
            origin_tenure=int(to),
            target_tenure=tn,
            market=int(i),
            parity=int(nn),
            child_state=int(cs),
            P=P,
            price=price,
            shared=shared,
        )
        c_val = interp_policy_scalar(b_grid, c_pol[:, tn, i, j, zz, nn, cs], bw)
        bp_val = interp_policy_scalar(b_grid, bp_pol[:, tn, i, j, zz, nn, cs], bw)
        if tn <= 0:
            h_val = interp_policy_scalar(b_grid, hR_pol[:, 0, i, j, zz, nn, cs], bw)
            target_rooms = 0.0
            target_label = "renter"
        else:
            target_rooms = float(h_own[tn - 1])
            h_val = target_rooms
            target_label = f"owner_{target_rooms:g}"
        fert_choice, expected_children = deterministic_fertility_choice(fert_probs, bb, to, i, j, zz, nn)
        completed_children = completed_children_for_packet(P, int(nn), int(cs))
        current_rooms = 0.0 if to <= 0 else float(h_own[to - 1])
        rows.append(
            {
                "age": float(age_grid[j]),
                "age_index": int(j),
                "z_index": int(zz),
                "z": float(z_grid[zz]) if zz < z_grid.size else float(zz),
                "current_tenure_index": int(to),
                "current_tenure": "renter" if to <= 0 else f"owner_{current_rooms:g}",
                "current_tenure_rooms": current_rooms,
                "parity_state": int(nn),
                "child_state": int(cs),
                "completed_children": int(completed_children),
                "liquid_wealth": b,
                "ergodic_mass": mass,
                "population_share": mass / total_mass if total_mass > 1e-14 else math.nan,
                "target_tenure_index": tn,
                "target_tenure": target_label,
                "target_tenure_rooms": target_rooms,
                "deterministic_owner_choice": float(tn > 0),
                "branch_liquid_wealth_after_transaction": bw,
                "consumption_policy": c_val,
                "next_liquid_wealth_policy": bp_val,
                "housing_services_policy": h_val,
                "deterministic_fertility_choice": fert_choice,
                "expected_children_probability_weighted": expected_children,
            }
        )
    return rows


def deterministic_branch_wealth(
    *,
    b: float,
    origin_tenure: int,
    target_tenure: int,
    market: int,
    parity: int,
    child_state: int,
    P: Any,
    price: np.ndarray,
    shared: Any,
) -> float:
    to = int(origin_tenure)
    tn = int(target_tenure)
    if tn == to:
        return float(b)
    p = maybe_vector_value(price, int(market))
    if not math.isfinite(p):
        return math.nan
    h_own = np.asarray(P.H_own, dtype=float).reshape(-1)
    sale = 0.0 if to <= 0 else (1.0 - float(P.psi)) * p * float(h_own[to - 1])
    if tn <= 0:
        return float(max(b + sale, 0.0))
    hcost = p * float(h_own[tn - 1])
    if to <= 0:
        branch = float(b - hcost)
        if bool(shared.birth_dp[parity, child_state, to, tn]):
            phi = float(shared.phi_choice[market, tn, parity, child_state])
            branch = max(branch, -phi * hcost)
        else:
            grant = float(shared.birth_entry_grant[market, tn, parity, child_state])
            if grant > 0.0:
                branch += grant
        return branch
    return float(b + sale - hcost)


def deterministic_fertility_choice(
    fert_probs: np.ndarray,
    bb: int,
    to: int,
    i: int,
    j: int,
    zz: int,
    nn: int,
) -> tuple[float, float]:
    try:
        if fert_probs.ndim == 6:
            probs = np.asarray(fert_probs[bb, to, i, j, zz, :], dtype=float).reshape(-1)
        else:
            probs = np.asarray(fert_probs[bb, to, i, j, zz, nn, :], dtype=float).reshape(-1)
    except IndexError:
        return math.nan, math.nan
    if probs.size == 0 or not np.any(np.isfinite(probs)):
        return math.nan, math.nan
    probs = np.nan_to_num(probs, nan=0.0, posinf=0.0, neginf=0.0)
    total = float(np.sum(probs))
    if total > 0.0:
        probs = probs / total
    grid = np.arange(probs.size, dtype=float)
    return float(np.argmax(probs)), float(np.dot(probs, grid))


def completed_children_for_packet(P: Any, nn: int, cs: int) -> int:
    k = int(getattr(P, "n_child_stages", 1))
    if cs == 0:
        return 0
    if cs == k + 1:
        return 1
    if cs == k + 2:
        return 2
    return int(nn)


def plot_ergodic_deterministic_policy_slices(rows: list[dict[str, Any]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    selected = selected_childless_renter_policy_slices(rows, max_slices=4)
    if not selected:
        return
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 8.4), sharex=True)
    specs = [
        (axes[0, 0], "consumption_policy", "Consumption", "consumption"),
        (axes[0, 1], "housing_services_policy", "Housing services", "housing services"),
        (axes[1, 0], "target_tenure_rooms", "Deterministic tenure/product", "0=renter, owner rooms"),
        (axes[1, 1], "deterministic_fertility_choice", "Deterministic fertility choice", "children"),
    ]
    colors = plt.cm.viridis(np.linspace(0.12, 0.88, len(selected)))
    for color, ((age, z), panel) in zip(colors, selected):
        panel = sorted(panel, key=lambda r: maybe_float(r["liquid_wealth"]))
        x = [maybe_float(r["liquid_wealth"]) for r in panel]
        mass = sum(maybe_float(r["population_share"]) for r in panel)
        label = f"age={age:g}, z={z:g}, mass={100.0 * mass:.2f}%"
        for ax, key, _, _ in specs:
            ax.plot(x, [maybe_float(r[key]) for r in panel], marker="o", ms=3.0, lw=1.8, color=color, label=label)
    for ax, _, title, ylabel in specs:
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.2)
    for ax in axes[1, :]:
        ax.set_xlabel("liquid wealth b")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center", ncols=2, frameon=False, fontsize=8)
    fig.suptitle("Deterministic argmax policies on common occupied childless-renter states", y=0.985)
    fig.tight_layout(rect=(0.0, 0.08, 1.0, 0.95))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_ergodic_deterministic_policy_top_slice(rows: list[dict[str, Any]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    selected = selected_childless_renter_policy_slices(rows, max_slices=12)
    if not selected:
        return
    key, panel = max(
        selected,
        key=lambda item: (len(item[1]) >= 8, sum(maybe_float(r["ergodic_mass"]) for r in item[1])),
    )
    panel = sorted(panel, key=lambda r: maybe_float(r["liquid_wealth"]))
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.8), sharex=True)
    specs = [
        (axes[0, 0], "consumption_policy", "Consumption", "consumption"),
        (axes[0, 1], "housing_services_policy", "Housing services", "housing services"),
        (axes[1, 0], "target_tenure_rooms", "Deterministic tenure/product", "0=renter, owner rooms"),
        (axes[1, 1], "deterministic_fertility_choice", "Deterministic fertility choice", "children"),
    ]
    x = [maybe_float(r["liquid_wealth"]) for r in panel]
    for ax, field, title, ylabel in specs:
        ax.plot(x, [maybe_float(r[field]) for r in panel], color="tab:blue", marker="o", ms=3.5, lw=1.8)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.2)
    for ax in axes[1, :]:
        ax.set_xlabel("liquid wealth b")
    mass = sum(maybe_float(r["population_share"]) for r in panel)
    age, z = key
    fig.suptitle(
        f"Deterministic argmax policy for one common occupied state: childless renter, age={age:g}, z={z:g}, mass={100.0 * mass:.2f}%",
        y=0.985,
    )
    fig.tight_layout(rect=(0.0, 0.03, 1.0, 0.93))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def selected_childless_renter_policy_slices(
    rows: list[dict[str, Any]],
    *,
    max_slices: int,
) -> list[tuple[tuple[float, float], list[dict[str, Any]]]]:
    groups: dict[tuple[float, float], list[dict[str, Any]]] = {}
    for row in rows:
        if int(row["current_tenure_index"]) != 0:
            continue
        if int(row["completed_children"]) != 0 or int(row["child_state"]) != 0:
            continue
        key = (float(row["age"]), float(row["z"]))
        groups.setdefault(key, []).append(row)
    ranked = []
    for key, panel in groups.items():
        if len(panel) < 3:
            continue
        mass = sum(maybe_float(r["ergodic_mass"]) for r in panel)
        ranked.append((mass, key, panel))
    ranked.sort(reverse=True, key=lambda x: x[0])
    return [(key, panel) for _, key, panel in ranked[:max_slices]]


def plot_ergodic_deterministic_policy_bins(rows: list[dict[str, Any]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    binned = deterministic_policy_by_liquid_wealth(rows)
    if not binned:
        return
    x = np.asarray([r["liquid_wealth"] for r in binned], dtype=float)
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 8.4), sharex=True)
    specs = [
        (axes[0, 0], "population_share", "Ergodic mass", "population share"),
        (axes[0, 1], "mean_consumption", "Mass-weighted consumption", "consumption"),
        (axes[1, 0], "mean_housing_services", "Mass-weighted housing", "housing services"),
        (axes[1, 1], "owner_choice_share", "Argmax owner choice share", "share"),
    ]
    for ax, key, title, ylabel in specs:
        ax.plot(x, [maybe_float(r[key]) for r in binned], marker="o", ms=2.5, lw=1.6)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.2)
    for ax in axes[1, :]:
        ax.set_xlabel("liquid wealth b")
    fig.suptitle("Mass-weighted deterministic argmax policies over all occupied states", y=0.985)
    fig.tight_layout(rect=(0.0, 0.03, 1.0, 0.95))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def deterministic_policy_by_liquid_wealth(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_b: dict[float, list[dict[str, Any]]] = {}
    for row in rows:
        by_b.setdefault(float(row["liquid_wealth"]), []).append(row)
    out: list[dict[str, Any]] = []
    total_mass = sum(maybe_float(r["ergodic_mass"]) for r in rows)
    for b, panel in sorted(by_b.items()):
        mass = sum(maybe_float(r["ergodic_mass"]) for r in panel)
        if mass <= 0.0:
            continue
        out.append(
            {
                "liquid_wealth": float(b),
                "population_share": mass / total_mass if total_mass > 1e-14 else math.nan,
                "mean_consumption": weighted_row_mean(panel, "consumption_policy", mass),
                "mean_housing_services": weighted_row_mean(panel, "housing_services_policy", mass),
                "owner_choice_share": weighted_row_mean(panel, "deterministic_owner_choice", mass),
                "mean_deterministic_fertility_choice": weighted_row_mean(
                    panel, "deterministic_fertility_choice", mass
                ),
            }
        )
    return out


def weighted_row_mean(rows: list[dict[str, Any]], key: str, mass: float) -> float:
    if mass <= 0.0:
        return math.nan
    total = 0.0
    used_mass = 0.0
    for row in rows:
        value = maybe_float(row.get(key))
        weight = maybe_float(row.get("ergodic_mass"))
        if math.isfinite(value) and math.isfinite(weight) and weight > 0.0:
            total += weight * value
            used_mass += weight
    return total / used_mass if used_mass > 0.0 else math.nan


def write_policy_cases(
    outdir: Path,
    *,
    source: dict[str, Any],
    grid: dict[str, Any],
    baseline: dict[str, Any],
    targets: dict[str, float],
    weights: dict[str, float],
    write_case_diagnostics: bool,
) -> list[dict[str, Any]]:
    outdir.mkdir(parents=True, exist_ok=True)
    cases = [
        {
            "case": "baseline",
            "label": "Baseline",
            "note": "Same solved theta as the mechanics packet baseline.",
            "overrides": {},
            "result": baseline,
        },
        {
            "case": "parent_ltv95",
            "label": "Parent LTV relief",
            "note": "Sets financed share to 0.95 for new-parent owner choices in the birth child-state.",
            "overrides": {
                "parent_dp_waiver": True,
                "parent_dp_waiver_phi": 0.95,
                "parent_dp_waiver_birth_state_only": True,
            },
        },
        {
            "case": "property_tax_up_1pp",
            "label": "Property tax +1pp annual",
            "note": "Raises annual property tax from 1 percent to 2 percent and re-clears the housing market.",
            "overrides": {"tau_H": 0.02 * 4.0},
        },
        {
            "case": "estate_tax_30pct",
            "label": "Estate tax wedge 30 percent",
            "note": "Terminal bequest-tax wedge only; no government budget, rebate, or inheritance kernel.",
            "overrides": {"estate_tax_rate": 0.30},
        },
    ]
    records: list[dict[str, Any]] = []
    baseline_moments = dict(baseline["moments"])
    for case in cases:
        result = case.get("result")
        if result is None:
            result = solve_candidate(
                theta=source["theta"],
                grid=grid,
                extra_overrides=dict(case["overrides"]),
                targets=targets,
                weights=weights,
                label=str(case["case"]),
            )
        case_dir = outdir / str(case["case"])
        case_dir.mkdir(parents=True, exist_ok=True)
        if write_case_diagnostics:
            write_diagnostics(result["sol"], result["P"], case_dir / "diagnostics")
        record = policy_record(case, result, targets, baseline_moments)
        write_json(case_dir / "record.json", record)
        records.append(record)
        write_csv(outdir / "policy_cases_partial.csv", records)
        write_json(outdir / "policy_cases_partial.json", records)
        print(
            f"policy {case['case']}: loss={record['rank_loss']:.6g}, "
            f"resid={record['market_residual']:.2e}, elapsed={record['elapsed_sec']:.1f}s",
            flush=True,
        )
    write_csv(outdir / "policy_cases.csv", records)
    write_json(outdir / "policy_cases.json", records)
    plot_policy_moment_comparison(records, targets, outdir / "policy_moment_comparison.png")
    return records


def policy_record(
    case: dict[str, Any],
    result: dict[str, Any],
    targets: dict[str, float],
    baseline_moments: dict[str, float],
) -> dict[str, Any]:
    moments = dict(result["moments"])
    key_moments = sorted(
        set(targets)
        | {
            "mean_age_first_birth",
            "aggregate_own_rate",
            "prime30_55_childless_owner_minus_renter_mean_rooms",
            "renter25_45_all_cap_share",
        }
    )
    return {
        "case": str(case["case"]),
        "label": str(case["label"]),
        "note": str(case["note"]),
        "status": "diagnostic_mechanics_not_policy_estimate",
        "overrides": jsonable(case["overrides"]),
        "rank_loss": float(result["rank_loss"]),
        "market_residual": float(result["market_residual"]),
        "elapsed_sec": float(result["elapsed_sec"]),
        "p_eq": jsonable(result["p_eq"]),
        "H_own": jsonable(np.asarray(result["P"].H_own, dtype=float)),
        "hR_max": float(result["P"].hR_max),
        "moments": {k: maybe_float(moments.get(k)) for k in key_moments},
        "delta_from_baseline": {
            k: maybe_float(moments.get(k)) - maybe_float(baseline_moments.get(k))
            for k in key_moments
            if math.isfinite(maybe_float(moments.get(k))) and math.isfinite(maybe_float(baseline_moments.get(k)))
        },
    }


def plot_first_look(
    policy_rows: list[dict[str, Any]],
    market_rows: list[dict[str, Any]],
    path: Path,
    *,
    mode: str,
    xlim: tuple[float, float] | None = None,
    wealth_key: str = "liquid_wealth",
    wealth_label: str = "liquid wealth",
    mass_filter_min: float | None = None,
    fast: bool = False,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not policy_rows and not market_rows:
        return
    plot_rows = select_first_look_policy_rows(policy_rows, mode=mode)
    ages = sorted({float(r["age"]) for r in plot_rows})
    z_values = sorted({float(r["z"]) for r in plot_rows})
    colors = plt.cm.viridis(np.linspace(0.12, 0.88, max(len(z_values), 1)))
    color_for_z = {z: colors[i] for i, z in enumerate(z_values)}
    age_styles = ["-", "--", ":"]
    style_for_age = {age: age_styles[i % len(age_styles)] for i, age in enumerate(ages)}

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.2) if fast else (13.0, 8.8))
    ax_c, ax_h, ax_m, ax_f = axes.ravel()
    policy_specs = [
        (ax_c, "consumption", "Consumption policy", "consumption"),
        (ax_f, "expected_children", "Fertility policy", "expected children"),
    ]
    for ax, key, title, ylabel in policy_specs:
        for z in z_values:
            for age in ages:
                panel = [r for r in plot_rows if float(r["z"]) == z and float(r["age"]) == age]
                panel.sort(key=lambda r: maybe_float(r.get(wealth_key)))
                panel = filter_policy_plot_rows(
                    panel,
                    xlim,
                    wealth_key=wealth_key,
                    mass_filter_min=mass_filter_min,
                )
                if not panel:
                    continue
                ax.plot(
                    [maybe_float(r.get(wealth_key)) for r in panel],
                    [maybe_float(r[key]) for r in panel],
                    color=color_for_z[z],
                    linestyle=style_for_age[age],
                    lw=1.7,
                    label=f"z={z:g}, age={age:g}",
                )
        ax.set_title(title)
        ax.set_xlabel(wealth_label)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.2)
        if xlim is not None:
            ax.set_xlim(*xlim)

    plot_housing_by_tenure(
        ax_h,
        plot_rows,
        z_values,
        ages,
        color_for_z,
        style_for_age,
        xlim=xlim,
        wealth_key=wealth_key,
        wealth_label=wealth_label,
        mass_filter_min=mass_filter_min,
    )
    if xlim is not None:
        ax_h.set_xlim(*xlim)
    plot_market_panel(ax_m, market_rows)
    handles, labels = ax_c.get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center", ncols=min(3, len(labels)), frameon=False, fontsize=8)
    mode_label = "simple" if mode == "simple" else "full"
    fig.suptitle(
        f"First-look model mechanics ({mode_label}, x-axis: {wealth_label}): policies and housing-market equilibrium",
        y=0.985,
    )
    fig.tight_layout(rect=(0.0, 0.07, 1.0, 0.96))
    fig.savefig(path, dpi=110 if fast else 180)
    plt.close(fig)


def filter_policy_plot_rows(
    rows: list[dict[str, Any]],
    xlim: tuple[float, float] | None,
    *,
    wealth_key: str = "liquid_wealth",
    mass_filter_min: float | None = None,
) -> list[dict[str, Any]]:
    filtered = filter_rows_to_xlim(rows, xlim, wealth_key=wealth_key)
    if mass_filter_min is None:
        return filtered
    return [
        r
        for r in filtered
        if maybe_float(r.get("ergodic_mass")) > float(mass_filter_min)
        or maybe_float(r.get("mass_within_age_z_childless_renter")) > float(mass_filter_min)
    ]


def filter_rows_to_xlim(
    rows: list[dict[str, Any]],
    xlim: tuple[float, float] | None,
    *,
    wealth_key: str = "liquid_wealth",
) -> list[dict[str, Any]]:
    if xlim is None:
        return rows
    lo, hi = xlim
    return [r for r in rows if lo <= maybe_float(r.get(wealth_key)) <= hi]


def select_first_look_policy_rows(rows: list[dict[str, Any]], *, mode: str) -> list[dict[str, Any]]:
    if mode != "simple":
        return rows
    ages = sorted({float(r["age"]) for r in rows})
    z_values = sorted({float(r["z"]) for r in rows})
    keep_ages = set(ages if len(ages) <= 2 else [ages[0], ages[-1]])
    keep_z = set(z_values if len(z_values) <= 2 else [z_values[0], z_values[-1]])
    return [r for r in rows if float(r["age"]) in keep_ages and float(r["z"]) in keep_z]


def wealth_density_rows(sol: Any, P: Any, *, wealth_measure: str = "liquid") -> list[dict[str, Any]]:
    g = np.asarray(getattr(sol, "g", np.zeros(0)), dtype=float)
    b_grid = np.asarray(getattr(sol, "b_grid", np.arange(g.shape[0] if g.ndim else 0)), dtype=float).reshape(-1)
    if g.ndim != 7 or g.shape[0] != len(b_grid):
        return []
    total_mass = float(np.sum(g))
    asset_price = owner_asset_price_vector(sol, P)
    measure = "total_wealth_net_liquidation" if wealth_measure == "total" else "liquid_wealth"
    group_specs = [
        ("all", "all", "all", range(g.shape[1]), slice(None)),
        ("childless_renter", "childless", "renter", range(0, min(g.shape[1], 1)), slice(0, 1)),
        ("childless_owner", "childless", "owner", range(1, g.shape[1]), slice(0, 1)),
        ("parent_renter", "parent", "renter", range(0, min(g.shape[1], 1)), slice(1, g.shape[5])),
        ("parent_owner", "parent", "owner", range(1, g.shape[1]), slice(1, g.shape[5])),
    ]
    out: list[dict[str, Any]] = []
    for group, parent_status, tenure_status, tenure_indices, parity_slice in group_specs:
        mass_by_wealth: dict[float, float] = {}
        for ten in tenure_indices:
            for market in range(g.shape[2]):
                arr = g[:, ten, market, :, :, parity_slice, :]
                mass_by_b = np.sum(arr, axis=tuple(range(1, arr.ndim))) if arr.size else np.zeros_like(b_grid)
                housing_shift = (
                    liquidated_housing_value(P, asset_price, int(ten), int(market))
                    if wealth_measure == "total"
                    else 0.0
                )
                if not math.isfinite(housing_shift):
                    continue
                wealth_values = b_grid + housing_shift
                for wealth, mass_value in zip(wealth_values, mass_by_b):
                    wealth_key = float(np.round(float(wealth), 12))
                    mass_by_wealth[wealth_key] = mass_by_wealth.get(wealth_key, 0.0) + float(mass_value)
        group_mass = float(sum(mass_by_wealth.values()))
        for wealth, mass_value in sorted(mass_by_wealth.items()):
            mass = float(mass_value)
            out.append(
                {
                    "wealth_measure": measure,
                    "group": group,
                    "parent_status": parent_status,
                    "tenure_status": tenure_status,
                    "wealth": float(wealth),
                    "liquid_wealth": float(wealth) if measure == "liquid_wealth" else math.nan,
                    "total_wealth": float(wealth) if measure != "liquid_wealth" else math.nan,
                    "ergodic_mass": mass,
                    "population_share": float(mass / total_mass) if total_mass > 1e-14 else math.nan,
                    "within_group_share": float(mass / group_mass) if group_mass > 1e-14 else math.nan,
                    "group_population_share": float(group_mass / total_mass) if total_mass > 1e-14 else math.nan,
                    "total_population_mass": total_mass,
                    "liquidated_housing_discount_1_minus_psi": 1.0 - float(getattr(P, "psi", 0.0)),
                }
            )
    return out


def row_wealth(row: dict[str, Any]) -> float:
    wealth = maybe_float(row.get("wealth"))
    if math.isfinite(wealth):
        return wealth
    return maybe_float(row.get("liquid_wealth"))


def plot_first_look_wealth_density(
    rows: list[dict[str, Any]],
    path: Path,
    *,
    wealth_label: str,
    title: str,
    fast: bool = False,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not rows:
        return
    aggregate = sorted([r for r in rows if str(r.get("group")) == "all"], key=row_wealth)
    if not aggregate:
        return
    xs = np.asarray([row_wealth(r) for r in aggregate], dtype=float)
    ys = np.asarray([maybe_float(r["population_share"]) for r in aggregate], dtype=float)
    valid = np.isfinite(xs) & np.isfinite(ys)
    if not np.any(valid):
        return
    xs = xs[valid]
    ys = ys[valid]
    if xs.size > 1:
        diffs = np.diff(np.sort(xs))
        width = 0.8 * float(np.nanmin(diffs[diffs > 0])) if np.any(diffs > 0) else 0.8
    else:
        width = 0.8

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.6) if fast else (14.5, 4.4), sharex=True, sharey=True)
    ax_all, ax_tenure, ax_parent = axes

    ax_all.bar(xs, ys, width=width, color="0.35", alpha=0.82, align="center")
    ax_all.plot(xs, ys, color="0.05", lw=1.7)
    ax_all.set_title("Aggregate")
    ax_all.set_ylabel("population share")

    renter = density_series(rows, xs, ["childless_renter", "parent_renter"])
    owner = density_series(rows, xs, ["childless_owner", "parent_owner"])
    ax_tenure.plot(xs, renter, color="tab:blue", lw=2.0, label=f"renters ({100.0 * np.nansum(renter):.1f}%)")
    ax_tenure.plot(xs, owner, color="tab:orange", lw=2.0, label=f"owners ({100.0 * np.nansum(owner):.1f}%)")
    ax_tenure.set_title("By tenure")
    ax_tenure.legend(frameon=False, fontsize=8)

    childless = density_series(rows, xs, ["childless_renter", "childless_owner"])
    parents = density_series(rows, xs, ["parent_renter", "parent_owner"])
    ax_parent.plot(
        xs,
        childless,
        color="tab:purple",
        lw=2.0,
        label=f"childless ({100.0 * np.nansum(childless):.1f}%)",
    )
    ax_parent.plot(xs, parents, color="tab:green", lw=2.0, label=f"parents ({100.0 * np.nansum(parents):.1f}%)")
    ax_parent.set_title("By fertility status")
    ax_parent.legend(frameon=False, fontsize=8)

    for ax in axes:
        ax.set_xlabel(wealth_label)
        ax.grid(axis="y", alpha=0.2)
    y_max = max(float(np.nanmax(ys)) * 1.15, 1e-12)
    for arr in [renter, owner, childless, parents]:
        if arr.size:
            y_max = max(y_max, float(np.nanmax(arr)) * 1.15)
    for ax in axes:
        ax.set_ylim(0.0, y_max)
    xlim = density_plot_xlim(xs, ys)
    if xlim is not None:
        for ax in axes:
            ax.set_xlim(*xlim)
    fig.suptitle(title, y=0.98)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))
    fig.savefig(path, dpi=110 if fast else 180)
    plt.close(fig)


def density_series(rows: list[dict[str, Any]], xs: np.ndarray, groups: list[str]) -> np.ndarray:
    values = {float(x): 0.0 for x in xs}
    keep = set(groups)
    for row in rows:
        if str(row.get("group")) not in keep:
            continue
        wealth = row_wealth(row)
        share = maybe_float(row.get("population_share"))
        if math.isfinite(wealth) and math.isfinite(share):
            values[float(wealth)] = values.get(float(wealth), 0.0) + share
    return np.asarray([values.get(float(x), 0.0) for x in xs], dtype=float)


def density_xlim_from_rows(rows: list[dict[str, Any]]) -> tuple[float, float] | None:
    aggregate = sorted([r for r in rows if str(r.get("group")) == "all"], key=row_wealth)
    if not aggregate:
        return None
    xs = np.asarray([row_wealth(r) for r in aggregate], dtype=float)
    ys = np.asarray([maybe_float(r["population_share"]) for r in aggregate], dtype=float)
    return density_plot_xlim(xs, ys)


def policy_xlim_from_rows(
    rows: list[dict[str, Any]],
    *,
    wealth_key: str,
    base_xlim: tuple[float, float] | None = None,
) -> tuple[float, float] | None:
    vals: list[float] = []
    for row in rows:
        base_x = maybe_float(row.get("liquid_wealth"))
        if base_xlim is not None and not (base_xlim[0] <= base_x <= base_xlim[1]):
            continue
        x = maybe_float(row.get(wealth_key))
        if math.isfinite(x):
            vals.append(x)
    if not vals:
        return base_xlim
    lo = float(np.nanmin(vals))
    hi = float(np.nanmax(vals))
    pad = max(0.5, 0.06 * max(hi - lo, 1.0))
    return lo - pad, hi + pad


def density_plot_xlim(xs: np.ndarray, ys: np.ndarray) -> tuple[float, float] | None:
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    if not np.any(valid):
        return None
    x = x[valid]
    y = y[valid]
    threshold = 1e-3
    support = x[y > threshold]
    if support.size == 0:
        return None
    lo = float(np.nanmin(support))
    hi = float(np.nanmax(support))
    pad = max(0.5, 0.06 * max(hi - lo, 1.0))
    return lo - pad, hi + pad


def wealth_grid_coverage_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    aggregate = sorted([r for r in rows if str(r.get("group")) == "all"], key=row_wealth)
    if not aggregate:
        return []
    xs = np.asarray([row_wealth(r) for r in aggregate], dtype=float)
    ys = np.asarray([maybe_float(r["population_share"]) for r in aggregate], dtype=float)
    valid = np.isfinite(xs) & np.isfinite(ys)
    xs = xs[valid]
    ys = ys[valid]
    total_share = float(np.sum(ys))
    out: list[dict[str, Any]] = []
    n_grid = int(xs.size)
    grid_min = float(np.nanmin(xs)) if xs.size else math.nan
    grid_max = float(np.nanmax(xs)) if xs.size else math.nan
    wealth_measure = str(aggregate[0].get("wealth_measure", "liquid_wealth"))
    for threshold in [1e-7, 1e-6, 1e-5, 1e-4, 1e-3]:
        keep = ys > threshold
        support = xs[keep]
        out.append(
            {
                "wealth_measure": wealth_measure,
                "diagnostic": "per_gridpoint_mass_cutoff",
                "cutoff_population_share": float(threshold),
                "grid_points": n_grid,
                "grid_min": grid_min,
                "grid_max": grid_max,
                "occupied_grid_points": int(np.sum(keep)),
                "occupied_grid_share": float(np.sum(keep) / max(n_grid, 1)),
                "occupied_mass_share": float(np.sum(ys[keep]) / max(total_share, 1e-14)),
                "support_min": float(np.nanmin(support)) if support.size else math.nan,
                "support_max": float(np.nanmax(support)) if support.size else math.nan,
            }
        )
    intervals = [
        ("negative_b" if wealth_measure == "liquid_wealth" else "negative_total_wealth", -math.inf, 0.0),
        ("b_0_to_6" if wealth_measure == "liquid_wealth" else "total_wealth_0_to_6", 0.0, 6.0),
        ("b_6_to_10" if wealth_measure == "liquid_wealth" else "total_wealth_6_to_10", 6.0, 10.0),
        ("b_ge_10" if wealth_measure == "liquid_wealth" else "total_wealth_ge_10", 10.0, math.inf),
    ]
    for name, lo, hi in intervals:
        keep = (xs >= lo) & (xs < hi)
        occupied = keep & (ys > 0.0)
        out.append(
            {
                "wealth_measure": wealth_measure,
                "diagnostic": name,
                "cutoff_population_share": math.nan,
                "grid_points": n_grid,
                "grid_min": grid_min,
                "grid_max": grid_max,
                "occupied_grid_points": int(np.sum(keep)),
                "occupied_grid_share": float(np.sum(keep) / max(n_grid, 1)),
                "occupied_mass_share": float(np.sum(ys[keep]) / max(total_share, 1e-14)),
                "support_min": float(np.nanmin(xs[occupied])) if np.any(occupied) else math.nan,
                "support_max": float(np.nanmax(xs[occupied])) if np.any(occupied) else math.nan,
            }
        )
    return out


def plot_housing_by_tenure(
    ax: Any,
    rows: list[dict[str, Any]],
    z_values: list[float],
    ages: list[float],
    color_for_z: dict[float, Any],
    style_for_age: dict[float, str],
    xlim: tuple[float, float] | None = None,
    wealth_key: str = "liquid_wealth",
    wealth_label: str = "liquid wealth",
    mass_filter_min: float | None = None,
) -> None:
    for z in z_values:
        for age in ages:
            panel = [r for r in rows if float(r["z"]) == z and float(r["age"]) == age]
            panel.sort(key=lambda r: maybe_float(r.get(wealth_key)))
            panel = filter_policy_plot_rows(
                panel,
                xlim,
                wealth_key=wealth_key,
                mass_filter_min=mass_filter_min,
            )
            if not panel:
                continue
            xs = [maybe_float(r.get(wealth_key)) for r in panel]
            ys = [maybe_float(r["housing_services_after_tenure"]) for r in panel]
            ax.plot(xs, ys, color=color_for_z[z], linestyle=style_for_age[age], lw=1.8, alpha=0.9)
            owner_points = [r for r in panel if str(r.get("chosen_tenure")) == "owner"]
            if owner_points:
                ax.scatter(
                    [maybe_float(r.get(wealth_key)) for r in owner_points],
                    [maybe_float(r["housing_services_after_tenure"]) for r in owner_points],
                    color=color_for_z[z],
                    s=10,
                    alpha=0.9,
                )
    ax.set_title("Expected housing by tenure choice")
    ax.set_xlabel(wealth_label)
    ax.set_ylabel("housing services")
    ax.grid(alpha=0.2)
    ax.text(
        0.02,
        0.98,
        "line: expected over tenure shocks\nmarker: modal owner branch",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "0.85", "alpha": 0.85, "pad": 3},
    )


def plot_market_panel(ax: Any, rows: list[dict[str, Any]]) -> None:
    if not rows:
        ax.axis("off")
        return
    rows = sorted(rows, key=lambda r: int(r["market"]))
    markets = np.arange(len(rows))
    rental = np.asarray([maybe_float(r["rental_demand"]) for r in rows], dtype=float)
    owner = np.asarray([maybe_float(r["owner_demand"]) for r in rows], dtype=float)
    total = rental + owner
    supply = np.asarray([maybe_float(r["supply"]) for r in rows], dtype=float)
    residual = (total - supply) / np.maximum(supply, 1e-12)
    max_residual = float(np.nanmax(np.abs(residual))) if residual.size else math.nan
    width = 0.32
    ax.bar(markets - width / 2, rental, width, label="renter demand", color="tab:blue", alpha=0.85)
    ax.bar(markets + width / 2, owner, width, label="owner demand", color="tab:orange", alpha=0.85)
    if math.isfinite(max_residual) and max_residual > 1e-3:
        ax.scatter(markets, supply, marker="x", s=70, color="0.15", label="supply")
    ax.set_xticks(markets, [str(int(r["market"])) for r in rows])
    ax.set_xlabel("housing-services market")
    ax.set_ylabel("quantity demanded")
    ax.set_title(f"Equilibrium prices, rents, and quantities (max residual={max_residual:.1e})")
    ax.grid(axis="y", alpha=0.2)

    ax2 = ax.twinx()
    asset = np.asarray([maybe_float(r.get("house_price", r["asset_price"])) for r in rows], dtype=float)
    rent = np.asarray([maybe_float(r["flow_rent_or_user_cost"]) for r in rows], dtype=float)
    ax2.plot(markets, asset, color="tab:red", marker="o", lw=1.8, label="house price")
    ax2.plot(markets, rent, color="tab:green", marker="s", lw=1.8, label="flow rent/user cost")
    ax2.set_ylabel("price")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8, loc="upper left")


def plot_room_bins(rows: list[dict[str, Any]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0), sharey=True)
    x = np.arange(len(ROOM_BIN_ORDER))
    width = 0.36
    for ax, tenure in zip(axes, TENURE_ORDER):
        panel = [r for r in rows if r["tenure"] == tenure]
        data = [float(next(r["data_share"] for r in panel if r["bin"] == b)) for b in ROOM_BIN_ORDER]
        model = [float(next(r["model_share"] for r in panel if r["bin"] == b)) for b in ROOM_BIN_ORDER]
        ax.bar(x - width / 2, data, width, label="data")
        ax.bar(x + width / 2, model, width, label="model")
        ax.set_xticks(x, ROOM_BIN_ORDER)
        ax.set_ylim(0.0, 1.0)
        ax.set_title(f"{tenure} room-bin shares")
        ax.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("share within tenure")
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_owner_rungs(rows: list[dict[str, Any]], path: Path, *, ylabel: str = "share of owners") -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    ax.bar(x, [float(r["share"]) for r in rows])
    ax.set_xticks(x, [f"{float(r['rooms']):g}" for r in rows])
    ax.set_xlabel("owner rung rooms")
    ax.set_ylabel(ylabel)
    ax.set_ylim(0.0, 1.0)
    ax.set_title("Owner rung shares")
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_age_profiles(rows: list[dict[str, Any]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ages = np.asarray([float(r["age"]) for r in rows])
    specs = [
        ("ownership_rate", "ownership rate"),
        ("fertility", "expected births"),
        ("mean_housing_services", "housing services"),
        ("mean_liquid_wealth", "liquid wealth"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.0), sharex=True)
    for ax, (key, ylabel) in zip(axes.ravel(), specs):
        ax.plot(ages, [float(r[key]) for r in rows], lw=2.0)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.2)
    for ax in axes[-1, :]:
        ax.set_xlabel("age")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_tenure_by_age(rows: list[dict[str, Any]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ages = np.asarray([float(r["age"]) for r in rows])
    owner = np.asarray([float(r["owner_share"]) for r in rows])
    renter = np.asarray([float(r["renter_share"]) for r in rows])
    young_mask = (ages >= 25.0) & (ages <= 34.0)
    old_mask = (ages >= 65.0) & (ages <= 75.0)
    young_owner = float(np.mean(owner[young_mask])) if np.any(young_mask) else math.nan
    old_owner = float(np.mean(owner[old_mask])) if np.any(old_mask) else math.nan

    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    ax.axvspan(25.0, 34.0, color="0.92", label="young target ages 25-34")
    ax.plot(ages, owner, marker="o", lw=2.2, color="#4C78A8", label="owners")
    ax.plot(ages, renter, marker="s", lw=2.0, color="#F58518", label="renters")
    ax.set_xlabel("age")
    ax.set_ylabel("population share")
    ax.set_ylim(0.0, 1.03)
    ax.set_title(f"Tenure by age: young owner share={young_owner:.3f}, old owner share={old_owner:.3f}")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, loc="center right")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_owner_entry_thresholds(rows: list[dict[str, Any]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    z_values = sorted({float(r["z"]) for r in rows})
    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    for z in z_values:
        panel = [r for r in rows if float(r["z"]) == z]
        panel.sort(key=lambda r: float(r["age"]))
        ax.plot(
            [float(r["age"]) for r in panel],
            [maybe_float(r["threshold_prob_0p50"]) for r in panel],
            lw=1.8,
            marker="o",
            ms=3,
            label=f"z={z:g}",
        )
    ax.set_xlabel("age")
    ax.set_ylabel("liquid wealth threshold for owner prob. >= 0.5")
    ax.set_title("Owner-entry thresholds for childless renters")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_policy_moment_comparison(records: list[dict[str, Any]], targets: dict[str, float], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    keys = [
        "tfr",
        "childless_rate",
        "own_rate",
        "own_rate_2534",
        "old_age_own_rate",
        "prime30_55_childless_owner_minus_renter_mean_rooms",
        "prime30_55_childless_renter_mean_rooms",
        "old_parent_childless_nonhousing_wealth_to_income_gap_6575",
    ]
    labels = [r["case"] for r in records]
    x = np.arange(len(records))
    fig, axes = plt.subplots(4, 2, figsize=(12.0, 10.5))
    for ax, key in zip(axes.ravel(), keys):
        ax.bar(x, [float(r["moments"].get(key, math.nan)) for r in records])
        if key in targets:
            ax.axhline(float(targets[key]), color="0.25", linestyle="--", lw=1.0)
        ax.set_title(key)
        ax.set_xticks(x, labels, rotation=30, ha="right", fontsize=8)
        ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_contact_sheet(outdir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    candidates = [
        outdir / "ergodic_deterministic_policy_top_slice.png",
        outdir / "ergodic_deterministic_policy_slices.png",
        outdir / "ergodic_deterministic_policy_bins.png",
        outdir / "first_look_policies_markets_on_path.png",
        outdir / "first_look_policies_markets.png",
        outdir / "first_look_wealth_density.png",
        outdir / "first_look_policies_markets_total_wealth_on_path.png",
        outdir / "first_look_policies_markets_total_wealth.png",
        outdir / "first_look_total_wealth_density.png",
        outdir / "first_look_policies_markets_full.png",
        outdir / "first_look_policies_markets_total_wealth_full.png",
        outdir / "tenure_by_age.png",
        outdir / "diagnostics/ownership_by_age.png",
        outdir / "diagnostics/policy_childless_renter_age30.png",
        outdir / "diagnostics/policy_childless_renter_age42.png",
        outdir / "room_bin_shares_prime30_55_childless.png",
        outdir / "owner_rung_shares_all_owners.png",
        outdir / "owner_rung_shares_prime30_55_childless.png",
        outdir / "owner_entry_thresholds.png",
    ]
    paths = [p for p in candidates if p.exists()]
    if not paths:
        return
    ncols = 2
    nrows = int(math.ceil(len(paths) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(11.0, 4.4 * nrows))
    axes_arr = np.asarray(axes).reshape(-1)
    for ax, img_path in zip(axes_arr, paths):
        img = plt.imread(img_path)
        ax.imshow(img)
        ax.set_title(img_path.name)
        ax.axis("off")
    for ax in axes_arr[len(paths) :]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(outdir / "contact_sheet.png", dpi=180)
    plt.close(fig)


def write_readme(
    outdir: Path,
    source: dict[str, Any],
    baseline: dict[str, Any],
    target_set: str,
    targets: dict[str, float],
    policy_records: list[dict[str, Any]],
) -> None:
    moments = baseline["moments"]
    lines = [
        "# Intergen Mechanics Packet",
        "",
        "Diagnostic mechanics packet for the June 2026 one-market intergenerational model. This is not a production calibration or quantitative policy estimate.",
        "",
        f"Source: `{source['source_meta'].get('path')}`.",
        f"Target set: `{target_set}`.",
        "",
        "## Headline Fit",
        "",
        "| Moment | Model | Target |",
        "|---|---:|---:|",
    ]
    for key in [
        "tfr",
        "childless_rate",
        "own_rate",
        "own_rate_2534",
        "old_age_own_rate",
        "prime30_55_childless_renter_mean_rooms",
        "prime30_55_childless_owner_mean_rooms",
        "prime30_55_childless_owner_minus_renter_mean_rooms",
    ]:
        target = targets.get(key, math.nan)
        lines.append(f"| `{key}` | {fmt(moments.get(key))} | {fmt(target)} |")
    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `diagnostics/`: standard plots from `intergen_housing_fertility.diagnostics.write_diagnostics`.",
            "- `first_look_policies_markets_on_path.png`: simpler 2-by-2 inspection panel using tenure-probability expected policies and dropping zero-mass points in each displayed childless-renter age/income slice.",
            "- `first_look_policies_markets.png`: same simple panel on the displayed wealth support without the on-path mass filter; use it for numerical grid inspection, not first-pass economics.",
            "- `first_look_policies_markets_full.png`: fuller version with all selected income states and ages.",
            "- `first_look_policies_markets_total_wealth_on_path.png`: on-path first-look panel against expected total wealth after tenure choice.",
            "- `first_look_policies_markets_total_wealth.png` and `_full.png`: same first-look policy panels with the policy x-axis set to expected liquid wealth plus net liquidated housing value after tenure choice.",
            "- `first_look_wealth_density.png` and `.csv`: aggregate ergodic mass over liquid wealth `b`, plus non-stacked renter/owner and childless/parent density panels.",
            "- `first_look_total_wealth_density.png` and `.csv`: same density diagnostics over total wealth `b + (1 - psi) p H` for owners and `b` for renters.",
            "- `wealth_grid_coverage.csv`: how much of the liquid-wealth grid is economically occupied under several mass thresholds.",
            "- `total_wealth_grid_coverage.csv`: occupied support under the total-wealth definition.",
            "- `first_look_policy_lines.csv` and `first_look_market_summary.csv`: source data for the first-look panels, including chosen tenure, ergodic mass, and house-price/user-cost accounting.",
            "- `ergodic_deterministic_policy_states.csv`: deterministic argmax policy rows for every occupied state in the stationary distribution.",
            "- `ergodic_deterministic_policy_top_slice.png`: one literal deterministic policy function for a common occupied childless-renter state.",
            "- `ergodic_deterministic_policy_slices.png`: readable deterministic policy slices for the most common occupied childless-renter age/income cells.",
            "- `ergodic_deterministic_policy_bins.png`: mass-weighted deterministic policies over all occupied states by liquid wealth.",
            "- `target_fit.csv` and `target_fit.md`: full target/model/gap table with loss contributions.",
            "- `room_bin_fit_prime30_55_childless.csv` and `room_bin_shares_prime30_55_childless.png`: model versus ACS room-bin shares.",
            "- `owner_rung_shares_all_owners.csv` and `.png`: realized owner mass across the full owner room ladder.",
            "- `owner_rung_shares_prime30_55_childless.csv` and `.png`: owner rung concentration among prime-age childless owners.",
            "- `age_profiles.csv` and `.png`: lifecycle ownership, fertility, housing, and liquid wealth profiles.",
            "- `tenure_by_age.csv` and `.png`: simple renter/owner population shares by age.",
            "- `owner_entry_thresholds.csv` and `.png`: childless-renter owner-entry probability thresholds by age and income state.",
            "- `owner_entry_policy_childless_renter_age30_42.csv`: owner-entry probability lines used for age-30 and age-42 inspection.",
            "- `solution_cache.pkl`: local trusted Python pickle with the solved `sol`, `P`, and `p_eq` objects. Use `--refresh-plots-from-cache` to rebuild figures and CSVs from it without re-solving.",
            "- `contact_sheet.png`: quick visual index when standard diagnostics are enabled.",
        ]
    )
    if policy_records:
        lines.extend(
            [
                "- `policy_cases/`: diagnostic policy proof-of-concept cases from the same theta.",
                "",
                "## Policy Caveat",
                "",
                "Policy cases are fixed-theta diagnostic mechanics only. The estate-tax case is a terminal bequest wedge without revenue rebates or inheritance kernels.",
            ]
        )
    lines.append("")
    (outdir / "README.md").write_text("\n".join(lines))


def write_markdown_table(path: Path, rows: list[dict[str, Any]], title: str) -> None:
    lines = [
        f"# {title}",
        "",
        "| Moment | Target | Model | Gap | Weight | Contribution |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| `{row['moment']}` | {fmt(row['target'])} | {fmt(row['model'])} | "
            f"{fmt(row['gap_model_minus_target'])} | {fmt(row['weight'])} | {fmt(row['loss_contribution'])} |"
        )
    lines.append("")
    path.write_text("\n".join(lines))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    flat = [flatten_row(row) for row in rows]
    fields: list[str] = []
    for row in flat:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(flat)


def flatten_row(row: dict[str, Any], prefix: str = "") -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in row.items():
        name = f"{prefix}{key}"
        if isinstance(value, dict):
            out.update(flatten_row(value, f"{name}."))
        elif isinstance(value, (list, tuple, np.ndarray)):
            out[name] = json.dumps(jsonable(value))
        else:
            out[name] = value
    return out


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(jsonable(obj), indent=2, sort_keys=True))


def maybe_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def safe_nanmin(values: np.ndarray) -> float:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    return float(np.min(vals)) if vals.size else math.nan


def safe_nanmax(values: np.ndarray) -> float:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    return float(np.max(vals)) if vals.size else math.nan


def fmt(value: Any) -> str:
    val = maybe_float(value)
    if not math.isfinite(val):
        return ""
    return f"{val:.3f}"


if __name__ == "__main__":
    main()
