#!/usr/bin/env python3
"""Reduced-grid smoke run for the competitive housing-product market branch."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

MODEL_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[3]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from dt_cp_model.parameters import asdict, build_calibration_setup  # noqa: E402
from dt_cp_model.solver import run_model_cp_dt  # noqa: E402


def main() -> None:
    args = parse_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    branch = git_text(["git", "branch", "--show-current"])
    commit = git_text(["git", "rev-parse", "HEAD"])
    dirty = bool(git_text(["git", "status", "--porcelain"]))

    setup = build_calibration_setup(args.setup)
    P = SimpleNamespace(**asdict(setup.P_base))
    apply_smoke_overrides(P, args)
    P.housing_product_market = True

    t0 = time.perf_counter()
    sol, P_out, Lambda = run_model_cp_dt(P, verbose=not args.quiet)
    elapsed = time.perf_counter() - t0

    old_ok, old_error = run_old_path_smoke(setup, args)

    price_records = sol.product_price_table
    demand_records = sol.product_demand_table
    capacity_records = sol.capacity_table
    moment_records = build_moment_records(sol, Lambda, elapsed, old_ok)
    room_records = (
        [{"tenure": "R", **row} for row in getattr(sol, "prime_age_renter_room_distribution", [])]
        + [{"tenure": "O", **row} for row in getattr(sol, "prime_age_owner_room_distribution", [])]
    )

    write_csv(out_dir / "product_price_table.csv", price_records)
    write_csv(out_dir / "product_demand_table.csv", demand_records)
    write_csv(out_dir / "capacity_residuals.csv", capacity_records)
    write_csv(out_dir / "moment_summary.csv", moment_records)
    write_csv(out_dir / "prime_age_room_distributions.csv", room_records)
    (out_dir / "solver_trace.json").write_text(json.dumps(getattr(sol, "product_market_trace", []), indent=2))

    summary = build_summary(
        branch=branch,
        commit=commit,
        dirty=dirty,
        setup=args.setup,
        args=args,
        sol=sol,
        Lambda=Lambda,
        elapsed=elapsed,
        old_ok=old_ok,
        old_error=old_error,
        demand_records=demand_records,
        room_records=room_records,
    )
    (out_dir / "summary.md").write_text(summary)

    print(summary)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Smoke-test the housing-product market branch.")
    parser.add_argument("--setup", choices=["fast", "benchmark"], default="benchmark")
    parser.add_argument("--J", type=int, default=None)
    parser.add_argument("--Nb", type=int, default=None)
    parser.add_argument("--capacity-max-iter", type=int, default=None)
    parser.add_argument("--capacity-tol", type=float, default=None)
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "output/model/housing_product_market_smoke")
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args()


def apply_smoke_overrides(P: SimpleNamespace, args: argparse.Namespace) -> None:
    if args.J is not None:
        P.J = int(args.J)
        P.J_R = max(2, P.J - 1)
        P.A_f_end = min(max(2, P.J - 2), getattr(P, "A_f_end", P.J - 2))
    if args.Nb is not None:
        P.Nb = int(args.Nb)
    if args.capacity_max_iter is not None:
        P.capacity_max_iter = int(args.capacity_max_iter)
    if args.capacity_tol is not None:
        P.capacity_tol = float(args.capacity_tol)
    P.force_full_bellman = True


def run_old_path_smoke(setup: SimpleNamespace, args: argparse.Namespace) -> tuple[bool, str | None]:
    P = SimpleNamespace(**asdict(setup.P_base))
    P.housing_product_market = False
    if args.J is not None:
        P.J = int(args.J)
        P.J_R = max(2, P.J - 1)
        P.A_f_end = min(max(2, P.J - 2), getattr(P, "A_f_end", P.J - 2))
    if args.Nb is not None:
        P.Nb = int(args.Nb)
    P.max_iter_eq = 1
    P.solve_mode = "pe"
    P.force_full_bellman = True
    P.p_fixed = np.array([1.0, 1.2])
    P.w_fixed = np.ones(2)
    P.entry_shares_fixed = np.array([0.55, 0.45])
    try:
        run_model_cp_dt(P, verbose=False)
    except Exception as exc:  # pragma: no cover - smoke diagnostic
        return False, repr(exc)
    return True, None


def build_moment_records(sol: SimpleNamespace, Lambda: np.ndarray, elapsed: float, old_ok: bool) -> list[dict[str, Any]]:
    pairs = {
        "Lambda_P": float(Lambda[0]),
        "Lambda_C": float(Lambda[1]),
        "max_capacity_residual": float(getattr(sol, "max_capacity_residual", np.nan)),
        "product_market_converged": float(bool(getattr(sol, "product_market_converged", False))),
        "old_path_smoke_ok": float(old_ok),
        "tfr": float(2 * sol.mean_parity),
        "childless_rate": float(sol.parity_dist[0]),
        "own_rate": float(sol.own_rate_3055),
        "own_gradient": float(sol.own_gradient_3055),
        "own_family_gap": float(sol.own_gap_newparent_nonparent_3055),
        "prime_childless_renter_median_rooms": float(sol.prime_childless_renter_median_rooms),
        "prime_childless_owner_median_rooms": float(sol.prime_childless_owner_median_rooms),
        "H01": float(getattr(sol, "housing_increment_0to1_eventstudy_t3", np.nan)),
        "H12": float(getattr(sol, "housing_increment_1to2_proxy_t3", np.nan)),
        "center_share_nonparents": float(sol.center_share_nonparents_2245),
        "center_share_newparents": float(sol.center_share_newparents_2245),
        "elapsed_sec": float(elapsed),
    }
    return [{"moment": key, "value": value} for key, value in pairs.items()]


def build_summary(
    *,
    branch: str,
    commit: str,
    dirty: bool,
    setup: str,
    args: argparse.Namespace,
    sol: SimpleNamespace,
    Lambda: np.ndarray,
    elapsed: float,
    old_ok: bool,
    old_error: str | None,
    demand_records: list[dict[str, Any]],
    room_records: list[dict[str, Any]],
) -> str:
    demand_lines = [
        "| loc | h | tenure | rooms | demand | capacity use |",
        "|---|---:|---|---:|---:|---:|",
    ]
    for row in demand_records:
        if abs(float(row["demand"])) < 1e-12:
            continue
        demand_lines.append(
            f"| {row['location']} | {row['h']} | {row['tenure']} | {float(row['size_h']):.2f} | "
            f"{float(row['demand']):.6f} | {float(row['capacity_use']):.6f} |"
        )
    room_lines = [
        "| tenure | rooms | share | mass |",
        "|---|---:|---:|---:|",
    ]
    for row in room_records:
        room_lines.append(
            f"| {row['tenure']} | {float(row['size_h']):.2f} | {float(row['share']):.6f} | {float(row['mass']):.6f} |"
        )
    checks = getattr(sol, "product_market_checks", {})
    check_lines = [f"- `{key}`: `{value}`" for key, value in checks.items()]
    old_line = "yes" if old_ok else f"no: `{old_error}`"
    dirty_note = "dirty worktree" if dirty else "clean worktree"
    return "\n".join(
        [
            "# Housing Product Market Smoke",
            "",
            f"- Branch: `{branch}`",
            f"- Commit: `{commit}` ({dirty_note})",
            f"- Setup loaded: `{setup}` with `J={args.J if args.J is not None else 'benchmark'}`, `Nb={args.Nb if args.Nb is not None else 'benchmark'}`",
            f"- Old benchmark path still runs: {old_line}",
            f"- Housing-product-market path converged: `{bool(getattr(sol, 'product_market_converged', False))}`",
            f"- Lambda_P: `{float(Lambda[0]):.6f}`",
            f"- Lambda_C: `{float(Lambda[1]):.6f}`",
            f"- Max capacity residual: `{float(getattr(sol, 'max_capacity_residual', np.nan)):.6g}`",
            f"- Elapsed seconds: `{elapsed:.2f}`",
            "",
            "## Product Demand",
            "",
            *demand_lines,
            "",
            "## Prime-Age Room Distributions",
            "",
            *room_lines,
            "",
            "## Sanity Checks",
            "",
            *check_lines,
            "",
            "## Known Issues",
            "",
            "- This smoke validates mechanics and is not a calibrated benchmark result.",
            "- The product-market Bellman uses compiled kernels; the forward distribution is still a first-pass Python implementation and should be profiled before long calibration sweeps.",
        ]
    )


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def git_text(cmd: list[str]) -> str:
    try:
        out = subprocess.check_output(cmd, cwd=REPO_ROOT, text=True, stderr=subprocess.DEVNULL)
    except Exception:
        return "unknown"
    return out.strip()


if __name__ == "__main__":
    main()
