"""Command line interface for the intergenerational housing/fertility model."""

from __future__ import annotations

import argparse
from pathlib import Path

from .diagnostics import write_diagnostics
from .solver import solve_model
from .utils import dumps_json


def main() -> None:
    parser = argparse.ArgumentParser(description="No-geography intergenerational housing/fertility model")
    sub = parser.add_subparsers(dest="cmd", required=True)

    smoke = sub.add_parser("smoke", help="Run a small fixed-price smoke solve")
    smoke.add_argument("--quiet", action="store_true")
    smoke.add_argument("--json", type=Path, default=None)

    ge = sub.add_parser("solve", help="Run owner-size price iteration")
    ge.add_argument("--mode", choices=["smoke", "benchmark"], default="smoke")
    ge.add_argument("--max-iter-eq", type=int, default=None)
    ge.add_argument("--quiet", action="store_true")
    ge.add_argument("--json", type=Path, default=None)

    diag = sub.add_parser("diagnostics", help="Solve and write a diagnostic packet")
    diag.add_argument("--mode", choices=["smoke", "benchmark"], default="smoke")
    diag.add_argument("--max-iter-eq", type=int, default=20)
    diag.add_argument("--fixed-prices", action="store_true")
    diag.add_argument("--outdir", type=Path, default=Path("../../output/model/intergen_housing_fertility_smoke"))
    diag.add_argument("--quiet", action="store_true")

    args = parser.parse_args()
    if args.cmd == "smoke":
        sol, P = solve_model(mode="smoke", solve_prices=False, verbose=not args.quiet)
        report = summarize(sol, P, solve_prices=False)
    elif args.cmd == "solve":
        overrides = {}
        if args.max_iter_eq is not None:
            overrides["max_iter_eq"] = args.max_iter_eq
        sol, P = solve_model(overrides, mode=args.mode, solve_prices=True, verbose=not args.quiet)
        report = summarize(sol, P, solve_prices=True)
    elif args.cmd == "diagnostics":
        overrides = {"max_iter_eq": args.max_iter_eq}
        sol, P = solve_model(overrides, mode=args.mode, solve_prices=not args.fixed_prices, verbose=not args.quiet)
        write_diagnostics(sol, P, args.outdir)
        report = summarize(sol, P, solve_prices=not args.fixed_prices)
        report["outdir"] = str(args.outdir)
    else:
        raise ValueError(args.cmd)

    text = dumps_json(report)
    print(text)
    if getattr(args, "json", None) is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text)


def summarize(sol, P, *, solve_prices: bool) -> dict:
    return {
        "mode": P.mode,
        "solve_prices": solve_prices,
        "converged": bool(getattr(sol, "converged", True)),
        "elapsed_sec": float(getattr(sol, "elapsed_sec", 0.0)),
        "owner_user_cost": sol.owner_user_cost,
        "owner_asset_price": sol.owner_asset_price,
        "owner_demand_by_size": sol.owner_demand_by_size,
        "owner_supply": sol.owner_supply,
        "owner_excess_by_size": sol.owner_excess_by_size,
        "aggregate_owner_demand": sol.aggregate_owner_demand,
        "aggregate_owner_supply": sol.aggregate_owner_supply,
        "aggregate_owner_excess": sol.aggregate_owner_excess,
        "own_rate": sol.own_rate,
        "young_owner_rate": sol.young_owner_rate,
        "old_owner_rate": sol.old_owner_rate,
        "mean_completed_fertility": sol.mean_completed_fertility,
        "childless_rate": sol.childless_rate,
        "price_iterations": len(getattr(sol, "price_trace", [])),
        "best_market_metric": float(getattr(sol, "best_market_metric", 0.0)),
        "best_max_abs_rel_excess": float(getattr(sol, "best_max_abs_rel_excess", 0.0)),
        "last_price_trace": getattr(sol, "price_trace", [])[-5:],
    }


if __name__ == "__main__":
    main()
