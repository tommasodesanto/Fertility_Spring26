"""Single normal-path entry point (local one core or Torch). Results are PRELIMINARY.

    python -m experiments.stationary_single_market.run fixed-price --credit reference --reference-root R --bundle B --bundle-sha S --out O
    python -m experiments.stationary_single_market.run fixed-price --credit corrected --d-bar 0 ...      (feasibility diagnostic)
    python -m experiments.stationary_single_market.run equilibrium --credit reference --initial-price-factor 1.05 ...

Only `experiments.stationary_single_market.engine` is imported: no frozen tools, no monkey-patching,
no calibration or psi normalization. A successful run is a preliminary
solve; certification is `verification/acceptance_oracle.py` (frozen reference gates,
14/31 tables, the 17 dated standard plots). Corrected credit is limited to a
fixed-price feasibility diagnostic: with d_bar=0 the known age-18 entrant
cells are Bellman-dead and the run stops with the engine's census. There is
no corrected GE until the author decides the entry contract.

In-process timers only; external wall time comes from verification/budget_run.py.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import json
import math
import os
import pickle
import time
from pathlib import Path

import numpy as np

from . import LABEL, credit, inputs
from .contract import validate_contract

STATUS = "preliminary_solve_only_not_certified"


def jsonable(value):
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist() if value.size <= 64 else dict(shape=list(value.shape), dtype=str(value.dtype))
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return repr(value)
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    return repr(value)


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(jsonable(value), indent=1, sort_keys=True) + "\n")


def prepare(args):
    loaded = inputs.load_inputs(args.bundle, args.reference_root, args.bundle_sha)
    P = copy.deepcopy(loaded.parameters)
    for flag in ("native_due_stayer_credit", "native_exact_inherited_distribution", "native_fixed_reference_entry",
                 "native_explicit_transaction_grid", "native_purchase_income", "native_exact_allocation_output"):
        if getattr(P, flag) is not True:
            raise RuntimeError("Reference contract flag not set in checkpoint: " + flag)
    validate_contract(P)
    credit.bind_engine_credit(P, args.credit, args.d_bar)
    # I/O location only (checkpoint value names the old run directory).
    P.native_inherited_distribution_evidence_dir = str(args.out / "inherited_state_failures")
    return loaded, P


def feasibility_report(exc, out: Path, d_bar) -> None:
    """Write the engine's own dead-mass census; nothing is repaired."""
    write_json(out / "feasibility_report.json", dict(
        status="infeasible_corrected_credit", d_bar=d_bar, stage=exc.stage, dead_mass=exc.dead_mass,
        census=exc.census, message=str(exc),
        action="Author decision required on entrant debt under the corrected contract; "
               "no asset truncation, mass deletion, transfer, positive d_bar or recalibration applied."))


def save_stage(out: Path, sol, sd, P, grid, price) -> dict:
    """Portable arrays (solution + lab precompute_shared as 'shared.*') and the
    trusted acceptance pickle (types.SimpleNamespace, NumPy, experiments.stationary_single_market.engine.*)."""
    out.mkdir(parents=True, exist_ok=False)
    arrays = {k: v for k, v in vars(sol).items() if isinstance(v, np.ndarray) and v.dtype != object}
    arrays.update({"shared." + k: v for k, v in vars(sd).items() if isinstance(v, np.ndarray) and v.dtype != object})
    np.savez(out / "solution_arrays.npz", **arrays)
    with gzip.open(out / "acceptance_solution.pkl.gz", "wb", compresslevel=1) as stream:
        pickle.dump(dict(solution=sol, shared=sd, parameters=P, b_grid=grid, price=np.asarray(price)), stream, protocol=5)
    return dict(arrays=len(arrays), path=str(out))


def _platform() -> dict:
    import platform
    import numba
    cache = os.environ.get("NUMBA_CACHE_DIR")
    files = sum(len(f) for _, _, f in os.walk(cache)) if cache and os.path.isdir(cache) else 0
    return dict(machine=platform.machine(), system=platform.platform(), python=platform.python_version(),
                numpy=np.__version__, numba=numba.__version__, numba_threads=numba.get_num_threads(),
                numba_cache_dir=cache, numba_cache_files_at_start=files,
                cache_state="cold" if files == 0 else "warm (pre-existing cache files)")


def ge_counters(timings: dict) -> dict:
    refine = timings.get("scalar_market_refine") or {}
    return dict(unique_fast_lifecycle_solves=timings.get("unique_fast_price_evaluations"),
                refinement_price_evaluations=refine.get("price_evaluations"),
                refinement_expansions=refine.get("expansions"),
                refinement_seconds=refine.get("price_evaluation_time_sec"),
                price_cache_hits=timings.get("price_cache_hits"),
                fast_solve_seconds_initial=timings.get("markov_income_solve_time"),
                final_full_solution="see profiled.final_solution (requires --count-calls)",
                note="iterations_completed is the initial-evaluation count, not total solves")


def main() -> None:
    t_main = time.perf_counter()
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("check-inputs", "fixed-price", "equilibrium"))
    ap.add_argument("--reference-root", type=Path, required=True)
    ap.add_argument("--bundle", type=Path, required=True)
    ap.add_argument("--bundle-sha", required=True)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--credit", choices=("reference", "corrected"))
    ap.add_argument("--d-bar", type=float)
    ap.add_argument("--repeat", type=int, default=1)
    ap.add_argument("--initial-price-factor", type=float)
    ap.add_argument("--skip-plots", action="store_true")
    ap.add_argument("--max-lifecycle", type=int, help="hard budget: at-price call N+1 is never started")
    ap.add_argument("--solve-deadline", type=float, help="hard solve-stage seconds (excludes import/output)")
    ap.add_argument("--count-calls", action="store_true",
                    help="verification-only sys.setprofile counter (use for both engines or neither)")
    args = ap.parse_args()
    if args.command == "check-inputs":
        print(json.dumps(inputs.load_inputs(args.bundle, args.reference_root, args.bundle_sha).identity, indent=1))
        return
    if args.credit is None or args.out is None:
        raise SystemExit("--credit and --out are required; there are no defaults")
    args.out.mkdir(parents=True, exist_ok=False)
    timing, counts = {}, {}
    t = time.perf_counter()
    loaded, P = prepare(args)
    manifest = inputs.load_manifest(args.reference_root)
    timing["inputs"] = time.perf_counter() - t
    t = time.perf_counter()
    from .engine import diagnostics, e5f_stationary_paygo as paygo, solver
    import numba
    timing["import"] = time.perf_counter() - t
    grid = loaded.b_grid.copy()
    from .engine.utils import make_grid
    np.testing.assert_array_equal(make_grid(P), grid)
    receipt = dict(status=STATUS, label=LABEL, command=args.command, credit=args.credit, d_bar=args.d_bar,
                   identity=loaded.identity, numba_threads=int(numba.get_num_threads()),
                   omp=os.environ.get("OMP_NUM_THREADS"), platform=_platform())
    if args.command == "fixed-price":
        price = loaded.reference_price.copy()
        solves, stages = [], []
        for k in range(args.repeat):
            Pk = copy.deepcopy(P)
            t = time.perf_counter()
            sd = solver.precompute_shared(Pk, grid)
            t_sd = time.perf_counter() - t
            try:
                sol = solver.solve_markov_income_at_prices(price, Pk, grid, SD=sd, verbose=False, fast_stats=False)
            except solver.InfeasibleThetaError as exc:
                if args.credit != "corrected":
                    raise
                feasibility_report(exc, args.out, args.d_bar)
                raise SystemExit(3)
            t_solve = time.perf_counter() - t - t_sd
            t = time.perf_counter()
            stages.append(save_stage(args.out / f"rep{k + 1}", sol, sd, Pk, grid, price))
            solves.append(dict(repetition=k + 1, precompute_shared=t_sd, lifecycle=t_solve,
                               serialization=time.perf_counter() - t))
            write_json(args.out / "progress.json", dict(status=STATUS, completed_repetitions=k + 1,
                                                        planned=args.repeat, solves=solves))
        timing["lifecycle_solves"] = solves
        counts["lifecycle_solves"] = args.repeat
        receipt["stages"] = stages
        P_used = Pk
    else:
        if args.credit != "reference":
            raise SystemExit("Corrected-credit GE is blocked by the known entrant infeasibility (author decision)")
        if P.markov_equilibrium_method != "direct_brent" or int(P.I) != 1:
            raise RuntimeError("Reference GE contract is the one-market direct Brent solver")
        if args.initial_price_factor is None:
            raise SystemExit("--initial-price-factor is required (lead-selected start)")
        price = args.initial_price_factor * loaded.reference_price
        start_price = price
        from .verification.callcount import BudgetExceeded, CallCounter
        write_json(args.out / "parameters_initial.json", inputs.serialized(vars(P)))
        counter = CallCounter(solver, max_lifecycle=args.max_lifecycle, deadline_seconds=args.solve_deadline) \
            if args.count_calls else None
        t = time.perf_counter()
        try:
            if counter is not None:
                with counter:
                    sol, P_used, prices, fiscal = paygo.solve_balanced_initial_equilibrium(
                        model=solver, parameters=copy.deepcopy(P), b_grid=grid, initial_prices=price,
                        payroll_tax=float(P.tau_pay), marginal_tolerance=1e-9, fiscal_tolerance=1e-6)
            else:
                sol, P_used, prices, fiscal = paygo.solve_balanced_initial_equilibrium(
                    model=solver, parameters=copy.deepcopy(P), b_grid=grid, initial_prices=price,
                    payroll_tax=float(P.tau_pay), marginal_tolerance=1e-9, fiscal_tolerance=1e-6)
        except BudgetExceeded as exc:
            write_json(args.out / "budget_failure.json", dict(
                status="failed_budget_not_a_certificate", reason=str(exc), elapsed_solve_stage=time.perf_counter() - t,
                profiled=counter.summary() if counter else None))
            raise SystemExit(4)
        timing["equilibrium_solve_stage"] = time.perf_counter() - t
        timing["equilibrium_solve_stage_profiled"] = counter is not None
        write_json(args.out / "parameters_final.json", inputs.serialized(vars(P_used)))
        counts.update(ge_counters(sol.timings))
        if counter is not None:
            counts["profiled"] = counter.summary()
            counts["final_full_solution"] = counts["profiled"]["final_solution"]
        price = np.asarray(prices, dtype=float).reshape(-1)   # accepted GE price for all artifacts
        t = time.perf_counter()
        sd = solver.precompute_shared(P_used, grid)          # lab shared object for the oracle
        receipt["stages"] = [save_stage(args.out / "ge", sol, sd, P_used, grid, price)]
        timing["serialization_ge_stage"] = time.perf_counter() - t
        receipt.update(initial_price=start_price.tolist(), price=price.tolist(), fiscal=fiscal,
                       engine_timings=sol.timings, strict_converged=bool(sol.timings.get("strict_converged")))
        if not receipt["strict_converged"]:
            raise SystemExit("GE strict gate failed")
    receipt["adult_entry_relative_gap"] = float(getattr(sol, "adult_entry_stationary_relative_gap", float("nan")))
    if not args.skip_plots:
        t = time.perf_counter()
        target = args.out / "package_diagnostics"
        diagnostics.write_diagnostics(sol, P_used, target)
        produced = sorted(p.name for p in target.glob("*.png"))
        receipt["package_diagnostics"] = dict(
            produced=produced, missing_standard=sorted(set(manifest["standard_diagnostic_names"]) - set(produced)),
            note="package write_diagnostics only; the 17 standard dated plots are produced and hash-checked by acceptance_oracle")
        timing["reporting"] = time.perf_counter() - t
    timing["main_in_process"] = time.perf_counter() - t_main
    receipt.update(timing=timing, counts=counts)
    write_json(args.out / "receipt.json", receipt)
    print(json.dumps(jsonable(dict(status=STATUS, timing=timing, counts=counts)), indent=1))


if __name__ == "__main__":
    main()
