"""Native full-GE benchmark: ORIGINAL or LAB engine, identical call and gates.

VERIFICATION ONLY. Status label: native_GE_benchmark_only. This measures the
authenticated numerical engine through its own native GE gates (strict market
convergence, returned PAYGO fiscal certificate). It does NOT run the frozen
calendar/observer stack, so it certifies no 113-path packet, 14/31 tables or
17 standard plots; that certificate is separate (Torch, intact tree).

    PY -m refactor_lab.verification.native_ge_benchmark run --engine original|lab --reference-root ROOT \
        --bundle B --bundle-sha S --initial-price-factor 1.05 --credit reference --out OUT
    PY -m refactor_lab.verification.native_ge_benchmark compare-params A/parameters_final.json B/parameters_final.json

--engine original imports the genuine `intergen_eqscale_seq_optimized.solver`
and `code/model/tools/e5f_stationary_paygo.py` from ROOT after authenticating
the 13 pinned sources (engine_receipt_pass2_reference.json). Every model/tool
module imported from ROOT is then audited (after import and after the solve)
against those pins or the frozen source_manifest; any unpinned or changed
file stops the run. No frozen calibration runtime is imported.
--engine lab imports `refactor_lab.engine` from PYTHONPATH (the stage chosen by
the caller) and records every imported lab file hash; importing any original
package module is a failure.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import json
import math
import os
import pickle
import sys
import time
from pathlib import Path

import numpy as np

T_PROCESS = time.perf_counter()
STATUS = "native_GE_benchmark_only"
RUN_FIELDS = {"native_inherited_distribution_evidence_dir", "eq_iter"}
PINS = Path(__file__).resolve().parent / "engine_receipt_pass2_reference.json"
MANIFEST_REL = "output/model/fertility_identification_20260928/fixed_reference_manifest.json"
PKG = "intergen_eqscale_seq_optimized"


def sha(path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as s:
        for b in iter(lambda: s.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def jsonable(v):
    if isinstance(v, dict):
        return {str(k): jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [jsonable(x) for x in v]
    if isinstance(v, np.ndarray):
        return v.tolist() if v.size <= 64 else dict(shape=list(v.shape), dtype=str(v.dtype))
    if isinstance(v, np.generic):
        v = v.item()
    if isinstance(v, float) and not math.isfinite(v):
        return repr(v)
    return v if v is None or isinstance(v, (str, bool, int, float)) else repr(v)


def write(path: Path, value) -> None:
    path.write_text(json.dumps(jsonable(value), indent=1, sort_keys=True) + "\n")


def original_pins(root: Path) -> dict:
    """13 original sources: path -> pinned sha (verified BEFORE import)."""
    receipt = json.loads(PINS.read_text())
    pins = {str((root / m["source"]).resolve()): m["source_sha256"] for m in receipt["modules"].values()
            if not Path(m["source"]).is_absolute()}
    bad = [p for p, h in pins.items() if not Path(p).exists() or sha(p) != h]
    if len(pins) != 13 or bad:
        raise SystemExit(f"original source authentication failed ({len(pins)} pins): {bad}")
    return pins


def audit_imports(root: Path, pins: dict, frozen: dict) -> dict:
    """Every module file under ROOT/code/model must be pinned or match the frozen source_manifest."""
    model_dir = (root / "code/model").resolve()
    seen, problems = {}, []
    for name, mod in sorted(sys.modules.items()):
        f = getattr(mod, "__file__", None)
        if not f:
            continue
        p = Path(f).resolve()
        original_name = (name == "intergen_eqscale_seq_optimized"
                         or name.startswith("intergen_eqscale_seq_optimized.")
                         or name.startswith("e5f_"))
        if original_name and not p.is_relative_to(model_dir):
            raise SystemExit(f"original model/helper module resolved outside ROOT/code/model: {name} -> {p}")
        if not p.is_relative_to(model_dir) or "refactor_lab" in p.parts:
            continue
        rel, digest = str(p.relative_to(root.resolve())), sha(p)
        if str(p) in pins:
            ok, basis = digest == pins[str(p)], "pass2_pin"
        elif rel in frozen:
            ok, basis = digest == frozen[rel], "frozen_source_manifest"
        else:
            ok, basis = False, "UNPINNED"
        seen[rel] = dict(module=name, sha256=digest, basis=basis, ok=ok)
        if not ok:
            problems.append(rel)
    if problems:
        raise SystemExit("imported original model/tool source not authenticated: " + ", ".join(problems))
    return seen


def lab_imports() -> dict:
    out = {}
    for name, mod in sorted(sys.modules.items()):
        f = getattr(mod, "__file__", None)
        if name.startswith(PKG):
            raise SystemExit("lab engine run imported original package module: " + name)
        if f and name.startswith("refactor_lab"):
            out[name] = dict(file=str(Path(f).resolve()), sha256=sha(f))
    return out


def public(P) -> dict:
    from refactor_lab.inputs import serialized
    return serialized({k: v for k, v in vars(P).items() if not k.startswith("_")})


def run(a) -> int:
    if a.credit != "reference" or a.initial_price_factor != 1.05:
        raise SystemExit("this benchmark supports --credit reference and --initial-price-factor 1.05 only")
    root = a.reference_root.resolve()
    cache = os.environ.get("NUMBA_CACHE_DIR")
    if not cache or (Path(cache).exists() and any(Path(cache).iterdir())):
        raise SystemExit("a fresh, empty per-role NUMBA_CACHE_DIR is required (cold benchmark)")
    if any(os.environ.get(v) != "1" for v in ("OMP_NUM_THREADS", "NUMBA_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                               "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")):
        raise SystemExit("one thread required: set OMP/NUMBA/OPENBLAS/MKL/VECLIB thread variables to 1")
    a.out.mkdir(parents=True, exist_ok=False)
    timing = {}
    t = time.perf_counter()
    from refactor_lab import inputs, credit
    loaded = inputs.load_inputs(a.bundle, root, a.bundle_sha)
    P = copy.deepcopy(loaded.parameters)
    for flag in ("native_due_stayer_credit", "native_exact_inherited_distribution", "native_fixed_reference_entry",
                 "native_explicit_transaction_grid", "native_purchase_income", "native_exact_allocation_output"):
        if getattr(P, flag) is not True:
            raise SystemExit("reference contract flag absent: " + flag)
    credit.bind_engine_credit(P, "reference")
    P.native_inherited_distribution_evidence_dir = str(a.out / "inherited_state_failures")
    grid = loaded.b_grid.copy()
    start = a.initial_price_factor * loaded.reference_price
    initial_public = public(P)
    write(a.out / "parameters_initial.json", initial_public)
    timing["inputs"] = time.perf_counter() - t
    identity = dict(engine=a.engine, bundle=loaded.identity)
    t = time.perf_counter()
    if a.engine == "original":
        pins = original_pins(root)                                  # before any original import
        manifest = json.loads((root / MANIFEST_REL).read_text())
        sm_path = Path(manifest["source_manifest"]["path"])
        if sha(sm_path) != manifest["source_manifest"]["sha256"]:
            raise SystemExit("frozen source_manifest hash differs")
        frozen = json.loads(sm_path.read_text())["files"]
        sys.path[:0] = [str(root / "code/model"), str(root / "code/model/tools")]
        from intergen_eqscale_seq_optimized import solver
        from intergen_eqscale_seq_optimized import utils as solver_utils
        import e5f_stationary_paygo as paygo
        identity["imports_after_import"] = audit_imports(root, pins, frozen)
    else:
        from refactor_lab.engine import solver, e5f_stationary_paygo as paygo
        from refactor_lab.engine import utils as solver_utils
        identity["imports_after_import"] = lab_imports()
    import numba
    timing["import"] = time.perf_counter() - t
    np.testing.assert_array_equal(solver_utils.make_grid(P), grid)
    from refactor_lab.verification.callcount import BudgetExceeded, CallCounter
    counter = CallCounter(solver, max_lifecycle=18, deadline_seconds=900)
    t = time.perf_counter()
    try:
        with counter:
            sol, P_final, prices, fiscal = paygo.solve_balanced_initial_equilibrium(
                model=solver, parameters=copy.deepcopy(P), b_grid=grid, initial_prices=start,
                payroll_tax=float(P.tau_pay), marginal_tolerance=1e-9, fiscal_tolerance=1e-6)
    except BudgetExceeded as exc:
        write(a.out / "budget_failure.json", dict(status="failed_budget_not_a_benchmark", reason=str(exc),
                                                  profiled=counter.summary()))
        return 4
    timing["solve_stage_profiled"] = time.perf_counter() - t
    identity["imports_after_solve"] = (audit_imports(root, pins, frozen) if a.engine == "original" else lab_imports())
    # Native gates (must hold; the full observer certificate is separate).
    tim = sol.timings
    gates = dict(strict_converged=bool(tim.get("strict_converged")), final_eq_error=tim.get("final_eq_error"),
                 tol_eq=float(P.tol_eq), marginal_gate=fiscal.get("marginal_gate"), fiscal_gate=fiscal.get("fiscal_gate"),
                 normalized_age_income_max_gap=fiscal.get("normalized_age_income_max_gap"))
    # Parameter identity: only run/output fields and stationary-PAYGO-derived values may change.
    derived, _ = paygo.bind_initial_balanced_pension(copy.deepcopy(P), payroll_tax=float(P.tau_pay))
    derived_public, final_public = public(derived), public(P_final)
    diffs, unexplained = {}, []
    for k in sorted(set(initial_public) | set(final_public)):
        if initial_public.get(k) == final_public.get(k):
            continue
        if k in RUN_FIELDS:
            diffs[k] = "run_or_output_field"
        elif final_public.get(k) == derived_public.get(k):
            diffs[k] = "derived_stationary_paygo"
        else:
            unexplained.append(k)
    write(a.out / "parameters_final.json", final_public)
    # Renewal: diagnostic only, versus the frozen reference's own recorded gate.
    ref = json.loads((root / MANIFEST_REL).read_text())["inherited_gates"]["adult_entry_gate"]
    ref_gap = abs(ref["entry_residual"]) / max(ref["entry_E"], ref["potential_B"], 1e-12)
    gap = float(sol.adult_entry_stationary_relative_gap)
    renewal = dict(relative_gap=gap, reference_relative_gap_from_adult_entry_gate=ref_gap, diagnostic_threshold=1e-6,
                   classification="consistent" if abs(gap - ref_gap) <= 1e-6 else "deviates",
                   note="lead-selected diagnostic; no psi normalization or closure change")
    t = time.perf_counter()
    price = np.asarray(prices, dtype=float).reshape(-1)
    sd = solver.precompute_shared(P_final, grid)
    arrays = {k: v for k, v in vars(sol).items() if isinstance(v, np.ndarray) and v.dtype != object}
    arrays.update({"shared." + k: v for k, v in vars(sd).items() if isinstance(v, np.ndarray) and v.dtype != object})
    arrays.update({"benchmark.b_grid": grid, "benchmark.price": price, "benchmark.start_price": np.asarray(start)})
    nonfinite = {k: int(np.count_nonzero(~np.isfinite(v))) for k, v in arrays.items()
                 if v.dtype.kind in "fc" and not np.isfinite(v).all()}
    np.savez(a.out / "solution_arrays.npz", **arrays)
    with gzip.open(a.out / "stage.pkl.gz", "wb", compresslevel=1) as stream:
        pickle.dump(dict(solution=sol, shared=sd, parameters=P_final, b_grid=grid, price=price), stream, protocol=5)
    timing["serialization"] = time.perf_counter() - t
    timing["process_total_in_python"] = time.perf_counter() - T_PROCESS
    cache = os.environ.get("NUMBA_CACHE_DIR")
    receipt = dict(status=STATUS, engine=a.engine, start_price=start, price=price, gates=gates, fiscal=fiscal,
                   renewal=renewal, parameter_differences=diffs, unexplained_parameter_differences=unexplained,
                   nonfinite_arrays=nonfinite, arrays=len(arrays), profiled=counter.summary(),
                   engine_timings=tim, timing=timing, identity=identity,
                   runtime=dict(__import__("refactor_lab.verification.baseline_identity", fromlist=["x"]).runtime_identity(),
                                numba_threads=numba.get_num_threads(), numba_cache_dir=cache,
                                pythonpath=os.environ.get("PYTHONPATH"), executable=os.path.realpath(sys.executable)),
                   scope="native GE gates only; frozen calendar/observer certificate (113 paths, 14/31, 17 plots) is separate")
    write(a.out / "receipt.json", receipt)
    failed = [k for k, v in (("strict_market", gates["strict_converged"]), ("marginal_gate", gates["marginal_gate"] is True),
                              ("fiscal_gate", gates["fiscal_gate"] is True), ("parameter_identity", not unexplained),
                              ("price_finite", np.isfinite(price).all()), ("distribution_finite", "g" not in nonfinite),
                              ("solution_arrays_finite", not nonfinite))
              if not v]
    print(json.dumps(dict(status=STATUS if not failed else "failed_native_gates", failed=failed,
                          price=price.tolist(), solve_stage_profiled=timing["solve_stage_profiled"],
                          lifecycle_evaluations=counter.summary()["lifecycle_evaluations"],
                          final_solution=counter.summary()["final_solution"]), indent=1))
    return 1 if failed else 0


def compare_params(a) -> int:
    x, y = json.loads(a.a.read_text()), json.loads(a.b.read_text())
    diff = sorted(k for k in set(x) | set(y) if k not in RUN_FIELDS and x.get(k) != y.get(k))
    print(json.dumps(dict(identical_except_run_fields=not diff, differing=diff, excluded=sorted(RUN_FIELDS)), indent=1))
    return 1 if diff else 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--engine", choices=("original", "lab"), required=True)
    r.add_argument("--reference-root", type=Path, required=True)
    r.add_argument("--bundle", type=Path, required=True)
    r.add_argument("--bundle-sha", required=True)
    r.add_argument("--initial-price-factor", type=float, required=True)
    r.add_argument("--credit", choices=("reference", "corrected"), required=True)
    r.add_argument("--out", type=Path, required=True)
    c = sub.add_parser("compare-params")
    c.add_argument("a", type=Path)
    c.add_argument("b", type=Path)
    a = ap.parse_args()
    return run(a) if a.cmd == "run" else compare_params(a)


if __name__ == "__main__":
    sys.exit(main())
