"""Fixed-price, permanent-low-productivity comparison at the chain-13 point.

This experiment constructs a separate, consistent native input object. It does
not change the production loader, source snapshot, or stationary GE closure.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

for name in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
             "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
MODEL = ROOT / "code/model"
sys.path.insert(0, str(MODEL))
from production.inputs import load_inputs, validate_inputs
from production.parameter_files import load_parameter_file
from production.engine.e5f_stationary_paygo import (
    bind_initial_balanced_pension, stationary_age_income_mass,
)
from production.engine.e5f_social_security import fiscal_accounts
from production.engine.diagnostics import write_diagnostics
from production.engine.shared import income_at_state
from production.storage import StoredResult, save_case, load_case
from production import equilibrium

HERE = Path(__file__).resolve().parent
BEST = MODEL / "parameters/best_params.py"
PRICE = 0.7760569760205563


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def production_hashes() -> dict[str, str]:
    files = sorted((MODEL / "production").rglob("*.py"))
    files += [BEST, MODEL / "tools/model_playground.py", Path(__file__)]
    return {str(path.relative_to(ROOT)): sha(path) for path in files}


def prepare(phi: float):
    spec = load_parameter_file(BEST)
    if spec["price_guess"] != PRICE:
        raise RuntimeError("Chain-13 price changed")
    P, grid = load_inputs(spec["parameters"], spec["external_inputs"], spec["native_overrides"])
    if int(P.Nz) != 9 or not bool(P.native_fixed_reference_entry):
        raise RuntimeError("Expected nine distinct income nodes and fixed reference entry")
    old_z = P.z_grid.copy()
    old_w = P.z_weights.copy()
    old_cond = P.fixed_reference_entry_conditional.copy()
    low = int(np.argmin(old_z))
    if low != 0 or np.unique(old_z).size != 9:
        raise RuntimeError("Lowest income node or distinct-node contract changed")
    entry_marginal = np.sum(old_cond * old_w[None, :], axis=1)
    if not np.isclose(entry_marginal.sum(), 1, rtol=0, atol=2e-12):
        raise RuntimeError("Original entrant wealth marginal is invalid")
    P.z_weights = np.eye(9, dtype=float)[low]
    P.Pi_z = np.tile(P.z_weights, (9, 1))
    P.fixed_reference_entry_conditional[:, low] = entry_marginal
    P.phi = np.full_like(P.phi, phi)
    P, pension_receipt = bind_initial_balanced_pension(P, payroll_tax=P.tau_pay)
    validate_inputs(P, grid)
    counterfactual_marginal = np.sum(
        P.fixed_reference_entry_conditional * P.z_weights[None, :], axis=1)
    if not np.array_equal(counterfactual_marginal, entry_marginal):
        raise RuntimeError("Entrant wealth marginal changed")
    mass = stationary_age_income_mass(P)
    if np.any(mass[:, :, 1:] != 0):
        raise RuntimeError("Nonlow income states have anticipated mass")
    gross_annual = P.w_hat[0] * P.income_age_profile[:P.J_R] * P.z_grid[low]
    baseline_gross_annual = P.w_hat[0] * P.income_age_profile[:P.J_R]
    receipt = {
        "reference": "post-interest soft chain 13, no Estate A",
        "counterfactual": "everyone at original lowest productivity node throughout life",
        "price_held": PRICE, "phi": phi, "payroll_tax_held": float(P.tau_pay),
        "lowest_z": float(P.z_grid[low]), "lowest_index": low,
        "original_z_grid": old_z.tolist(), "original_z_weights": old_w.tolist(),
        "counterfactual_z_weights": P.z_weights.tolist(),
        "counterfactual_Pi_z": P.Pi_z.tolist(),
        "baseline_gross_annual_earnings_by_working_age": baseline_gross_annual.tolist(),
        "counterfactual_gross_annual_earnings_by_working_age": gross_annual.tolist(),
        "counterfactual_after_tax_period_income_by_age": [
            income_at_state(P, 0, j, float(P.z_grid[low])) for j in range(int(P.J))],
        "balanced_pension_period": float(P.pension),
        "predicted_fiscal_accounts": pension_receipt["predicted_accounts"],
        "baseline_entry_wealth_marginal": entry_marginal.tolist(),
        "counterfactual_entry_wealth_marginal": counterfactual_marginal.tolist(),
        "original_conditional_entry_sha256": hashlib.sha256(old_cond.tobytes()).hexdigest(),
        "best_params_sha256": sha(BEST),
        "source_sha256": production_hashes(),
    }
    return P, grid, spec, receipt


def run(phi: float, max_seconds: int):
    import numba
    numba.set_num_threads(1)
    P, grid, spec, receipt = prepare(phi)
    case = HERE / ("phi_08" if phi == 0.8 else "phi_10")
    if case.exists():
        raise FileExistsError(f"Refusing to overwrite {case}")
    case.mkdir(parents=True)
    (case / "input_contract.json").write_text(json.dumps(receipt, indent=2) + "\n")

    def expired(_signum, _frame):
        raise TimeoutError(f"Fixed-price solve exceeded {max_seconds} seconds")

    previous = signal.signal(signal.SIGALRM, expired)
    signal.alarm(max_seconds)
    start = time.monotonic()
    try:
        outcome = equilibrium.solve_at_price(P, grid, PRICE)
        sol, P, shared = outcome["solution"], outcome["P"], outcome["shared"]
    except Exception as exc:
        (case / "failure.json").write_text(json.dumps({"type": type(exc).__name__, "message": str(exc)}) + "\n")
        raise
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)
    elapsed = time.monotonic() - start
    fiscal = fiscal_accounts(sol.g, P)
    state_mass = np.asarray(sol.g).sum(axis=(0, 1, 2, 3, 5, 6))
    nonlow_mass = float(state_mass.sum() - state_mass[0])
    if (not np.isfinite(sol.g).all() or np.min(sol.g) < -1e-14
            or not np.isclose(sol.g.sum(), 1.0, rtol=0, atol=1e-8)):
        raise RuntimeError("Forward distribution is nonfinite, negative or unnormalized")
    for name in ("fert_probs", "tenure_probs", "loc_probs"):
        probabilities = getattr(sol, name, None)
        if probabilities is not None:
            probabilities = np.asarray(probabilities)
            if (not np.isfinite(probabilities).all() or np.min(probabilities) < -1e-12
                    or np.max(probabilities) > 1 + 1e-12):
                raise RuntimeError(f"Invalid {name}")
    if nonlow_mass > 1e-12 or abs(fiscal["scaled_pension_budget_residual"]) > 1e-6:
        raise RuntimeError(f"Income support or PAYGO failed: {nonlow_mass}, {fiscal}")
    result = StoredResult(sol, P, grid, PRICE, parameters=spec["parameters"],
                          label="fixed-price, all at lowest productivity")
    save_case(result, case, metadata={"experiment": "low_productivity_credit_v1",
                                       "input_contract_file": "input_contract.json",
                                       "fixed_price": True, "phi": phi})
    reopened, _ = load_case(case)
    if not np.array_equal(reopened.solution.g, sol.g):
        raise RuntimeError("Saved distribution roundtrip failed")
    arrays = {}
    for prefix, obj in (("solution", sol), ("shared", shared), ("inputs", P)):
        for name, value in vars(obj).items():
            if isinstance(value, np.ndarray) and not value.dtype.hasobject:
                arrays[f"{prefix}__{name}"] = value
    np.savez_compressed(case / "solution_arrays.npz", **arrays)
    write_diagnostics(sol, P, case / "standard_diagnostics")
    figures = sorted((case / "standard_diagnostics").glob("*.png"))
    if len(figures) != 17:
        raise RuntimeError(f"Expected 17 standard plots, found {len(figures)}")
    aggregates = result.aggregates()
    after_hashes = production_hashes()
    if after_hashes != receipt["source_sha256"]:
        raise RuntimeError("Production source or experiment driver changed during solve")
    summary = {"status": "complete_fixed_price_diagnostic", "phi": phi,
               "solve_seconds": elapsed, "price": PRICE, "lowest_z": receipt["lowest_z"],
               "income_state_mass": state_mass.tolist(), "nonlow_mass": nonlow_mass,
               "fiscal_accounts": fiscal, "aggregates": aggregates,
               "standard_diagnostic_pngs": [p.name for p in figures],
               "native_result_sha256": sha(case / "native_result.npz"),
               "solution_arrays_sha256": sha(case / "solution_arrays.npz"),
               "source_sha256_after": after_hashes}
    (case / "summary.json").write_text(json.dumps(summary, indent=2, default=lambda x: x.tolist() if isinstance(x, np.ndarray) else float(x)) + "\n")
    print(json.dumps({"phi": phi, "solve_seconds": elapsed, "pension": P.pension,
                      "nonlow_mass": nonlow_mass, "paygo_residual": fiscal["scaled_pension_budget_residual"],
                      "saved_case": str(case)}, indent=2))


def compare_phi_inputs():
    low, _, _, _ = prepare(0.8)
    high, _, _, _ = prepare(1.0)
    differences = []
    for name in vars(low):
        left, right = getattr(low, name), getattr(high, name)
        if isinstance(left, np.ndarray):
            same = np.array_equal(left, right, equal_nan=True)
        else:
            same = left == right
        if not same:
            differences.append(name)
    if differences != ["phi"]:
        raise RuntimeError(f"Phi arms differ in additional native fields: {differences}")


def supervise(phi: float, max_seconds: int):
    compare_phi_inputs()
    before = production_hashes()
    name = "phi_08" if phi == 0.8 else "phi_10"
    heartbeat = HERE / f"{name}_supervisor.json"
    log = HERE / f"{name}_run.log"
    if (HERE / name).exists():
        raise FileExistsError(f"Refusing existing {HERE / name}")
    with log.open("w") as stream:
        process = subprocess.Popen([sys.executable, str(Path(__file__)), "--phi", str(phi),
                                    "--max-seconds", str(max_seconds), "--worker"],
                                   stdout=stream, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        start = time.monotonic()
        status = "running"
        while process.poll() is None:
            elapsed = time.monotonic() - start
            rss = subprocess.run(["ps", "-o", "rss=", "-p", str(process.pid)],
                                 capture_output=True, text=True, check=False).stdout.strip()
            rss_kib = int(rss) if rss.isdigit() else None
            if elapsed > max_seconds or (rss_kib is not None and rss_kib > 8 * 1024 * 1024):
                os.killpg(process.pid, signal.SIGKILL)
                status = "time_limit" if elapsed > max_seconds else "memory_limit"
                process.wait()
                break
            heartbeat.write_text(json.dumps({"status": "running", "elapsed_seconds": elapsed,
                                             "rss_kib": rss_kib, "pid": process.pid}) + "\n")
            time.sleep(5)
        if status == "running":
            status = "complete" if process.returncode == 0 else "failed"
    after = production_hashes()
    heartbeat.write_text(json.dumps({"status": status,
                                     "elapsed_seconds": time.monotonic() - start,
                                     "returncode": process.returncode,
                                     "source_unchanged": before == after,
                                     "source_sha256_before": before,
                                     "source_sha256_after": after}, indent=2) + "\n")
    if before != after:
        raise RuntimeError("Source changed during supervised run")
    if status != "complete":
        raise RuntimeError(f"{name} stopped: {status}; see {log}")
    print(log.read_text())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--phi", type=float, choices=(0.8, 1.0), required=True)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--max-seconds", type=int, default=600)
    args = parser.parse_args()
    if args.max_seconds <= 0 or args.max_seconds > 600:
        raise ValueError("Per-arm time budget must be at most 600 seconds")
    if args.preflight:
        compare_phi_inputs()
        _, _, _, receipt = prepare(args.phi)
        print(json.dumps({k: receipt[k] for k in (
            "phi", "lowest_z", "balanced_pension_period", "price_held",
            "payroll_tax_held", "predicted_fiscal_accounts")}, indent=2))
    elif args.worker:
        run(args.phi, args.max_seconds)
    else:
        supervise(args.phi, args.max_seconds)
