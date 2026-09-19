"""Replay one saved full-grid stationary solution; no calibration or price search."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import platform
import signal
import sys
import time


class CompatibleUnpickler(pickle.Unpickler):
    """Read NumPy 2 / Python 3.13 namespace names on the migrated runtime."""
    def find_class(self, module, name):
        if module == "pathlib._local":
            module = "pathlib"
        if module.startswith("numpy._core"):
            module = module.replace("numpy._core", "numpy.core", 1)
        return super().find_class(module, name)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    for key in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[key] = "1"
    os.environ.update(NUMBA_DISABLE_JIT="0", MPLBACKEND="Agg",
                      NUMBA_CACHE_DIR=str(args.output.resolve() / "numba_cache"))
    # A single full-grid replay is expected to take seconds to minutes, not hours.
    def timeout(*_):
        raise TimeoutError("Single-solve benchmark exceeded its 15-minute cap")
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(900)
    root = args.source_root.resolve()
    manifest = json.loads(args.manifest.read_text())
    pins = {k: v for k, v in manifest["source_files"].items()
            if k.startswith("code/model/intergen_eqscale_seq_optimized/")}
    for relative, expected in pins.items():
        if sha(root / relative) != expected:
            raise RuntimeError("Frozen model source mismatch: " + relative)
    sys.path[:0] = [str(root / "code/model"), str(root / "code/model/tools")]
    import numpy as np
    import numba
    from intergen_eqscale_seq_optimized import solver, diagnostics
    from e5f_stationary_paygo import certify_initial_pension
    with gzip.open(args.checkpoint, "rb") as stream:
        packet = CompatibleUnpickler(stream).load()
    P, grid, reference = packet["parameters"], packet["b_grid"], packet["solution"]
    metadata = dict(machine=platform.machine(), platform=platform.platform(),
                    python=sys.version, numpy=np.__version__, numba=numba.__version__,
                    threads=1, J=int(P.J), Nb=len(grid),
                    checkpoint_sha256=sha(args.checkpoint), source_sha256=pins,
                    timing_scope="one full household Bellman and stationary distribution at saved equilibrium prices",
                    cache_state="fresh empty Numba disk cache; compilation included",
                    historical_torch_phase_seconds=reference.timings)
    (args.output / "inputs.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print("Starting one full-grid solve: 17 ages, 120 wealth nodes, one thread", flush=True)
    started = time.perf_counter()
    solution = solver.solve_markov_income_at_prices(reference.p_eq, P, grid, verbose=True)
    elapsed = time.perf_counter() - started
    comparisons = {}
    for name in ("V", "c_pol", "hR_pol", "bp_pol", "tenure_choice", "tenure_probs",
                 "loc_probs", "fert_probs", "fert_value", "g"):
        actual, expected = getattr(solution, name), getattr(reference, name)
        if actual is None or expected is None:
            comparisons[name] = {"equal": actual is expected}
        else:
            a, b = np.asarray(actual), np.asarray(expected)
            comparisons[name] = dict(max_abs_gap=float(np.max(np.abs(a-b))),
                                    equal=bool(np.allclose(a, b, rtol=1e-9, atol=1e-9)))
    demand = np.asarray(solution.housing_demand)
    supply = np.asarray(solution.housing_supply)
    residual = float(np.max(np.abs(demand-supply)/np.maximum(np.abs(supply), 1e-12)))
    solution.best_max_abs_rel_excess = residual
    fiscal = certify_initial_pension(solution.g, P, marginal_tolerance=1e-9, fiscal_tolerance=1e-6)
    passed = all(row["equal"] for row in comparisons.values()) and residual <= P.tol_eq
    result = dict(status="passed" if passed else "comparison_failed", wall_seconds=elapsed,
                  timings=solution.timings, comparisons=comparisons,
                  market_relative_residual=residual, fiscal=fiscal,
                  full_equilibrium_search=False, calibration_run=False)
    (args.output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in ("status", "wall_seconds", "timings", "market_relative_residual")}), flush=True)
    diagnostics.write_diagnostics(solution, P, args.output / "diagnostics")
    signal.alarm(0)
    if not passed:
        raise RuntimeError("Saved stationary replay differs; inspect comparison before comparing speed")


if __name__ == "__main__":
    main()
