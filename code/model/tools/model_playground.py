"""Interactive access to the saved soft reference and its pinned native solver.

Importing this module loads no model state and performs no lifecycle solve.
"""
from __future__ import annotations

import hashlib
import importlib
import json
import os
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

for _thread_var in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_thread_var] = "1"

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
EXPLORER_CASES = ROOT / "output/model/fixed_reference_economics_20260928/soft_timing_review_v1/explorer_cases.json"
SELECTION = ROOT / "output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_selected.json"
REFERENCE_CLOSURE = ROOT / "output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_postcheck/selected_postcheck/phase_b_ge/selected_repeat_final/closure.json"
V2_DIR = ROOT / "output/model/fixed_reference_economics_20260928/normalized_calibration_v2"
ORIGINAL_ENGINE = ROOT / "output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed/source"


def load_saved_solution(case: str = "soft") -> SimpleNamespace:
    """Load a saved explorer case after checking its recorded SHA-256; no solve."""
    config = json.loads(EXPLORER_CASES.read_text())
    spec = next((item for item in config["cases"] if item["id"] == case), None)
    if spec is None:
        raise KeyError(f"Unknown saved case {case!r}; available: {[x['id'] for x in config['cases']]}")
    path = Path(spec["arrays"])
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != spec["sha256"]:
        raise RuntimeError(f"Saved array fingerprint differs for {case}: {path}")
    with np.load(path, allow_pickle=False) as data:
        values = {key: data[key].copy() for key in data.files}
    values.update(price=float(spec["price"]), timing=spec["timing"], case_id=case)
    return SimpleNamespace(**values)


def _install_read_only_overlay() -> None:
    """Install the existing hash-checked, read-only local source mapping once."""
    marker = "_fertility_model_playground_overlay"
    overlay = ROOT / "output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/local_runtime/bootstrap.py"
    installed = getattr(sys, marker, None)
    if installed:
        if installed != str(overlay):
            raise RuntimeError("A different frozen-source overlay is already installed")
        return
    source = overlay.read_text()
    prefix = source.split("if '--preflight-context' in sys.argv:", 1)[0]
    if "DIGESTS=" not in prefix or "Frozen overlay write forbidden" not in prefix:
        raise RuntimeError("The authenticated read-only overlay has an unexpected shape")
    scope = {"__file__": str(overlay), "__name__": "authenticated_local_overlay"}
    exec(compile(prefix, str(overlay), "exec"), scope)
    setattr(sys, marker, str(overlay))


def load_reference_model() -> tuple[SimpleNamespace, np.ndarray, SimpleNamespace]:
    """Build the authenticated soft reference inputs and return ``P, b_grid, solver``.

    This only constructs inputs and shared model objects. It does not solve.
    The returned price is in ``P.reference_price``; the closure's derived H0 is
    retained in ``P.H0``.
    """
    _install_read_only_overlay()
    sys.path.insert(0, str(V2_DIR))
    import run_psi as v2
    import numba
    numba.set_num_threads(1)

    v2.native.verify_sources()
    selected = json.loads(SELECTION.read_text())
    source = ROOT / selected["source"]
    if hashlib.sha256(source.read_bytes()).hexdigest() != selected["source_sha256"]:
        raise RuntimeError("Selected soft parameter source hash differs")
    source_record = json.loads(source.read_text())
    if source_record[selected["source_key"]]["best"] != selected["selected"]:
        raise RuntimeError("Selected soft parameter record differs from its source")
    point = selected["selected"]["parameters"]
    if selected["selected"]["weight_fingerprint"] != v2.weight_fingerprint({}):
        raise RuntimeError("Soft reference weight fingerprint differs")
    if v2.native.target_identity(selected["selected"]["target_fit"]) != v2.CONFIG["base_target_contract"]:
        raise RuntimeError("Soft reference target contract differs")

    import copy
    original_lane = copy.deepcopy(v2.inputs.LANES["floor_s0"])
    _, bounds, _ = v2.inputs.seed_and_bounds("floor_s0")
    bounds = {name: tuple(value) for name, value in bounds.items()}
    if bounds["h_P"] != (0.1, 2.3):
        raise RuntimeError("Expected base h_P bound changed")
    bounds.update(h_P=(0.1, 2.6), psi_child=tuple(v2.CONFIG["psi_bounds"]))
    try:
        v2.inputs.LANES["floor_s0"].update(seed=dict(point), bounds=bounds, free_coordinates=list(point))
        P, b_grid = v2.inputs.proposal("floor_s0")
        P, _ = v2.inputs.entry(P, b_grid, "nonnegative_mean")
        with tempfile.TemporaryDirectory(prefix="model-playground-check-") as temp:
            P = v2.native.utility_checks(P, b_grid, "floor_s0", Path(temp))
    finally:
        v2.inputs.LANES["floor_s0"] = original_lane

    # Match the fixed-price experiment: corrected zero-credit contract.
    sys.path.insert(0, str(ORIGINAL_ENGINE))
    from small_credit_lab import credit
    from small_credit_lab.engine import solver
    expected_package = (ORIGINAL_ENGINE / "small_credit_lab").resolve()
    if Path(credit.__file__).resolve() != expected_package / "credit.py":
        raise RuntimeError(f"Unexpected credit module origin: {credit.__file__}")
    if Path(solver.__file__).resolve() != expected_package / "engine/solver.py":
        raise RuntimeError(f"Unexpected solver module origin: {solver.__file__}")
    credit.bind_engine_credit(P, "corrected", 0.0)
    closure = json.loads(REFERENCE_CLOSURE.read_text())
    if closure["normalized_population"] != 1.0 or P.N_target != 1.0:
        raise RuntimeError("Soft reference population normalization changed")
    P.H0 = np.array([closure["H0_derived"]])
    P.reference_price = float(closure["price"])
    P.native_inherited_distribution_evidence_dir = tempfile.mkdtemp(prefix="model-playground-inherited-state-")
    return P, b_grid, solver


def solve_at_price(P: SimpleNamespace, b_grid: np.ndarray, solver: SimpleNamespace, price: float) -> SimpleNamespace:
    """Run one explicit one-core stationary solve at ``price``; no market root."""
    if not np.isfinite(price) or price <= 0:
        raise ValueError("price must be finite and positive")
    shared = solver.precompute_shared(P, b_grid)
    return solver.solve_markov_income_at_prices(
        np.array([float(price)]), P, b_grid, SD=shared, fast_stats=False
    )


def main() -> None:
    global sol
    if len(sys.argv) > 1 and sys.argv[1] in {"-h", "--help"}:
        print(__doc__)
        print("Interactive: code/model/.venv/bin/python -i code/model/tools/model_playground.py")
        return
    sol = load_saved_solution()
    print(f"Saved case: {sol.case_id}; price={sol.price:.15g}; timing={sol.timing}")
    print(f"V shape={sol.V.shape}; wealth grid nodes={len(sol.b_grid)}; arrays loaded with no solve")
    print("Open Python: code/model/.venv/bin/python -i code/model/tools/model_playground.py")


if __name__ == "__main__":
    main()
