"""Bounded saved-population replay diagnostic; no model solve."""
from __future__ import annotations
import argparse, gzip, hashlib, importlib, json, pickle, sys, time
from pathlib import Path
import numpy as np

POLICY_SOURCE = "run_e5f_native_rental_access_diagnostic"

def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

def load(path: Path):
    with gzip.open(path, "rb") as f:
        return pickle.load(f)

def install(source_model: Path, adapter_tools: Path, local_tools: Path) -> None:
    # Source and source tools must precede local tools before unpickling.
    source_tools = source_model / "tools"
    for p in (local_tools, adapter_tools, source_tools, source_model):
        s = str(p.resolve())
        if s in sys.path: sys.path.remove(s)
        sys.path.insert(0, s)

def array_metrics(a, b):
    a, b = np.asarray(a), np.asarray(b)
    d = a - b
    return {"shape": list(a.shape), "l1": float(np.abs(d).sum()),
            "linf": float(np.max(np.abs(d))),
            "exact_equal": bool(np.array_equal(a, b)),
            "allclose_atol_1e-10_rtol_0": bool(np.allclose(a, b, atol=1e-10, rtol=0)),
            "nonzero_changed_count": int(np.count_nonzero(d)),
            "sum_a": float(a.sum()), "sum_b": float(b.sum()),
            "mass_delta": float(a.sum() - b.sum()),
            "relative_l1": float(np.abs(d).sum() / max(np.abs(b).sum(), 1e-300))}

def axis_sums(a):
    a = np.asarray(a)
    out = {}
    for axis in range(a.ndim):
        other = tuple(i for i in range(a.ndim) if i != axis)
        out[f"axis_{axis}"] = np.asarray(a.sum(axis=other)).tolist()
    return out

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--checkpoint", type=Path, required=True)
    p.add_argument("--source-root", type=Path, required=True)
    p.add_argument("--adapter-tools", type=Path, required=True)
    p.add_argument("--plan", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--expected-checkpoint-sha", required=True)
    p.add_argument("--source-sha", required=True)
    args = p.parse_args()
    started = time.time(); args.output.mkdir(parents=True, exist_ok=True)
    local_tools = Path(__file__).resolve().parent
    # Validate the staged adapter plan before loading the pickle.
    install(args.source_root, args.adapter_tools, local_tools)
    adapter = importlib.import_module("run_e5f_income_candidate_calibration")
    plan = adapter.read(args.plan)
    adapter.validate_plan(plan, require_source=True)
    if sha(args.checkpoint) != args.expected_checkpoint_sha:
        raise ValueError("checkpoint sha mismatch")
    src_manifest = Path(plan["source_manifest_path"])
    if not src_manifest.exists(): raise FileNotFoundError(src_manifest)
    if sha(src_manifest) != args.source_sha:
        raise ValueError("source manifest sha mismatch")
    x = load(args.checkpoint)
    stored = np.asarray(x["stationary_g_pre"])
    evaluation = x["evaluation"]
    stored_eval = np.asarray(evaluation.g_pre)
    rec = {"status": "started", "checkpoint_sha256": sha(args.checkpoint),
           "source_manifest_sha256": sha(src_manifest), "source_root": str(args.source_root),
           "plan": str(args.plan), "checkpoint": str(args.checkpoint),
           "stationary_g_pre_shape": list(stored.shape), "evaluation_g_pre_shape": list(stored_eval.shape),
           "stored_stationary_vs_evaluation_g_pre": array_metrics(stored, stored_eval),
           "stored_stationary_by_axis": axis_sums(stored), "evaluation_g_pre_by_axis": axis_sums(stored_eval),
           "stored_feasibility_projection_mass": float(getattr(evaluation, "feasibility_projection_mass", np.nan))}
    # Replay one period through the native continuation mapper, supplying saved policy.
    primitive = importlib.import_module("run_e5f_matched_pf_smoke")
    _, model = primitive.pf.transition.configure_sequential_model()
    primitive.pf.calendar.apply_fertility = primitive.pf.transition.apply_sequential_fertility
    primitive.pf.calendar.advance_calendar_distribution = primitive.pf.transition.advance_sequential_calendar_distribution
    P = x["parameters"]; grid = x["b_grid"]
    policy = evaluation.policy
    P._fert2_probs = np.asarray(policy.fert2_probs).copy()
    shared = model.precompute_shared(P, grid)
    price = np.asarray(policy.price)
    observed = primitive.pf.calendar.evaluate_period(price, stored, P, grid, shared,
        primitive.pf.calendar.SolveCounter(), supply_rule=x["supply_rule"], supplied_policy=policy)
    rec["observed_vs_stored_stationary_g_pre"] = array_metrics(observed.g_pre, stored)
    rec["observed_vs_stored_evaluation_g_pre"] = array_metrics(observed.g_pre, stored_eval)
    rec["observed_current_vs_saved"] = array_metrics(observed.g_current, evaluation.g_current)
    rec["observed_births_vs_saved"] = array_metrics(observed.births, evaluation.births)
    rec["observed_g_pre_by_axis"] = axis_sums(observed.g_pre)
    rec["observed_current_by_axis"] = axis_sums(observed.g_current)
    rec["observed_births_by_axis"] = axis_sums(observed.births)
    rec["observed_current_sum"] = float(np.asarray(observed.g_current).sum())
    rec["observed_births_sum"] = float(np.asarray(observed.births).sum())
    rec["saved_current_sum"] = float(np.asarray(evaluation.g_current).sum())
    rec["saved_births_sum"] = float(np.asarray(evaluation.births).sum())
    rec["elapsed_seconds"] = time.time() - started
    rec["status"] = "completed"
    (args.output / "report.json").write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n")

if __name__ == "__main__": main()
