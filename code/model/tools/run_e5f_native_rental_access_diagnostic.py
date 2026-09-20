"""Native 2x2 rental-access by mortgage-financing diagnostic.

This is a fixed-price, partial-equilibrium diagnostic.  It solves three new
arms (an exact baseline control, open rental access at ``phi=.8``, and open
rental access at ``phi=1``) and reuses the saved capped baseline and
mortgage-only arrays for the fourth cell.  It makes no GE, recalibration, or
mediation claim.
"""
from __future__ import annotations

import argparse, copy, gzip, hashlib, importlib, json, os, signal, subprocess, sys, time
from pathlib import Path
from typing import Any, Mapping

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "code/model/tools"
DEFAULT_CHECKPOINT = ROOT / "output/model/paper_baseline_sep14/replay_20260917/native_output/raw/repetition_02/initial_state.pkl.gz"
DEFAULT_REPLAY = ROOT / "output/model/paper_baseline_sep14/replay_20260917"
DEFAULT_SOURCE_ROOT = ROOT / "tmp/paper_baseline_sep14/code/model"
DEFAULT_OUTPUT = ROOT / "output/model/native_financing_diagnostic_20260919/rental_access"
DEFAULT_PRIOR_INPUT = ROOT / "output/model/native_financing_diagnostic_20260919/run"
CHECKPOINT_SHA256 = "3322a61994fb3654d67f4b1d6cf2d0f7cacbb3668d06a417e192ee363c174993"
SOLVER_SHA256 = "2992412586b81cef3a3e58d92191bb51f54d3f9cc600d7675bbadaed7d1682da"
POLICY_NAMES = ("V", "c_pol", "hR_pol", "bp_pol", "tenure_choice", "tenure_probs", "loc_probs", "fert_probs", "fert_value", "fert2_probs", "price")
ARMS = {
    "capped_phi08": ("capped", 0.8),
    "capped_phi1": ("capped", 1.0),
    "open_phi08": ("open", 0.8),
    "open_phi1": ("open", 1.0),
}
NEW_CASES = {"baseline_control": ("capped", 0.8), "rental_access": ("open", 0.8), "rental_access_mortgage": ("open", 1.0)}
CASE_ORDER = tuple(NEW_CASES)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def write(path: Path, obj: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str) + "\n")
    tmp.replace(path)


def packet(path: Path) -> Any:
    with gzip.open(path, "rb") as f: return __import__("pickle").load(f)


def install_paths(source: Path) -> None:
    for p in (source, source / "tools", BASE):
        if str(p.resolve()) not in sys.path: sys.path.insert(0, str(p.resolve()))


def validate_contract(checkpoint: Path, replay: Path, source: Path) -> dict[str, str]:
    selected = replay / "native_output/selected_checkpoint.json"
    if not selected.exists(): raise ValueError("checkpoint contract hash mismatch: selected checkpoint record missing")
    actual = sha(checkpoint)
    if actual != CHECKPOINT_SHA256 or json.loads(selected.read_text()).get("checkpoint_sha256") != CHECKPOINT_SHA256:
        raise ValueError("checkpoint contract hash mismatch")
    solver = source / "intergen_eqscale_seq_optimized/solver.py"
    if sha(solver) != SOLVER_SHA256: raise ValueError("frozen solver hash mismatch")
    return {"checkpoint_sha256": actual, "solver_sha256": sha(solver)}


def changed(a: Any, b: Any) -> list[str]:
    def eq(x: Any, y: Any) -> bool:
        try: return np.array_equal(x, y) if isinstance(x, np.ndarray) or isinstance(y, np.ndarray) else bool(x == y)
        except Exception: return False
    return sorted(k for k in set(vars(a)) | set(vars(b)) if not eq(getattr(a, k, None), getattr(b, k, None)))


def policy_arrays(policy: Any) -> dict[str, np.ndarray]:
    out = {n: np.asarray(getattr(policy, n, None)) for n in POLICY_NAMES}
    missing = [n for n in POLICY_NAMES if getattr(policy, n, None) is None]
    if missing: raise ValueError(f"missing mandatory policy arrays: {missing}")
    return out


def compare_exact(a: Any, b: Any) -> None:
    for n, x in policy_arrays(a).items(): np.testing.assert_allclose(policy_arrays(b)[n], x, atol=1e-10, rtol=0, err_msg=n)


def refresh_finance(P: Any) -> None:
    importlib.import_module("intergen_eqscale_seq_optimized.parameters").build_debt_caps(P)
    if not np.isfinite(P.debt_taper_weights).all() or not np.isfinite(P.debt_caps).all(): raise ValueError("nonfinite rebuilt finance arrays")


def arm_parameters(base: Any, cap: str, phi: float) -> Any:
    P = copy.deepcopy(base)
    P.phi = np.full_like(np.asarray(base.phi, dtype=float), phi)
    if cap == "open": P.hR_max = float(np.max(np.asarray(base.H_own, dtype=float)))
    if cap == "open" and P.hR_max < float(base.hR_max): raise ValueError("expanded renter cap is below the native baseline cap")
    if phi != 0.8: refresh_finance(P)
    allowed = {"phi", "hR_max", "debt_taper_weights", "debt_caps"}
    altered = set(changed(base, P))
    if altered - allowed: raise ValueError(f"unaccounted parameter changes: {sorted(altered - allowed)}")
    return P


def gates(ev: Any) -> dict[str, float]:
    out = {"pre_mass": float(ev.g_pre.sum()), "post_fertility_mass": float(ev.g_post_fertility.sum()), "current_mass": float(ev.g_current.sum()), "birth_mass": float(np.asarray(ev.births).sum())}
    if not np.isfinite(list(out.values())).all(): raise ValueError(f"nonfinite mass: {out}")
    if any(np.min(a) < -1e-12 for a in (ev.g_pre, ev.g_post_fertility, ev.g_current)): raise ValueError("negative household mass")
    if max(abs(out["pre_mass"] - out["post_fertility_mass"]), abs(out["pre_mass"] - out["current_mass"])) > 1e-10: raise ValueError(f"period mass gate failed: {out}")
    for n, a in policy_arrays(ev.policy).items():
        if not np.isfinite(a).all(): raise ValueError(f"nonfinite policy: {n}")
    for n in ("tenure_probs", "loc_probs", "fert_probs", "fert2_probs"):
        a = policy_arrays(ev.policy)[n]
        if np.min(a) < -1e-12 or np.max(a) > 1 + 1e-12: raise ValueError(f"probability gate failed: {n}")
    return out


def native_solve(x: Mapping[str, Any], P: Any) -> tuple[Any, Any, Any, Any]:
    model = importlib.import_module("intergen_eqscale_seq_optimized.solver")
    primitive = importlib.import_module("run_e5f_matched_pf_smoke")
    price, grid = np.asarray(x["evaluation"].policy.price), x["b_grid"]
    _, model = primitive.pf.transition.configure_sequential_model()
    primitive.pf.calendar.apply_fertility = primitive.pf.transition.apply_sequential_fertility
    primitive.pf.calendar.advance_calendar_distribution = primitive.pf.transition.advance_sequential_calendar_distribution
    shared = model.precompute_shared(P, grid)
    sol = model.solve_markov_income_at_prices(price, P, grid, verbose=False, fast_stats=False)
    # Required by the native continuation mapper; keep this mutation inside the arm copy.
    P._fert2_probs = np.asarray(sol.fert2_probs).copy()
    policy = primitive.pf.calendar.policy_from_solution(sol, price, P, grid, shared)
    ev = primitive.pf.calendar.evaluate_period(price, x["stationary_g_pre"], P, grid, shared, primitive.pf.calendar.SolveCounter(), supply_rule=x["supply_rule"], supplied_policy=policy)
    budget = primitive.dated_budget(ev, P, shared, grid, float(P.user_cost_rate * price[0]))
    return ev, budget, shared, model


def standard_graphs(x: Mapping[str, Any], P: Any, ev: Any, shared: Any, model: Any, out: Path) -> dict[str, Any]:
    try:
        audit = importlib.import_module("run_e5f_independent_numerical_audit")
        audit.standard_diagnostics({"parameters": P, "b_grid": x["b_grid"], "evaluation": ev, "shared": shared}, out, validate_production_young=False)
        graphs = sorted((out / "standard_diagnostics").glob("*.png"))
        if len(graphs) != 17: raise RuntimeError(f"expected 17 standard PNGs, found {len(graphs)}")
        return {"status": "completed", "count": len(graphs), "paths": [str(p) for p in graphs]}
    except Exception as exc:
        return {"status": "unavailable", "reason": f"{type(exc).__name__}: {exc}", "route": "reconstruct the audit packet inputs and rerun standard_diagnostics"}


def build_comparison(args: argparse.Namespace) -> dict[str, Any]:
    """Merge saved capped cells and new open cells using the validated report formulas."""
    install_paths(args.source_root)
    report = importlib.import_module("build_e5f_native_financing_report")
    meta, _ = report._load_metadata(args.checkpoint, args.source_root)
    sources = {
        "capped_phi08": args.prior_input / "cases/baseline_01",
        "capped_phi1": args.prior_input / "cases/mortgage_only_03",
        "open_phi08": args.output / "cases/rental_access",
        "open_phi1": args.output / "cases/rental_access_mortgage",
    }
    rows = []
    missing = {}
    for arm, path in sources.items():
        try:
            row = report._measure(arm, report._load_case(path), meta)
            row.update({"arm": arm, "access": "open" if arm.startswith("open") else "capped", "phi": 1.0 if arm.endswith("phi1") else 0.8})
            rows.append(row)
        except Exception as exc:
            missing[arm] = f"{type(exc).__name__}: {exc}"
    by = {r["arm"]: r for r in rows}
    if len(by) == 4:
        by["interaction_births"] = ((by["open_phi1"]["births_per_initial_household"] - by["open_phi08"]["births_per_initial_household"]) - (by["capped_phi1"]["births_per_initial_household"] - by["capped_phi08"]["births_per_initial_household"]))
    args.output.mkdir(parents=True, exist_ok=True)
    fields = sorted({k for r in rows for k in r} | {"arm", "status", "reason"})
    with (args.output / "comparisons.csv").open("w", newline="") as f:
        import csv
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n"); w.writeheader()
        for r in rows: w.writerow({**r, "status": "available", "reason": ""})
        for arm, reason in missing.items(): w.writerow({"arm": arm, "status": "missing", "reason": reason})
    lines = ["# Native rental-access × mortgage-financing diagnostic", "", "Fixed-price partial-equilibrium comparison using the same saved beginning-of-period mass `g_pre` for all four cells.", "", "| Arm | Access | Phi | Births / initial household | First births / initial household | Mean physical rooms |", "|---|---|---:|---:|---:|---:|"]
    for r in rows: lines.append("| {arm} | {access} | {phi:.1f} | {births_per_initial_household:.8g} | {first_birth_flow_initial_n0:.8g} | {mean_physical_rooms:.8g} |".format(**r))
    for arm, reason in missing.items(): lines.append(f"| {arm} | missing | missing | missing | missing | missing |\nReason: {reason}")
    if "interaction_births" in by:
        lines += ["", f"Interaction in births per initial household: `{by['interaction_births']:.12g}` = (open phi=1 − open phi=.8) − (capped phi=1 − capped phi=.8).", "This is a descriptive interaction at frozen prices and population; it is not identified mediation, a causal decomposition, or a recalibration."]
    (args.output / "comparison_report.md").write_text("\n".join(lines) + "\n")
    write(args.output / "comparison_metadata.json", {"status": "complete" if not missing else "incomplete", "arms": list(sources), "missing": missing, "interaction": by.get("interaction_births"), "normalization": "validated native full flows, first-birth n=0 mass loss, and physical rooms", "finite_status_checks": "build_e5f_native_financing_report._measure completed for every available arm; it includes probability postchecks and finite numeric checks"})
    if missing: raise RuntimeError(f"comparison incomplete; see comparison_metadata.json: {missing}")
    return {"status": "complete", "interaction_births": by["interaction_births"], "rows": len(rows)}


def run_case(args: argparse.Namespace) -> dict[str, Any]:
    contract = validate_contract(args.checkpoint, args.replay, args.source_root); install_paths(args.source_root)
    x = packet(args.checkpoint); base = x["parameters"]; cap, phi = NEW_CASES[args.case]
    P = arm_parameters(base, cap, phi); parameter_changes = changed(base, P)
    if args.case == "baseline_control" and parameter_changes: raise ValueError(f"baseline changed fields: {parameter_changes}")
    start = time.monotonic(); ev, budget, shared, model = native_solve(x, P)
    if float(budget.get("budget_excess_mass", np.inf)) > 2e-10 or float(budget.get("maximum_occupied_excess", np.inf)) > 1e-9: raise ValueError(f"native budget gate failed: {budget}")
    if not np.array_equal(ev.g_pre, x["stationary_g_pre"]): raise ValueError("period mapper changed saved beginning-of-period mass")
    if args.case == "baseline_control":
        compare_exact(x["evaluation"].policy, ev.policy)
        for n in ("g_current", "births"): np.testing.assert_allclose(getattr(ev, n), getattr(x["evaluation"], n), atol=1e-10, rtol=0, err_msg=n)
    out = args.output / "cases" / args.case
    if out.exists(): raise FileExistsError(f"refusing to overwrite existing case: {out}")
    out.mkdir(parents=True)
    import run_e5f_independent_numerical_audit as audit
    audit_result = audit.policy_array_audit({"evaluation": ev, "parameters": P, "b_grid": x["b_grid"]}, out)
    if audit_result["occupied_negative_steps"]: raise ValueError("occupied value monotonicity gate failed")
    np.savez_compressed(out / "arrays.npz", g_pre=ev.g_pre, g_post_fertility=ev.g_post_fertility, g_current=ev.g_current, births=ev.births, **policy_arrays(ev.policy))
    graphs = standard_graphs(x, P, ev, shared, model, out)
    result = {"status": "completed", "case": args.case, "cap": cap, "phi": phi, "hR_max": float(P.hR_max), "elapsed_seconds": time.monotonic() - start, "contract": contract, "parameter_changes": parameter_changes, "mass": gates(ev), "budget": budget, "relative_market_residual": float(ev.relative_market_residual), "policy_audit": audit_result, "standard_graphs": graphs, "interpretation": "fixed-price PE interaction diagnostic; not identified mediation or recalibration"}
    write(out / "receipt.json", result); return result


def inspect(args: argparse.Namespace) -> dict[str, Any]:
    contract = validate_contract(args.checkpoint, args.replay, args.source_root); install_paths(args.source_root); x = packet(args.checkpoint); base = x["parameters"]
    rows = []
    for case, (cap, phi) in NEW_CASES.items():
        P = arm_parameters(base, cap, phi); rows.append({"case": case, "cap": cap, "phi": phi, "hR_max": float(P.hR_max), "changed_fields": changed(base, P)})
    return {"status": "inspect_only", "contract": contract, "cases": rows, "exact_loop": list(CASE_ORDER), "scope": "fixed-price full-lifecycle PE; all contract and overrides checked before solve"}


def drive_case_order(case_runner: Any, progress_writer: Any, comparison_runner: Any) -> None:
    """Small orchestration seam used by the no-solve sequence smoke test."""
    for i, case in enumerate(CASE_ORDER, 1):
        case_runner(case)
        progress_writer(i, case)
    comparison_runner()


def sequence(args: argparse.Namespace) -> None:
    started = time.monotonic(); write(args.output / "inspect.json", inspect(args))
    for i, case in enumerate(CASE_ORDER, 1):
        left = args.total_budget_seconds - (time.monotonic() - started)
        if left <= 0: raise TimeoutError("total diagnostic budget exhausted")
        cmd = [sys.executable, str(Path(__file__).resolve()), "--mode", "case", "--case", case, "--checkpoint", str(args.checkpoint), "--replay", str(args.replay), "--source-root", str(args.source_root), "--output", str(args.output)]
        child = subprocess.Popen(cmd, start_new_session=True); limit = min(args.case_budget_seconds, left); case_start = time.monotonic()
        try:
            while child.poll() is None:
                if time.monotonic() - case_start > limit: raise subprocess.TimeoutExpired(cmd, limit)
                write(args.output / "heartbeat.json", {"case": case, "completed": i - 1, "elapsed_seconds": time.monotonic() - started})
                time.sleep(min(30., max(0.1, limit - (time.monotonic() - case_start))))
        except subprocess.TimeoutExpired:
            write(args.output / "timeout.json", {"status": "timeout", "case": case, "limit_seconds": limit}); os.killpg(child.pid, signal.SIGTERM)
            try: child.wait(timeout=15)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL); child.wait(timeout=15)
            raise
        if child.returncode: raise RuntimeError(f"{case} child exited {child.returncode}")
        write(args.output / "latest_completed.json", {"status": "progress", "completed": i, "of": 3, "last_case": case, "elapsed_seconds": time.monotonic() - started})
    build_comparison(args)


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(); p.add_argument("--mode", choices=("inspect-only", "case", "sequence"), default="inspect-only"); p.add_argument("--case", choices=tuple(NEW_CASES)); p.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT); p.add_argument("--replay", type=Path, default=DEFAULT_REPLAY); p.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT); p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT); p.add_argument("--prior-input", type=Path, default=DEFAULT_PRIOR_INPUT); p.add_argument("--case-budget-seconds", type=float, default=600); p.add_argument("--total-budget-seconds", type=float, default=1500); a = p.parse_args(argv)
    try:
        if a.mode == "inspect-only": write(a.output / "inspect.json", inspect(a))
        elif a.mode == "case": run_case(a)
        else: sequence(a)
    except Exception as exc:
        write(a.output / f"failure_{a.case or a.mode}.json", {"status": "failed", "type": type(exc).__name__, "error": str(exc)}); raise


if __name__ == "__main__": main()
