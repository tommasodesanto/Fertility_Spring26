"""Credit-policy retention replay: three household solves against the frozen refit contract.

Replays three arms of the completed ``finance_dose_refit_v2`` refit_new_income
factorial (phi=.8, cap=6, lambda in {0, 1, 5}) through the tested factorial
runner copied from that runtime, retains the policy and distribution arrays the
old run discarded, reproduces the old scalar receipts within 1e-10, and reports
whether the arm policies coincide.  Coincidence is a finding, not a gate.
Fixed-price partial equilibrium; no new target, parameter, entry, population,
GE or calibration change.  Runs only from a staged experiment root whose
launch manifest pins every input by SHA256.
"""
from __future__ import annotations
import argparse, hashlib, importlib, json, os, signal, sys, time
from pathlib import Path
from typing import Any
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
CASES = ((0.8, 0.0, 6.0, "case_0_baseline_lambda0"), (0.8, 1.0, 6.0, "case_1_credit_lambda1"), (0.8, 5.0, 6.0, "case_2_credit_lambda5"))
CORE = "intergen_eqscale_seq_optimized"
RUNTIME_HELPERS = ("run_e5f_financing_factorial", "run_e5f_native_rental_access_diagnostic", "run_e5f_native_financing_diagnostic",
                   "run_e5f_native_income_cohort_diagnostic", "build_e5f_native_financing_report")
EXTRA_ARRAYS = ("g_pre", "g_post_fertility", "g_current", "births")
METRIC_KEYS = ("birth_flow", "first_birth_flow", "mean_rooms", "ownership", "renter_mass", "owner_mass", "market_residual", "pre_mass")
COHORT_KEYS = ("cumulative_explicit_births_per_initial_household", "first_births_per_initial_household", "first_birth_mean_age", "initial_mass", "rows")
AUDIT_KEYS = ("maximum_occupied_value_drop", "occupied_negative_steps", "lower_node_mass_at_negative_steps", "share_pre_choice_mass_at_negative_steps")
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS", "NUMBA_DISABLE_JIT", "NUMBA_CACHE_DIR")
TOL = 1e-10


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def read_json(path: Path) -> Any: return json.loads(Path(path).read_text())


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str) + "\n"); tmp.replace(path)


def mod(name: str) -> Any: return importlib.import_module(name)


def under(p: Path, root: Path) -> bool:
    try: p.relative_to(root); return True
    except ValueError: return False


def manifest_checks(m: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Fail-closed SHA256 checks of every pinned input, in both the copied and the original locations."""
    exp, rows = Path(m["experiment_root"]), {}
    def check(key: str, path: Any, expected: str) -> None:
        p = Path(path); actual = sha(p) if p.is_file() else None
        rows[key] = {"path": str(p), "expected": expected, "actual": actual, "ok": actual == expected}
    check("checkpoint", m["checkpoint"]["path"], m["checkpoint"]["sha256"])
    check("summary_runtime", m["summary"]["remote_path"], m["summary"]["sha256"])
    check("plan_runtime", m["plan"]["remote_path"], m["plan"]["sha256"]); check("plan_copy", m["plan"]["copy"], m["plan"]["sha256"])
    for name, pin in m["core"].items():
        check(f"core_copy/{name}", Path(m["frozen_copy"]) / "code/model" / CORE / name, pin)
        check(f"core_frozen/{name}", Path(m["frozen_root"]) / "code/model" / CORE / name, pin)
    for name, pin in m["runtime_helpers"].items():
        check(f"helper_copy/{name}", exp / "code/model/tools" / name, pin)
        check(f"helper_runtime/{name}", Path(m["runtime_root"]) / "code/model/tools" / name, pin)
    check("wrapper", m["wrapper"]["path"], m["wrapper"]["sha256"])
    check("frozen_source_manifest", m["source_manifest"]["path"], m["source_manifest"]["sha256"])
    for rel, pin in read_json(m["source_manifest"]["path"]).items():
        check(f"frozen_python_copy/{rel}", Path(m["frozen_copy"]) / "code/model" / rel, pin)
        check(f"frozen_python_original/{rel}", Path(m["frozen_root"]) / "code/model" / rel, pin)
    if Path(__file__).resolve() != Path(m["wrapper"]["path"]).resolve():
        raise RuntimeError(f"running wrapper {Path(__file__).resolve()} is not the pinned wrapper {m['wrapper']['path']}")
    copied = sorted(p.name for p in (exp / "code/model/tools").glob("*.py"))
    expected = sorted(set(m["runtime_helpers"]) | {Path(m["wrapper"]["path"]).name})
    if copied != expected: raise RuntimeError(f"copied helper set {copied} differs from allowed {expected}")
    bad = [k for k, r in rows.items() if not r["ok"]]
    if bad: raise RuntimeError("fail-closed hash check failed: " + ", ".join(bad))
    return rows


def check_sys_path(exp: Path) -> None:
    for entry in sys.path:
        p = (Path(entry) if entry else Path.cwd()).resolve()
        if ("Fertility_Spring26" in str(p) or entry == "") and not under(p, exp):
            raise RuntimeError(f"foreign sys.path entry: {entry!r}")


def install(m: dict[str, Any]) -> tuple[Any, Any, Any]:
    """Import the copied runtime helpers and point every old install function at the copied frozen core."""
    exp = Path(m["experiment_root"]).resolve(); frozen = Path(m["frozen_copy"]).resolve() / "code/model"
    if ROOT.resolve() != exp: raise RuntimeError(f"wrapper ROOT {ROOT} is not the experiment root {exp}")
    if not (frozen / CORE / "solver.py").is_file(): raise RuntimeError(f"copied frozen core missing under {frozen}")
    check_sys_path(exp)
    f = mod("run_e5f_financing_factorial")
    if Path(f.ROOT).resolve() != exp: raise RuntimeError(f"factorial ROOT resolved to {f.ROOT}, expected {exp}")
    f.install(frozen)
    rental = mod("run_e5f_native_rental_access_diagnostic"); rental.install_paths(frozen)
    native = mod("run_e5f_native_financing_diagnostic"); native.install_paths(frozen)
    mod("run_e5f_native_income_cohort_diagnostic").install_paths(frozen)
    check_sys_path(exp)
    for name in (CORE + ".parameters", CORE + ".solver", "run_e5f_matched_pf_smoke", "run_e5f_independent_numerical_audit", "build_e5f_native_financing_report"): mod(name)
    check_sys_path(exp)
    return f, rental, native


def import_origins(m: dict[str, Any]) -> dict[str, str]:
    """Every project module must come from the copied frozen core or the copied allowed runtime helpers."""
    exp = Path(m["experiment_root"]).resolve(); frozen = Path(m["frozen_copy"]).resolve() / "code/model"
    helpers, wrapper, rows = exp / "code/model/tools", Path(m["wrapper"]["path"]).resolve(), {}
    for name, module in list(sys.modules.items()):
        file = getattr(module, "__file__", None)
        if not file: continue
        p = Path(file).resolve()
        if "Fertility_Spring26" not in str(p) and not under(p, exp): continue
        rows[name] = str(p); top = name.split(".")[0]
        if name in ("__main__", "__mp_main__") or top == wrapper.stem: ok = p == wrapper
        elif top == CORE: ok = under(p, frozen / CORE)
        elif top in RUNTIME_HELPERS: ok = p == helpers / (top + ".py")
        else: ok = under(p, frozen)
        if not ok: raise RuntimeError(f"import origin violation: {name} loaded from {p}")
    for req in (CORE + ".solver", CORE + ".parameters", *RUNTIME_HELPERS):
        if req not in rows: raise RuntimeError(f"required module not loaded: {req}")
    check_sys_path(exp)
    return rows


def old_rows(summary: dict[str, Any], phi: float, lam: float, cap: float) -> list[dict[str, Any]]:
    if summary.get("status") != "complete" or summary.get("design") != "dose": raise RuntimeError("old summary is not a complete dose summary")
    rows = [r for r in summary["cases"] if (float(r["phi"]), float(r["lambda"]), float(r["rental_cap"])) == (phi, lam, cap)]
    if not rows: raise RuntimeError(f"no old receipt for phi={phi} lambda={lam} cap={cap}")
    return rows


def reconcile_contract(row: dict[str, Any], x: Any, m: dict[str, Any], plan: dict[str, Any]) -> dict[str, dict[str, Any]]:
    c, P = row["contract"], x["parameters"]; fp = getattr(P, "income_candidate_fingerprint", None)
    checks = {"checkpoint_path": (c["checkpoint"], m["checkpoint"]["path"]), "checkpoint_sha256": (c["checkpoint_sha256"], m["checkpoint"]["sha256"]),
              "source_root": (c["source_root"], m["frozen_root"]), "population_source": (c["population_source"], "saved_evaluation"),
              "family": (c["family"], "refit_new_income"), "price": (c["price"], np.asarray(x["evaluation"].policy.price).tolist()),
              "initial_population_shape": (c["initial_population_shape"], list(np.asarray(x["stationary_g_pre"]).shape)),
              "candidate_fingerprint": (c["candidate_fingerprint"], fp), "plan_candidate_fingerprint": (plan["candidate_payload_fingerprint"], fp),
              "plan_source_root": (plan["source_root"], m["frozen_root"]), "preferences": (row["preferences"], "family checkpoint unchanged"),
              "population": (row["population"], "common contract pre-choice mass"), "entry_rule": (row["cohort"]["entry_rule"], "native conditional entrant for each process"),
              "mortgage_access_label": (row["mortgage_access_label"], "native joint deposit-and-collateral access")}
    out = {k: {"old": a, "expected": b, "ok": a == b} for k, (a, b) in checks.items()}
    bad = [k for k, v in out.items() if not v["ok"]]
    if bad: raise RuntimeError(f"contract reconciliation failed for {row['label']}: " + ", ".join(bad))
    return out


def compare_scalars(new: dict[str, Any], old: dict[str, Any], keys: tuple[str, ...]) -> dict[str, dict[str, Any]]:
    rows = {}
    for k in keys:
        a, b = new.get(k), old.get(k)
        if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool) and not isinstance(b, bool):
            d = abs(float(a) - float(b)); rows[k] = {"new": a, "old": b, "abs_diff": d, "within": bool(d <= TOL)}
        else: rows[k] = {"new": a, "old": b, "abs_diff": None, "within": a == b}
    return rows


def reproduce(new: dict[str, Any], rows: list[dict[str, Any]]) -> dict[str, Any]:
    per, passed = [], True
    for r in rows:
        metrics, cohort = compare_scalars(new["metrics"], r["metrics"], METRIC_KEYS), compare_scalars(new["cohort"], r["cohort"], COHORT_KEYS)
        ok = all(v["within"] for v in {**metrics, **cohort}.values()); passed &= ok
        per.append({"old_label": r["label"], "metrics": metrics, "cohort": cohort, "metrics_and_cohort_within_tolerance": ok,
                    "budget_informational": compare_scalars(new["budget"], r["budget"], tuple(r["budget"])),
                    "policy_audit_informational": compare_scalars(new["policy_audit"], r["policy_audit"], AUDIT_KEYS)})
    return {"tolerance": TOL, "gate": "scalar metrics and cohort summaries versus every matching old receipt", "passed": bool(passed), "rows": per}


def hash_tree(root: Path) -> dict[str, str]:
    return {str(p.relative_to(root)): sha(p) for p in sorted(root.rglob("*")) if p.is_file()}


def run_case(a: argparse.Namespace, m: dict[str, Any]) -> None:
    phi, lam, cap, label = CASES[a.case]; out = Path(a.results) / label
    if out.exists(): raise RuntimeError(f"case output already exists (run_case requires a fresh directory): {out}")
    if a.deadline_epoch and time.time() > a.deadline_epoch: raise RuntimeError("total time budget exhausted before case start")
    hashes = manifest_checks(m); f, rental, _ = install(m); origins = import_origins(m)
    summary, plan = read_json(m["summary"]["remote_path"]), read_json(m["plan"]["copy"])
    rows = old_rows(summary, phi, lam, cap)
    x = f.prepare_population(f.packet(Path(m["checkpoint"]["path"])), "saved_evaluation")
    contract = {r["label"]: reconcile_contract(r, x, m, plan) for r in rows}
    captured: dict[str, np.ndarray] = {}; calls: list[float] = []; original = rental.native_solve
    def observed(xx: Any, P: Any) -> Any:
        result = original(xx, P); ev = result[0]
        captured.update({n: np.array(getattr(ev.policy, n), copy=True) for n in f.POLICY_NAMES})
        captured.update({n: np.array(getattr(ev, n), copy=True) for n in EXTRA_ARRAYS})
        calls.append(time.time()); return result
    rental.native_solve = observed; started = time.time()
    try: receipt = f.run_case(x, phi, lam, cap, out, True)
    finally: rental.native_solve = original
    elapsed = time.time() - started
    if len(calls) != 1: raise RuntimeError(f"expected exactly one household solve, observed {len(calls)}")
    if receipt.get("status") != "completed" or receipt.get("standard_graphs", {}).get("count") != 17: raise RuntimeError(f"run_case receipt incomplete: {receipt.get('standard_graphs')}")
    if receipt.get("cohort", {}).get("status") != "completed": raise RuntimeError("cohort receipt incomplete")
    arrays = out / "policy_arrays.npz"; np.savez_compressed(arrays, **captured)
    repro = reproduce(receipt, rows)
    receipt.update({"status": "completed" if repro["passed"] else "failed_reproduction", "case": label, "case_index": a.case,
                    "old_labels": [r["label"] for r in rows], "reproduction": repro, "contract_reconciliation": contract, "hashes": hashes,
                    "import_origins": import_origins(m), "import_origins_before_solve": origins, "household_solves": len(calls), "solve_wall_seconds": elapsed,
                    "policy_arrays": {"path": str(arrays), "sha256": sha(arrays), "names": sorted(captured), "shapes": {n: list(v.shape) for n, v in captured.items()}},
                    "threads": {k: os.environ.get(k) for k in THREAD_VARS}, "job": os.environ.get("SLURM_JOB_ID"), "policy_identity_is_gate": False,
                    "manifest_sha256": sha(a.manifest), "scope": "fixed-price partial-equilibrium replay of an existing arm; no new target, parameter, entry, population, GE or calibration change"})
    receipt["file_hashes"] = hash_tree(out); write_json(out / "receipt.json", receipt)
    if not repro["passed"]: raise RuntimeError(f"{label} does not reproduce the old receipt within {TOL}; see {out / 'receipt.json'}")


def check_smoke(a: argparse.Namespace, m: dict[str, Any]) -> None:
    label = CASES[0][3]; out = Path(a.results) / label; r = read_json(out / "receipt.json"); problems = []
    if r.get("status") != "completed" or r.get("case") != label: problems.append("status")
    g = r.get("standard_graphs", {}); paths = [Path(p) for p in g.get("paths", [])]
    if g.get("count") != 17 or len(paths) != 17 or not all(p.is_file() for p in paths): problems.append("standard17")
    if len(list((out / "standard_diagnostics").glob("*.png"))) != 17: problems.append("png_count")
    arrays = out / "policy_arrays.npz"; pa = r.get("policy_arrays", {})
    if not arrays.is_file() or sha(arrays) != pa.get("sha256"): problems.append("policy_arrays")
    else:
        names = sorted(np.load(arrays).files)
        if names != pa.get("names") or not set(EXTRA_ARRAYS) <= set(names): problems.append("policy_array_names")
    if not r.get("hashes") or not all(v.get("ok") for v in r["hashes"].values()): problems.append("hashes")
    if not r.get("reproduction", {}).get("passed"): problems.append("reproduction")
    if r.get("cohort", {}).get("rows") != 17 or r.get("household_solves") != 1: problems.append("cohort_or_solve_count")
    if r.get("manifest_sha256") != sha(a.manifest): problems.append("manifest")
    stale = [rel for rel, h in r.get("file_hashes", {}).items() if not (out / rel).is_file() or sha(out / rel) != h]
    if stale: problems.append("file_hashes:" + ",".join(stale[:5]))
    if problems: raise RuntimeError("smoke receipt check failed: " + ", ".join(problems))
    write_json(Path(a.results) / "smoke_check.json", {"status": "verified", "case": label, "receipt_sha256": sha(out / "receipt.json"), "job": os.environ.get("SLURM_JOB_ID"), "checked": time.time()})


def compare(a: argparse.Namespace, m: dict[str, Any]) -> None:
    f, _, _ = install(m); results, labels = Path(a.results), [c[3] for c in CASES]; loaded = {}
    for label in labels:
        r = read_json(results / label / "receipt.json")
        if r.get("status") != "completed" or not r.get("reproduction", {}).get("passed"): raise RuntimeError(f"{label} is not a completed, reproduced case")
        array_path = results / label / "policy_arrays.npz"
        if sha(array_path) != r['policy_arrays']['sha256']: raise RuntimeError(f'array hash mismatch: {label}')
        loaded[label] = dict(np.load(array_path))
    names = sorted(set(f.POLICY_NAMES) | set(EXTRA_ARRAYS)); pairs = {}
    for label, arrays in loaded.items():
        if sorted(arrays) != names: raise RuntimeError(f'mandatory array set mismatch: {label}')
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            A, B, per = loaded[labels[i]], loaded[labels[j]], {}
            for n in names:
                if A[n].shape != B[n].shape: per[n] = {"shape_mismatch": [list(A[n].shape), list(B[n].shape)], "within_tolerance": False}; continue
                finite = np.isfinite(A[n]) & np.isfinite(B[n])
                d = np.abs(A[n][finite].astype(float) - B[n][finite].astype(float)); mx = float(d.max(initial=0.0))
                if not np.array_equal(A[n][~finite], B[n][~finite]): mx = float('inf')
                per[n] = {"max_abs_diff": mx, "count_above_tolerance": int(np.count_nonzero(d > TOL)), "identical": bool(np.array_equal(A[n], B[n])), "within_tolerance": bool(mx <= TOL)}
            pairs[f"{labels[i]}__vs__{labels[j]}"] = {"arrays": per, "policies_coincide_within_tolerance": all(per[n]["within_tolerance"] for n in names if n in f.POLICY_NAMES),
                                                      "distributions_coincide_within_tolerance": all(per[n]["within_tolerance"] for n in names if n in EXTRA_ARRAYS)}
    write_json(results / "policy_comparison.json", {"status": "completed", "question": "do household policies coincide across lambda in {0,1,5} at phi=.8, cap=6?",
                                                    "identity_is_gate": False, "tolerance": TOL, "array_names": names, "pairs": pairs,
                                                    "receipts": {l: sha(results / l / "receipt.json") for l in labels}, "job": os.environ.get("SLURM_JOB_ID")})


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=("verify", "case", "check-smoke", "compare"), required=True)
    p.add_argument("--manifest", type=Path, required=True); p.add_argument("--results", type=Path, required=True)
    p.add_argument("--case", type=int, choices=(0, 1, 2)); p.add_argument("--deadline-epoch", type=float); p.add_argument("--receipt", type=Path)
    a = p.parse_args(argv); m = read_json(a.manifest)
    def stop(signum: int, frame: Any) -> None: raise KeyboardInterrupt(f"terminated by signal {signum}")
    signal.signal(signal.SIGTERM, stop)
    a.results.mkdir(parents=True, exist_ok=True)
    try:
        if a.mode == "verify":
            hashes = manifest_checks(m); install(m)
            write_json(a.receipt or a.results / "verify.json", {"status": "verified", "hashes": hashes, "import_origins": import_origins(m), "sys_path": list(sys.path), "job": os.environ.get("SLURM_JOB_ID")})
        elif a.mode == "case":
            if a.case is None: raise RuntimeError("--case required")
            run_case(a, m)
        elif a.mode == "check-smoke": check_smoke(a, m)
        else: compare(a, m)
    except BaseException as exc:
        write_json(a.results / f"failure_{a.mode}_{'' if a.case is None else a.case}_{int(time.time())}.json", {"status": "failed", "mode": a.mode, "case": a.case, "error": f"{type(exc).__name__}: {exc}"}); raise
    return 0


if __name__ == "__main__": raise SystemExit(main())
