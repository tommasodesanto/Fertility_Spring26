"""Run the bounded 5x3, 9x3, and 15x3 native-entry income-grid comparison.

This driver is intentionally a launchable diagnostic, not a calibration or GE
driver.  It holds the checkpoint prices, wealth grid, fiscal objects, and all
non-income parameters fixed, explicitly rebuilds the income process with the
requested persistent-state count, then uses the existing native solver and
entry-cohort/17-graph diagnostics.  The launcher that accompanies this file is
dry-run only until the lead reviews its manifest.
"""
from __future__ import annotations

import argparse, copy, gzip, hashlib, importlib, inspect, json, pickle, signal, sys, time
from pathlib import Path
from typing import Any, Mapping
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
CHECKPOINT_SHA256 = "b3491eedcee6250cf94833067646d3e6463496cbf5a64bdcabe7b13bc7e89eb2"
SOLVER_SHA256 = "2992412586b81cef3a3e58d92191bb51f54d3f9cc600d7675bbadaed7d1682da"
PARAMS_SHA256 = "c0c1c18500fba069152659eaf588c3c895993cdcecb104cee7d6edca6bfae6a5"
FROZEN_REFIT_SOURCE_ROOT = "/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/corrected_initial_source_v2"
DEFAULT_CHECKPOINT = ROOT / "output/model/paper_baseline_sep14/replay_20260917/native_output/raw/repetition_02/initial_state.pkl.gz"
DEFAULT_SOURCE = ROOT / "tmp/paper_baseline_sep14/code/model"
DEFAULT_CANDIDATE = ROOT / "output/model/native_financing_diagnostic_20260919/earnings_candidate/candidate.json"
DEFAULT_OUTPUT = ROOT / "output/model/native_financing_diagnostic_20260919/specification_followup/income_grid_cohort_v1"
DEFAULT_BASELINE_SUMMARY = ROOT / "output/model/native_financing_diagnostic_20260919/overnight/final_mechanisms/refit_smoke/summary.json"
CASES = {"5x3": 5, "9x3": 9, "15x3": 15}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def install_paths(source: Path) -> None:
    for p in (source, source / "tools"):
        if str(p.resolve()) not in sys.path:
            sys.path.insert(0, str(p.resolve()))


def packet(path: Path) -> Any:
    with gzip.open(path, "rb") as f:
        return pickle.load(f)


def load_candidate(path: Path) -> dict[str, Any]:
    x = json.loads(path.read_text())
    if x.get("period_years") != 4 or x.get("permanent_types") is not False:
        raise ValueError("candidate must be the four-year no-permanent-type artifact")
    return x


def load_baseline_summary(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    summary = json.loads(path.read_text())
    for case in summary.get("cases", []):
        if (case.get("lambda") == 0.0 and case.get("phi") == 0.8 and case.get("rental_cap") == 6.0
                and case.get("cohort", {}).get("status") == "completed"):
            return case["contract"], case["cohort"]
    raise ValueError("baseline summary lacks retained refit lambda=0, phi=.8, cap=6 cohort")


def apply_grid_income(P: Any, candidate: Mapping[str, Any], persistent_states: int) -> tuple[Any, dict[str, Any]]:
    if persistent_states not in CASES.values():
        raise ValueError("persistent state count must be one of 5, 9, 15")
    from build_persistent_transitory_income_candidate import build_persistent_transitory_income_candidate
    annual = candidate["annual_coefficients_recovered_from_nested_fitted_covariances"]
    mapping = candidate["four_year_diagnostic_mapping"]
    rho = float(annual["rho_annual"]); vp = float(annual["persistent_variance"])
    overrides, meta = build_persistent_transitory_income_candidate(
        rho_annual=rho,
        persistent_innovation_sd_annual=float(np.sqrt((1.0-rho*rho)*vp)),
        transitory_log_sd_period=float(mapping["transitory_log_sd_period"]),
        period_years=float(candidate["period_years"]),
        persistent_states=int(persistent_states),
    )
    if meta["transitory_states"] != 3 or meta["joint_states"] != 3*persistent_states:
        raise ValueError("constructor did not return the required persistent x 3 state shape")
    for key in ("z_grid", "z_weights", "Pi_z", "income_shock_persistence", "use_income_types",
                "income_type_transition", "permanent_income_levels_enabled", "permanent_income_log_variance"):
        setattr(P, key, copy.deepcopy(overrides[key]))
    importlib.import_module("intergen_eqscale_seq_optimized.parameters").build_debt_caps(P)
    return P, meta


def shape_for_process(x: Mapping[str, Any], joint_states: int) -> dict[str, Any]:
    out = dict(x)
    shape = list(np.asarray(x["stationary_g_pre"]).shape)
    if len(shape) != 7 or shape[4] != 15:
        raise ValueError(f"checkpoint stationary_g_pre has unexpected shape {shape}")
    shape[4] = int(joint_states)
    out["stationary_g_pre"] = np.zeros(tuple(shape), dtype=float)
    return out


def module_hashes(driver: Path, source: Path, candidate: Path) -> dict[str, str]:
    local_panel = source / "intergen_eqscale_seq_optimized" / "local_panel.py"
    return {"driver": sha256(driver), "constructor": sha256(ROOT / "code/model/tools/build_persistent_transitory_income_candidate.py"),
            "candidate": sha256(candidate), "local_panel": sha256(local_panel)}


def verify_frozen_source(source: Path) -> dict[str, str]:
    root = source.parent.parent
    launch = json.loads((root / "launch_manifest.json").read_text())
    manifest_path = root / "source_manifest.json"
    if sha256(manifest_path) != launch["source_manifest_sha256"]:
        raise ValueError("source manifest checksum differs from launch")
    checked = json.loads(manifest_path.read_text())
    for rel, expected in checked.items():
        path = source / rel
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"staged source hash mismatch: {rel}")
    for rel, expected in (("intergen_eqscale_seq_optimized/solver.py", SOLVER_SHA256),
                          ("intergen_eqscale_seq_optimized/parameters.py", PARAMS_SHA256)):
        if sha256(source / rel) != expected:
            raise ValueError(f"frozen core hash mismatch: {rel}")
    for name, pin in launch["inputs"].items():
        if sha256(root / pin["copy"]) != pin["sha256"]:
            raise ValueError(f"input checksum differs from launch: {name}")
    if sha256(Path(__file__)) != launch["driver_sha256"]:
        raise ValueError("driver checksum mismatch")
    return checked


def income_payload_fingerprint(P: Any) -> str:
    h = hashlib.sha256()
    for key in ("z_grid", "z_weights", "Pi_z"):
        value = np.ascontiguousarray(np.asarray(getattr(P, key)))
        h.update(key.encode()); h.update(str(value.dtype).encode()); h.update(str(value.shape).encode()); h.update(value.tobytes())
    return h.hexdigest()


def module_origins() -> dict[str, str]:
    names = ("run_e5f_native_income_cohort_diagnostic", "run_e5f_native_rental_access_diagnostic",
             "run_e5f_matched_pf_smoke", "intergen_eqscale_seq_optimized.local_panel", "intergen_eqscale_seq_optimized.solver",
             "intergen_eqscale_seq_optimized.parameters", "build_persistent_transitory_income_candidate")
    out = {}
    for name in names:
        try:
            out[name] = str(Path(inspect.getfile(importlib.import_module(name))).resolve())
        except (ImportError, TypeError):
            out[name] = "unavailable"
    return out


def parameter_fingerprint(P: Any) -> str:
    """Stable receipt hash for the complete copied parameter object."""
    rows = {}
    for key, value in sorted(vars(P).items()):
        if isinstance(value, np.ndarray):
            rows[key] = {"dtype": str(value.dtype), "shape": list(value.shape), "sha256": hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()}
        else:
            try:
                json.dumps(value)
                rows[key] = value
            except TypeError:
                rows[key] = repr(value)
    return hashlib.sha256(json.dumps(rows, sort_keys=True, default=repr).encode()).hexdigest()


def run_with_alarm(fn: Any, seconds: float) -> Any:
    def alarm(_signum: int, _frame: Any) -> None:
        raise TimeoutError(f"case exceeded {seconds:.0f}-second wall-clock alarm")
    previous = signal.signal(signal.SIGALRM, alarm)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        return fn()
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, previous)


def run_case(x: Mapping[str, Any], candidate: Mapping[str, Any], case: str, out: Path,
             source: Path, checkpoint_policy: Any, baseline_cohort: Mapping[str, Any], deadline: float) -> dict[str, Any]:
    import run_e5f_native_income_cohort_diagnostic as base
    rental = importlib.import_module("run_e5f_native_rental_access_diagnostic")
    primitive = importlib.import_module("run_e5f_matched_pf_smoke")
    model = importlib.import_module("intergen_eqscale_seq_optimized.solver")
    P, meta = apply_grid_income(copy.deepcopy(x["parameters"]), candidate, CASES[case])
    x_case = shape_for_process(x, meta["joint_states"])
    x_case["stationary_g_pre"] = base.native_entry_cohort(P, np.asarray(x["b_grid"]), tuple(np.asarray(x_case["stationary_g_pre"]).shape))
    ev, budget, shared, solver_model = rental.native_solve(x_case, P)
    rental.gates(ev)
    if float(budget.get("budget_excess_mass", np.inf)) > 2e-10 or float(budget.get("maximum_occupied_excess", np.inf)) > 1e-9:
        raise ValueError(f"budget gate failed for {case}: {budget}")
    if time.monotonic() > deadline:
        raise TimeoutError(f"case {case} exceeded 900-second alarm")
    baseline_policy_gap = None
    if case == "5x3":
        from run_e5f_native_income_cohort_diagnostic import policy_arrays
        gaps = {}
        for name, got in policy_arrays(ev.policy).items():
            np.testing.assert_allclose(np.asarray(got), np.asarray(policy_arrays(checkpoint_policy)[name]), atol=1e-10, rtol=0)
            left, right = np.asarray(got, dtype=float), np.asarray(policy_arrays(checkpoint_policy)[name], dtype=float)
            finite = np.isfinite(left) & np.isfinite(right)
            if not np.array_equal(left[~finite], right[~finite], equal_nan=True):
                raise ValueError(f"nonfinite baseline policy mismatch: {name}")
            gaps[name] = float(np.abs(left[finite] - right[finite]).max(initial=0.0))
        baseline_policy_gap = max(gaps.values())
    case_out = out / case; case_out.mkdir(parents=True, exist_ok=True)
    base.save_full_arrays(case_out / "full_policy_arrays.npz", ev, x_case["stationary_g_pre"])
    receipt = base.run_cohort(x_case, P, ev.policy, case, case_out)
    # Standard graphs consume the normalized lifetime cohort distribution.
    cohort = np.load(case_out / "cohort_arrays.npz")
    x_graph = dict(x_case); x_graph["stationary_g_pre"] = cohort["lifetime_g_pre"]
    graph_shared = solver_model.precompute_shared(P, np.asarray(x["b_grid"]))
    graph_policy = primitive.pf.calendar.policy_from_solution(ev.policy, np.asarray(ev.policy.price), P, np.asarray(x["b_grid"]), graph_shared)
    graph_ev = primitive.pf.calendar.evaluate_period(np.asarray(ev.policy.price), x_graph["stationary_g_pre"], P, np.asarray(x["b_grid"]), graph_shared, primitive.pf.calendar.SolveCounter(), supply_rule=x["supply_rule"], supplied_policy=graph_policy)
    rental.gates(graph_ev)
    graphs = rental.standard_graphs(x_graph, P, graph_ev, graph_shared, solver_model, case_out)
    if graphs.get("status") != "completed" or graphs.get("count") != 17:
        raise RuntimeError(f"standard 17-graph gate failed for {case}: {graphs}")
    np.testing.assert_allclose(np.asarray(graph_ev.g_pre), np.asarray(x_graph["stationary_g_pre"]), atol=1e-10, rtol=0)
    entry = np.load(case_out / "initial_native_entry.npz")
    g = np.asarray(entry["g_pre"])
    # Conditional distributions are explicit because entry wealth can change with the grid.
    income_mass = g.sum(axis=(0, 1, 2, 3, 5, 6))
    wealth_marginal = g.sum(axis=(1, 2, 3, 4, 5, 6))
    wealth_by_income = g.sum(axis=(1, 2, 3, 5, 6))
    wealth_given_income = wealth_by_income / np.maximum(wealth_by_income.sum(axis=0, keepdims=True), 1e-300)
    income_given_wealth = wealth_by_income / np.maximum(wealth_by_income.sum(axis=1, keepdims=True), 1e-300)
    np.savez_compressed(case_out / "entry_marginal_conditional_arrays.npz", income_marginal=income_mass, wealth_marginal=wealth_marginal,
                        wealth_by_income=wealth_by_income, wealth_given_income=wealth_given_income,
                        income_given_wealth=income_given_wealth)
    (case_out / "entry_distribution_summary.json").write_text(json.dumps({
        "income_marginal": income_mass.tolist(), "wealth_marginal": wealth_marginal.tolist(),
        "wealth_given_income": wealth_given_income.tolist(), "income_given_wealth": income_given_wealth.tolist(),
        "interpretation": "native conditional entry distribution; grid changes can alter both marginal and conditional wealth/income composition"
    }, indent=2) + "\n")
    receipt.update({"persistent_states": CASES[case], "transitory_states": 3, "joint_income_states": meta["joint_states"],
                    "price": np.asarray(ev.policy.price).tolist(), "budget": budget,
                    "entry_distribution": "native conditional wealth/income rule rebuilt for each grid; no cross-grid bitwise equality asserted",
                    "entry_distribution_summary": "entry_distribution_summary.json (marginal and conditional wealth/income arrays)",
                    "full_parameter_fingerprint": parameter_fingerprint(P),
                    "income_payload_fingerprint": income_payload_fingerprint(P),
                    "interpretation": "fixed-price entry-cohort specification comparison; marginal and conditional entry distributions may change; not a pure policy effect",
                    "baseline_policy_max_abs_gap": baseline_policy_gap, "standard_graphs": graphs, "household_solves": 1})
    if case == "5x3":
        baseline_cohort_gap = 0.0
        for key in ("cumulative_explicit_births_per_initial_household", "first_births_per_initial_household", "first_birth_mean_age"):
            gap = abs(float(receipt[key]) - float(baseline_cohort[key]))
            np.testing.assert_allclose(float(receipt[key]), float(baseline_cohort[key]), atol=1e-10, rtol=0)
            baseline_cohort_gap = max(baseline_cohort_gap, gap)
        receipt["baseline_cohort_max_abs_gap"] = baseline_cohort_gap
    if time.monotonic() > deadline:
        raise TimeoutError(f"case {case} exceeded 900-second alarm")
    (case_out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=("smoke", "production"), default="smoke")
    ap.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    ap.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE)
    ap.add_argument("--candidate-json", type=Path, default=DEFAULT_CANDIDATE)
    ap.add_argument("--baseline-summary", type=Path, required=True)
    ap.add_argument("--baseline-summary-sha256", required=True)
    ap.add_argument("--smoke-receipt", type=Path)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    a = ap.parse_args(argv)
    if a.mode == "production" and a.smoke_receipt is None:
        raise SystemExit("production requires --smoke-receipt from the completed 5x3 smoke")
    if a.source_root.name != "model" or a.source_root.parent.name != "code":
        raise SystemExit("--source-root must be the normalized code/model tree")
    if a.output.exists() and any(a.output.iterdir()):
        raise SystemExit(f"refusing nonempty output directory: {a.output}")
    if sha256(a.checkpoint) != CHECKPOINT_SHA256:
        raise SystemExit("checkpoint hash is not the frozen b3491e artifact")
    install_paths(a.source_root)
    frozen_source_hashes = verify_frozen_source(a.source_root)
    before = module_hashes(Path(__file__), a.source_root, a.candidate_json)
    candidate = load_candidate(a.candidate_json)
    baseline_contract, baseline_cohort = load_baseline_summary(a.baseline_summary)
    if sha256(a.baseline_summary) != a.baseline_summary_sha256:
        raise SystemExit("baseline summary hash mismatch")
    if baseline_contract.get("checkpoint_sha256") != CHECKPOINT_SHA256:
        raise SystemExit("retained refit summary is not pinned to b3491e checkpoint")
    if int(candidate.get("state_count", -1)) != 15 or int(candidate.get("transitory_states", -1)) != 3:
        raise SystemExit("candidate metadata is not the retained 15-state, 3-transitory artifact")
    if candidate.get("source_sha256", {}).get("constructor.py") != before["constructor"]:
        raise SystemExit("candidate constructor fingerprint mismatch")
    if a.mode == "production":
        smoke = json.loads(a.smoke_receipt.read_text())
        if smoke.get("status") != "completed" or smoke.get("mode") != "smoke" or smoke.get("checkpoint_sha256") != CHECKPOINT_SHA256:
            raise SystemExit("smoke receipt is incomplete or has mismatched checkpoint")
        if smoke.get("source_manifest_hashes") != frozen_source_hashes or smoke.get("baseline_summary_sha256") != a.baseline_summary_sha256:
            raise SystemExit("smoke source/summary contract mismatch")
        if smoke.get("source_hashes_before") != before or smoke.get("source_hashes_after") != before:
            raise SystemExit("smoke receipt source hashes do not match current immutable inputs")
        if len(smoke.get("results", [])) != 1:
            raise SystemExit("smoke must contain exactly one completed baseline")
        baseline_result = smoke["results"][0]
        if baseline_result.get("status") != "completed" or baseline_result.get("household_solves") != 1:
            raise SystemExit("smoke baseline result incomplete")
        for key in ("baseline_policy_max_abs_gap", "baseline_cohort_max_abs_gap"):
            if baseline_result.get(key) is None or baseline_result[key] > 1e-10:
                raise SystemExit("smoke baseline reproduction failed: " + key)
        graph = baseline_result.get("standard_graphs", {})
        if graph.get("count") != 17 or len(graph.get("paths", [])) != 17 or not all(Path(p).is_file() for p in graph["paths"]):
            raise SystemExit("smoke standard plot packet missing")
    import run_e5f_native_income_cohort_diagnostic as base
    primitive = importlib.import_module("run_e5f_matched_pf_smoke"); primitive.pf.transition.configure_sequential_model()
    x = packet(a.checkpoint)
    if list(np.asarray(x["stationary_g_pre"]).shape) != list(baseline_contract.get("initial_population_shape", [])):
        raise SystemExit("checkpoint population shape differs from retained refit contract")
    if not np.isfinite(np.asarray(x["evaluation"].policy.price)).all() or abs(float(np.asarray(x["evaluation"].policy.price)[0]) - float(baseline_contract["price"][0])) > 1e-10:
        raise SystemExit("checkpoint price differs from retained refit contract")
    contract = baseline_contract
    a.output.mkdir(parents=True, exist_ok=True)
    cases = ["5x3"] if a.mode == "smoke" else ["9x3", "15x3"]
    manifest = {"status": "planned", "mode": a.mode, "cases": cases,
                "total_household_solves": 1 if a.mode == "smoke" else 2, "chain_household_solves": 3, "per_case_seconds": 900, "total_seconds": 1050 if a.mode == "smoke" else 2250,
                "checkpoint_sha256": sha256(a.checkpoint), "frozen_refit_source_root": FROZEN_REFIT_SOURCE_ROOT,
                "baseline_summary_sha256": a.baseline_summary_sha256, "source_hashes_before": before, "source_manifest_hashes": frozen_source_hashes, "module_origins": module_origins(), "contract": contract,
                "fixed": ["checkpoint prices", "wealth grid", "all non-income parameters", "native entry rule"], "no": ["GE", "psi normalization", "calibration", "target fit"]}
    (a.output / "launch_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (a.output / "latest_progress.json").write_text(json.dumps({"status": "planned", "cases_completed": []}, indent=2) + "\n")
    origins = module_origins()
    for name in ("run_e5f_native_income_cohort_diagnostic", "run_e5f_native_rental_access_diagnostic", "run_e5f_matched_pf_smoke", "intergen_eqscale_seq_optimized.local_panel", "intergen_eqscale_seq_optimized.solver",
             "intergen_eqscale_seq_optimized.parameters", "build_persistent_transitory_income_candidate"):
        if origins.get(name, "unavailable") == "unavailable" or not origins[name].startswith(str(a.source_root.resolve())):
            raise SystemExit(f"helper module origin is not inside frozen source tree: {name} -> {origins.get(name)}")
    manifest["module_origins"] = origins
    checkpoint_policy = getattr(x["evaluation"], "policy")
    results = []
    for case in manifest["cases"]:
        started = time.monotonic()
        results.append(run_with_alarm(lambda: run_case(x, candidate, case, a.output, a.source_root, checkpoint_policy, baseline_cohort, started + 900.0), 900.0))
        for module_name, module in list(sys.modules.items()):
            file = getattr(module, "__file__", None)
            if file and "Fertility_Spring26" in str(Path(file).resolve()):
                if not Path(file).resolve().is_relative_to(a.source_root.resolve()):
                    raise RuntimeError(f"foreign project import: {module_name} -> {file}")
        if time.monotonic() - started > 900.0:
            raise TimeoutError(f"case {case} exceeded 900-second wall clock")
        (a.output / "latest_progress.json").write_text(json.dumps({"status": "running", "cases_completed": [r["arm"] for r in results]}, indent=2) + "\n")
    if a.mode == "smoke":
        manifest["baseline_gate"] = {"cohort_metric_max_abs_gap": results[0].get("baseline_cohort_max_abs_gap"), "tolerance": 1e-10}
    after = module_hashes(Path(__file__), a.source_root, a.candidate_json)
    if verify_frozen_source(a.source_root) != frozen_source_hashes:
        raise RuntimeError("staged source changed during computation")
    if before != after:
        raise RuntimeError("driver/constructor/candidate/local_panel changed during computation")
    manifest.update({"status": "completed", "source_hashes_after": after, "results": results})
    (a.output / "launch_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
