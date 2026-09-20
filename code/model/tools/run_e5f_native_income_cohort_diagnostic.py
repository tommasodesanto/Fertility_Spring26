"""Fixed-price entry-cohort comparison for the no-permanent-type income candidate.

This is a partial-equilibrium diagnostic.  It uses the saved original policy
arrays for the two native arms and solves only the two candidate-income arms.
All arms start from the same age-zero entry cohort with income integrated out
and then redrawn from each process's stationary weights.
"""
from __future__ import annotations

import argparse, copy, gzip, hashlib, importlib, json, pickle, sys, time
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Mapping
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CHECKPOINT = ROOT / "output/model/paper_baseline_sep14/replay_20260917/native_output/raw/repetition_02/initial_state.pkl.gz"
DEFAULT_SOURCE = ROOT / "tmp/paper_baseline_sep14/code/model"
DEFAULT_PRIOR = ROOT / "output/model/native_financing_diagnostic_20260919"
DEFAULT_CANDIDATE = DEFAULT_PRIOR / "earnings_candidate/candidate.json"
DEFAULT_OUTPUT = DEFAULT_PRIOR / "income_cohort"
DEFAULT_REPLAY = ROOT / "output/model/paper_baseline_sep14/replay_20260917"
POLICY_NAMES = ("V", "c_pol", "hR_pol", "bp_pol", "tenure_choice", "tenure_probs", "loc_probs", "fert_probs", "fert_value", "fert2_probs", "price")


def install_paths(source: Path) -> None:
    for p in (source, source / "tools"):
        if str(p.resolve()) not in sys.path:
            sys.path.insert(0, str(p.resolve()))


def packet(path: Path) -> Any:
    with gzip.open(path, "rb") as f:
        return pickle.load(f)


def load_json(path: Path) -> dict[str, Any]:
    x = json.loads(path.read_text())
    if x.get("status") != "diagnostic_candidate_only" or x.get("permanent_types") is not False:
        raise ValueError("candidate JSON is not the verified no-permanent-type artifact")
    builder = ROOT / "code/model/tools/build_persistent_transitory_income_candidate.py"
    h = hashlib.sha256(builder.read_bytes()).hexdigest()
    if x.get("source_sha256", {}).get("constructor.py") != h:
        raise ValueError("candidate constructor source hash mismatch")
    return x


def policy_arrays(policy: Any) -> dict[str, np.ndarray]:
    out = {}
    for name in POLICY_NAMES:
        value = getattr(policy, name, None) if not isinstance(policy, Mapping) else policy.get(name)
        if value is None:
            raise ValueError(f"missing policy array: {name}")
        out[name] = np.asarray(value)
    return out


def policy_from_arrays(path: Path) -> SimpleNamespace:
    z = np.load(path)
    return SimpleNamespace(**{name: np.asarray(z[name]) for name in POLICY_NAMES})


def apply_candidate_income(P: Any, candidate: Mapping[str, Any]) -> Any:
    from build_persistent_transitory_income_candidate import build_persistent_transitory_income_candidate
    rho = float(candidate["annual_coefficients_recovered_from_nested_fitted_covariances"]["rho_annual"])
    vp = float(candidate["annual_coefficients_recovered_from_nested_fitted_covariances"]["persistent_variance"])
    sd_eta = float(np.sqrt((1.0 - rho * rho) * vp))
    sd_eps = float(candidate["four_year_diagnostic_mapping"]["transitory_log_sd_period"])
    overrides, _ = build_persistent_transitory_income_candidate(
        rho_annual=rho, persistent_innovation_sd_annual=sd_eta,
        transitory_log_sd_period=sd_eps, period_years=float(candidate["period_years"]),
    )
    for key in ("z_grid", "z_weights", "Pi_z", "income_shock_persistence", "use_income_types", "income_type_transition", "permanent_income_levels_enabled", "permanent_income_log_variance"):
        setattr(P, key, copy.deepcopy(overrides[key]))
    importlib.import_module("intergen_eqscale_seq_optimized.parameters").build_debt_caps(P)
    return P


def standardized_entry_cohort(P: Any, grid: np.ndarray, z_weights: np.ndarray, base_P: Any, shape: tuple[int, ...]) -> np.ndarray:
    """Build age-zero mass with common non-income marginal and supplied z weights."""
    calendar = importlib.import_module("run_e5f_matched_pf_smoke").pf.calendar
    base = np.asarray(calendar.entrant_cohort(np.ones(int(P.I)), base_P, grid), dtype=float)
    if base.ndim != 6 or len(shape) != 7:
        raise ValueError("native entrant cohort must be 6D and stationary mass 7D")
    # entrant_cohort axes are (b, tenure, location, income, n, m); insert age.
    marginal = base.sum(axis=3, keepdims=True)
    seed6 = np.repeat(marginal, len(z_weights), axis=3)
    seed6 *= np.asarray(z_weights, dtype=float).reshape((1, 1, 1, -1, 1, 1))
    seed = np.zeros(shape, dtype=float)
    seed[:, :, :, 0, :, :, :] = seed6
    return seed


def synthetic_mass_flow(pre: np.ndarray, post: np.ndarray, deaths: np.ndarray, entry: np.ndarray) -> float:
    arrays = [np.asarray(v, dtype=float) for v in (pre, post, deaths, entry)]
    if any(np.any(v < -1e-12) for v in arrays):
        raise ValueError("negative mass in cohort flow")
    return float(arrays[0].sum() - arrays[1].sum() - arrays[2].sum() + arrays[3].sum())


def save_full_arrays(path: Path, ev: Any, g_pre: np.ndarray) -> None:
    vals = {"g_pre": np.asarray(g_pre), "g_post_fertility": np.asarray(ev.g_post_fertility),
            "g_current": np.asarray(ev.g_current), "births": np.asarray(ev.births)}
    vals.update(policy_arrays(ev.policy))
    np.savez_compressed(path, **vals)


def housing_stats(ev: Any, P: Any) -> tuple[float, float]:
    from build_e5f_native_financing_report import _physical_housing
    owner, rooms, _ = _physical_housing(
        {"g_current": np.asarray(ev.g_current), "hR_pol": np.asarray(ev.policy.hR_pol)},
        {"H_own": np.asarray(P.H_own)})
    return rooms, owner


def run_cohort(x: Mapping[str, Any], P: Any, policy: Any, label: str, out: Path) -> dict[str, Any]:
    out.mkdir(parents=True, exist_ok=True)
    primitive = importlib.import_module("run_e5f_matched_pf_smoke")
    transition = primitive.pf.transition
    calendar = primitive.pf.calendar
    calendar.apply_fertility = transition.apply_sequential_fertility
    calendar.advance_calendar_distribution = transition.advance_sequential_calendar_distribution
    grid = np.asarray(x["b_grid"])
    model = importlib.import_module("intergen_eqscale_seq_optimized.solver")
    shared = model.precompute_shared(P, grid)
    P._fert2_probs = np.asarray(policy.fert2_probs).copy()
    policy = calendar.policy_from_solution(policy, np.asarray(policy.price), P, grid, shared)
    g = standardized_entry_cohort(P, grid, np.asarray(P.z_weights), x["parameters"], tuple(np.asarray(x["stationary_g_pre"]).shape))
    price = np.asarray(x["evaluation"].policy.price)
    initial_mass = float(g.sum())
    if abs(initial_mass - 1.0) > 1e-10: raise ValueError("entry mass must equal one")
    rental = importlib.import_module("run_e5f_native_rental_access_diagnostic")
    audit = importlib.import_module("run_e5f_independent_numerical_audit")
    rows, births_total, first_total = [], 0.0, 0.0
    for age in range(int(P.J)):
        ev = calendar.evaluate_period(price, g, P, grid, shared, calendar.SolveCounter(), supply_rule=x["supply_rule"], supplied_policy=policy)
        if not np.allclose(ev.g_pre, g, atol=1e-10, rtol=0):
            # evaluate_period returns the calendar's gated copy.  Its only
            # allowed mutation is the native dead-node relocation; retaining
            # this as a hard failure keeps the common-entry comparison honest.
            delta = float(np.abs(np.asarray(ev.g_pre) - np.asarray(g)).sum())
            projected = float(getattr(ev, "projected_mass", np.nan))
            raise ValueError(
                f"standardized cohort is infeasible under {label} age {age}: "
                f"native dead-node projection mass={projected:.17g}, L1 change={delta:.17g}; "
                "common non-income entry marginal cannot be silently projected"
            )
        if abs(float(g.sum()) - float(g[:, :, :, age, ...].sum())) > 1e-10:
            raise ValueError("cohort occupies unexpected ages")
        rental.gates(ev)
        budget = primitive.dated_budget(ev, P, shared, grid, float(P.user_cost_rate * price[0]))
        if float(budget.get("budget_excess_mass", np.inf)) > 2e-10 or float(budget.get("maximum_occupied_excess", np.inf)) > 1e-9:
            raise ValueError(f"cohort budget gate failed at age {age}: {budget}")
        audit_result = audit.policy_array_audit({"evaluation": ev, "parameters": P, "b_grid": grid}, out)
        if audit_result["occupied_negative_steps"]: raise ValueError(f"occupied value gate failed at age {age}")
        first = float(np.asarray(ev.g_pre)[:, :, :, age, :, 0, :].sum() - np.asarray(ev.g_post_fertility)[:, :, :, age, :, 0, :].sum())
        births = float(np.asarray(ev.births).sum())
        current = float(np.asarray(ev.g_current).sum()); post_childless=float(np.asarray(ev.g_post_fertility)[:, :, :, age, :, 0, :].sum())
        rooms, owner = housing_stats(ev, P)
        rows.append({"arm": label, "age_index": age, "age_years": float(P.age_start + age * P.da), "surviving_mass": current, "explicit_birth_flow": births, "exact_first_birth_flow": first, "childless_share": post_childless / max(current, 1e-30), "mean_rooms": rooms, "ownership": owner})
        births_total += births; first_total += first
        nxt, _, deaths, residual = transition.advance_sequential_calendar_distribution(ev, np.zeros(int(P.I)), P, grid, shared)
        if abs(float(residual)) > 1e-10 or abs(synthetic_mass_flow(g, nxt, deaths, np.zeros(int(P.I)))) > 1e-10:
            raise ValueError(f"cohort mass gate failed at age {age}")
        (out / "heartbeat.json").write_text(json.dumps({"arm": label, "age": age, "status": "running"}) + "\n")
        if abs(float(nxt[:, :, :, 0, ...].sum())) > 1e-12: raise ValueError("unexpected new entry")
        g = nxt
        import csv
        with (out / "cohort_by_age.csv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    out.mkdir(parents=True, exist_ok=True)
    import csv
    with (out / "cohort_by_age.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    np.savez_compressed(out / "cohort_arrays.npz", g_final=g)
    receipt = {"status": "completed", "arm": label, "initial_mass": 1.0, "cumulative_explicit_births_per_initial_household": births_total, "first_births_per_initial_household": first_total, "first_birth_mean_age": (sum(r["age_years"] * r["exact_first_birth_flow"] for r in rows) / first_total if first_total else None), "rows": len(rows), "scope": "fixed-price entry-cohort partial equilibrium; no GE, empirical fit, or production closure"}
    (out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT); p.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE)
    p.add_argument("--prior-input", type=Path, default=DEFAULT_PRIOR); p.add_argument("--candidate-json", type=Path, default=DEFAULT_CANDIDATE); p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT); p.add_argument("--replay", type=Path, default=DEFAULT_REPLAY)
    p.add_argument("--mode", choices=("inspect", "candidate"), default="inspect"); a = p.parse_args(argv)
    install_paths(a.source_root)
    contract = importlib.import_module("run_e5f_native_financing_diagnostic").validate_contract(a.checkpoint, a.replay, a.source_root)
    primitive = importlib.import_module("run_e5f_matched_pf_smoke")
    primitive.pf.transition.configure_sequential_model()
    x = packet(a.checkpoint); candidate = load_json(a.candidate_json)
    if a.mode == "inspect":
        (a.output / "inspect.json").parent.mkdir(parents=True, exist_ok=True); (a.output / "inspect.json").write_text(json.dumps({"status": "ready", **contract, "candidate_states": candidate["state_count"], "saved_arms": ["baseline_01", "mortgage_only_03"], "candidate_arms": ["candidate_phi08", "candidate_phi1"]}, indent=2) + "\n"); return
    a.output.mkdir(parents=True, exist_ok=True)
    receipts = []
    for label, phi in (("baseline", .8), ("mortgage_only", 1.0), ("candidate_phi08", .8), ("candidate_phi1", 1.0)):
        (a.output / label).mkdir(parents=True, exist_ok=True)
        if label in {"baseline", "mortgage_only"}:
            P = copy.deepcopy(x["parameters"]); P.phi = np.full_like(np.asarray(P.phi, dtype=float), phi)
        else:
            P = apply_candidate_income(copy.deepcopy(x["parameters"]), candidate); P.phi = np.full_like(np.asarray(P.phi, dtype=float), phi)
        importlib.import_module("intergen_eqscale_seq_optimized.parameters").build_debt_caps(P)
        if label in {"baseline", "mortgage_only"}:
            policy_path = a.prior_input / "cases" / ("baseline_01" if label == "baseline" else "mortgage_only_03") / "arrays.npz"
            policy = policy_from_arrays(policy_path)
            np.savez_compressed(a.output / label / "saved_policy_arrays.npz", **policy_arrays(policy))
            receipts.append(run_cohort(x, P, policy, label, a.output / label))
            (a.output / "latest_completed.json").write_text(json.dumps(receipts, indent=2))
            continue
        rental = importlib.import_module("run_e5f_native_rental_access_diagnostic")
        x2 = dict(x); x2["stationary_g_pre"] = standardized_entry_cohort(P, np.asarray(x["b_grid"]), np.asarray(P.z_weights), x["parameters"], tuple(np.asarray(x["stationary_g_pre"]).shape))
        ev, budget, shared, model = rental.native_solve(x2, P)
        rental.gates(ev)
        policy = ev.policy
        (a.output / label).mkdir(parents=True, exist_ok=True)
        save_full_arrays(a.output / label / "full_policy_arrays.npz", ev, x2["stationary_g_pre"])
        graphs = rental.standard_graphs(x2, P, ev, shared, model, a.output / label)
        (a.output / label / "standard_graphs.json").write_text(json.dumps(graphs, indent=2))
        if graphs["status"] != "completed": raise RuntimeError(graphs)
        # Keep the original non-income entry marginal for every arm.
        receipts.append(run_cohort(x, P, policy, label, a.output / label))
        (a.output / "latest_completed.json").write_text(json.dumps(receipts, indent=2))
    import csv
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for receipt in receipts:
        label = receipt["arm"]
        with (a.output / label / "cohort_by_age.csv").open() as f: rows = list(csv.DictReader(f))
        ages = [float(r["age_years"]) for r in rows]
        for ax, key in zip(axes, ("explicit_birth_flow", "mean_rooms", "ownership")):
            ax.plot(ages, [float(r[key]) for r in rows], label=label)
            ax.set_title(key.replace("_", " ")); ax.set_xlabel("Age")
    axes[0].legend(fontsize=7); fig.tight_layout(); fig.savefig(a.output / "cohort_comparison.png", dpi=160); plt.close(fig)
    (a.output / "receipt.json").write_text(json.dumps({"status": "completed", "arms": receipts, "contract": contract, "candidate_sha256": hashlib.sha256(a.candidate_json.read_bytes()).hexdigest(), "entry": "common non-income native entrant marginal; income independently redrawn in every arm; no offspring entry", "scope": "fixed-price diagnostic; no recalibration or GE"}, indent=2))


if __name__ == "__main__": main()
