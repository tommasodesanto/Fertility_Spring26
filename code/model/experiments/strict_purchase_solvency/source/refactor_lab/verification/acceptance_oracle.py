"""Acceptance oracle: frozen reference gates applied to a lab solution (Torch).

This is the ONLY lab file that executes the frozen September 28 reporting and
audit stack. The normal path (`run.py`) never imports it. It loads the actual
fixed-reference driver `run_fixed_price.py` by path (SHA pinned to its plan's
`driver_sha256`) and repeats that driver's post-solve sequence verbatim, with
the single change that the lifecycle solution comes from the lab engine:

  stationary reconstruction of g_pre from the lab policy (checked against the
  saved `stationary_g_pre` inside exact_control) -> dated evaluation ->
  gates() [budget, purchase/transaction wealth, estate funding, fiscal
  certificate, policy arrays, operator nesting, mass accounting] -> fertility,
  housing and recent-parent observers -> 14 target rows -> 31 parameter rows ->
  exact_control [every numeric array vs the checkpoint packet, all 14 rows
  exact] -> baseline-state impact gates -> dated standard_diagnostics (17 PNGs,
  hash-equal to the reference) -> public-parameter immutability -> source pins.

Any failed requirement raises; there is no partial certificate.

Phases
  fixed-price  --lab-solution rep1/acceptance_solution.pkl.gz rep2/...   (each fully certified)
  ge-certify   --lab-solution ge/acceptance_solution.pkl.gz              (gates/tables/plots, no equivalence)
  ge-old       old engine GE at the same explicit start, then ge certificate
The lab `shared` object (its own precompute_shared output) is used in the
packet, so exact_control compares lab shared arrays with the checkpoint and
the nested path set must equal the retained control_1 set (113 paths).

Class resolution for the lab pickle: `refactor_lab.engine.*` (on PYTHONPATH).
The pickle is a trusted local acceptance artifact, never a normal input.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import importlib.util
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np

FROZEN_DRIVER_REL = "output/model/fixed_reference_economics_20260928/sources/fixed_price_v1/run_fixed_price.py"
FROZEN_DRIVER_SHA256 = "96d6923a252f57bc4d8c44fd6479b13f48ba217d74edf8ef629d120428b03b44"
MUTABLE_IO_FIELDS = {"native_inherited_distribution_evidence_dir"}


def sha(path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as s:
        for b in iter(lambda: s.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def load_frozen(root: Path):
    path = root / FROZEN_DRIVER_REL
    if sha(path) != FROZEN_DRIVER_SHA256:
        raise RuntimeError("Frozen fixed-reference driver differs from its plan pin")
    spec = importlib.util.spec_from_file_location("frozen_fixed_price", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["frozen_fixed_price"] = module
    spec.loader.exec_module(module)
    return module


def platform_info() -> dict:
    import platform, numba, numpy, scipy
    return dict(machine=platform.machine(), system=platform.platform(), python=platform.python_version(),
                numpy=numpy.__version__, numba=numba.__version__, scipy=scipy.__version__,
                numba_threads=numba.get_num_threads(), numba_cache_dir=__import__("os").environ.get("NUMBA_CACHE_DIR"))


RUN_FIELDS = {"native_inherited_distribution_evidence_dir", "eq_iter"}   # output path; solver iteration counter


def parameter_identity(fp, final: dict, reference_P, label: str) -> dict:
    """Every public field must equal the checkpoint, except RUN_FIELDS, or equal
    the value the existing stationary PAYGO rule derives from the checkpoint
    parameters (bind_initial_balanced_pension). Anything else fails."""
    import e5f_stationary_paygo as paygo_frozen   # frozen tools path (after authenticate)
    ref = {k: v for k, v in vars(reference_P).items() if not k.startswith("_")}
    derived, _ = paygo_frozen.bind_initial_balanced_pension(copy.deepcopy(reference_P), payroll_tax=float(reference_P.tau_pay))
    derived = vars(derived)
    final = {k: v for k, v in final.items() if not k.startswith("_")}
    rows, unexplained = {}, []
    for k in sorted(set(ref) | set(final)):
        if k in RUN_FIELDS:
            rows[k] = "run_or_output_field"
            continue
        a, b = fp.serialized(final.get(k)), fp.serialized(ref.get(k))
        if k in final and k in ref and a == b:
            continue
        if k in final and a == fp.serialized(derived.get(k)):
            rows[k] = "derived_stationary_paygo"
        else:
            unexplained.append(k)
    fp.require(not unexplained, f"{label}: parameter drift vs checkpoint: {unexplained}")
    return dict(label=label, identical_except=rows, unexplained=unexplained)


RETAINED_CONTROL_REL = "output/model/fixed_reference_economics_20260928/fixed_price_v1/control_1/control_arrays.json"


def authenticate(args):
    fp = load_frozen(args.root)
    args.out.mkdir(parents=True, exist_ok=False)
    return (fp,) + tuple(fp.authenticate(args.out))


def load_stage(path: Path) -> dict:
    with gzip.open(path, "rb") as stream:
        return pickle.load(stream)


def postsolve(fp, manifest, contract, objective, runtime, prepared, reference, stage, out: Path, mode: str,
              renewal_tolerance: float | None = None, comparison: dict | None = None) -> dict:
    """run_fixed_price.child post-solve sequence on a supplied solution.

    ORACLE BOUNDARY: policy_from_solution, reconstruct_stationary_pre_fertility,
    evaluate_period, gates, observers, scoring and standard_diagnostics are the
    frozen calendar/tool implementations; the lab supplies solution, shared
    (its own precompute_shared output), parameters, grid and price.
    mode='fixed-price': exact control of every nested path vs the checkpoint.
    mode='ge': same gates/tables/plots, market+fiscal required, no checkpoint equality.
    """
    require, rt = fp.require, prepared.rt
    cal = rt["primitive"].pf.calendar
    out.mkdir(parents=True, exist_ok=False)
    sol, sd, P = stage["solution"], stage["shared"], stage["parameters"]
    grid, price = np.asarray(stage["b_grid"]), np.asarray(stage["price"], dtype=float).reshape(-1)
    require(np.array_equal(grid, reference["b_grid"]), "Grid differs")
    pub = lambda Q: {k: v for k, v in vars(Q).items() if not k.startswith("_") and k not in MUTABLE_IO_FIELDS}
    baseline = mode == "baseline"      # original engine on this machine: gates required, equality recorded
    fixed = mode == "fixed-price" or baseline
    if fixed:
        require(fp.serialized(pub(P)) == fp.serialized(pub(reference["parameters"])), "Public parameters differ from checkpoint")
        require(np.array_equal(price, np.asarray(reference["solution"].p_eq).reshape(-1)), "Replay price differs")
    public_before = fp.serialized({k: v for k, v in vars(P).items() if not k.startswith("_")})
    require(P.native_due_stayer_credit and P.native_exact_inherited_distribution, "Reference DUE/inherited-state contract absent")
    actual = fp.actual_parameters(prepared, P, grid)
    require(len(actual) == 31, "Complete parameter table unavailable")
    params = copy.deepcopy(manifest["full_parameter_table"])
    for row in params:
        if fixed:
            require(float(actual[row["parameter"]]) == float(row["estimate"]), "Parameter mismatch: " + row["parameter"])
        else:
            row["reference_estimate"], row["estimate"] = row["estimate"], str(float(actual[row["parameter"]]))
    P._fert2_probs = sol.fert2_probs.copy()
    policy = cal.policy_from_solution(sol, price, P, grid, sd)
    pre, reconstruction = cal.reconstruct_stationary_pre_fertility(sol, policy, P, grid, sd)
    runtime.require_abs_gate(reconstruction["stationary_post_fertility_nesting_l1"], 5e-9, "Stationary reconstruction")
    runtime.require_abs_gate(reconstruction["stationary_feasibility_projection_mass"], 0., "Stationary projection")
    if fixed and not baseline and comparison is None:   # otherwise compared inside the nested census
        require(np.array_equal(pre, reference["stationary_g_pre"]), "stationary_g_pre differs from checkpoint")
    supply = cal.HousingSupplyRule("static-elastic", float(price[0]),
        float(P.H0[0]*(P.user_cost_rate*price[0]/P.r_bar[0])**P.xi_supply[0]), float(P.xi_supply[0]))
    ev = cal.evaluate_period(price, pre, P, grid, sd, cal.SolveCounter(), supply_rule=supply, supplied_policy=policy)
    packet = dict(parameters=P, b_grid=grid, shared=sd, solution=sol, evaluation=ev,
        stationary_g_pre=pre, supply_rule=supply, demographic_seed=reference.get("demographic_seed"))
    LAST_PACKET[0] = packet
    gates = fp.gates(packet, prepared, out, stationary=True)
    market = float(ev.relative_market_residual)
    if not fixed:
        runtime.require_abs_gate(market, 2e-4, "housing market")
        require(bool(sol.timings.get("strict_converged")), "GE strict gate failed")
        cert = gates.get("fiscal_certificate") or {}   # certify_initial_pension raises on failure
        require(cert.get("marginal_gate") is True and cert.get("fiscal_gate") is True, "Fiscal certificate failed")
    fertility = {q: rt["observe_initial_fertility"](ev, P, age_projection=q) for q in ("uniform_birth_time", "constant_post_cell")}
    housing = rt["observe_initial_housing_wealth"](ev, P, grid, sd, diagnostic_enabled=True,
        age_projection="uniform_within_age_cell", diagnostic_allow_family_proxies=True,
        include_wealth=True, include_birth_response=True)
    recent = rt["observe_recent_parent_flow"](ev, P, diagnostic_enabled=True, snapshot=rt["SNAPSHOT"],
        age_projection=rt["AGE_PROJECTION"], diagnostic_allow_residence_proxy=True,
        input_provenance=dict(case_id="refactor_lab_" + mode, reference_checkpoint_sha256=manifest["checkpoint"]["sha256"]))
    completed = float(rt["chain"].extract_moments(sol, P)["tfr"])
    fits = runtime.score_targets(objective, fertility, housing, recent["model_value"], completed)
    require(len(fits) == 14 and len(params) == 31, "All 14 fit rows and 31 parameter rows required")
    control = None
    native = sys.modules["e5f_current_transition_runtime"]
    retained = json.loads((Path(ORIG_ROOT[0]) / RETAINED_CONTROL_REL).read_text())
    if fixed and comparison is None and not baseline:
        control = fp.exact_control(reference, packet, fits, manifest, prepared, out)
        mine = json.loads((out / "control_arrays.json").read_text())
        require(set(mine["arrays"]) == set(retained["arrays"]) and len(mine["arrays"]) == 113,
                "Nested control path set differs from retained control_1 (113 paths)")
    elif fixed:
        against = comparison["packet"] if comparison else reference
        check = native.compare_arrays(against, packet)
        fp.write(out / "control_arrays.json", check)
        require(set(check["arrays"]) == set(retained["arrays"]) and len(check["arrays"]) == 113,
                "Nested control path set differs from retained control_1 (113 paths)")
        inexact = {k: v for k, v in check["arrays"].items()
                   if v.get("status") != "compared" or not v.get("exact") or not v.get("finite")}
        if comparison:   # lab vs SAME-MACHINE original-engine baseline: exact, finite, complete
            require(not inexact, "Lab differs from same-machine old baseline: " + ", ".join(list(inexact)[:15]))
            for name, rows_new, rows_old in (("fit", fits, comparison["fits"]), ("parameter", params, comparison["params"])):
                require(len(rows_new) == len(rows_old) == (14 if name == "fit" else 31), f"{name} row count differs")
                key = "moment" if name == "fit" else "parameter"
                require([r[key] for r in rows_new] == [r[key] for r in rows_old], f"{name} row names/order differ")
                for x, y in zip(rows_new, rows_old):
                    for col in y:
                        xa, ya = str(x.get(col, "")), str(y[col])
                        same = xa == ya
                        if not same:
                            try:
                                same = float(xa) == float(ya)
                            except ValueError:
                                same = False
                        require(same, f"{name} row differs from baseline: {y[key]}/{col}")
            control = dict(status="passed_vs_same_machine_old_baseline", arrays=len(check["arrays"]),
                           baseline=comparison["label"], torch_checkpoint_comparison="diagnostic only (see driver)")
        else:            # baseline vs Torch checkpoint: platform differences RECORDED, not relaxed or hidden
            control = dict(status="recorded_vs_torch_checkpoint", arrays=len(check["arrays"]),
                           nonexact_paths={k: {q: v.get(q) for q in ("status", "max_abs", "max_rel", "exact", "finite") if q in v}
                                           for k, v in inexact.items()},
                           fit_rows_vs_manifest=[dict(moment=x["moment"], model=x["model"], canonical=y["model"])
                                                 for x, y in zip(fits, manifest["full_target_table"])
                                                 if str(x["model"]) != str(y["model"])])
    if fixed and not baseline and comparison is None:
        impact_P = copy.deepcopy(P)
        impact_sd = rt["model"].precompute_shared(impact_P, grid)
        impact = cal.evaluate_period(price, reference["stationary_g_pre"], impact_P, grid, impact_sd,
            cal.SolveCounter(), supply_rule=supply, supplied_policy=policy)
        (out / "baseline_state_impact").mkdir()
        fp.gates(dict(parameters=impact_P, b_grid=grid, shared=impact_sd, evaluation=impact,
                      stationary_g_pre=reference["stationary_g_pre"]), prepared, out / "baseline_state_impact", stationary=False)
        for name in ("g_pre", "g_post_fertility", "g_current", "g_stay_distribution"):
            require(np.array_equal(getattr(impact, name), getattr(reference["evaluation"], name)), "Impact array differs: " + name)
    fp.table(out / "target_fit.csv", fits)
    fp.table(out / "parameters.csv", params)
    fp.write(out / "observers.json", cal.jsonable(dict(fertility=fertility, housing_wealth=housing, recent_parent=recent)))
    rt["audit"].standard_diagnostics(packet, out, validate_production_young=False)
    plots = sorted(q.name for q in (out / "standard_diagnostics").glob("*.png"))
    require(plots == sorted(manifest["standard_diagnostic_names"]), "Standard 17-plot set differs")
    plot_hashes = {name: sha(out / "standard_diagnostics" / name) for name in plots}
    if fixed and comparison is not None:
        require(sorted(comparison["plot_hashes"]) == plots, "Plot set differs from same-machine baseline")
        for name in plots:
            require(plot_hashes[name] == comparison["plot_hashes"][name], "Diagnostic differs from same-machine baseline: " + name)
    elif fixed and not baseline:
        for name in plots:
            require(plot_hashes[name] == manifest["artifact_hashes"]["standard_diagnostics/" + name],
                    "Diagnostic hash differs: " + name)
    elif baseline:
        control["plots_vs_torch_canonical"] = {n: plot_hashes[n] == manifest["artifact_hashes"]["standard_diagnostics/" + n] for n in plots}
    require(fp.serialized({k: v for k, v in vars(P).items() if not k.startswith("_")}) == public_before, "Public parameter mutated")
    runtime.verify_sources(dict(contract, objective=manifest["objective"]))
    ref_gap = float(reference["solution"].adult_entry_stationary_relative_gap)
    gap = float(sol.adult_entry_stationary_relative_gap)
    if fixed:
        renewal_class = "checked_by_exact_control"
    elif renewal_tolerance is None:
        renewal_class = "unclassified_no_tolerance_supplied"
    else:
        renewal_class = "consistent_with_reference" if abs(gap - ref_gap) <= renewal_tolerance else "deviates_from_reference"
    if baseline:
        status = "old_engine_same_machine_baseline_gates_passed"
    elif fixed and comparison is not None:
        status = "certified_fixed_price_replay_vs_same_machine_old_baseline"
    elif fixed:
        status = "certified_fixed_price_replay"
    elif renewal_class == "consistent_with_reference":
        status = "ge_gates_passed_renewal_consistent_not_equivalence"
    else:
        status = "ge_limited_renewal_" + renewal_class
    receipt = dict(status=status, plot_hashes=plot_hashes, platform=platform_info(),
        mode=mode, oracle_boundary=postsolve.__doc__.split("\n\n")[1], price=price.tolist(),
        market_relative_residual=market, completed_fertility=completed, rows_fit=len(fits), rows_parameters=len(params),
        target_comparison_loss=sum(float(r["loss_contribution"]) for r in fits if r["loss_contribution"] != ""),
        renewal=dict(closure=f"{P.adult_entry_clock}/{P.population_closure}, unchanged; no psi normalization",
                     relative_gap=gap, signed_residual=float(sol.adult_entry_stationary_residual),
                     reference_relative_gap=ref_gap, tolerance=renewal_tolerance,
                     tolerance_source="lead-selected DIAGNOSTIC threshold for this comparison; the frozen reference encodes no renewal tolerance (adult_entry_gate records entry_E and residual only)",
                     classification=renewal_class,
                     inherited_reference_gate=manifest["inherited_gates"]["adult_entry_gate"]),
        control=control, gates=gates, standard_plot_count=len(plots), plots_hash_checked=fixed)
    fp.write(out / "receipt.json", cal.jsonable(receipt))
    return receipt


ORIG_ROOT = [None]


def load_baseline(root: Path, directory: Path, pin: str) -> dict:
    """Accept a same-machine baseline only via baseline_identity.verify (pin, original
    engine bytes, checkpoint/manifest/driver hashes, artifact hashes, current runtime)."""
    import csv
    from refactor_lab.verification import baseline_identity
    r = baseline_identity.verify(root, directory, pin)
    rows = lambda name: list(csv.DictReader((directory / "certificate" / name).open()))
    return dict(packet=load_stage(directory / "baseline_packet.pkl.gz"), fits=rows("target_fit.csv"),
                params=rows("parameters.csv"), plot_hashes=r["plots"], label=str(directory), receipt=r)


def fixed_price(args) -> dict:
    fp, *auth = authenticate(args)
    ORIG_ROOT[0] = str(args.root)
    comparison = load_baseline(args.root, args.comparison_reference, args.comparison_pin) if args.comparison_reference else None
    receipts = {}
    for k, path in enumerate(args.lab_solution, start=1):
        receipts[f"rep{k}"] = postsolve(fp, *auth, load_stage(path), args.out / f"rep{k}", "fixed-price",
                                        comparison=comparison)
    fp.write(args.out / "summary.json", dict(status="certified_fixed_price_replay_all_repetitions",
                                             repetitions=len(receipts), lab_solutions=[str(p) for p in args.lab_solution],
                                             lab_solution_sha256=[sha(p) for p in args.lab_solution]))
    return receipts


def old_fixed_price(args) -> dict:
    """ORIGINAL engine at the reference price and parameters on this machine (baseline)."""
    fp, manifest, contract, objective, runtime, prepared, reference = authenticate(args)
    ORIG_ROOT[0] = str(args.root)
    model = prepared.rt["model"]
    P = copy.deepcopy(reference["parameters"])
    P.native_inherited_distribution_evidence_dir = str(args.out / "inherited_state_failures")
    grid = np.asarray(reference["b_grid"]).copy()
    price = np.asarray(reference["solution"].p_eq, dtype=float).reshape(-1).copy()
    sd = model.precompute_shared(P, grid)
    t = time.perf_counter()
    sol = model.solve_markov_income_at_prices(price, P, grid, SD=sd, verbose=False, fast_stats=False)
    solve_s = time.perf_counter() - t
    stage = dict(solution=sol, shared=sd, parameters=P, b_grid=grid, price=price)
    cert = postsolve(fp, manifest, contract, objective, runtime, prepared, reference, stage, args.out / "certificate", "baseline")
    # The 113-path census compares evaluation/stationary_g_pre too; store the certified packet objects.
    with gzip.open(args.out / "baseline_packet.pkl.gz", "wb", compresslevel=1) as stream:
        pickle.dump(LAST_PACKET[0], stream, protocol=5)
    flat = {k: v for k, v in vars(sol).items() if isinstance(v, np.ndarray) and v.dtype != object}
    flat.update({"shared." + k: v for k, v in vars(sd).items() if isinstance(v, np.ndarray) and v.dtype != object})
    np.savez(args.out / "solution_arrays.npz", **flat)          # same layout as run.save_stage, for compare.py
    from refactor_lab.verification import baseline_identity
    receipt = dict(baseline_identity.record(args.root, args.out), lifecycle_solve_seconds=solve_s,
                   certificate_status=cert["status"], vs_torch_checkpoint=cert["control"],
                   note="Torch-checkpoint differences are recorded platform differences, not relaxed gates")
    fp.write(args.out / "baseline_receipt.json", prepared.rt["primitive"].pf.calendar.jsonable(receipt))
    print("baseline_receipt_sha256", sha(args.out / "baseline_receipt.json"))
    return receipt


LAST_PACKET = [None]


def ge_certify(args) -> dict:
    fp, *auth = authenticate(args)
    ORIG_ROOT[0] = str(args.root)
    reference = auth[-1]
    stage = load_stage(args.lab_solution[0])
    initial = json.loads((Path(args.lab_solution[0]).parents[1] / "parameters_initial.json").read_text())
    ident = dict(initial=_serialized_identity(fp, initial, reference["parameters"], "lab initial"),
                 final=parameter_identity(fp, vars(stage["parameters"]), reference["parameters"], "lab final"))
    fp.write(args.out / "parameter_identity.json", ident)
    return postsolve(fp, *auth, stage, args.out / "certificate", "ge", args.renewal_tolerance)


def _serialized_identity(fp, initial_serialized: dict, reference_P, label: str) -> dict:
    """Initial inputs: already-serialized public fields must equal the checkpoint (RUN_FIELDS aside)."""
    ref = fp.serialized({k: v for k, v in vars(reference_P).items() if not k.startswith("_") and k not in RUN_FIELDS})
    mine = {k: v for k, v in initial_serialized.items() if not k.startswith("_") and k not in RUN_FIELDS}
    diff = sorted(k for k in set(ref) | set(mine) if ref.get(k) != mine.get(k))
    fp.require(not diff, f"{label}: input drift vs checkpoint: {diff}")
    return dict(label=label, identical=True, excluded=sorted(RUN_FIELDS))


def ge_old(args) -> dict:
    """Old-engine GE at the same explicit start, then the same GE certificate."""
    fp, manifest, contract, objective, runtime, prepared, reference = authenticate(args)
    ORIG_ROOT[0] = str(args.root)
    rt = prepared.rt
    model = rt["model"]
    P = copy.deepcopy(reference["parameters"])
    P.native_inherited_distribution_evidence_dir = str(args.out / "inherited_state_failures")
    grid = np.asarray(reference["b_grid"]).copy()
    start = float(args.price_factor) * np.asarray(reference["solution"].p_eq, dtype=float).reshape(-1)
    initial = _serialized_identity(fp, fp.serialized(vars(P)), reference["parameters"], "old initial")
    solve = lambda: rt["solve_balanced_initial_equilibrium"](model=model, parameters=P, b_grid=grid,
        initial_prices=start, payroll_tax=float(P.tau_pay), marginal_tolerance=1e-9, fiscal_tolerance=1e-6)
    from refactor_lab.verification.callcount import BudgetExceeded, CallCounter
    counter = CallCounter(model, max_lifecycle=args.max_lifecycle, deadline_seconds=args.solve_deadline) \
        if args.count_calls else None
    t = time.perf_counter()
    try:
        if counter is not None:
            with counter:
                sol, P2, prices, fiscal = solve()
        else:
            sol, P2, prices, fiscal = solve()
    except BudgetExceeded as exc:
        fp.write(args.out / "budget_failure.json", dict(status="failed_budget_not_a_certificate", reason=str(exc),
                 profiled=counter.summary() if counter else None))
        raise SystemExit(4)
    solve_s = time.perf_counter() - t
    price = np.asarray(prices, dtype=float).reshape(-1)
    sd = model.precompute_shared(P2, grid)
    arrays = {k: v for k, v in vars(sol).items() if isinstance(v, np.ndarray) and v.dtype != object}
    arrays.update({"shared." + k: v for k, v in vars(sd).items() if isinstance(v, np.ndarray) and v.dtype != object})
    np.savez(args.out / "solution_arrays.npz", **arrays)
    stage = dict(solution=sol, shared=sd, parameters=P2, b_grid=grid, price=price)
    ident = dict(initial=initial, final=parameter_identity(fp, vars(P2), reference["parameters"], "old final"))
    if args.lab_solution:   # final effective parameters: lab and old must agree field by field
        lab_P = load_stage(args.lab_solution[0])["parameters"]
        pick = lambda Q: fp.serialized({k: v for k, v in vars(Q).items() if not k.startswith("_") and k not in RUN_FIELDS})
        a, b = pick(lab_P), pick(P2)
        diff = sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))
        fp.require(not diff, f"lab/old final parameter mismatch: {diff}")
        ident["lab_vs_old_final"] = dict(identical=True, excluded=sorted(RUN_FIELDS))
    fp.write(args.out / "parameter_identity.json", ident)
    receipt = dict(engine="old_frozen", start=start.tolist(), price=price.tolist(), equilibrium_solve_stage_seconds=solve_s,
                   equilibrium_solve_stage_profiled=counter is not None,
                   engine_timings=sol.timings, fiscal=fiscal, profiled=counter.summary() if counter else None,
                   note="process includes frozen authentication/imports; not comparable to the lab process total")
    fp.write(args.out / "solve_receipt.json", rt["primitive"].pf.calendar.jsonable(receipt))
    receipt["certificate"] = postsolve(fp, manifest, contract, objective, runtime, prepared, reference,
                                       stage, args.out / "certificate", "ge", args.renewal_tolerance)["status"]
    return receipt


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("phase", choices=("fixed-price", "old-fixed-price", "ge-certify", "ge-old"))
    ap.add_argument("--comparison-reference", type=Path, help="same-machine old baseline directory")
    ap.add_argument("--comparison-pin", help="sha256 of its baseline_receipt.json")
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--lab-solution", type=Path, nargs="+")
    ap.add_argument("--count-calls", action="store_true")
    ap.add_argument("--max-lifecycle", type=int)
    ap.add_argument("--solve-deadline", type=float)
    ap.add_argument("--renewal-tolerance", type=float)
    ap.add_argument("--price-factor", type=float)
    a = ap.parse_args()
    if a.phase == "old-fixed-price":
        old_fixed_price(a)
        return
    if a.comparison_reference and not a.comparison_pin:
        raise SystemExit("--comparison-pin required with --comparison-reference")
    if a.phase in ("fixed-price", "ge-certify"):
        if not a.lab_solution:
            raise SystemExit("--lab-solution required")
        (fixed_price if a.phase == "fixed-price" else ge_certify)(a)
    else:
        if a.price_factor is None:
            raise SystemExit("--price-factor required (lead: 1.05)")
        ge_old(a)


if __name__ == "__main__":
    main()
