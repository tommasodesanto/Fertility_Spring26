#!/usr/bin/env python3
"""Six prescribed-price cases on one frozen candidate; no GE root is solved."""
from __future__ import annotations

import copy
import csv
import hashlib
import importlib.util
import json
import os
import signal
import sys
import time
import traceback
from pathlib import Path
from types import SimpleNamespace

for _thread_var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_thread_var] = "1"

import runtime_overlay
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
PSI = ROOT / "output/model/fixed_reference_economics_20260928/utility_floor_psi_v1"
OLD = ROOT / "output/model/fixed_reference_economics_20260928/utility_floor_round2_v1"
BASE = ROOT / "output/model/publication_refactor_20260929/grid_resolution_v1/credit053_v2/runner"
CONTRACT = json.loads((HERE / "contract.json").read_text())
BINDING = json.loads((HERE / "source_binding.json").read_text())


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False, default=str) + "\n")
    tmp.replace(path)


def table(path, rows):
    if not rows:
        raise RuntimeError(f"Cannot write empty table: {path}")
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


class CaseDeadline(TimeoutError):
    """Only the explicitly declared per-cell alarm may continue beyond q0."""


def _case_alarm(_signum, _frame):
    raise CaseDeadline("300-second fixed-price case deadline reached")


def authenticate_candidate(runtime_out=None):
    require(CONTRACT["candidate_root"].endswith("verified_global_20261001T0941NY_chain7_0173/ROOT"),
            "Frozen candidate identity changed")
    for row in BINDING["files"]:
        path = ROOT / row["path"]
        require(path.is_file() and sha(path) == row["sha256"], "Pinned input drift: " + str(path))
    source_pins_path = ROOT / BINDING["source_pins_path"]
    require(sha(source_pins_path) == BINDING["source_pins_sha256"], "Original source pins changed")
    pins = json.loads(source_pins_path.read_text())
    for rel, digest in pins.items():
        require(sha(ROOT / rel) == digest, "Original engine source changed: " + rel)

    sys.path.insert(0, str(OLD))
    import inputs
    import runner as native
    case = ROOT / CONTRACT["candidate_root"]
    params_rows = native.readtable(case / "parameters.csv")
    params = {r["parameter"]: float(r["estimate"]) for r in params_rows}
    require(len(params_rows) == 31 and params["psi_child"] == float(BINDING["candidate_psi"]),
            "Winner parameter count or psi differs from the binding")
    closure = json.loads((case / "closure.json").read_text())
    require(float(closure["price"]) == float(BINDING["candidate_price"]),
            "Verified winner price differs from the binding")
    fits = native.readtable(case / "target_fit.csv")
    require(len(fits) == 14 and abs(sum(float(row["loss_contribution"] or 0.0)
            for row in fits) - float(BINDING["candidate_loss"])) < 1e-8,
            "Winner target rows or base loss differ from the binding")
    winner_input = json.loads((ROOT / "output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/deployment/monitor_snapshot/verified_global_20261001T0941NY_chain7_0173/input_contract.json").read_text())
    require(winner_input["base_target_contract_sha256"] == BINDING["winner_target_contract_sha256"]
            and winner_input["entry"]["conditional_sha256"] == BINDING["winner_entry_conditional_sha256"],
            "Verified winner target or entry fingerprint differs")
    seed, bounds, _ = inputs.seed_and_bounds("floor_s0")
    bounds = dict(bounds)
    psi_config = json.loads((PSI / "plan.json").read_text())
    bounds["psi_child"] = tuple(psi_config["psi_bounds"])
    point = {k: params[k] for k in inputs.parameters("floor_s0")}
    point["psi_child"] = params["psi_child"]
    P, grid = inputs.proposal("floor_s0")
    P, entry = inputs.entry(P, grid, "nonnegative_mean")
    P = inputs.bind(P, point, bounds, "floor")
    sys.path.insert(0, str(ROOT / "code/model"))
    credit_source = ROOT / "output/model/fixed_reference_economics_20260928/credit_no_taper_v1/small_credit_v1/source"
    sys.path.insert(0, str(credit_source))
    from small_credit_lab import credit
    from refactor_lab.engine import solver
    credit.bind_engine_credit(P, "corrected", 0.0)
    natural = copy.deepcopy(P)
    natural.native_due_stayer_credit = False
    natural.native_solvency_credit = True
    natural.unsecured_credit_limit = None
    require(solver.validate_native_solvency_mode(natural), "Native lifetime-repayment-only mode failed its gate")

    # Reuse the candidate runner's frozen observer/reporting stack and require
    # that its executed preference/solver code matches the active refactor.
    sys.path.insert(0, str(BASE))
    spec = importlib.util.spec_from_file_location("mechanism_base_workflow", BASE / "run_comparison.py")
    base = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)
    sys.path.insert(0, str(OLD))
    import phase_b_pilot as ge
    runtime_out = Path(runtime_out) if runtime_out is not None else HERE / "runtime_context" / ("check_"+str(time.time_ns()))
    runtime_out.mkdir(parents=True, exist_ok=True)
    ctx = base.authored.context_from_bundle(SimpleNamespace(
        bundle=ROOT / "output/model/publication_refactor_20260929/local_export_v1/inputs",
        reference_root=ROOT, out=runtime_out))
    native.install_reporter_on_authored(base.authored)
    base.authored.authenticate_frozen(ctx)
    dims = (int(P.Nb), int(P.Nz))
    expected = native.expected_parameters(point, dimensions=dims, arm="floor")
    ctx.update(P=P, b_grid=grid, selected_d_bar=0.0, reference_psi=float(P.psi_child),
               expected_dimensions={"wealth_grid_nodes": dims[0], "income_states": dims[1]},
               expected_parameters=expected, free_coordinates=[],
               price_start=float(BINDING["candidate_price"]))

    actual = ctx["fp"].actual_parameters(ctx["prepared"], P, grid)
    ge.validate_parameter_estimates(ctx, native.PLAN["reference_parameter_table"], actual)
    candidate_rows = native.readtable(ROOT / CONTRACT["candidate_root"] / "parameters.csv")
    candidate_values = {r["parameter"]: float(r["estimate"]) for r in candidate_rows}
    require(set(actual) == set(candidate_values) and len(actual) == 31,
            "Native effective proposal does not match the full candidate table")
    for key, value in actual.items():
        require(abs(float(value) - candidate_values[key]) <= 2e-12 * max(1.0, abs(candidate_values[key])),
                "Actual effective parameter differs from frozen candidate: " + key)

    checked = solver.precompute_shared(P, grid)
    from small_credit_lab.engine import solver as live_solver
    live = live_solver.precompute_shared(P, grid)
    for key in ("h_bar", "c_bar", "g_bar", "alpha_flat", "psi_v", "escale_flat"):
        np.testing.assert_array_equal(getattr(checked, key), getattr(live, key))
    from small_credit_lab.engine import shared as es, child_preferences as ec, kernels as ek
    from refactor_lab.engine import shared as cs, child_preferences as cc, kernels as ck
    for current, executed in ((cs, es), (cc, ec), (ck, ek)):
        require(sha(current.__file__) == sha(executed.__file__), "Frozen overlay differs from active refactor: " + current.__name__)
    require(entry["conditional_sha256"] == json.loads((ROOT / "output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/deployment/monitor_snapshot/verified_global_20261001T0941NY_chain7_0173/input_contract.json").read_text())["entry"]["conditional_sha256"],
            "Candidate entry fingerprint changed")
    require(np.asarray(grid).size == 120 and int(P.Nz) == 9, "Candidate 120x9 grid changed")
    helpers_dir=ROOT / "output/model/fixed_reference_economics_20260928/elasticity_v1/source_v2"
    helper_spec=importlib.util.spec_from_file_location("mechanism_credit_report_helpers",helpers_dir/"run_credit.py")
    report_helpers=importlib.util.module_from_spec(helper_spec);helper_spec.loader.exec_module(report_helpers)
    return dict(inputs=inputs, native=native, credit=credit, solver=solver, ge=ge, report_helpers=report_helpers,
                base=base, context=ctx, P=P, natural=natural, grid=grid,
                point=point, entry=entry, params_rows=params_rows, actual_parameters=actual)


def make_price_cell(auth, regime, factor, out, deadline):
    native = auth["native"]
    P0 = auth["P"] if regime == "reference" else auth["natural"]
    P = copy.deepcopy(P0)
    grid = auth["grid"]
    q0 = float(BINDING["candidate_price"])
    q = q0 * float(factor)
    P.native_inherited_distribution_evidence_dir=str(out/"inherited_state_diagnostics")
    sd = auth["solver"].precompute_shared(P, grid)
    start = time.monotonic()
    auth["lifecycle_solves"]=auth.get("lifecycle_solves",0)+1
    require(auth["lifecycle_solves"]<=6,"Six actual lifecycle-attempt cap exceeded")
    write(out/"latest.json",dict(status="lifecycle_claimed",regime=regime,price=q,lifecycle_used=auth["lifecycle_solves"],deadline_epoch=deadline,pid=os.getpid()))
    if regime == "lifetime_repayment_only":
        import refactor_lab.engine.household as household
        from credit_mode import trace_native_support, lower_grid_diagnostic
        with trace_native_support(household) as support_calls:
            sol = auth["solver"].solve_markov_income_at_prices(np.asarray([q]), P, grid,
                SD=sd, verbose=False, fast_stats=False)
        support = None
    else:
        sol = auth["solver"].solve_markov_income_at_prices(np.asarray([q]), P, grid,
            SD=sd, verbose=False, fast_stats=False)
        support = {"status": "not_applicable_reference_credit"}
    require(float(getattr(P,"_entry_censored_mass",0.0)) <= auth["credit"].DEAD_MASS_TOL,
            "Inherited entry censoring would remove occupied mass")
    require(time.time() < deadline,"Case/global deadline exceeded during solve")
    elapsed = time.monotonic() - start
    require(elapsed <= float(CONTRACT["price_response_budget"]["maximum_seconds_per_case"]),
            "Per-case lifecycle/report deadline exceeded")

    cal = auth["context"]["prepared"].rt["primitive"].pf.calendar
    P._fert2_probs = sol.fert2_probs.copy()
    policy = cal.policy_from_solution(sol, np.asarray([q]), P, grid, sd)
    cohort_pre, reconstruction = cal.reconstruct_stationary_pre_fertility(sol, policy, P, grid, sd)
    auth["context"]["runtime"].require_abs_gate(reconstruction["stationary_post_fertility_nesting_l1"],5e-9,"Cohort reconstruction")
    auth["context"]["runtime"].require_abs_gate(reconstruction["stationary_feasibility_projection_mass"],0.,"Cohort projection")
    supply = cal.HousingSupplyRule("static-elastic", q,
        float(P.H0[0] * (P.user_cost_rate * q / P.r_bar[0]) ** P.xi_supply[0]), float(P.xi_supply[0]))
    ev = cal.evaluate_period(np.asarray([q]), cohort_pre, P, grid, sd,
        cal.SolveCounter(), supply_rule=supply, supplied_policy=policy)
    mass_rows = {name: float(getattr(ev, name).sum()) for name in ("g_pre", "g_post_fertility", "g_current")}
    require(max(abs(value - mass_rows["g_pre"]) for value in mass_rows.values()) <= 2e-10,
            "Fixed-price within-period household mass does not conserve")
    if regime=="lifetime_repayment_only":
        support=lower_grid_diagnostic(sol,grid,support_calls,P,realized_distribution=ev.g_current)
    packet = dict(parameters=P,b_grid=grid,shared=sd,solution=sol,evaluation=ev,
                  stationary_g_pre=cohort_pre,supply_rule=supply,demographic_seed=None)
    cohort_gates=case_gates(auth,packet,out,regime,stationary=True)
    fiscal_certificate=cohort_gates["fiscal_certificate"]

    fertility = {k: auth["context"]["prepared"].rt["observe_initial_fertility"](ev, P, age_projection=k)
                 for k in ("uniform_birth_time", "constant_post_cell")}
    housing = auth["context"]["prepared"].rt["observe_initial_housing_wealth"](ev, P, grid, sd,
        diagnostic_enabled=True, age_projection="uniform_within_age_cell",
        diagnostic_allow_family_proxies=True, include_wealth=True, include_birth_response=True)
    recent = auth["context"]["prepared"].rt["observe_recent_parent_flow"](ev, P,
        diagnostic_enabled=True, snapshot=auth["context"]["prepared"].rt["SNAPSHOT"],
        age_projection=auth["context"]["prepared"].rt["AGE_PROJECTION"],
        diagnostic_allow_residence_proxy=True,
        input_provenance={"case_id": f"{regime}_{factor:+.2%}", "candidate_case": CONTRACT["candidate_case"]})
    completed = float(auth["context"]["prepared"].rt["chain"].extract_moments(sol, P)["tfr"])
    runtime = auth["context"]["runtime"]
    fits = runtime.score_targets(auth["context"]["objective"], fertility, housing,
                                 recent["model_value"], completed)
    # The unchanged native comparator expects the calibration CSV representation.
    native.table(out / "target_fit.csv", fits)
    fits = native.readtable(out / "target_fit.csv")
    require(len(fits)==14 and native.target_identity(fits)==native.PLAN["target_contract"],"All 14 original targets/weights required")
    native.residual(fits)
    actual=auth["context"]["fp"].actual_parameters(auth["context"]["prepared"],P,grid)
    auth["ge"].validate_parameter_estimates(auth["context"],auth["params_rows"],actual)
    require(actual==auth["actual_parameters"],"Fixed-price case changed an effective candidate parameter")

    # The first cell is the sole source of q0 inherited states. Every impact
    # cell reuses that exact PRE-fertility distribution and its own policy.
    if regime == "reference" and float(factor) == 1.0:
        impact_pre = np.asarray(cohort_pre).copy()
        require(np.isfinite(impact_pre).all() and float(impact_pre.sum()) > 0.0,
                "q0 reference pre-fertility distribution is missing or invalid")
        np.savez_compressed(out.parent / "q0_reference_inherited_states.npz", g_pre=impact_pre)
        write(out.parent / "q0_reference_inherited_states.json", {
            "source_case": "reference_price_1.0", "shape": list(impact_pre.shape),
            "timing": "reconstructed_stationary_pre_fertility", "sum": float(impact_pre.sum()), "sha256": hashlib.sha256(impact_pre.tobytes()).hexdigest()})
    else:
        p = out.parent / "q0_reference_inherited_states.npz"
        require(p.is_file(), "q0 reference impact distribution not yet checkpointed")
        with np.load(p,allow_pickle=False) as saved: impact_pre=saved["g_pre"].copy()
    impact_P=copy.deepcopy(P)
    impact_sd = auth["solver"].precompute_shared(impact_P, grid)
    impact = cal.evaluate_period(np.asarray([q]), impact_pre, impact_P, grid, impact_sd,
        cal.SolveCounter(), supply_rule=supply, supplied_policy=policy)

    impact_out=out/"baseline_state_impact";impact_out.mkdir()
    impact_packet=dict(packet,parameters=impact_P,shared=impact_sd,evaluation=impact,stationary_g_pre=impact_pre)
    impact_gates=case_gates(auth,impact_packet,impact_out,regime,stationary=False)
    if regime=="reference" and float(factor)==1.0:
        for name in ("g_pre","g_post_fertility","g_current"):
            require(np.array_equal(getattr(impact,name),getattr(ev,name)),"q0 pre-fertility impact replay differs: "+name)
    # Existing stable 17-panel audit packet; market/renewal closure is recorded,
    # never imposed, at these prescribed prices.
    packet = dict(parameters=P, b_grid=grid, shared=sd, solution=sol,
                  evaluation=ev, stationary_g_pre=cohort_pre, supply_rule=supply,
                  demographic_seed=None)
    auth["context"]["prepared"].rt["audit"].standard_diagnostics(packet, out,
        validate_production_young=False)
    plots = sorted(x.name for x in (out / "standard_diagnostics").glob("*.png"))
    require(plots==sorted(auth["context"]["manifest"]["standard_diagnostic_names"]),"Standard 17 plot identities differ")

    params = copy.deepcopy(auth["params_rows"])
    for row in params:
        row["reference_estimate"] = row["estimate"]
        row["status"] = "Fixed current-candidate value during prescribed-price response"
        if row["parameter"] == "psi_child":
            row["status"] = "Fixed at current candidate psi; no refit or fertility normalization"
    native.table(out / "target_fit.csv", fits)
    native.table(out / "parameters.csv", params)
    native.write(out / "observers.json", cal.jsonable(dict(fertility=fertility,
        housing_wealth=housing, recent_parent=recent)))
    np.savez_compressed(out / "solution_arrays.npz", **{
        k: v for k, v in vars(sol).items() if isinstance(v, np.ndarray) and v.dtype != object})
    report_helpers=auth["report_helpers"]
    cohort_summary=report_helpers.aggregates(auth["context"]["fp"],ev,P,grid,auth["context"]["prepared"].rt["model"])
    cohort_summary["reference_label"]=CONTRACT["candidate_label"]
    cohort_summary.update(completed_fertility=completed,mean_first_birth_age=next(float(x["model"]) for x in fits if x["moment"]=="nchs_mean_age"))
    impact_summary=report_helpers.aggregates(auth["context"]["fp"],impact,impact_P,grid,auth["context"]["prepared"].rt["model"])
    impact_summary["reference_label"]=CONTRACT["candidate_label"]
    impact_summary["inherited_distribution_sha256"]=hashlib.sha256(impact_pre.tobytes()).hexdigest()
    write(impact_out/"summary.json",cal.jsonable(impact_summary))
    transition=auth["context"]["prepared"].rt["primitive"].pf.transition
    adjusted_births=float(transition.calendar_topcode_birth_accounting(ev.g_pre,ev.g_post_fertility,float(ev.births),P)["topcode_adjusted_birth_children"])
    entry_rate=float(sol.entry_rate)
    require(np.isfinite(adjusted_births) and np.isfinite(entry_rate) and entry_rate>0.,"Invalid prescribed-price renewal accounting")
    closure = {"status": "prescribed_price_no_market_or_renewal_root", "price": q,
        "adjusted_births":adjusted_births,"actual_entry_rate":entry_rate,
        "price_factor": float(factor), "mapped_rent": float(P.user_cost_rate * q),
        "renewal_residual_reported_not_imposed": float(adjusted_births / (2.1 * entry_rate) - 1.0),
        "relative_market_residual_reported_not_imposed": float(ev.relative_market_residual),
        "housing_demand": np.asarray(ev.demand_by_loc).tolist(), "housing_supply": np.asarray(ev.supply_by_loc).tolist(),
        "mass_conservation": mass_rows, "pension_paygo_certificate_reported_not_imposed": fiscal_certificate,
        "reconstruction": reconstruction, "support_diagnostic": support, "lifecycle_seconds": elapsed,
        "lifecycle_solves": 1, "production_adoption": False, "grid_nodes": int(len(grid)),
        "candidate_base_loss": float(BINDING["candidate_loss"]), "candidate_price_q0": q0,
        "candidate_psi_child": float(P.psi_child), "standard_plot_count": len(plots),
        "baseline_state_impact": impact_summary, "cohort_summary":cohort_summary,
        "cohort_gates":cohort_gates,"impact_gates":impact_gates,"entry_censored_mass":float(getattr(P,"_entry_censored_mass",0.)),
        "natural_support_certified":False,"market_clearing_certified":False,"demographic_renewal_certified":False,
        "economic_changes": [] if regime == "reference" else ["DUE stayer rule disabled", "fixed unsecured limit cleared", "native solvency credit enabled"],
        "support_limitation": "Occupied support guard does not certify unoccupied continuation alternatives or exact full natural support."}
    write(out / "closure.json", closure)
    write(out / "receipt.json", {"status": "completed_experimental_case" if regime=="reference" else "completed_support_limited_diagnostic",
        "natural_support_certified":False, "regime": regime,
        "price_factor": float(factor), "elapsed_seconds": elapsed,
        "driver_sha256": sha(__file__), "plot_names": plots,
        "target_fit_sha256": sha(out / "target_fit.csv"), "parameters_sha256": sha(out / "parameters.csv"),
        "support_status": support["status"], "source_binding_sha256": sha(HERE / "source_binding.json")})
    require(time.time()<deadline,"Case/global deadline exceeded during reporting")
    return closure


def case_gates(auth,packet,out,regime,*,stationary):
    """Reuse universal native gates; natural credit replaces only the artificial LTV ledger."""
    prepared=auth["context"]["prepared"];fp=auth["context"]["fp"]
    if regime=="reference":return fp.gates(packet,prepared,out,stationary=stationary)
    sys.path.insert(0,str(ROOT/"code/model/tools"))
    from e5f_solvency_credit_benchmark import audit_purchase_accounting
    original=prepared.rt["accounting"]
    prepared.rt["accounting"]=SimpleNamespace(audit_purchase_accounting=audit_purchase_accounting)
    try:result=fp.gates(packet,prepared,out,stationary=stationary)
    finally:prepared.rt["accounting"]=original
    result.update(natural_support_certified=False,credit_audit_scope="exact transactions and possible-death net estates; universal household/forward gates retained; full natural support unverified")
    write(Path(out)/"gates.json",prepared.rt["primitive"].pf.calendar.jsonable(result))
    return result


def write_comparison(output,records):
    """The existing elasticity experiment's central/one-sided log differences at ±1%."""
    import math
    rows=[];outcomes=("births_per_household","first_births","second_births","third_bin_entries","completed_fertility","mean_first_birth_age")
    for regime in ("reference","lifetime_repayment_only"):
        completed=[r for r in records if r.get("regime")==regime and r.get("status") in ("completed_experimental_case","completed_support_limited_diagnostic")]
        if len(completed)!=3:
            write(output/(regime+"_comparison_unresolved.json"),dict(status="incomplete_three_price_cells",cases=completed));continue
        receipts={float(r["price_factor"]):json.loads(Path(r["closure_path"]).read_text()) for r in completed}
        support_statuses=[x["support_diagnostic"]["status"] for x in receipts.values()]
        occupied_pass=all(x=="occupied_support_pass_unoccupied_alternatives_unverified" for x in support_statuses)
        if regime!="reference" and not occupied_pass:
            write(output/(regime+"_comparison_unresolved.json"),dict(status="support_uncertified",support_statuses=support_statuses,natural_support_certified=False));continue
        for scope,key in (("impact","baseline_state_impact"),("cohort","cohort_summary")):
            for outcome in outcomes if scope=="cohort" else outcomes[:4]:
                if not all(outcome in cell[key] for cell in receipts.values()):continue
                values={factor:float(cell[key][outcome])/float(cell[key]["household_mass"]) if outcome in ("first_births","second_births","third_bin_entries") else float(cell[key][outcome]) for factor,cell in receipts.items()}
                def ratio(a,b,fa,fb):return (math.log(b)-math.log(a))/(math.log(fb)-math.log(fa)) if a>0 and b>0 else None
                rows.append(dict(regime=regime,scope=scope,outcome=outcome,price_099=values[.99],price_100=values[1.],price_101=values[1.01],central_log_elasticity=ratio(values[.99],values[1.01],.99,1.01),lower_log_elasticity=ratio(values[.99],values[1.],.99,1.),upper_log_elasticity=ratio(values[1.],values[1.01],1.,1.01),status="prescribed_price_response" if regime=="reference" else "support_limited_diagnostic",natural_support_certified=False,market_clearing_certified=False))
    if rows:table(output/"elasticities.csv",rows)


def run(output, deadline_epoch, *, auth=None, evaluate=None, mock=False):
    started=time.time();deadline_epoch=min(float(deadline_epoch),started+1200.)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    try:auth = authenticate_candidate(output/"runtime_preparation") if auth is None else auth
    except BaseException as exc:
        write(output/"failure.json",dict(status="fatal_source_or_initializer_failure",error_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc(),lifecycle_solves=0,no_auto_retry=True))
        raise
    evaluate = make_price_cell if evaluate is None else evaluate
    attempted=0
    write(output / "launch.json", {"status": "running", "contract_sha256": sha(HERE / "contract.json"),
        "source_binding_sha256": sha(HERE / "source_binding.json"), "driver_sha256": sha(__file__),
        "started_epoch": started, "deadline_epoch": deadline_epoch, "lifecycle_solves": 0,
        "candidate": CONTRACT["candidate_case"], "case_budget_seconds": 300,
        "global_budget_seconds": 1200, "lifecycle_solve_cap": 6,
        "cpu_threads": 1, "memory_gib_cap": 24})
    cases = [("reference", 1.0), ("reference", .99), ("reference", 1.01),
             ("lifetime_repayment_only", .99), ("lifetime_repayment_only", 1.0),
             ("lifetime_repayment_only", 1.01)]
    records = []
    write(output / "latest_completed.json", {"completed": records, "lifecycle_solves": 0})
    for index, (regime, factor) in enumerate(cases):
        if time.time() >= deadline_epoch:
            records.append({"regime": regime, "price_factor": factor, "status": "not_started_global_deadline"})
            break
        label = f"{index:02d}_{regime}_p{factor:.2f}"
        out = output / label
        out.mkdir(parents=True, exist_ok=False)
        case_deadline = min(deadline_epoch, time.time() + 300.0)
        write(out / "latest.json", {"status": "running", "regime": regime,
            "price_factor": factor, "started_epoch": time.time(), "case_deadline_epoch": case_deadline})
        old_handler = signal.getsignal(signal.SIGALRM)
        signal.signal(signal.SIGALRM, _case_alarm)
        old_timer = signal.setitimer(signal.ITIMER_REAL, max(.001, case_deadline - time.time()))
        try:
            attempted+=1
            require(attempted<=6,"Six lifecycle-attempt cap exceeded")
            result = evaluate(auth, regime, factor, out, case_deadline)
            row = {"label": label, "regime": regime, "price_factor": factor,
                "status": "completed_experimental_case" if regime=="reference" else "completed_support_limited_diagnostic", "lifecycle_seconds": result["lifecycle_seconds"],
                "support_status": result["support_diagnostic"]["status"], "closure_path": str(out / "closure.json")}
        except Exception as exc:
            row = {"label": label, "regime": regime, "price_factor": factor,
                "status": "failed_or_unresolved", "error": repr(exc), "traceback": traceback.format_exc()}
            row["failure_class"]="declared_case_timeout" if isinstance(exc,CaseDeadline) else "fatal_accounting_source_target_shape_or_unexpected"
            write(out / "failure.json", row)
            if index==0 or not isinstance(exc,CaseDeadline):
                records.append(row)
                failure=dict(status="fatal_q0_reference_failure" if index==0 else "fatal_case_failure",failed_case=row,cases=records,lifecycle_solves=0 if mock else auth.get("lifecycle_solves",0),case_attempts=attempted,no_auto_retry=True,dependent_cases_blocked=True)
                write(output/"failure.json",failure)
                write(output/"latest_completed.json",failure)
                write(output/"completed.json",failure)
                raise
        finally:
            signal.setitimer(signal.ITIMER_REAL, *old_timer)
            signal.signal(signal.SIGALRM, old_handler)
        records.append(row)
        write(output / "latest_completed.json", {"completed": records,
            "lifecycle_solves":0 if mock else auth.get("lifecycle_solves",0),"case_attempts":attempted,
            "elapsed_seconds": time.time() - float(json.loads((output / "launch.json").read_text())["started_epoch"])})
        write(output / "best_so_far.json", {"status": "mechanism_responses_not_adopted",
            "completed": records, "candidate_loss": float(BINDING["candidate_loss"])})
    if not mock:write_comparison(output,records)
    write(output / "completed.json", {"status": "finished_with_case_statuses", "cases": records,
        "lifecycle_solves":0 if mock else auth.get("lifecycle_solves",0),"case_attempts":attempted,
        "production_adoption": False})
    return records


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--deadline-epoch", type=float)
    parser.add_argument("--smoke-loop", action="store_true")
    args = parser.parse_args()
    if args.smoke_loop:
        def fake(auth,regime,factor,out,end):
            write(out/"mock_case.json",dict(regime=regime,factor=factor,lifecycle_solves=0))
            return dict(lifecycle_seconds=0.,support_diagnostic=dict(status="mock_no_support_certificate"))
        records=run(args.out,time.time()+1200,auth={},evaluate=fake,mock=True)
        require(len(records)==6,"Exact mock loop did not visit all six cells")
        print(json.dumps(dict(status="mock_loop_passed",lifecycle_solves=0,cases=records)))
    else:
        if args.deadline_epoch is None:
            parser.error("--deadline-epoch is required for a numerical dispatch")
        run(args.out, args.deadline_epoch)
