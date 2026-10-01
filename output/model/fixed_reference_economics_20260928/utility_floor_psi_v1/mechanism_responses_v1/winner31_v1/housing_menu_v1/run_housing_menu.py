#!/usr/bin/env python3
"""Two authorized fixed-q Stone–Geary owned-room menu experiments; no GE/refit."""
import copy,json,os,signal,sys,time,traceback
from pathlib import Path
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'): os.environ[k]='1'
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import fixed_price_responses as d
from fixed_price_responses import (np,require,write,case_gates,CONTRACT,BINDING,sha)
BASELINE=HERE.parent/'purchase_ltv_v1/local_run/retry5/results/baseline_80_80'
MENUS=[('add3',[2.,3.,4.,6.,8.,10.]),('add1_3_5_7',[1.,2.,3.,4.,5.,6.,7.,8.,10.])]
def alarm(*_): raise TimeoutError('120-second cell or 300-second global deadline')
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
    require(auth["lifecycle_solves"]<=2,"Two actual lifecycle-attempt cap exceeded")
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

    impact_summary = None
    impact_gates = None
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
        "economic_changes": ["Experimental owned-room menu expanded; fixed tenure logit scale implies added variety value"],
        "owned_room_menu": P.H_own.tolist(), "n_owned_options": P.n_house, "n_tenure_options": 1+P.n_house,
        "impact_scope": "stationary cohort only; no fixed-PRE state embedding or replay",
        "support_limitation": "Occupied support guard does not certify unoccupied continuation alternatives or exact full natural support."}
    from refactor_lab.engine.parameters import children_at_home_count
    current=np.asarray(ev.g_current).sum(axis=(0,2,4))  # tenure, age, ever-born, child-state
    parents=np.asarray([[children_at_home_count(n,c,P)>0 for c in range(P.n_child_states)] for n in range(P.n_parity)])
    ages=float(P.age_start)+np.arange(P.J)*float(P.da)
    young=(ages>=18)&(ages<=42)
    menu_rows=[]
    for t,H in enumerate([None]+P.H_own.tolist()):
        mass=float(current[t].sum()); pmass=float((current[t]*parents[None,:,:]).sum())
        youngmass=float(current[t,young].sum()); yparent=float((current[t,young]*parents[None,:,:]).sum())
        menu_rows.append(dict(tenure_index=t,owned_rooms=H,current_household_mass=mass,current_parent_mass=pmass,age18_42_mass=youngmass,age18_42_parent_mass=yparent))
    write(out/'menu_use.json',dict(state_age_clock='model age-cell start; supplemental diagnostic, no target clock change',rows=menu_rows,household_mass=float(current.sum()),current_parent_mass=float((current*parents[None,None,:,:]).sum()),age18_42_mass=float(current[:,young].sum()),age18_42_parent_mass=float((current[:,young]*parents[None,None,:,:]).sum()),first_birth_conditional_menu_mass='not extracted; stationary realized tenure rather than event-conditioned causal response'))
    closure['menu_use_path']=str(out/'menu_use.json')
    write(out / "closure.json", closure)
    write(out / "receipt.json", {"status": "completed_experimental_case" if regime=="reference" else "completed_support_limited_diagnostic",
        "natural_support_certified":False, "regime": regime,
        "price_factor": float(factor), "elapsed_seconds": elapsed,
        "driver_sha256": sha(__file__), "plot_names": plots,
        "target_fit_sha256": sha(out / "target_fit.csv"), "parameters_sha256": sha(out / "parameters.csv"),
        "support_status": support["status"], "source_binding_sha256": sha(d.HERE / "source_binding.json")})
    require(time.time()<deadline,"Case/global deadline exceeded during reporting")
    return closure

def main(output,initialize_only=False):
 started=time.time();deadline=started+300.;output=Path(output);output.mkdir(parents=True,exist_ok=False)
 auth=d.authenticate_candidate(output/'runtime_preparation')
 P=auth['P'];require(np.array_equal(P.H_own,[2,4,6,8,10]) and P.n_house==5,'Reference menu drift')
 require(np.all(np.asarray(P.phi)==.8) and P.native_due_stayer_credit and not getattr(P,"native_solvency_credit",False) and P.native_purchase_income and not P.use_pti_constraint and not P.joint_nested_choice,'Reference credit drift')
 require(abs(float(BINDING['candidate_price'])-.719168368828958)<1e-15,'Reference q drift')
 require(BASELINE.joinpath('target_fit.csv').is_file(),'Saved baseline unavailable')
 original=copy.deepcopy(P);checks=[]
 for label,menu in MENUS:
  trial=copy.deepcopy(original);trial.H_own=np.asarray(menu);trial.n_house=len(menu)
  sd=auth['solver'].precompute_shared(trial,auth['grid'])
  require(sd.birth_dp.shape[-2:]==(1+len(menu),1+len(menu)),'Birth choice axes did not update')
  require(auth['context']['fp'].actual_parameters(auth['context']['prepared'],trial,auth['grid'])==auth['actual_parameters'],'31 native parameters changed')
  checks.append(dict(label=label,H_own=menu,n_house=trial.n_house,nt=1+trial.n_house,birth_dp_shape=list(sd.birth_dp.shape),phi_choice_shape=list(sd.phi_choice.shape) if hasattr(sd,'phi_choice') else None))
 write(output/'initializer.json',dict(status='zero_lifecycle_initializer_passed',cells=checks,baseline_path=str(BASELINE),baseline_target_sha256=sha(BASELINE/'target_fit.csv'),baseline_parameters_sha256=sha(BASELINE/'parameters.csv'),driver_sha256=sha(__file__),all31_parameters_unchanged=True,price=float(BINDING['candidate_price']),no_fixed_PRE_replay=True,deadline_epoch=deadline))
 if initialize_only:return
 records=[]
 for label,menu in MENUS:
  end=min(deadline,time.time()+120);require(end>time.time(),'Global deadline')
  auth['P']=copy.deepcopy(original);auth['P'].H_own=np.asarray(menu);auth['P'].n_house=len(menu)
  case=output/label;case.mkdir();write(output/'latest_completed.json',dict(completed=records,running=label,deadline_epoch=deadline))
  old=signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,end-time.time())
  try:closure=make_price_cell(auth,'reference',1.,case,end)
  except BaseException as e:
   write(output/'failure.json',dict(error=str(e),traceback=traceback.format_exc(),completed=records,lifecycle_attempts=auth.get('lifecycle_solves',0)));raise
  finally:signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)
  records.append(dict(label=label,H_own=menu,closure_path=str(case/'closure.json'),lifecycle_seconds=closure['lifecycle_seconds'],cohort_summary=closure['cohort_summary']))
  write(output/'latest_completed.json',dict(completed=records,lifecycle_attempts=auth['lifecycle_solves']))
 write(output/'completed.json',dict(status='two_completed_fixed_price_menu_experiments',completed=records,lifecycle_attempts=auth['lifecycle_solves'],total_seconds=time.time()-started,production_adoption=False))
if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--initialize-only',action='store_true');a=ap.parse_args();main(a.out,a.initialize_only)
