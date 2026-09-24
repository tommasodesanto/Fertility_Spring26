#!/usr/bin/env python3
"""Experimental stationary calibration from the verified B-floor checkpoint.

The native model, selected checkpoint, objective ancestry, and PAYGO solver are
pinned. No demographic birth-entry queue is implemented in this solve.
"""
from __future__ import annotations
import argparse, copy, csv, fcntl, gzip, hashlib, importlib.util, json, math, os, pickle, random, signal, sys, time
from pathlib import Path
from types import SimpleNamespace
import numpy as np

BASE=Path("/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a")
BUNDLE=BASE/"utility_overnight_20260923_v1"
WORK=BASE/"commute_calibration_20260924_v1"
TAX_DRIVER=BASE/"paygo_tax_comparison_20260924/run_paygo_two_rate.py"
PLAN=BUNDLE/"results/production/B_floor/worker09_proposal16/plan.json"
CHECKPOINT=BUNDLE/"results/production/B_floor/worker09_proposal16/result/evaluation/raw/repetition_01/initial_state.pkl.gz"
SOURCE=BUNDLE/"source"
OBJECTIVE=WORK/"inputs/objective.json"
OBJECTIVE_SHA="c28e5d620d463dc3a592c038ca2771b47bf37ed1a9cf5ba31f7dc793b58b61db"
TAX=0.08751017424959717
THETA1=0.008193084126995582
ANNUAL_DEP=0.01416143718381309
ANNUAL_PROPERTY_TAX=0.010598360773872594
START_PSI=0.14245465246024056
FREE=("H0","beta_annual","chi","first_birth_fixed_cost","h_P","kappa_fert","kappa_fert_continuation","theta0")


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""): h.update(block)
    return h.hexdigest()


def write(path, value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+".tmp")
    tmp.write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+"\n")
    tmp.replace(path)


def table(path, rows):
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open("w",newline="") as f:
        writer=csv.DictWriter(f,fields); writer.writeheader(); writer.writerows(rows)


def load_tax_driver():
    spec=importlib.util.spec_from_file_location("verified_tax_driver", TAX_DRIVER)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def prepare_selected():
    tax=load_tax_driver()
    plan=tax.source_preflight(SOURCE,PLAN)
    assert sha(CHECKPOINT)==tax.CHECKPOINT_SHA
    assert sha(OBJECTIVE)==OBJECTIVE_SHA
    sys.path[:0]=[str(SOURCE/"code/model/tools"),str(SOURCE/"code/model")]
    sys.dont_write_bytecode=True
    with gzip.open(CHECKPOINT,"rb") as f: selected=pickle.load(f)
    tax.check_checkpoint(selected)
    objective=json.loads(OBJECTIVE.read_text())
    assert len(objective["parameter_restrictions"])==8
    assert tuple(sorted(r["parameter"] for r in objective["parameter_restrictions"]))==tuple(sorted(FREE))
    assert len(objective["target_rows"])==13
    assert sum(r["actual_weight"] is not None for r in objective["target_rows"])==12
    assert objective["external_restrictions"]["theta1"]==THETA1
    return tax,plan,selected,objective


def apply_point(selected, point):
    """Bind eight free coordinates with reviewed parenthood utility mapping."""
    assert set(point)==set(FREE)
    from e5f_parenthood_utility import bind_parenthood_utility,validate_parenthood_utility
    P=bind_parenthood_utility(copy.deepcopy(selected["parameters"]),point)
    P.theta1=THETA1  # external B15 wealth/earnings restriction; below old searched bound
    P.delta=1.0-(1.0-ANNUAL_DEP)**int(P.period_years)
    P.tau_H=ANNUAL_PROPERTY_TAX*int(P.period_years)
    P.user_cost_rate=P.q+P.delta+P.tau_H
    validate_parenthood_utility(P)
    assert P.theta1==THETA1 and P.sigma==2.0 and P.psi==0.06
    assert np.all(np.asarray(P.phi)==0.8) and P.xi_supply[0]==0.63
    assert P.user_cost_rate==P.q+P.delta+P.tau_H
    return P


def normalize_point(*, selected, model, chain, native_solver, calibration, point, deadline_epoch, ledger_path):
    """Call the unchanged reviewed 2.1 normalizer through a scoped GE bridge."""
    base=apply_point(selected,point)
    original=calibration.closure
    solve_ledger=[]
    def solve_ge(_chain, overrides):
        if time.time()>=deadline_epoch: raise TimeoutError("shared deadline before native solve")
        P=copy.deepcopy(base)
        P.psi_child=float(overrides["psi_child"])
        record=dict(index=len(solve_ledger),psi_child=float(P.psi_child),status="started",started_epoch=time.time())
        solve_ledger.append(record);write(ledger_path,solve_ledger)
        start=time.monotonic()
        try:
            sol,P,price,fiscal=native_solver(model=model,parameters=P,
                b_grid=selected["b_grid"],initial_prices=selected["solution"].p_eq,
                payroll_tax=TAX,marginal_tolerance=1e-9,fiscal_tolerance=1e-6)
        except Exception as exc:
            record.update(status="failed",error_type=type(exc).__name__,error=str(exc),seconds=time.monotonic()-start)
            write(ledger_path,solve_ledger);raise
        record.update(status="completed",seconds=time.monotonic()-start,
                      price=float(price[0]),market_residual=float(sol.timings["best_eq_error"]))
        write(ledger_path,solve_ledger)
        if time.time()>=deadline_epoch: raise TimeoutError("shared deadline after native solve")
        return sol,P,price,time.monotonic()-start
    calibration.closure=SimpleNamespace(solve_ge=solve_ge)
    try:
        return calibration.solve_old_steady_state(chain,{},initial_psi=START_PSI,
            completed_fertility_target=2.1,completed_fertility_tolerance=5e-4,normalize=True)
    finally:
        calibration.closure=original


def setup_runtime(tax, plan, selected, output):
    """Install the same pinned native runtime as the completed PAYGO comparison."""
    bundle_tools=SOURCE.parent/"tools"
    assert sha(bundle_tools/"e5f_earnings_wealth_contract.py")==plan["files"]["accounting"]["sha256"]
    assert sha(bundle_tools/"run_e5f_preference_share_candidate.py")==plan["files"]["adapter"]["sha256"]
    sys.path[:0]=[str(SOURCE/"code/model/tools"),str(SOURCE/"code/model")]
    import run_e5f_matched_pf_smoke as primitive
    import run_e5f_transition_calibration as calibration
    sys.path.insert(0,str(bundle_tools))
    import e5f_earnings_wealth_contract as accounting
    import run_e5f_independent_numerical_audit as audit
    from e5f_stationary_paygo import solve_balanced_initial_equilibrium,certify_initial_pension
    from e5f_social_security import fiscal_accounts
    from e5f_initial_fertility_observer import observe_initial_fertility
    from e5f_initial_housing_observer import observe_initial_housing_wealth
    from e5f_recent_parent_flow_observer import observe_recent_parent_flow,SNAPSHOT,AGE_PROJECTION
    chain,model=primitive.pf.transition.configure_sequential_model()
    accounting.install_fixed_entry(model)
    accounting.install_explicit_transaction_grid(model)
    runtime_diff=output/"reviewed_purchase_runtime.diff"
    accounting.install_purchase_income(model,runtime_diff)
    expected={
        ".diff":"592a3c90bbc187ba6aecd5fe0b6491ef2f90aa965f4043eae0101baceac23e4f",
        ".generated.py":"900b8bac964a2f91e8f9f47c566cd89c25e316e5db0afbe0975b58cf1a52bd28",
        ".tenure.py":"01ff6c5402cf8317d4b927916f0ac1233026fece2716138163a44990b13403a0",
        ".allocation.py":"86c56d2666f69f757dd69a5b87647686bbd5c4205a8a07cf20651317065e4d8a",
    }
    for suffix,pin in expected.items(): assert sha(runtime_diff.with_suffix(suffix))==pin
    primitive.pf.calendar.apply_fertility=primitive.pf.transition.apply_sequential_fertility
    primitive.pf.calendar.advance_calendar_distribution=primitive.pf.transition.advance_sequential_calendar_distribution
    primitive.pf.transition.calendar.model=model
    return locals()


def evaluate_point(*, tax, objective, selected, runtime, point, output, deadline_epoch, graphs=False):
    bounds={r["parameter"]:(float(r["lower"]),float(r["upper"]))
            for r in objective["parameter_restrictions"]}
    assert set(bounds)==set(point)==set(FREE)
    assert all(math.isfinite(float(point[name])) and low<=float(point[name])<=high
               for name,(low,high) in bounds.items())
    output.mkdir(parents=True,exist_ok=False)
    chain=runtime["chain"]; model=runtime["model"]; primitive=runtime["primitive"]
    calibration=runtime["calibration"]
    sol,P,price,chosen_solve_seconds,normalization=normalize_point(selected=selected,model=model,
        chain=chain,native_solver=runtime["solve_balanced_initial_equilibrium"],
        calibration=calibration,point=point,deadline_epoch=deadline_epoch,
        ledger_path=output/"stationary_solves.json")
    assert abs(normalization["completed_fertility"]-2.1)<=5e-4
    assert P.theta1==THETA1 and P.tau_pay==TAX
    grid=selected["b_grid"]
    shared=model.precompute_shared(P,grid)
    P._fert2_probs=sol.fert2_probs.copy()
    policy=primitive.pf.calendar.policy_from_solution(sol,price,P,grid,shared)
    pre,reconstruction=primitive.pf.calendar.reconstruct_stationary_pre_fertility(sol,policy,P,grid,shared)
    operator=primitive.pf.transition.operator_gates(sol,policy,pre,P,grid,shared)
    operator.update(reconstruction)
    for name in ("stationary_post_fertility_nesting_l1","one_step_constant_path_nesting_l1",
                 "mature_flow_abs_error","birth_flow_abs_error","topcode_adjusted_birth_flow_abs_error"):
        assert abs(operator[name])<=5e-9,name
    assert abs(operator["zero_entry_mass_accounting_residual"])<=2e-8
    assert operator["stationary_feasibility_projection_mass"]<=1e-6
    supply=primitive.pf.calendar.HousingSupplyRule("static-elastic",float(price[0]),
        float(P.H0[0]*(P.user_cost_rate*price[0]/P.r_bar[0])**P.xi_supply[0]),float(P.xi_supply[0]))
    evaluation=primitive.pf.calendar.evaluate_period(price,pre,P,grid,shared,
        primitive.pf.calendar.SolveCounter(),supply_rule=supply,supplied_policy=policy)
    assert evaluation.relative_market_residual<=2e-4
    budget=primitive.dated_budget(evaluation,P,shared,grid,float(P.user_cost_rate*price[0]))
    purchase=runtime["accounting"].audit_purchase_accounting(evaluation,P,shared,grid,model)
    fiscal=runtime["certify_initial_pension"](evaluation.g_current,P,
        marginal_tolerance=1e-9,fiscal_tolerance=1e-6)
    packet=dict(parameters=P,b_grid=grid,evaluation=evaluation,shared=shared,supply_rule=supply,
        solution=sol,stationary_g_pre=pre,demographic_seed=selected.get("demographic_seed"),
        contract_sha256=OBJECTIVE_SHA,ancestry_contract_sha256=selected.get("contract_sha256"))
    arrays=runtime["audit"].policy_array_audit(packet,output)
    assert arrays["occupied_negative_steps"]==0
    assert all(not x["nonfinite"] and x["minimum"]>=0 and x["maximum"]<=1
               for x in arrays["probabilities"].values())
    fertility={projection:runtime["observe_initial_fertility"](evaluation,P,age_projection=projection)
               for projection in ("uniform_birth_time","constant_post_cell")}
    housing=runtime["observe_initial_housing_wealth"](evaluation,P,grid,shared,
        diagnostic_enabled=True,age_projection="uniform_within_age_cell",
        diagnostic_allow_family_proxies=True,include_wealth=True,include_birth_response=True)
    early=dict(fertility=fertility,housing_wealth=housing)
    checkpoint=output/"initial_state.pkl.gz"
    with gzip.open(checkpoint,"wb",compresslevel=1) as f: pickle.dump(packet,f,protocol=5)
    checkpoint_sha=sha(checkpoint)
    recent=runtime["observe_recent_parent_flow"](evaluation,P,diagnostic_enabled=True,
        snapshot=runtime["SNAPSHOT"],age_projection=runtime["AGE_PROJECTION"],
        diagnostic_allow_residence_proxy=True,
        input_provenance=dict(case_id=output.name,checkpoint_sha256=checkpoint_sha))
    moments=chain.extract_moments(sol,P)
    rows=tax.target_rows(objective,early,recent["model_value"],float(moments["tfr"]))
    for row in rows: row.pop("descriptive_fixed_psi",None)
    assert len(rows)==13 and sum(r["weight"]!="" for r in rows)==12
    for row in rows:
        assert all(math.isfinite(float(row[key])) for key in ("target","model","gap"))
        if row["weight"]!="":
            assert math.isfinite(float(row["weight"])) and float(row["weight"])>0
            assert math.isfinite(float(row["loss_contribution"])) and float(row["loss_contribution"])>=0
        else:
            assert row["moment"]=="initial_normalization" and row["loss_contribution"]==""
    loss=sum(float(r["loss_contribution"]) for r in rows if r["loss_contribution"]!="")
    assert math.isfinite(loss) and loss>=0
    params=[dict(parameter=name,estimate=tax.actual_parameters(P)[name],
        lower=next(r["lower"] for r in objective["parameter_restrictions"] if r["parameter"]==name),
        upper=next(r["upper"] for r in objective["parameter_restrictions"] if r["parameter"]==name),
        status="experimental free coordinate") for name in FREE]
    for row in params:
        low=float(row["lower"]);high=float(row["upper"]);v=float(row["estimate"])
        row["near_bound"]=min(v-low,high-v)<=.01*(high-low)
    params += [dict(parameter=name,estimate=value,lower="",upper="",near_bound="",status=status)
        for name,value,status in (
          ("theta1",P.theta1,"externally fixed B15"),
          ("psi_child",P.psi_child,"normalized to 2.1"),
          ("payroll_tax",P.tau_pay,"experimental PAYGO rate"),
          ("pension_period",P.pension,"endogenous balanced PAYGO"),
          ("housing_supply_elasticity",P.xi_supply[0],"retained external setting"),
          ("tenure_choice_kappa",P.tenure_choice_kappa,"retained external setting"),
          ("alpha_cons",P.alpha_cons,"retained external setting"),
          ("sigma",P.sigma,"retained external setting"),
          ("selling_cost",P.psi,"retained external setting"),
          ("financed_share",P.phi[0],"retained external setting"),
          ("annual_depreciation",ANNUAL_DEP,"author adopted input"),
          ("period_depreciation",P.delta,"four-year compounded"),
          ("annual_property_tax",ANNUAL_PROPERTY_TAX,"author adopted input"),
          ("period_property_tax",P.tau_H,"four-year linear source convention"),
          ("income_process",15,"retained B15 persistent-state count"),
          ("entrant_conversion_factor",P.entrant_conversion_factor,"retained stationary law"))]
    table(output/"target_fit.csv",rows);table(output/"parameters.csv",params)
    receipt=dict(status="verified_experimental_point",loss=loss,point=point,normalization=normalization,
        target_system_sha256=OBJECTIVE_SHA,source_manifest_sha256=tax.SOURCE_SHA,
        selected_checkpoint_sha256=tax.CHECKPOINT_SHA,case_checkpoint_sha256=checkpoint_sha,
        free_count=8,weighted_count=12,display_count=13,theta1_external=THETA1,
        payroll_tax=TAX,annual_depreciation=ANNUAL_DEP,annual_property_tax=ANNUAL_PROPERTY_TAX,
        period_depreciation=P.delta,period_property_tax=P.tau_H,user_cost_rate=P.user_cost_rate,
        price=float(price[0]),market_residual=evaluation.relative_market_residual,
        fiscal=fiscal,operator_gates=operator,household_budget=budget,purchase_accounting=purchase,
        policy_array_gates=arrays,chosen_solve_seconds=chosen_solve_seconds,
        objective_stationary_solve_seconds=normalization["stationary_solve_seconds"],
        objective_stationary_solves=normalization["stationary_solves"],
        model_observer_warning="First-birth housing is an unmatched stationary proxy for the PSID panel +3/+4 versus -3/-2 estimate",
        demographic_queue="Adopted 16/20 entry queue not implemented; existing stationary law retained")
    write(output/"receipt.json",tax.finite_json(primitive.pf.calendar.jsonable(receipt)))
    if graphs:
        runtime["audit"].standard_diagnostics(packet,output,validate_production_young=False)
        assert len(list((output/"standard_diagnostics").glob("*.png")))==17
    return receipt


def core_preflight():
    tax,plan,selected,objective=prepare_selected()
    sys.path[:0]=[str(SOURCE/"code/model/tools"),str(SOURCE/"code/model")]
    point={k:tax.PARAMETERS[k] for k in FREE}
    P=apply_point(selected,point)
    bounds={r["parameter"]:(float(r["lower"]),float(r["upper"]))
            for r in objective["parameter_restrictions"]}
    widths=plan["search_config"]["coordinate_widths"]
    assert set(bounds)==set(FREE) and set(FREE).issubset(widths)
    for worker_id in range(1,9):
        candidate,_=proposal(point,bounds,widths,random.Random(20260924+worker_id),1,worker_id)
        assert all(bounds[k][0]<=candidate[k]<=bounds[k][1] for k in FREE)
    return {"status":"core_preflight_passed_no_solve","source_sha256":tax.SOURCE_SHA,
        "checkpoint_sha256":tax.CHECKPOINT_SHA,"objective_sha256":OBJECTIVE_SHA,
        "free_coordinates":FREE,"theta1_external":P.theta1,"psi_start":START_PSI,
        "annual_depreciation":ANNUAL_DEP,"period_depreciation":P.delta,
        "annual_property_tax":ANNUAL_PROPERTY_TAX,"period_property_tax":P.tau_H,
        "payroll_tax":TAX,"user_cost_rate":P.user_cost_rate,"weighted_moments":12,
        "display_rows":13,"normalizer_source":"run_e5f_transition_calibration.py:765-868",
        "demographic_queue":"adopted but unimplemented; existing stationary law retained"}


def record_global(run_root, receipt, case_path):
    lock=run_root/"summary.lock"
    with lock.open("a") as f:
        fcntl.flock(f,fcntl.LOCK_EX)
        item=dict(status=receipt["status"],loss=receipt["loss"],point=receipt["point"],
                  case_path=str(case_path),updated_epoch=time.time(),
                  target_system_sha256=OBJECTIVE_SHA)
        write(run_root/"latest_completed.json",item)
        best_path=run_root/"best_so_far.json"
        if not best_path.exists() or item["loss"]<json.loads(best_path.read_text())["loss"]:
            write(best_path,item)
        fcntl.flock(f,fcntl.LOCK_UN)


def proposal(center, bounds, widths, rng, index, worker_id):
    name=FREE[(index+worker_id-1)%len(FREE)]
    lower,upper=bounds[name]
    scale=(0.8 if worker_id<=3 else 1.6 if worker_id<=6 else 3.2)
    for _ in range(100):
        value=center[name]+rng.gauss(0,widths[name]*scale)
        if lower<=value<=upper:
            candidate=dict(center);candidate[name]=value
            return candidate,dict(changed_coordinate=name,scale=scale)
    raise RuntimeError("bounded proposal generation failed: "+name)


def run_stage(stage, worker_id, run_root):
    run_root.mkdir(parents=True,exist_ok=True)
    if stage=="smoke":
        deadline_path=run_root/"deadline.json"
        assert not deadline_path.exists(),"smoke already started"
        deadline_epoch=time.time()+3600
        write(deadline_path,dict(first_smoke_start_epoch=time.time(),deadline_epoch=deadline_epoch,
            total_seconds=3600,export_reserve_seconds=240,objective_sha256=OBJECTIVE_SHA))
    else:
        data=json.loads((run_root/"deadline.json").read_text())
        assert data["objective_sha256"]==OBJECTIVE_SHA
        deadline_epoch=float(data["deadline_epoch"])
        assert (run_root/"best_so_far.json").exists(),"exact-loop smoke did not complete"
    # A hard process alarm ends native calls before the shared deadline and
    # reserves four minutes for selected-case export. Slurm/timeout also cap jobs.
    remaining=deadline_epoch-time.time()-240
    if remaining<=0: raise TimeoutError("no shared budget remains before export reserve")
    def alarm(_signum,_frame): raise TimeoutError("shared calibration deadline reached")
    signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,remaining)
    try:
        tax,plan,selected,objective=prepare_selected()
        sys.dont_write_bytecode=True
        case_root=run_root/("smoke" if stage=="smoke" else f"worker_{worker_id:02d}")
        case_root.mkdir(parents=True,exist_ok=False)
        os.environ["NUMBA_CACHE_DIR"]=str(case_root/"numba_cache")
        (case_root/"numba_cache").mkdir()
        runtime=setup_runtime(tax,plan,selected,case_root)
        initial={k:tax.PARAMETERS[k] for k in FREE}
        bounds={r["parameter"]:(float(r["lower"]),float(r["upper"]))
                for r in objective["parameter_restrictions"]}
        assert set(bounds)==set(FREE)
        if stage=="smoke":
            case_path=case_root/"point_00"
            receipt=evaluate_point(tax=tax,objective=objective,selected=selected,runtime=runtime,
                point=initial,output=case_path,deadline_epoch=deadline_epoch,graphs=True)
            record_global(run_root,receipt,case_path)
            write(case_root/"complete.json",dict(status="exact_loop_smoke_passed",loss=receipt["loss"],
                stationary_solves=receipt["objective_stationary_solves"],graph_count=17))
            return
        rng=random.Random(20260924+worker_id)
        start_best=json.loads((run_root/"best_so_far.json").read_text())
        center=dict(start_best["point"]);best=float(start_best["loss"])
        widths=plan["search_config"]["coordinate_widths"]
        for index in range(1,6):
            if deadline_epoch-time.time()<900: break
            point,metadata=proposal(center,bounds,widths,rng,index,worker_id)
            case_path=case_root/f"point_{index:02d}"
            try:
                receipt=evaluate_point(tax=tax,objective=objective,selected=selected,runtime=runtime,
                    point=point,output=case_path,deadline_epoch=deadline_epoch,graphs=False)
            except TimeoutError: raise
            except Exception as exc:
                write(case_path/"failure.json",dict(status="numerical_proposal_rejected",
                    error_type=type(exc).__name__,error=str(exc),point=point,metadata=metadata))
                continue
            saved=json.loads((case_path/"receipt.json").read_text())
            saved["proposal_metadata"]=metadata
            write(case_path/"receipt.json",saved)
            record_global(run_root,receipt,case_path)
            if receipt["loss"]<best:
                best=receipt["loss"];center=point
                write(case_root/"best_so_far.json",dict(loss=best,point=center,case_path=str(case_path)))
            write(case_root/"latest_completed.json",dict(loss=receipt["loss"],point=point,
                case_path=str(case_path),completed=index))
        write(case_root/"complete.json",dict(status="bounded_worker_complete",worker=worker_id,
            proposals_max=5,completed=len(list(case_root.glob("point_*/receipt.json"))),
            deadline_epoch=deadline_epoch))
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)


def export_selected(run_root):
    """Render selected diagnostics and complete tables on Torch before cutoff."""
    deadline=json.loads((run_root/"deadline.json").read_text())["deadline_epoch"]
    if time.time()>=deadline: raise TimeoutError("one-hour calibration/export deadline passed")
    def alarm(_signum,_frame): raise TimeoutError("shared export deadline reached")
    signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,deadline-time.time())
    try:
        best=json.loads((run_root/"best_so_far.json").read_text())
        case=Path(best["case_path"])
        receipt=json.loads((case/"receipt.json").read_text())
        assert receipt["target_system_sha256"]==OBJECTIVE_SHA
        assert receipt["case_checkpoint_sha256"]==sha(case/"initial_state.pkl.gz")
        fits=list(csv.DictReader((case/"target_fit.csv").open()))
        params=list(csv.DictReader((case/"parameters.csv").open()))
        assert len(fits)==13 and sum(r["weight"]!="" for r in fits)==12
        assert len(params)>=17
        out=run_root/"selected_export";out.mkdir(parents=True,exist_ok=False)
        graphs=case/"standard_diagnostics"
        if len(list(graphs.glob("*.png")))!=17:
            tax,plan,selected,objective=prepare_selected()
            os.environ["NUMBA_CACHE_DIR"]=str(out/"numba_cache")
            (out/"numba_cache").mkdir()
            runtime=setup_runtime(tax,plan,selected,out)
            with gzip.open(case/"initial_state.pkl.gz","rb") as f: packet=pickle.load(f)
            assert packet["contract_sha256"]==OBJECTIVE_SHA
            runtime["audit"].standard_diagnostics(packet,out,validate_production_young=False)
            graphs=out/"standard_diagnostics"
        assert len(list(graphs.glob("*.png")))==17
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.image as mpimg
        from matplotlib.backends.backend_pdf import PdfPages
        def fmt(value):
            if value=="": return ""
            try: return f"{float(value):.3g}"
            except (ValueError,TypeError): return str(value)
        def table_page(pdf,title,headers,rows):
            fig,ax=plt.subplots(figsize=(11,8.5));ax.axis("off")
            ax.set_title(title,fontsize=14,loc="left",pad=20)
            cell=[[fmt(v) for v in row] for row in rows]
            t=ax.table(cellText=cell,colLabels=headers,loc="center",cellLoc="left")
            t.auto_set_font_size(False);t.set_fontsize(7);t.scale(1,1.35)
            pdf.savefig(fig,bbox_inches="tight");plt.close(fig)
        pdf_path=out/"experimental_commute_calibration.pdf"
        with PdfPages(pdf_path) as pdf:
            fig,ax=plt.subplots(figsize=(11,8.5));ax.axis("off")
            text=("Experimental one-hour steady-state calibration\n"
                  f"Loss: {receipt['loss']:.3g}; 12 weighted moments; 8 free coordinates\n"
                  "PAYGO payroll tax 0.08751017424959717; fertility normalized to 2.1\n"
                  "New first-birth target 1.465 rooms (head/spouse at baseline).\n"
                  "Model housing observer is an unmatched stationary proxy for the PSID panel estimate.\n"
                  "Model bequest observer counts gross positive estates, including childless households;\n"
                  "the SCF target is child-directed. Recipient mapping remains unresolved.\n"
                  "Adopted 16/20 birth-entry queue is unimplemented; existing stationary law retained.\n"
                  f"Objective SHA256: {OBJECTIVE_SHA}\nSelected case: {case}")
            ax.text(.04,.95,text,va="top",fontsize=11,linespacing=1.5)
            pdf.savefig(fig,bbox_inches="tight");plt.close(fig)
            table_page(pdf,"Complete target fit",("Moment","Target","Model","Gap","Weight","Loss"),
                [(r["moment"],r["target"],r["model"],r["gap"],r["weight"],r["loss_contribution"]) for r in fits])
            for start in (0,14):
                chunk=params[start:start+14]
                if chunk:
                    table_page(pdf,"Parameters and restrictions",("Parameter","Estimate","Lower","Upper","Near bound","Status"),
                        [(r["parameter"],r["estimate"],r["lower"],r["upper"],r["near_bound"],r["status"]) for r in chunk])
            for path in sorted(graphs.glob("*.png")):
                fig,ax=plt.subplots(figsize=(11,8.5));ax.axis("off")
                ax.set_title(path.stem,fontsize=10)
                ax.imshow(mpimg.imread(path));pdf.savefig(fig,bbox_inches="tight");plt.close(fig)
        assert time.time()<deadline
        final=dict(status="selected_export_complete",loss=receipt["loss"],case_path=str(case),
            pdf_path=str(pdf_path),pdf_sha256=sha(pdf_path),graph_count=17,
            target_rows=13,weighted_rows=12,parameter_rows=len(params),
            objective_sha256=OBJECTIVE_SHA,source_sha256=receipt["source_manifest_sha256"],
            selected_checkpoint_sha256=receipt["selected_checkpoint_sha256"],
            elapsed_from_first_smoke=time.time()-json.loads((run_root/"deadline.json").read_text())["first_smoke_start_epoch"])
        write(out/"final_receipt.json",final)
        print(json.dumps(final,sort_keys=True))
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--stage",choices=("preflight","smoke","worker","export"),required=True)
    parser.add_argument("--worker-id",type=int,default=0)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    if args.stage=="preflight": print(json.dumps(core_preflight(),sort_keys=True,indent=2))
    elif args.stage=="export":
        if args.output is None: raise SystemExit("--output required")
        export_selected(args.output.resolve())
    else:
        if args.output is None: raise SystemExit("--output required")
        if args.stage=="worker" and not 1<=args.worker_id<=8: raise SystemExit("worker id 1..8 required")
        run_stage(args.stage,args.worker_id,args.output.resolve())
