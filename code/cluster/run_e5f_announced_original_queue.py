"""One isolated, announced four-level preference path under the original queue.

All four preference levels are known at 2007.  This wrapper patches only the
in-memory queue evaluator; frozen scientific source and its numerical gates are
left unchanged.  It is intended for the immutable manifest made by
``prepare_e5f_announced_original_queue.py``.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import json
import os
from contextlib import contextmanager
from pathlib import Path
import pickle
import subprocess
import sys
import threading
import time
from types import SimpleNamespace as NS
from unittest.mock import patch

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_key] = "1"


FINAL_PSI = 0.09221854783921073
NATIVE_TOL = 2e-10


def read(path): return json.loads(Path(path).read_text())
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(driver, path, value):
    """Use the frozen atomic serializer, including on failure receipts."""
    driver.save(Path(path), value)


def exact_gap(a, b):
    import numpy as np
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("Nonfinite or differently shaped reproduction objects")
    return float(np.max(np.abs(a-b))) if a.size else 0.


def validate_manifest(m, manifest_path):
    required = {"spec", "endpoint_pickle", "endpoint_receipt", "output", "psi_levels",
                "periods_after_four", "absolute_deadline_unix", "seconds", "smoke_seconds",
                "max_root_evaluations", "file_sha256"}
    missing = required - set(m)
    if missing: raise ValueError("Manifest misses: " + ", ".join(sorted(missing)))
    if list(map(float, m["psi_levels"])) != [0.12891531457859182, 0.11696608375682901,
                                                0.10564290922456478, FINAL_PSI]:
        raise ValueError("The approved announced psi levels are required exactly")
    if int(m["periods_after_four"]) != 100 or int(m["max_root_evaluations"]) != 8:
        raise ValueError("Authorized horizon/root budget is T=104, eight mappings")
    for path, digest in m["file_sha256"].items():
        if sha(path) != digest: raise ValueError("Changed pinned input: " + str(path))
    if not Path(m["output"]).is_absolute(): raise ValueError("Output must be absolute")
    return sha(manifest_path)


def endpoint_from_manifest(m):
    receipt = read(m["endpoint_receipt"])
    with gzip.open(m["endpoint_pickle"], "rb") as f: endpoint = pickle.load(f)
    if not getattr(endpoint, "verified", False) or not receipt.get("verified", False):
        raise ValueError("Endpoint is not verified")
    if endpoint.receipt != receipt: raise ValueError("Endpoint pickle and receipt differ")
    if float(endpoint.parameters.psi_child) != FINAL_PSI: raise ValueError("Endpoint final psi differs")
    return endpoint, receipt


def fresh_endpoint_check(c, endpoint, receipt, folder, deadline):
    """Fresh same-native terminal solve and one-step audit before path injection."""
    import numpy as np
    import e5f_original_queue_terminal as terminal
    candidate = terminal._evaluate_trial(old=c.old, psi=FINAL_PSI,
        coordinates=np.asarray(endpoint.coordinates, dtype=float), audit=c.audit,
        deadline=deadline, trial=1)
    audit = terminal._one_step_audit(candidate, c.old, FINAL_PSI)
    gates = audit.get("checks", {}) if isinstance(audit, dict) else {}
    passed = bool(isinstance(audit, dict) and audit.get("status") == "passed" and gates and all(gates.values()))
    # Policy and endogenous state must belong to the current frozen primitives.
    comparisons = {}
    for name, a, b in (("asset_price", candidate.asset_price, endpoint.asset_price),
                       ("pension", candidate.parameters.pension, endpoint.parameters.pension),
                       ("rebate", candidate.parameters.property_tax_lump_sum_transfer, endpoint.parameters.property_tax_lump_sum_transfer)):
        comparisons[name] = exact_gap(a, b)
        if comparisons[name] > NATIVE_TOL: raise ValueError("Fresh endpoint differs: " + name)
    for name, a, b in (("policy_V", candidate.policy.V, endpoint.policy.V),
                       ("g_pre", candidate.state.g_pre, endpoint.state.g_pre),
                       ("scheduled_entries", candidate.state.scheduled_entries, endpoint.state.scheduled_entries),
                       ("scheduled_raw_entries", candidate.state.scheduled_raw_entries, endpoint.state.scheduled_raw_entries)):
        gap=exact_gap(a,b)
        comparisons[name]=gap
        if gap > NATIVE_TOL: raise ValueError("Fresh endpoint differs: " + name)
    saved_final = receipt.get("final", {}).get("residual")
    if saved_final is None or exact_gap(candidate.residual,saved_final) > NATIVE_TOL:
        raise ValueError("Fresh endpoint residual differs from saved final mapping")
    if not passed: raise ValueError("Fresh one-step stationary audit did not pass")
    save(c.driver, Path(folder)/"fresh_terminal_check.json", dict(passed=True, comparisons=comparisons,
         endpoint_receipt_verified=True, final_psi=FINAL_PSI, one_step_audit=audit,
         household_mass=float(candidate.state.g_pre.sum()), renewal_ratio=float(candidate.renewal_ratio),
         housing_relative_gap=float(candidate.housing_relative_gap),
         pension_relative_gap=float(candidate.pension_relative_gap),rebate_relative_gap=float(candidate.rebate_relative_gap),
         raw_queue=candidate.state.scheduled_raw_entries,adjusted_queue=candidate.state.scheduled_entries))
    return candidate


def announced_queue_path(c, psi_path):
    """Queue evaluator with a full known vector passed directly to native PF."""
    import numpy as np
    rebated, queue, pf = c.rebated, c.queue, c.joined.pf
    psi_path = np.asarray(psi_path, dtype=float)
    if psi_path.ndim != 1 or not len(psi_path) or not np.isfinite(psi_path).all():
        raise ValueError("Finite one-dimensional announced preferences required")
    def evaluate(*, inherited, old_state, prices, pensions, transfers, psi, terminal,
                 observer=None, demographics=None, demographic_evaluator=None):
        del demographics, demographic_evaluator
        p, benefits, rebates = queue._validated_paths(prices, pensions, transfers)
        if len(p) != len(psi_path): raise ValueError("Announced vector/horizon mismatch")
        active_psi = psi_path.copy()
        terminal_P, terminal_policy, terminal_price, _ = rebated._terminal_parts(terminal)
        if not np.isclose(float(terminal_P.psi_child), float(active_psi[-1]), rtol=0, atol=1e-14):
            raise ValueError("Terminal must use final announced preference")
        observed=0; native_period=pf.calendar.evaluate_period
        def dated(*args, **kwargs):
            nonlocal observed
            if observed >= len(p) or not np.isclose(float(args[2].psi_child), active_psi[observed], rtol=0, atol=1e-14):
                raise RuntimeError("Native household preference differs from announced date")
            e=native_period(*args, **kwargs)
            if observer is not None: observer(observed,e,args[2],args[3],args[4])
            observed += 1; return e
        with patch.object(pf.calendar, "evaluate_period", dated):
            native=pf.evaluate_path_at_prices(prices=p, psi_path=active_psi, transfer_path=rebates,
                terminal_price=float(terminal_price), terminal_V=terminal_policy.V, base_parameters=old_state.parameters,
                b_grid=old_state.b_grid, initial_state=inherited.households, supply_rule=old_state.supply_rule,
                birth_to_entry_conversion=queue.BIRTH_TO_ENTRY_CONVERSION, historical_conditioning=None,
                pension_path=benefits, payroll_tax_path=np.full(len(p), queue.PAYROLL_TAX))
        if observed != len(p): raise RuntimeError("Native announced evaluation skipped a date")
        if len(native.rows) != len(p) or len(native.values) != len(p)+1:
            raise RuntimeError("Native announced rows/values horizon mismatch")
        if any(not np.isfinite(v) or v > tol for v,tol in (
                (native.maximum_mass_accounting_error,2e-8),
                (native.maximum_policy_reproduction_error,2e-10),
                (native.maximum_feasibility_projection_mass,1e-6))):
            raise RuntimeError("Announced native path failed original queue numerical gates")
        for i,row in enumerate(native.rows):
            row.update(period=i, calendar_year=int(inherited.year)+4*i, announced_psi=float(active_psi[i]),
                       annual_net_migration_over_period=0., net_migrant_heads_over_period=0.,
                       pension_period=float(benefits[i]), pension_period_units=float(benefits[i]))
            if "psi_child" not in row or not np.isclose(float(row["psi_child"]), active_psi[i], rtol=0, atol=1e-14):
                raise RuntimeError("Native row psi differs from announced schedule")
            row.update(queue.annotate_original_queue_metadata())
        return NS(history=NS(rows=[],values=[],bellman_solves=0), person_tail=native, rows=native.rows,
            values=native.values, bellman_solves=native.bellman_solves,
            maximum_market_residual=native.maximum_market_residual,
            maximum_mass_accounting_error=native.maximum_mass_accounting_error,
            maximum_policy_reproduction_error=native.maximum_policy_reproduction_error,
            maximum_feasibility_projection_mass=native.maximum_feasibility_projection_mass,
            elapsed_seconds=native.elapsed_seconds, pension_period=benefits.copy(), initial_2023_age_head_gap=0.,
            metadata=queue.annotate_original_queue_metadata())
    def first(*, inherited, old_state, demographics, path, prices, pensions, transfers, psi, demographic_evaluator=None):
        # First accepted date must replay only psi_0, while its V reflects all later news.
        p,b,t=queue._validated_paths(prices,pensions,transfers)
        if len(p)<2 or len(path.values)<2 or not path.rows:
            raise ValueError("First replay needs two prices/values and a dated row")
        one_evaluate,_=announced_queue_path(c, np.asarray(psi_path[:1]))
        replay_terminal=NS(parameters=NS(psi_child=float(psi_path[0])),policy=NS(V=np.asarray(path.values[1])),asset_price=float(p[1]))
        one=one_evaluate(inherited=inherited,old_state=old_state,prices=p[:1],pensions=b[:1],transfers=t[:1],
                     psi=float(psi_path[0]),terminal=replay_terminal,observer=None)
        if exact_gap(one.values[0],path.values[0]) > NATIVE_TOL: raise RuntimeError("First announced replay V differs")
        for key in ("asset_price","renter_price","housing_demand","housing_supply","owner_rate",
                    "birth_children_topcode_adjusted","pension_period_units","payroll_tax_revenue",
                    "pension_outlays","property_tax_revenue","equal_transfer_outlays",
                    "effective_mature_entrant_flow_B","raw_state_scheduled_mature_entrant_flow_B",
                    "entrant_flow_next","mass_accounting_residual","psi_child"):
            expected,actual=path.rows[0].get(key),one.rows[0].get(key)
            if ((expected is None)!=(actual is None) or
                    (expected is not None and exact_gap(expected,actual)>NATIVE_TOL)):
                raise RuntimeError("First announced replay row differs: "+key)
        return c.rebated.InheritedState(int(inherited.year)+4,one.person_tail.terminal_state)
    return evaluate, first


def initial_guess(c, endpoint, m, count):
    import numpy as np
    if count != 104: raise ValueError("Saved seed extension requires exactly 104 dates")
    seed=read(m["initial_seed_json"]); x=np.asarray(seed["prices"],dtype=float)
    if x.size != 300: raise ValueError("Seed must contain 300 coordinates")
    base=x.reshape(3,100); target=np.asarray(endpoint.coordinates,dtype=float)
    ext=np.empty((3,count)); ext[:,:100]=base
    for j,w in enumerate((.2,.4,.6,.8),start=100): ext[:,j]=(1-w)*base[:,-1]+w*target
    return ext


def terminal_distance(actual, target):
    import numpy as np
    mass=max(float(target.g_pre.sum()),1e-15)
    return dict(distribution_relative_l1=float(np.abs(actual.g_pre-target.g_pre).sum()/mass),
      population_relative_gap=float(abs(actual.g_pre.sum()/mass-1)),
      queue_relative_max=float(np.max(np.abs(np.asarray(actual.scheduled_entries)/np.asarray(target.scheduled_entries)-1))),
      raw_queue_relative_max=float(np.max(np.abs(np.asarray(actual.scheduled_raw_entries)/np.asarray(target.scheduled_raw_entries)-1))))


def make_references(c, endpoint):
    import numpy as np
    import run_e5f_transition_calibration as fertility
    import run_e5f_original_queue_experiments as runner
    calendar=c.primitive.calendar
    def one(end):
        P=end.parameters; shared=calendar.model.precompute_shared(P,c.old.b_grid)
        e=calendar.evaluate_period(np.array([end.asset_price]),end.state.g_pre,P,c.old.b_grid,
            shared,calendar.SolveCounter(),supply_rule=c.old.supply_rule,supplied_policy=end.policy)
        return dict(fertility=fertility.period_fertility_diagnostics(e,P),quantities=runner.reference(c,e,P))
    initial=NS(parameters=c.old.parameters,policy=c.old.policy,
               asset_price=float(c.old.policy.price[0]),state=c.old.initial_state)
    return dict(initial=one(initial),terminal=one(endpoint),stationary_endpoint_verified=True)


def plot_packet(c,m,folder,references,*,status):
    save(c.driver,Path(folder)/"stationary_reference.json",references)
    save(c.driver,Path(folder)/"irf_contract.json",dict(label="Announced four-step preference path",
        shock_calendar_year=2007,status_label=status,
        shock_description="Four levels known in 2007; final level permanent.",
        production_eligible=False))
    subprocess.run([sys.executable,m["plotter"],"--case-dir",str(folder),
                    "--fertility-rate","--include-pre-shock"],check=True,timeout=120)


def run_path(c, endpoint, m, folder, deadline, *, psi_path=None, guess=None):
    import numpy as np
    psi_path=np.r_[np.asarray(m["psi_levels"],float),np.full(100,FINAL_PSI)] if psi_path is None else np.asarray(psi_path,float)
    count=len(psi_path)
    if count < 2: raise ValueError("At least two announced dates required")
    controls=dict(c.controls); [controls.pop(k,None) for k in ("automatic_fiscal_polish","fiscal_tolerance","fiscal_slope","initial_jacobian")]
    controls["slope"]=controls.pop("market_slope",1.63); controls["max_evaluations"]=8
    guess=initial_guess(c,endpoint,m,count) if guess is None else np.asarray(guess,float); observations=[]; snapshot={}
    references=make_references(c,endpoint)
    def observe(i,e,P,grid,shared):
        import run_e5f_transition_calibration as fertility
        observations.append(dict(period=i,calendar_year=2007+4*i,announced_psi=float(psi_path[i]),**fertility.period_fertility_diagnostics(e,P)))
        if i==0:snapshot.update(parameters=P,b_grid=grid,evaluation=e,shared=shared,supply_rule=c.old.supply_rule)
    def progress(row): save(c.driver,Path(folder)/"latest_completed.json",row); save(c.driver,Path(folder)/"best_so_far.json",row) if row.get("new_best") else None
    mapping=[0]; native=c.rebated.evaluate_forecast
    def capture(**kwargs):
        observations.clear(); snapshot.clear()
        result=native(**kwargs); mapping[0]+=1
        packet=Path(folder)/"mappings"/f"mapping_{mapping[0]:02d}"
        save(c.driver,packet/"rows.json",result.rows)
        save(c.driver,packet/"fertility.json",observations)
        save(c.driver,packet/"terminal_distance.json",terminal_distance(result.person_tail.terminal_state,endpoint.state))
        save(c.driver,packet/"native_queue_gates.json",dict(mass=result.maximum_mass_accounting_error,
             policy=result.maximum_policy_reproduction_error, feasibility=result.maximum_feasibility_projection_mass,
             rows=len(result.rows), values=len(result.values)))
        save(c.driver,Path(folder)/"latest_completed_mapping.json",dict(mapping=mapping[0],path=str(packet)))
        if count==104:
            plot_packet(c,m,packet,references,status="Unconverged mapping; see market and terminal-distance checks.")
            if snapshot:
                from run_e5f_successive_surprises_overnight import standard_graphs
                standard_graphs(snapshot,NS(path=result),packet/"graphs")
        return result
    with patch.object(c.rebated,"evaluate_forecast",capture):
      result=c.rebated.solve_rebated_forecast(inherited=c.rebated.InheritedState(2007,c.old.initial_state),psi=float(psi_path[-1]),
        old_state=c.old,terminal=NS(parameters=endpoint.parameters,policy=endpoint.policy,asset_price=endpoint.asset_price),
        demographic_primitives=None,count=count,initial_prices=guess[0],initial_pensions=guess[1],initial_transfers=guess[2],
        audit_controls=c.audit,root_controls=controls,deadline_monotonic=deadline,callback=progress,observer=observe)
    receipt=result.root_receipt; receipt.update(schema="e5f_announced_original_queue_v1", expected_psi_path=psi_path.tolist(),
      information="All four announced preference levels are known from 2007; no later surprise occurs.",
      fiscal_information="PAYGO payroll tax 0.179 and equal 1% property-tax rebate each date.")
    receipt.update(c.queue.annotate_original_queue_metadata())
    save(c.driver,Path(folder)/"root_receipt.json",receipt)
    if result.path is not None:
        save(c.driver,Path(folder)/"rows.json",result.path.rows); save(c.driver,Path(folder)/"fertility.json",observations)
        save(c.driver,Path(folder)/"terminal_distance.json",terminal_distance(result.path.person_tail.terminal_state,endpoint.state))
        if count==104:
            finite=bool(receipt.get("finite_horizon_market_fiscal_converged"))
            plot_packet(c,m,folder,references,status=("Market/fiscal root passed; terminal approach unverified."
                if finite else "Unconverged path; see root and terminal-distance receipts."))
            if snapshot:
                from run_e5f_successive_surprises_overnight import standard_graphs
                standard_graphs(snapshot,result,Path(folder)/"graphs")
    return result


def compare_native_paths(a,b):
    gaps={"values":exact_gap(a.values,b.values)}
    for key in ("g_pre","scheduled_entries","scheduled_raw_entries"):
        gaps[key]=exact_gap(getattr(a.person_tail.terminal_state,key),getattr(b.person_tail.terminal_state,key))
    if len(a.rows)!=len(b.rows): raise ValueError("Constant-path row counts differ")
    row_gap=0.
    for left,right in zip(a.rows,b.rows):
        for key,value in left.items():
            if key not in right: raise ValueError("Native constant row missing "+key)
            if isinstance(value,(int,float)) and not isinstance(value,bool):
                row_gap=max(row_gap,exact_gap(value,right[key]))
            elif value != right[key]: raise ValueError("Native constant metadata differs: "+key)
    gaps["rows"]=row_gap
    if max(gaps.values())>NATIVE_TOL: raise ValueError("Constant-vector/native mapping mismatch")
    return gaps


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--manifest",type=Path,required=True); ap.add_argument("--mode",choices=("smoke","run"),required=True); a=ap.parse_args()
    m=read(a.manifest); manifest_sha=validate_manifest(m,a.manifest); out=Path(m["output"]); folder=out/a.mode
    spec=read(m["spec"]); sys.path.insert(0,str(Path(spec["batch"])/"source"))
    import run_e5f_original_queue_experiments as runner
    # Import only after immutable hashes, then load_context installs frozen imports/pins.
    c=runner.load_context(m["spec"]); c.spec_path=Path(m["spec"]); endpoint,receipt=endpoint_from_manifest(m)
    if folder.exists(): raise ValueError("Refusing to overwrite a started mode")
    folder.mkdir(parents=True); end=min(float(m["absolute_deadline_unix"]),time.time()+float(m["smoke_seconds"] if a.mode=="smoke" else m["seconds"]))
    deadline=time.monotonic()+end-time.time(); save(c.driver,folder/"startup_receipt.json",dict(manifest_sha256=manifest_sha,mode=a.mode,deadline_unix=end))
    stop=threading.Event()
    def heartbeat():
        while not stop.wait(60):
            save(c.driver,folder/"controller_heartbeat.json",dict(remaining_seconds=max(0.,deadline-time.monotonic()),mode=a.mode))
            if time.monotonic() >= deadline:
                save(c.driver,folder/"controller_failure.json",dict(error="Authorized deadline reached"))
                os._exit(124)
    threading.Thread(target=heartbeat,daemon=True).start()
    try:
      # This is deliberately outside announced injection.
      fresh_endpoint_check(c,endpoint,receipt,folder,min(deadline,time.monotonic()+600))
      if a.mode=="run":
        smoke=read(out/"smoke_summary.json")
        if not (smoke.get("passed") and smoke.get("manifest_sha256")==manifest_sha): raise ValueError("Matching native smoke required")
      psi_path=__import__("numpy").r_[m["psi_levels"],__import__("numpy").full(100,FINAL_PSI)]
      evaluate,first=announced_queue_path(c,psi_path)
      with c.queue.original_queue_adapter(), patch.object(c.rebated,"evaluate_forecast",evaluate), patch.object(c.rebated,"first_period_state",first), c.cache.policy_cache(c.joined.pf,max_bytes=12*1024**3):
        if a.mode=="smoke":
          # (a,b) Six-date exact stationary mapping and joint root, with a constant
          # announced vector.  It checks vector routing and accepted-date replay.
          P=c.old.parameters; q=float(c.packet["evaluation"].policy.price[0])
          stationary=NS(parameters=P,policy=c.packet["evaluation"].policy,asset_price=q,state=c.old.initial_state)
          constant=__import__("numpy").full(6,float(P.psi_child))
          ev0, fp0=announced_queue_path(c,constant)
          g=__import__("numpy").repeat([q,float(P.pension),float(P.property_tax_lump_sum_transfer)],6).reshape(3,6)
          inputs=dict(inherited=c.rebated.InheritedState(2007,c.old.initial_state),old_state=c.old,
              prices=g[0],pensions=g[1],transfers=g[2],psi=float(P.psi_child),
              terminal=NS(parameters=P,policy=stationary.policy,asset_price=q))
          baseline=c.queue.queue_path(**inputs)
          routed=ev0(**inputs)
          save(c.driver,folder/"constant_native_comparison.json",dict(passed=True,**compare_native_paths(baseline,routed)))
          with patch.object(c.rebated,"evaluate_forecast",ev0), patch.object(c.rebated,"first_period_state",fp0):
            stationary_result=run_path(c,stationary,m,folder/"constant_root",deadline,psi_path=constant,guess=g)
          if stationary_result.next_state is None: raise RuntimeError("Constant-vector six-date root/replay failed")
          drift=terminal_distance(stationary_result.path.person_tail.terminal_state,c.old.initial_state)
          if max(drift.values())>1e-5: raise RuntimeError("Constant-vector stationary drift failed")
          save(c.driver,folder/"constant_stationary_drift.json",dict(passed=True,**drift))
          # (c) One direct six-date mapping containing the announced four-level
          # transition and two final-psi dates; no equilibrium claim is made.
          announced=__import__("numpy").r_[m["psi_levels"],FINAL_PSI,FINAL_PSI]
          g=__import__("numpy").repeat(__import__("numpy").asarray(endpoint.coordinates,float)[:,None],6,axis=1)
          evA,fpA=announced_queue_path(c,announced)
          from e5f_balanced_terminal import _household_checks
          import run_e5f_transition_calibration as fertility
          audits=[]; observations=[]
          rents=c.joined.pf.rents_from_asset_prices(g[0],endpoint.asset_price,c.old.parameters)
          def audit_announced(i,e,P,grid,shared):
              diagnostics,gates=_household_checks(e,P,shared,grid,float(rents[i]),c.primitive,c.audit)
              if not gates or not all(gates.values()): raise RuntimeError("Announced dated household audit failed")
              if P.pension!=g[1,i] or P.property_tax_lump_sum_transfer!=g[2,i] or P.tau_pay!=.179:
                  raise RuntimeError("Announced dated fiscal parameters failed")
              audits.append(dict(period=i,psi=float(P.psi_child),diagnostics=diagnostics,gates=gates))
              observations.append(dict(period=i,calendar_year=2007+4*i,**fertility.period_fertility_diagnostics(e,P)))
          direct=evA(inherited=c.rebated.InheritedState(2007,c.old.initial_state),old_state=c.old,prices=g[0],pensions=g[1],transfers=g[2],psi=FINAL_PSI,
             terminal=NS(parameters=endpoint.parameters,policy=endpoint.policy,asset_price=endpoint.asset_price),observer=audit_announced)
          fpA(inherited=c.rebated.InheritedState(2007,c.old.initial_state),old_state=c.old,demographics=None,
              path=direct,prices=g[0],pensions=g[1],transfers=g[2],psi=FINAL_PSI)
          ddir=folder/"announced_direct"; save(c.driver,ddir/"rows.json",direct.rows); save(c.driver,ddir/"native_values_shape.json",dict(values=len(direct.values),rows=len(direct.rows)))
          save(c.driver,ddir/"household_audits.json",audits); save(c.driver,ddir/"fertility.json",observations)
          save(c.driver,ddir/"terminal_distance.json",terminal_distance(direct.person_tail.terminal_state,endpoint.state))
          if len(direct.rows)!=6 or len(direct.values)!=7: raise RuntimeError("Six-date announced mapping shape failed")
          if [r["calendar_year"] for r in direct.rows] != [2007+4*i for i in range(6)]: raise RuntimeError("Announced dates wrong")
          if [r["announced_psi"] for r in direct.rows] != announced.tolist(): raise RuntimeError("Announced psi rows wrong")
          plot_packet(c,m,ddir,make_references(c,endpoint),status="Six-date smoke only; equilibrium not solved.")
        else: result=run_path(c,endpoint,m,folder,deadline)
      c.driver.verify_pins(m["file_sha256"])
      if a.mode=="smoke": save(c.driver,out/"smoke_summary.json",dict(passed=True,manifest_sha256=manifest_sha,native_gates=True,final_endpoint_fresh_check=True))
      save(c.driver,folder/"controller_complete.json",dict(completed=True,mode=a.mode,
          finite_horizon_market_fiscal_converged=(bool(result.root_receipt.get("finite_horizon_market_fiscal_converged")) if a.mode=="run" else None),
          stationary_endpoint_verified=True,horizon_verified=False,production_eligible=False))
    except BaseException as exc:
      save(c.driver,folder/"controller_failure.json",dict(error_type=type(exc).__name__,error=str(exc)))
      raise
    finally: stop.set()

if __name__=="__main__": main()
