#!/usr/bin/env python3
"""Diagnostic global-saving comparison, conditional on the certified 2023 state.

The production source files are immutable. This process explicitly substitutes
two saving kernels at runtime, then uses the existing fertility, tenure,
distribution, and market operators. This is not a re-estimated calibration.
"""
from __future__ import annotations
import argparse
import copy
import gzip
import json
import pickle
import threading
import time
from pathlib import Path
import numpy as np
import run_e5f_independent_numerical_audit as audit
from intergen_eqscale_seq_optimized.kernels import njit

model,baseline,calendar,transition=audit.model,audit.baseline,audit.calendar,audit.transition
LOCAL_RENTER=model.full_renter_block_kernel
LOCAL_OWNER=model.full_owner_block_kernel


@njit(cache=True)
def global_renter(Rv,Rvt,Vc,previous,has_prev,bg,cb,hb,psi,gb,alphas,scales,
                  rent,hmax,cmin,cb0,hb0,alpha,oms,beta,snext,Dnext,a1,a2,tol):
    if has_prev: raise ValueError("Global diagnostic requires the full feasible interval")
    nb,nc=Vc.shape
    values=np.empty((nb,nc));savings=np.empty((nb,nc));cons=np.empty((nb,nc));housing=np.empty((nb,nc))
    for c in range(nc):
        al=alphas[c]
        for b in range(nb):
            R=Rv[b]+min(max(gb[c]-Rvt[b],0),gb[c])
            floor=min(snext*min(bg[b],0),-Dnext)
            lo=max(floor,bg[0]);hi=max(R-cb[c]-rent*hb[c]-1e-6,lo)
            bp,val=audit.segment_oracle(lo,hi,R,Vc[:,c],bg,rent,hb[c],cb[c],psi[c],hmax,al,oms,beta,scales[c],0.,0.,False)
            values[b,c]=val;savings[b,c]=bp
            surplus=R-cb[c]-rent*hb[c]-bp
            if surplus<=1e-10:
                cons[b,c]=cb0+cmin;housing[b,c]=hb0+.01
            else:
                ht=(1-al)*surplus/rent
                if hb[c]+ht>hmax:
                    ct=max(R-cb[c]-rent*hmax-bp,1e-10)
                    cons[b,c]=cb[c]+max(ct,cmin)
                    housing[b,c]=hb[c]+max(hmax-hb[c],.01)
                else:
                    cons[b,c]=cb[c]+max(al*surplus,cmin)
                    housing[b,c]=hb[c]+max(ht,.01)
    return values,savings,cons,housing


@njit(cache=True)
def global_owner(Rv,Rvt,Vc,previous,has_prev,bg,cb,hb,psi,gb,alphas,scales,
                 floors,cost,housing,hbscale,premium,cmin,alpha,oms,beta,snext,Dnext,a1,a2,tol,strict=0):
    if has_prev: raise ValueError("Global diagnostic requires the full feasible interval")
    nb,nc=Vc.shape
    values=np.empty((nb,nc));savings=np.empty((nb,nc));cons=np.empty((nb,nc))
    for c in range(nc):
        residual=housing-hbscale*hb[c]
        if strict and residual<=0:
            for b in range(nb):
                values[b,c]=-1e10;savings[b,c]=bg[0];cons[b,c]=cb[c]+cmin
            continue
        Ko=(premium*max(residual,1e-10))**((1-alphas[c])*oms)
        for b in range(nb):
            R=Rv[b]+min(max(gb[c]-Rvt[b],0),gb[c])
            unsecured=min(snext*min(bg[b]-floors[c],0),-Dnext)
            lo=max(floors[c]+unsecured,bg[0]);hi=max(R-cost-cb[c]-1e-6,lo)
            # Rent/housing-cap arguments are unused on the owner branch.
            bp,val=audit.segment_oracle(lo,hi,R,Vc[:,c],bg,1.,0.,cb[c],psi[c],1.,alphas[c],oms,beta,scales[c],cost,Ko,True)
            values[b,c]=val;savings[b,c]=bp
            cons[b,c]=cb[c]+max(R-cost-cb[c]-bp,cmin)
    return values,savings,cons


def set_method(method):
    model.full_renter_block_kernel=GLOBAL_RENTER if method=="global" else LOCAL_RENTER
    model.full_owner_block_kernel=GLOBAL_OWNER if method=="global" else LOCAL_OWNER


GLOBAL_RENTER=global_renter
GLOBAL_OWNER=global_owner


def quantities(e,P):
    g=e.g_current;mass=float(e.g_post_fertility.sum())
    births=transition.calendar_topcode_birth_accounting(e.g_pre,e.g_post_fertility,e.births,P)
    rental=g[:,0];rental_h=e.policy.hR_pol[:,0]
    parent_rental=g[:,0,...,1:]
    rows=dict(adult_households=mass,explicit_births=float(e.births),
        adjusted_births=float(births["topcode_adjusted_birth_children"]),
        births_per_household=float(births["topcode_adjusted_birth_children"])/mass,
        ownership=float(g[:,1:].sum()/g.sum()),
        rooms_per_household=float(e.demand_by_loc.sum()/mass),
        asset_price=float(e.policy.price[0]),market_residual=float(e.relative_market_residual),
        feasibility_projection=float(e.feasibility_projection_mass),
        renter_mass=float(rental.sum()),renter_at_cap_share=float(rental[rental_h>=P.hR_max-1e-8].sum()/rental.sum()),
        dependent_renter_at_cap_share=float(parent_rental[rental_h[...,1:]>=P.hR_max-1e-8].sum()/max(parent_rental.sum(),1e-30)))
    return rows


def run_case(packet,method,name,clear,out,*,fixed_price=None,market_tol=2e-4):
    set_method(method)
    P=copy.deepcopy(packet["parameters"])
    audit.policy.apply_policy(P,audit.policy.POLICIES[name])
    rule=audit.policy.policy_supply_rule(packet["supply_rule"],audit.policy.POLICIES[name])
    state=copy.deepcopy(packet["state"])
    state.initial_policy=None
    state.price_guess=float(packet["evaluation"].policy.price[0] if fixed_price is None else fixed_price)
    counter=calendar.SolveCounter();started=time.perf_counter()
    if clear:
        e,shared,fallback=baseline.evaluate_state(state,P,packet["b_grid"],rule,counter,market_tol,30)
    else:
        shared=model.precompute_shared(P,packet["b_grid"])
        e=calendar.evaluate_period(np.asarray([state.price_guess]),state.g_pre,P,packet["b_grid"],shared,counter,rule)
        fallback=False
    health=calendar.distribution_health(dict(pre=e.g_pre,post=e.g_post_fertility,current=e.g_current))
    if health["nonfinite_distribution_count"] or health["min_distribution_mass"] < -1e-14:
        raise RuntimeError("Distribution health gate failed")
    if abs(e.g_post_fertility.sum()-packet["evaluation"].g_post_fertility.sum())>2e-10:
        raise RuntimeError("Impact common population gate failed")
    if e.feasibility_projection_mass>1e-6: raise RuntimeError("Production feasibility gate failed")
    if clear and e.relative_market_residual>2e-4: raise RuntimeError("Production market gate failed")
    pricing="separately_cleared" if clear else ("fixed_local_policy_price" if fixed_price is not None else "fixed_production_price")
    case_out=out/"cases"/f"{method}_{name}_{pricing}"
    case_out.mkdir(parents=True)
    case_packet=dict(packet,parameters=P,evaluation=e,shared=shared,supply_rule=rule,state=state,
        diagnostic_case=dict(method=method,policy=name,pricing=pricing))
    budget=audit.budget_audit(case_packet,case_out)
    arrays=audit.policy_array_audit(case_packet,case_out)
    audit.standard_diagnostics(case_packet,case_out,validate_production_young=False)
    checkpoint=case_out/"case_state.pkl.gz"
    with gzip.open(checkpoint,"wb",compresslevel=1) as stream:
        pickle.dump(case_packet,stream,protocol=5)
    row=dict(method=method,policy=name,market_cleared=clear,pricing=pricing,**quantities(e,P),
        budget_excess_share=float(budget["budget_excess_share"]),
        weighted_positive_budget_gap=float(budget["weighted_positive_budget_gap"]),
        maximum_occupied_budget_excess=float(budget["maximum_occupied_excess"]),
        reported_budget_exact=bool(budget["budget_excess_mass"]==0),
        requested_market_tolerance=market_tol,diagnostic_directory=str(case_out),
        case_checkpoint_sha256=audit.digest(checkpoint),
        bellman_solves=counter.bellman,elapsed_seconds=time.perf_counter()-started,grid_fallback=fallback)
    if name=="baseline":
        delta=e.policy.V-packet["evaluation"].policy.V
        occupied=packet["evaluation"].g_pre>1e-12
        row.update(min_occupied_value_change=float(delta[occupied].min()),
                   max_occupied_value_change=float(delta[occupied].max()))
    else:
        row.update(min_occupied_value_change=None,max_occupied_value_change=None)
    audit.save_json(out/"latest_completed_case.json",row)
    audit.save_json(case_out/"case_receipt.json",dict(row=row,budget=budget,policy_arrays=arrays,
        inherited_state_checkpoint_sha256=packet["checkpoint_sha256"],
        interpretation="Reported floors are inherited from production; nonzero excess is measured, not declared feasible. Common inherited pre-choice population; prices are described by the pricing field."))
    audit.progress(out,"case_complete",**row)
    return row


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint",type=Path,required=True)
    parser.add_argument("--checkpoint-sha256",required=True)
    parser.add_argument("--outdir",type=Path,required=True)
    parser.add_argument("--mode",choices=("smoke","equilibrium"),required=True)
    parser.add_argument("--smoke-summary",type=Path)
    parser.add_argument("--smoke-summary-sha256")
    parser.add_argument("--market-tol",type=float,default=2e-4)
    args=parser.parse_args();out=args.outdir.resolve()
    if not 0 < args.market_tol <= 2e-4: raise ValueError("Only production or tighter market tolerance is allowed")
    if out.exists() and any(out.iterdir()): raise FileExistsError(out)
    out.mkdir(parents=True,exist_ok=True)
    if audit.digest(args.checkpoint)!=args.checkpoint_sha256: raise RuntimeError("Checkpoint hash failed")
    if args.mode=="equilibrium":
        if not args.smoke_summary or audit.digest(args.smoke_summary)!=args.smoke_summary_sha256:
            raise RuntimeError("Exact-loop smoke receipt missing or mismatched")
        receipt=json.loads(args.smoke_summary.read_text())
        if receipt.get("status")!="complete" or receipt.get("driver_sha256")!=audit.digest(__file__):
            raise RuntimeError("Exact-loop smoke did not pass this driver")
        if not receipt.get("global_dominance_gate_passed"): raise RuntimeError("Global method failed smoke dominance gate")
        if receipt.get("oracle_driver_sha256")!=audit.digest(audit.__file__):
            raise RuntimeError("Oracle or diagnostic helper differs from the smoke")
    transition.configure_sequential_model()
    packet=audit.load_checkpoint(args.checkpoint)
    packet["checkpoint_sha256"]=args.checkpoint_sha256
    _,configured=transition.configure_sequential_model()
    live=baseline.calibration.code_fingerprint_contract(configured)["bundle_sha256"]
    if live!=audit.BUNDLE or packet["input_hashes"]["code_bundle_sha256"]!=live:
        raise RuntimeError("Frozen production source differs")
    started=time.perf_counter();stop=threading.Event()
    def heartbeat():
        while not stop.wait(60):
            audit.save_json(out/"runtime_heartbeat.json",dict(elapsed_seconds=time.perf_counter()-started,
                 utc=time.strftime("%Y-%m-%d %H:%M:%S UTC",time.gmtime())))
    threading.Thread(target=heartbeat,daemon=True).start()
    audit.save_json(out/"latest_completed_case.json",dict(status="no_case_complete"))
    audit.save_json(out/"best_so_far.json",dict(status="diagnostic_comparison_no_calibration_selection"))
    original_solve=calendar.solve_policy
    def traced_solve(price,P,bg,shared,counter):
        start=time.perf_counter()
        result=original_solve(price,P,bg,shared,counter)
        audit.progress(out,"bellman_complete",price=np.asarray(price).tolist(),
                       seconds=time.perf_counter()-start,count=counter.bellman)
        return result
    calendar.solve_policy=traced_solve
    rows=[]
    policies=("baseline",) if args.mode=="smoke" else ("baseline","supply-plus-20","dependent-child-ltv95")
    for method in ("local","global"):
        for name in policies:
            rows.append(run_case(packet,method,name,args.mode=="equilibrium",out,market_tol=args.market_tol))
            baseline.write_csv(out/"case_summary.csv",rows)
    if args.mode=="smoke":
        local,global_=rows
        # At identical prices and unchanged continuation architecture, the
        # exhaustive method must weakly improve every occupied value. This
        # tests the full backward loop, including continuation propagation.
        dominance=global_["min_occupied_value_change"]>=-1e-7
        if not dominance: raise RuntimeError("Global Bellman dominance failed")
        if max(abs(local[k]-quantities(packet["evaluation"],packet["parameters"])[k]) for k in ("births_per_household","ownership","rooms_per_household","asset_price"))>2e-10:
            raise RuntimeError("Local fixed-price identity reproduction failed")
        # Exercise the full market loop with the substituted kernels before
        # authorizing the larger policy panel, in addition to the two
        # common-price Bellman checks above.
        rows.append(run_case(packet,"global","baseline",True,out,market_tol=args.market_tol))
        baseline.write_csv(out/"case_summary.csv",rows)
    else:
        dominance=None  # Values at different equilibrium prices do not have a dominance ordering.
        comparisons=[]
        for method in ("local","global"):
            base=next(r for r in rows if r["method"]==method and r["policy"]=="baseline")
            for r in rows:
                if r["method"]==method and r["policy"]!="baseline":
                    comparisons.append(dict(method=method,policy=r["policy"],
                       births_percent=100*(r["births_per_household"]/base["births_per_household"]-1),
                       rooms_percent=100*(r["rooms_per_household"]/base["rooms_per_household"]-1),
                       ownership_pp=100*(r["ownership"]-base["ownership"]),
                       price_percent=100*(r["asset_price"]/base["asset_price"]-1)))
        baseline.write_csv(out/"policy_comparison.csv",comparisons)
        # Isolate the direct optimizer difference at each local-method policy
        # price. These extra evaluations are conditional-price comparisons;
        # their housing markets are deliberately not cleared again.
        fixed_price_comparisons=[]
        for name in policies:
            local=next(r for r in rows if r["method"]=="local" and r["policy"]==name)
            r=run_case(packet,"global",name,False,out,fixed_price=local["asset_price"],market_tol=args.market_tol)
            rows.append(r)
            baseline.write_csv(out/"case_summary.csv",rows)
            fixed_price_comparisons.append(dict(policy=name,price=local["asset_price"],
                births_global_minus_local_percent=100*(r["births_per_household"]/local["births_per_household"]-1),
                rooms_global_minus_local_percent=100*(r["rooms_per_household"]/local["rooms_per_household"]-1),
                ownership_global_minus_local_pp=100*(r["ownership"]-local["ownership"])))
        baseline.write_csv(out/"fixed_price_method_comparison.csv",fixed_price_comparisons)
    summary=dict(status="complete",mode=args.mode,elapsed_seconds=time.perf_counter()-started,
        global_dominance_gate_passed=dominance,rows=rows,driver_sha256=audit.digest(__file__),
        oracle_driver_sha256=audit.digest(audit.__file__),input_hashes=packet["input_hashes"],
        checkpoint_sha256=args.checkpoint_sha256,production_changed=False,
        market_tolerance=args.market_tol,
        interpretation="Conditional on the original 2023 distribution and inherited history. Only the saving maximization method differs; current prices clear separately in equilibrium mode. This is numerical sensitivity, not a new calibration or confidence interval.")
    audit.save_json(out/"summary.json",summary)
    audit.progress(out,"complete",elapsed_seconds=summary["elapsed_seconds"])
    stop.set();set_method("local")


if __name__=="__main__":main()
