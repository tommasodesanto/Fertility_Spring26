#!/usr/bin/env python3
"""Read-only dated-state reconstruction and independent saving-choice audit.

Run from the frozen production source snapshot named in the run plan. The
oracle maximizes the same fixed-continuation objective on every interpolation
segment; it does not change the Bellman solver or recompute an equilibrium.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import json
import math
import pickle
import sys
import threading
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "code/model"), str(ROOT / "code/model/tools")]
import run_e5f_post2023_policy_mechanisms as policy
from intergen_eqscale_seq_optimized import solver as model
from intergen_eqscale_seq_optimized.diagnostics import write_diagnostics
from intergen_eqscale_seq_optimized.kernels import (
    eval_owner_scalar, eval_renter_scalar, njit,
)

baseline, calendar, transition = policy.baseline, policy.calendar, policy.transition
PRODUCTION = "e5f_transition_calibration_eta063_kappa005_rooms_repair_gauss_newton_coord_20260904a"
BUNDLE = "630ba20bca6a1b54eb4c46aca904c4a087afb8c808b9c7f4660d5fcd316a970e"


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_checkpoint(path):
    opener=gzip.open if Path(path).suffix==".gz" else open
    with opener(path,"rb") as stream:
        return pickle.load(stream)


def save_json(path, value):
    calendar.write_json_atomic(Path(path), calendar.jsonable(value))


def progress(out, phase, **extra):
    record = dict(phase=phase, utc=time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()), **extra)
    save_json(out / "heartbeat.json", record)
    print(json.dumps(record), flush=True)


def reconstruct(production_root, out):
    chain, configured = transition.configure_sequential_model()
    base = production_root / "output/model" / PRODUCTION
    task = base / "task_010"
    case = task / "cases/task_010"
    source = production_root / "output/model/intergen_e5f_child_room_floor_psinneg_extended_20260806/report/results.json"
    contracts = baseline.validate_input_contracts(
        report_path=base / "report/summary.json", task_summary_path=task / "summary.json",
        case_dir=case, case_transition_path=case / "transition_path.csv", source_path=source,
        expected_report_sha256="4d9bf21f679f1bd63d2156228f73c672cc39e229489c11e289c319c958df590b",
        expected_task_summary_sha256="7f4ad73e05b5bed319532fae406456aa0f232dd2e47e2d99d8044ac7f5f55fdd",
        expected_case_transition_sha256="1f1e8d5246690a0d7e72ae9085ebc08297820f37dcdafa71dad17f6a1e537863",
        expected_source_sha256="0afcb82d4735bd15aaa143ea04e3105a5d43df152122d02b983372102f20eef6",
        expected_target_fingerprint="3726c17e62c8233ce62d5f4c95f44fd2cc2ea6cfa3d2492795461b4569300497",
        expected_code_bundle_sha256=BUNDLE,
        expected_renewal_contract_sha256="be4df40c258cc92758c75c0f99f32089677c48938d6c5d6507eb77e6dfb61fb3",
        expected_scientific_contract_sha256=None, expected_selection_sha256=None,
        chain=chain, model=configured,
    )
    progress(out, "source_contract_passed", input_hashes=contracts["hashes"])
    prepared = baseline.prepare_model(contracts, source, chain, configured)
    progress(out, "stationary_reconstruction_passed", gates=prepared.stationary_gates)
    history, state, P, mass, reproduction = policy.reconstruct_matched_2023(
        prepared, contracts, market_tol=2e-4, market_max_iter=30, progress_dir=out,
    )
    evaluation, shared, fallback = baseline.evaluate_state(
        state, P, prepared.b_grid, prepared.supply_rule, calendar.SolveCounter(), 2e-4, 30,
    )
    if evaluation.relative_market_residual > 2e-4 or evaluation.feasibility_projection_mass > 2e-10:
        raise RuntimeError("Reconstructed 2023 market/feasibility gate failed")
    progress(out, "dated_reconstruction_passed", reproduction=reproduction,
             market_residual=evaluation.relative_market_residual)
    packet = dict(parameters=P, b_grid=prepared.b_grid, evaluation=evaluation,
                  shared=shared, state=state, supply_rule=prepared.supply_rule,
                  history=history, reproduction=reproduction, input_hashes=contracts["hashes"],
                  driver_sha256=digest(__file__), production_task=str(task))
    checkpoint = out / "dated_state.pkl"
    with checkpoint.open("wb") as stream:
        pickle.dump(packet, stream, protocol=5)
    save_json(out / "checkpoint_manifest.json", dict(path=str(checkpoint), sha256=digest(checkpoint),
              input_hashes=contracts["hashes"], driver_sha256=digest(__file__)))
    return packet


def standard_diagnostics(packet, out, *, validate_production_young=True):
    P, bg, e = packet["parameters"], packet["b_grid"], packet["evaluation"]
    p, g = e.policy, e.g_current
    stats = model.compute_markov_statistics(
        g, p.fert_probs, p.loc_probs, P, bg, p.price, p.hR_pol,
        asset_g=e.g_post_fertility, bequest_g=g, bp_pol=p.bp_pol,
    )
    for key in ("V", "c_pol", "hR_pol", "bp_pol", "tenure_choice", "tenure_probs", "loc_probs", "fert_probs", "fert_value"):
        setattr(stats, key, getattr(p, key))
    stats.fert2_probs = P._fert2_probs
    stats.g, stats.b_grid = g, bg
    stats.type_values,stats.type_weights,_=model.income_transition_values(P)
    stats.owner_asset_price = p.price
    stats.owner_user_cost = P.user_cost_rate * p.price
    stats.housing_demand, stats.housing_supply = e.demand_by_loc, e.supply_by_loc
    stats.aggregate_housing_demand = float(e.demand_by_loc.sum())
    stats.aggregate_housing_supply = float(e.supply_by_loc.sum())
    stats.aggregate_housing_excess = stats.aggregate_housing_demand - stats.aggregate_housing_supply
    stats.best_max_abs_rel_excess = e.relative_market_residual
    stats.rental_demand_by_market = np.sum(g[:, 0] * p.hR_pol[:, 0], axis=(0, 2, 3, 4, 5))
    stats.owner_demand_by_size = np.asarray([float(g[:, t].sum()) * float(P.H_own[t-1]) for t in range(1, g.shape[1])])
    stats.aggregate_rental_demand = float(stats.rental_demand_by_market.sum())
    stats.aggregate_owner_demand = float(stats.owner_demand_by_size.sum())
    if not math.isclose(stats.aggregate_rental_demand + stats.aggregate_owner_demand,
                        stats.aggregate_housing_demand, abs_tol=2e-10, rel_tol=0):
        raise RuntimeError("Diagnostic demand decomposition does not add up")
    quantity_fields=("housing_demand","housing_supply","aggregate_housing_demand",
        "aggregate_housing_supply","aggregate_housing_excess","rental_demand_by_market",
        "owner_demand_by_size","aggregate_rental_demand","aggregate_owner_demand")
    mass=float(g.sum())
    save_json(out/"market_quantity_units.json",dict(adult_households=mass,
        aggregate_quantities={key:getattr(stats,key) for key in quantity_fields},
        plotted_quantity_unit="rooms per adult-household unit; aggregate quantities divided by current household mass"))
    for key in quantity_fields:
        setattr(stats,key,getattr(stats,key)/mass)
    # Preserve the established graph set and layouts. Its older per-adult
    # captions assumed unit household mass; the dated distribution is not
    # normalized to one. Correct both that conversion and its unit label.
    from matplotlib.axes import Axes
    from unittest.mock import patch
    original_ylabel=Axes.set_ylabel
    def dated_ylabel(self,label,*args,**kwargs):
        return original_ylabel(self,"rooms per household" if label=="service units per adult" else label,*args,**kwargs)
    with patch.object(Axes,"set_ylabel",dated_ylabel):
        write_diagnostics(stats, P, out / "standard_diagnostics")
    # Save the standard graph set intact. Supplemental dated flows and exposure
    # definitions are explicit, since inherited graph labels are historical.
    rows = []
    for j in range(P.J):
        gj, pre = g[:, :, :, j], e.g_pre[:, :, :, j]
        mj = float(gj.sum())
        housing = sum(float(np.sum(gj[:, t] * (p.hR_pol[:, t, :, j] if t == 0 else P.H_own[t-1]))) for t in range(g.shape[1]))
        rows.append(dict(age_node=float(P.age_start+j*P.da), mass=mj,
                         owner_rate=float(gj[:, 1:].sum()/mj),
                         mean_rooms=housing/mj,
                         mean_liquid_wealth=float(np.sum(gj * bg[:, None, None, None, None, None])/mj),
                         childless_rate=float(gj[..., 0, :].sum()/mj),
                         pre_fertility_mass=float(pre.sum())))
    baseline.write_csv(out / "lifecycle_2023.csv", rows)
    age_rows = []
    # Uniform-within-cell alternatives are measurements, not changed targets.
    for name, starts in (("cells_start_at_node", np.asarray([r["age_node"] for r in rows])),
                         ("cells_centered_at_node", np.asarray([r["age_node"] for r in rows])-P.da/2)):
        overlap = np.maximum(0, np.minimum(starts+P.da, 35)-np.maximum(starts, 25))/P.da
        den = sum(overlap[j]*rows[j]["mass"] for j in range(P.J))
        own = sum(overlap[j]*rows[j]["mass"]*rows[j]["owner_rate"] for j in range(P.J))/den
        age_rows.append(dict(measurement=name, own_rate_2534=own, interpretation="diagnostic uniform density within each four-year cell"))
    lo, hi = model.age_to_index(P,25), model.age_to_index(P,34)
    den = sum(rows[j]["mass"] for j in range(lo,hi+1))
    own = sum(rows[j]["mass"]*rows[j]["owner_rate"] for j in range(lo,hi+1))/den
    age_rows.insert(0,dict(measurement="production_whole_nodes",own_rate_2534=own,interpretation="unchanged model age_to_index aggregation"))
    if validate_production_young and abs(own - .31098385018672514) > 2e-8:
        raise RuntimeError(f"2023 young ownership reproduction mismatch: {own}")
    baseline.write_csv(out / "young_ownership_age_measurement.csv", age_rows)
    return stats


def budget_audit(packet, out):
    P, bg, e, SD = packet["parameters"], packet["b_grid"], packet["evaluation"], packet["shared"]
    p, g = e.policy, e.g_current
    rows, worst = [], []
    maximum = 0.0
    for j in range(P.J):
        for ten in range(g.shape[1]):
            mass = positive_gap = bad_mass = floor_mass = 0.0
            max_gap = 0.0
            for zz, z in enumerate(P.z_grid):
                for nn in range(P.n_parity):
                    for cs in range(P.n_child_states):
                        idx = (slice(None), ten, 0, j, zz, nn, cs)
                        gm = g[idx]
                        if gm.sum() <= 0:
                            continue
                        cflat = nn + P.n_parity*cs
                        income = model.income_at_state(P,0,j,float(z))
                        Rv = P.R_gross*bg+income
                        Rvt = P.R_gross*np.maximum(bg,0)+income
                        grant = float(SD.gb_flat.reshape(-1)[cflat])
                        Rv += np.clip(grant-Rvt,0,grant)
                        if ten == 0:
                            housing_cost = P.user_cost_rate*p.price[0]*p.hR_pol[idx]
                        else:
                            h = P.H_own[ten-1]
                            housing_cost = (P.delta+P.tau_H)*p.price[0]*h
                            housing_cost += getattr(P,"owner_size_cost",0)*p.price[0]*max(h-getattr(P,"owner_size_cost_ref",6),0)**getattr(P,"owner_size_cost_power",2)
                        gap = p.c_pol[idx] + housing_cost + p.bp_pol[idx] - Rv
                        cb = float(SD.cb_flat.reshape(-1)[cflat])
                        binding = p.c_pol[idx] <= cb+P.c_min+1e-12
                        occupied = gm > 1e-12
                        mass += float(gm.sum())
                        bad_mass += float(gm[gap>1e-9].sum())
                        floor_mass += float(gm[binding].sum())
                        positive_gap += float(np.sum(gm*np.maximum(gap,0)))
                        if np.any(occupied):
                            m = float(np.max(gap[occupied])); max_gap = max(max_gap,m)
                            if m > 1e-9:
                                b = int(np.argmax(np.where(occupied,gap,-np.inf)))
                                worst.append(dict(age=float(P.age_start+j*P.da),tenure=ten,income_index=zz,
                                    parity=nn,dependent_children=cs,wealth=float(bg[b]),mass=float(gm[b]),
                                    consumption=float(p.c_pol[idx][b]),saving=float(p.bp_pol[idx][b]),
                                    excess_expenditure=float(gap[b]),resources=float(Rv[b])))
            maximum=max(maximum,max_gap)
            rows.append(dict(age=float(P.age_start+j*P.da),tenure=ten,mass=mass,
                 budget_excess_mass=bad_mass,budget_excess_share=bad_mass/max(mass,1e-30),
                 consumption_floor_mass=floor_mass,weighted_positive_budget_gap=positive_gap,
                 maximum_occupied_excess=max_gap))
    baseline.write_csv(out / "budget_by_age_tenure.csv",rows)
    if worst:
        baseline.write_csv(out / "worst_occupied_budget_states.csv",sorted(worst,key=lambda r:-r["excess_expenditure"])[:100])
    summary=dict(household_mass=float(g.sum()),
        budget_excess_mass=sum(r["budget_excess_mass"] for r in rows),
        budget_excess_share=sum(r["budget_excess_mass"] for r in rows)/g.sum(),
        weighted_positive_budget_gap=sum(r["weighted_positive_budget_gap"] for r in rows),
        maximum_occupied_excess=maximum,occupancy_threshold_for_maximum=1e-12,budget_tolerance=1e-9,
        scope="exact post-tenure grid mass; no changed policies or equilibrium")
    save_json(out / "budget_summary.json",summary)
    return summary


def policy_array_audit(packet, out):
    """Occupied pre-choice value monotonicity and probability range receipts."""
    e,P=packet["evaluation"],packet["parameters"]
    p=e.policy
    dV=np.diff(p.V,axis=0)
    occupied=e.g_pre[:-1]>1e-12
    drops=occupied & (dV < -1e-7)
    indices=np.argwhere(drops)
    rows=[]
    if len(indices):
        order=np.argsort(dV[drops])[:100]
        for idx in indices[order]:
            b,ten,i,j,z,n,m=map(int,idx)
            rows.append(dict(wealth=float(packet["b_grid"][b]),next_wealth=float(packet["b_grid"][b+1]),
              age=float(P.age_start+P.da*j),tenure=ten,income_index=z,parity=n,dependent_children=m,
              mass=float(e.g_pre[tuple(idx)]),value_change=float(dV[tuple(idx)])))
        baseline.write_csv(out/"value_monotonicity_worst.csv",rows)
    probabilities={}
    for name,arr in (("tenure",p.tenure_probs),("first_birth",p.fert_probs),
                     ("continuation_birth",P._fert2_probs),("location",p.loc_probs)):
        if arr is not None:
            probabilities[name]=dict(minimum=float(np.min(arr)),maximum=float(np.max(arr)),
                         nonfinite=int(np.count_nonzero(~np.isfinite(arr))))
    result=dict(occupied_negative_steps=int(drops.sum()),
        lower_node_mass_at_negative_steps=float(e.g_pre[:-1][drops].sum()),
        share_pre_choice_mass_at_negative_steps=float(e.g_pre[:-1][drops].sum()/e.g_pre.sum()),
        maximum_occupied_value_drop=float(-min(np.min(dV[occupied]),0)),
        value_drop_tolerance=1e-7,probabilities=probabilities,
        interpretation="Adjacent-wealth screen on occupied pre-choice states; fixed-age/income/tenure/family comparisons. A flagged drop requires diagnosis, not automatic attribution to an economic mechanism.")
    save_json(out/"policy_array_summary.json",result)
    return result


@njit(cache=True)
def segment_oracle(lo, hi, resources, continuation, bg, rent, hb, cb, pc,
                   hmax, alpha, oms, beta, es, owner_cost, owner_K, owner):
    """Exhaustive endpoints and stationary points of each concave segment."""
    dc=cb+rent*hb
    cap=rent*(hmax-hb)/(1-alpha)
    hcap=max(hmax-hb,1e-10)
    Kr=(alpha**alpha*((1-alpha)/rent)**(1-alpha))**oms
    points=np.empty(bg.size+3)
    n=0
    points[n]=lo; n+=1
    for x in bg:
        if lo < x < hi:
            points[n]=x; n+=1
    if not owner:
        kink=resources-dc-cap
        if lo < kink < hi:
            points[n]=kink; n+=1
    points[n]=hi; n+=1
    points=np.sort(points[:n])
    best=-1e300; bpbest=lo
    for k in range(n):
        x=points[k]
        if owner:
            value=eval_owner_scalar(x,resources,continuation,bg,owner_cost,cb,pc,owner_K,alpha,oms,beta,es)
        else:
            value=eval_renter_scalar(x,resources,continuation,bg,dc,pc,cap,cb,rent,hmax,hcap,Kr,alpha,oms,beta,es)
        if value>best:
            best=value;bpbest=x
        if k==n-1 or points[k+1]-x<1e-14:
            continue
        mid=(x+points[k+1])/2
        if mid<=bg[0] or mid>=bg[-1]:
            slope=0.0
        else:
            ix=np.searchsorted(bg,mid)-1
            slope=(continuation[ix+1]-continuation[ix])/(bg[ix+1]-bg[ix])
        if slope<=0:
            continue
        if owner:
            optimal_c=(beta*slope/(es*owner_K*alpha))**(1/(alpha*oms-1))
            candidate=resources-owner_cost-cb-optimal_c
        elif resources-dc-mid>cap:
            K=hcap**((1-alpha)*oms)
            optimal_c=(beta*slope/(es*K*alpha))**(1/(alpha*oms-1))
            candidate=resources-cb-rent*hmax-optimal_c
        else:
            optimal_surplus=(beta*slope/(es*Kr))**(1/(oms-1))
            candidate=resources-dc-optimal_surplus
        if x<candidate<points[k+1]:
            if owner:
                value=eval_owner_scalar(candidate,resources,continuation,bg,owner_cost,cb,pc,owner_K,alpha,oms,beta,es)
            else:
                value=eval_renter_scalar(candidate,resources,continuation,bg,dc,pc,cap,cb,rent,hmax,hcap,Kr,alpha,oms,beta,es)
            if value>best:
                best=value;bpbest=candidate
    return bpbest,best


def saving_audit(packet,out,draws):
    P,bg,e,SD=packet["parameters"],packet["b_grid"],packet["evaluation"],packet["shared"]
    p,g=e.policy,e.g_current
    rng=np.random.default_rng(2026090501)
    selected,counts=np.unique(rng.choice(g.size,draws,p=g.reshape(-1)/g.sum()),return_counts=True)
    count_map=dict(zip(selected.tolist(),counts.tolist()))
    # Add deterministic mass-bearing age/tenure boundaries, separate from the
    # probability sample so their unequal inclusion cannot bias its readout.
    for j in range(P.J):
        for ten in range(g.shape[1]):
            gj=g[:,ten,0,j]
            occ=np.argwhere(gj>1e-12)
            if not len(occ): continue
            for r in (occ[0],occ[-1],np.asarray(np.unravel_index(gj.argmax(),gj.shape))):
                b,zz,nn,cs=map(int,r)
                count_map.setdefault(int(np.ravel_multi_index((b,ten,0,j,zz,nn,cs),g.shape)),0)
    coords=[(np.unravel_index(k,g.shape),v) for k,v in count_map.items()]
    coords.sort(key=lambda x:(x[0][3],x[0][4]))
    z_values,_,Pi=model.income_transition_values(P)
    nb,nt,I,_,Nz,npar,ncs=g.shape
    Vbq=np.zeros((nb,nt,I,npar,ncs))
    for ten in range(nt):
        hv=p.price[0]*P.H_own[ten-1] if ten else 0
        for nn in range(npar):
            for cs in range(ncs):
                nk=model.get_completed_fertility(nn,cs,P)
                Vbq[:,ten,0,nn,cs]=model.bequest_utility_vec(bg+hv,nk,P)
    cache_key=None;rows=[];max_gain=0.0
    for k,(idx,count) in enumerate(coords):
        b,ten,i,j,zz,nn,cs=map(int,idx)
        if cache_key!=(j,zz):
            if j==P.J-1: Vnr=Vbq.copy()
            else:
                Vnr=np.zeros_like(Vbq)
                for zn in range(Nz):
                    if Pi[zz,zn]>0:
                        Vnr+=Pi[zz,zn]*p.V[:,:,:,j+1,zn]
                if P.use_age_survival:
                    Vnr=P.survival_probs[j]*Vnr+(1-P.survival_probs[j])*Vbq
            Vc=model.apply_child_aging(Vnr,P,nb,nt,I,npar,ncs,age_index=j)
            cache_key=(j,zz)
        cflat=nn+npar*cs
        cb,hb,pc,gb,al,es=[float(getattr(SD,name).reshape(-1)[cflat]) for name in ("cb_flat","hb_flat","psi_flat","gb_flat","alpha_flat","escale_flat")]
        income=model.income_at_state(P,i,j,float(z_values[zz]))
        resources=P.R_gross*bg[b]+income+np.clip(gb-(P.R_gross*max(bg[b],0)+income),0,gb)
        rent=P.user_cost_rate*p.price[i]
        oms=1-P.sigma
        oc=Ko=0.0
        if ten:
            h=P.H_own[ten-1]
            oc=(P.delta+P.tau_H)*p.price[i]*h+getattr(P,"owner_size_cost",0)*p.price[i]*max(h-getattr(P,"owner_size_cost_ref",6),0)**getattr(P,"owner_size_cost_power",2)
            rawfloor=-SD.phi_choice[i,ten,nn,cs]*p.price[i]*h
            lo=max(float(model.owner_borrowing_floor(P,bg[b],rawfloor,j)),bg[0])
            hi=max(resources-oc-cb-1e-6,lo)
            residual_h=h-getattr(P,"owner_h_bar_scale",1)*hb
            if residual_h<=0:
                raise RuntimeError("Occupied owner state has infeasible housing floor")
            Ko=(max(P.chi,1e-8)*residual_h)**((1-al)*oms)
        else:
            lo=max(float(model.renter_borrowing_floor(P,bg[b],j)),bg[0])
            hi=max(resources-cb-rent*hb-1e-6,lo)
        cont=np.ascontiguousarray(Vc[:,ten,i,nn,cs])
        bp=float(p.bp_pol[idx])
        if bp<lo-1e-9 or bp>hi+1e-9: raise RuntimeError("Saved saving policy outside production feasible interval")
        optimum,value=segment_oracle(lo,hi,resources,cont,bg,rent,hb,cb,pc,P.hR_max,al,oms,P.beta,es,oc,Ko,bool(ten))
        if ten:
            actual=float(eval_owner_scalar(bp,resources,cont,bg,oc,cb,pc,Ko,al,oms,P.beta,es))
        else:
            Kr=(al**al*((1-al)/rent)**(1-al))**oms
            actual=float(eval_renter_scalar(bp,resources,cont,bg,cb+rent*hb,pc,rent*(P.hR_max-hb)/(1-al),cb,rent,P.hR_max,max(P.hR_max-hb,1e-10),Kr,al,oms,P.beta,es))
        gain=value-actual
        if gain < -1e-8: raise RuntimeError("Independent global oracle worse than saved feasible choice")
        max_gain=max(max_gain,gain)
        rows.append(dict(age=float(P.age_start+j*P.da),tenure=ten,income_index=zz,parity=nn,
          dependent_children=cs,wealth=float(bg[b]),mass=float(g[idx]),sample_multiplicity=count,
          saved_saving=bp,global_saving=optimum,saving_gap=optimum-bp,
          saved_objective=actual,global_objective=value,value_gain=gain,
          relative_value_gain=gain/max(abs(actual),1e-12),lo=lo,hi=hi))
        if (k+1)%250==0:
            baseline.write_csv(out / "saving_oracle_latest.csv",rows)
            save_json(out / "best_so_far.json",max(rows,key=lambda r:r["value_gain"]))
            progress(out,"saving_oracle",completed=k+1,total=len(coords),max_value_gain=max_gain)
    baseline.write_csv(out / "saving_oracle_states.csv",rows)
    baseline.write_csv(out / "saving_oracle_worst.csv",sorted(rows,key=lambda r:-r["value_gain"])[:100])
    result=dict(draws=draws,unique_states=len(rows),seed=2026090501,
        maximum_value_gain=max_gain,maximum_absolute_saving_gap=max(abs(r["saving_gap"]) for r in rows),
        sample_share_value_gain_above_1e_6=sum(r["sample_multiplicity"] for r in rows if r["value_gain"]>1e-6)/draws,
        sample_share_value_gain_above_1e_4=sum(r["sample_multiplicity"] for r in rows if r["value_gain"]>1e-4)/draws,
        sample_share_value_gain_above_1e_3=sum(r["sample_multiplicity"] for r in rows if r["value_gain"]>1e-3)/draws,
        sample_share_saving_gap_above_001=sum(r["sample_multiplicity"] for r in rows if abs(r["saving_gap"])>.01)/draws,
        scope="fixed production continuation, exact per-segment maximization; weighted random sample plus separately weighted boundary states; not a globally re-solved equilibrium")
    save_json(out / "saving_oracle_summary.json",result)
    save_json(out / "latest_completed_case.json",dict(status="complete",**result))
    save_json(out / "best_so_far.json",max(rows,key=lambda r:r["value_gain"]))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--production-root",type=Path,required=True)
    parser.add_argument("--outdir",type=Path,required=True)
    parser.add_argument("--checkpoint",type=Path)
    parser.add_argument("--checkpoint-sha256")
    parser.add_argument("--draws",type=int,default=240)
    parser.add_argument("--graphs-only",action="store_true",help="Regenerate the standard packet from a saved checkpoint without solving or sampling.")
    args=parser.parse_args()
    if args.graphs_only and not args.checkpoint:
        raise ValueError("--graphs-only requires an existing checkpoint")
    if not 100<=args.draws<=50000: raise ValueError("draws must be in [100,50000]")
    out=args.outdir.resolve()
    if out.exists() and any(out.iterdir()): raise FileExistsError(out)
    out.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter()
    progress(out,"started",draws=args.draws,driver_sha256=digest(__file__))
    save_json(out / "latest_completed_case.json",dict(status="started",completed=0))
    save_json(out / "best_so_far.json",dict(status="no_audited_saving_states_yet"))
    stop_heartbeat=threading.Event()
    def heartbeat():
        while not stop_heartbeat.wait(60):
            save_json(out / "runtime_heartbeat.json",dict(elapsed_seconds=time.perf_counter()-started,
                      utc=time.strftime("%Y-%m-%d %H:%M:%S UTC",time.gmtime())))
    threading.Thread(target=heartbeat,daemon=True).start()
    if args.checkpoint:
        if digest(args.checkpoint)!=args.checkpoint_sha256: raise RuntimeError("Checkpoint hash mismatch")
        # Only the locally generated, hash-pinned checkpoint from this driver.
        _, configured = transition.configure_sequential_model()
        packet=load_checkpoint(args.checkpoint)
        live = baseline.calibration.code_fingerprint_contract(configured)["bundle_sha256"]
        if live != BUNDLE or packet["input_hashes"]["code_bundle_sha256"] != live:
            raise RuntimeError("Checkpoint or live source contract differs")
    else:
        packet=reconstruct(args.production_root.resolve(),out)
    if args.graphs_only:
        standard_diagnostics(packet,out,validate_production_young=not bool(packet.get("diagnostic_case")))
        save_json(out / "summary.json",dict(status="complete_standard_diagnostics",input_hashes=packet["input_hashes"],
                  checkpoint_sha256=args.checkpoint_sha256,driver_sha256=digest(__file__),model_solves=0))
        stop_heartbeat.set()
        return
    budget=budget_audit(packet,out)
    arrays=policy_array_audit(packet,out)
    progress(out,"budget_audit_complete",summary=budget)
    oracle=saving_audit(packet,out,args.draws)
    standard_diagnostics(packet,out)
    result=dict(status="complete_independent_numerical_audit",elapsed_seconds=time.perf_counter()-started,
        budget=budget,oracle=oracle,policy_arrays=arrays,reproduction=packet["reproduction"],input_hashes=packet["input_hashes"],
        driver_sha256=digest(__file__),production_changed=False)
    save_json(out / "summary.json",result)
    progress(out,"complete",elapsed_seconds=result["elapsed_seconds"])
    stop_heartbeat.set()


if __name__=="__main__": main()
