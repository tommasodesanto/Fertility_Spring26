#!/usr/bin/env python3
"""Two-node fixed-continuation branch audit; no Bellman/GE/population solve.

Independent scalar objective plus derivative bisection on every linear-value
segment. Run only after lead review, with one CPU/16GB and a five-minute cap.
"""
from __future__ import annotations
import os
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_key] = "1"
import argparse
import bisect
import gzip
import hashlib
import json
import math
from pathlib import Path
import pickle
import sys
import time


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def verify(path, expected):
    if digest(path) != expected:
        raise RuntimeError(f"Input SHA256 mismatch: {path}")


def interp(grid, values, x):
    if x <= grid[0]:
        return values[0]
    if x >= grid[-1]:
        return values[-1]
    i = bisect.bisect_right(grid, x) - 1
    w = (x-grid[i])/(grid[i+1]-grid[i])
    return (1-w)*values[i]+w*values[i+1]


def allocation(bp, c):
    if c['owner']:
        consumption = c['resources']-c['owner_cost']-c['cb']-bp
        housing = c['owner_effective_housing']
        capped = True
        surplus = consumption
    else:
        surplus = c['resources']-c['cb']-c['rent']*c['hb']-bp
        housing = min((1-c['alpha'])*surplus/c['rent'], c['hmax']-c['hb'])
        capped = (1-c['alpha'])*surplus/c['rent'] > c['hmax']-c['hb']
        consumption = (c['resources']-c['cb']-c['rent']*c['hmax']-bp
                       if capped else c['alpha']*surplus)
    return consumption, housing, surplus, capped


def objective(bp, c):
    consumption, housing, surplus, _ = allocation(bp, c)
    if consumption <= 1e-10 or housing <= 0 or surplus <= 1e-10:
        return -1e10
    composite = consumption**c['alpha'] * housing**(1-c['alpha'])
    return (c['es']*composite**c['oms']/c['oms'] + c['pc']
            + c['beta']*interp(c['grid'], c['continuation'], bp))


def derivative(bp, slope, c):
    consumption, housing, surplus, capped = allocation(bp, c)
    if consumption <= 0 or housing <= 0 or surplus <= 0:
        return -math.inf
    scaled_power = c['es']*(consumption**c['alpha']*housing**(1-c['alpha']))**c['oms']
    marginal = c['alpha']*scaled_power/consumption if capped else scaled_power/surplus
    return -marginal+c['beta']*slope


def maximize(c):
    """All segment endpoints plus roots of the decreasing segment derivative."""
    knots = [c['lo'], c['hi']] + [b for b in c['grid'] if c['lo'] < b < c['hi']]
    if not c['owner']:
        kink = c['resources']-c['cb']-c['rent']*c['hb']-c['rent']*(c['hmax']-c['hb'])/(1-c['alpha'])
        if c['lo'] < kink < c['hi']:
            knots.append(kink)
    knots = sorted(set(knots))
    candidates = [(b, objective(b, c), 'endpoint') for b in knots]
    for a, b in zip(knots[:-1], knots[1:]):
        midpoint = (a+b)/2
        ix = bisect.bisect_right(c['grid'], midpoint)-1
        slope = (0. if midpoint <= c['grid'][0] or midpoint >= c['grid'][-1]
                 else (c['continuation'][ix+1]-c['continuation'][ix])/(c['grid'][ix+1]-c['grid'][ix]))
        if derivative(a, slope, c) <= 0 or derivative(b, slope, c) >= 0:
            continue
        left, right = a, b
        for _ in range(70):
            mid = (left+right)/2
            if derivative(mid, slope, c) > 0:
                left = mid
            else:
                right = mid
        root = (left+right)/2
        candidates.append((root, objective(root, c), 'stationary'))
    best = max(candidates, key=lambda x:x[1])
    return dict(saving=best[0], value=best[1], candidate_kind=best[2],
                segment_count=max(len(knots)-1,0), candidate_count=len(candidates))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--contract', type=Path, required=True)
    ap.add_argument('--contract-sha256', required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    verify(args.contract, args.contract_sha256)
    contract = json.loads(args.contract.read_text())
    verify(__file__, contract['driver_sha256'])
    root = Path(contract['snapshot'])
    verify(contract['parent_contract'], contract['parent_contract_sha256'])
    parent = json.loads(Path(contract['parent_contract']).read_text())
    for relative, pin in parent['source_sha256'].items():
        path = (root/relative).resolve()
        if not path.is_relative_to(root):
            raise RuntimeError('Source path escapes snapshot')
        verify(path, pin)
    verify(contract['checkpoint'], contract['checkpoint_sha256'])
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    sys.path[:0] = [str(root/'code/model'), str(root/'code/model/tools')]
    import numpy as np
    import run_e5f_open_population_transition as transition
    _, model = transition.configure_sequential_model()
    with gzip.open(contract['checkpoint'], 'rb') as stream:
        packet = pickle.load(stream)
    P, e, SD = packet['parameters'], packet['evaluation'], packet['shared']
    p = e.policy
    grid = np.asarray(packet['b_grid']).reshape(-1)
    j, zz, nn, cs = 3, 8, 0, 0
    assert P.age_start+j*P.da == 30 and int(P.I)==1
    assert not P.joint_nested_choice and P.exhaustive_saving_control
    assert not P.use_pti_constraint and str(P.interp_method)=='linear'
    assert float(P.hbar_child_rooms)==0 and float(P.sigma)==2
    assert abs(float(P.pension)-2.0463613896121218)<1e-12
    assert abs(float(P.psi_child)-.2900515293650047)<1e-12
    assert abs(float(P.hbar_first_child_jump)-.7168500088318734)<1e-12
    assert list(np.asarray(P.H_own)) == [2,4,6,8,10]
    assert np.all(np.asarray(SD.birth_entry_grant)==0) and not np.any(SD.birth_dp)
    z, _, Pi = model.income_transition_values(P)
    assert abs(float(z[zz])-1.1432279552449633)<1e-12
    nb, nt, I, _, nz, npar, ncs = p.V.shape
    # Reproduce the solver's expectation, death bequest and child-aging order.
    Vbq = np.zeros((nb,nt,I,npar,ncs))
    for ten in range(nt):
        hv = float(p.price[0])*float(P.H_own[ten-1]) if ten else 0.
        for n in range(npar):
            for m in range(ncs):
                Vbq[:,ten,0,n,m] = model.bequest_utility_vec(grid+hv,model.get_completed_fertility(n,m,P),P)
    expected = np.zeros_like(Vbq)
    for zn in range(nz):
        if Pi[zz,zn]>0:
            expected += Pi[zz,zn]*p.V[:,:,:,j+1,zn,:,:]
    survival = float(P.survival_probs[j]) if P.use_age_survival else 1.
    continuation = model.apply_child_aging(survival*expected+(1-survival)*Vbq,
        P,nb,nt,I,npar,ncs,age_index=j)
    cb,hb,pc,grant,alpha,es = [float(getattr(SD,k).reshape(-1)[0]) for k in
        ('cb_flat','hb_flat','psi_flat','gb_flat','alpha_flat','escale_flat')]
    assert cb==0 and hb==0 and pc==0 and es==1
    price, rent = float(p.price[0]), float(P.user_cost_rate*p.price[0])
    income = float(model.income_at_state(P,0,j,float(z[zz])))
    base = dict(grid=grid.tolist(),cb=cb,hb=hb,pc=pc,alpha=alpha,es=es,
        beta=float(P.beta),oms=1-float(P.sigma),rent=rent,hmax=float(P.hR_max))
    contexts, records, lookups = {}, [], {}
    count = 0
    def audit(ten, wealth, grid_index=None):
        nonlocal count
        key = (ten,float(wealth))
        if key in lookups:
            return lookups[key]
        count += 1
        if count > contract['maximum_scalar_branch_optimizations']:
            raise RuntimeError('Scalar branch budget exceeded')
        c = dict(base,owner=ten>0,continuation=continuation[:,ten,0,0,0].tolist())
        test_resources = float(P.R_gross)*max(wealth,0)+income
        c['resources'] = float(P.R_gross)*wealth+income+min(max(grant-test_resources,0),grant)
        c['owner_cost'], c['owner_effective_housing'] = 0.,0.
        if ten:
            h = float(P.H_own[ten-1]); phi=float(SD.phi_choice[0,ten,0,0])
            c['owner_cost'] = ((float(P.delta)+float(P.tau_H))*price*h
                +float(getattr(P,'owner_size_cost',0))*price*max(h-float(getattr(P,'owner_size_cost_ref',6)),0)**float(getattr(P,'owner_size_cost_power',2)))
            c['owner_effective_housing'] = max(float(P.chi),1e-8)*(h-float(getattr(P,'owner_h_bar_scale',1))*hb)
            c['lo'] = max(float(model.owner_borrowing_floor(P,wealth,-phi*price*h,j)),float(grid[0]))
            c['hi'] = max(c['resources']-c['owner_cost']-cb-1e-6,c['lo'])
        else:
            c['lo'] = max(float(model.renter_borrowing_floor(P,wealth,j)),float(grid[0]))
            c['hi'] = max(c['resources']-cb-rent*hb-1e-6,c['lo'])
        optimum = maximize(c)
        record = dict(tenure=ten,current_branch_wealth=wealth,grid_index=grid_index,
            optimum=optimum,lo=c['lo'],hi=c['hi'],resources=c['resources'])
        if grid_index is not None:
            saved = float(p.bp_pol[grid_index,ten,0,j,zz,0,0])
            if not c['lo']-1e-9 <= saved <= c['hi']+1e-9:
                raise RuntimeError('Stored saving lies outside exact branch domain')
            record.update(saved_saving=saved,saved_action_value=objective(saved,c),
                optimum_minus_saved=optimum['value']-objective(saved,c),
                saved_consumption=float(p.c_pol[grid_index,ten,0,j,zz,0,0]))
            co,ho,_,_=allocation(saved,c)
            record['reconstructed_consumption']=cb+co
            record['saved_consumption_gap']=record['saved_consumption']-cb-co
            if not ten:
                record['saved_housing_gap']=float(p.hR_pol[grid_index,ten,0,j,zz,0,0])-hb-ho
        records.append(record);lookups[key]=record;contexts[key]=c
        return record
    nodes=[]
    for ib in (52,53):
        wealth=float(grid[ib]); assert abs(wealth-contract['wealth_values'][ib-52])<1e-12
        renter=audit(0,wealth,ib)
        saved_values=[renter['saved_action_value']]; best_values=[renter['optimum']['value']]
        products=[]
        for ten,h in enumerate(P.H_own,1):
            phi=float(SD.phi_choice[0,ten,0,0]);post=wealth-price*float(h)
            dp=(1-phi)*price*float(h);floor=-phi*price*float(h)
            if wealth<dp or post<floor:
                saved_values.append(-1e10);best_values.append(-1e10);continue
            ix,w=model.interp_indices(grid,np.array([post])); ix=int(ix[0]);w=float(w[0])
            lower=audit(ten,float(grid[ix]),ix);upper=audit(ten,float(grid[ix+1]),ix+1)
            sv=(1-w)*lower['saved_action_value']+w*upper['saved_action_value']
            bv=(1-w)*lower['optimum']['value']+w*upper['optimum']['value']
            exact=audit(ten,post)
            saved_values.append(sv);best_values.append(bv)
            products.append(dict(tenure=ten,rooms=float(h),post_purchase_wealth=post,
                downpayment=dp,collateral_floor=floor,lower_index=ix,upper_index=ix+1,
                upper_weight=w,interpolated_saved_action_value=sv,interpolated_optimal_value=bv,
                direct_offgrid_optimal_value=exact['optimum']['value'],
                direct_minus_interpolated=exact['optimum']['value']-bv))
        kappa=float(P.tenure_choice_kappa)
        predicted=np.exp((np.asarray(saved_values)-max(saved_values))/kappa);predicted/=predicted.sum()
        observed=np.asarray(p.tenure_probs[ib,0,0,j,zz,0,0,:],dtype=float)
        loggap=kappa*math.log(float(observed[3]/observed[0]))-(saved_values[3]-saved_values[0])
        nodes.append(dict(wealth_index=ib,wealth=wealth,products=products,
            saved_branch_values=saved_values,optimal_branch_values=best_values,
            observed_probabilities=observed.tolist(),reconstructed_probabilities=predicted.tolist(),
            maximum_probability_gap=float(np.max(np.abs(observed-predicted))),
            sixroom_minus_renter_logodds_value_gap=loggap,
            post_fertility_mass=float(e.g_post_fertility[ib,0,0,j,zz,0,0])))
    crosses=[]
    for a in records:
        for b in records:
            if a is b or a['tenure']!=b['tenure'] or 'saved_saving' not in a or 'saved_saving' not in b:
                continue
            c=contexts[(a['tenure'],float(a['current_branch_wealth']))];action=b['saved_saving']
            if c['lo']<=action<=c['hi']:
                crosses.append(dict(target_tenure=a['tenure'],target_grid_index=a['grid_index'],
                    source_grid_index=b['grid_index'],saving=action,value=objective(action,c),
                    gain_over_target_saved=objective(action,c)-a['saved_action_value']))
    result=dict(status='completed_fixed_continuation_diagnostic',certified_global_solution=False,
        full_Bellman_solves=0,GE_solves=0,scalar_branch_optimizations=count,
        survival=survival,pension=float(P.pension),psi=float(P.psi_child),price=price,
        utility=dict(cb=cb,hb=hb,pc=pc,alpha=alpha,equivalence_scale=es,chi=float(P.chi),beta=float(P.beta)),
        branch_records=records,tenure_nodes=nodes,crossed_feasible_actions=crosses,
        maximum_saved_action_objective_gain=max(r.get('optimum_minus_saved',0) for r in records),
        maximum_logodds_value_gap=max(abs(n['sixroom_minus_renter_logodds_value_gap']) for n in nodes),
        maximum_probability_gap=max(n['maximum_probability_gap'] for n in nodes),
        checkpoint_sha256=contract['checkpoint_sha256'],seconds=time.monotonic()-started,
        interpretation='Direct off-grid optima are supplemental interpolation diagnostics, not the production choice values. Inspect gains and reproduction gaps before certifying either economics or numerics.')
    (args.output/'audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['status','scalar_branch_optimizations','maximum_saved_action_objective_gain','maximum_logodds_value_gap','maximum_probability_gap','seconds']}))


if __name__=='__main__':
    main()
