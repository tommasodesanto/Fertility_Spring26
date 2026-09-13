#!/usr/bin/env python3
"""Bounded A0/A+ successive surprises with an endogenous finite boundary.

Both migration comparisons use the same finite constant-conditions truncation.
No open-population stationary endpoint is imported or claimed. A successful
finite root is provisional until separate fixed-shock horizon comparisons.
"""
from __future__ import annotations
import argparse
import copy
import csv
from dataclasses import replace
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import resource
import sys
import threading
import time
from types import SimpleNamespace as NS
from unittest.mock import patch

for _thread_key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_thread_key] = '1'
import numpy as np


class _AutomaticFiscalPolish(BaseException):
    def __init__(self, record):self.record=record


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def clean(value):
    if isinstance(value, np.ndarray): return clean(value.tolist())
    if isinstance(value, np.generic): return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value): return None
    if isinstance(value, dict): return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)): return [clean(v) for v in value]
    if isinstance(value, Path): return str(value)
    return value


def save(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(clean(value), indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def checkpoint(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    with gzip.open(temporary, 'wb', compresslevel=1) as f:
        pickle.dump(value, f, protocol=5)
    temporary.replace(path)
    return dict(path=str(path), sha256=sha(path))


def verify_pins(pins):
    if not isinstance(pins, dict) or not pins:
        raise ValueError('Nonempty explicit source/input SHA256 mapping required')
    for path, digest in pins.items():
        if sha(path) != digest:
            raise ValueError('Pinned source/input changed: '+str(path))


def migration_case(demographics, case):
    if case == 'A+': return demographics
    if case != 'A0': raise ValueError('This driver supports A0/A+ only')
    zero = {year: np.zeros_like(values, dtype=float)
            for year, values in demographics.net_migration.items()}
    result = replace(demographics, net_migration=zero)
    if any(np.count_nonzero(values) for values in result.net_migration.values()):
        raise RuntimeError('Zero migration must hold in every age-sex cell')
    return result


def unpack_coordinates(values, count):
    x = np.asarray(values, dtype=float)
    if x.shape != (3*(count+1),) or not np.isfinite(x).all() or np.any(x <= 0):
        raise ValueError('Positive complete N+1 price/pension/transfer blocks required')
    return tuple(x.reshape(3, count+1))


def initial_seed_step(plan, manifest):
    if 'seed_step' in manifest:return float(manifest['seed_step'])
    seeds=plan.get('seed_steps',[-.01])
    return float(next(iter(seeds.values())) if isinstance(seeds,dict) else seeds[0])


def initial_coordinate_seed(manifest,case,count,default):
    """Load an optional pinned numerical guess without importing solved state."""
    if case not in ('A0','A+') or count not in (6,24,100):
        raise ValueError('Initial coordinate seed requires a supported history track')
    default=np.asarray(default,dtype=float);unpack_coordinates(default,count)
    profiles=manifest.get('initial_coordinate_seeds')
    key=('A0' if case=='A0' else 'Aplus')+f'_{count}'
    allowed={f'{name}_{n}' for name in ('A0','Aplus') for n in (6,24,100)}
    base=dict(mode='default',case=case,count=count,selection='default',seed_origin=None,
        default_coordinates=default.tolist(),selected_coordinates=default.tolist(),
        coordinates_changed=False,
        driver=dict(path=str(Path(__file__).resolve()),sha256=sha(Path(__file__).resolve())))
    if profiles is None:return default.copy(),base
    if not isinstance(profiles,dict) or any(name not in allowed for name in profiles):
        raise ValueError('Initial coordinate seeds require recognized A0/Aplus count keys')
    if key not in profiles:return default.copy(),base
    if 'resume_history' in manifest:
        raise ValueError('Initial coordinate seeds cannot be combined with history resume')
    entry=profiles[key]
    if not isinstance(entry,dict) or set(entry)!= {'path','sha256'}:
        raise ValueError('Initial coordinate seed entry must contain only path and sha256')
    seed_path=Path(entry['path'])
    if not seed_path.is_absolute() or str(seed_path) not in manifest.get('file_sha256',{}):
        raise ValueError('Initial coordinate seed JSON must be an absolute manifest-pinned input')
    if manifest['file_sha256'][str(seed_path)]!=entry['sha256'] or sha(seed_path)!=entry['sha256']:
        raise ValueError('Initial coordinate seed JSON hash mismatch')
    seed=json.loads(seed_path.read_text())
    fields={'case','count','start_year','coordinates','label','source_root_receipt','selection'}
    extension_fields={'source_count','rule'}
    if not isinstance(seed,dict) or not fields.issubset(seed) or seed['label']!='numerical_guess_only':
        raise ValueError('Initial coordinate seed JSON schema mismatch')
    extension=set(seed)-fields
    if extension and extension != extension_fields:
        raise ValueError('Initial coordinate seed JSON schema mismatch')
    extended=bool(extension)
    if extended and (seed['source_count']!=24 or seed['rule']!='hold_last' or count!=100):
        raise ValueError('Initial coordinate seed extension requires source_count 24, hold_last, and target count 100')
    if seed['case']!=case or seed['count']!=count or seed['start_year']!=2007:
        raise ValueError('Initial coordinate seed track mismatch')
    origin=seed['source_root_receipt']
    if not isinstance(origin,dict) or set(origin)!= {'path','sha256'}:
        raise ValueError('Seed source root receipt requires path and sha256')
    root_path=Path(origin['path'])
    if not root_path.is_absolute() or sha(root_path)!=origin['sha256']:
        raise ValueError('Seed source root receipt hash mismatch')
    root=json.loads(root_path.read_text())
    source_count=seed['source_count'] if extended else count
    if root.get('case')!=case or root.get('count')!=source_count or root.get('start_year')!=2007:
        raise ValueError('Seed source root receipt track mismatch')
    selection=seed['selection']
    if selection not in ('final','best'):
        raise ValueError('Seed price selection must be final or best')
    selected=root.get(selection)
    if not isinstance(selected,dict) or selected.get('mapping_valid') is not True:
        raise ValueError('Selected seed root mapping is not valid')
    prices=np.asarray(selected.get('prices'),dtype=float);source_prices=unpack_coordinates(prices,source_count)
    if extended:
        prices=np.concatenate([np.pad(block,(0,count-source_count),'edge') for block in source_prices])
    coordinates=np.asarray(seed['coordinates'],dtype=float);unpack_coordinates(coordinates,count)
    if not np.array_equal(coordinates,prices):
        raise ValueError('Seed coordinates do not exactly match selected root prices')
    receipt=dict(base,mode='pinned_numerical_guess',selection=selection,
        seed_origin=dict(path=str(seed_path),sha256=entry['sha256'],
            source_root_receipt=dict(path=str(root_path),sha256=origin['sha256'])),
        selected_coordinates=coordinates.tolist(),coordinates_changed=not np.array_equal(default,coordinates))
    if extended:
        receipt['seed_origin'].update(source_count=source_count,rule='hold_last')
    return coordinates.copy(),receipt


def shift_forecast_coordinates(values,count):
    """Advance the numerical starting guess by one realized four-year date."""
    unpack_coordinates(values,count)
    blocks=np.asarray(values,dtype=float).reshape(3,count+1)
    return np.column_stack((blocks[:,1:],blocks[:,-1])).ravel()


def cached_boundary(template, actual_g, old, audit, runtime):
    """Re-account the same lifetime policy on actual carried terminal heads.

    The lifetime Bellman value depends on boundary coordinates, not aggregate
    population. This avoids a second household solve without substituting an
    inherited or normalized population into the boundary equilibrium equations.
    """
    _, _, primitive, _, _ = runtime
    import e5f_balanced_terminal as balanced
    from e5f_social_security import fiscal_accounts
    import e5f_closed_finite_boundary as closed
    P, policy, grid = template.parameters, template.policy, template.b_grid
    calendar = primitive.calendar
    shared = primitive.model.precompute_shared(P, grid)
    actual = calendar.evaluate_period(policy.price, actual_g.copy(), P, grid,
        shared, calendar.SolveCounter(), supply_rule=old.supply_rule,
        supplied_policy=policy)
    diagnostics, gates = balanced._household_checks(actual, P, shared, grid,
        float(P.user_cost_rate)*float(policy.price[0]), primitive, audit)
    heads = float(actual.g_current.sum())
    supply, demand = float(actual.supply_by_loc[0]), float(actual.demand_by_loc[0])
    accounts = fiscal_accounts(actual.g_current, P)
    revenue = float(primitive.model.property_tax_revenue_from_distribution(
        actual.g_current, actual.policy.hR_pol, actual.policy.price, P))
    outlays = float(P.property_tax_lump_sum_transfer)*heads
    if not np.isfinite([heads,supply,demand,revenue,outlays]).all() or heads <= 0 or supply <= 0:
        raise ValueError('Invalid actual boundary mass or market/fiscal accounting')
    gates['actual_head_mass_preserved'] = abs(heads-float(actual_g.sum())) <= 2e-9*max(1.,float(actual_g.sum()))
    accounts.update(household_heads=heads, property_tax_revenue=revenue,
        equal_transfer_outlays=outlays, property_tax_budget_residual=revenue-outlays)
    residuals = dict(housing_relative=(demand-supply)/supply,
        pension_relative=accounts['scaled_pension_budget_residual'],
        rebate_relative=closed._scaled_difference(revenue,outlays))
    return closed.FiniteBoundaryEvaluation(P,grid,actual_g.copy(),actual.policy,
        residuals,accounts,dict(diagnostics, future_fiscal_consistency_verified=False,
        horizon_status='unverified_finite_truncation'),gates)


def validated_fixed_asset_prices(values, count, price_bounds, pf, parameters):
    prices=np.asarray(values,dtype=float);lo,hi=price_bounds
    if (prices.shape!=(count+1,) or not np.isfinite(prices).all()
            or np.any(prices<=0) or np.any(prices<lo) or np.any(prices>hi)):
        raise ValueError('Fixed asset prices must be a positive finite in-bound N+1 path')
    rents=np.asarray(pf.rents_from_asset_prices(prices[:-1],float(prices[-1]),parameters),dtype=float)
    if rents.shape!=(count,) or not np.isfinite(rents).all() or np.any(rents<=0):
        raise ValueError('Fixed asset prices must imply positive finite PF rents')
    return prices.copy()


def fiscal_polish_switch_prices(record, width, market_tolerance, maximum_evaluations):
    """Return the price block only after an eligible coupled iterate mapping."""
    if (record.get('phase')!='iterate' or record.get('mapping_valid') is not True
            or not isinstance(record.get('evaluation'),int)
            or maximum_evaluations-record['evaluation']<2):
        return None
    residual=np.asarray(record.get('residual'),dtype=float)
    coordinates=np.asarray(record.get('prices'),dtype=float)
    if (residual.shape!=(3*width,) or coordinates.shape!=(3*width,)
            or not np.isfinite(residual).all() or not np.isfinite(coordinates).all()
            or np.any(np.abs(residual[:width])>=market_tolerance)
            or np.max(np.abs(residual[width:]))<market_tolerance):
        return None
    return coordinates[:width].copy()


def solve_forecast(*, inherited, old, demographics, psi, count, initial,
                   controls, audit, deadline, folder, case, initial_jacobian=None,
                   fixed_asset_prices=None):
    import e5f_rebated_surprises as rebated
    import e5f_closed_finite_boundary as closed
    from e5f_balanced_terminal import _household_checks
    from e5f_social_security import fiscal_accounts
    from e5f_matched_pf_path_root import solve_price_path
    import run_e5f_transition_calibration as fertility
    runtime = rebated._runtime()
    _, joined, primitive, _, rent_domain = runtime
    folder = Path(folder); folder.mkdir(parents=True, exist_ok=True)
    width=count+1; rc=dict(controls)
    bounds=[tuple(rc.pop(name)) for name in ('price_bounds','pension_bounds','transfer_bounds')]
    for lo,hi in bounds:
        if not np.isfinite([lo,hi]).all() or not 0 < lo < hi:
            raise ValueError('Joint log-root needs explicit positive bounds')
    fixed_prices=None
    if fixed_asset_prices is not None:
        fixed_prices=validated_fixed_asset_prices(
            fixed_asset_prices,count,bounds[0],joined.pf,old.parameters)
    automatic_fiscal_polish=rc.pop('automatic_fiscal_polish',False)
    if type(automatic_fiscal_polish) is not bool:
        raise ValueError('automatic_fiscal_polish must be Boolean')
    automatic_fiscal_polish=bool(automatic_fiscal_polish and fixed_prices is None)
    if not 2 <= rc['max_evaluations'] <= 24 or not 0 < rc['market_tolerance'] <= 2e-4:
        raise ValueError('Retained mapping budget/market gate required')
    if not 0 <= rc['final_reproduction_tolerance'] <= 2e-10:
        raise ValueError('Retained exact replay gate required')
    rc.pop('fiscal_tolerance',None); rc.pop('fiscal_slope',None)
    rc.setdefault('slope',rc.pop('market_slope',1.))
    allowed={'slope','market_tolerance','max_log_step','damping','max_evaluations',
        'max_condition_number','worsening_factor','final_reproduction_tolerance'}
    rc={k:v for k,v in rc.items() if k in allowed}
    latest={}; mapping=0
    def progress(record):
        save(folder/'latest_completed.json',record)
        if record.get('new_best'): save(folder/'best_so_far.json',record)
    def project(raw):
        x=np.asarray(raw,dtype=float).reshape(3,width).copy()
        for j,(lo,hi) in enumerate(bounds): x[j]=np.clip(x[j],lo,hi)
        if fixed_prices is None:
            endpoint=NS(parameters=old.parameters,asset_price=float(x[0,-1]))
            x[0,:-1]=rent_domain.project_price_path_to_positive_rents(
                x[0,:-1],terminal=endpoint,minimum_rent_share=1e-6)[0]
            if np.any(x[0]>bounds[0][1]): raise ValueError('Positive-rent projection exceeds price bound')
        else:
            x[0]=fixed_prices
        return x.ravel()
    def evaluate(raw):
        nonlocal mapping,latest
        mapping+=1; began=time.monotonic(); latest={}
        prices,pensions,transfers=unpack_coordinates(raw,count)
        Q=copy.deepcopy(old.parameters); Q.psi_child=float(psi)
        # Only the policy from this seed call is used. Actual boundary budgets
        # below use the state carried by THIS complete trial path.
        terminal=closed.boundary_evaluation(parameters=Q,g_pre=inherited.households.g_pre,
            grid=old.b_grid,supply_rule=old.supply_rule,price=float(prices[-1]),
            pension=float(pensions[-1]),transfer=float(transfers[-1]),audit_controls=audit,
            deadline_monotonic=deadline,callback=lambda row:save(folder/'latest_phase.json',dict(row,mapping=mapping)))
        rows=[];audits=[];observations=[];snapshot={};snapshot2023={}
        rents=joined.pf.rents_from_asset_prices(prices[:-1],float(prices[-1]),old.parameters)
        def observe(i,e,P,grid,shared):
            if time.monotonic()>=deadline: raise TimeoutError('Forecast deadline reached during dated sweep')
            if i != len(rows): raise RuntimeError('Nonsequential dated observer')
            if P.pension != pensions[i] or P.property_tax_lump_sum_transfer != transfers[i] or P.tau_pay != .179:
                raise RuntimeError('Dated fiscal values differ from root coordinates')
            d,gates=_household_checks(e,P,shared,grid,float(rents[i]),primitive,audit)
            if not all(gates.values()): raise RuntimeError('Dated household audit failed')
            a=fiscal_accounts(e.g_current,P)
            revenue=float(primitive.model.property_tax_revenue_from_distribution(e.g_current,e.policy.hR_pol,e.policy.price,P))
            tax=rebated.rebated_tax_accounts(property_tax_revenue=revenue,
                transfer_per_head=transfers[i],head_mass=float(e.g_current.sum()))
            rows.append(rebated.dated_residual(demand=e.demand_by_loc[0],supply=e.supply_by_loc[0],payroll_accounts=a,tax_accounts=tax))
            audits.append(dict(year=inherited.year+4*i,diagnostics=d,gates=gates,payroll=a,rebate=tax))
            observations.append(dict(calendar_year=inherited.year+4*i,**fertility.period_fertility_diagnostics(e,P)))
            pack=dict(parameters=P,b_grid=grid,evaluation=e,shared=shared,supply_rule=old.supply_rule)
            if i==0:snapshot.update(pack)
            if inherited.year+4*i==2023:snapshot2023.update(pack)
            save(folder/'latest_date.json',dict(mapping=mapping,year=inherited.year+4*i,elapsed_seconds=time.monotonic()-began))
        solves=0;native=joined.pf.solve_date_policy
        def dated(*args,**kwargs):
            nonlocal solves
            if time.monotonic()>=deadline:raise TimeoutError('Forecast deadline reached during backward sweep')
            solves+=1
            save(folder/'latest_backward.json',dict(mapping=mapping,household_date_solve=solves))
            return native(*args,**kwargs)
        with patch.object(joined.pf,'solve_date_policy',dated):
            path=rebated.evaluate_forecast(inherited=inherited,old_state=old,
                demographics=demographics,prices=prices[:-1],pensions=pensions[:-1],
                transfers=transfers[:-1],psi=psi,terminal=terminal,observer=observe)
        if len(rows)!=count:raise RuntimeError('Missing dated market/fiscal observation')
        if case=='A0':
            for row in path.person_tail.rows:
                if row['annual_net_migration_over_period'] != 0. or row['net_migrant_heads_over_period'] != 0.:
                    raise RuntimeError('Nonzero migration in closed forecast')
        actual=path.person_tail.terminal_state.g_pre
        boundary=cached_boundary(terminal,actual,old,audit,runtime)
        rows.append(np.array([boundary.residuals['housing_relative'],
            200.*boundary.residuals['pension_relative'],200.*boundary.residuals['rebate_relative']]))
        latest=dict(path=path,coordinates=np.asarray(raw).copy(),mapping=mapping,
            boundary=boundary,observations=observations,snapshot=snapshot,snapshot2023=snapshot2023)
        save(folder/'mapping_receipt.json',dict(mapping=mapping,elapsed_seconds=time.monotonic()-began,
            explicit_dates=count,boundary_dates=1,bellman_date_solves=solves,
            peak_rss_native_units=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            peak_rss_unit='bytes' if sys.platform=='darwin' else 'KiB',
            boundary_accounts=boundary.actual_accounts,boundary_residuals=boundary.residuals,
            boundary_household_gates=boundary.gates,finite_horizon_only=True))
        return dict(residual=rebated.stack_dated_residuals(rows),mapping_valid=boundary.mapping_valid,
            payload=dict(mapping=mapping,dated_audits=audits,boundary_accounts=boundary.actual_accounts,
                boundary_residuals=boundary.residuals,boundary_household_gates=boundary.gates))
    default_jacobian=np.diag(np.r_[np.full(width,-float(rc['slope'])),np.full(2*width,-200.)])
    root_started=time.monotonic()
    root_arguments=dict(initial_prices=initial,evaluate=evaluate,project=project,
        deadline_monotonic=deadline,default_jacobian=default_jacobian,**rc)
    if not automatic_fiscal_polish:
        receipt=solve_price_path(callback=progress,initial_jacobian=initial_jacobian,**root_arguments)
    else:
        coupled_events=[]
        def coupled_progress(record):
            enriched=dict(record,root_phase='coupled',
                total_evaluation=record.get('evaluation'),forecast_elapsed_seconds=time.monotonic()-root_started)
            coupled_events.append(copy.deepcopy(enriched));progress(enriched)
            prices=fiscal_polish_switch_prices(record,width,rc['market_tolerance'],rc['max_evaluations'])
            if prices is not None:raise _AutomaticFiscalPolish(record)
        try:
            receipt=solve_price_path(callback=coupled_progress,initial_jacobian=initial_jacobian,**root_arguments)
        except _AutomaticFiscalPolish as switch:
            switch_evaluation=int(switch.record['evaluation'])
            remaining=rc['max_evaluations']-switch_evaluation
            polish_budget=min(6,remaining)
            fixed_prices=validated_fixed_asset_prices(
                np.asarray(switch.record['prices'],dtype=float)[:width],count,bounds[0],joined.pf,old.parameters)
            polish_events=[];polish_started_offset=time.monotonic()-root_started
            def polish_progress(record):
                local=record.get('evaluation')
                enriched=dict(record,root_phase='fiscal_polish',
                    total_evaluation=None if local is None else switch_evaluation+local,
                    forecast_elapsed_seconds=time.monotonic()-root_started)
                polish_events.append(copy.deepcopy(enriched));progress(enriched)
            polish_controls=dict(rc,max_evaluations=polish_budget)
            receipt=solve_price_path(initial_prices=np.asarray(switch.record['prices'],dtype=float),
                evaluate=evaluate,project=project,deadline_monotonic=deadline,callback=polish_progress,
                initial_jacobian=None,default_jacobian=default_jacobian,**polish_controls)
            coupled_by_evaluation={};coupled_order=[]
            for row in coupled_events:
                evaluation=row.get('evaluation')
                if evaluation is not None:
                    if evaluation not in coupled_by_evaluation:coupled_order.append(evaluation)
                    coupled_by_evaluation[evaluation]=row
            coupled_history=[coupled_by_evaluation[evaluation] for evaluation in coupled_order]
            polish_history=[]
            for row in receipt.get('history',[]):
                value=copy.deepcopy(row);value['root_phase']='fiscal_polish'
                value['phase_evaluation']=value['evaluation'];value['evaluation']=switch_evaluation+value['evaluation']
                value['forecast_elapsed_seconds']=polish_started_offset+float(value.get('elapsed_seconds',0.))
                polish_history.append(value)
            local_evaluations=int(receipt.get('evaluations',0))
            receipt['history']=coupled_history+polish_history
            receipt['evaluations']=switch_evaluation+local_evaluations
            receipt['elapsed_seconds']=time.monotonic()-root_started
            receipt['root_phase_ledger']=[dict(root_phase=row.get('root_phase'),
                evaluation=row.get('evaluation'),phase=row.get('phase'),score=row.get('score'),
                mapping_valid=row.get('mapping_valid'),forecast_elapsed_seconds=row.get('forecast_elapsed_seconds'))
                for row in receipt['history']]
            receipt['automatic_fiscal_polish']=dict(switched=True,switch_evaluation=switch_evaluation,
                coupled_evaluations=switch_evaluation,fiscal_polish_evaluations=local_evaluations,
                total_actual_mappings=mapping,original_maximum_evaluations=rc['max_evaluations'],
                fiscal_polish_maximum_evaluations=polish_budget,fixed_asset_prices=fixed_prices)
            if receipt['evaluations']>rc['max_evaluations'] or mapping!=receipt['evaluations']:
                raise RuntimeError('Automatic fiscal polish mapping accounting violated the original budget')
        else:
            if mapping!=receipt.get('evaluations'):
                raise RuntimeError('Automatic fiscal-polish coupled mapping accounting differs from root receipt')
            history=[]
            for row in receipt.get('history',[]):
                value=copy.deepcopy(row);value['root_phase']='coupled'
                value['total_evaluation']=value['evaluation'];history.append(value)
            receipt['history']=history
            receipt['root_phase_ledger']=[dict(root_phase='coupled',evaluation=row['evaluation'],
                phase=row.get('phase'),score=row.get('score'),mapping_valid=row.get('mapping_valid'),
                forecast_elapsed_seconds=row.get('elapsed_seconds')) for row in history]
            receipt['automatic_fiscal_polish']=dict(switched=False,
                coupled_evaluations=receipt.get('evaluations',0),fiscal_polish_evaluations=0,
                total_actual_mappings=mapping,
                original_maximum_evaluations=rc['max_evaluations'])
    final=receipt.get('final')
    finite=bool(receipt.get('converged') and final is not None and latest
        and final['payload']['mapping']==latest['mapping']
        and np.array_equal(final['prices'],latest['coordinates']))
    receipt.update(finite_horizon_market_fiscal_converged=finite,start_year=inherited.year,
        psi=float(psi),case=case,count=count,boundary_status='unverified_finite_truncation',
        terminal_distance_passed=False,horizon_verified=False,production_eligible=False)
    if fixed_prices is not None:receipt['fixed_asset_prices']=fixed_prices
    save(folder/'root_receipt.json',receipt)
    if latest:
        save(folder/'rows.json',latest['path'].rows);save(folder/'fertility.json',latest['observations'])
    if not finite:return NS(path=latest.get('path'),root_receipt=receipt,next_state=None,realized_row=None),latest
    p,b,t=unpack_coordinates(final['prices'],count)
    next_state=rebated.first_period_state(inherited=inherited,old_state=old,demographics=demographics,
        path=latest['path'],prices=p[:-1],pensions=b[:-1],transfers=t[:-1],psi=psi)
    result=NS(path=latest['path'],root_receipt=receipt,next_state=next_state,realized_row=latest['path'].rows[0])
    checkpoint(folder/'accepted_forecast.pkl.gz',dict(result=result,boundary=latest['boundary'],coordinates=final['prices']))
    checkpoint(folder/'first_period_diagnostics.pkl.gz',latest['snapshot'])
    if latest['snapshot2023']:
        checkpoint(folder/'native_2023_snapshot.pkl.gz',latest['snapshot2023'])
    return result,latest


def reusable_forecast_jacobian(result, count, enabled):
    """Copy an accepted root Jacobian for the immediately following psi trial."""
    if not enabled or result.next_state is None:
        return None
    receipt=result.root_receipt
    if (not receipt.get('finite_horizon_market_fiscal_converged')
            or receipt.get('fixed_asset_prices') is not None):
        return None
    width=3*(count+1)
    matrix=np.asarray(receipt.get('final_jacobian'),dtype=float)
    if matrix.shape!=(width,width) or not np.isfinite(matrix).all():
        return None
    return matrix.copy()


def _artifact(entry, label):
    if not isinstance(entry,dict) or set(entry)!= {'path','sha256'}:
        raise ValueError(label+' requires exactly path and sha256')
    path=Path(entry['path']).resolve()
    if sha(path)!=entry['sha256']:raise ValueError('Pinned resume artifact changed: '+label)
    return path


def _economic_manifest(value):
    """Remove only resume, numerical budget, and relocated history-driver metadata."""
    result=copy.deepcopy(value)
    for key in ('resume_history','reuse_forecast_jacobian','policy_reserve_seconds','forecast_seconds'):
        result.pop(key,None)
    pins=result.get('file_sha256',{})
    result['file_sha256']={path:digest for path,digest in pins.items()
        if Path(path).name!='run_e5f_final_rebated_history.py'}
    return result


def _same_state(left, right, seen=None):
    """Exact recursive comparison of every serialized inherited-state field."""
    if seen is None:seen=set()
    pair=(id(left),id(right))
    if pair in seen:return True
    seen.add(pair)
    if type(left) is not type(right):return False
    if isinstance(left,np.ndarray):
        return left.dtype==right.dtype and left.shape==right.shape and np.array_equal(left,right,equal_nan=True)
    if isinstance(left,np.generic):return left.dtype==right.dtype and left.tobytes()==right.tobytes()
    if isinstance(left,dict):
        return left.keys()==right.keys() and all(_same_state(left[k],right[k],seen) for k in left)
    if isinstance(left,(tuple,list)):
        return len(left)==len(right) and all(_same_state(a,b,seen) for a,b in zip(left,right))
    if hasattr(left,'__dict__'):
        return _same_state(vars(left),vars(right),seen)
    if isinstance(left,float):
        return left==right or (np.isnan(left) and np.isnan(right))
    return left==right


def load_resume_history(spec, *, current_manifest, targets, tolerance, case, count,
                        initial_checkpoint_sha256, source_root):
    """Validate and load one explicitly pinned contiguous accepted history prefix."""
    required={'source_manifest','source_contract','realized_fit','last_realized_state','windows'}
    if not isinstance(spec,dict) or set(spec)!=required:
        raise ValueError('resume_history has an incomplete explicit artifact contract')
    source_manifest_path=_artifact(spec['source_manifest'],'source manifest')
    source_manifest=json.loads(source_manifest_path.read_text())
    verify_pins(source_manifest['file_sha256'])
    if _economic_manifest(source_manifest)!=_economic_manifest(current_manifest):
        raise ValueError('Resume and current manifests differ outside driver, reuse, or budget metadata')
    source_contract_path=_artifact(spec['source_contract'],'source contract')
    source_contract=json.loads(source_contract_path.read_text())
    if (source_contract.get('manifest_sha256')!=spec['source_manifest']['sha256']
            or source_contract.get('case')!=case or source_contract.get('count')!=count
            or source_contract.get('initial_checkpoint_sha256')!=initial_checkpoint_sha256
            or Path(source_contract.get('source_root','')).resolve()!=Path(source_root).resolve()):
        raise ValueError('Resume source contract differs in case, count, initial state, or model source')
    realized_path=_artifact(spec['realized_fit'],'realized fit')
    realized=json.loads(realized_path.read_text())
    windows=spec['windows'];years=[int(row['decision_year']) for row in targets]
    if (not isinstance(windows,list) or not windows or len(windows)>len(years)
            or [entry.get('year') for entry in windows]!=years[:len(windows)]
            or len(realized)!=len(windows)):
        raise ValueError('Resume windows are not a nonempty contiguous historical prefix')
    last_result=None;last_root=None
    for index,(entry,fit,target) in enumerate(zip(windows,realized,targets)):
        if set(entry)!= {'year','fit','root_receipt','accepted_forecast'}:
            raise ValueError('Each resumed window needs exact fit/root/forecast pins')
        fit_path=_artifact(entry['fit'],f'window {entry["year"]} fit')
        root_path=_artifact(entry['root_receipt'],f'window {entry["year"]} root')
        accepted_path=_artifact(entry['accepted_forecast'],f'window {entry["year"]} accepted forecast')
        if (root_path.parent!=accepted_path.parent
                or root_path.parent not in (fit_path.parent,fit_path.parent/'alternative')):
            raise ValueError('Resumed root and accepted forecast do not belong to the selected trial')
        saved_fit=json.loads(fit_path.read_text());root=json.loads(root_path.read_text())
        if saved_fit!=fit or Path(fit['folder']).resolve()!=fit_path.parent:
            raise ValueError('Resumed selected fit differs from the realized-fit ledger')
        year=years[index];desired=float(target['period_tfr_arithmetic_mean'])
        if (fit.get('year')!=year or float(fit['target'])!=desired
                or not np.isfinite([float(fit['psi']),float(fit['model']),float(fit['gap'])]).all()
                or float(fit['gap'])!=float(fit['model'])-desired
                or abs(float(fit['gap']))>tolerance):
            raise ValueError('Resumed fit fails its pinned empirical target or tolerance')
        final=root.get('final');reproduction=float(root.get('final_reproduction_max_abs',np.nan))
        if (root.get('converged') is not True or root.get('status')!='converged'
                or root.get('finite_horizon_market_fiscal_converged') is not True
                or root.get('start_year')!=year or root.get('case')!=case or root.get('count')!=count
                or float(root.get('psi',np.nan))!=float(fit['psi']) or not isinstance(final,dict)
                or final.get('mapping_valid') is not True or not np.isfinite(reproduction)
                or reproduction>min(2e-10,float(tolerance))):
            raise ValueError('Resumed selected root lacks its finite exact acceptance certificate')
        with gzip.open(accepted_path,'rb') as stream:accepted=pickle.load(stream)
        result=accepted.get('result');coordinates=np.asarray(accepted.get('coordinates'),dtype=float)
        final_prices=np.asarray(final.get('prices'),dtype=float);width=3*(count+1)
        if (result is None or result.next_state is None or result.next_state.year!=year+4
                or clean(result.root_receipt)!=root or coordinates.shape!=(width,)
                or final_prices.shape!=(width,) or not np.array_equal(coordinates,final_prices)):
            raise ValueError('Accepted forecast packet differs from its selected root or next year')
        last_result=result;last_root=root
    state_path=_artifact(spec['last_realized_state'],'last realized state')
    with gzip.open(state_path,'rb') as stream:inherited=pickle.load(stream)
    if inherited.year!=years[len(windows)-1]+4 or not _same_state(inherited,last_result.next_state):
        raise ValueError('Last accepted forecast and resumed inherited state differ')
    return dict(realized=realized,inherited=inherited,
        initial=shift_forecast_coordinates(last_root['final']['prices'],count),
        completed_years=years[:len(windows)],source_manifest=str(source_manifest_path),
        source_contract=str(source_contract_path),last_realized_state=str(state_path))


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',required=True);parser.add_argument('--case',choices=['A0','A+'],required=True)
    parser.add_argument('--count',type=int,choices=[6,24,100],required=True)
    parser.add_argument('--output',required=True);parser.add_argument('--seconds',type=float,required=True)
    args=parser.parse_args(argv)
    if not np.isfinite(args.seconds) or not 0<args.seconds<=12*3600:parser.error('--seconds must be in (0,43200]')
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    started=time.monotonic();deadline=started+args.seconds;realized=[]
    stop=threading.Event()
    def heartbeat():
        while not stop.wait(min(60.,args.seconds)):
            save(out/'heartbeat.json',dict(elapsed_seconds=time.monotonic()-started,
                case=args.case,count=args.count,completed_windows=len(realized)))
            if time.monotonic()>=deadline:
                save(out/'failure.json',dict(error='Hard total controller deadline exhausted',realized=realized))
                os._exit(124)
    threading.Thread(target=heartbeat,daemon=True).start()
    try:
        manifest=json.loads(Path(args.manifest).read_text());verify_pins(manifest['file_sha256'])
        if manifest['prior_plan'] not in manifest['file_sha256']:
            raise ValueError('The complete prior plan must itself be SHA256 pinned')
        plan=json.loads(Path(manifest['prior_plan']).read_text());verify_pins(plan['file_sha256'])
        if plan['empirical_blocks'] not in {**plan['file_sha256'],**manifest['file_sha256']}:
            raise ValueError('The empirical fertility-window input must be SHA256 pinned')
        root=Path(plan['source_root']);sys.path[:0]=[str(Path(__file__).resolve().parent),str(root/'code/model/tools'),str(root/'code/model')]
        import e5f_rebated_surprises as rebated
        from e5f_rebated_initial_bridge import build_rebated_initial_state
        from e5f_balanced_terminal import TerminalAuditControls
        from run_e5f_successive_surprises_overnight import next_psi,standard_graphs
        rebated._runtime()  # Configure the pinned sequential package before unpickling.
        summary=json.loads(Path(manifest['initial_summary']).read_text())
        if summary['status']!='verified_rebated_initial_smoke':raise ValueError('Initial rebated prerequisite has not passed')
        initial_root=Path(summary['source_root']).resolve()
        if initial_root!=root.resolve():
            if str(initial_root)!=manifest.get('initial_source_root'):
                raise ValueError('Unpinned initial scientific source root')
            equivalence=Path(manifest['kernel_equivalence'])
            if str(equivalence) not in manifest['file_sha256']:
                raise ValueError('Kernel equivalence receipt must be pinned')
            pairs=json.loads(equivalence.read_text())['pairs']
            if not pairs:raise ValueError('Missing initial/history kernel comparison')
            for pair in pairs:
                if sha(pair['initial'])!=pair['sha256'] or sha(pair['history'])!=pair['sha256']:
                    raise ValueError('Initial/history kernel equivalence failed')
        item=summary['checkpoint'];cp=item.get('path',item.get('checkpoint'))
        digest=item.get('sha256',item.get('checkpoint_sha256'))
        if not cp or not digest or sha(cp)!=digest:raise ValueError('Rebated initial checkpoint hash mismatch')
        with gzip.open(cp,'rb') as f:packet=pickle.load(f)
        initial_raw=json.loads(Path(manifest.get('initial_raw_summary',str(Path(cp).parent/'summary.json'))).read_text())
        old=build_rebated_initial_state(packet=packet,normalization=initial_raw['normalization'],
            outside_origin_entry_share=plan['outside_origin_entry_share'],preference_change_2023=0.)
        demographics=migration_case(packet['demographic_seed'],args.case)
        inherited=rebated.InheritedState(2007,old.initial_state)
        initial_psi=float(old.parameters.psi_child)
        audit=TerminalAuditControls(**plan['terminal_template']['audit_controls'])
        controls=dict(plan['history_root_controls']);controls.update(manifest.get('root_controls',{}))
        controls.setdefault('transfer_bounds',[1e-10,10.]);controls['max_evaluations']=min(24,int(controls['max_evaluations']))
        reuse_forecast_jacobian=manifest.get('reuse_forecast_jacobian',False)
        if type(reuse_forecast_jacobian) is not bool:
            raise ValueError('reuse_forecast_jacobian must be Boolean')
        targets=list(csv.DictReader(Path(plan['empirical_blocks']).open()))
        if [int(t['decision_year']) for t in targets]!=[2007,2011,2015,2019]:raise ValueError('Four pinned fertility windows required')
        tolerance=float(plan['fertility_fit_tolerance'])
        if not 0<tolerance<=.005:raise ValueError('Retained fertility fit tolerance required')
        save(out/'contract_receipt.json',dict(case=args.case,count=args.count,manifest_sha256=sha(args.manifest),
            initial_checkpoint_sha256=digest,source_root=str(root),migration_zero_by_cell=args.case=='A0',
            boundary='Actual-carried-state finite snapshot for both migration cases; replaces prior stationary-boundary comparison',
            reuse_forecast_jacobian=reuse_forecast_jacobian,
            resume_history_requested='resume_history' in manifest,
            horizon_verified=False,production_eligible=False,full_2023_table_status='native snapshot saved for authoritative observer replay'))
        width=args.count+1
        initial=np.r_[np.full(width,float(packet['evaluation'].policy.price[0])),
            np.full(width,float(old.parameters.pension)),np.full(width,float(old.parameters.property_tax_lump_sum_transfer))]
        initial,seed_receipt=initial_coordinate_seed(manifest,args.case,args.count,initial)
        save(out/'initial_coordinate_seed_receipt.json',seed_receipt)
        completed_years=[]
        if 'resume_history' in manifest:
            resumed=load_resume_history(manifest['resume_history'],current_manifest=manifest,
                targets=targets,tolerance=tolerance,case=args.case,count=args.count,
                initial_checkpoint_sha256=digest,source_root=root)
            realized=resumed['realized'];inherited=resumed['inherited'];initial=resumed['initial']
            completed_years=resumed['completed_years']
            save(out/'resume_receipt.json',dict(status='verified_completed_history_prefix',
                completed_years=completed_years,source_manifest=resumed['source_manifest'],
                source_contract=resumed['source_contract'],last_realized_state=resumed['last_realized_state'],
                inherited_year=inherited.year,horizon_verified=False,production_eligible=False))
            save(out/'realized_fit.json',realized)
        fit_deadline=deadline-float(manifest.get('policy_reserve_seconds',min(3*3600,args.seconds*.25)))
        alternative_used=False;smoked=False
        for target in targets:
            year=int(target['decision_year']);desired=float(target['period_tfr_arithmetic_mean'])
            if year in completed_years:continue
            if inherited.year!=year:raise RuntimeError('Inherited historical clock mismatch')
            trials=[];seen=[];winner=None;forecast_jacobian=None
            center=initial_psi if not realized else realized[-1]['psi']
            seed=initial_seed_step(plan,manifest)
            for attempt in range(min(6,int(plan.get('maximum_trials_per_window',6)))):
                if time.monotonic()>=fit_deadline:raise TimeoutError('Historical fit deadline reached')
                psi=float(np.clip(next_psi(trials,center,seed,(initial_psi-.20,initial_psi+.02)),initial_psi-.20,initial_psi+.02))
                if any(abs(psi-q)<1e-8 for q in seen):
                    psi=float(np.clip(center+(attempt+1)*seed,initial_psi-.20,initial_psi+.02))
                    if any(abs(psi-q)<1e-8 for q in seen):break
                seen.append(psi);folder=out/f'window_{year}'/f'trial_{attempt:02d}'
                try:
                    result,detail=solve_forecast(inherited=inherited,old=old,demographics=demographics,
                        psi=psi,count=args.count,initial=initial,controls=controls,audit=audit,
                        deadline=min(fit_deadline,time.monotonic()+float(manifest.get('forecast_seconds',fit_deadline-time.monotonic()))),folder=folder,case=args.case,
                        initial_jacobian=forecast_jacobian)
                    if result.next_state is None and not alternative_used and result.root_receipt.get('best'):
                        forecast_jacobian=None
                        alternative_used=True;initial=np.asarray(result.root_receipt['best']['prices'])
                        alt=dict(controls,damping=float(controls['damping'])*.5,
                            automatic_fiscal_polish=False)
                        save(folder/'alternative_start.json',dict(reason='Best admissible root coordinates with half damping; sole track alternative'))
                        result,detail=solve_forecast(inherited=inherited,old=old,demographics=demographics,
                            psi=psi,count=args.count,initial=initial,controls=alt,audit=audit,
                            deadline=fit_deadline,folder=folder/'alternative',case=args.case,
                            initial_jacobian=forecast_jacobian)
                    if result.next_state is None:
                        forecast_jacobian=None
                        save(folder/'rejected.json',dict(reason='Finite market/fiscal root failed'));continue
                    forecast_jacobian=reusable_forecast_jacobian(
                        result,args.count,reuse_forecast_jacobian)
                    initial=np.asarray(result.root_receipt['final']['prices'])
                    if not smoked:
                        standard_graphs(detail['snapshot'],result,folder/'native_graphs')
                        save(out/'native_smoke.json',dict(passed=True,count=args.count,seconds=time.monotonic()-started,
                            folder=str(folder),horizon_verified=False,boundary_status='unverified_finite_truncation'))
                        smoked=True
                    value=float(detail['observations'][0]['period_tfr_topcode_adjusted']);gap=value-desired
                    trials.append((psi,gap));fit=dict(year=year,psi=psi,target=desired,model=value,gap=gap,folder=str(folder))
                    save(folder/'fit.json',fit)
                    if winner is None or abs(gap)<abs(winner[0]['gap']):winner=(fit,result,detail)
                    save(out/'best_so_far.json',dict(realized=realized,current=winner[0]))
                    if abs(gap)<=tolerance:break
                except (RuntimeError,ValueError,TimeoutError) as exc:
                    forecast_jacobian=None
                    save(folder/'rejected.json',dict(error_type=type(exc).__name__,error=str(exc)))
                    if isinstance(exc,TimeoutError):raise
            if winner is None or abs(winner[0]['gap'])>tolerance:
                raise RuntimeError('No fitted admissible forecast; failed window is not inherited')
            fit,result,detail=winner
            standard_graphs(detail['snapshot'],result,out/f'window_{year}'/'selected_graphs')
            realized.append(fit);inherited=result.next_state
            initial=shift_forecast_coordinates(result.root_receipt['final']['prices'],args.count)
            checkpoint(out/f'realized_state_{inherited.year}.pkl.gz',inherited)
            save(out/'realized_fit.json',realized)
        save(out/'finite_history_complete.json',dict(realized=realized,horizon_verified=False,production_eligible=False))
        for annual_tax in (.01,.02):
            folder=out/'policies'/('baseline_rebate' if annual_tax==.01 else 'tax2_rebate')
            if time.monotonic()>=deadline:break
            policy_old=copy.copy(old);policy_old.parameters=copy.deepcopy(old.parameters)
            P=policy_old.parameters;P.tau_H=4.*annual_tax
            P.user_cost_rate=float(P.R_gross)+float(P.delta)+float(P.tau_H)-1.
            try:
                result,detail=solve_forecast(inherited=inherited,old=policy_old,demographics=demographics,
                    psi=realized[-1]['psi'],count=args.count,initial=initial,controls=controls,audit=audit,
                    deadline=deadline,folder=folder,case=args.case)
                if result.next_state is not None:standard_graphs(detail['snapshot'],result,folder/'graphs')
                save(folder/'summary.json',dict(finite_converged=result.next_state is not None,
                    annual_tax=annual_tax,horizon_verified=False,production_eligible=False))
            except (RuntimeError,ValueError,TimeoutError) as exc:
                save(folder/'failure.json',dict(error_type=type(exc).__name__,error=str(exc)))
        save(out/'summary.json',dict(status='finite_history_complete_policy_attempts_finished',case=args.case,
            count=args.count,realized=realized,horizon_verified=False,production_eligible=False))
    except Exception as exc:
        save(out/'failure.json',dict(error_type=type(exc).__name__,error=str(exc),realized=realized,
            elapsed_seconds=time.monotonic()-started));raise
    finally:
        stop.set()


if __name__=='__main__':main()
