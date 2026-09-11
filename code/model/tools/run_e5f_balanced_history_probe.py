"""Pinned announced-history replay or explicit joint price/pension root; no objective.

With exactly two root evaluations, the existing root samples the supplied path
then repeats it. It does not take a price/pension optimization step. Conditional
mapping/replay success is reported separately from market/fiscal convergence.
The separate actual-root schema permits six or 28 dates and bounded updates.
"""
from __future__ import annotations
import argparse
import copy
from dataclasses import fields
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import pickle
import resource
import sys
import threading
import time

for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[_key]='1'
os.environ.setdefault('MPLBACKEND','Agg')

import run_e5f_balanced_terminal_probe as terminal_driver

ROOT=Path(__file__).resolve().parents[3]
SCHEMA='e5f_balanced_history_probe_v1'
ROOT_SCHEMA='e5f_balanced_history_root_v1'
PIN_NAMES=('initial_checkpoint','initial_summary','initial_contract','terminal_checkpoint',
           'terminal_summary','terminal_contract','terminal_root_receipt')
ROOT_KEYS={'price_bounds','pension_bounds','market_tolerance','fiscal_tolerance',
    'market_slope','fiscal_slope','max_log_step','damping','max_evaluations',
    'max_condition_number','worsening_factor','final_reproduction_tolerance','initial_jacobian'}
digest,verify,pinned_json=terminal_driver.digest,terminal_driver.verify,terminal_driver.pinned_json


def validate_contract(c):
    required={'schema','source_sha256',*PIN_NAMES,'originating_source_commit','count',
        'outside_origin_entry_share','outside_entry_status','psi_change_from_initial','preference_rule',
        'demographic_seed_mode','terminal_demographic_year','seconds','root_seconds',
        'identity_tolerance','initial_fertility_tolerance','audit_controls','root_controls',
        'initial_prices','initial_pensions','standard_graph_count'}
    actual=c.get('schema')==ROOT_SCHEMA
    if actual:required.add('mode')
    if (set(c)!=required or c['schema'] not in (SCHEMA,ROOT_SCHEMA)
            or (actual and c['mode']!='actual_root')):
        raise ValueError('Exact explicit history replay or actual-root schema and fields required')
    if (c['originating_source_commit']!='c6dd3508' or type(c['count']) is not int or c['count'] not in ((6,28) if actual else (6,))
            or c['outside_origin_entry_share']!=.169 or c['outside_entry_status']!='diagnostic_outstanding_not_estimated'
            or c['psi_change_from_initial']!=-.25
            or c['preference_rule']!='announced_linear_2007_2023_then_constant_diagnostic'
            or c['demographic_seed_mode']!='serialized_frozen_2023_primitives_without_realignment'
            or c['terminal_demographic_year']!=2100 or c['standard_graph_count']!=17):
        raise ValueError('Require approved diagnostic timing, explicit entry/preference and demographic declarations')
    ceiling=7200 if actual and c['count']==28 else 1800
    if (type(c['seconds']) is not int or type(c['root_seconds']) is not int
            or not 1<=c['seconds']<=ceiling or not 1<=c['root_seconds']<=c['seconds']-120
            or (actual and c['count']==6 and c['root_seconds']>1620)):
        raise ValueError('Explicit bounded root/watchdog and120 seconds reporting reserve required')
    for name,ceiling in [('identity_tolerance',2e-9),('initial_fertility_tolerance',5e-4)]:
        if not math.isfinite(c[name]) or not 0<c[name]<=ceiling:
            raise ValueError(f'Invalid {name}')
    for name in PIN_NAMES:
        spec=c[name]
        if (set(spec)!={'path','sha256'} or not Path(spec['path']).is_absolute()
                or not isinstance(spec['sha256'],str) or len(spec['sha256'])!=64
                or any(ch not in '0123456789abcdef' for ch in spec['sha256'])):
            raise ValueError(f'Explicit absolute path and completed SHA256 required: {name}')
    if (c['initial_checkpoint']['sha256']!=terminal_driver.INITIAL_SHA256
            or Path(c['initial_checkpoint']['path']).parts[-3:]!=('new_balanced_smoke','repetition_02','initial_state.pkl.gz')):
        raise ValueError('Use the verified17358109 repetition_02 initial checkpoint')
    initial_parent=Path(c['initial_checkpoint']['path']).parent
    terminal_parent=Path(c['terminal_checkpoint']['path']).parent
    if (Path(c['initial_summary']['path'])!=initial_parent/'summary.json'
            or Path(c['terminal_checkpoint']['path']).name!='terminal_state.pkl.gz'
            or Path(c['terminal_summary']['path'])!=terminal_parent/'summary.json'
            or Path(c['terminal_root_receipt']['path'])!=terminal_parent/'root_receipt.json'):
        raise ValueError('Each result packet must use its own canonical checkpoint/summary/receipt directory')
    r=c['root_controls']
    if set(r)!=ROOT_KEYS or type(r['max_evaluations']) is not int or r['max_evaluations'] not in (range(3,9) if actual else (2,)) or r['initial_jacobian'] is not None:
        raise ValueError('Explicit bounded evaluation count and no supplied Jacobian required')
    for name,bound in [('initial_prices','price_bounds'),('initial_pensions','pension_bounds')]:
        values=c[name]
        lo,hi=r[bound]
        if (not isinstance(values,list) or len(values)!=c['count'] or not all(math.isfinite(v) and v>0 for v in values)
                or not all(math.isfinite(v) for v in (lo,hi)) or not 0<lo<hi
                or not all(lo<=v<=hi for v in values)):
            raise ValueError(f'All dated {name} must be explicit and within declared bounds')


def validate_parents(c, initial, initial_summary, initial_contract, terminal, terminal_summary, terminal_contract, root):
    if (initial['contract_sha256']!=c['initial_contract']['sha256']
            or initial_contract.get('schema')!='e5f_parenthood_initial_probe_v1'
            or initial_contract.get('normalize') is not True
            or initial_summary.get('status')!='passed_initial_diagnostic'
            or initial_summary.get('checkpoint_sha256')!=c['initial_checkpoint']['sha256']):
        raise ValueError('Initial packet requires its ORIGINAL input contract and verified summary')
    if (terminal.get('schema')!='e5f_balanced_terminal_probe_checkpoint_v1'
            or terminal['contract_sha256']!=c['terminal_contract']['sha256']
            or terminal_summary.get('status')!='passed_terminal_root_diagnostic'
            or terminal_summary.get('endpoint_numerically_verified') is not True
            or not terminal_summary.get('final_checks') or not all(terminal_summary['final_checks'].values())
            or terminal_summary['checkpoint_sha256']!=c['terminal_checkpoint']['sha256']
            or terminal_contract.get('schema')!=terminal_driver.SCHEMA):
        raise ValueError('Terminal packet must be an accepted balanced root with its original input contract')
    for name in ('initial_checkpoint','initial_summary','initial_contract'):
        if terminal_contract[name]!=c[name]:
            raise ValueError(f'Terminal uses another normalized initial input: {name}')
    if (terminal_contract['psi_change_from_initial']!=c['psi_change_from_initial']
            or terminal_contract['terminal_demographic_year']!=c['terminal_demographic_year']
            or terminal_contract['demographic_seed_mode']!=c['demographic_seed_mode']):
        raise ValueError('Terminal preference or demographic continuation differs')
    prior_sources=terminal_contract.get('source_sha256',{})
    # No shared-file changes are authorized by this first probe. New driver,
    # adapter and test files may be added; every prior source remains identical.
    # A future diagnostic-only exception must enumerate exact reviewed paths.
    if not prior_sources or any(c['source_sha256'].get(name)!=pin for name,pin in prior_sources.items()):
        raise ValueError('Terminal and history economic source pins differ')
    if (root.get('schema')!='e5f_balanced_terminal_v1' or not root.get('converged')
            or not root.get('endpoint_production_eligible') or not root.get('fresh_endpoint_matches_final')
            or not root.get('gates') or not all(root['gates'].values())):
        raise ValueError('Converged and freshly repeated terminal root required')
    return [dict(path=name,sha256=pin) for name,pin in sorted(prior_sources.items())]


def validate_replay(c, receipt, signatures, records):
    """Associate the root-selected point with its own fresh final path."""
    def values(v):return v.tolist() if hasattr(v,'tolist') else v
    count=c['count'];actual=c['schema']==ROOT_SCHEMA
    evaluations=receipt['evaluations']
    if (type(evaluations) is not int or not 2<=evaluations<=c['root_controls']['max_evaluations']
            or (not actual and evaluations!=2) or len(signatures)!=evaluations
            or len(records)!=evaluations or len(receipt['history'])!=evaluations
            or not receipt['fresh_path_matches_final']):
        raise RuntimeError('Missing complete mappings or fresh final association')
    years=list(range(2007,2007+4*count,4))
    if any(len(rows)!=count or [row['year'] for row in rows]!=years for rows in signatures):
        raise RuntimeError('Incomplete dated policy/distribution signatures')
    if [row['evaluation'] for row in records]!=list(range(1,evaluations+1)):
        raise RuntimeError('Mapping records must count each evaluation exactly once')
    best,final=receipt['best'],receipt['final']
    if best is None or final is None:
        raise RuntimeError('No selected mapping and reserved fresh final replay')
    chosen=best['payload']['trial'];repeated=final['payload']['trial']
    if (type(chosen) is not int or type(repeated) is not int
            or not 1<=chosen<repeated or repeated!=evaluations
            or records[-1]['phase']!='final'):
        raise RuntimeError('Selected/fresh trial identifiers are inconsistent')
    for point,trial in ((best,chosen),(final,repeated)):
        recorded=records[trial-1];historical=receipt['history'][trial-1]
        if (not point['mapping_valid'] or point['payload']['bellman_solves']!=2*count
                or not all(a['passed'] for a in point['payload']['mapping_gates'].values())
                or len(point['payload']['dated_household_audits'])!=count
                or not all(all(a['gates'].values()) for a in point['payload']['dated_household_audits'])):
            raise RuntimeError('Selected/final mapping audits or Bellman count failed')
        for name in ('prices','fiscal_values','market_residual','fiscal_residual'):
            if values(point[name])!=values(recorded[name]) or values(point[name])!=values(historical[name]):
                raise RuntimeError('Saved mapping coordinates differ from selected trial')
    if (any(values(best[k])!=values(final[k]) for k in ('prices','fiscal_values'))
            or signatures[chosen-1]!=signatures[repeated-1]
            or not all(receipt['gates'][k] for k in ('mapping','market_replay','fiscal_replay'))):
        raise RuntimeError('Selected whole-path policy/distribution/fiscal replay failed')
    tolerance=c['root_controls']['final_reproduction_tolerance']
    for name in ('market_residual','fiscal_residual'):
        left,right=values(best[name]),values(final[name])
        if (len(left)!=count or len(right)!=count
                or any(not math.isfinite(a) or not math.isfinite(b) or abs(a-b)>tolerance
                       for a,b in zip(left,right))):
            raise RuntimeError('Selected/final numerical residual reproduction failed')
    for row in records:
        if row['bellman_solves']!=2*count:
            raise RuntimeError('Completed mapping Bellman count differs from announced loop')
    return dict(verified=True,selected_trial=chosen,fresh_final_trial=repeated,
        completed_mappings=evaluations,bellman_solves=sum(row['bellman_solves'] for row in records))


def run(args):
    verify(args.contract,args.contract_sha256)
    c=json.loads(args.contract.read_text());validate_contract(c)
    actual=c['schema']==ROOT_SCHEMA;count=c['count']
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    started=time.monotonic();finished=threading.Event();lock=threading.Lock()
    state=dict(phase='preflight',completed_mappings=0)
    def jsonable(v):
        if isinstance(v,dict):return {str(k):jsonable(x) for k,x in v.items()}
        if isinstance(v,(list,tuple)):return [jsonable(x) for x in v]
        if hasattr(v,'tolist'):return jsonable(v.tolist())
        if isinstance(v,Path):return str(v)
        return v
    def save(name,value):
        with lock:
            p=out/name;q=p.with_suffix(p.suffix+'.tmp')
            q.write_text(json.dumps(jsonable(value),indent=2,sort_keys=True)+'\n');os.replace(q,p)
    def watchdog():
        while not finished.wait(min(30.,max(.01,c['seconds']-(time.monotonic()-started)))):
            elapsed=time.monotonic()-started;save('heartbeat.json',dict(state,elapsed_seconds=elapsed))
            if elapsed>=c['seconds']:
                save('failure.json',dict(state,status='timeout',elapsed_seconds=elapsed));os._exit(124)
    threading.Thread(target=watchdog,daemon=True).start()
    save('contract.json',dict(c,contract_sha256=args.contract_sha256));save('heartbeat.json',state)
    try:
        verified_sources=terminal_driver.verify_sources(c)
        initial_summary=pinned_json(c['initial_summary']);initial_contract=pinned_json(c['initial_contract'])
        terminal_summary=pinned_json(c['terminal_summary']);terminal_contract=pinned_json(c['terminal_contract'])
        root=pinned_json(c['terminal_root_receipt'])
        for name in ('initial_checkpoint','terminal_checkpoint'):verify(c[name]['path'],c[name]['sha256'])
        sys.path[:0]=[str(ROOT/'code/model'),str(ROOT/'code/model/tools')]
        import numpy as np
        import e5f_balanced_history as adapter
        import e5f_balanced_terminal as terminal_adapter
        from e5f_approved_initial_state import build_approved_initial_state
        import run_e5f_matched_pf_smoke as primitive
        import run_e5f_perfect_foresight_person_demography as person
        from e5f_social_security import fiscal_accounts
        adapter._runtime()  # Activate the identical sequential hooks before bridge arithmetic.
        with gzip.open(c['initial_checkpoint']['path'],'rb') as f:initial=pickle.load(f)
        with gzip.open(c['terminal_checkpoint']['path'],'rb') as f:terminal_packet=pickle.load(f)
        terminal_source_comparison=validate_parents(c,initial,initial_summary,initial_contract,
            terminal_packet,terminal_summary,terminal_contract,root)
        if jsonable(terminal_packet['root_receipt'])!=root:
            raise ValueError('Serialized terminal numerical receipt differs from its pinned JSON')
        old=build_approved_initial_state(packet=initial,normalization=initial_summary['normalization'],
            outside_origin_entry_share=c['outside_origin_entry_share'],
            preference_change_2023=c['psi_change_from_initial'],fertility_tolerance=c['initial_fertility_tolerance'])
        demographics=initial['demographic_seed'];terminal_demographics=terminal_packet['demographic_seed']
        demography_receipt=terminal_driver.validate_demographics(demographics,old.parameters,c,person,np)
        payload=root['final']['payload']
        terminal=terminal_adapter.BalancedTerminalEndpoint(terminal_packet['parameters'],terminal_packet['b_grid'],
            terminal_packet['policy'],terminal_packet['endpoint'],terminal_packet['social_security'],
            payload['diagnostics'],payload['household_gates'])
        audit_controls=terminal_adapter.TerminalAuditControls(**c['audit_controls'])
        if set(c['audit_controls'])!={f.name for f in fields(terminal_adapter.TerminalAuditControls)}:
            raise ValueError('Every household audit control must be explicit')
        save('preflight.json',dict(status='passed',source_files_verified=verified_sources,
            terminal_source_comparison=terminal_source_comparison,allowed_shared_source_differences=[],
            initial_bridge=old.diagnostics,demographics=demography_receipt,
            initial_supply=vars(old.supply_rule),outside_entry_status=c['outside_entry_status'],
            targets_loaded=False,objective_computed=False,production_eligible=False))
        save('sizing.json',dict(dates=count,backward_and_forward_calls_per_mapping=2*count,
            maximum_mappings=c['root_controls']['max_evaluations'],
            maximum_bellman_calls=2*count*c['root_controls']['max_evaluations'],extra_reporting_bellman_calls=0,prior_initial_ge_seconds_approximate=70.,
            six_date_mapping_seconds='unknown: this exact smoke measures it; GE time is not a dated Bellman timing',
            root_budget_seconds=c['root_seconds'],total_budget_seconds=c['seconds'],
            reporting_reserve_seconds=c['seconds']-c['root_seconds']))
        del initial,terminal_packet
        dated_signatures=[];snapshot={};records=[]
        save('latest_completed.json',dict(status='no_completed_mapping'))
        save('best_so_far.json',dict(status='no_valid_mapping'))
        def array_hash(a):
            a=np.asarray(a);h=hashlib.sha256(str((a.shape,a.dtype.str)).encode())
            h.update(np.ascontiguousarray(a).tobytes());return h.hexdigest()
        def observe(index,e,P,grid,shared):
            if index==0:dated_signatures.append([])
            record=dict(year=2007+4*index,policy={k:array_hash(v) for k,v in primitive.policy_arrays(e.policy).items()},
                distributions={k:array_hash(getattr(e,k)) for k in ('g_pre','g_post_fertility','g_current')},
                fiscal=fiscal_accounts(e.g_current,P))
            dated_signatures[-1].append(record)
            state.update(phase='historical_forward',mapping=len(dated_signatures),year=record['year'])
            save('latest_date.json',record)
            if index==4:
                snapshot.clear();snapshot.update(parameters=P,b_grid=grid,evaluation=e,shared=shared,
                    supply_rule=old.supply_rule,demographic_seed=demographics,calendar_year=2023,
                    state_timing='pre-choice g_pre and corresponding current2023 choices retained separately')
        def callback(record):
            if record.get('event')=='complete':save('root_completion.json',record);return
            record=dict(record,bellman_solves=2*len(dated_signatures[record['evaluation']-1]))
            if record.get('event')=='safeguard':
                records[-1]=record;save('root_evaluations.json',records);return
            records.append(record);state.update(phase='between_mappings',completed_mappings=record['evaluation'])
            save('latest_completed.json',record);save('root_evaluations.json',records)
            save('dated_reproduction.json',dict(signatures=dated_signatures))
            if record.get('new_best'):save('best_so_far.json',record)
        result=adapter.solve_balanced_history(old_state=old,terminal=terminal,terminal_root_receipt=root,
            demographic_primitives=demographics,terminal_demographic_primitives=terminal_demographics,
            count=count,initial_prices=c['initial_prices'],initial_pensions=c['initial_pensions'],
            audit_controls=audit_controls,deadline_monotonic=started+c['root_seconds'],
            callback=callback,observer=observe,**c['root_controls'])
        receipt=result.root_receipt;save('root_receipt.json',receipt)
        save('dated_reproduction.json',dict(signatures=dated_signatures))
        replay=validate_replay(c,receipt,dated_signatures,records)
        save('replay_verification.json',replay)
        if result.path is None:raise RuntimeError('No fresh final path to serialize')
        primitive.pf.write_csv(out/'transition_path.csv',result.path.rows)
        save('dated_household_audits.json',receipt['final']['payload']['dated_household_audits'])
        rents=primitive.pf.rents_from_asset_prices(np.asarray(receipt['final']['prices']),
            float(terminal.policy.price[0]),old.parameters)
        snapshot.update(schema='e5f_balanced_history_probe_checkpoint_v1',contract_sha256=args.contract_sha256,
            root_receipt=receipt,actual_renter_price=float(rents[4]),
            historical_terminal_state=result.path.history.terminal_state,
            person_terminal_state=result.path.person_tail.terminal_state,
            production_eligible=False,calibrated_smm=False,horizon_verified=False)
        state['phase']='checkpoint_and_graphs'
        checkpoint=out/'dated_2023.pkl.gz'
        with gzip.open(checkpoint,'wb',compresslevel=1) as f:pickle.dump(snapshot,f,protocol=5)
        with gzip.open(checkpoint,'rb') as f:loaded=pickle.load(f)
        for name,value in primitive.policy_arrays(snapshot['evaluation'].policy).items():
            np.testing.assert_array_equal(value,primitive.policy_arrays(loaded['evaluation'].policy)[name])
        for name in ('g_pre','g_post_fertility','g_current'):
            np.testing.assert_array_equal(getattr(snapshot['evaluation'],name),getattr(loaded['evaluation'],name))
        for name in ('persons','heads'):
            np.testing.assert_array_equal(getattr(loaded['person_terminal_state'].persons,name),
                getattr(snapshot['person_terminal_state'].persons,name))
        np.testing.assert_array_equal(loaded['person_terminal_state'].g_pre,snapshot['person_terminal_state'].g_pre)
        np.testing.assert_array_equal(loaded['historical_terminal_state'].g_pre,snapshot['historical_terminal_state'].g_pre)
        if (loaded['contract_sha256']!=args.contract_sha256 or loaded['calendar_year']!=2023
                or fiscal_accounts(loaded['evaluation'].g_current,loaded['parameters'])!=dated_signatures[-1][4]['fiscal']):
            raise RuntimeError('Dated checkpoint reload fiscal/timing/contract mismatch')
        del loaded
        import run_e5f_independent_numerical_audit as audit
        from unittest.mock import patch
        original_writer=audit.write_diagnostics
        def dated_writer(stats,P,destination):
            stats.owner_user_cost=np.array([snapshot['actual_renter_price']])
            return original_writer(stats,P,destination)
        with patch.object(audit,'write_diagnostics',dated_writer):
            audit.standard_diagnostics(snapshot,out,validate_production_young=False)
        graphs=sorted((out/'standard_diagnostics').glob('*.png'))
        if len(graphs)!=17:raise RuntimeError('Stable17-graph packet incomplete')
        save('diagnostics_receipt.json',dict(standard_graph_count=17,
            graph_sha256={str(p.relative_to(out)):digest(p) for p in graphs},
            actual_dated_renter_price=snapshot['actual_renter_price'],
            stationary_carrying_cost_reference=float(snapshot['parameters'].user_cost_rate)*float(snapshot['evaluation'].policy.price[0]),
            price_panel='existing user-cost field displays actual PF date rent; no parameter modification',
            supplemental_tables='legacy diagnostic statistics only; not a loaded target system or objective'))
        maxrss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        summary=dict(status=('converged_finite_history_root_diagnostic' if receipt['finite_horizon_market_fiscal_converged']
                else 'incomplete_finite_history_root_diagnostic') if actual else 'reproduced_conditional_history_mapping',
            mapping_replay_verified=True,
            finite_horizon_market_fiscal_converged=receipt['finite_horizon_market_fiscal_converged'],
            market_gate=receipt['gates']['housing'],social_security_gate=receipt['gates']['social_security'],
            root_status=receipt['status'],maximum_market_residual=float(np.max(np.abs(receipt['final']['market_residual']))),
            maximum_fiscal_residual=float(np.max(np.abs(receipt['final']['fiscal_residual']))),
            root_evaluations=receipt['evaluations'],bellman_solves=replay['bellman_solves'],standard_graph_count=17,checkpoint_reload_verified=True,
            checkpoint_sha256=digest(checkpoint),elapsed_seconds=time.monotonic()-started,
            maximum_rss_bytes=int(maxrss if sys.platform=='darwin' else maxrss*1024),
            historical_path_fitted=False,horizon_verified=False,production_eligible=False,calibrated_smm=False,
            terminal_distance=receipt['terminal_distance'],
            interpretation=('bounded announced joint price/pension root; finite-path gates distinct from uncertified horizon'
                if actual else 'announced six-date prescribed path with exact replay; fiscal/market imbalance remains explicit'))
        save('summary.json',summary)
        return 0
    except Exception as exc:
        save('failure.json',dict(state,status='failed_history_probe',error_type=type(exc).__name__,
            error=str(exc),elapsed_seconds=time.monotonic()-started,production_eligible=False))
        raise
    finally:finished.set()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--contract',type=Path,required=True);p.add_argument('--contract-sha256',required=True)
    p.add_argument('--output',type=Path,required=True)
    return run(p.parse_args())


if __name__=='__main__':raise SystemExit(main())
