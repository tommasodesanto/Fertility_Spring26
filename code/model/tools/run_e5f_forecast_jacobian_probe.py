#!/usr/bin/env python3
"""Run one bounded forecast from a pinned accepted root and its learned Jacobian."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import time
from unittest.mock import patch

import numpy as np


class _ProbeFinished(BaseException):
    pass


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def load(name, path):
    spec=importlib.util.spec_from_file_location(name,Path(path).resolve())
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


def accepted_seed(receipt, case, count, reproduction_tolerance):
    width=3*(count+1);final=receipt.get('final');best=receipt.get('best')
    if (receipt.get('converged') is not True or receipt.get('status')!='converged'
            or receipt.get('finite_horizon_market_fiscal_converged') is not True
            or receipt.get('start_year')!=2007 or receipt.get('case')!=case
            or receipt.get('count')!=count or not isinstance(final,dict) or not isinstance(best,dict)
            or final.get('mapping_valid') is not True):
        raise ValueError('Seed is not an accepted same-case/count 2007 forecast root')
    prices=np.asarray(final.get('prices'),dtype=float)
    best_prices=np.asarray(best.get('prices'),dtype=float)
    jacobian=np.asarray(receipt.get('final_jacobian'),dtype=float)
    reproduction=float(receipt.get('final_reproduction_max_abs',math.nan))
    if (prices.shape!=(width,) or best_prices.shape!=(width,)
            or not np.isfinite(prices).all() or np.any(prices<=0)
            or not np.array_equal(prices,best_prices)
            or jacobian.shape!=(width,width) or not np.isfinite(jacobian).all()
            or not np.isfinite(reproduction) or reproduction>reproduction_tolerance):
        raise ValueError('Seed root lacks a finite exact replay, prices, or learned Jacobian')
    damping=float(receipt.get('final_damping',math.nan))
    psi=float(receipt.get('psi',math.nan))
    if not np.isfinite([damping,psi]).all() or not 0<damping<=1:
        raise ValueError('Seed root damping or psi is invalid')
    return prices.copy(),jacobian.copy(),psi


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--manifest-sha256',required=True)
    parser.add_argument('--driver',type=Path,required=True)
    parser.add_argument('--driver-sha256',required=True)
    parser.add_argument('--cache',type=Path,required=True)
    parser.add_argument('--cache-sha256',required=True)
    parser.add_argument('--cache-proof',type=Path,required=True)
    parser.add_argument('--cache-proof-sha256',required=True)
    parser.add_argument('--seed-root',type=Path,required=True)
    parser.add_argument('--seed-root-sha256',required=True)
    parser.add_argument('--seed-contract-sha256',required=True)
    parser.add_argument('--case',choices=['A0','A+'],required=True)
    parser.add_argument('--count',type=int,choices=[6,24,100],required=True)
    parser.add_argument('--psi',type=float)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--seconds',type=float,default=3600.)
    args=parser.parse_args(argv)
    if not np.isfinite(args.seconds) or not 0<args.seconds<=3600:
        parser.error('--seconds must be in (0,3600]')
    pins=((args.manifest,args.manifest_sha256,'manifest'),
          (args.driver,args.driver_sha256,'history driver'),
          (args.cache,args.cache_sha256,'policy cache'),
          (args.cache_proof,args.cache_proof_sha256,'cache proof'),
          (args.seed_root,args.seed_root_sha256,'seed root'))
    for path,digest,label in pins:
        if sha(path)!=digest:raise ValueError(f'Pinned {label} changed')
    out=args.output.resolve()
    if out.exists():raise FileExistsError(f'Probe output already exists: {out}')
    out.mkdir(parents=True)

    manifest=read(args.manifest);proof=read(args.cache_proof);seed=read(args.seed_root)
    driver_path=args.driver.resolve();cache_path=args.cache.resolve()
    if (proof.get('status')!='verified' or proof.get('cache_sha256')!=args.cache_sha256
            or proof.get('exact_mapping_equal') is not True
            or proof.get('household_and_accounting_mapping_valid') is not True):
        raise ValueError('Native exact-cache proof does not match its source')
    seed_contract_path=args.seed_root.resolve().parents[2]/'contract_receipt.json'
    if sha(seed_contract_path)!=args.seed_contract_sha256:
        raise ValueError('Pinned seed history contract changed')
    seed_contract=read(seed_contract_path)
    if (seed_contract.get('manifest_sha256')!=args.manifest_sha256
            or seed_contract.get('case')!=args.case or seed_contract.get('count')!=args.count):
        raise ValueError('Seed root was not run under this manifest, case, and count')
    plan=read(manifest['prior_plan']);summary=read(manifest['initial_summary'])
    checkpoint=summary['checkpoint'];checkpoint_sha=checkpoint.get('sha256',checkpoint.get('checkpoint_sha256'))
    if (seed_contract.get('source_root')!=plan['source_root']
            or seed_contract.get('initial_checkpoint_sha256')!=checkpoint_sha):
        raise ValueError('Seed model source or initial checkpoint differs from this manifest')
    controls=dict(plan['history_root_controls']);controls.update(manifest.get('root_controls',{}))
    prices,jacobian,seed_psi=accepted_seed(
        seed,args.case,args.count,float(controls['final_reproduction_tolerance']))
    psi=seed_psi-.01 if args.psi is None else float(args.psi)
    if not np.isfinite(psi):raise ValueError('Probe psi must be finite')

    driver=load('run_e5f_final_rebated_history',driver_path)
    cache=load('e5f_exact_policy_cache',cache_path)
    probe_manifest=dict(manifest)
    probe_manifest['policy_reserve_seconds']=0
    probe_manifest['forecast_seconds']=args.seconds
    probe_manifest['forecast_jacobian_probe']=dict(
        canonical_manifest=str(args.manifest.resolve()),canonical_manifest_sha256=args.manifest_sha256,
        seed_root=str(args.seed_root.resolve()),seed_root_sha256=args.seed_root_sha256,
        seed_contract=str(seed_contract_path),seed_contract_sha256=args.seed_contract_sha256,
        wrapper_only_budget_override=True)
    probe_manifest_path=out/'probe_manifest.json';driver.save(probe_manifest_path,probe_manifest)
    driver.save(out/'probe_contract.json',dict(case=args.case,count=args.count,seconds=args.seconds,
        seed_psi=seed_psi,probe_psi=psi,initial_prices_source='accepted seed final exact replay',
        initial_jacobian_source='accepted seed final learned Jacobian',source_pins={str(p):d for p,d,_ in pins},
        seed_contract_sha256=args.seed_contract_sha256,canonical_manifest_sha256=args.manifest_sha256,
        cache_native_proof_sha256=args.cache_proof_sha256,
        probe_driver_sha256=sha(Path(__file__).resolve()),maximum_forecast_solves=1))

    captured={};original=driver.solve_forecast;started=time.monotonic()
    def one_forecast(**kwargs):
        if captured.get('called'):raise RuntimeError('Probe attempted more than one forecast')
        captured['called']=True;captured['folder']=Path(kwargs['folder'])
        if (kwargs['case']!=args.case or kwargs['count']!=args.count
                or kwargs['inherited'].year!=2007):
            captured['error']=ValueError('First intercepted forecast differs from seed case/count/year')
            raise _ProbeFinished()
        call=dict(kwargs,initial=prices.copy(),initial_jacobian=jacobian.copy(),psi=psi,
                  deadline=min(kwargs['deadline'],time.monotonic()+args.seconds))
        import e5f_rebated_surprises as rebated
        _,joined,*_=rebated._runtime()
        try:
            cache_gib=24 if args.count==100 else 6
            with cache.policy_cache(joined.pf,max_bytes=cache_gib*1024**3) as stats:
                try:result,detail=original(**call)
                finally:captured['cache']=stats.snapshot()
            captured.update(result=result,detail=detail)
        except Exception as exc:
            captured['error']=exc
        raise _ProbeFinished()

    try:
        with patch.object(driver,'solve_forecast',one_forecast):
            driver.main(['--manifest',str(probe_manifest_path),'--case',args.case,
                         '--count',str(args.count),'--output',str(out),'--seconds',str(args.seconds)])
    except _ProbeFinished:
        pass
    elapsed=time.monotonic()-started
    driver.save(out/'elapsed.json',dict(elapsed_seconds=elapsed,maximum_seconds=args.seconds,
        forecast_calls=1 if captured.get('called') else 0))
    if not captured.get('called'):
        raise RuntimeError('History driver reached no forecast solve')
    if 'error' in captured:
        exc=captured['error'];driver.save(out/'failure.json',dict(status='failed_forecast_jacobian_probe',
            error_type=type(exc).__name__,error=str(exc),elapsed_seconds=elapsed,cache=captured.get('cache')))
        raise RuntimeError('Single forecast Jacobian probe failed') from exc
    result,detail=captured['result'],captured['detail'];receipt=result.root_receipt
    valid=bool(result.next_state is not None and receipt.get('finite_horizon_market_fiscal_converged'))
    if not valid:
        driver.save(out/'failure.json',dict(status='failed_forecast_jacobian_probe',
            error='Forecast did not pass the unchanged finite root and replay gates',elapsed_seconds=elapsed,
            root_receipt=str(Path(detail.get('folder',''))/'root_receipt.json'),cache=captured.get('cache')))
        raise RuntimeError('Single forecast Jacobian probe did not converge')
    from run_e5f_successive_surprises_overnight import standard_graphs
    graph_folder=out/'graphs';standard_graphs(detail['snapshot'],result,graph_folder)
    graphs=sorted(str(path) for path in (graph_folder/'standard_diagnostics').glob('*.png'))
    if len(graphs)!=17:raise RuntimeError('Probe did not write all 17 standard graphs')
    root_receipt_path=captured['folder']/'root_receipt.json'
    root_summary={key:driver.clean(value) for key,value in receipt.items()
                  if key not in ('history','best','final','final_jacobian')}
    root_summary.update(root_receipt=str(root_receipt_path),root_receipt_sha256=driver.sha(root_receipt_path),
        seed_final_jacobian_reused=True)
    driver.save(out/'probe_root_summary.json',root_summary)
    observation=detail['observations'][0]
    driver.save(out/'probe_evaluation.json',dict(calendar_year=2007,psi=psi,
        period_fertility=observation,root_evaluations=receipt['evaluations'],
        mapping_receipt=str(captured['folder']/'mapping_receipt.json'),
        cache=captured['cache'],graphs=graphs))
    driver.save(out/'valid.json',dict(status='verified',finite_horizon_market_fiscal_converged=True,
        exact_replay_verified=True,historical_state_carried=False,root_receipt=str(root_receipt_path)))
    summary_out=dict(status='verified_forecast_jacobian_probe',case=args.case,count=args.count,
        elapsed_seconds=elapsed,seed_psi=seed_psi,probe_psi=psi,
        finite_horizon_market_fiscal_converged=True,exact_replay_verified=True,
        evaluations=receipt['evaluations'],final_reproduction_max_abs=receipt['final_reproduction_max_abs'],
        final_damping=receipt['final_damping'],standard_graph_count=17,
        root_receipt=str(root_receipt_path),
        evaluation=str(out/'probe_evaluation.json'),cache=captured['cache'],
        historical_state_carried=False,production_eligible=False)
    driver.save(out/'summary.json',summary_out);print(json.dumps(summary_out),flush=True)


if __name__=='__main__':main()
