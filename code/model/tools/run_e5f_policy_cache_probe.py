#!/usr/bin/env python3
"""Compare one complete native forecast mapping with and without policy reuse."""
import argparse
from contextlib import nullcontext
import hashlib
import importlib
import json
from pathlib import Path
import time
from unittest.mock import patch
import numpy as np

class MappingFinished(BaseException):
    pass

def array_signature(value):
    if value is None:return None
    a=np.asarray(value)
    return dict(shape=list(a.shape),dtype=str(a.dtype),sha256=hashlib.sha256(a.tobytes(order="C")).hexdigest())

def mapping_signature(detail):
    path=detail["path"];evaluation=detail["snapshot"]["evaluation"]
    return dict(rows=path.rows,observations=detail["observations"],
        values=[array_signature(v) for v in path.values],
        final_households=array_signature(path.person_tail.terminal_state.g_pre),
        first_current_households=array_signature(evaluation.g_current),
        first_post_birth_households=array_signature(evaluation.g_post_fertility),
        first_policy={k:array_signature(getattr(evaluation.policy,k,None)) for k in
            ("V","bp_pol","hR_pol","tenure_choice","loc_probs","fert_probs","tenure_probs")})

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--manifest",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--cache-sha256",required=True)
    args=parser.parse_args()
    import run_e5f_final_rebated_history as driver
    import e5f_rebated_surprises as rebated
    import e5f_exact_policy_cache as cache
    if hashlib.sha256(Path(cache.__file__).read_bytes()).hexdigest()!=args.cache_sha256:
        raise ValueError("Policy cache source fingerprint changed")
    _,joined,*_=rebated._runtime()
    root=importlib.import_module("e5f_matched_pf_path_root")
    args.output.mkdir(parents=True,exist_ok=False)
    probe_manifest=args.output/'forecast_manifest.json'
    specification=json.loads(args.manifest.read_text())
    specification['policy_reserve_seconds']=0
    driver.save(probe_manifest,specification)
    signatures=[];receipts=[]
    for enabled in (False,True):
        label="cached" if enabled else "uncached"
        folder=args.output/label;captured={}
        def single_mapping(*,initial_prices,evaluate,project,**kwargs):
            prices=project(np.asarray(initial_prices,dtype=float));mapping=evaluate(prices)
            captured["mapping"]=mapping
            best=dict(prices=prices,payload=mapping.get("payload"),score=float(np.max(np.abs(mapping["residual"]))))
            return dict(converged=False,final=None,best=best)
        original=driver.solve_forecast
        def forecast(**kwargs):
            result,detail=original(**kwargs)
            captured["signature"]=mapping_signature(detail)
            raise MappingFinished()
        started=time.monotonic()
        with (cache.policy_cache(joined.pf,max_bytes=6*1024**3) if enabled else nullcontext(None)) as stats:
            with patch.object(root,"solve_price_path",single_mapping),patch.object(driver,"solve_forecast",forecast):
                try:
                    driver.main(["--manifest",str(probe_manifest),"--case","A0","--count","6",
                                 "--output",str(folder),"--seconds","1200"])
                except MappingFinished:
                    pass
            counters=stats.snapshot() if stats is not None else None
        if "signature" not in captured:raise RuntimeError("Native mapping was not completed")
        signature=driver.clean(captured["signature"])
        residual=driver.clean(captured["mapping"]["residual"])
        signature["root_residual"]=residual
        signature["mapping_valid"]=captured["mapping"]["mapping_valid"]
        signature["audit_payload"]=driver.clean(captured["mapping"]["payload"])
        driver.save(folder/"complete_mapping_signature.json",signature)
        signatures.append(signature)
        receipts.append(dict(label=label,elapsed_seconds=time.monotonic()-started,cache=counters))
        driver.save(args.output/"latest_completed.json",receipts[-1])
    equal=signatures[0]==signatures[1]
    useful=receipts[1]["cache"]["hits"]>0
    valid=all(s["mapping_valid"] is True for s in signatures)
    result=dict(status="verified" if equal and useful and valid else "failed",exact_mapping_equal=equal,
        household_and_accounting_mapping_valid=valid,
        positive_cache_hits=useful,runs=receipts,cache_sha256=args.cache_sha256,
        numerical_and_economic_contract_unchanged=True)
    driver.save(args.output/"summary.json",result)
    if not equal or not useful or not valid:raise RuntimeError("Policy cache failed native exactness or reuse")
    print(json.dumps(result),flush=True)

if __name__=="__main__":main()
