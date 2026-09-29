#!/usr/bin/env python3
"""Torch/Slurm-only native readiness smoke for preference-estimation staging.

It is deliberately a code/native baseline proof: it never fits historical
fertility, dispatches work, retries, or certifies a long shocked equilibrium.
"""
from __future__ import annotations
import argparse, gzip, os, pickle, subprocess, sys, threading, time
from pathlib import Path


DEFAULT_ENDPOINT = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_transition_20260928/four_shock_v1/estimation_tests/18759951/native/endpoints/endpoint_1/receipt.json")
DEFAULT_ENDPOINT_CHECKPOINT_SHA256 = "f58ecc3b0e72aba4964b0733d7fe5c46b20cb711cb8d9249bc0142f239f9709a"


def main(a):
    if sys.platform != "linux" or not os.environ.get("SLURM_JOB_ID", "").isdigit():
        raise RuntimeError("Torch Slurm only")
    source, out = a.source_dir.resolve(), a.output.resolve()
    if type(a.cache_max_bytes) is not int or a.cache_max_bytes < 0:
        raise ValueError("Cache budget must be a nonnegative integer number of bytes")
    if out.exists(): raise RuntimeError("output already exists")
    out.mkdir(parents=True); started=time.monotonic(); stop=threading.Event()
    sys.path[:0]=[str(source), "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools"]
    import numpy as np
    import run_e5f_preference_estimation as estimator
    engine=estimator.inner
    pins={name:engine.sha(source/name) for name in estimator.SOURCES}
    contract=estimator.target_contract(a.empirical_blocks,a.annual)
    def write(name,value): engine.write(out/name,value)
    def beat():
        while not stop.wait(60): write("heartbeat.json",dict(epoch=time.time(),phase="native_smoke_backward_work"))
    thread=threading.Thread(target=beat,daemon=True);thread.start()
    try:
        write("latest_completed.json",dict(phase="created",epoch=time.time()))
        tests=["test_e5f_preference_estimation","test_e5f_preference_estimation_batch",
               "test_e5f_preference_transition","test_e5f_four_shock_acceleration",
               "test_e5f_exact_policy_cache"]
        with (out/"tests.log").open("w") as log:
            result=subprocess.run([sys.executable,"-m","unittest","-v",*tests],cwd=source,
              env=dict(os.environ,PYTHONPATH=os.pathsep.join(sys.path[:2])),stdout=log,stderr=subprocess.STDOUT,timeout=120)
        engine.require(result.returncode==0,"Pure estimator tests failed")
        manifest,packet,runtime=engine.load_reference(out/"reference")
        psi=float(packet["parameters"].psi_child); q=float(packet["solution"].p_eq[0]); pension=float(packet["parameters"].pension); H=a.seed_horizon
        endpoint_receipt=Path(a.endpoint_receipt) if a.endpoint_receipt else DEFAULT_ENDPOINT
        fresh_endpoint=a.fresh_endpoint or not endpoint_receipt.exists()
        maximum_policy_calls=10*H+24+2+10+(4 if fresh_endpoint else 0)
        plan=estimator.draft_plan("four_successive");plan.update(source_pins=pins,target_contract=contract,horizons=[6,8])
        reused_seed_pin=None
        if a.jacobian_receipt is not None:
            reused_seed_pin=dict(path=str(a.jacobian_receipt.resolve()),sha256=engine.sha(a.jacobian_receipt))
        plan["acceleration"].update(seed_horizon=H,perturbed_date=a.perturbed_date,seed_receipt=reused_seed_pin)
        plan["budget"].update(total_seconds=a.total_seconds,candidate_seconds=a.total_seconds,
          mapping_seconds=a.mapping_seconds,jacobian_seconds=a.jacobian_seconds,path_seconds=a.jacobian_seconds,
          endpoint_seconds=a.jacobian_seconds,maximum_policy_calls=maximum_policy_calls)
        plan["endpoint"]["max_evaluations"]=2;plan["path"]["max_evaluations"]=2
        plan["path"]["cache_max_bytes"]=a.cache_max_bytes
        write("plan.json",dict(plan,smoke_maximum_policy_calls=maximum_policy_calls,fresh_endpoint=fresh_endpoint))
        write("best_so_far.json",dict(status="PENDING",phase="native_smoke_not_started",maximum_policy_calls=maximum_policy_calls))
        native=estimator.NativeEstimator(plan,out/"native",manifest,packet,runtime);native.deadline=started+a.total_seconds;native.candidate_deadline=native.deadline
        write("latest_completed.json",dict(phase="authenticated_jacobian" if reused_seed_pin else "measured_jacobian",epoch=time.time()))
        native.prepare_jacobian()
        seed=engine.pinned(reused_seed_pin) if reused_seed_pin else out/"native/jacobian_seed/measured/receipt.json"
        engine.require(seed.exists(),"Authenticated measured seed receipt missing")
        reused=False
        endpoint_receipt_pin=None
        if not fresh_endpoint:
            endpoint_receipt_pin=dict(path=str(endpoint_receipt),sha256=engine.sha(endpoint_receipt))
            endpoint=engine.read(engine.pinned(endpoint_receipt_pin))
            if a.endpoint_receipt is None:
                engine.require(endpoint["checkpoint"]["sha256"]==DEFAULT_ENDPOINT_CHECKPOINT_SHA256,"Default endpoint checkpoint SHA-256 changed")
            checkpoint=engine.pinned(endpoint["checkpoint"])
            with gzip.open(checkpoint,"rb") as stream: terminal=pickle.load(stream)
            engine.check_endpoint_primitives(packet["parameters"],terminal["parameters"])
            engine.require(np.array_equal(packet["b_grid"],terminal["b_grid"]) and endpoint["repeat_verified"] is True and
              endpoint["native_one_step_verified"] is True and endpoint["terminal"]["all_checks_pass"] is True and endpoint["housing"]=="fixed_stock" and
              endpoint["psi_child"]==terminal["parameters"].psi_child==psi and endpoint["price"]==terminal["solution"].p_eq[0]==q and
              endpoint["pension"]==terminal["parameters"].pension==pension and endpoint["reference_manifest_sha256"]==engine.MANIFEST_SHA and
              endpoint["source_manifest_sha256"]==manifest["source_manifest"]["sha256"],"Unchanged endpoint evidence failed identity/primitives/grid/terminal authentication")
            endpoint=dict(endpoint,checkpoint=dict(path=str(checkpoint),sha256=engine.sha(checkpoint)));reused=True
        else:
            terminal,endpoint=native.endpoint(psi)
        write("latest_completed.json",dict(phase="six_date_path",epoch=time.time(),reused_endpoint=reused))
        receipt,latest=native.path(psi,terminal,endpoint,6,out/"native/forecast")
        observed=float(receipt["fertility"][0]["period_tfr_topcode_adjusted"])
        native.latest=dict(psi=psi,latest=latest,terminal=terminal,endpoint=endpoint,summary={"payload":{"candidate":0,"models":[observed]}})
        synthetic=dict(converged=True,root={"final":dict(mapping_valid=True,prices=[psi],payload={"candidate":0})},parameter=dict(estimate=psi,lower=.01*psi,upper=2*psi,near_bound=False))
        native.advance(0,synthetic,diagnostics=False)
        from types import SimpleNamespace
        boundary=dict(evaluation=SimpleNamespace(policy=SimpleNamespace(V=latest["result"].values[6])))
        carried,record=native.guarded(a.mapping_seconds,lambda:engine.mapping(packet,runtime,boundary,endpoint,latest["prices"][1:],latest["pensions"][1:],np.full(5,psi),"fixed_stock",out/"native/carried",plan["path"]["cache_max_bytes"],initial_state=native.inherited,start_year=2011,measure_fertility=True))
        gaps=[]
        engine.require(len(record["rows"])==5 and len(receipt["rows"])==6,"Carried and full paths must have five and six rows")
        engine.require(len(record["fertility"])==5 and len(receipt["fertility"])==6,"Carried and full paths must have five and six fertility rows")
        for actual,expected in zip(record["rows"],receipt["rows"][1:]):
            engine.require(actual["period"]+1==expected["period"],"Wrong carried period offset")
            gaps += [abs(float(actual[k])-float(v))/max(1.,abs(float(v))) for k,v in expected.items() if k!="period" and isinstance(v,(int,float))]
        fertility_gaps=[abs(float(actual["period_tfr_topcode_adjusted"])-float(expected["period_tfr_topcode_adjusted"])) for actual,expected in zip(record["fertility"],receipt["fertility"][1:])]
        pf=runtime.rt["primitive"].pf; g=float(np.max(np.abs(carried.terminal_state.g_pre-latest["result"].terminal_state.g_pre)))
        queues={n:float(np.max(np.abs(pf.birth_queue_values(getattr(carried.terminal_state,n))-pf.birth_queue_values(getattr(latest["result"].terminal_state,n))))) for n in ("scheduled_entries","scheduled_raw_entries")}
        engine.require(all(record["gates"].values()) and max(gaps)<=2e-10 and max(fertility_gaps)<=1e-12 and g<=1e-12 and max(queues.values())<=1e-12,"Carried-state five-date remainder differs")
        write("readiness.json",dict(status="PASS",tests_passed=True,native_endpoint_and_fertility_smoke_passed=True,estimator_sources=pins,historical_fit=False,preference_changes=False,target_contract=contract,seed_receipt=dict(path=str(seed),sha256=engine.sha(seed)),reused_seed=bool(reused_seed_pin),reused_endpoint=dict(receipt=endpoint_receipt_pin,reused=reused,identity=endpoint),gates=dict(native_path=receipt["root_and_terminal_pass"],carried=all(record["gates"].values()),rows=max(gaps),fertility=max(fertility_gaps),g_pre=g,queues=queues),timings=dict(seconds=time.monotonic()-started),cache=record["cache"],cache_max_bytes=a.cache_max_bytes,maximum_policy_calls=maximum_policy_calls,execution_enabled=False,production_horizon_verified=False))
        write("best_so_far.json",dict(status="PASS",native_fertility=observed))
        write("latest_completed.json",dict(status="PASS",phase="complete",epoch=time.time(),reused_endpoint=reused))
    except BaseException as exc:
        write("failure.json",dict(status="FAIL",error_type=type(exc).__name__,error=str(exc),seconds=time.monotonic()-started,historical_fit=False,preference_changes=False,estimator_sources=pins))
        raise
    finally: stop.set();thread.join(timeout=1)

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--source-dir",type=Path,required=True);p.add_argument("--output",type=Path,required=True);p.add_argument("--empirical-blocks",type=Path,required=True);p.add_argument("--annual",type=Path,required=True);p.add_argument("--seed-horizon",type=int,default=10);p.add_argument("--perturbed-date",type=int,default=5);p.add_argument("--mapping-seconds",type=float,required=True);p.add_argument("--jacobian-seconds",type=float,required=True);p.add_argument("--total-seconds",type=float,required=True);p.add_argument("--cache-max-bytes",type=int,default=2*1024**3);p.add_argument("--endpoint-receipt",type=Path);p.add_argument("--jacobian-receipt",type=Path);p.add_argument("--fresh-endpoint",action="store_true");main(p.parse_args())
