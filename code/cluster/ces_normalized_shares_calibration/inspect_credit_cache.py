"""Postprocess a saved CES credit diagnostic after a rejected native gate.

This reader performs no lifecycle solve and never calls the production gate
stack. It reconstructs reports from the saved native arrays and retains the
failed estate audit as a rejected result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np

PRICE = 0.6408361332017416
H0 = 7.306962620836552
CHAIN = Path("/work/deployment/followup_tools/credit_reference_chain1.json")
CHAIN_SHA256 = "26bd0fbb6901f9080d4e15472029c48847f117b21cb95f6d0de0be4cc3976e85"
EXPECTED_NATIVE_FAILURE = "Negative-estate production gate failed; saved ledger retained"


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Not JSON serializable: {type(value).__name__}")


def write(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    safe, paths = _sanitize(value)
    safe["nonfinite_paths"] = paths
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(safe, indent=2, sort_keys=True, default=_jsonable,
                              allow_nan=False) + "\n")
    tmp.replace(path)


def _sanitize(value, path=""):
    paths=[]
    def visit(node, current):
        if isinstance(node, np.ndarray):
            return visit(node.tolist(), current)
        if isinstance(node, dict):
            out={}
            for key,item in node.items():
                child=f"{current}.{key}" if current else str(key)
                out[str(key)]=visit(item,child)
            return out
        if isinstance(node,(list,tuple)):
            return [visit(item,f"{current}[{i}]") for i,item in enumerate(node)]
        if isinstance(node,(float,np.floating)) and not math.isfinite(float(node)):
            paths.append(current)
            return None
        if isinstance(node,np.generic):
            return visit(node.item(),current)
        return node
    result=visit(value,path)
    # Gather paths once across the recursive traversal.
    def collect(node,current=""):
        found=[]
        if isinstance(node,dict):
            for key,item in node.items(): found.extend(collect(item,f"{current}.{key}" if current else str(key)))
        elif isinstance(node,(list,tuple)):
            for i,item in enumerate(node): found.extend(collect(item,f"{current}[{i}]"))
        elif isinstance(node,(float,np.floating)) and not math.isfinite(float(node)):
            found.append(current)
        return found
    return result,collect(value,path)


def sha256(path: Path):
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(1<<20),b""):
            h.update(block)
    return h.hexdigest()


def compare_named(control, case):
    """Compare only explicit observer moment names and rows keyed by moment."""
    def moments(left,right):
        names=sorted(set(left or {})|set(right or {}))
        out={}
        for name in names:
            a=(left or {}).get(name); b=(right or {}).get(name)
            row={"control":a,"phi_095":b}
            if isinstance(a,(int,float)) and isinstance(b,(int,float)):
                if math.isfinite(float(a)) and math.isfinite(float(b)):
                    row["change"]=float(b)-float(a)
            out[name]=row
        return out
    def rows(left,right):
        l={str(x["moment"]):x for x in (left or []) if isinstance(x,dict) and "moment" in x}
        r={str(x["moment"]):x for x in (right or []) if isinstance(x,dict) and "moment" in x}
        return moments({k:v.get("model_value") for k,v in l.items()},
                       {k:v.get("model_value") for k,v in r.items()})
    return {
        "fertility_moments":moments(control.get("fertility",{}).get("moments"),
                                     case.get("fertility",{}).get("moments")),
        "housing_moments":moments(control.get("rooms_and_ownership_by_young_and_30_55",{}).get("moments"),
                                   case.get("housing_wealth",{}).get("moments")),
        "housing_rows_by_moment":rows(control.get("rooms_and_ownership_by_young_and_30_55",{}).get("rows"),
                                       case.get("housing_wealth",{}).get("rows")),
        "control_status":control.get("status"),
        "phi_095_status":"native_estate_gate_rejected",
    }


def run(input_root: Path, out: Path):
    if not CHAIN.is_file() or sha256(CHAIN)!=CHAIN_SHA256:
        raise RuntimeError("Pinned chain-1 source missing or SHA-256 mismatch")
    arm=input_root/"phi_095"
    stage=arm/"stage"
    arrays_path=stage/"solution_arrays.npz"
    stage_summary=json.loads((stage/"summary.json").read_text())
    failure=json.loads((arm/"failure.json").read_text())
    if failure.get("error")!=EXPECTED_NATIVE_FAILURE:
        raise RuntimeError("Cached arm is not the expected negative-estate-gate rejection")
    control_path=input_root/"phi_080"/"summary.json"
    control=json.loads(control_path.read_text())
    if stage_summary.get("label")!="selected_diagnostic_phi_095":
        raise RuntimeError("Cached native stage label mismatch")
    if float(stage_summary["q"])!=PRICE or not arrays_path.is_file():
        raise RuntimeError("Cached price or native arrays missing/mismatched")
    source=json.loads(CHAIN.read_text())
    selected=source["selected"]
    theta={k:float(v) for k,v in selected["parameters"].items()}

    from experiments.ces_normalized_shares import adapter
    from production import equilibrium
    from production.credit import bind_engine_credit
    from production.engine import solver

    started=time.time()
    out.mkdir(parents=True,exist_ok=False)
    write(out/"begun.json",dict(status="cached_postprocessing_started",started_epoch=started,
        no_lifecycle_solve=True,no_native_gate_call=True,native_gate_status="rejected",
        native_gate_error=failure["error"],input_arrays_sha256=sha256(arrays_path),
        control_summary_sha256=sha256(control_path),stage_summary=stage_summary))
    try:
        P,grid=adapter.load_inputs(theta)
        P.phi=np.full_like(P.phi,.95,dtype=float)
        P.H0=np.full_like(P.H0,H0,dtype=float)
        bind_engine_credit(P,"corrected",float(P.unsecured_credit_limit))
        with np.load(arrays_path,allow_pickle=False) as cached:
            sol_arrays={k:cached[k].copy() for k in cached.files if not k.startswith("shared.")}
            shared_arrays={k[len("shared."):]:cached[k].copy() for k in cached.files if k.startswith("shared.")}
        sol=SimpleNamespace(**sol_arrays)
        P._fert2_probs=sol.fert2_probs.copy()
        out_context=out/"context"
        with adapter.install():
            context=equilibrium.build_context(P,grid,out_context,price_start=PRICE,
                deadline=time.time()+280.,max_lifecycle=1,closure="fixed_h0")
            shared=solver.precompute_shared(P,grid)
            mismatches=[]
            for name,cached_value in shared_arrays.items():
                if not hasattr(shared,name) or not np.array_equal(np.asarray(getattr(shared,name)),cached_value):
                    mismatches.append(name)
            if mismatches:
                raise RuntimeError("Cached native shared inputs differ on reconstruction: "+", ".join(mismatches))
            cal=context["prepared"].rt["primitive"].pf.calendar
            price=np.asarray([PRICE],dtype=float)
            policy=cal.policy_from_solution(sol,price,P,grid,shared)
            pre,reconstruction=cal.reconstruct_stationary_pre_fertility(sol,policy,P,grid,shared)
            supply=cal.HousingSupplyRule("static-elastic",PRICE,
                float(P.H0[0]*(P.user_cost_rate*PRICE/P.r_bar[0])**P.xi_supply[0]),float(P.xi_supply[0]))
            ev=cal.evaluate_period(price,pre,P,grid,shared,cal.SolveCounter(),
                supply_rule=supply,supplied_policy=policy)

            # Audit and preserve the rejection; deliberately do not invoke the
            # native gate or downgrade its status.
            estate=context["prepared"].estate.audit(ev,P,grid)
            write(out/"estate_ledger.json",dict(status="diagnostic_ledger_native_gate_rejected",
                native_failure=failure,ledger=estate,source_stage_summary=stage_summary))
            fertility=context["prepared"].rt["observe_initial_fertility"](
                ev,P,age_projection="uniform_birth_time")
            accounting=fertility["accounting"]
            flow_matrix=np.asarray(accounting["parity_birth_flows_by_age"],dtype=float)
            if flow_matrix.shape!=(int(P.J),3):
                raise RuntimeError(f"Native parity birth-flow observer has unexpected shape {flow_matrix.shape}")
            names=("first_births","second_births","third_births")
            flow_totals={name:float(flow_matrix[:,i].sum()) for i,name in enumerate(names)}
            flows={name:flow_matrix[:,i].tolist() for i,name in enumerate(names)}
            if abs(flow_totals["first_births"]-float(accounting["first_birth_flow"]))>1e-10:
                raise RuntimeError("Native first-birth flow identity fails")
            if abs(sum(flow_totals.values())-float(accounting["explicit_birth_flow"]))>1e-10:
                raise RuntimeError("Native explicit birth-flow identity fails")
            control_flows=control["birth_flow_totals"]
            control_matrix=np.asarray(control["fertility"]["accounting"]["parity_birth_flows_by_age"],dtype=float)
            control_matrix_totals={name:float(control_matrix[:,i].sum()) for i,name in enumerate(names)}
            control_flow_gaps={name:control_matrix_totals[name]-float(control_flows[name]) for name in names}
            if max(abs(v) for v in control_flow_gaps.values())>1e-10:
                raise RuntimeError(f"Baseline observer flows differ from saved native flows: {control_flow_gaps}")
            flow_gaps={name:flow_totals[name]-float(control_flows[name]) for name in names}
            housing=context["prepared"].rt["observe_initial_housing_wealth"](
                ev,P,grid,shared,diagnostic_enabled=True,age_projection="uniform_within_age_cell",
                diagnostic_allow_family_proxies=True,include_wealth=True,include_birth_response=True)
            packet=dict(parameters=P,b_grid=grid,shared=shared,solution=sol,evaluation=ev,
                stationary_g_pre=pre,supply_rule=supply,
                demographic_seed=context["reference"].get("demographic_seed"))
            plots_dir=out/"standard_diagnostics"
            plot_failure=None
            try:
                context["prepared"].rt["audit"].standard_diagnostics(
                    packet,out,validate_production_young=False)
                plots=sorted(p.name for p in plots_dir.glob("*.png"))
                if plots!=sorted(context["manifest"]["standard_diagnostic_names"]):
                    raise RuntimeError(f"Expected native 17 plots; found {len(plots)}")
            except Exception as exc:
                plots=[]
                plot_failure=dict(error_type=type(exc).__name__,error=str(exc),
                    missing_solution_statistics="see exact error; cached scalar statistics are not imputed")
            case=dict(fertility=fertility,housing_wealth=housing)
            control_tfr=control.get("native_moments",{}).get("tfr")
            control_own=control.get("native_moments",{}).get("aggregate_own_rate")
            paired_core={
                "native_tfr":{"phi_080":control_tfr,"phi_095":stage_summary.get("tfr"),
                    "change":(float(stage_summary["tfr"])-float(control_tfr)) if control_tfr is not None else None},
                "aggregate_ownership":{"phi_080":control_own,"phi_095":stage_summary.get("own_rate"),
                    "change":(float(stage_summary["own_rate"])-float(control_own)) if control_own is not None else None},
                "birth_flow_change":{"phi_080":control_flows,"phi_095":flow_totals,"change":flow_gaps},
            }
            write(out/"cached_observations.json",dict(status="descriptive_cached_arrays_native_gate_rejected",
                native_gate="rejected",native_gate_error=failure["error"],phi=.95,price=PRICE,H0=H0,
                source_parameters=theta,flow_totals=flow_totals,flows_by_age=flows,
                paired_core=paired_core,estate_ledger=estate,
                observer_birth_accounting={k:accounting[k] for k in ("first_birth_flow","explicit_birth_flow",
                    "maximum_parity_flow_error")},
                fertility=fertility,housing_wealth=housing,
                reconstruction_checks=reconstruction,plot_names=plots,plot_failure=plot_failure,
                control_comparison=compare_named(control,case),no_lifecycle_solve=True,
                no_gate_waiver=True,interpretation="diagnostic only; rejected native estate gate"))
            if plot_failure:
                write(out/"plots_incomplete.json",plot_failure)
            write(out/"completed.json",dict(status="cached_postprocessing_completed_with_rejected_native_gate",
                lifecycle_solves=0,estate_ledger=str(out/"estate_ledger.json"),
                observations=str(out/"cached_observations.json"),standard_plot_count=len(plots),
                native_gate="rejected",elapsed_seconds=time.time()-started))
    except Exception as exc:
        write(out/"failure.json",dict(status="cached_postprocessing_failed_no_retry",
            error_type=type(exc).__name__,error=str(exc),lifecycle_solves=0,
            native_gate="rejected",time_epoch=time.time()))
        raise


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root",type=Path,default=Path("/work/results/credit_diagnostic"))
    parser.add_argument("--out",type=Path,default=Path("/work/results/cached_inspection_v2"))
    args=parser.parse_args()
    run(args.input_root,args.out)
