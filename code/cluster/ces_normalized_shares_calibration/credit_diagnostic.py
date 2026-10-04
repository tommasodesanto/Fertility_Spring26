"""Matched fixed-price phi=.80/.95 diagnostic for the immutable CES v5 stage."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import time
import math

import numpy as np

PRICE = 0.6408361332017416
H0 = 7.306962620836552
CHAIN = Path("/work/deployment/followup_tools/credit_reference_chain1.json")
CHAIN_SHA256 = "26bd0fbb6901f9080d4e15472029c48847f117b21cb95f6d0de0be4cc3976e85"
THREADS = ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")


def write(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, default=_jsonable, allow_nan=False) + "\n")
    tmp.replace(path)


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return str(value)


def _public(P):
    values = json.loads(json.dumps(vars(P), default=_jsonable))
    return {k: v for k, v in values.items() if not k.startswith("_")}


def _summary_json(value):
    """Make undefined optional observer fields explicit without hiding their paths."""
    paths=[]
    def visit(node,path):
        if isinstance(node,np.ndarray):
            return visit(node.tolist(),path)
        if isinstance(node,dict):
            return {str(k):visit(v,f"{path}.{k}" if path else str(k)) for k,v in node.items()}
        if isinstance(node,(list,tuple)):
            return [visit(v,f"{path}[{i}]") for i,v in enumerate(node)]
        if isinstance(node,(float,np.floating)) and not math.isfinite(float(node)):
            paths.append(path)
            return None
        if isinstance(node,np.generic):
            return visit(node.item(),path)
        return node
    safe=visit(value,"")
    mandatory=("observed_fixed_price", "birth_flow_totals", "birth_flows_by_age",
               "utility_contract", "parameters", "price", "h0", "native_moments.tfr")
    fatal=[]
    for path in paths:
        low=path.lower()
        if low.startswith(mandatory):
            fatal.append(path)
    return safe,paths,fatal


def write_summary(path: Path, value):
    safe,paths,fatal=_summary_json(value)
    write(path.with_name("summary_nonfinite_fields.json"),
          dict(nonfinite_paths=paths, fatal_paths=fatal,
               note="Optional nonfinite observer fields are null in summary.json; raw solution arrays are retained."))
    if fatal:
        raise ValueError("Nonfinite mandatory closure, birth-flow or target-useful observer field(s): " + ", ".join(fatal))
    safe["nonfinite_paths"]=paths
    safe["observer_completeness"]="complete" if not paths else "missing_optional_observer_fields"
    safe["target_interpretation"]="descriptive_only" if not paths else "defer_until_missing_fields_are_reviewed"
    write(path,safe)


def _utility_receipt(P, sd):
    """Check the reached alpha and material multiplier arrays, cell by cell."""
    from experiments.ces_normalized_shares import adapter
    alpha = np.asarray(sd.alpha_flat).reshape(-1)
    multiplier = np.asarray(sd.escale_flat).reshape(-1)
    expected_a = np.empty((int(P.n_parity), int(P.n_child_states)))
    expected_m = np.empty_like(expected_a)
    for n in range(int(P.n_parity)):
        for cs in range(int(P.n_child_states)):
            m = adapter._children_at_home(P, n, cs)
            a = .733 if m == 0 else float(np.clip(.733-float(P.delta_alpha_jump)-float(P.delta_alpha)*m,.05,.95))
            e = ((2.+.7*m)/2.)**.7
            expected_a[n,cs] = a
            expected_m[n,cs] = (e*a**a*(1.-a)**(1.-a))**(float(P.sigma)-1.)
    expected_a=expected_a.reshape(-1,order="F")
    expected_m=expected_m.reshape(-1,order="F")
    if not np.allclose(alpha, expected_a, rtol=0., atol=1e-14):
        raise RuntimeError("Reached alpha array differs from normalized CES share contract")
    if not np.allclose(multiplier, expected_m, rtol=0., atol=1e-14):
        raise RuntimeError("Reached utility multiplier differs from normalized CES contract")
    return dict(alpha=alpha.tolist(), material_multiplier=multiplier.tolist(),
                delta_alpha_jump=float(P.delta_alpha_jump), delta_alpha=float(P.delta_alpha),
                alpha0=.733, normalized_ces_limit_shares=bool(P.normalized_ces_limit_shares),
                housing_floor=bool(P.child_room_floor), reference_rent=None)


def _standard_packet(context, live, destination):
    cal = context["prepared"].rt["primitive"].pf.calendar
    P, grid, sd, sol = (live[k] for k in ("P", "b_grid", "sd", "sol"))
    price = np.asarray(live["price"])
    policy = cal.policy_from_solution(sol, price, P, grid, sd)
    pre, _ = cal.reconstruct_stationary_pre_fertility(sol, policy, P, grid, sd)
    supply = cal.HousingSupplyRule("static-elastic", float(price[0]),
        float(P.H0[0]*(P.user_cost_rate*price[0]/P.r_bar[0])**P.xi_supply[0]), float(P.xi_supply[0]))
    ev = cal.evaluate_period(price, pre, P, grid, sd, cal.SolveCounter(),
                             supply_rule=supply, supplied_policy=policy)
    packet = dict(parameters=P, b_grid=grid, shared=sd, solution=sol,
                  evaluation=ev, stationary_g_pre=pre, supply_rule=supply,
                  demographic_seed=context["reference"].get("demographic_seed"))
    context["prepared"].rt["audit"].standard_diagnostics(packet, destination,
                                                           validate_production_young=False)
    names = sorted(p.name for p in (destination/"standard_diagnostics").glob("*.png"))
    if len(names) != 17 or names != sorted(context["manifest"]["standard_diagnostic_names"]):
        raise RuntimeError("Expected the exact native 17-plot packet")
    return ev, policy, names


def run(out: Path, prepare_only=False):
    if os.environ.get("CES_NORMALIZED_SHARES_STAGED_CONTEXT") != "1":
        raise RuntimeError("Requires the authenticated CES v5 staged container")
    bad = {k: os.environ.get(k) for k in THREADS if os.environ.get(k) != "1"}
    if bad:
        raise RuntimeError(f"Single CPU thread caps must all be 1: {bad}")
    source = json.loads(CHAIN.read_text())
    selected = source["selected"]
    actual_chain_sha=hashlib.sha256(CHAIN.read_bytes()).hexdigest()
    if actual_chain_sha != CHAIN_SHA256 or os.environ.get("CES_CREDIT_CHAIN_SHA256") != CHAIN_SHA256:
        raise RuntimeError("Chain-1 source SHA-256 pin missing or mismatched")
    theta = {k: float(v) for k, v in selected["parameters"].items()}
    if len(theta) != 11 or selected["price"] != PRICE or selected["H0_derived"] != H0:
        raise RuntimeError("Pinned chain-1 vector, price or derived H0 changed")
    from experiments.ces_normalized_shares import adapter
    from production import equilibrium, native_price, native_phase_b

    started = time.time()
    end = started + 1800.
    out.mkdir(parents=True, exist_ok=False)
    cases = []
    write(out/"begun.json", dict(status="running_fixed_price_diagnostic", started_epoch=started,
        deadline_epoch=end, source=str(CHAIN), source_sha256=actual_chain_sha,
        source_selected_fingerprint=selected["effective_input_fingerprint"], price=PRICE, H0=H0,
        fixed_parameters=theta, phi_values=[.8,.95], no_GE=True, no_recalibration=True))
    try:
        inputs = []
        for phi in (.8,.95):
            P, grid = adapter.load_inputs(theta)
            P.phi = np.full_like(P.phi, phi, dtype=float)
            P.H0 = np.full_like(P.H0, H0, dtype=float)
            inputs.append((P, grid))
        P0, grid0 = inputs[0]
        P1, grid1 = inputs[1]
        changed = sorted(k for k in set(_public(P0)) | set(_public(P1))
                         if _public(P0).get(k) != _public(P1).get(k))
        if changed != ["phi"] or not np.array_equal(grid0,grid1):
            raise RuntimeError(f"Paired public inputs differ beyond phi: {changed}")
        if not np.allclose(P0.H0, H0, rtol=0., atol=0.) or not np.allclose(P1.H0,H0,rtol=0.,atol=0.):
            raise RuntimeError("H0 must be held at chain-1 derived value in both arms")
        from production.engine import solver
        utility_preflight=[]
        with adapter.install():
            for phi,(P,grid) in zip((.8,.95),inputs):
                P.H0=np.full_like(P.H0,H0,dtype=float)
                sd=solver.precompute_shared(P,grid)
                utility_preflight.append(dict(phi=phi,**_utility_receipt(P,sd)))
        write(out/"prepared.json", dict(status="passed", changed_public_fields=changed,
            grid_sha256=hashlib.sha256(np.asarray(grid0,dtype=np.float64).tobytes()).hexdigest(),
            chain1_source_sha256=hashlib.sha256(CHAIN.read_bytes()).hexdigest(),
            price=PRICE, H0=H0, phi=[.8,.95], fixed_parameters=theta,
            utility_preflight=utility_preflight,
            preparation_only=prepare_only))
        if prepare_only:
            write(out/"completed.json", dict(status="preparation_passed_zero_solves"))
            return
        with adapter.install():
            for arm, (P, grid) in zip(("phi_080","phi_095"),inputs):
                case_start=time.time(); case_end=min(end,case_start+600.)
                if case_start >= end: raise TimeoutError("30-minute overall diagnostic budget exhausted")
                case=out/arm; case.mkdir()
                write(case/"begun.json",dict(status="started",phi=float(P.phi[0]),started_epoch=case_start,deadline_epoch=case_end))
                try:
                    context=equilibrium.build_context(P,grid,case,price_start=PRICE,
                        deadline=case_end,max_lifecycle=1,closure="fixed_h0")
                    budget=equilibrium.Budget(case,case_end,1)
                    live=native_price.solve_fixed_price(context,float(P.unsecured_credit_limit),PRICE,
                        budget, "selected_diagnostic_"+arm,case/"stage")
                    utility=_utility_receipt(live["P"],live["sd"])
                    observed=native_phase_b._observe_with_deadline(context,live,arm,final=False)
                    ev,policy,plots=_standard_packet(context,live,case/"diagnostic_packet")
                    Pactual=live["P"]
                    # Native distribution writes these exact event flows during evaluate_period.
                    flows={k:np.asarray(getattr(Pactual,k),dtype=float).tolist() for k in
                           ("_first_births_by_age","_second_births_by_age","_third_births_by_age")}
                    flow_totals={k.removeprefix("_").removesuffix("_by_age"):float(np.sum(v))
                                 for k,v in flows.items()}
                    moments=context["prepared"].rt["chain"].extract_moments(live["sol"],Pactual)
                    fertility=context["prepared"].rt["observe_initial_fertility"](
                        ev,Pactual,age_projection="uniform_birth_time")
                    housing=context["prepared"].rt["observe_initial_housing_wealth"](
                        ev,Pactual,grid,live["sd"],diagnostic_enabled=True,
                        age_projection="uniform_within_age_cell",diagnostic_allow_family_proxies=True,
                        include_wealth=True,include_birth_response=True)
                    arrays={k:v for k,v in vars(live["sol"]).items()
                            if isinstance(v,np.ndarray) and v.dtype!=object}
                    arrays.update({"shared."+k:v for k,v in vars(live["sd"]).items()
                                   if isinstance(v,np.ndarray) and v.dtype!=object})
                    np.savez_compressed(case/"solution_arrays.npz",**arrays)
                    write(case/"executed_P.json",vars(Pactual))
                    write_summary(case/"summary.json",dict(status="fixed_price_diagnostic_only",phi=float(Pactual.phi[0]),
                        price=PRICE,H0=float(Pactual.H0[0]),parameters={k:theta[k] for k in theta},
                        utility_contract=utility,observed_fixed_price=observed,native_moments=moments,
                        fertility=fertility,completed_childlessness_and_exactly_one={
                            k:v for k,v in fertility.items() if "childless" in k.lower() or "exactly_one" in k.lower()},
                        first_birth_age=next((v for k,v in fertility.items() if "first_birth_age" in k.lower()),None),
                        fertility_observer_fields=list(fertility),birth_flow_totals=flow_totals,birth_flows_by_age=flows,
                        rooms_and_ownership_by_young_and_30_55=housing,standard_plots=plots,
                        lifecycle_solves=budget.used_lifecycle,no_GE_certification=True))
                    row=dict(case=arm,summary=str(case/"summary.json"),plots=plots,elapsed_seconds=time.time()-case_start)
                    write(case/"completed.json",dict(status="completed",**row))
                    cases.append(row); write(out/"latest_completed.json",dict(cases=cases))
                except Exception as exc:
                    write(case/"failure.json",dict(status="failed_no_retry",error_type=type(exc).__name__,error=str(exc),time_epoch=time.time()))
                    raise
        write(out/"completed.json",dict(status="two_fixed_price_cases_completed",cases=cases,
            price=PRICE,H0=H0,only_economic_change="uniform phi .80 to .95",no_GE=True,no_recalibration=True))
    except Exception as exc:
        write(out/"failure.json",dict(status="stopped_no_retry",error_type=type(exc).__name__,error=str(exc),completed_cases=cases,time_epoch=time.time()))
        raise


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=Path("/work/results/credit_diagnostic"))
    parser.add_argument("--prepare-only",action="store_true")
    args=parser.parse_args()
    run(args.out,args.prepare_only)
