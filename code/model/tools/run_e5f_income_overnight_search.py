"""Finite, provenance-pinned multivariate overnight income search controller.

The native evaluator and score gates remain owned by the existing income
candidate helper.  This file only constructs a deterministic bounded proposal
set, dispatches it, and records reviewable receipts.
"""
from __future__ import annotations

import argparse, csv, hashlib, json, math, os, signal, subprocess, sys, time, random, shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Callable

PARAMETERS = ("beta_annual", "kappa_fert", "kappa_fert_continuation", "chi", "H0",
              "theta0", "theta1", "first_birth_fixed_cost", "h_P")
POSITIVE = {"kappa_fert", "kappa_fert_continuation", "chi", "H0", "theta1", "h_P"}
OBJECTIVE = "4440ea07f4de957740ca6c04961d2806d9b9ef782c7a0e7dad4ce73e1db651b1"
FIXED_SEED = 20260919
STOP_DISPATCH_SECONDS = 900.0

from run_e5f_income_candidate_search import ContractError, STOP

def read(path: Path): return json.loads(Path(path).read_text())
def write(path: Path, value: Any):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(path) + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    tmp.replace(path)
def fingerprint(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()
def sha(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""): h.update(block)
    return h.hexdigest()

def old_helper():
    import run_e5f_income_candidate_search as old
    return old

def validate_score_contract(receipt: dict, plan: dict):
    old_helper().validate_score_contract(receipt, plan)
def validate_parameters(receipt: dict, expected: dict[str, float]):
    old_helper().validate_parameters(receipt, expected)
def extract_score(receipt: dict) -> float:
    return old_helper().extract_score(receipt)

def _bounds(plan):
    bounds = plan.get("parameter_bounds") or plan.get("pilot_box")
    if not bounds or set(bounds) != set(PARAMETERS): raise ContractError("complete parameter bounds required")
    return {n: tuple(map(float, bounds[n])) for n in PARAMETERS}

def _clip(point, bounds):
    return {n: min(bounds[n][1], max(bounds[n][0], float(point[n]))) for n in PARAMETERS}

def _seed_from(row):
    if not isinstance(row, dict): return None
    point = row.get("parameters") or row.get("initial_psi_parameters") or row.get("proposal")
    if isinstance(point, dict) and all(n in point for n in PARAMETERS):
        return {n: float(point[n]) for n in PARAMETERS}
    return None

def load_seeds(plan: dict) -> list[dict[str, Any]]:
    """Load pilot plus four pinned prior seeds; reject absent/malformed provenance."""
    pilot = plan.get("pilot_parameters")
    if not isinstance(pilot, dict) or set(pilot) != set(PARAMETERS): raise ContractError("pilot parameter pin missing")
    seeds = [{"label": "originalseed", "parameters": {n: float(pilot[n]) for n in PARAMETERS}}]
    rows = plan.get("prior_successful_cases") or plan.get("seed_cases") or []
    if isinstance(rows, str): rows = read(Path(rows)).get("cases", [])
    for row in rows:
        p = _seed_from(row)
        if p is None: continue
        label = str(row.get("label") or row.get("case_id") or f"prior_{len(seeds)}")
        seeds.append({"label": label, "parameters": p, "source": row.get("source")})
    if len(seeds) < 5:
        raise ContractError("plan must pin original seed plus four prior successful low-loss seeds")
    # Distinct and bounded; provenance is retained verbatim in each proposal.
    bounds = _bounds(plan); out=[]; seen=set()
    for item in seeds[:5]:
        p = _clip(item["parameters"], bounds); key=fingerprint(p)
        if key in seen: raise ContractError("duplicate pinned search seed")
        seen.add(key); out.append({**item, "parameters": p})
    return out

def prepare_plan(plan: dict, prior_summary: Path | None) -> dict:
    if prior_summary is None: return plan
    summary = read(prior_summary)
    if summary.get("status") != "verified_selection" or summary.get("verification",{}).get("status") != "verified":
        raise ContractError("prior summary lacks verified selection")
    rows=[]; seen={fingerprint(plan["pilot_parameters"])}; known=[plan["pilot_parameters"]]
    candidates=sorted([x for x in summary.get("cases",[]) if x.get("status")=="completed"],key=lambda x:(x["objective"],x["case"]))
    for row in candidates:
        receipt=old_helper().load_verified_score(Path(row["receipt"]["_evaluation"])/"scored_repetition_01/score.json",plan)
        validate_parameters(receipt,row["parameters"])
        if receipt["loss"] != row["objective"]: raise ContractError("prior row loss mismatch")
        point=_seed_from(row); known.append(point); key=fingerprint(point)
        if key in seen: continue
        seen.add(key)
        if len(rows)<4: rows.append({"label":f"prior_case_{row['case']}","parameters":point,"source":str(prior_summary)})
    if len(rows)<4: raise ContractError("fewer than four distinct successful seeds")
    return {**plan,"prior_successful_cases":rows,"known_parameters":known,
            "prior_summary_sha256":sha(prior_summary),"prior_summary":str(prior_summary)}

def proposals(plan: dict, max_proposals: int = 96) -> list[dict[str, Any]]:
    """Two deterministic log-local scales around five pinned seeds (<=96 points).

    Every proposal moves all nine coordinates.  Sign vectors are fixed and
    antithetic, so this is a local multivariate experiment rather than a
    collection of one-coordinate probes.
    """
    bounds = _bounds(plan); seeds = load_seeds(plan); rows=[]; seen=set()
    rng = random.Random(FIXED_SEED)
    vectors = [tuple(rng.uniform(-1.,1.) for _ in PARAMETERS) for _ in range(6)]
    vectors = [u for v in vectors for u in (v,tuple(-x for x in v))]
    widths = {"beta_annual":.012,"kappa_fert":.55,"kappa_fert_continuation":.55,"chi":.20,
              "H0":.35,"theta0":.20,"theta1":.55,"first_birth_fixed_cost":.20,"h_P":.30}
    seen={fingerprint(x) for x in plan.get("known_parameters",[])} | {fingerprint(x["parameters"]) for x in seeds}
    for j,vector in enumerate(vectors):
        for seed in seeds:
            for scale_label,scale in (("near",.4),("wide",1.0)):
                base=seed["parameters"]; point={}
                for name,u in zip(PARAMETERS,vector):
                    point[name]=base[name]*math.exp(scale*widths[name]*u) if name in POSITIVE else base[name]+scale*widths[name]*u
                point=_clip(point,bounds); key=fingerprint(point)
                if key in seen: continue
                seen.add(key);rows.append({"parameters":point,"seed_label":seed["label"],"scale":scale_label,"direction":j,"proposal_seed":FIXED_SEED})
                if len(rows)>=max_proposals:return rows
    return rows

def _numeric_fit(r):
    return {"loss": r.get("loss"), "price": r.get("_native_price"),
            "parameters": [(x.get("parameter"), x.get("estimate")) for x in r.get("parameters", [])],
            "target_fit": [(x.get("restriction_id", x.get("target")), x.get("target"), x.get("model"),
                            x.get("gap"), x.get("actual_weight"), x.get("loss_contribution")) for x in r.get("target_fit", [])],
            "normalization": (r.get("normalization", {}).get("completed_fertility"),
                              r.get("normalization", {}).get("psi_child"))}

def adapter_evaluator(plan_path: Path, output: Path, parameters: dict[str, float], *, psi: float,
                      repetitions: int, timeout: float, case_id: str) -> dict:
    """Use the existing evaluator, including its heartbeat and psutil cleanup."""
    return old_helper().adapter_evaluator(plan_path, output, parameters, psi=psi,
                                          repetitions=repetitions, timeout=timeout,
                                          case_id=case_id)

def run_search(plan: dict, output: Path, evaluator: Callable[..., dict], *, incumbent: dict,
               search_seconds: float = 16200, max_proposals: int = 96, workers: int = 8,
               case_timeout: float = 900, verification_seconds: float = 1800) -> dict:
    try: validate_score_contract(incumbent, plan)
    except ContractError:
        STOP.set(); raise
    seeds = load_seeds(plan); all_props = proposals(plan, max_proposals)
    actual = {row["parameter"]: float(row["estimate"]) for row in incumbent.get("parameters", [])
              if row.get("structural_coordinate")}
    if set(actual) != set(PARAMETERS): raise ContractError("incumbent lacks all structural coordinates")
    output.mkdir(parents=True, exist_ok=False); write(output / "incumbent_score.json", incumbent)
    best = {"status":"incumbent", "objective":extract_score(incumbent), "parameters":actual,
            "psi": float(incumbent.get("_initial_psi", plan.get("pilot_initial_psi", 0.1489153145785918))),
            "seed_label":"originalseed", "receipt":incumbent}
    write(output / "best_so_far.json", best); write(output / "latest_completed.json", {"status":"awaiting_first_case"}); completed=[]; failed_batches=0; started=time.monotonic()
    for batch_start in range(0, len(all_props), workers):
        remaining = search_seconds - (time.monotonic() - started)
        if remaining < STOP_DISPATCH_SECONDS: write(output/"stop.json", {"status":"stop_dispatch_budget","completed":len(completed),"remaining":remaining}); break
        batch = all_props[batch_start:batch_start+workers]
        def one(item):
            idx, proposal = item; case=output/"cases"/f"case_{idx:03d}"
            try:
                receipt=evaluator(proposal["parameters"],case,psi=float(plan.get("pilot_initial_psi",0.1489153145785918)),repetitions=1,
                                  timeout=min(case_timeout,remaining),case_id=f"income_overnight_v1_{idx:03d}")
                validate_score_contract(receipt,plan); validate_parameters(receipt,proposal["parameters"])
                return {"case":idx,"status":"completed","objective":extract_score(receipt),**proposal,"receipt":receipt}
            except Exception as exc:
                # Old helper's ContractError is intentionally fatal; numerical
                # and solver failures are reviewable case rows.
                if isinstance(exc,ContractError):
                    STOP.set()
                    raise
                return {"case":idx,"status":"numerical_failure","error":str(exc),**proposal}
        with ThreadPoolExecutor(max_workers=min(workers,len(batch))) as pool:
            futures=[pool.submit(one,(batch_start+j+1,p)) for j,p in enumerate(batch)]
            batch_rows=[]
            for fut in as_completed(futures):
                row=fut.result(); batch_rows.append(row); completed.append(row); write(output/"latest_completed.json",row); write(output/"cases.json",completed)
                if row["status"]=="completed" and (row["objective"],row["case"]) < (best["objective"],best.get("case",0)):
                    best={"status":"improved","objective":row["objective"],"parameters":row["parameters"],"psi":float(plan.get("pilot_initial_psi",0.1489153145785918)),"seed_label":row["seed_label"],"receipt":row["receipt"],"case":row["case"]}; write(output/"best_so_far.json",best)
        write(output/"cases.json",completed)
        if not any(x["status"]=="completed" for x in batch_rows): failed_batches += 1
        else: failed_batches=0
        if failed_batches >= 6: write(output/"stop.json", {"status":"six_failed_batches","completed":len(completed)}); break
    verification={"status":"failed"}
    try:
        check=evaluator(best["parameters"],output/"selected_verification",psi=best["psi"],repetitions=2,timeout=verification_seconds,case_id="income_overnight_v1_selected_repetition")
        validate_score_contract(check,plan); validate_parameters(check,best["parameters"])
        exact=extract_score(check)==float(best["objective"]); summary=check.get("_summary",{})
        numeric=_numeric_fit(check)==_numeric_fit(best["receipt"])
        verification={"status":"verified" if exact and summary.get("repetitions")==2 and summary.get("exact_loss_equality") is True and numeric else "mismatch","exact_objective":exact,"numeric_fit_equal":numeric,"receipt":check}
    except Exception as exc: verification={"status":"failed","error":str(exc)}
    chosen=verification.get("receipt",best.get("receipt"));
    if chosen:
        for name, rows in (("selected_target_fit",chosen.get("target_fit",[])),("selected_parameters",chosen.get("parameters",[]))):
            fields=list(dict.fromkeys(k for r in rows for k in r));
            with (output/(name+".csv")).open("w",newline="") as f: csv.DictWriter(f,fieldnames=fields).writeheader(); csv.DictWriter(f,fieldnames=fields).writerows(rows)
    if chosen and chosen.get("_evaluation"):
        rep="repetition_02" if verification["status"]=="verified" else "repetition_01"
        source=Path(chosen["_evaluation"])/"raw"/rep/"standard_diagnostics"
        if len(list(source.glob("*.png")))!=17:raise ContractError("missing 17 selected diagnostics")
        shutil.copytree(source,output/"selected_standard_diagnostics")
    proposal_success = any(x.get("status")=="completed" for x in completed)
    verified = verification["status"]=="verified" and proposal_success
    result={"schema":"e5f_income_overnight_v1","status":"verified_selection" if verified else "requires_review","selected":best,"verification":verification,"cases":completed,"proposal_count":len(completed),"proposal_success":proposal_success,"max_proposals":max_proposals,"poll_complete":len(completed)==len(all_props),"objective_canonical_sha256":OBJECTIVE,"incumbent_objective":extract_score(incumbent),"fixed_seed":FIXED_SEED,"old_helper_controller_sha256":plan.get("controller_sha256"),"overnight_controller_sha256":sha(Path(__file__))}
    write(output/"summary.json",result); return result

def main():
    p=argparse.ArgumentParser();p.add_argument("--mode",choices=("search","smoke"),required=True)
    p.add_argument("--plan",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
    p.add_argument("--prior-summary",type=Path,required=True);p.add_argument("--smoke-summary",type=Path)
    a=p.parse_args();plan=read(a.plan)
    from run_e5f_income_candidate_calibration import validate_plan
    validate_plan(plan,require_source=True)
    if sha(Path(__file__))!=plan["overnight_controller_sha256"]:raise ContractError("overnight controller hash mismatch")
    if sha(Path(old_helper().__file__))!=plan["controller_sha256"]:raise ContractError("helper hash mismatch")
    contract=read(Path(plan["source_manifest_path"]))
    for r in contract["parameter_restrictions"]:
        expected=[r["lower"],.99 if r["parameter"]=="beta_annual" else r["upper"]]
        if plan["parameter_bounds"][r["parameter"]]!=expected:raise ContractError("parameter bounds mismatch")
    def stop(signum,frame):
        STOP.set();raise KeyboardInterrupt("overnight terminated")
    signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
    plan=prepare_plan(plan,a.prior_summary)
    prior=read(a.prior_summary)
    def verified_selected(summary):
        if summary.get("status")!="verified_selection" or summary.get("verification",{}).get("status")!="verified":raise ContractError("unverified dependency")
        r=summary["verification"]["receipt"]
        actual=old_helper().load_verified_score(Path(r["_evaluation"])/"scored_repetition_01/score.json",plan)
        validate_parameters(actual,summary["selected"]["parameters"])
        if _numeric_fit(r)!=_numeric_fit(actual):raise ContractError("dependency score changed")
        return actual
    incumbent=verified_selected(prior)
    if a.smoke_summary:
        sm=read(a.smoke_summary);candidate=verified_selected(sm)
        if not sm.get("proposal_success"):raise ContractError("smoke had no successful new proposals")
        plan["smoke_summary_sha256"]=sha(a.smoke_summary)
        plan["known_parameters"] += [r["parameters"] for r in sm["cases"]]
        if candidate["loss"]<incumbent["loss"]:incumbent=candidate
    elif a.mode=="search":raise ContractError("production requires verified smoke")
    # Keep the original pilot as an exploratory center; warm normalization only.
    plan["pilot_initial_psi"]=float(incumbent["normalization"]["psi_child"])
    runtime=a.output.parent/(a.output.name+".plan.json");write(runtime,plan)
    result=run_search(plan,a.output,lambda params,case,**kw:adapter_evaluator(runtime,case,params,**kw),
        incumbent=incumbent,max_proposals=2 if a.mode=="smoke" else 96,
        workers=2 if a.mode=="smoke" else 8,search_seconds=1800 if a.mode=="smoke" else 16200,
        verification_seconds=1500 if a.mode=="smoke" else 1800)
    if result["status"]!="verified_selection":raise SystemExit(2)
if __name__=="__main__":main()
