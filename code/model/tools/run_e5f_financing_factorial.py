"""Bounded native financing x rental-cap factorial diagnostic.

This controller composes the existing native fixed-price primitives.  It never
changes prices, preferences, entry rules, KFE masses, or the target contract.
Families are separate checkpoint contracts; comparisons across families are
descriptive because their income processes and entry distributions can differ.
"""
from __future__ import annotations
import argparse, copy, gzip, hashlib, importlib, json, pickle, sys, time, subprocess, signal, csv
from pathlib import Path
from typing import Any, Mapping
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
POLICY_NAMES = ("V", "c_pol", "hR_pol", "bp_pol", "tenure_choice", "tenure_probs", "loc_probs", "fert_probs", "fert_value", "fert2_probs", "price")
PHIS = (0.8, 1.0)
LAMBDAS = (0.0, 5.0)
CAPS = (6.0, 10.0)
FAMILY_NAMES = ("original", "stationary_new_income", "refit_new_income")

def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

def write(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    tmp.replace(path)

def packet(path: Path) -> Any:
    with gzip.open(path, "rb") as f: return pickle.load(f)

def install(source: Path) -> None:
    for p in (source, source / "tools", ROOT / "code/model/tools"):
        if str(p.resolve()) not in sys.path: sys.path.insert(0, str(p.resolve()))

def validate(args: argparse.Namespace) -> dict[str, Any]:
    from run_e5f_income_candidate_calibration import validate_plan, read
    plan=read(args.plan);validate_plan(plan,require_source=True)
    if sha(Path(__file__))!=plan["factorial_controller_sha256"]:raise ValueError("factorial controller hash mismatch")
    for name,pin in plan["factorial_helper_sha256"].items():
        if sha(Path(__file__).with_name(name))!=pin:raise ValueError("factorial helper hash mismatch: "+name)
    args.source_root=Path(plan["source_root"])/"code/model"
    install(args.source_root)
    if args.selection_summary:
        if args.family!="refit_new_income":raise ValueError("selection only allowed for refit family")
        from run_e5f_income_candidate_search import load_verified_score,validate_parameters
        summary=read(args.selection_summary)
        if summary.get("status")!="verified_selection" or summary.get("verification",{}).get("status")!="verified":raise ValueError("refit selection unverified")
        evaluation=Path(summary["verification"]["receipt"]["_evaluation"])
        score=load_verified_score(evaluation/"scored_repetition_01/score.json",plan)
        validate_parameters(score,summary["selected"]["parameters"])
        if score["loss"]!=summary["selected"]["objective"]:raise ValueError("selection score mismatch")
        args.checkpoint=evaluation/"raw/repetition_02/initial_state.pkl.gz"
        args.expected_hash=read(args.checkpoint.parent/"summary.json")["checkpoint_sha256"]
    if not args.checkpoint or not args.expected_hash:raise ValueError("checkpoint and expected hash required")
    actual=sha(args.checkpoint)
    if actual!=args.expected_hash:raise ValueError("checkpoint hash mismatch")
    x=packet(args.checkpoint);P=x["parameters"]
    if args.family=="original":
        if actual!="3322a61994fb3654d67f4b1d6cf2d0f7cacbb3668d06a417e192ee363c174993":raise ValueError("wrong original checkpoint")
    else:
        if getattr(P,"income_candidate_fingerprint",None)!=plan["candidate_payload_fingerprint"] or P.permanent_income_levels_enabled:raise ValueError("wrong candidate income")
        if args.family=="stationary_new_income" and actual!=plan["incumbent_checkpoint_sha256"]:raise ValueError("wrong stationary pilot checkpoint")
    return {"family":args.family,"checkpoint":str(args.checkpoint),"checkpoint_sha256":actual,
            "source_manifest":str(plan["source_manifest_path"]),"source_root":plan["source_root"],
            "candidate_fingerprint":getattr(P,"income_candidate_fingerprint",None),
            "price":np.asarray(x["evaluation"].policy.price).tolist(),
            "initial_population_shape":list(np.asarray(x["stationary_g_pre"]).shape),
            "scope":"fixed-price partial equilibrium; within-family preferences and initial population fixed",
            "closure":"checkpoint population, fiscal objects, geography, prices and entry law held fixed; no GE claim"}

def arm_parameters(base: Any, phi: float, lam: float, cap: float) -> Any:
    P = copy.deepcopy(base)
    P.phi = np.full_like(np.asarray(base.phi, dtype=float), phi)
    P.lambda_d = float(lam)
    if lam:
        P.debt_taper_start_age, P.debt_taper_end_age = 82.0, 86.0
    P.hR_max = float(cap)
    importlib.import_module("intergen_eqscale_seq_optimized.parameters").build_debt_caps(P)
    native=importlib.import_module("run_e5f_native_financing_diagnostic")
    altered=native.changed(base,P)
    allowed={"phi","lambda_d","hR_max","debt_taper_start_age","debt_taper_end_age","debt_taper_weights","debt_caps"}
    if set(altered)-allowed:raise ValueError("unexpected parameter changes: "+str(set(altered)-allowed))
    if phi==.8 and lam==0 and cap==6 and altered:raise ValueError("baseline changed parameters")
    return P

def policy_arrays(policy: Any) -> dict[str, np.ndarray]:
    out = {n: np.asarray(getattr(policy, n, None)) for n in POLICY_NAMES}
    missing = [n for n in POLICY_NAMES if getattr(policy, n, None) is None]
    if missing: raise ValueError(f"missing policy arrays: {missing}")
    return out

def metrics(ev: Any, P: Any) -> dict[str, float]:
    gpre, gpost, gc = map(np.asarray, (ev.g_pre, ev.g_post_fertility, ev.g_current))
    first = float(gpre[:, :, :, :, :, 0, :].sum() - gpost[:, :, :, :, :, 0, :].sum())
    from build_e5f_native_financing_report import _physical_housing
    owner, rooms, _ = _physical_housing({"g_current": gc, "hR_pol": np.asarray(ev.policy.hR_pol)}, {"H_own": np.asarray(P.H_own)})
    return {"birth_flow": float(np.asarray(ev.births).sum()), "first_birth_flow": first,
            "mean_rooms": float(rooms), "ownership": float(owner),
            "renter_mass": float(gc[:, 0, ...].sum()), "owner_mass": float(gc[:, 1:, ...].sum()),
            "market_residual": float(ev.relative_market_residual), "pre_mass": float(gpre.sum())}

def run_cohort(*args):
    return importlib.import_module("run_e5f_native_income_cohort_diagnostic").run_cohort(*args)

def run_case(x: Mapping[str, Any], phi: float, lam: float, cap: float, out: Path, make_graphs: bool) -> dict[str, Any]:
    native = importlib.import_module("run_e5f_native_financing_diagnostic")
    rental = importlib.import_module("run_e5f_native_rental_access_diagnostic")
    out.mkdir(parents=True, exist_ok=False)
    P = arm_parameters(x["parameters"], phi, lam, cap)
    ev, budget, shared, model = rental.native_solve(x, P)
    grid=x["b_grid"]
    if not np.array_equal(ev.g_pre, x["stationary_g_pre"]): raise ValueError("all arms changed saved initial population")
    rental.gates(ev)
    if float(budget.get("budget_excess_mass", np.inf)) > 2e-10 or float(budget.get("maximum_occupied_excess", np.inf)) > 1e-9:
        raise ValueError(f"budget gate failed: {budget}")
    audit = importlib.import_module("run_e5f_independent_numerical_audit")
    pa = audit.policy_array_audit({"evaluation": ev, "parameters": P, "b_grid": grid}, out)
    if pa["occupied_negative_steps"]: raise ValueError("occupied value monotonicity gate failed")
    if phi == 0.8 and lam == 0.0 and cap == 6.0:
        native.compare_exact(x["evaluation"].policy, ev.policy)
        np.testing.assert_allclose(ev.g_current, x["evaluation"].g_current, atol=1e-10, rtol=0)
        np.testing.assert_allclose(ev.births, x["evaluation"].births, atol=1e-10, rtol=0)
    graph = rental.standard_graphs(x, P, ev, shared, model, out) if make_graphs else {"status": "deferred"}
    if make_graphs and (graph.get("status") != "completed" or graph.get("count") != 17):
        raise RuntimeError(f"stable 17-plot packet failed: {graph}")
    cohort=run_cohort(x,P,ev.policy,out.name,out/"cohort")
    return {"status": "completed", "cohort":cohort,"phi": phi, "lambda": lam, "rental_cap": cap,
            "mortgage_access_label": "native joint deposit-and-collateral access", "metrics": metrics(ev, P),
            "budget": budget, "policy_audit": pa, "standard_graphs": graph,
            "population": "common saved stationary pre-choice mass", "preferences": "family checkpoint unchanged"}

def case_order(mode):
    controls=[(.8,0.,6.),(.8,0.,6.)]
    return controls+([(1.,5.,10.)] if mode=="smoke" else [(p,l,c) for p in PHIS for l in LAMBDAS for c in CAPS])

def run(args):
    contract=validate(args);out=args.output/args.family;out.mkdir(parents=True,exist_ok=False)
    write(out/"contract.json",contract);write(out/"latest_completed.json",{"status":"awaiting_first_case"})
    rows=[];started=time.monotonic();order=case_order(args.mode)
    for i,(phi,lam,cap) in enumerate(order,1):
        label=f"arm_{i:02d}_phi{phi:g}_lambda{lam:g}_cap{cap:g}";case=out/"cases"/label
        cmd=[sys.executable,str(Path(__file__)),"--mode","case","--plan",str(args.plan),"--family",args.family,
             "--checkpoint",str(args.checkpoint),"--expected-hash",args.expected_hash,"--output",str(case),
             "--phi",str(phi),"--lam",str(lam),"--cap",str(cap)]
        case.parent.mkdir(parents=True,exist_ok=True)
        with (case.parent/(label+".log")).open("w") as f:
            proc=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
            try:
                start=time.monotonic();last=0.
                while proc.poll() is None:
                    elapsed=time.monotonic()-start
                    if elapsed>600 or time.monotonic()-started>9000:raise TimeoutError("factorial time budget exceeded")
                    if time.monotonic()-last>30:
                        write(out/"heartbeat.json",{"case":label,"elapsed":elapsed,"updated":time.time()});last=time.monotonic()
                    time.sleep(.5)
            except BaseException:
                from run_e5f_income_candidate_search import _kill_process_group
                _kill_process_group(proc);raise
        if proc.returncode:raise RuntimeError(f"{label} failed; see {label}.log")
        row=json.loads((case/"receipt.json").read_text());row.update(label=label)
        rows.append(row);write(out/"latest_completed.json",row)
        write(out/"completed_summary.json",{"status":"partial","contract":contract,"cases":rows})
    fields=["label","phi","lambda","rental_cap","birth_flow","first_birth_flow","mean_rooms","ownership","cumulative_explicit_births_per_initial_household","first_births_per_initial_household","first_birth_mean_age"]
    with (out/"comparisons.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore");w.writeheader()
        for row in rows:w.writerow({**row,**row["metrics"],**row["cohort"]})
    write(out/"summary.json",{"status":"complete","contract":contract,"cases":rows,"elapsed_seconds":time.monotonic()-started})

def main(argv=None):
    p=argparse.ArgumentParser();p.add_argument("--plan",type=Path,required=True);p.add_argument("--family",choices=FAMILY_NAMES,required=True)
    p.add_argument("--checkpoint",type=Path);p.add_argument("--expected-hash");p.add_argument("--selection-summary",type=Path)
    p.add_argument("--output",type=Path,required=True);p.add_argument("--mode",choices=("inspect","smoke","run","case"),required=True)
    p.add_argument("--phi",type=float);p.add_argument("--lam",type=float);p.add_argument("--cap",type=float)
    a=p.parse_args(argv)
    def stop(signum,frame):raise KeyboardInterrupt("factorial terminated")
    signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
    try:
        if a.mode in ("inspect","case"):
            contract=validate(a)
            if a.mode=="inspect":write(a.output/"inspect.json",contract)
            else:
                if (a.phi,a.lam,a.cap) not in case_order("run"):raise ValueError("unapproved factorial arm")
                row=run_case(packet(a.checkpoint),a.phi,a.lam,a.cap,a.output,True);row["contract"]=contract;write(a.output/"receipt.json",row)
        else:run(a)
    except Exception as exc:
        write(a.output/"failure.json",{"status":"failed","error":str(exc)});raise
if __name__=="__main__":main()
