#!/usr/bin/env python3
"""Prepared, unsubmitted paired stationary search; fail closed on reviewed lock."""
from __future__ import annotations
import argparse, copy, csv, fcntl, gzip, importlib.util, json, math, os, pickle, signal, subprocess, sys, tempfile, time
from pathlib import Path
from types import SimpleNamespace

BASE = Path("/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a")
WORK = BASE / "nightpair_20260925_v1"
SOURCE = WORK / "source"
PLAN = BASE / "utility_overnight_20260923_v1/results/production/B_floor/worker09_proposal16/plan.json"
CHECKPOINT = BASE / "utility_overnight_20260923_v1/results/production/B_floor/worker09_proposal16/result/evaluation/raw/repetition_01/initial_state.pkl.gz"
OBJECTIVE = WORK / "inputs/objective.json"
BANK = WORK / "inputs/proposal_bank.json"
SOURCE_MANIFEST = WORK / "inputs/source_manifest.json"
LOCK = WORK / "inputs/launch_lock.json"
ANCESTOR = WORK / "ancestor_commute.py"
RATES = {"greaney_179": 0.179, "oasi_087510": 0.08751017424959717}
FIRST_BIRTH_EXACT = 1.465
BUDGET_SECONDS = 21600
SEARCH_RESERVE_SECONDS = 4500
EXPORT_RESERVE_SECONDS = 900
WORKERS_PER_ARM = 20
POINTS_PER_WORKER_MAX = 18
OBJECTIVE_CAP_SECONDS = 3100


def load_ancestor():
    spec = importlib.util.spec_from_file_location("paired_commute_ancestor", ANCESTOR)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def read_contract(old, allow_preflight=False):
    if sys.flags.optimize:
        raise RuntimeError("contract assertions require an unoptimized Python runtime")
    expected = os.environ.get("EXPECTED_PAIR_LOCK_SHA256")
    if not expected or len(expected) != 64:
        raise RuntimeError("reviewed pair lock SHA256 environment is required")
    if not LOCK.is_file() or old.sha(LOCK) != expected:
        raise RuntimeError("paired launch lock missing or hash differs")
    lock = json.loads(LOCK.read_text())
    if lock["status"] != "reviewed_ready_to_launch":
        if not (allow_preflight and lock["status"]=="preflight_pending_review"):
            raise RuntimeError("paired launch lock has not passed lead review")
    assert lock["run_budget_seconds"] == BUDGET_SECONDS
    assert lock["search_reserve_seconds"] == SEARCH_RESERVE_SECONDS
    assert lock["export_reserve_seconds"] == EXPORT_RESERVE_SECONDS
    assert lock["workers_per_arm"] == WORKERS_PER_ARM
    assert lock["points_per_worker_max"] == POINTS_PER_WORKER_MAX
    assert lock["rates"] == RATES
    assert lock["first_birth_target_exact"] == FIRST_BIRTH_EXACT
    for path,key in ((OBJECTIVE,"objective_sha256"),(BANK,"proposal_bank_sha256"),
                     (SOURCE_MANIFEST,"source_manifest_sha256"),(ANCESTOR,"ancestor_sha256")):
        if old.sha(path) != lock[key]: raise RuntimeError("frozen pair input differs: " + str(path))
    required_runtime = {"run_pair.py","run_task.sh","submit.sh","render_pair.py",
                        "prepare_bank.py","tools/e5f_earnings_wealth_contract.py",
                        "tools/run_e5f_preference_share_candidate.py"}
    if not required_runtime.issubset(lock["runtime_file_sha256"]):
        raise RuntimeError("runtime file pins incomplete")
    for rel,pin in lock["runtime_file_sha256"].items():
        path=(WORK/rel).resolve()
        if not path.is_relative_to(WORK.resolve()) or old.sha(path)!=pin:
            raise RuntimeError("runtime file differs: "+rel)
    if old.sha(old.TAX_DRIVER)!=lock["tax_driver_sha256"]:
        raise RuntimeError("external verified tax driver differs")
    manifest = json.loads(SOURCE_MANIFEST.read_text())
    files = manifest["files"]
    assert len(files) >= 641 and manifest["source_root"] == str(SOURCE)
    for rel,expected_hash in files.items():
        path = (SOURCE/rel).resolve()
        if not path.is_relative_to(SOURCE.resolve()) or not path.is_file() or old.sha(path) != expected_hash:
            raise RuntimeError("paired source differs: " + rel)
    objective = json.loads(OBJECTIVE.read_text())
    rows = {r["restriction_id"]:r for r in objective["target_rows"]}
    assert len(rows) == 13 and len(objective["target_rows"]) == 13
    assert sum(r["actual_weight"] is not None for r in rows.values()) == 12
    assert rows["first_birth_rooms"]["target"] == FIRST_BIRTH_EXACT
    assert rows["first_birth_rooms"]["actual_weight"] == 137.5652749002964
    assert objective["normalization_tolerance"] == 0.0005
    bounds = {r["parameter"]:(float(r["lower"]),float(r["upper"]))
              for r in objective["parameter_restrictions"]}
    assert set(bounds) == set(old.FREE)
    assert bounds["beta_annual"] == (0.94,0.99)
    assert bounds["h_P"][1] == 2.3
    assert objective["external_restrictions"]["theta1"] == old.THETA1
    bank = json.loads(BANK.read_text())
    assert bank["schema"] == "paired_fixed_proposal_bank_v1"
    assert len(bank["seeds"]) == 2 and len(bank["workers"]) == WORKERS_PER_ARM
    for point in list(bank["seeds"].values()) + [p for pts in bank["workers"].values() for p in pts]:
        assert set(point) == set(old.FREE)
        assert all(math.isfinite(float(point[k])) and low <= float(point[k]) <= high
                   for k,(low,high) in bounds.items())
    assert all(len(bank["workers"][str(i)]) == POINTS_PER_WORKER_MAX
               for i in range(1,WORKERS_PER_ARM+1))
    return lock,objective,bank


def configure(old, arm, lock):
    if arm not in RATES: raise RuntimeError("unknown tax arm")
    old.SOURCE = SOURCE
    old.WORK = WORK
    old.OBJECTIVE = OBJECTIVE
    old.OBJECTIVE_SHA = lock["objective_sha256"]
    old.TAX = RATES[arm]
    original_apply = old.apply_point
    original_normalize = old.normalize_point
    original_evaluate = old.evaluate_point
    last = {}
    def apply_with_clock(selected, point):
        P = original_apply(selected,point)
        P.adult_entry_clock = "split_birth_vintage"
        return P
    def normalize_with_gate(**kwargs):
        result = original_normalize(**kwargs)
        sol,P,price,seconds,normalization = result
        if P.adult_entry_clock != "split_birth_vintage":
            raise RuntimeError("split birth-vintage adult entry was not applied")
        from intergen_eqscale_seq_optimized.adult_entry import require_closed_stationary_renewal
        tolerance = float(json.loads(OBJECTIVE.read_text())["normalization_tolerance"])
        gate = require_closed_stationary_renewal(
            sol.entry_rate,sol.adult_entry_adjusted_birth_children,tolerance)
        gap = float(sol.adult_entry_stationary_relative_gap)
        if not math.isfinite(gap):
            raise RuntimeError("nonfinite adult-entry stationary gap")
        last.clear()
        last.update(adult_entry_stationary_relative_gap=gap,
                    adult_entry_adjusted_birth_children=float(sol.adult_entry_adjusted_birth_children),
                    adult_entry_entry_rate=float(sol.entry_rate),
                    adult_entry_gate=gate)
        return result
    def evaluate_with_gate(**kwargs):
        receipt = original_evaluate(**kwargs)
        receipt["demographic_queue"] = "split 16/20 birth-vintage entry, single /2.1 conversion; final normalized closed-renewal gate passed"
        receipt.update(last)
        # The ancestor table describes the old departure conversion. It is
        # retained only as a diagnostic under the split birth-vintage clock.
        parameter_path=kwargs["output"]/"parameters.csv"
        with parameter_path.open(newline="") as f:
            rows=list(csv.DictReader(f))
        for row in rows:
            if row["parameter"]=="entrant_conversion_factor":
                row["status"]="legacy child-departure diagnostic; inactive in split-birth entry"
        rows.append(dict(parameter="adult_entry_birth_to_household_conversion",
                         estimate=1.0/2.1,lower="",upper="",near_bound="",
                         status="effective closed stationary birth conversion"))
        old.table(parameter_path,rows)
        with gzip.open(kwargs["output"]/"initial_state.pkl.gz","rb") as f:
            saved=pickle.load(f)
        if saved["parameters"].adult_entry_clock!="split_birth_vintage":
            raise RuntimeError("case checkpoint lost split entry clock")
        receipt["adult_entry_birth_to_household_conversion"]=1.0/2.1
        old.write(kwargs["output"]/"receipt.json",receipt)
        return receipt
    old.apply_point = apply_with_clock
    old.normalize_point = normalize_with_gate
    old.evaluate_point = evaluate_with_gate


def prepare(old, lock):
    tax = old.load_tax_driver()
    if old.sha(PLAN) != tax.PLAN_SHA or old.sha(CHECKPOINT) != tax.CHECKPOINT_SHA:
        raise RuntimeError("selected reference plan/checkpoint differs")
    plan = json.loads(PLAN.read_text())
    assert plan["case_id"] == "worker09_proposal16"
    sys.path[:0] = [str(SOURCE/"code/model/tools"),str(SOURCE/"code/model")]
    sys.dont_write_bytecode = True
    with gzip.open(CHECKPOINT,"rb") as f: selected = pickle.load(f)
    tax.check_checkpoint(selected)
    objective = json.loads(OBJECTIVE.read_text())
    tax.SOURCE_SHA = lock["source_manifest_sha256"]
    def selected_preflight():
        return tax,plan,selected,objective
    old.prepare_selected = selected_preflight
    return tax,plan,selected,objective


def deadline(run_root, old):
    run_root.mkdir(parents=True,exist_ok=True)
    mutex = run_root/"deadline.lock"
    with mutex.open("a") as f:
        fcntl.flock(f,fcntl.LOCK_EX)
        path=run_root/"deadline.json"
        if not path.exists():
            now=time.time()
            old.write(path,dict(first_smoke_start_epoch=now,
                deadline_epoch=now+BUDGET_SECONDS,search_cutoff_epoch=now+BUDGET_SECONDS-SEARCH_RESERVE_SECONDS,
                export_cutoff_epoch=now+BUDGET_SECONDS-EXPORT_RESERVE_SECONDS,
                total_seconds=BUDGET_SECONDS))
        data=json.loads(path.read_text())
        fcntl.flock(f,fcntl.LOCK_UN)
    return data


def record_arm(old, run_root, arm, receipt, case):
    arm_root=run_root/arm
    mutex=arm_root/"summary.lock"
    with mutex.open("a") as f:
        fcntl.flock(f,fcntl.LOCK_EX)
        item=dict(case_path=str(case),loss=receipt["loss"],point=receipt["point"],
                  target_system_sha256=receipt["target_system_sha256"],updated_epoch=time.time())
        old.write(arm_root/"latest_completed.json",item)
        best=arm_root/"best_so_far.json"
        if not best.exists() or item["loss"]<json.loads(best.read_text())["loss"]:
            old.write(best,item)
        fcntl.flock(f,fcntl.LOCK_UN)


def run(stage,arm,slot,repeat_id,run_root):
    old=load_ancestor()
    lock,objective,bank=read_contract(old)
    configure(old,arm,lock)
    tax,plan,selected,objective=prepare(old,lock)
    cutoffs=deadline(run_root,old) if stage=="smoke" else json.loads((run_root/"deadline.json").read_text())
    cutoff=cutoffs["search_cutoff_epoch"] if stage in ("smoke","worker") else cutoffs["export_cutoff_epoch"]
    remaining=cutoff-time.time()
    if remaining<=0: raise TimeoutError("shared stage cutoff exhausted")
    def alarm(*_): raise TimeoutError("shared paired-stage cutoff reached")
    signal.signal(signal.SIGALRM,alarm)
    signal.setitimer(signal.ITIMER_REAL,remaining)
    def evaluate_budgeted(**kwargs):
        point_cutoff=min(cutoff,time.time()+OBJECTIVE_CAP_SECONDS)
        kwargs["deadline_epoch"]=point_cutoff
        signal.setitimer(signal.ITIMER_REAL,max(0.001,point_cutoff-time.time()))
        try:
            return old.evaluate_point(**kwargs)
        finally:
            signal.setitimer(signal.ITIMER_REAL,max(0.001,cutoff-time.time()))
    try:
        arm_root=run_root/arm
        arm_root.mkdir(parents=True,exist_ok=True)
        name="smoke" if stage=="smoke" else f"worker_{slot:02d}" if stage=="worker" else f"repeat_{repeat_id:02d}"
        task_root=arm_root/name
        task_root.mkdir(parents=True,exist_ok=False)
        os.environ["NUMBA_CACHE_DIR"]=str(task_root/"numba_cache")
        (task_root/"numba_cache").mkdir()
        runtime=old.setup_runtime(tax,plan,selected,task_root)
        if stage=="smoke":
            for index,(label,point) in enumerate(bank["seeds"].items(),1):
                case=task_root/f"seed_{index:02d}_{label}"
                rec=evaluate_budgeted(tax=tax,objective=objective,selected=selected,runtime=runtime,
                    point=point,output=case,graphs=True)
                record_arm(old,run_root,arm,rec,case)
            old.write(task_root/"complete.json",dict(status="paired_exact_loop_smoke_passed",
                scored=2,graph_count_per_case=17))
        elif stage=="worker":
            assert 1<=slot<=WORKERS_PER_ARM
            assert (arm_root/"smoke/complete.json").exists()
            scored=rejected=attempted=0
            for index,point in enumerate(bank["workers"][str(slot)],1):
                if time.time()>=cutoff: break
                attempted+=1
                case=task_root/f"point_{index:02d}"
                try:
                    rec=evaluate_budgeted(tax=tax,objective=objective,selected=selected,runtime=runtime,
                        point=point,output=case,graphs=False)
                except TimeoutError: raise
                except Exception as exc:
                    case.mkdir(parents=True,exist_ok=True)
                    old.write(case/"failure.json",dict(status="unclassified_point_failure_stop",
                        error_type=type(exc).__name__,error=str(exc),point=point))
                    raise
                scored+=1
                record_arm(old,run_root,arm,rec,case)
                old.write(task_root/"latest_completed.json",dict(index=index,scored=scored,
                    rejected=rejected,loss=rec["loss"],case_path=str(case)))
            old.write(task_root/"complete.json",dict(status="bounded_worker_complete",attempted=attempted,
                scored=scored,rejected=rejected,unrun=POINTS_PER_WORKER_MAX-attempted))
        else:
            assert repeat_id in (1,2)
            selected_best=json.loads((arm_root/"best_so_far.json").read_text())
            case=task_root/"selected_repeat"
            rec=evaluate_budgeted(tax=tax,objective=objective,selected=selected,runtime=runtime,
                point=selected_best["point"],output=case,graphs=False)
            old.write(task_root/"complete.json",dict(status="exact_selected_repeat_complete",
                selected_case=selected_best["case_path"],selected_loss=selected_best["loss"],
                repeated_loss=rec["loss"],repeat_id=repeat_id,case_path=str(case)))
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)


def zero_solve_preflight(old,lock,objective,bank,arm):
    """Install each native runtime and bind both common points without a GE call."""
    records=[]
    (WORK/"results").mkdir(exist_ok=True)
    for arm in (arm,):
        arm_old=load_ancestor()
        configure(arm_old,arm,lock)
        tax,plan,selected,_=prepare(arm_old,lock)
        with tempfile.TemporaryDirectory(prefix="preflight_",dir=WORK/"results") as scratch:
            runtime=arm_old.setup_runtime(tax,plan,selected,Path(scratch))
            if not callable(runtime["solve_balanced_initial_equilibrium"]):
                raise RuntimeError("native PAYGO solve callback missing")
            for label,point in bank["seeds"].items():
                P=arm_old.apply_point(selected,point)
                if (P.adult_entry_clock!="split_birth_vintage" or P.theta1!=arm_old.THETA1
                    or P.delta!=1-(1-arm_old.ANNUAL_DEP)**int(P.period_years)
                    or P.tau_H!=arm_old.ANNUAL_PROPERTY_TAX*int(P.period_years)
                    or not 0.94<=float(point["beta_annual"])<=0.99
                    or not float(point["h_P"])<=2.3 or arm_old.TAX!=RATES[arm]):
                    raise RuntimeError("preflight parameter binding differs")
                records.append(dict(arm=arm,seed=label,tax_fiscal_override=arm_old.TAX,
                    adult_entry_clock=P.adult_entry_clock,theta1=P.theta1,
                    period_depreciation=P.delta,period_property_tax=P.tau_H,
                    beta_annual=point["beta_annual"],h_P=point["h_P"]))
    return records


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--stage",required=True,choices=("preflight","preflight-one","smoke","worker","repeat"))
    p.add_argument("--arm",choices=tuple(RATES),default="greaney_179")
    p.add_argument("--slot",type=int,default=0)
    p.add_argument("--repeat-id",type=int,default=0)
    p.add_argument("--run-root",type=Path)
    a=p.parse_args()
    old=load_ancestor()
    lock,objective,bank=read_contract(old,allow_preflight=a.stage in ("preflight","preflight-one"))
    if a.stage=="preflight":
        records=[]
        for arm in RATES:
            child=subprocess.run([sys.executable,str(Path(__file__).resolve()),"--stage","preflight-one",
                                  "--arm",arm],check=True,capture_output=True,text=True)
            records.extend(json.loads(child.stdout)["bound_seed_runtime_records"])
        print(json.dumps(dict(status="zero_solve_native_preflight_passed",lock_sha256=old.sha(LOCK),
            objective_sha256=lock["objective_sha256"],source_manifest_sha256=lock["source_manifest_sha256"],
            proposal_bank_sha256=lock["proposal_bank_sha256"],seed_count=len(bank["seeds"]),
            workers_per_arm=WORKERS_PER_ARM,points_per_worker=POINTS_PER_WORKER_MAX,
            first_birth_target=FIRST_BIRTH_EXACT,beta_upper=.99,rates=RATES,
            bound_seed_runtime_records=records,native_solve_count=0),sort_keys=True))
        return
    if a.stage=="preflight-one":
        records=zero_solve_preflight(old,lock,objective,bank,a.arm)
        print(json.dumps(dict(status="zero_solve_arm_passed",arm=a.arm,
                              bound_seed_runtime_records=records,native_solve_count=0),sort_keys=True))
        return
    if a.run_root is None: raise SystemExit("--run-root required")
    run(a.stage,a.arm,a.slot,a.repeat_id,a.run_root.resolve())

if __name__=="__main__": main()
