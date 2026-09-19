"""Bounded local mini-calibration controller (no changes to model sources)."""
from __future__ import annotations
import argparse, csv, hashlib, json, math, os, random, signal, subprocess, sys, threading, time, concurrent.futures
from pathlib import Path

OBJECTIVE = "4440ea07f4de957740ca6c04961d2806d9b9ef782c7a0e7dad4ce73e1db651b1"
PARAMETERS = ("H0", "beta_annual", "chi", "first_birth_fixed_cost", "h_P",
              "kappa_fert", "kappa_fert_continuation", "theta0", "theta1")
CAPS = {"beta_annual": (0.94, .99), "h_P": (.1, 2.3)}
CASE_SECONDS, WORKERS, RESERVE, RSS_LIMIT = 420, 6, 240, 24 * 1024**3

def load(path): return json.loads(Path(path).read_text())
def atomic(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    tmp.replace(path)
def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

def _find_rows(obj):
    if isinstance(obj, dict):
        for key in ("parameters", "parameter_restrictions", "restrictions"):
            if isinstance(obj.get(key), list): yield from obj[key]
        for v in obj.values(): yield from _find_rows(v)
    elif isinstance(obj, list):
        for v in obj: yield from _find_rows(v)

def restrictions(objective, reference=None):
    out = {}
    for row in _find_rows(objective):
        if not isinstance(row, dict): continue
        name = row.get("parameter") or row.get("name")
        if name in PARAMETERS and "lower" in row and "upper" in row:
            try: out[name] = (float(row["lower"]), float(row["upper"]))
            except (TypeError, ValueError): pass
    for name, bounds in CAPS.items(): out[name] = bounds
    if reference:
        for row in reference.get("score", {}).get("parameters", []):
            name = row.get("parameter")
            if name in PARAMETERS and name not in out:
                try: out[name] = (float(row["lower"]), float(row["upper"]))
                except (TypeError, ValueError): pass
    missing = [p for p in PARAMETERS if p not in out]
    if missing: raise ValueError("objective contract lacks bounds: " + ",".join(missing))
    return out

def reference_values(reference):
    score = reference.get("score", reference)
    pars = {r["parameter"]: float(r["estimate"]) for r in score.get("parameters", [])
            if r.get("parameter") in PARAMETERS}
    norm = score.get("normalization", {})
    if len(pars) != len(PARAMETERS) or not norm: raise ValueError("reference lacks full score")
    return pars, norm

def proposals(original, bounds, seed=20260919):
    p = dict(original); seen = set(); result = []
    def add(label, values):
        vals = dict(values); key = tuple(round(float(vals[x]), 14) for x in PARAMETERS)
        if key in seen: raise ValueError("duplicate proposal")
        seen.add(key); result.append({"case_id": label, "initial_psi": p.get("initial_psi"),
                                      "parameters": vals, "repetitions": 1})
    for name in PARAMETERS:
        lo, hi = bounds[name]; x = float(p[name]); step = .001 if name == "beta_annual" else (.01 * abs(x) if x else .01 * (hi-lo))
        y = x + step
        if y > hi: y = x - step
        if y < lo: y = x + step
        if not lo <= y <= hi: raise ValueError(f"outbound perturbation {name}")
        q = dict(p); q[name] = y; add("coord_%s" % name, q)
    rng = random.Random(seed)
    for j in range(3):
        q = dict(p)
        for name in PARAMETERS:
            lo, hi = bounds[name]
            if name == "beta_annual": delta = -.0005
            else: delta = (1 if rng.random() >= .5 else -1) * .005 * abs(float(q[name]))
            q[name] = min(hi, max(lo, float(q[name]) + delta))
        add(f"joint_{j+1}", q)
    return result

def _numeric(a, b, path=""):
    if isinstance(a, bool) or isinstance(b, bool): return a == b
    if isinstance(a, (int,float)) and isinstance(b, (int,float)):
        return math.isclose(float(a), float(b), rel_tol=1e-10, abs_tol=1e-8)
    if isinstance(a, dict) and isinstance(b, dict):
        ignored = ("hash", "sha", "path", "seconds", "timing", "runtime", "checkpoint")
        keys = [k for k in a if not any(x in k.lower() for x in ignored)]
        return set(keys) == {k for k in b if not any(x in k.lower() for x in ignored)} and all(_numeric(a[k], b[k], path+"."+k) for k in keys)
    if isinstance(a, list) and isinstance(b, list): return len(a)==len(b) and all(_numeric(x,y,path) for x,y in zip(a,b))
    return a == b

def compare_score(actual, reference):
    aa = actual.get("score", actual); rr = reference.get("score", reference)
    fields = ("loss", "target_fit", "parameters", "normalization")
    return all(_numeric(aa.get(k), rr.get(k), k) for k in fields)

def process_tree(roots):
    rows=[list(map(int,r.split())) for r in subprocess.check_output(["ps","-axo","pid=,ppid=,pgid=,rss="],text=True).splitlines()]
    ids=set(roots)
    while True:
        children={r[0] for r in rows if r[1] in ids}
        if children<=ids: break
        ids|=children
    return [r for r in rows if r[0] in ids]

def kill_group(proc):
    groups={r[2] for r in process_tree({proc.pid})}
    for group in groups:
        try: os.killpg(group,signal.SIGKILL)
        except ProcessLookupError: pass
    proc.wait()

def run_process(command, folder, timeout=CASE_SECONDS, env=None, active=None, lock=None):
    folder = Path(folder); folder.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic(); log = (folder.parent / (folder.name + ".log")).open("w")
    proc = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, env=env,
                            start_new_session=True)
    if active is not None:
        with lock: active[proc.pid] = proc
    try: proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        kill_group(proc); return {"status":"timeout", "elapsed_seconds":time.monotonic()-started}
    finally:
        log.close()
        if active is not None:
            with lock: active.pop(proc.pid, None)
    result = {"status":"exited", "returncode":proc.returncode,
              "elapsed_seconds":time.monotonic()-started}
    receipt = folder / "candidate_result.json"
    if receipt.exists():
        try:
            payload = json.loads(receipt.read_text()); result.update(payload)
            result["status"] = payload.get("status", result["status"])
        except json.JSONDecodeError: result["status"] = "failed"
    elif proc.returncode != 0: result["status"] = "failed"
    return result

def rss_groups(groups):
    return sum(r[3]*1024 for r in process_tree(groups)) if groups else 0

def write_tables(out, records):
    for name, key in (("all_target_fits.csv", "target_fit"), ("all_parameters.csv", "parameters")):
        rows=[]
        for r in records:
            for row in r.get("score", {}).get(key, []):
                row=dict(row)
                if key=="parameters" and row.get("parameter")=="beta_annual":
                    row.update(upper=.99, near_bound=abs(float(row["estimate"])-.99)<.0005)
                rows.append(dict(case_id=r.get("case_id"), **row))
        if rows:
            with (Path(out)/name).open("w", newline="") as f:
                w=csv.DictWriter(f, lineterminator="\n", fieldnames=list(dict.fromkeys(k for x in rows for k in x))); w.writeheader(); w.writerows(rows)

def write_comparison(out, reference, selected):
    for key, identifier in (("target_fit","restriction_id"),("parameters","parameter")):
        base={r[identifier]:r for r in reference["score"][key]}
        rows=[]
        for row in selected["score"][key]:
            r=dict(row); old=base[r[identifier]]
            if key=="target_fit":
                for field in ("model","gap","loss_contribution"): r["baseline_"+field]=old[field]
            else:
                r["baseline_estimate"]=old["estimate"]
                if r[identifier]=="beta_annual": r.update(upper=.99,near_bound=abs(float(r["estimate"])-.99)<.0005)
            rows.append(r)
        with (Path(out)/("comparison_"+key+".csv")).open("w",newline="") as f:
            writer=csv.DictWriter(f,lineterminator="\n",fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
            writer.writeheader(); writer.writerows(rows)

def main(argv=None):
    ap=argparse.ArgumentParser(); ap.add_argument("--output",type=Path,required=True); ap.add_argument("--template",type=Path,required=True)
    ap.add_argument("--runtime-dir",type=Path,required=True); ap.add_argument("--objective",type=Path,required=True); ap.add_argument("--smoke-only",action="store_true")
    ap.add_argument("--validated-baseline",type=Path)
    args=ap.parse_args(argv); out=args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    recipe=Path("output/model/local_mini_calibration_20260919/staged/recipe").resolve(); driver=recipe/"run_e5f_joint_rebated_initial_scored.py"
    helper=recipe/"run_e5f_rebated_initial_overnight.py"; joint=recipe/"run_e5f_joint_rebated_initial_probe.py"; original=load(recipe/"proposal.json")
    refpath=recipe/"native_output/candidate_result.json"
    if not refpath.exists(): refpath=Path("output/model/paper_baseline_sep14/replay_20260917/native_output/candidate_result.json")
    if sha(helper)!="d9aa97b890442d45971ec622b4b41687da10ffa26eedaf0198f352c7e6ecb790" or sha(joint)!="9eee3bca39f2a98f4a58cf18196d695b6a9db3e93ef18f8eaa2bf4dbb1243bbb":
        raise ValueError("frozen recipe fingerprint mismatch")
    ref=load(refpath); obj=load(args.objective); bounds=restrictions(obj,ref); pars,_=reference_values(ref)
    canonical = hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode()).hexdigest()
    if canonical != OBJECTIVE: raise ValueError("objective fingerprint mismatch")
    all_records=[]; active={}; lock=threading.Lock(); stop=threading.Event(); finished=threading.Event(); peak=0; started=time.monotonic()
    deadline=started+1500
    def monitor():
        nonlocal peak
        while not finished.is_set():
            with lock: processes=list(active.values())
            rss=rss_groups({p.pid for p in processes}); peak=max(peak,rss)
            reason="memory_limit" if rss>RSS_LIMIT else ("time_limit" if time.monotonic()>deadline else None)
            atomic(out/"heartbeat.json",dict(elapsed_seconds=time.monotonic()-started,active_processes=len(processes),rss_bytes=rss,peak_rss_bytes=peak,stop_reason=reason))
            if reason:
                stop.set()
                for proc in processes: kill_group(proc)
            finished.wait(5)
    threading.Thread(target=monitor,daemon=True).start()
    env=dict(os.environ, PYTHONPATH=str(args.runtime_dir), NUMBA_CACHE_DIR=str(out/"cache"), OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", NUMBA_NUM_THREADS="1")
    def run(item, repetitions=1):
        nonlocal peak
        case=out/"cases"/item["case_id"]; proposal=dict(item, repetitions=repetitions); atomic(case.parent/(item["case_id"]+".json"),proposal)
        cmd=[sys.executable,"-B",str(driver),"--helper",str(helper),"--helper-sha256",sha(helper),"--joint",str(joint),"--joint-sha256",sha(joint),"--template",str(args.template),"--output",str(case),"--proposal",str(case.parent/(item["case_id"]+".json"))]
        res=run_process(cmd,case,timeout=max(1,min(CASE_SECONDS,deadline-time.monotonic())),env=env,active=active,lock=lock); return {**res, "case_id":item["case_id"]}
    baseline=dict(original,case_id="baseline_smoke",repetitions=1); rec=load(args.validated_baseline) if args.validated_baseline else run(baseline); rec["case_id"]="baseline_smoke"; all_records.append(rec); atomic(out/"latest_completed.json",rec); atomic(out/"best_so_far.json",rec)
    if rec["status"] != "verified" or not rec.get("score") or not compare_score(rec, ref): atomic(out/"summary.json",{"status":"failed_smoke", "peak_rss_bytes":peak}); return 2
    if args.smoke_only: atomic(out/"summary.json",{"status":"smoke_complete", "peak_rss_bytes":peak}); return 0
    items = proposals(pars,bounds)
    for item in items: item["initial_psi"] = original.get("initial_psi")
    def save():
        valid=[x for x in all_records if x.get("status")=="verified"]
        atomic(out/"latest_completed.json",all_records[-1])
        if valid: atomic(out/"best_so_far.json",min(valid,key=lambda x:x["score"]["loss"]))
        write_tables(out,all_records)
    for offset in range(0,len(items),WORKERS):
        if stop.is_set() or time.monotonic()+CASE_SECONDS+RESERVE>deadline: break
        with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as pool:
            futures=[pool.submit(run,item) for item in items[offset:offset+WORKERS]]
            for future in concurrent.futures.as_completed(futures):
                rec=future.result(); all_records.append(rec); save()
        if any(x.get("status")!="verified" for x in all_records): break
    selected=min((x for x in all_records if x.get("status")=="verified"),key=lambda x:x["score"]["loss"])
    repeat_status="not_run_budget"
    if not stop.is_set() and deadline-time.monotonic()>30:
        proposal=baseline if selected["case_id"]=="baseline_smoke" else next(x for x in items if x["case_id"]==selected["case_id"])
        final=run(dict(proposal,case_id="selected_final_repeat")); all_records.append(final); save()
        repeat_status="verified" if final.get("status")=="verified" and compare_score(final,selected) else "failed"
    finished.set()
    atomic(out/"selected.json",selected)
    selected_dir=out/"selected_tables"; selected_dir.mkdir(); write_tables(selected_dir,[selected])
    write_comparison(out,ref,selected)
    atomic(out/"summary.json",dict(status="completed_diagnostic" if repeat_status=="verified" else "provisional",completed_cases=len(all_records),peak_rss_bytes=peak,elapsed_seconds=time.monotonic()-started,objective_sha256=OBJECTIVE,final_repeat=repeat_status,baseline_loss=all_records[0]["score"]["loss"],selected_case=selected["case_id"],selected_loss=selected["score"]["loss"]))
    return 0
if __name__ == "__main__": raise SystemExit(main())
