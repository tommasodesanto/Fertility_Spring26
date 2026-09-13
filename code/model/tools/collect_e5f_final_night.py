#!/usr/bin/env python3
"""Collect bounded overnight receipts; optional cluster-side monitoring loop."""
import argparse
import json
from pathlib import Path
import subprocess
import time

NAMES = ("latest_completed.json", "best_so_far.json", "final_summary.json", "summary.json",
         "failure.json", "heartbeat.json", "realized_fit.json", "native_smoke.json",
         "finite_history_complete.json", "joint_progress.json", "search_contract.json")

def read(path):
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        return {"unreadable": str(exc)}

def collect(root, jobs):
    stages = []
    locations = [root / name for name in ("joint_initial_smoke", "initial_search_joint", "age_pilot_joint")]
    for group in ("histories_joint", "histories_refit"):
        locations += sorted((root / group).glob("A*"))
    for folder in locations:
        row = {"stage": str(folder.relative_to(root)), "exists": folder.is_dir()}
        for name in NAMES:
            path = folder / name
            if path.is_file(): row[name] = read(path)
        for path in sorted(folder.glob("policies/*/summary.json")):
            row[str(path.relative_to(folder))] = read(path)
        # Bounded latest numerical progress, not a traversal of binary checkpoints.
        phases = list(folder.glob("window_*/trial_*/latest_*.json"))
        if phases:
            latest = max(phases, key=lambda p: p.stat().st_mtime)
            row["latest_numerical_progress"] = dict(path=str(latest), age_seconds=time.time()-latest.stat().st_mtime, receipt=read(latest))
        row["full_fit_tables"] = [str(p) for p in folder.glob("*target*csv")]
        row["parameter_tables"] = [str(p) for p in folder.glob("*parameter*csv")]
        stages.append(row)
    job_ids=set(jobs.split(','))
    for path in root.glob('histories*/submission.json'):
        job=read(path).get('job','')
        if str(job).isdigit():job_ids.add(str(job))
    jobs=','.join(sorted(job_ids,key=int))
    queue = subprocess.run(["squeue", "-h", "-j", jobs, "-o", "%i|%T|%M|%R"], capture_output=True, text=True, timeout=30)
    return dict(collected_unix=time.time(), jobs=jobs, queue=queue.stdout.splitlines(),
                queue_error=queue.stderr.strip(), stages=stages,
                outstanding=["B household formation unit mapping", "B+ signed migrant allocation", "finite-boundary horizon verification"])

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("runroot",type=Path);parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--jobs",required=True);parser.add_argument("--watch-seconds",type=int,default=0)
    parser.add_argument("--interval",type=int,default=300)
    args=parser.parse_args()
    if not 0<=args.watch_seconds<=43200 or args.interval<60:raise ValueError("Bounded, nonbusy monitoring required")
    if not all(job.isdigit() for job in args.jobs.split(',')):raise ValueError("Numeric job ids required")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    deadline=time.monotonic()+args.watch_seconds
    while True:
        result=collect(args.runroot,args.jobs)
        temporary=args.output.with_suffix('.tmp');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(args.output)
        if time.monotonic()>=deadline or (not result['queue'] and not result['queue_error']):break
        time.sleep(min(args.interval,max(0,deadline-time.monotonic())))

if __name__=="__main__":main()
