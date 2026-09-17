"""Prepare and dispatch the bounded September 13 original-rule experiment.

Run on the Torch login node after copying the four experiment source files to
the declared new batch's source directory. Preparation submits only a native
smoke. That job dispatches independent experiments if and only if it passes.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import time

NUMERICAL_PYTHON="/share/apps/anaconda3/2025.06/bin/python"
BASE=Path("/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913")


def read(path):return json.loads(Path(path).read_text())
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,data):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(p.suffix+".tmp");tmp.write_text(json.dumps(data,indent=2)+"\n");tmp.replace(p)


def submit(spec_path,mode,count,seconds,name,*,dispatch_after=False):
    spec=read(spec_path);batch=Path(spec["batch"]);source=batch/"source"
    args=[NUMERICAL_PYTHON,str(source/"run_e5f_original_queue_experiments.py"),
          "--spec",str(spec_path),"--output",str(batch/name),"--mode",mode,
          "--count",str(count),"--seconds",str(seconds)]
    minutes=int(seconds//60)+10
    script="\n".join(["#!/bin/bash",f"#SBATCH --job-name=e5f_orig_{name}",
        "#SBATCH --account=torch_pr_570_general","#SBATCH --cpus-per-task=1",
        "#SBATCH --mem=32G",f"#SBATCH --time={minutes}",
        f"#SBATCH --output={batch}/logs/{name}_%j.log", "set -euo pipefail",
        "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1",
        shlex.join(args)])+"\n"
    if dispatch_after:
        script+=shlex.join(["/usr/bin/python3",str(source/Path(__file__).name),"--dispatch",str(spec_path)])+"\n"
    path=batch/f"{name}.sbatch";path.write_text(script)
    done=subprocess.run(["sbatch","--parsable",str(path)],text=True,capture_output=True,check=True)
    job=done.stdout.strip().split(";")[0]
    if not job.isdigit():raise RuntimeError("Unexpected sbatch response: "+done.stdout)
    return dict(job_id=int(job),mode=mode,count=count,seconds=seconds,output=str(batch/name),script_sha256=sha(path))


def dispatch(spec_path):
    spec=read(spec_path);batch=Path(spec["batch"])
    receipt=batch/"dispatch.json"
    if receipt.exists():raise ValueError("Already dispatched; no duplicate jobs")
    smoke=read(spec["smoke_summary"])
    if smoke.get("status")!="passed" or smoke.get("spec_sha256")!=sha(spec_path):
        raise ValueError("Native exact-loop smoke has not passed")
    if spec["absolute_deadline_unix"]-time.time()<1800:
        raise TimeoutError("Insufficient afternoon time to dispatch")
    result=dict(spec_sha256=sha(spec_path),smoke_sha256=sha(spec["smoke_summary"]),jobs=[],failures=[])
    # A failure to submit one experiment does not prevent independent arms.
    for mode,count,seconds,name in [
        ("finite",24,14400,"irf_finite_24"),("finite",100,18000,"irf_finite_100"),
        ("terminal",100,21600,"irf_terminal_100"),
        ("history",6,21600,"history_6"),("history",24,21600,"history_24")]:
        try:result["jobs"].append(submit(spec_path,mode,count,seconds,name))
        except Exception as exc:result["failures"].append(dict(name=name,error=str(exc)))
        save(receipt,result)
    result["finished_submission"]=True;save(receipt,result)


def prepare(batch):
    batch=batch.resolve();source=batch/"source";spec_path=batch/"spec.json"
    if spec_path.exists():raise ValueError("Refusing to replace a prepared batch")
    (batch/"logs").mkdir(parents=True,exist_ok=True)
    names=["e5f_original_queue_experiment.py","e5f_original_queue_terminal.py",
        "run_e5f_original_queue_experiments.py","build_e5f_stationary_shock_figures.py",Path(__file__).name]
    for name in names:
        if not (source/name).is_file():raise ValueError("Missing frozen source: "+name)
    base_spec=BASE/"cutoff_horizon_v1/spec.json"
    if not base_spec.exists():base_spec=BASE/"cutoff_horizon_source_v1/spec.json"
    # The previous morning spec is normally stored in cutoff_horizon_v1.
    if not base_spec.exists():raise ValueError("Pinned morning spec must be supplied at its recorded location")
    base=read(base_spec);manifest=read(base["runtime_manifest"])
    manifest.pop("initial_coordinate_seeds",None);manifest.pop("resume_history",None)
    manifest.update(policy_reserve_seconds=3600,forecast_seconds=7200,
        disclosure="Original stationary household distribution and births/2.1 queue at every date; no observed age bridge, no person/headship reset, no migration; equal rebates and balanced PAYGO. Finite boundary remains provisional.",
        population_law="original_household_birth_vintage_queue")
    # Pin new adapter alongside the untouched scientific and empirical pins.
    manifest["file_sha256"].update({str(source/name):sha(source/name) for name in names})
    history_manifest=batch/"history_manifest.json";save(history_manifest,manifest)
    fits_path=Path(base["source_case"])/"realized_fit.json"
    fits=read(fits_path)
    if [row["year"] for row in fits]!=[2007,2011,2015,2019]:raise ValueError("Expected original four-shock source")
    spec=dict(batch=str(batch),base_spec=str(base_spec),history_manifest=str(history_manifest),
        permanent_psi=float(fits[-1]["psi"]),permanent_shock_source=str(fits_path),
        permanent_shock_source_sha256=sha(fits_path),
        smoke_summary=str(batch/"smoke/summary.json"),
        absolute_deadline_unix=time.time()+7*3600,
        numerical_worker_cap=5,structural_parameters_reestimated=False,
        history_shocks_reestimated=True,presentation_results_replaced=False,
        population_law="original_household_birth_vintage_queue",birth_to_entry_conversion=1/2.1,
        historical_age_conditioning=False,person_headship_transition=False,
        file_sha256={str(source/name):sha(source/name) for name in names})
    spec["file_sha256"].update({str(base_spec):sha(base_spec),str(history_manifest):sha(history_manifest),str(fits_path):sha(fits_path)})
    save(spec_path,spec)
    job=submit(spec_path,"smoke",6,2400,"smoke",dispatch_after=True)
    save(batch/"submission.json",dict(smoke=job,spec_sha256=sha(spec_path),long_jobs_submitted=False,
        next="Smoke dispatches five independent jobs only on successful exact-loop verification."))
    print(json.dumps(dict(batch=str(batch),smoke_job=job["job_id"],spec=str(spec_path))))


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument("--prepare",type=Path);g.add_argument("--dispatch",type=Path)
    a=p.parse_args()
    if a.prepare:prepare(a.prepare)
    else:dispatch(a.dispatch)

if __name__=="__main__":main()
