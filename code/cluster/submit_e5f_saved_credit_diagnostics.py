#!/usr/bin/env python3
"""Stage and optionally submit the bounded saved-cohort credit diagnostic.

Preparation is read-only locally and submission is opt-in.  The remote stage is
fresh by construction, and all inputs are copied from the pinned local summaries
and analyzer before either Slurm job is submitted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shlex
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ANALYZER = ROOT / "code/model/tools/analyze_e5f_saved_credit_diagnostics.py"
SUMMARY_BASE = ROOT / "output/model/native_financing_diagnostic_20260919/overnight/final_mechanisms"
REMOTE_HOST = "torch"
REMOTE_STAGE = "/scratch/td2248/projects/Fertility_Spring26_specification_20260920/saved_credit_v1"
FAMILIES = ("original", "stationary_new_income", "refit_new_income")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def run(argv: list[str], *, dry_run: bool, input_text: str | None = None) -> str:
    if argv[0] == 'ssh' and argv[2:4] == ['bash', '-lc']:
        argv = ['ssh', '-o', 'BatchMode=yes', argv[1], 'bash -lc ' + shlex.quote(argv[4])]
    print("$", shlex.join(argv))
    if dry_run:
        return ""
    return subprocess.run(argv, input=input_text, text=True, check=True,
                          capture_output=True).stdout.strip()


def summary_paths() -> dict[str, Path]:
    return {"original": SUMMARY_BASE / "original/summary.json",
            "stationary_new_income": SUMMARY_BASE / "stationary_new_income/summary.json",
            "refit_new_income": SUMMARY_BASE / "refit_new_income/summary.json"}


def command_argv(family_args: list[str], output: str, smoke: bool, wall: int) -> list[str]:
    args = ["python", "-u", f"{REMOTE_STAGE}/analyze_e5f_saved_credit_diagnostics.py",
            "--summary-original", f"{REMOTE_STAGE}/summary_original.json",
            "--summary-pilot", f"{REMOTE_STAGE}/summary_stationary_new_income.json",
            "--summary-refit", f"{REMOTE_STAGE}/summary_refit_new_income.json",
            "--output", f"{REMOTE_STAGE}/{output}", "--strict", "--wall-time-seconds", str(wall),
            "--per-case-seconds", "120"]
    if smoke:
        args.append("--smoke")
    args += ["--families", *family_args]
    return args


def receipt_check(output: str, expected_cases: int, expected_families: int) -> str:
    code = ("import json; from pathlib import Path; "
            f"d=json.loads(Path('{REMOTE_STAGE}/{output}/receipt.json').read_text()); "
            f"assert d.get('status')=='complete' and len(d.get('families',{{}}))=={expected_families}; "
            "cases=[c for f in d['families'].values() for c in f.get('cases',[])]; "
            f"assert len(cases)=={expected_cases}; "
            "assert all(len(c.get('supplemental_plots',{}).get('paths',[]))==2 for c in cases)")
    return shlex.join(["python", "-c", code])


def main() -> None:
    global REMOTE_STAGE
    parser = argparse.ArgumentParser()
    parser.add_argument("--submit", action="store_true", help="stage and submit; default is dry-run")
    parser.add_argument("--families", nargs="+", choices=FAMILIES, default=list(FAMILIES),
                        help="families for the full run; smoke remains refit_new_income")
    parser.add_argument("--remote-stage", default=REMOTE_STAGE)
    args = parser.parse_args()
    REMOTE_STAGE = args.remote_stage
    if not REMOTE_STAGE.startswith('/scratch/td2248/projects/Fertility_Spring26_specification_20260920/') or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_/-' for c in REMOTE_STAGE):
        raise ValueError('invalid remote stage')
    local = ROOT / 'output/model/native_financing_diagnostic_20260919/specification_followup/saved_credit_v1'
    if args.submit:
        local.mkdir(exist_ok=False)

    summaries = summary_paths()
    files = {"analyzer": ANALYZER, **{f"summary_{family}": path for family, path in summaries.items()}}
    missing = [str(path) for path in files.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("missing local input(s): " + ", ".join(missing))
    manifest = {"remote_stage": REMOTE_STAGE, "files": {key: {"local": str(path), "sha256": digest(path)} for key, path in files.items()},
                "families": args.families, "smoke_family": ["refit_new_income"],
                "smoke": {"wall_seconds": 240, "outer_timeout": 270, "slurm_time": "00:05:00"},
                "full": {"wall_seconds": 720, "outer_timeout": 840, "slurm_time": "00:15:00"},
                "environment": {"threads": "1", "python_module": "anaconda3/2025.06"}}
    manifest.update(household_solves=0, hypothesis='Measure retained cohort debt and statutory/grid support; treatment policy arrays were not retained.', stop_rule='Stop on failed source, cohort, flow, entry, plot or time gate; no automatic retry.')
    if args.submit:
        (local/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest, indent=2, sort_keys=True))
    remote = REMOTE_HOST
    run(["ssh", remote, "bash", "-lc", f"test ! -e {shlex.quote(REMOTE_STAGE)} && mkdir -p {shlex.quote(REMOTE_STAGE)}"], dry_run=not args.submit)
    if args.submit:
        for key, path in files.items():
            run(["scp", str(path), f"{remote}:{REMOTE_STAGE}/{key}.json" if key.startswith("summary_") else f"{remote}:{REMOTE_STAGE}/analyze_e5f_saved_credit_diagnostics.py"], dry_run=False)
        run(["ssh", remote, "bash", "-lc", f"test ! -e {shlex.quote(REMOTE_STAGE + '/smoke')} && test ! -e {shlex.quote(REMOTE_STAGE + '/full')}"], dry_run=False)
        run(["ssh", remote, "bash", "-lc", f"cat > {shlex.quote(REMOTE_STAGE + '/manifest.json')}"], dry_run=False,
            input_text=json.dumps(manifest, indent=2) + "\n")
        verify_code = ("import hashlib,json; from pathlib import Path; "
                       f"m=json.loads(Path('{REMOTE_STAGE}/manifest.json').read_text()); "
                       "[(lambda p,e: (_ for _ in ()).throw(SystemExit('hash mismatch '+str(p))) if hashlib.sha256(p.read_bytes()).hexdigest()!=e else None)(Path('" + REMOTE_STAGE + "/'+k+'.json') if k.startswith('summary_') else Path('" + REMOTE_STAGE + "/analyze_e5f_saved_credit_diagnostics.py'), v['sha256']) for k,v in m['files'].items()]")
        run(["ssh", remote, "bash", "-lc", 'module load anaconda3/2025.06; '+shlex.join(["python", "-c", verify_code])], dry_run=False)

    env = "set -euo pipefail; module load anaconda3/2025.06; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONUNBUFFERED=1; "
    hash_cmd = shlex.join(['python','-c',verify_code])+'; ' if args.submit else 'true; '
    smoke_cmd = 'bash -lc ' + shlex.quote(env + hash_cmd + "timeout 270s " + shlex.join(command_argv(["refit_new_income"], "smoke", True, 240)) + '; ' + receipt_check('smoke',4,1))
    full_cmd = 'bash -lc ' + shlex.quote(env + hash_cmd + receipt_check("smoke", 4, 1) + " && timeout 840s " + shlex.join(command_argv(list(args.families), "full", False, 720)) + " && " + receipt_check("full", 4 * len(args.families), len(args.families)))
    smoke_sbatch = ["sbatch", "--parsable", "--account=torch_pr_570_general", "--job-name=e5f_credit_smoke",
                    "--cpus-per-task=1", "--mem=24G", "--time=00:05:00", f"--output={REMOTE_STAGE}/smoke_%j.out",
                    "--wrap", smoke_cmd]
    smoke_job = run(["ssh", remote, "bash", "-lc", shlex.join(smoke_sbatch)], dry_run=not args.submit)
    if args.submit:
        (local/'submission.json').write_text(json.dumps({'smoke_job':smoke_job,'remote_root':REMOTE_STAGE,'full_job':None},indent=2)+'\n')
    full_sbatch = ["sbatch", "--parsable", "--account=torch_pr_570_general", "--dependency=afterok:" + (smoke_job or "SMOKE_JOB_ID"),
                   "--job-name=e5f_credit_full", "--cpus-per-task=1", "--mem=24G", "--time=00:15:00",
                   f"--output={REMOTE_STAGE}/full_%j.out", "--wrap", full_cmd]
    full_job = run(["ssh", remote, "bash", "-lc", shlex.join(full_sbatch)], dry_run=not args.submit)
    job_ids = {"smoke": smoke_job or None, "full": full_job or None, "manifest_sha256": digest(ANALYZER)}
    print(json.dumps(job_ids, indent=2, sort_keys=True))
    if args.submit:
        (local/'submission.json').write_text(json.dumps({'smoke_job':smoke_job,'full_job':full_job,'remote_root':REMOTE_STAGE},indent=2)+'\n')
        run(["ssh", remote, "bash", "-lc", f"cat > {shlex.quote(REMOTE_STAGE + '/job_ids.json')}"], dry_run=False,
            input_text=json.dumps(job_ids, indent=2) + "\n")


if __name__ == "__main__":
    main()
