"""Read-only, one-shot status for the fresh 80% purchase-rule calibration.

Run: python3 monitor.py
The command reads Torch receipts through SSH and writes only monitor/ in this
packet. It never invokes a model evaluator or changes Slurm state.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_ROOT = "/scratch/td2248/projects/purchase_fresh_calibration_v1"
REMOTE_CODE = r'''
import hashlib,json,math,pathlib,re,subprocess,time
ROOT=pathlib.Path(ROOT_TEXT)
DESIGN_HASH=EXPECTED_HASH
def load(path):
 try:return json.loads(path.read_text()) if path.exists() else None
 except Exception as e:return {"_read_error":str(e)}
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None
out={"root_exists":ROOT.is_dir(),"remote_root":str(ROOT),"remote_design_sha256":None,
 "remote_design_path":None,"submission":None,"scheduler":{},"slots":[],"warnings":[],"winners":{}}
if not ROOT.is_dir():print(json.dumps(out));raise SystemExit(0)
for path in (ROOT/"design.json",ROOT/"source/design.json",ROOT/"source/fresh_calibration_v1/design.json",
 ROOT/"source/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/fresh_calibration_v1/design.json"):
 if path.exists():out["remote_design_sha256"]=sha(path);out["remote_design_path"]=str(path);break
design_ok=out["remote_design_sha256"]==DESIGN_HASH
receipts=sorted(ROOT.glob("submission_receipt*.json"),key=lambda p:p.stat().st_mtime)
if receipts:out["submission"]=load(receipts[-1])
job_id=str((out["submission"] or {}).get("job_id") or (out["submission"] or {}).get("production_array_job_id") or "")
if re.fullmatch(r"[0-9]+",job_id):
 p=subprocess.run(["squeue","-h","-j",job_id,"-o","%i|%T|%M|%R"],capture_output=True,text=True)
 out["scheduler"]={"job_id":job_id,"returncode":p.returncode,"rows":p.stdout.strip().splitlines(),"stderr":p.stderr.strip()[:300]}
now=time.time();winner={"hard":{"provisional":None,"verified":None},"quarter":{"provisional":None,"verified":None}}
for i in range(24):
 d=ROOT/"results"/f"slot_{i}";c=load(d/"search/cases.json");b=load(d/"search/best_so_far.json");p=load(d/"postcheck/completed.json");t=load(d/"launcher_terminal.json")
 q=d/"search/latest_completed.json";contracts=[load(d/f"{stage}_fresh_contract.json") for stage in ("init","search","postcheck")]
 counts={}
 if isinstance(c,list):
  for case in c:
   key=str(case.get("status","missing_status"))
   if key=="passed":
    try:valid=math.isfinite(float(case.get("base_loss")))
    except (TypeError,ValueError):valid=False
    if not valid:key="passed_without_computed_loss"
   counts[key]=counts.get(key,0)+1
 candidate=b.get("best") if isinstance(b,dict) else None
 selected=p.get("selected") if isinstance(p,dict) else None
 checked=p.get("selected_postcheck") if isinstance(p,dict) else None
 post_ok=(isinstance(selected,dict) and isinstance(checked,dict) and
          p.get("status")=="selected_numerically_verified" and checked.get("status")=="passed" and
          selected.get("base_loss")==checked.get("base_loss") and
          selected.get("target_fit")==checked.get("base_target_fit"))
 search_hash=contracts[1].get("design_sha256") if isinstance(contracts[1],dict) else None
 post_hash=contracts[2].get("design_sha256") if isinstance(contracts[2],dict) else None
 search_authenticated=design_ok and search_hash==DESIGN_HASH
 post_authenticated=search_authenticated and post_hash==DESIGN_HASH
 def valid_point(point):
  if not isinstance(point,dict) or point.get("status")!="passed":return False
  try:return math.isfinite(float(point.get("base_loss")))
  except (TypeError,ValueError):return False
 age=round(now-q.stat().st_mtime) if q.exists() else None
 row={"slot":i,"result_dir_exists":d.is_dir(),"counts":counts,"checkpoint_age_seconds":age,
      "terminal":t,"postcheck_status":p.get("status") if isinstance(p,dict) else None,
      "postcheck_exact":post_ok,"provisional_loss":candidate.get("base_loss") if isinstance(candidate,dict) else None,
      "provisional_case":candidate.get("label") if isinstance(candidate,dict) else None,
      "verified_loss":checked.get("base_loss") if post_ok else None,
      "search_authenticated":search_authenticated,"postcheck_authenticated":post_authenticated,
      "contracts":{stage:(x.get("design_sha256") if isinstance(x,dict) else None) for stage,x in zip(("init","search","postcheck"),contracts)}}
 out["slots"].append(row)
 for stage,value,point,authenticated in (("provisional",row["provisional_loss"],candidate,search_authenticated),
                                         ("verified",row["verified_loss"],selected,post_authenticated and post_ok)):
  if authenticated and valid_point(point) and isinstance(value,(int,float)) and math.isfinite(float(value)):
   arm="hard" if i<12 else "quarter";prior=winner[arm][stage]
   if prior is None or value<prior["loss"]:winner[arm][stage]={"slot":i,"loss":value,"point":point,"postcheck":checked if stage=="verified" else None}
out["winners"]=winner
print(json.dumps(out))
'''


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    keys = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--remote-root", default=DEFAULT_ROOT)
    parser.add_argument("--ssh-host", default="torch")
    parser.add_argument("--output-dir", type=Path, default=HERE / "monitor")
    args = parser.parse_args()

    design_path = HERE / "design.json"
    raw_design = design_path.read_bytes()
    design = json.loads(raw_design)
    starts = design.get("starts", [])
    if len(starts) != 24 or any(start.get("slot") != i or start.get("arm") != ("hard" if i < 12 else "quarter") for i, start in enumerate(starts)):
        raise SystemExit("Design must map slots 0–11 to hard and 12–23 to quarter")
    design_hash = hashlib.sha256(raw_design).hexdigest()
    remote = "ROOT_TEXT=" + repr(args.remote_root) + "\nEXPECTED_HASH=" + repr(design_hash) + "\n" + REMOTE_CODE
    ssh_command = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
                   "-o", "ConnectionAttempts=1", "-o", "ServerAliveInterval=15",
                   "-o", "ServerAliveCountMax=1", args.ssh_host, "python3", "-"]
    try:
        proc = subprocess.run(ssh_command, input=remote, text=True,
                              capture_output=True, timeout=45)
    except subprocess.TimeoutExpired as exc:
        failure = f"Remote read timed out after {exc.timeout} seconds"
    else:
        failure = "Remote read failed: " + proc.stderr.strip()[:500] if proc.returncode else None
    if failure:
        out = args.output_dir.resolve()
        out.mkdir(parents=True, exist_ok=True)
        (out / "failed_attempt.json").write_text(json.dumps({
            "as_of_utc": datetime.now(timezone.utc).isoformat(),
            "remote_root": args.remote_root,
            "ssh_host": args.ssh_host,
            "error": failure,
        }, indent=2) + "\n")
        raise SystemExit(failure)
    snapshot = json.loads(proc.stdout)
    snapshot.update(as_of_utc=datetime.now(timezone.utc).isoformat(), local_design_sha256=design_hash,
                    target_fingerprint=design.get("target_fingerprint"), weight_fingerprint=design.get("weight_fingerprint"))
    warnings = snapshot["warnings"]
    remote_hash = snapshot.get("remote_design_sha256")
    if remote_hash is None:
        warnings.append("Remote design file unavailable: design identity is unverified")
    elif remote_hash != design_hash:
        warnings.append("Remote design SHA-256 differs from local design: do not compare or use results")
    for row in snapshot["slots"]:
        expected_arm = starts[row["slot"]]["arm"]
        row["arm"] = expected_arm
        for stage, digest in row["contracts"].items():
            if digest is not None and digest != design_hash:
                warnings.append(f"slot {row['slot']} {stage} contract design hash mismatch")
        if row["provisional_loss"] is not None and not row["search_authenticated"]:
            warnings.append(f"slot {row['slot']} provisional result suppressed: search design identity unavailable or mismatched")
        if row["postcheck_status"] is not None and not row["postcheck_authenticated"]:
            warnings.append(f"slot {row['slot']} postcheck result suppressed: postcheck design identity unavailable or mismatched")
        if row["terminal"] is not None and row["terminal"].get("exit_code") != 0:
            warnings.append(f"slot {row['slot']} launcher exit {row['terminal'].get('exit_code')}")
        if row["terminal"] is None and row["checkpoint_age_seconds"] is not None and row["checkpoint_age_seconds"] > 1800:
            warnings.append(f"slot {row['slot']} checkpoint older than 30 minutes")

    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    lines = ["# Fresh 80% calibration monitor", "", f"Snapshot {snapshot['as_of_utc']}. Remote root: `{args.remote_root}`.",
             f"Local design SHA-256: `{design_hash}`. Remote design SHA-256: `{remote_hash or 'unavailable'}`.", ""]
    if not snapshot["root_exists"]:
        lines += ["The remote campaign has not been staged yet.", ""]
    else:
        scheduler = snapshot.get("scheduler") or {}
        lines += [f"Slurm job: `{scheduler.get('job_id', 'not submitted')}`; active rows: {len(scheduler.get('rows', []))}.", "",
                  "| Arm | Valid scored cases | Budget-uncomputed | Numerically inadmissible | Other uncomputed | Terminal / 12 | Fresh postchecks / 12 |",
                  "|---|---:|---:|---:|---:|---:|---:|"]
        for arm in ("hard", "quarter"):
            rows = [r for r in snapshot["slots"] if r["arm"] == arm]
            total = lambda key: sum(r["counts"].get(key, 0) for r in rows)
            other = sum(sum(v for k, v in r["counts"].items() if k not in ("passed", "budget_exhausted", "inadmissible_numerical")) for r in rows)
            lines.append(f"| {arm} | {total('passed')} | {total('budget_exhausted')} | {total('inadmissible_numerical')} | {other} | {sum(r['terminal'] is not None for r in rows)} | {sum(r['postcheck_exact'] for r in rows)} |")
        lines += ["", "A completed model attempt is counted as valid only when its case status is `passed` and it has a computed loss. The other statuses are listed separately.", ""]
        for arm in ("hard", "quarter"):
            lines += [f"## {arm.capitalize()}", ""]
            for stage in ("verified", "provisional"):
                win = snapshot["winners"][arm][stage]
                if win is None:
                    lines.append(f"{stage.capitalize()} candidate: unavailable.")
                    continue
                point = win.pop("point")
                win.pop("postcheck", None)
                prefix = f"{arm}_{stage}_slot{win['slot']}"
                fit = point.get("target_fit") or point.get("base_target_fit") or []
                params = point.get("effective_parameters") or []
                if isinstance(fit, list) and isinstance(params, list) and len(fit) == 14 and len(params) == 31:
                    fit_path, par_path = out / f"{prefix}_target_fit.csv", out / f"{prefix}_parameters.csv"
                    write_csv(fit_path, fit)
                    write_csv(par_path, params)
                    win.update(target_fit_path=str(fit_path), parameters_path=str(par_path))
                    lines.append(f"{stage.capitalize()} candidate: loss **{win['loss']:.6f}**, slot {win['slot']}; [14-row target fit]({fit_path}), [31-parameter bounds]({par_path}).")
                else:
                    warnings.append(f"{arm} {stage} slot {win['slot']} missing complete 14/31 tables")
                    lines.append(f"{stage.capitalize()} candidate: loss **{win['loss']:.6f}**, slot {win['slot']}; complete 14/31 tables unavailable.")
            lines.append("")
        if any(snapshot["winners"][arm]["verified"] is not None for arm in ("hard", "quarter")):
            lines += ["Fresh selected-point verification confirms the quoted numerical point. It does not certify optimizer convergence or adopt the experimental economics.", ""]
    if warnings:
        lines += ["## Health and identity warnings", ""] + [f"- {warning}" for warning in warnings] + [""]
    (out / "STATUS.md").write_text("\n".join(lines))
    (out / "latest_status.json").write_text(json.dumps(snapshot, indent=2) + "\n")
    print(f"Wrote {out / 'STATUS.md'}; warnings={len(warnings)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
