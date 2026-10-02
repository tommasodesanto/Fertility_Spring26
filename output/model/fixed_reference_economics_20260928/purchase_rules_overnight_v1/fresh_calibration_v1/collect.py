"""Collect the finished fresh 80% calibration without invoking the model.

The remote validation runs read-only. Files are copied only into this packet's
collection directory; the earlier calibration selection is never modified.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = "/scratch/td2248/projects/purchase_fresh_calibration_v1"
OUT = HERE / "collection"
RUNNER_REL = "output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/run_psi.py"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


REMOTE = r'''
import csv,hashlib,json,math,pathlib
ROOT=pathlib.Path(ROOT_TEXT)
DESIGN=json.loads(DESIGN_TEXT)
SUBMIT=json.loads(SUBMIT_TEXT)
PINS=json.loads(PINS_TEXT)
RUNNER_SHA=RUNNER_SHA_TEXT
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def rows(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def require(ok,why):
 if not ok:raise ValueError(why)
def close(a,b,eps=1e-8):return math.isclose(float(a),float(b),rel_tol=0,abs_tol=eps)
require(ROOT.is_dir(),'missing fresh root')
require(read(ROOT/'source/design.json')==DESIGN,'remote design content drift')
require(read(ROOT/'submission_receipt.json')==SUBMIT,'remote submission drift')
require(read(ROOT/'source_sha256.json')==PINS,'remote source pins drift')
for name,expected in PINS.items():require(sha(ROOT/'source'/name)==expected,'source hash drift: '+name)
require(SUBMIT['design_sha256']==sha(ROOT/'source/design.json'),'submission design hash drift')
require(SUBMIT['source_sha256']==PINS and SUBMIT['job_id']=='19040483','submission identity drift')
require(SUBMIT['economic_contract_unchanged'] is True,'economic contract drift')
require(len(DESIGN['starts'])==24 and DESIGN['no_economic_or_target_change'] is True,'design drift')
require(DESIGN['maximum_objective_calls_per_chain']==250 and DESIGN['final_reserve_seconds']==900,'search budget drift')
cal=ROOT.parent/'purchase_rules_overnight_v1'
runner_rel='output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/run_psi.py'
inventory=read(cal/'inventory.json')['files']
require(runner_rel in inventory,'missing runner pin')
require(inventory[runner_rel]==RUNNER_SHA and sha(cal/'source'/runner_rel)==RUNNER_SHA,'runner source pin mismatch')
out={'schema':'fresh_calibration_native_collection_v1','remote_root':str(ROOT),'design_sha256':SUBMIT['design_sha256'],
 'job_id':SUBMIT['job_id'],'slots':[],'winners':{},'files':{}}
for slot,start in enumerate(DESIGN['starts']):
 arm='hard' if slot<12 else 'quarter';d=ROOT/'results'/f'slot_{slot}'
 require(start['slot']==slot and start['arm']==arm,'slot mapping drift')
 terminal=read(d/'launcher_terminal.json');search=read(d/'search/search_completed.json')
 post=read(d/'postcheck/completed.json');input_contract=read(d/'postcheck/input_contract.json')
 contracts={stage:read(d/f'{stage}_fresh_contract.json') for stage in ('search','postcheck')}
 require(terminal['slot']==slot and terminal['exit_code']==0,'launcher not successful: '+str(slot))
 require(terminal['maximum_objective_calls']==250 and terminal['final_reserve_seconds']==900,'launcher budget drift')
 require(terminal['deadline_epoch']<=DESIGN['absolute_deadline_epoch'],'deadline drift')
 require(search['objective_calls']<=250 and search['selected'] is not None,'search budget or selection failure')
 for stage,c in contracts.items():
  require(c['slot']==slot and c['arm']==arm and c['stage']==stage,'stage identity drift')
  require(c['design_sha256']==SUBMIT['design_sha256'],'stage design hash drift')
  require(c['runner_sha256']==RUNNER_SHA,'runner hash drift')
  require(c['reviewed_censor_sha256']==DESIGN['reviewed_censor_source_sha256'],'censor hash drift')
  require(c['original_plan_sha256']==DESIGN['original_plan_sha256'],'plan hash drift')
  require(c['target_contract_sha256']==DESIGN['target_fingerprint'] and c['weight_fingerprint']==DESIGN['weight_fingerprint'],'target or weight drift')
  require(c['parameters']==start['parameters'] and c['max_objective_calls']==250,'start or budget drift')
  require(c['deadline_epoch']==terminal['deadline_epoch'],'stage deadline drift')
 require(post['status']=='selected_numerically_verified','postcheck status failure')
 selected=post['selected'];checked=post['selected_postcheck']
 require(selected==search['selected'] and checked['status']=='passed','selected point drift')
 require(selected['base_loss']==checked['base_loss'] and selected['target_fit']==checked['base_target_fit'],'selected/postcheck mismatch')
 require(input_contract['purchase_rule']==arm and input_contract['owner_financed_share']==.8,'purchase rule drift')
 require(input_contract['base_target_contract_sha256']==DESIGN['target_fingerprint'] and input_contract['weight_contract_sha256']==DESIGN['weight_fingerprint'],'input fingerprint drift')
 base=d/'postcheck/selected_postcheck/phase_b_ge';full=base/'selected_root'
 target=rows(full/'target_fit.csv');params=rows(full/'parameters.csv');closure=read(full/'closure.json')
 require(len(target)==14 and len(params)==31,'14/31 table count failure')
 require(target==checked['base_target_fit'],'native target table mismatch')
 require(sum(x['role']=='scored' for x in target)==10,'scored target count failure')
 loss=sum(float(x['loss_contribution']) for x in target if x['loss_contribution'])
 require(close(loss,checked['weighted_loss']) and close(loss,selected['loss']),'loss arithmetic failure')
 bounds=DESIGN['full_approved_bounds'];free=set()
 for row in params:
  key=row['parameter']
  if key in bounds:
   free.add(key);require(close(row['lower'],bounds[key][0],1e-12) and close(row['upper'],bounds[key][1],1e-12),'bound drift: '+key)
   require(float(bounds[key][0])-1e-12<=float(row['estimate'])<=float(bounds[key][1])+1e-12,'estimate outside bound: '+key)
 require(free==set(bounds),'free parameter count failure')
 plots=list((full/'standard_diagnostics').glob('*.png'))
 require(closure['standard_plot_count']==17 and len(plots)==17,'17 standard plots missing')
 require(closure['population_scale']==1. and closure['normalized_population']==1.,'population normalization drift')
 require(abs(float(closure['renewal_residual']))<=1e-6 and abs(float(closure['absolute_housing_residual']))<=1e-6,'equilibrium closure failure')
 require((base/'selected_repeat/stage/solution_arrays.npz').is_file(),'selected native arrays missing')
 out['slots'].append({'slot':slot,'arm':arm,'loss':loss,'terminal':terminal,'search_completed':search,'postcheck_completed':post,
  'input_contract':input_contract,'search_contract':contracts['search'],'postcheck_contract':contracts['postcheck']})
 winner=out['winners'].get(arm)
 if winner is None or loss<winner['loss']:out['winners'][arm]={'slot':slot,'loss':loss}
for arm,w in out['winners'].items():
 slot=w['slot'];base=ROOT/'results'/f'slot_{slot}'/'postcheck/selected_postcheck/phase_b_ge'
 manifest={}
 for name in ('selected_root','selected_repeat','selected_repeat_final'):
  folder=base/name;require(folder.is_dir(),'missing selected packet '+name)
  for p in folder.rglob('*'):
   if p.is_file():manifest[str(p.relative_to(base))]={'sha256':sha(p),'bytes':p.stat().st_size}
 out['files'][arm]=manifest
print(json.dumps(out,allow_nan=False))
'''


def main() -> None:
    design = json.loads((HERE / "design.json").read_text())
    submit = json.loads((HERE / "submission_receipt.json").read_text())
    pins = json.loads((HERE / "source_sha256.json").read_text())
    for name, expected in pins.items():
        if digest(HERE / name) != expected:
            raise RuntimeError(f"local source pin drift: {name}")
    if submit["design_sha256"] != digest(HERE / "design.json"):
        raise RuntimeError("local submission/design drift")
    runner_path = HERE.parent / "run_psi.py"
    if not runner_path.is_file():
        raise RuntimeError("missing local runner pin")
    runner_sha = digest(runner_path)
    script = "ROOT_TEXT=" + repr(ROOT) + "\nDESIGN_TEXT=" + repr(json.dumps(design)) + "\nSUBMIT_TEXT=" + repr(json.dumps(submit)) + "\nPINS_TEXT=" + repr(json.dumps(pins)) + "\nRUNNER_SHA_TEXT=" + repr(runner_sha) + "\n" + REMOTE
    proc = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", "torch", "python3", "-"],
                          input=script, text=True, capture_output=True, timeout=120)
    if proc.returncode:
        raise RuntimeError("remote validation failed: " + proc.stderr[-1200:])
    snapshot = json.loads(proc.stdout)
    if len(snapshot["slots"]) != 24 or set(snapshot["winners"]) != {"hard", "quarter"}:
        raise RuntimeError("incomplete remote receipt")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "all24_postchecks.json").write_text(json.dumps(snapshot["slots"], indent=2) + "\n")
    for arm, winner in snapshot["winners"].items():
        slot = winner["slot"]
        base = f"{ROOT}/results/slot_{slot}/postcheck/selected_postcheck/phase_b_ge"
        destination = OUT / arm / f"slot_{slot}"
        destination.mkdir(parents=True, exist_ok=True)
        for name in ("selected_root", "selected_repeat", "selected_repeat_final"):
            subprocess.run(["rsync", "-a", "--timeout=120", f"torch:{base}/{name}/", str(destination / name) + "/"],
                           check=True, timeout=300)
        manifest = snapshot["files"][arm]
        for rel, expected in manifest.items():
            path = destination / rel
            if not path.is_file() or path.stat().st_size != expected["bytes"] or digest(path) != expected["sha256"]:
                raise RuntimeError(f"{arm} packet hash mismatch: {rel}")
        with (destination / "selected_root/target_fit.csv").open(newline="") as stream:
            if len(list(csv.DictReader(stream))) != 14:
                raise RuntimeError(f"{arm} collected target row count drift")
        with (destination / "selected_root/parameters.csv").open(newline="") as stream:
            if len(list(csv.DictReader(stream))) != 31:
                raise RuntimeError(f"{arm} collected parameter row count drift")
    receipt = {k: snapshot[k] for k in ("schema", "remote_root", "design_sha256", "job_id", "winners", "files")}
    receipt["postchecked_slots"] = len(snapshot["slots"])
    receipt["verified_local_hashes"] = True
    (OUT / "completed.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"postchecked_slots": 24, "winners": snapshot["winners"], "receipt": str(OUT / "completed.json")}))


if __name__ == "__main__":
    main()
