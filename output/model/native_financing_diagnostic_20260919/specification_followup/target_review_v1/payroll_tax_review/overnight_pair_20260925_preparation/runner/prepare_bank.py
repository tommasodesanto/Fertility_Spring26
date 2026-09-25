#!/usr/bin/env python3
"""Build identical, static paired proposal coverage on Torch; no model imports."""
import hashlib,json,math,random
from pathlib import Path

BASE=Path("/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a")
WORK=BASE/"nightpair_20260925_v1"
PLAN=BASE/"utility_overnight_20260923_v1/results/production/B_floor/worker09_proposal16/plan.json"
COMMUTE_BEST=BASE/"commute_calibration_20260924_v1/results/run_001/best_so_far.json"
OBJECTIVE=WORK/"inputs/objective.json"
OUTPUT=WORK/"inputs/proposal_bank.json"
FREE=("H0","beta_annual","chi","first_birth_fixed_cost","h_P","kappa_fert","kappa_fert_continuation","theta0")

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):h.update(block)
    return h.hexdigest()

plan=json.loads(PLAN.read_text())
obj=json.loads(OBJECTIVE.read_text())
commute=json.loads(COMMUTE_BEST.read_text())
bounds={r["parameter"]:(float(r["lower"]),float(r["upper"])) for r in obj["parameter_restrictions"]}
widths=plan["search_config"]["coordinate_widths"]
assert set(bounds)==set(FREE) and set(FREE).issubset(widths)
assert bounds["beta_annual"]==(0.94,0.99) and bounds["h_P"][1]==2.3
seed_bfloor={k:plan["structural_parameters"][k] for k in FREE}
seed_commute={k:commute["point"][k] for k in FREE}
seeds={"b_floor_sep23":seed_bfloor,"commute_selected_sep24":seed_commute}
assert all(bounds[k][0]<=float(p[k])<=bounds[k][1] for p in seeds.values() for k in FREE)
workers={}
for slot in range(1,21):
    rng=random.Random(20260925+slot)
    seq=[]
    for index in range(1,19):
        # Each center receives both joint and one-coordinate draws.
        center=seed_bfloor if (index-1)//2%2==0 else seed_commute
        joint=index%2==0
        start=(index+slot-1)%len(FREE)
        names=[FREE[start]] if not joint else [FREE[(start+j)%len(FREE)] for j in range(3)]
        scale=.8 if slot<=7 else 1.6 if slot<=14 else 3.2
        point=dict(center)
        for name in names:
            lower,upper=bounds[name]
            for attempt in range(1000):
                value=center[name]+rng.gauss(0,float(widths[name])*scale)
                if lower<=value<=upper:break
            else: raise RuntimeError(f"bounded draw failed slot={slot} index={index} coordinate={name}")
            point[name]=value
        seq.append(point)
    workers[str(slot)]=seq
out=dict(schema="paired_fixed_proposal_bank_v1",seeds=seeds,workers=workers,
    proposal_rule="Same raw-coordinate Gaussian bank in both tax arms. Each worker has nine joint three-coordinate and nine one-coordinate draws, with both types at both centers; inherited Sep23 coordinate widths at scales 0.8/1.6/3.2 by worker stratum; bounded rejection, no transforms.",
    worker_rng_seed="20260925+slot",points_per_worker=18,workers_per_arm=20,
    plan_sha256=sha(PLAN),commute_best_sha256=sha(COMMUTE_BEST),objective_sha256=sha(OBJECTIVE))
OUTPUT.parent.mkdir(parents=True,exist_ok=True)
if OUTPUT.exists(): raise RuntimeError("proposal bank already exists; immutable")
OUTPUT.write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+"\n")
print(json.dumps(dict(path=str(OUTPUT),sha256=sha(OUTPUT),seeds=list(seeds),
    proposals_per_arm=sum(map(len,workers.values())),beta_max=max(p["beta_annual"] for ps in workers.values() for p in ps)),sort_keys=True))
