"""Prepare a reproducible, full-bound Sobol exploration plan (zero solves)."""
import hashlib, json, math, sys
from pathlib import Path
from scipy.stats import qmc

ROOT=Path('/scratch/td2248/projects/estate_birth_global_search_20261004_v1')
PARENT=Path('/scratch/td2248/projects/estate_birth_binary_continuation_20261004_v2')
LOG={'chi','kappa_fert','kappa_fert_continuation','tenure_choice_kappa'}
SEED=20261004
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canonical(x): return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def main():
    parent=json.loads((PARENT/'control/starts.json').read_text())
    inventory=json.loads((PARENT/'inventory.json').read_text())
    assert sha(PARENT/'control/starts.json')==ROOT.joinpath('parent_starts.sha256').read_text().strip()
    assert parent['target_fingerprint']==inventory['target_fingerprint'] and parent['weight_fingerprint']==inventory['weight_fingerprint']
    bounds=parent['bounds']; names=sorted(bounds)
    assert len(names)==10 and LOG<=set(names)
    unit=qmc.Sobol(d=10,scramble=True,seed=SEED).random_base2(m=6)
    points=[]
    for u in unit:
        row={}
        for k,z in zip(names,u):
            lo,hi=map(float,bounds[k]);row[k]=math.exp(math.log(lo)+(math.log(hi)-math.log(lo))*float(z)) if k in LOG else lo+(hi-lo)*float(z)
        points.append(row)
    assert len(points)==64 and len({canonical(p) for p in points})==64
    assert all(0<min(unit[:,j]) and max(unit[:,j])<1 for j in range(10))
    assert all(len(set(int(v*64) for v in unit[:,j]))==64 for j in range(10)), 'Sobol marginal strata drift'
    incumbent=json.loads((ROOT/'incumbent.json').read_text())
    assert incumbent['status']=='selected_numerically_verified' and incumbent['native_loss']==21.275413361071312
    plan=dict(stage='global_exploration_v1',seed=SEED,sobol_m=6,unit_coordinates=unit.tolist(),points=points,
        transforms={k:('log' if k in LOG else 'linear') for k in names},bounds=bounds,
        incumbent_control=dict(parameters=incumbent['selected']['parameters'],native_loss=incumbent['native_loss'],source='recovery_chain_1_verified'),
        target_contract=parent['target_contract'],target_fingerprint=parent['target_fingerprint'],weight_fingerprint=parent['weight_fingerprint'],
        parent_inventory_sha256=sha(PARENT/'inventory.json'),parent_starts_sha256=sha(PARENT/'control/starts.json'),
        task_count=16,cases_per_task=4,task_wall_seconds=5400,case_budget_seconds=1200,stop_loss=13.,
        no_auto_extension=True,no_scientific_adoption=True)
    out=ROOT/'control/plan.json';assert not out.exists();out.write_text(json.dumps(plan,sort_keys=True,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status='prepared_zero_solves',plan_sha256=sha(out),points=64,task_count=16,seed=SEED,
                          transformed=sorted(LOG),marginal_strata_per_coordinate=64)))
if __name__=='__main__':main()
