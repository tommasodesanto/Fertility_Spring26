"""Derive a bounded 36-point completion stage from the pinned global stage."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OLD = ROOT / 'output/model/experiments/birth_count_choice/estate_a_global_search_20261004_v1/deployment'
OUT = ROOT / 'output/model/experiments/birth_count_choice/estate_a_global_completion_20261005_v1/deployment'
ORIG_SHA = '44ffc96620deab7577b8ecaaa744e41782b9294e643d04d2d5bafa7a7ce7e01b'
ORIG_MANIFEST_SHA = '0e6801e3126b0c613be733705327bab5c519cd9ec93f7bc80fda441a66df1cf1'
CHUNKS = [[1,2,3,5],[6,7,9,10],[11,13,14,15],[18,19,22,23],
          [29,30,31,34],[35,38,39,42],[43,46,47,49],[50,51,53,54],[55,61,62,63]]
COMPLETED = [16,20,24,25,26,27,32,36,40,44,56,57,58]
FATAL_ATTEMPTED = [0,4,8,12,17,21,28,33,37,41,45,48,52,59,60]

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def replace_one(text, old, new):
    assert text.count(old) == 1, (old, text.count(old))
    return text.replace(old, new)
def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value)

def main():
    assert sha(OLD/'plan.json') == ORIG_SHA
    assert sha(OLD/'stage_manifest.json') == ORIG_MANIFEST_SHA
    assert len(COMPLETED)==13 and len(FATAL_ATTEMPTED)==15
    assert sorted(COMPLETED+FATAL_ATTEMPTED+sum(CHUNKS,[])) == list(range(64))
    oldplan=json.loads((OLD/'plan.json').read_text())
    assert oldplan['stage']=='global_exploration_v1' and len(oldplan['points'])==64
    plan=dict(oldplan)
    plan.update(stage='global_completion_v1',task_count=9,cases_per_task=4,
                chunks=CHUNKS,original_plan_sha256=ORIG_SHA,
                original_stage_manifest_sha256=ORIG_MANIFEST_SHA,
                previous_completed_indices=COMPLETED,previous_fatal_attempted_indices=FATAL_ATTEMPTED,
                absolute_cutoff_utc='2026-10-05T18:00:00Z',no_auto_retry=True)
    write(OUT/'control/original_plan.json',(OLD/'plan.json').read_text())
    write(OUT/'control/plan.json',json.dumps(plan,indent=2,sort_keys=True)+'\n')
    for name in ('incumbent.json','parent_starts.sha256'):
        (OUT/name).write_bytes((OLD/name).read_bytes())
    driver=(HERE.parent/'global_search_v1/explore.py').read_text()
    driver=replace_one(driver,"require(0<=z.task<16,'Task index outside 0..15')",
                       "require(0<=z.task<9,'Task index outside 0..8')")
    driver=replace_one(driver,"require(plan['stage']=='global_exploration_v1' and len(plan['points'])==64 and plan['task_count']==16 and plan['cases_per_task']==4,'Plan cardinality drift')",
                       "require(plan['stage']=='global_completion_v1' and len(plan['points'])==64 and plan['task_count']==9 and plan['cases_per_task']==4 and plan['chunks']=="+repr(CHUNKS)+",'Completion plan drift')\n    require(sha(DEPLOY/'control/original_plan.json')==plan['original_plan_sha256']=='"+ORIG_SHA+"','Original Sobol plan drift')")
    driver=replace_one(driver,"first_case_index=4*z.task,case_count=2 if z.mode=='smoke' else 1 if z.mode=='preflight' else 4,",
                       "original_sobol_indices=plan['chunks'][z.task],case_count=2 if z.mode=='smoke' else 1 if z.mode=='preflight' else 4,")
    driver=replace_one(driver,"indexes=[0,1] if z.mode=='smoke' else [-1] if z.mode=='preflight' else list(range(4*z.task,4*z.task+4))",
                       "indexes=[1,2] if z.mode=='smoke' else [-1] if z.mode=='preflight' else plan['chunks'][z.task]")
    classify="""def classify_native_runtime(exc):
    message=str(exc)
    base=dict(reason=message,error_type=type(exc).__name__,lifecycle_solves=None)
    if message=='native GE acceptance failed: uncomputed_bounded_budget':
        return dict(base,status='budget_exhausted',rejection_kind='native_solve_cap')
    if message=='native GE acceptance failed: uncomputed_price_unbracketed':
        return dict(base,status='inadmissible_numerical',rejection_kind='price_unbracketed_diagnostic_caps')
    if type(exc).__name__=='InheritedDistributionInfeasible' and getattr(exc,'classification',None)=='inherited_distribution_infeasible' and hasattr(exc,'audit'):
        return dict(base,status='inadmissible_numerical',rejection_kind='typed_inherited_distribution_gate')
    return None

"""
    driver=replace_one(driver,'def main():\n',classify+'def main():\n')
    driver=replace_one(driver,"if str(exc)=='native GE acceptance failed: uncomputed_bounded_budget':\n                    result=dict(status='budget_exhausted',reason=str(exc),error_type=type(exc).__name__,\n                                rejection_kind='native_solve_cap',lifecycle_solves=None)\n                else:raise",
                       "result=classify_native_runtime(exc)\n                if result is None:raise")
    write(OUT/'explore.py',driver)
    oldlaunch=(HERE.parent/'global_search_v1/launch_torch.sh').read_text()
    launch=oldlaunch.replace('estate_birth_global_search_20261004_v1','estate_birth_global_completion_20261005_v1')
    launch=replace_one(launch,'^(0|[1-9]|1[0-5])$','^(0|[1-8])$')
    launch=replace_one(launch,'deadline_epoch=$((start_epoch+wall_seconds))',
                       'deadline_epoch=$((start_epoch+wall_seconds))\nif [[ "$mode" == production && "$deadline_epoch" -gt 1791223200 ]]; then deadline_epoch=1791223200; fi')
    gate='''# Keep recovery successors from overlapping the dependent global stage.
if [[ "$mode" == production ]]; then
  while true; do
    queue_snapshot=$(squeue -h -u "$USER" -o '%j %T') || { echo 'Recovery queue query failed; refusing native cases' >&2; exit 4; }
    if ! awk '$1 ~ /^(estatebirth|softtiming)$/ && $2 ~ /^(PENDING|RUNNING|COMPLETING)$/ {found=1} END {exit !found}' <<< "$queue_snapshot"; then break; fi
    "$python" - "$out" <<'PYWAIT'
import json,sys,time
from pathlib import Path
Path(sys.argv[1],'queue_wait.json').write_text(json.dumps(dict(status='waiting_for_overnight_recovery_successors',epoch=time.time()))+'\\n')
PYWAIT
    if (( $(date +%s) >= deadline_epoch - 300 )); then
      "$python" - "$out" <<'PYDEFER'
import json,sys,time
from pathlib import Path
Path(sys.argv[1],'deferred_concurrency.json').write_text(json.dumps(dict(status='deferred_no_native_cases',epoch=time.time()))+'\\n')
PYDEFER
      exit 0
    fi
    sleep 60
  done
fi
'''
    launch=replace_one(launch,'export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib\n',
                       gate+'export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib\n')
    write(OUT/'launch_torch.sh',launch)
    verify=(HERE.parent/'global_search_v1/verify_stage.py').read_text().replace('estate_birth_global_search_20261004_v1','estate_birth_global_completion_20261005_v1')
    verify=replace_one(verify,"inv=json.loads((root/'inventory.json').read_text())",
                       "inv=json.loads((root/'inventory.json').read_text())\n assert sha(root/'control/plan.json')==m['completion_plan_sha256']\n assert sha(root/'control/original_plan.json')==m['original_plan_sha256']")
    write(OUT/'verify_stage.py',verify)
    oldmanifest=json.loads((OLD/'stage_manifest.json').read_text())
    oldmanifest.update(stage='estate_a_global_completion_20261005_v1',
                       entrypoints={n:sha(OUT/n) for n in ['explore.py','verify_stage.py','launch_torch.sh']},
                       original_plan_sha256=ORIG_SHA,original_stage_manifest_sha256=ORIG_MANIFEST_SHA,
                       completion_plan_sha256=sha(OUT/'control/plan.json'),
                       no_auto_retry=True,no_auto_extension=True)
    write(OUT/'stage_manifest.json',json.dumps(oldmanifest,indent=2,sort_keys=True)+'\n')
    print(json.dumps(dict(plan_sha256=sha(OUT/'control/plan.json'),manifest_sha256=sha(OUT/'stage_manifest.json'),
                          indexes=sum(CHUNKS,[]),original_stage_manifest_sha256=sha(OLD/'stage_manifest.json'))))
if __name__=='__main__': main()
