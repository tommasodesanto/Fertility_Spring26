"""Execute on Torch; version the stationary and calendar probability correction."""
import ast
import copy
import hashlib
import json
import shutil
from pathlib import Path

B=Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913')
REL='code/model/intergen_eqscale_seq_optimized/solver.py'
OLD='3bd6782e2f3c27ebe533d54eb562bf74b0794a0e094c293c4766e82c39e4c6ae'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def canonical(d):return hashlib.sha256(json.dumps(d,sort_keys=True,separators=(',',':'),ensure_ascii=True,allow_nan=False).encode()).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')

for stem in ('corrected_initial_source','corrected_history_source'):
    previous=B/stem;current=B/(stem+'_v2');p=current/REL
    assert sha(previous/REL)==OLD
    source=(previous/REL).read_text();lines=source.splitlines(keepends=True)
    function=next(n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef) and n.name=='forward_distribution_markov_income')
    original=''.join(lines[function.lineno-1:function.end_lineno])
    assert 'normalized_probs' not in original
    token='                        for tn in range(nt):\n'
    assert original.count(token)==1
    insertion='''                        if tenure_probs is not None:
                            all_probs = np.asarray(
                                tenure_probs[:, to, id_, j, zz, nn, :, :], dtype=float
                            )
                            prob_sum = np.sum(all_probs, axis=-1)
                            normalized_probs = np.divide(
                                all_probs,
                                prob_sum[:, :, None],
                                out=np.zeros_like(all_probs),
                                where=prob_sum[:, :, None] > 0,
                            )
'''
    updated=original.replace(token,insertion+token)
    raw='pr = tenure_probs[:, to, id_, j, zz, nn, :, tn]'
    assert updated.count(raw)==1
    updated=updated.replace(raw,'pr = normalized_probs[:, :, tn]')
    replacement=''.join(lines[:function.lineno-1])+updated+''.join(lines[function.end_lineno:])
    if p.read_text()!=replacement:p.write_text(replacement)
    compile(p.read_text(),str(p),'exec')
    files=[x for x in (current/'code/model').rglob('*') if x.is_file() and '__pycache__' not in x.parts and x.suffix not in ('.pyc','.nbc','.nbi')]
    differences=[str(x.relative_to(current)) for x in files if sha(x)!=sha(previous/x.relative_to(current))]
    assert differences==[REL],differences
    save(B/(stem+'_v2_receipt.json'),dict(source_root=str(current),base_source_root=str(previous),only_changed_source=REL,source_files=len(files),source_sha256={str(x):sha(x) for x in files}))

NEW=sha(B/'corrected_initial_source_v2'/REL)
assert NEW==sha(B/'corrected_history_source_v2'/REL)
eq=read(B/'corrected_kernel_equivalence.json')
for pair in eq['pairs']:
    pair['initial']=pair['initial'].replace('/corrected_initial_source/','/corrected_initial_source_v2/')
    pair['history']=pair['history'].replace('/corrected_history_source/','/corrected_history_source_v2/')
    assert sha(pair['initial'])==sha(pair['history'])
    pair['sha256']=sha(pair['initial'])
save(B/'corrected_kernel_equivalence_v2.json',eq)
for stem in ('corrected_initial_source','corrected_history_source'):
    receipt_path=B/(stem+'_v2_receipt.json');receipt=read(receipt_path)
    receipt.update(kernel_equivalence=str(B/'corrected_kernel_equivalence_v2.json'),
                   kernel_equivalence_sha256=sha(B/'corrected_kernel_equivalence_v2.json'))
    save(receipt_path,receipt)

T=B/'corrected_initial_template_v6';OLDT=B/'corrected_initial_template_v4'
if not T.exists():shutil.copytree(OLDT,T,symlinks=True)
def relink(s):
    return s.replace(str(OLDT),str(T)).replace(str(B/'corrected_initial_source')+'/',str(B/'corrected_initial_source_v2')+'/').replace('"'+str(B/'corrected_initial_source')+'"','"'+str(B/'corrected_initial_source_v2')+'"').replace(OLD,NEW).replace('tenure_probability_mass_fix_v1','tenure_probability_mass_fix_v2').replace('tenure_probability_mass_conservation_v1','tenure_probability_mass_conservation_v2')
for rel in ['plan_capped_beta_099.json','run_capped_beta.py','run_profile.py','seed_case/initial_contract.json','seed_case/run_contract.json','inputs/working_contract.json','inputs/economic_source.json','inputs/observation_snapshot.json','inputs/run_scored_candidate.py']:
    (T/rel).write_text(relink((OLDT/rel).read_text()))
objective=read(T/'inputs/working_contract.json');prior=read(OLDT/'inputs/working_contract.json')
excluded={'contract_id','source_fingerprints','source_provenance'}
assert canonical({k:v for k,v in prior.items() if k not in excluded})=='b8400562a5a9ceedac83b9ad103515eeb9e0d594ce65c26fa93adde2be86ed4b'
objective['contract_id']=objective['contract_id'].replace('mass_fix_v1','mass_fix_v2')
for name,file in [('economic_source_manifest_7e872053','economic_source.json'),('observation_snapshot_manifest_70abd4a8','observation_snapshot.json')]:
    objective['source_fingerprints'][name]=canonical(read(T/'inputs'/file))
assert {k:v for k,v in objective.items() if k not in excluded}=={k:v for k,v in prior.items() if k not in excluded}
save(T/'inputs/working_contract.json',objective)
previous_fingerprint=canonical(prior);fingerprint=canonical(objective)
for rel in ('run_profile.py','run_capped_beta.py','inputs/run_scored_candidate.py'):
    p=T/rel;p.write_text(p.read_text().replace(previous_fingerprint,fingerprint));compile(p.read_text(),str(p),'exec')
initial=read(T/'seed_case/initial_contract.json')
initial['numerical_correction'].update(corrected_sha256=NEW,scope='Identical float64 tenure normalization in stationary Markov KFE and calendar Markov cohort advance',full_initial_replay_status='pending_v5')
save(T/'seed_case/initial_contract.json',initial)
run=read(T/'seed_case/run_contract.json')
run['working_objective']['canonical_sha256']=fingerprint
for key in ('working_objective','scorer','validator'):
    run[key]['sha256']=sha(run[key]['path'])
run['initial_solve_contract']=dict(path=str(T/'seed_case/initial_contract.json'),sha256=sha(T/'seed_case/initial_contract.json'))
run['wrapper_sha256']=sha(T/'inputs/run_scored_candidate.py')
save(T/'seed_case/run_contract.json',run)
plan=read(T/'plan_capped_beta_099.json')
plan['controller_sha256']=sha(T/'run_capped_beta.py');plan['resume_score_sha256']=sha(plan['resume_score_path'])
for name,digest in plan['file_sha256'].items():
    destination=Path(name) if Path(name).is_absolute() else T/name
    if not destination.exists():
        origin=Path(plan['seed_manifest_path']).parent/name
        assert sha(origin)==digest
        shutil.copyfile(origin,destination)
plan['file_sha256']={p:sha(Path(p) if Path(p).is_absolute() else T/p) for p in plan['file_sha256']}
save(T/'plan_capped_beta_099.json',plan)
script=(OLDT/'run.sh').read_text().replace(str(OLDT),str(T)).replace('corrected_initial_replay_v4','corrected_initial_replay_v6')
(T/'run.sh').write_text(script)
save(T/'source_objective_version.json',dict(status='explicit_corrected_source_objective_version',previous_objective_canonical_sha256=previous_fingerprint,corrected_objective_canonical_sha256=fingerprint,unchanged_target_weight_measurement_parameter_fields_sha256=canonical({k:v for k,v in objective.items() if k not in excluded}),changed_objective_fields=sorted(excluded),only_changed_numerical_source=REL,corrected_solver_sha256=NEW,canonical_json_method='ensure_ascii=True, allow_nan=False; sorted compact JSON',full_initial_replay_status='pending'))
print(json.dumps(dict(solver_sha256=NEW,objective_fingerprint=fingerprint,matching_kernel_files=len(eq['pairs']),template=str(T))))
