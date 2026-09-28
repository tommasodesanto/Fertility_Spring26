from pathlib import Path
import json,hashlib,subprocess,time,collections
root=Path('/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project');orig='/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26';b=root/'output/model/fertility_identification_20260928'
state=subprocess.check_output(['sacct','-j','18716710','-X','-n','--format=State'],text=True).strip();assert state=='FAILED',state
old=json.loads((b/'run_v1/complete.json').read_text());assert old['status']=='incomplete_or_fatal_stop';assert time.time()<old['clock']['search_cutoff'];assert old['clock']['end']==1790630594.4615145
paths=[b/'run_v1/complete.json',b/'approval_v1.json',root/'code/model/tools/run_e5f_fertility_identification_resume.py',root/'code/model/tools/test_run_e5f_fertility_identification_resume.py',b/'resume_v1.sh'];overlay={}
for r in old['records']:
 paths.append(Path(r['request_path'].replace(orig,str(root))))
 if r['status']=='fatal':
  f=paths[-1].parent/r['case']/'failure.json';paths.append(f);e=json.loads(f.read_text());assert e['error_type']=='NonpositiveNormalizedBenefit' and e['classifier_error']=='candidate identity/stage required' and e==r['error'];overlay[r['case']]={'status':'inadmissible','reason':'verified_nonpositive_benefit_stage_label_only'}
pins=[{'path':str(p).replace(str(root),orig),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths]
m={'status':'lead_reviewed_resume','contract_sha256':old['contract_sha256'],'old_job_terminal':True,'old_job':'18716710','old_complete':str(b/'run_v1/complete.json').replace(str(root),orig),'original_approval':str(b/'approval_v1.json').replace(str(root),orig),'pins':pins,'rejection_overlays':overlay,'clock':old['clock'],'original_status_counts':dict(collections.Counter(r['status'] for r in old['records'])),'tests':'10 synthetic tests18721444 PASS; actual runtime gate18721487 PASS; independent source review complete; original failures untouched'}
f=b/'resume_manifest_v1.json';assert not f.exists();f.write_text(json.dumps(m,indent=2));print('Manifest prepared',m['original_status_counts'])
