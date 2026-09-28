"""Prepare new bounded diagnostics from frozen overnight scientific inputs, Torch only."""
from pathlib import Path
import os,json,hashlib,copy,sys
assert sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit()
R=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26'); O=R/'output/model/fertility_identification_20260928/contract_v1';O.mkdir(parents=True,exist_ok=False)
read=lambda p:json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def pin(p):return dict(path=str(p),sha256=sha(p))
def canon(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')
src=R/'output/model/overnight_calibration_20260928/contract_v1/contract.json';c=read(src);assert sha(src)=='c83aaff1a90b1ba5bb0919151840e745e6816cb0a3a023ce746e750f77226c5d'
run=R/'output/model/overnight_calibration_20260928/gated_v1/search';selected=read(run/'selected.json')['selected']['block'];ep=read(run/'selected_export/block/export_receipt.json')
def ap(path):
 p=Path(path);return dict(case_path=str(p),**{n:pin(p/(n+'.json' if n=='receipt' else n+'.csv')) for n in ('receipt','target_fit','parameters')})
c['source_contract']=pin(src);c['schema']='e5f_fertility_identification_v1';c['files']['night_helper']=copy.deepcopy(c['files']['driver']);c['files']['driver']=pin(R/'code/model/tools/run_e5f_fertility_identification.py')
c['files']['identification_tests']=pin(R/'code/model/tools/test_run_e5f_fertility_identification.py');c['files']['identification_preparer']=pin(Path(__file__).resolve())
c['initial_point']=selected['point'];c['anchor']=ap(selected['case_path']);c['anchor']['repeats']=[ap(r['case_path']) for r in ep['repeats']]
base=read(c['lanes']['primary']['objective']['path']);lanes={}
for name,mult,factor in [('primary',1,None),('early10',10,None),('early100',100,None),('profile_half',10,.5),('profile_double',10,2),('profile_quadruple',10,4)]:
 obj=copy.deepcopy(base)
 for r in obj['target_rows']:
  if r['restriction_id']=='early_fertility':r['actual_weight']*=mult
 path=O/(name+'_objective.json');write(path,obj)
 lanes[name]=dict(objective=pin(path),canonical_sha256=canon(obj),target_weight_fingerprint=canon(obj['target_rows']),early_multiplier=mult,fixed={} if factor is None else {'kappa_fert_continuation':selected['point']['kappa_fert_continuation']*factor})
c['lanes']=lanes;c['jacobian_steps']={k:(.0005 if k=='beta_annual' else .02*abs(v)) for k,v in c['initial_point'].items()}
c['budget']=dict(workers=24,objective_cap_seconds=1800,max_objective_cases=538,max_search_cases=480,population=16,generations=4,total_seconds=21600,search_seconds=18000,repeat_seconds=21000)
c['seed']=2026092801;c['early_moment']='early_fertility'
alt=read(run/'initial_0487_identity/case/receipt.json')['point'];c['seed_points']=[selected['point'],alt];c['seed_artifacts']=[pin(Path(selected['case_path'])/'receipt.json'),pin(run/'initial_0487_identity/case/receipt.json')]
c['proposal_widths']=dict(H0=.4,beta_annual=.001,chi=.04,child_benefit_curvature=.08,delta_alpha_jump=.025,first_birth_fixed_cost=.15,kappa_fert=.05,kappa_fert_continuation=.15,theta0=.02,tenure_choice_kappa=.35)
c['author_scope']='Current full Jacobian, full early-weight reoptimization, and fixed continuation-scale profiles; no model or target changes. Diagnostic weighting only.'
write(O/'contract.json',c);write(O/'preparation_receipt.json',dict(contract_sha256=sha(O/'contract.json'),native_solves=0,maximum_objectives=538,scientific_source=pin(src)))
print(sha(O/'contract.json'))
