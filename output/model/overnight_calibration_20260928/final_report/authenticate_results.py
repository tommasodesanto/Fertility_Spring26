#!/usr/bin/env python3
"""Torch-only, stdlib final-artifact authentication. Never certifies visual review."""
import argparse, csv, hashlib, json, math, os
from collections import Counter
from pathlib import Path

EXPECTED='c83aaff1a90b1ba5bb0919151840e745e6816cb0a3a023ce746e750f77226c5d'
LOGICAL=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
BASE=Path('output/model/overnight_calibration_20260928')
def need(ok,why):
 if not ok:raise ValueError(why)
def read(p):return json.loads(Path(p).read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def canon(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def table(p,key,n):
 with Path(p).open(newline='') as f:rr=list(csv.DictReader(f))
 need(len(rr)==n and len({r[key] for r in rr})==n,'Wrong or duplicate row count: '+str(p));return rr
def csvwrite(p,rr):
 with Path(p).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')
def main(a):
 need(os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm only; checkpoint hashing is not local work')
 root=Path(a.project_root).resolve()
 def mapped(path):
  p=Path(path)
  if p.is_relative_to(LOGICAL):p=root/p.relative_to(LOGICAL)
  elif not p.is_absolute():p=root/p
  p=p.resolve();need(p.is_relative_to(root),'Path escapes physical project: '+str(p));return p
 run=mapped(BASE/'gated_v1/search');out=mapped(BASE/'gated_v1/final_review')
 need(not out.exists(),'Never overwrite review directory; move an explicitly failed prior review before retry')
 contract=mapped(BASE/'contract_v1/contract.json');need(sha(contract)==EXPECTED,'Contract SHA differs from approved night')
 c=read(contract);science=canon(c)
 for item in list(c['files'].values())+[c['source_manifest']]:need(sha(mapped(item['path']))==item['sha256'],'Source pin mismatch: '+item['path'])
 objs={}
 for lane,spec in c['lanes'].items():
  p=mapped(spec['objective']['path']);need(sha(p)==spec['objective']['sha256'],'Objective pin mismatch');obj=read(p)
  need(canon(obj)==spec['canonical_sha256'] and canon(obj['target_rows'])==spec['target_weight_fingerprint'],'Objective fingerprint mismatch');objs[lane]=obj
 need(set(objs)=={'primary','identity','block'},'Wrong lanes')
 final=read(run/'complete.json');need(final['status']=='bounded_search_complete','Search incomplete')
 need(final['contract_sha256']==EXPECTED and final['scientific_identity']==science,'Final identity mismatch')
 selected=final['selected'];key=final['common_primary_key'];best=selected[key]
 need(set(selected)>={'primary','identity','block'} and len(selected) in (3,4),'Wrong selection')
 records={r['case']:r for r in final['records']};need(len(records)==len(final['records']),'Duplicate cases')
 eligible=[r for r in records.values() if not r['case'].startswith('repeat_') and r['status']=='success']
 need(best['primary_rescore']==min(r['primary_rescore'] for r in eligible),'Not common-primary best')
 plotnames=sorted(c['standard_diagnostic_names']);need(len(plotnames)==len(set(plotnames))==17,'Plot contract mismatch')
 hashes={};lanes={};allfits=[];allparams=[];rescored_by_lane={};checkpoints={}
 primary={r['restriction_id']:r for r in objs['primary']['target_rows']}
 def check_case(row):
  need(row['status']=='success' and records.get(row['case'])==row,'Case not in successful final ledger')
  path=mapped(row['case_path']);need(path==run/row['case']/'case','Case path mismatch')
  receipt=read(path/'receipt.json');side=read(path.parent/'success.json');req=read(mapped(row['request_path']))
  need(req['id']==row['case'] and req['point']==row['point'] and req['lane']==row['lane'],'Request mismatch')
  need(req['contract_sha256']==EXPECTED and req['normalization_inputs']==c['normalization'],'Request contract mismatch')
  need(side['context']==req['context'],'Success context mismatch')
  ctx=req['context'];need(ctx['candidate_id']==row['case'] and ctx['point_sha256']==canon(row['point']) and ctx['contract_sha256']==EXPECTED and ctx['source_sha256']==c['files']['driver']['sha256'] and ctx['target_sha256']==c['lanes'][row['lane']]['objective']['sha256'],'Request source/target context mismatch')
  need(sha(path/'receipt.json')==row['receipt_sha256']==side['receipt_sha256'],'Receipt hash mismatch')
  lane=row['lane'];spec=c['lanes'][lane]
  need(receipt['status']=='verified_provisional_calibration_point','Case lacks passing receipt')
  need(receipt['contract_sha256']==EXPECTED and receipt['scientific_identity']==science,'Receipt identity mismatch')
  need(receipt['scientific_candidate_id']==canon(dict(science=science,lane=lane,point=row['point'],normalization=c['normalization'])),'Candidate identity mismatch')
  need(receipt['lane']==lane and receipt['point']==row['point'] and receipt['normalization_inputs']==c['normalization'],'Receipt parameters mismatch')
  need(receipt['target_system_sha256']==spec['objective']['sha256'] and receipt['target_weight_fingerprint']==spec['target_weight_fingerprint'] and receipt['source_manifest_sha256']==c['source_manifest']['sha256'],'Receipt source/target mismatch')
  h=sha(path/'initial_state.pkl.gz');checkpoints[row['case']]=h
  need(h==row['checkpoint_sha256']==side['checkpoint_sha256']==receipt['case_checkpoint_sha256'],'Checkpoint mismatch')
  norm=receipt['normalization'];need(norm['psi_child']>0 and norm['absolute_gap']<=5e-4 and norm['target']==2.1,'Normalization failed')
  ledger=read(path/'stationary_solves.json');need(len(ledger)==receipt['objective_stationary_solves']==norm['stationary_solves'] and 1<=len(ledger)<=23 and all(x['status']=='completed' for x in ledger),'Solve ledger failure')
  need(receipt['fiscal']['fiscal_gate'] and receipt['fiscal']['marginal_gate'] and receipt['household_budget']['budget_excess_mass']==0,'Fiscal/budget gate failure')
  fits=table(path/'target_fit.csv','moment',14);params=table(path/'parameters.csv','parameter',31)
  targets={r['restriction_id']:r for r in objs[lane]['target_rows']};need({r['moment'] for r in fits}==set(targets),'Target set mismatch')
  loss=0.;common=0.
  for r in fits:
   t=targets[r['moment']];gap=float(r['gap']);w=t['actual_weight'];need(float(r['target'])==t['target'] and math.isfinite(float(r['model'])) and math.isclose(float(r['model'])-float(r['target']),gap,rel_tol=1e-12,abs_tol=1e-12),'Target/gap mismatch')
   if w is None:need(r['weight']==r['loss_contribution']=='','Normalization weighted')
   else:need(float(r['weight'])==w and math.isclose(float(r['loss_contribution']),w*gap*gap,rel_tol=1e-12,abs_tol=1e-12),'Weight/contribution mismatch');loss+=float(r['loss_contribution'])
   common+=(primary[r['moment']]['actual_weight'] or 0)*gap*gap
  need(math.isclose(loss,receipt['loss'],rel_tol=1e-12,abs_tol=1e-12) and receipt['loss']==row['loss']==side['loss'] and math.isclose(common,row['primary_rescore'],rel_tol=1e-12,abs_tol=1e-12),'Loss mismatch')
  pm={r['parameter']:r for r in params};bounds={r['parameter']:r for r in objs[lane]['parameter_restrictions']}
  need(len(bounds)==10 and all(math.isfinite(float(r['estimate'])) for r in params),'Parameter table mismatch')
  for name,b in bounds.items():
   r=pm[name];v=float(r['estimate']);need(float(r['lower'])==b['lower']<=v<=b['upper']==float(r['upper']) and math.isclose(v,row['point'][name],rel_tol=2e-12,abs_tol=2e-12),'Parameter bounds/point mismatch')
  return path,fits,params
 for lane,original in selected.items():
  packet=run/'selected_export'/lane;export=read(packet/'export_receipt.json')
  need(export['contract_sha256']==EXPECTED and export['selected']==original,'Export source mismatch')
  reps=export['repeats'];need(len(reps)==2 and len({r['case'] for r in reps})==2,'Need two distinct repeats')
  source,fits,params=check_case(original)
  repeatpaths=[]
  for r in reps:
   need(r['case'].startswith('repeat_') and r['design']==lane and r['point']==original['point'] and r['lane']==original['lane'],'Repeat assignment mismatch')
   path,_,_=check_case(r);repeatpaths.append(path)
  need(mapped(export['graph_source'])==repeatpaths[0],'Graph source mismatch')
  for filename in ('target_fit.csv','parameters.csv'):
   h=sha(source/filename);need(all(sha(p/filename)==h for p in repeatpaths+[packet]),'Full table not byte identical: '+filename)
  need(sha(source/'receipt.json')==sha(packet/'receipt.json'),'Export receipt changed')
  for p in repeatpaths+[packet]:need(sorted(x.name for x in (p/'standard_diagnostics').glob('*.png'))==plotnames,'Missing/extra standard plots')
  for name in plotnames:
   h=sha(repeatpaths[0]/'standard_diagnostics'/name);need(all(sha(p/'standard_diagnostics'/name)==h for p in [repeatpaths[1],packet]),'Repeat/export plot mismatch: '+name)
  for path in [packet/n for n in ('target_fit.csv','parameters.csv','receipt.json','export_receipt.json')]+[packet/'standard_diagnostics'/n for n in plotnames]:hashes[str(path.relative_to(run))]=sha(path)
  rr=[]
  for r in fits:
   r=dict(r);w=primary[r['moment']]['actual_weight'];r['weight']='' if w is None else w;r['loss_contribution']='' if w is None else w*float(r['gap'])**2;rr.append(r)
  rescored_by_lane[lane]=rr;allfits.extend(dict(selected_lane=lane,**r) for r in rr);allparams.extend(dict(selected_lane=lane,**r) for r in params)
  lanes[lane]=dict(all14_31_byte_exact=True,all17_repeat_export_plot_hashes_exact=True,plots_visually_reviewed=False,selected_case=original['case'],repeat_cases=[r['case'] for r in reps])
 out.mkdir(parents=True)
 fitpath=out/'target_fit_primary_rescore.csv';csvwrite(fitpath,rescored_by_lane[key]);csvwrite(out/'all_lane_primary_fits.csv',allfits);csvwrite(out/'all_lane_parameters.csv',allparams)
 auth=dict(status='authenticated_final_export',visual_review_status='PENDING_LEAD_REVIEW',complete_sha256=sha(run/'complete.json'),contract_sha256=EXPECTED,common_primary_best=best,primary_fit_sha256=sha(fitpath),lanes=lanes,selected_files_sha256=hashes,checkpoint_sha256=checkpoints,model_solves=0)
 write(out/'authentication.proposed.json',auth)
 write(out/'comparison.json',dict(common_primary_key=key,selected=selected,search_counts=dict(Counter(r['status'] for r in records.values() if not r['case'].startswith('repeat_'))),repeat_count=sum(r['case'].startswith('repeat_') for r in records.values()),visual_review='pending',model_solves=0))
 print(json.dumps(dict(output=str(out),selected=list(selected),visual_review='pending',model_solves=0)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--project-root',required=True);main(p.parse_args())
