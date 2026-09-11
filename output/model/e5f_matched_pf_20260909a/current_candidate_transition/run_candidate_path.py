"""Bounded current-candidate diagnostic; terminal -> six-date smoke -> 28 dates.

No calibration, target changes or production promotion. Each stage retains its
own original hashed contract and complete numerical/diagnostic receipts.
"""
import argparse, csv, hashlib, json, os
from pathlib import Path
import subprocess, sys, time

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
INITIAL=Path('/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batches/extended_refinement_20260911/results/cases/selected_exact_repetitions')
PINS={
 'initial_checkpoint':{'path':str(INITIAL/'evaluation/raw/repetition_02/initial_state.pkl.gz'),'sha256':'738b9112f58d96f5bbbaf41ff31d7bf034927ddde71b430ebb020f53ff2b1c9c'},
 'initial_summary':{'path':str(INITIAL/'evaluation/raw/repetition_02/summary.json'),'sha256':'2026fb1eba211b2fb55a071ca3456668f8fa8fc68f1f8ffb2415d47d941084b3'},
 'initial_contract':{'path':str(INITIAL/'initial_contract.json'),'sha256':'d7e45f3ae2f694c4bcb2b209cb47f25b2653b57a4b357512793bf5057a78d552'}}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def pin(p):return dict(path=str(Path(p).resolve()),sha256=sha(p))
def read(p):return json.loads(Path(p).read_text())
def save(p,d):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');os.replace(q,p)

def run_stage(folder,name,c,driver):
 cp=folder/'contracts'/f'{name}.json';save(cp,c)
 out=folder/name
 command=[sys.executable,str(ROOT/'code/model/tools'/driver),'--contract',str(cp),'--contract-sha256',sha(cp),'--output',str(out)]
 save(folder/'latest_stage.json',dict(stage=name,status='running',command=command,start=time.time()))
 with (folder/f'{name}.log').open('w') as log:
  completed=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,timeout=c['seconds']+45)
 summary=read(out/'summary.json') if (out/'summary.json').exists() else {}
 save(folder/'latest_stage.json',dict(stage=name,exit_code=completed.returncode,summary=summary))
 if completed.returncode:raise RuntimeError(f'{name} did not pass; inspect saved numerical receipts before continuation')
 return out,cp,summary

def compare(out,folder):
 rows=list(csv.DictReader((out/'transition_path.csv').open()))
 data=list(csv.DictReader((HERE/'inputs/empirical_blocks.csv').open()))
 keyed={int(float(r['calendar_year'])):r for r in rows}
 first=float(keyed[2007]['birth_children_topcode_adjusted'])
 result=[]
 for d in data:
  year=int(d['decision_year']);m=float(keyed[year]['birth_children_topcode_adjusted'])/first;t=float(d['live_births_index_2008_2011'])
  result.append(dict(decision_year=year,birth_year_start=int(d['birth_year_start']),birth_year_end=int(d['birth_year_end']),target_birth_index=t,model_birth_index=m,gap_percentage_points=100*(m-t)))
 save(folder/f'{out.name}_birth_comparison.json',dict(rows=result,squared_index_error=sum((r['gap_percentage_points']/100)**2 for r in result[1:]),interpretation='Normalized birth-count path, not female TFR; national data versus model geography. Trial amplitude, not estimated shock.',production_eligible=False))

def main():
 a=argparse.ArgumentParser();a.add_argument('--delta',type=float,required=True);a.add_argument('--label',required=True);a.add_argument('--short-only',action='store_true');args=a.parse_args()
 if not (-.25<=args.delta<=0):raise ValueError('Explicit diagnostic decline bracket [-.25,0] required')
 folder=HERE/'results'/args.label;folder.mkdir(parents=True,exist_ok=False)
 for spec in PINS.values():
  if sha(spec['path'])!=spec['sha256']:raise ValueError('Initial parent hash mismatch')
 original=read(PINS['initial_contract']['path'])
 sources=dict(original['source_sha256'])
 for path in (ROOT/'code/model').rglob('*.py'):sources[str(path.relative_to(ROOT))]=sha(path)
 # The candidate driver independently rejects changes against original economic pins.
 common=dict(PINS,originating_source_commit=original['source_commit'],source_sha256=sources,psi_change_from_initial=args.delta)
 save(folder/'plan.json',dict(initial_candidate='r5_joint_09, exact repetition02',delta=args.delta,stages=['terminal8 mappings','six dates x8 mappings including replay','28 dates x8 mappings including replay'],maximum_bellman_calls=8+96+(0 if args.short_only else 448),budgets_seconds=[1800,1800]+([] if args.short_only else[7200]),prior_six_date_mapping_seconds_approximate=215,estimated_28_date_mapping_seconds_approximate=1003,stop='Fail stage on numerical/provenance gate; no downstream run from failed equilibrium',production_eligible=False,shock_estimated=False))
 try:
  c=read(HERE/'inputs/terminal_template.json');c.update(common,schema='e5f_candidate_terminal_v1')
  terminal,tc,ts=run_stage(folder,'terminal',c,'run_e5f_candidate_terminal.py')
  if not ts.get('endpoint_numerically_verified'):raise RuntimeError('Terminal not numerically verified')
  receipt=read(terminal/'root_receipt.json');pt=float(receipt['final']['prices'][0]);bt=float(receipt['final']['fiscal_values'][0])
  initial_price=float(read(PINS['initial_summary']['path'])['price'])
  for n in ([6] if args.short_only else [6,28]):
   c=read(HERE/'inputs/history_template.json');c.update(common,schema='e5f_candidate_history_root_v1',mode='actual_root',count=n,seconds=1800 if n==6 else 7200,root_seconds=1620 if n==6 else 7080)
   c['root_controls']['max_evaluations']=8
   c.update(terminal_checkpoint=pin(terminal/'terminal_state.pkl.gz'),terminal_summary=pin(terminal/'summary.json'),terminal_contract=pin(tc),terminal_root_receipt=pin(terminal/'root_receipt.json'))
   # Explicit numerical guesses only; all dates' pensions and prices are solved.
   c['initial_prices']=[initial_price+(pt-initial_price)*min(j/14,1) for j in range(n)]
   c['initial_pensions']=[3.51260051949617+(bt-3.51260051949617)*min(j/14,1) for j in range(n)]
   out,hc,hs=run_stage(folder,f'history_{n}',c,'run_e5f_candidate_history.py')
   compare(out,folder)
   if not hs.get('finite_horizon_market_fiscal_converged'):raise RuntimeError(f'{n}-date equilibrium incomplete; no horizon promotion')
  save(folder/'complete.json',dict(status='passed_diagnostic_stages',production_eligible=False,shock_estimated=False,horizon_verified=False))
 except Exception as exc:
  save(folder/'failure.json',dict(type=type(exc).__name__,message=str(exc),production_eligible=False));raise

if __name__=='__main__':main()
