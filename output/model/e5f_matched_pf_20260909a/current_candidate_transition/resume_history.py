"""One bounded continuation of a verified but unfinished six-date root.

Uses the exact saved point/Jacobian and original economic inputs. Requires the
parent to finish first; numerical gates/checkpoint associations remain in the
existing candidate driver. A passed six-date root then starts28dates.
"""
import copy, json, sys
from pathlib import Path
import run_candidate_path as pipeline

FILES=dict(contract='contract.json',root_receipt='root_receipt.json',summary='summary.json',checkpoint='dated_2023.pkl.gz',dated_reproduction='dated_reproduction.json',root_evaluations='root_evaluations.json')

def continuation_contract(prior):
 c=pipeline.read(prior/'contract.json');c.pop('contract_sha256')
 root=pipeline.read(prior/'root_receipt.json')
 c.update(schema='e5f_candidate_history_continuation_v1',continuation={k:pipeline.pin(prior/v) for k,v in FILES.items()},initial_prices=root['best']['prices'],initial_pensions=root['best']['fiscal_values'])
 c['root_controls'].update(initial_jacobian=root['final_jacobian'],damping=root['final_damping'],max_evaluations=8)
 c.update(seconds=1800,root_seconds=1620)
 return c

def long_contract(c,root,terminal_root):
 c=copy.deepcopy(c);c.pop('continuation');c.update(schema='e5f_candidate_history_root_v1',count=28,seconds=7200,root_seconds=7080)
 c['root_controls']['initial_jacobian']=None
 pt=terminal_root['final']['prices'][0];bt=terminal_root['final']['fiscal_values'][0]
 for name,key,end in [('initial_prices','prices',pt),('initial_pensions','fiscal_values',bt)]:
  values=root['final'][key];c[name]=values+[values[-1]+(end-values[-1])*min(j/14,1) for j in range(1,23)]
 return c

def main():
 parent=pipeline.HERE/'results/delta_m005/history_6';s=pipeline.read(parent/'summary.json')
 if s.get('finite_horizon_market_fiscal_converged'):
  print('Parent six-date root passed; original pipeline owns its28-date stage. Nothing to resume.');return
 dest=pipeline.HERE/'recovery/delta_m005';dest.mkdir(parents=True,exist_ok=False)
 pipeline.save(dest/'plan.json',dict(maximum_additional_mappings=16,maximum_additional_bellman_calls=544,stage_budgets_seconds=[1800,7200],prior_mapping_seconds=170,stop='One continued six-date root;28dates only if it converges; fail on any provenance/replay/household gate.',production_eligible=False))
 c=continuation_contract(parent)
 sys.path[:0]=[str(pipeline.ROOT/'code/model/tools'),str(pipeline.ROOT/'code/model')]
 import run_e5f_candidate_history as driver
 driver.validate_contract(c);driver.load_continuation(c)
 out,cp,s=pipeline.run_stage(dest,'history_6_resumed',c,'run_e5f_candidate_history.py');pipeline.compare(out,dest)
 if not s.get('finite_horizon_market_fiscal_converged'):raise RuntimeError('One continued root remains incomplete; preserved for review')
 pipeline.save(dest/'accepted_history_6.json',dict(summary=pipeline.pin(out/'summary.json'),contract=pipeline.pin(cp),root_receipt=pipeline.pin(out/'root_receipt.json')))
 lc=long_contract(c,pipeline.read(out/'root_receipt.json'),pipeline.read(c['terminal_root_receipt']['path']))
 driver.validate_contract(lc)
 out,cp,s=pipeline.run_stage(dest,'history_28',lc,'run_e5f_candidate_history.py');pipeline.compare(out,dest)
 pipeline.save(dest/'completion.json',dict(summary=s,production_eligible=False,horizon_verified=False,shock_fitted=False))

if __name__=='__main__':main()
