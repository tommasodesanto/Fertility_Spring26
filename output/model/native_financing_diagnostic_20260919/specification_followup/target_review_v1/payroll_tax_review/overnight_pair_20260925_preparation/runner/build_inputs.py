#!/usr/bin/env python3
"""Torch-only construction of a review-pending paired target and proposal bank."""
import json
from pathlib import Path

BASE=Path('/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a')
WORK=BASE/'nightpair_20260925_v1'
old=json.loads((BASE/'commute_calibration_20260924_v1/inputs/objective.json').read_text())
contract=json.loads((WORK/'inputs/proposed_common_contract.json').read_text())
old_rows={r['restriction_id']:r for r in old['target_rows']}
new_rows={r['id']:r for r in contract['target_rows']}
if set(old_rows)!=set(new_rows) or len(old_rows)!=13:
    raise RuntimeError('target row identity changed')
for rid,row in old_rows.items():
    new=new_rows[rid]
    row.update(target=new['target'],actual_weight=new['weight'],label=new['label'],
        definition=new['definition'],sample=new['sample'],
        model_observation=new['model_observation'],
        empirical_builder=new['source']['builder'],
        empirical_source_path=new['source']['path'],
        empirical_record_id=new['source']['record_id'],
        empirical_provenance_contract_id=new['source']['contract_id'],
        empirical_standard_error=new['standard_error'],
        uncertainty_status=new['uncertainty_status'],
        weight_status=new['weight_status'],
        warning=new['mapping_warning'])
old['parameter_restrictions']=contract['free_parameter_bounds']
old['objective_name']='paired_overnight_split_entry_20260925'
old['contract_id']='paired_overnight_split_entry_20260925_v1'
old['status']='new experimental target and source contract; pending code preflight and lead launch review'
old['source_provenance']={'paired_contract_path':str(WORK/'inputs/proposed_common_contract.json'),
                          'paired_contract_status':contract['status']}
old['experimental_changes']=contract['disclosures']
old['paired_contract']=contract
old['normalization_tolerance']=0.0005
assert old['external_restrictions']['theta1']==contract['fixed_settings']['theta1']['value']
assert old_rows['first_birth_rooms']['target']==1.465
assert sum(r['actual_weight'] is not None for r in old_rows.values())==12
assert {r['parameter']:(r['lower'],r['upper']) for r in old['parameter_restrictions']}['beta_annual']==(.94,.99)
out=WORK/'inputs/objective.json'
if out.exists(): raise RuntimeError('refuse overwrite of paired objective')
out.write_text(json.dumps(old,indent=2,sort_keys=True,allow_nan=False)+'\n')
print(out)
