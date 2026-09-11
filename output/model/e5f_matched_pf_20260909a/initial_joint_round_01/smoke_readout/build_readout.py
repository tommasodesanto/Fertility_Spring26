"""Deterministic readout of saved observations; no model or empirical solves."""
from pathlib import Path
import copy,csv,datetime,hashlib,importlib.util,json,shutil
ROOT=Path(__file__).resolve().parent
ROUND=ROOT.parent
BASE=ROUND.parent/'initial_fit_readout'
def read(p):return json.loads(p.read_text())
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,data):p.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
def csv_write(p,rows,fields=None):
 with p.open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=fields or list(rows[0]),extrasaction='ignore');writer.writeheader();writer.writerows(rows)
def lookup(data,key):
 for name in key.split('.'):data=data[name]
 return data
def fit_rows(prototype,summary,early,projection):
 out=[]
 for row in prototype:
  x=copy.deepcopy(row);where=x['model_observation'].replace('fertility.uniform_birth_time.','fertility.'+projection+'.')
  value=lookup(summary,where) if where.startswith('final.') else lookup(early,where)
  x.update(model=value,gap=None if value is None else value-x['target'],model_available=value is not None,model_observation=where,actual_weight=None,loss_contribution=None,calibrated_smm=False)
  x['model_source_path']=str(ROOT/'collected/smoke'/('summary.json' if where.startswith('final.') else 'repetition_02/early_measurement.json'))
  x['model_source_json_location']=where
  x['model_checkpoint_sha256']=summary['final']['checkpoint_sha256'];out.append(x)
 return out

def main():
 baseline=read(BASE/'summary.json');smoke=ROOT/'collected/smoke';summary=read(smoke/'summary.json');early=read(smoke/'repetition_02/early_measurement.json');contract=read(ROUND/'smoke_contract.json')
 assert digest(BASE/'target_fit.csv')==digest(ROUND/'inputs/target_fit.csv')
 assert digest(BASE/'parameters.csv')==digest(ROUND/'inputs/parameters.csv')
 assert digest(BASE/'reference_precision_options.csv')==digest(ROUND/'inputs/reference_precision_options.csv')
 spec=importlib.util.spec_from_file_location('panel_validator',ROUND/'inputs/panel_validator.py');v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
 for rep in (1,2):
  p=smoke/f'repetition_{rep:02d}';final=read(p/'summary.json');obs=read(p/'early_measurement.json');ges=read(p/'stationary_solves.json')
  v.validate_summary(dict(summary,repetitions=1,stationary_solves=len(ges),final=final),contract,obs,list(csv.DictReader((p/'parameters.csv').open())),ges)
 assert read(smoke/'repetition_01/early_measurement.json')==early
 market=v.validate_market(read(smoke/'repetition_02/market_quantity_units.json'))
 fits={projection:fit_rows(baseline['target_fit'],summary,early,projection) for projection in ('uniform_birth_time','constant_post_cell')}
 validation=fit_rows(baseline['validation_fit'],summary,early,'uniform_birth_time')
 assert len(fits['uniform_birth_time'])==13 and len(validation)==2
 assert [r['restriction_id'] for r in fits['uniform_birth_time'] if r['model'] is None]==['recent_parent_ownership']
 targets=list(csv.DictReader((ROUND/'inputs/target_fit.csv').open()))
 for prototype,row in zip(targets,baseline['target_fit']):
  assert prototype['restriction_id']==row['restriction_id'] and float(prototype['target'])==row['target'] and prototype['model_observation']==row['model_observation']
 fields=list(targets[0])+['model_checkpoint_sha256']
 csv_write(ROOT/'target_fit.csv',fits['uniform_birth_time'],fields)
 csv_write(ROOT/'target_fit_constant_post_cell.csv',fits['constant_post_cell'],fields)
 csv_write(ROOT/'validation_fit.csv',validation,fields)
 projection_rows=[]
 for row in baseline['cps_projection_sensitivity']:
  x=copy.deepcopy(row);source=next(r for r in fits[x['projection']] if r['restriction_id']==x['restriction_id'])
  x.update(model=source['model'],gap=source['gap'],actual_weight=None,loss_contribution=None,source_path=source['model_source_path']);projection_rows.append(x)
 csv_write(ROOT/'cps_projection_sensitivity.csv',projection_rows)
 compare=[]
 base_projection={(x['restriction_id'],x['projection']):x['model'] for x in baseline['cps_projection_sensitivity']}
 base_targets={x['restriction_id']:x for x in baseline['target_fit']}
 for projection,rows in fits.items():
  for row in rows:
   old=base_projection.get((row['restriction_id'],projection),base_targets[row['restriction_id']]['model']);new=row['model']
   compare.append({'restriction_id':row['restriction_id'],'label':row['label'],'projection':projection,'target':row['target'],'baseline_model':old,'baseline_gap':None if old is None else old-row['target'],'smoke_model':new,'smoke_gap':row['gap'],'model_change':None if new is None or old is None else new-old,'actual_weight':None,'loss_contribution':None,'model_observation':row['model_observation']})
 csv_write(ROOT/'comparison_target_fit.csv',compare)
 vcompare=[]
 for old,new in zip(baseline['validation_fit'],validation):
  vcompare.append({'restriction_id':new['restriction_id'],'target':new['target'],'baseline_model':old['model'],'baseline_gap':old['gap'],'smoke_model':new['model'],'smoke_gap':new['gap'],'model_change':new['model']-old['model'],'actual_weight':None,'loss_contribution':None})
 csv_write(ROOT/'comparison_validation_fit.csv',vcompare)
 shutil.copyfile(smoke/'repetition_02/parameters.csv',ROOT/'parameters.csv')
 parameters=list(csv.DictReader((ROOT/'parameters.csv').open()));assert len(parameters)==17
 oldparams={r['parameter']:r for r in csv.DictReader((ROUND/'inputs/parameters.csv').open())}
 pcompare=[dict(r,baseline_estimate=oldparams[r['parameter']]['estimate'],estimate_change=float(r['estimate'])-float(oldparams[r['parameter']]['estimate'])) for r in parameters]
 csv_write(ROOT/'comparison_parameters.csv',pcompare)
 shutil.copyfile(ROUND/'inputs/reference_precision_options.csv',ROOT/'reference_precision_options.csv')
 shutil.copyfile(BASE/'summary.json',ROOT/'baseline_readout_summary.json');shutil.copyfile(BASE/'validation_fit.csv',ROOT/'baseline_validation_fit.csv')
 result={'schema':'initial_joint_smoke_readout_v1','created_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'job_id':'17360699','status':'complete_saved_observation_readout','calibrated_smm':False,'actual_weight_profile':None,'summed_loss':None,'target_or_measurement_changes':False,'checkpoint_sha256':summary['final']['checkpoint_sha256'],'restriction_count':13,'validation_count':2,'parameter_count':17,'target_fit':fits['uniform_birth_time'],'target_fit_constant_post_cell':fits['constant_post_cell'],'validation_fit':validation,'cps_projection_sensitivity':projection_rows,'parameters':parameters,'comparison_target_fit':compare,'comparison_validation_fit':vcompare,'elapsed_seconds':summary['elapsed_seconds'],'stationary_solves':summary['stationary_solves'],'repetitions':2,'exact_early_equality':True,'source_contract_sha256':digest(ROUND/'smoke_contract.json'),'baseline_readout_sha256':digest(BASE/'summary.json'),'pinned_target_csv_sha256':digest(ROUND/'inputs/target_fit.csv'),'local_numeric_gates_verified':True,'market_residual_reconstructed':market,'normalization':summary['final']['normalization'],'fiscal':summary['final']['fiscal'],'household_budget':summary['final']['household_budget'],'operator_gates':summary['final']['operator_gates'],'policy_array_gates':summary['final']['policy_array_gates'],'unavailable_restrictions':['recent_parent_ownership'],'measurement_warnings_preserved':baseline['mapper_warnings'],'remote_gate_status':'awaiting_lead_receipt_path_correction','new_model_or_empirical_solves':0}
 gatepath=ROOT/'collected/control/smoke_gate_receipt.json'
 if gatepath.exists():
  gate=read(gatepath);assert gate['status']=='verified' and gate['source_pins_verified']==634
  result.update(remote_gate_status='verified',remote_gate_receipt=str(gatepath),remote_gate_receipt_sha256=digest(gatepath),reviewed_run_plan_sha256=gate['plan_sha256'],source_pins_verified=634,remote_checkpoint_graph_and_light_artifact_hashes_verified=len(gate['artifact_sha256']))
 write(ROOT/'summary.json',result)
 print(json.dumps({'status':result['status'],'rows':13,'validations':2,'parameters':17,'both_projections':True,'null_actual_weights':True,'numeric_gates_passed':True,'market_residual':market,'checkpoint_sha256':result['checkpoint_sha256']}))
if __name__=='__main__':main()
