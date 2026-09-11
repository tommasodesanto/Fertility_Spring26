"""Build the lead-specified frozen working contract, without scoring any case.

Reads only recorded decisions, target/provenance artifacts and saved observations.
Runs the existing pure scorer tests; never creates a model certification receipt,
changes an input packet, estimates weights or launches a model/cluster job.
"""
from pathlib import Path
import copy,csv,datetime,hashlib,json,math,subprocess,sys
import score_initial as score

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent

def read(path):return json.loads(path.read_text())
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def csv_rows(path):
 with path.open() as stream:return list(csv.DictReader(stream))
def require(ok,reason):
 if not ok:raise ValueError(reason)
def numeric(value):
 result=float(value);require(math.isfinite(result),'Nonfinite numeric input');return result
def lookup(packet,path):
 for key in path.split('.'):
  require(isinstance(packet,dict) and key in packet,'Missing model path: '+path);packet=packet[key]
 require(type(packet) in (int,float) and math.isfinite(packet),'Nonfinite model path: '+path)
 return packet

def build(*, replace_unreviewed=False):
 paths={
  'lead_decision':ROOT/'lead_methodological_decision.json',
  'working_weights':ROOT/'working_weights.csv',
  'target_fit':BASE/'initial_fit_readout/target_fit.csv',
  'parameter_table':BASE/'initial_fit_readout/parameters.csv',
  'fertility_provenance':BASE/'design_research/fertility_contract/fertility_target_contract.json',
  'housing_wealth_provenance':BASE/'design_research/observer_contract/observer_contract.json',
  'fertility_availability':BASE/'parameter_target_audit/fertility/fertility_availability.json',
  'economic_manifest_record':BASE/'initial_joint_round_01/preparation_metadata.json',
  'observer_manifest_record':BASE/'initial_fit_readout/recent_parent_probe/contract.json',
  'baseline_early_measurement':BASE/'initial_sensitivity_panel/collected/baseline/repetition_01/early_measurement.json',
  'joint_early_measurement':BASE/'initial_joint_round_01/smoke_readout/collected/smoke/repetition_02/early_measurement.json',
  'recent_panel_record':BASE/'initial_fit_readout/recent_parent_probe/completed_17362130/recent_parent_panel.json',
  'recent_joint_record':BASE/'initial_fit_readout/recent_parent_probe/completed_17362130/recent_parent_joint.json',
  'scorer':ROOT/'score_initial.py',
  'scorer_tests':ROOT/'test_score_initial.py',
  'builder':Path(__file__).resolve(),
  'inherited_parameter_gate':BASE/'initial_joint_round_01/inputs/panel_validator.py',
 }
 for name,path in paths.items():require(path.is_file(),'Missing input artifact: '+name+' '+str(path))
 decision=read(paths['lead_decision']);weights=csv_rows(paths['working_weights']);targets=csv_rows(paths['target_fit'])
 require(decision['primary_cps_projection']=='uniform_birth_time','Primary CPS decision changed')
 require(decision['structural_coordinates']==9 and decision['target_restrictions']==13 and decision['scored_restrictions']==12,'Lead design count changed')
 require(len(weights)==len(targets)==len(decision['rows'])==13,'Incomplete target/weight decision')
 weightmap={x['restriction_id']:x for x in weights};leadrows={x['restriction_id']:x for x in decision['rows']}
 require(len(weightmap)==len(leadrows)==13,'Duplicate target/weight decision')
 require([x['restriction_id'] for x in weights]==[x['restriction_id'] for x in targets],'Target/weight order mismatch')
 fert=read(paths['fertility_provenance']);housing=read(paths['housing_wealth_provenance']);availability=read(paths['fertility_availability'])
 fertility_records={x['id']:x for x in fert['candidate_moments']};housing_records={x['id']:x for x in housing['records']}
 fertility_indexes={x['id']:i for i,x in enumerate(fert['candidate_moments'])};housing_indexes={x['id']:i for i,x in enumerate(housing['records'])}
 cps_windows=[x['window'] for x in availability['cps_moments'] if x['window']=='2004+2006 pooled weighted records'];require(len(cps_windows)==1,'Authoritative pooled CPS window missing')
 require(len(fertility_records)==len(fert['candidate_moments']) and len(housing_records)==len(housing['records']),'Duplicate authoritative record IDs')
 economic=read(paths['economic_manifest_record']);observer=read(paths['observer_manifest_record'])
 economic_pins=economic['source_pins'];observer_pins=observer['source_sha256']
 require(economic['source_commit']=='7e872053' and len(economic_pins)==634,'Economic source identity changed')
 require(observer['source_commit']=='70abd4a80df979903b05d361f47ac2bf5aab2c23' and len(observer_pins)==641,'Observer snapshot identity changed')
 require(score.fingerprint(economic_pins)==economic['source_manifest_sha256'],'Economic manifest fingerprint mismatch')
 require(score.fingerprint(observer_pins)==observer['source_fingerprint'],'Observer manifest fingerprint mismatch')
 require(all(observer_pins.get(k)==v for k,v in economic_pins.items()),'Inherited economic source changed in observer snapshot')
 additions={k:v for k,v in observer_pins.items() if k not in economic_pins};require(len(additions)==7,'Observer addition count changed')
 source_fingerprints={'economic_source_manifest_7e872053':economic['source_manifest_sha256'],
                      'observation_snapshot_manifest_70abd4a8':observer['source_fingerprint']}
 for name in ('scorer','target_fit','working_weights','lead_decision','parameter_table','fertility_provenance','housing_wealth_provenance'):
  source_fingerprints[name+'_file_sha256']=digest(paths[name])
 baseline=read(paths['baseline_early_measurement']);joint=read(paths['joint_early_measurement'])
 panel=read(paths['recent_panel_record']);recentjoint=read(paths['recent_joint_record'])
 baseline_recent=[x for x in panel['rows'] if x['case_id']=='baseline'];require(len(baseline_recent)==1,'Missing/duplicate recent baseline record')
 recent_evidence=[]
 for name,row in (('baseline',baseline_recent[0]),('joint_smoke',recentjoint['row'])):
  path=Path(row['observation_path']);require(path.is_file(),'Missing recent-parent packet: '+str(path))
  require(digest(path)==row['observation_sha256'],'Recent-parent observation fingerprint mismatch: '+name)
  packet=read(path)
  require(packet['observer_id']==decision['recent_parent_observer']==score.RECENT_OBSERVER,'Recent observer name mismatch')
  require(packet['moment']==decision['recent_parent_model_moment']==score.RECENT_MOMENT,'Recent model name mismatch')
  require(packet['available'] is True and packet['model_value']==row['model_value'],'Recent observed value mismatch')
  require(packet['metadata']['policy_input_provenance']['checkpoint_sha256']==row['checkpoint_sha256'],'Recent checkpoint identity mismatch')
  require(packet['production_eligible'] is False and packet['target_contract_activated'] is False,'Original recent observer flags changed')
  paths['recent_'+name+'_observation']=path
  recent_evidence.append({'case_id':row['case_id'],'observation_path':str(path),'observation_sha256':digest(path),'checkpoint_sha256':row['checkpoint_sha256'],'observer_id':packet['observer_id'],'moment':packet['moment'],'original_packet_flags_preserved':True})
 approximation={'approximation_id':decision['name']+'__recent_parent_synchronized_proxy',
                'maintained':True,'exact_acs':False,'declaration':decision['recent_parent_approximation'],
                'observer_id':decision['recent_parent_observer'],'moment':decision['recent_parent_model_moment'],
                'snapshot':'synchronized_post_fertility_snapshot','age_projection':'uniform_within_age_cell',
                'diagnostic_allow_residence_proxy':True,'decision_source':str(paths['lead_decision']),
                'decision_sha256':digest(paths['lead_decision'])}
 for item in recent_evidence:
  packet=read(Path(item['observation_path']))
  for name in ('snapshot','age_projection','diagnostic_allow_residence_proxy'):
   require(packet['metadata'][name]==approximation[name],'Recent maintained convention mismatch: '+name)
 rows=[];mapping_checks=[];missing=[]
 numeric_metadata=('empirical_standard_error','synthetic_scale','reference_inverse_variance','reference_inverse_squared_synthetic_scale')
 for original in targets:
  row=copy.deepcopy(original);name=row['restriction_id'];w=weightmap[name];lead=leadrows[name]
  value=numeric(row['target']);require(value==numeric(w['target'])==lead['target'],'Target changed: '+name)
  for field in ('empirical_builder','definition','empirical_record_id','empirical_source_path','empirical_provenance_contract','role'):
   require(row[field]==w[field]==lead[field],'Target provenance/role changed: '+name+'.'+field)
  row['target']=value
  for field in numeric_metadata:row[field]=None if row.get(field,'')=='' else numeric(row[field])
  # Remove numerical baseline readout fields from a target-only contract.
  for field in ('model','gap','model_available','model_source_path','model_source_json_location','loss_contribution'):
   row.pop(field,None)
  row['calibrated_smm']=False;row['actual_weight']=lead['working_weight'];row['working_scale']=lead['working_scale']
  row['weight_status']='separate_unscored_normalization' if name=='initial_normalization' else 'frozen_lead_selected_working_minimum_distance_weight'
  row['weight_rationale']=lead['weight_rationale'];row['weight_decision_source']=str(paths['lead_decision'])
  if name=='initial_normalization':
   require(value==2.1 and w['working_scale']==w['working_weight']=='','Normalization must remain unweighted')
   row['sample']=row['definition'];row['sample_provenance']='Author-selected model normalization, as explicitly recorded in the preserved target definition'
   row['model_observation']='final.normalization.completed_fertility'
  else:
   require(numeric(w['working_scale'])==lead['working_scale'] and numeric(w['working_weight'])==lead['working_weight'],'Weight CSV/decision mismatch: '+name)
   require(row['actual_weight']>0 and math.isclose(row['actual_weight'],1/row['working_scale']**2,rel_tol=1e-14),'Working scale/weight arithmetic mismatch: '+name)
   recordid=row['empirical_record_id']
   if recordid in fertility_records:
    record=fertility_records[recordid];require(value==record['value'] and row['definition']==record['definition'],'Fertility authoritative target mismatch: '+name)
    if name.startswith('cps_'):
     row['sample']=availability['sample']+'; '+cps_windows[0]+'; '+record['definition']
     row['sample_provenance']=str(paths['fertility_availability'])+'#/sample and '+str(paths['fertility_provenance'])+'#/candidate_moments/'+str(fertility_indexes[recordid])+'/definition'
     row['sample_weights']=availability['weights']
     chosen=record['uncertainty']['pooled_candidates']['pooled_fixed_bases_max_positive_correlation_se']
    else:
     row['sample']=fert['nchs_source_cautions']['sample']+'; '+record['definition']
     row['sample_provenance']=str(paths['fertility_provenance'])+'#/nchs_source_cautions/sample and #/candidate_moments/'+str(fertility_indexes[recordid])+'/definition'
     row['residence_filter_caution']=fert['nchs_source_cautions']['resident_filter']
     chosen=record['uncertainty']['annual_2003_2006_sample_sd']
    require(chosen==row['working_scale'],'Selected fertility scale differs from authoritative record')
    row['authoritative_record']=copy.deepcopy(record)
   else:
    require(recordid in housing_records,'Missing authoritative housing/wealth record: '+recordid)
    record=housing_records[recordid];require(value==record['estimate'],'Housing/wealth authoritative target mismatch: '+name)
    row['sample']=record['sample'];row['sample_provenance']=str(paths['housing_wealth_provenance'])+'#/records/'+str(housing_indexes[recordid])+'/sample'
    for field in ('fixed_effects','clustering','geography'):row[field]=record[field]
    row['authoritative_record']=copy.deepcopy(record)
    if recordid.startswith('acs_'):row['common_acs_sample']=copy.deepcopy(housing['common_acs_sample'])
   row['model_observation']=score.MOMENT_PATHS[name].format(projection=decision['primary_cps_projection'])
   if name=='recent_parent_ownership':
    require(value==score.RECENT_TARGET,'Recent ACS target changed')
    row['maintained_model_approximation']=copy.deepcopy(approximation)
    mapping_checks.append({'restriction_id':name,'mapping':row['model_observation'],'baseline_and_joint_packets_verified':True,'legacy_unavailable_row_not_used':True})
   else:
    path=row['model_observation'];lookup(baseline,path);lookup(joint,path)
    mapping_checks.append({'restriction_id':name,'mapping':path,'baseline_and_joint_packets_verified':True})
  for field in score.PROVENANCE_FIELDS:
   if not isinstance(row.get(field),str) or not row[field].strip():missing.append(name+'.'+field)
  rows.append(row)
 require(not missing,'Missing provenance fields: '+', '.join(missing))
 params=csv_rows(paths['parameter_table']);restrictions=[]
 for row in params:
  if row['parameter'] in score.PARAMETERS:
   restrictions.append({k:numeric(row[k]) if k in ('lower','upper') else row[k] for k in ('parameter','lower','upper','transform')})
 require(len(restrictions)==9 and {r['parameter'] for r in restrictions}==set(score.PARAMETERS),'Nine explicit parameter restrictions required')
 require('.01*(high-low)' in paths['inherited_parameter_gate'].read_text(),'Inherited near-bound rule not located')
 artifacts={name:{'path':str(path),'sha256':digest(path)} for name,path in paths.items()}
 contract={'schema':score.SCHEMA,'contract_id':decision['name'],'objective_name':decision['name'],
           'status':'frozen_working_contract_pending_lead_review_and_real_evaluation_receipts',
           'cps_projection':decision['primary_cps_projection'],'cps_sensitivity_rule':decision['cps_sensitivity'],
           'target_rows':rows,'parameter_restrictions':restrictions,'near_bound_fraction':.01,
           'near_bound_rule_source':artifacts['inherited_parameter_gate'],'normalization_tolerance':5e-4,
           'source_fingerprints':source_fingerprints,'recent_parent_approximation':approximation,
           'source_provenance':{'economic_source':{'commit':economic['source_commit'],'file_count':634,'manifest_sha256':economic['source_manifest_sha256'],'source_sha256':economic_pins},
                                'observation_snapshot':{'commit':observer['source_commit'],'file_count':641,'manifest_sha256':observer['source_fingerprint'],'source_sha256':observer_pins},
                                'inherited_economic_files_identical':634,'added_files':additions,
                                'distinction':'The 70abd4a8 snapshot preserves all 634 economic-source files and adds seven files, including passive recent-parent observation and separate history tools. The additions are listed explicitly; they are not a change to the inherited economic code.'},
           'input_artifacts':artifacts,'mapping_verification':mapping_checks,'recent_mapping_evidence':recent_evidence,
           'weight_interpretation':decision['weight_interpretation'],'weight_rationales':{k:decision[k] for k in ('nchs_scale_rationale','cps_scale_rationale','other_scales_rationale')},
           'historical_provenance_status_note':'Original authoritative records retain their historical not-activated/unresolved flags verbatim. The supplied working-weight and maintained-approximation choices come from the separately pinned lead decision; original packets and empirical records are not relabelled.',
           'identification_note':decision['identification'],'calibrated_smm':False,'benchmark_certified':False,
           'evaluation_certification_fabricated':False,'numerical_score_computed':False}
 # Structural compatibility only; never construct a synthetic certificate for a real case.
 require(set(score._rows(rows,'restriction_id','target'))==set(score.MOMENT_PATHS)|{score.NORMALIZATION_ID},'Scorer row schema mismatch')
 score._parameter_table(params,restrictions,.01)
 tests=subprocess.run([sys.executable,'-B','-m','unittest','discover','-s',str(ROOT),'-p','test_score_initial.py','-q'],capture_output=True,text=True)
 require(tests.returncode==0,'Existing scorer tests failed: '+tests.stderr)
 data=json.dumps(contract,indent=2,sort_keys=True,allow_nan=False)+'\n';output=ROOT/'working_contract.json'
 previous_sha=digest(output) if output.exists() else None
 if output.exists() and output.read_text()!=data:
  require(replace_unreviewed,'Frozen working contract exists with different bytes; review before replacing it')
  output.write_text(data)
 elif not output.exists():output.write_text(data)
 receipt={'status':'built_frozen_contract_pending_lead_review','built_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'working_contract_path':str(output),'working_contract_file_sha256':digest(output),
          'predecessor_unreviewed_contract_file_sha256':previous_sha if replace_unreviewed else None,
          'working_contract_canonical_sha256':score.fingerprint(contract),
          'contract_fingerprint_convention':'Pass canonical_sha256 to score_initial; file_sha256 authenticates serialized artifact bytes.',
          'input_artifacts':artifacts,'economic_source_pins':634,'observer_snapshot_pins':641,'inherited_pins_identical':634,'added_files':list(additions),
          'restriction_count':13,'scored_restriction_count':12,'structural_parameter_count':9,'primary_cps_projection':'uniform_birth_time',
          'normalization_unscored':True,'model_paths_verified':mapping_checks,'recent_observation_evidence':recent_evidence,
          'missing_required_fields':missing,'scorer_tests':{'returncode':tests.returncode,'stdout':tests.stdout,'stderr':tests.stderr},
          'no_real_score_or_evaluation_receipt_created':True,'raw_observations_modified':False,'model_solves':0,'cluster_jobs_submitted':0,
          'remaining':'Lead reviews complete contract and creates independently verified evaluation receipts before numerical scoring.'}
 (ROOT/'build_receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True,allow_nan=False)+'\n')
 print(json.dumps({k:receipt[k] for k in ('status','working_contract_file_sha256','working_contract_canonical_sha256','restriction_count','scored_restriction_count','missing_required_fields')}))

if __name__=='__main__':
 try:
  require(sys.argv[1:] in ([],['--replace-unreviewed']),'Only --replace-unreviewed is supported')
  build(replace_unreviewed=sys.argv[1:]==['--replace-unreviewed'])
 except Exception as exc:
  (ROOT/'build_receipt.json').write_text(json.dumps({'status':'blocked','error_type':type(exc).__name__,'error':str(exc),'model_solves':0,'cluster_jobs_submitted':0},indent=2)+'\n')
  raise
