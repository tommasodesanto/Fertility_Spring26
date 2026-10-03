"""Generate five matched bounded starts; no model runtime or solve."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
PACKETS=ROOT/'output/model/fixed_reference_economics_20260928'
ANCHOR=PACKETS/'soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/completed.json'
OUT=ROOT/'output/model/experiments/birth_count_choice/estate_a_calibration_v1'
PAUSED=PACKETS/'alternative_wealth_cluster_20261003_v1/pause_20261003.json'
def canonical(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def main():
    anchor=json.loads(ANCHOR.read_text());base=dict(anchor['selected']['parameters'])
    assert anchor['status']=='selected_numerically_verified' and anchor['arm']=='alternative' and anchor['chain']==13
    bounds={r['parameter']:[float(r['lower']),float(r['upper'])] for r in anchor['parameters'] if r['parameter'] in base}
    assert bounds['beta_annual']==[.94,.99];bounds['beta_annual']=[.93,.99]
    paused=json.loads(PAUSED.read_text());chain=next(c for c in paused['chains'] if c['chain']==6);checkpoint=chain['files']['best_so_far']['best']
    assert paused['status']=='author_stopped_checkpoints_preserved' and paused['job_array']=='19111687'
    assert checkpoint['label']=='0060_nm' and checkpoint['loss']==22.141841386410267 and checkpoint['exploration_unverified'] is True
    starts=[base,dict(checkpoint['parameters'])];patterns=(
      {'beta_annual':-.0015,'kappa_fert':.04,'tenure_choice_kappa':-.035,'child_benefit_curvature':.04},
      {'theta0':.04,'kappa_fert_continuation':-.035,'h_P':.012,'psi_child':-.025},
      {'beta_annual':.001,'chi':-.025,'first_birth_fixed_cost':.035,'kappa_fert':-.035,'theta0':-.04})
    for pattern in patterns:
        row=dict(base)
        for k,shift in pattern.items():row[k]=min(bounds[k][1],max(bounds[k][0],base[k]*(1+shift)))
        starts.append(row)
    target=[{k:r[k] for k in ('moment','target','weight','role')} for r in anchor['target_fit']]
    wealth=next(r for r in target if r['moment']=='wealth_earnings')
    assert wealth['target']=='6.92658379107299' and wealth['weight']=='7.595098472533724'
    wealth['target']='4.45838713455674'
    assert [{k:r[k] for k in ('moment','target','weight','role')} for r in checkpoint['target_fit']]==target
    assert checkpoint['weight_fingerprint']==canonical(dict(base_contract=target,multipliers={}))
    assert set(checkpoint['parameters'])==set(bounds) and all(lo<=checkpoint['parameters'][k]<=hi for k,(lo,hi) in bounds.items())
    plan=dict(arms={'binary':1,'count3':3},starts=starts,bounds=bounds,target_contract=target,
       target_fingerprint=canonical(target),weight_fingerprint=canonical(dict(base_contract=target,multipliers={})),
       source_checkpoint=str(ANCHOR.relative_to(ROOT)),source_checkpoint_sha256=hashlib.sha256(ANCHOR.read_bytes()).hexdigest(),
       provisional_seed_source=str(PAUSED.relative_to(ROOT)),provisional_seed_source_sha256=hashlib.sha256(PAUSED.read_bytes()).hexdigest(),
       provisional_seed=dict(start_index=1,array='19111687',chain=6,label='0060_nm',saved_loss=22.141841386410267,verified=False,seed_only=True,not_adopted=True),
       deterministic_seed=20261003,perturbation_patterns=patterns,start_generation='chain13, preserved unverified new-target chain6 checkpoint, and three deterministic bounded chain13 relative perturbations; no RNG',
       coordinates=list(base),free_parameter_count=10,positive_weight_target_count=10,identification_rank_certified=False,
       per_chain=dict(cpus=1,memory_GiB=24,wall_seconds=21600,max_objective_calls=500,max_lifecycle_per_GE=32,final_native_reserve_seconds=1800),
       task_map=[dict(task=i,arm=arm,birth_cap=cap,chain=j) for i,(arm,cap,j) in enumerate((a,c,j) for a,c in [('binary',1),('count3',3)] for j in range(5))],
       solve_budget=dict(total_chains=10,search_max_GE_calls=5000,search_max_lifecycle_solves=160000,final_fresh_GE_calls=10,
          final_max_lifecycle_solves=320,smoke_search_GE_calls=4,smoke_final_fresh_GE_calls=2,smoke_max_lifecycle_solves=192,
          typical_fixed_price_seconds=6,typical_GE_seconds=90,rough_max_search_hours_at_90sec=125,
          rough_calls_per_chain_in_19800seconds_at_90sec=220,
          rough_parallel_wall_hours_for_500calls_at_90sec=12.5,actual_chain_wall_cap_hours=6,
          uncertainty='Prior unchanged economics timing is indicative only; estate and birth menu changes may alter convergence and solve times. Time cap will likely bind before 500calls.'),
       economic_changes=dict(estate_A='utility and deceased-account estate bp+(1-psi)*P*h_prime; no extra R',
          birth_menu='matched engine cap1 current binary economics versus cap3',wealth_target='experimental narrowerPSID wealth target with retained weight',
          beta_bound='both arms .93 to .99',unchanged='all other primitives, entry, timing, targets, weights and numerical gates'),
       scf_provisional='Retained target .007291... is not reconciled to recipient and wealth-scope conventions; no row dropped or reweighted.',
       storage_budget=dict(observed_previous_native_GE_MiB=82,per_objective_retained_cap_GiB=1,all_chains_search_max_GiB=5000,submission_required_free_GiB=5200,runtime_shared_free_floor_GiB=350,method='Preserve generated reports, arrays and plots; no pruning; stop on cap.'),
       no_auto_retry=True,no_parameter_promotion=True,production_release_required=True)
    assert len({canonical(s) for s in starts})==5
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'start_plan.json').write_text(json.dumps(plan,indent=2,sort_keys=True)+'\n')
    print(json.dumps(dict(status='prepared_zero_solves',starts=5,chains=10,target_fingerprint=plan['target_fingerprint'],weight_fingerprint=plan['weight_fingerprint'])))
if __name__=='__main__':main()
