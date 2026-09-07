import sys,json,hashlib,time,re,argparse,math
from pathlib import Path
root=Path(__file__).resolve().parents[3]
parser=argparse.ArgumentParser(description="Build an immutable experimental run contract after scientific source review")
parser.add_argument('--remote-root',required=True)
parser.add_argument('--seed-summary',type=Path,required=True)
parser.add_argument('--outdir',type=Path,required=True)
parser.add_argument('--finish-epoch',type=float,required=True)
parser.add_argument('--expected-history-seconds',type=float,required=True)
parser.add_argument('--runtime-estimate-status',choices=('measured','provisional'),required=True)
parser.add_argument('--profile',choices=('v1','wide32','parallel32'),required=True)
parser.add_argument('--policy-workers',type=int,choices=(1,4),default=1)
parser.add_argument('--parallel-policy-receipt',type=Path)
parser.add_argument('--imported-smoke-root',type=Path)
parser.add_argument('--imported-smoke-contract',type=Path)
parser.add_argument('--imported-smoke-contract-sha256')
parser.add_argument('--imported-smoke-verification-sha256')
parser.add_argument('--imported-policy-verification-sha256')
args=parser.parse_args()
if not 0 < args.expected_history_seconds <= 3600:raise ValueError('Expected history time must fit the per-case cap')
sys.path[:0]=[str(root/'code/model'),str(root/'code/model/tools')]
import run_e5f_transition_calibration as cal
from intergen_eqscale_seq_optimized import solver
bundle=cal.code_fingerprint_contract(solver)['bundle_sha256']

import run_e5f_joint_overnight_case as adapter
remote=args.remote_root.rstrip('/')
if bundle != adapter.BUNDLE:raise RuntimeError('Scientific bundle changed; review before pinning the case adapter')
out=args.outdir.resolve();out.mkdir(parents=True,exist_ok=True)
if (out/'contract.json').exists():raise RuntimeError('Refusing to replace an existing run contract')
profile={'v1':dict(max_workers=12,max_histories=360,population_size=32,max_generations=8,polish_rounds=2),
         'wide32':dict(max_workers=32,max_histories=640,population_size=64,max_generations=8,polish_rounds=2,final_reserve_seconds=16200),
         'parallel32':dict(max_workers=32,max_histories=640,population_size=32,max_generations=8,polish_rounds=2,final_reserve_seconds=12600)}[args.profile]
available_total_seconds=min(43200,max(0,args.finish_epoch-time.time()))
final_reserve=profile.get('final_reserve_seconds',10800)
available_search_seconds=min(32400,max(0,available_total_seconds-final_reserve))
if available_search_seconds < math.ceil(profile['population_size']/profile['max_workers'])*3600:
 raise RuntimeError('Deadline cannot fit the initial population and reserved final verification')
projected_search_histories=min(profile['max_histories']-28,
 int(available_search_seconds/args.expected_history_seconds)*profile['max_workers'])
imported_fields=(args.imported_smoke_root,args.imported_smoke_contract,args.imported_smoke_contract_sha256,
                 args.imported_smoke_verification_sha256,args.imported_policy_verification_sha256)
if (args.profile in ('wide32','parallel32') or any(imported_fields)) and not all(imported_fields):
 raise ValueError('Parallel search profiles require complete imported-smoke provenance')
original=None; imported_smoke=None
if all(imported_fields):
 smoke_root=args.imported_smoke_root.resolve(); original_contract=args.imported_smoke_contract.resolve()
 adapter.verify(original_contract,args.imported_smoke_contract_sha256)
 original=adapter.read_json(original_contract)
 # Check proof bytes before producing any launchable contract.
 adapter.verify(smoke_root/'smoke_verification.json',args.imported_smoke_verification_sha256)
 adapter.verify(smoke_root/'policy_loop_verification.json',args.imported_policy_verification_sha256)
 imported_smoke={'root':str(smoke_root),'original_contract':str(original_contract),'original_contract_sha256':args.imported_smoke_contract_sha256,
                 'smoke_verification_sha256':args.imported_smoke_verification_sha256,'policy_loop_verification_sha256':args.imported_policy_verification_sha256}
raw=args.seed_summary.resolve();s=adapter.read_json(raw)
seed={'status':'uncertified starting parameters; require fresh full-loop verification','old_psi_child':s['old_psi_child'],
      'best_candidate':{k:s['best_candidate'][k] for k in ('theta','new_psi_child')},'panel_design':{'unit_vector':s['panel_design']['unit_vector']}}
seed['best_candidate']['old_psi_child']=s['old_psi_child']
adapter.write_json(out/'seed_center.json',seed);(out/'seed_reference.json').write_bytes(raw.read_bytes())
helpers=['run_e5f_independent_numerical_audit.py','run_e5f_post2023_no_policy_continuations.py','run_e5f_post2023_policy_mechanisms.py']
base=dict(schema='e5f_joint_nested_overnight_v1',stage='template',source=s['source'],source_sha256=adapter.SOURCE,
 target_set=s['target_set'],target_fingerprint=adapter.TARGET,code_bundle_sha256=bundle,search_domain=adapter.SEARCH_DOMAIN,
 first_child_jump_upper=2.,adapter_sha256=adapter.digest(adapter.__file__),helper_sha256={h:adapter.digest(root/'code/model/tools'/h) for h in helpers},
 input_sha256={},production_promoted=False,cases=[],launch_deadline_epoch=args.finish_epoch)
c=dict(schema='e5f_joint_nested_long_v1',base_plan=base,base_plan_sha256=hashlib.sha256((json.dumps(base,indent=2,sort_keys=True)+'\n').encode()).hexdigest(),
 controller_sha256=adapter.digest(root/'code/model/tools/run_e5f_joint_nested_long_search.py'),adapter_sha256=adapter.digest(adapter.__file__),
 planner_sha256=adapter.digest(root/'code/model/tools/build_e5f_bounded_refinement_plan.py'),
 finalizer_driver=remote+'/code/model/tools/run_e5f_joint_nested_finalize.py',finalizer_sha256=adapter.digest(root/'code/model/tools/run_e5f_joint_nested_finalize.py'),
 seed_center=remote+'/output/model/joint_nested_overnight/seed_center.json',seed_center_sha256=adapter.digest(out/'seed_center.json'),
 seed_reference=remote+'/output/model/joint_nested_overnight/seed_reference.json',seed_reference_sha256=adapter.digest(out/'seed_reference.json'),
 output_root=remote+'/output/model/joint_nested_overnight',source_sha256=adapter.SOURCE,target_fingerprint=adapter.TARGET,
 code_bundle_sha256=bundle,search_domain=adapter.SEARCH_DOMAIN,run_profile=args.profile,case_timeout_seconds=3600,
 max_search_seconds=32400,max_total_seconds=43200,smoke_histories=4,**profile,
 random_seed=20260906,absolute_finish_epoch=args.finish_epoch,
 expected_history_seconds=args.expected_history_seconds,expected_history_solve_count_upper=160,
 runtime_estimate_status=args.runtime_estimate_status,
 run_size=f'120 wealth x6 housing x1 market x17 ages x15 income x4 parity x4 child counts; up to{profile["max_histories"]} attempted complete histories (including4 imported smoke histories), each five cleared historical dates and normalized old steady state; 22 final Jacobian probes and2 exact repeats are reserved',
 estimated_history_wall_hours=profile['max_histories']*args.expected_history_seconds/profile['max_workers']/3600,policy_path_dates=11,policy_path_cases=4,policy_workers=args.policy_workers,production_promoted=False,
 budget_estimate={'at_contract_creation_epoch':time.time(),'available_total_seconds':available_total_seconds,
  'available_search_seconds':available_search_seconds,'reserved_final_seconds':final_reserve,
  'projected_search_histories_at_measured_rate':projected_search_histories,
  'interpretation':'640/360 is a hard attempt ceiling, not a planned completion count; actual stages require their full timeout waves before the fixed cutoff'},
 imported_smoke=imported_smoke,
 closure={'expectations':'current-date prices treated as permanent; temporary equilibrium, not perfect foresight',
          'post2023_population':'maintained closed national: M=0,rho=1; inherited four-slot birth queue',
          'historical_population':'unchanged Census totals and ACS householder-age bridge 2007-2023',
          'fixed_supply_elasticity':.63,'replacement_conversion':1/2.1,'fiscal':'tax revenue discarded; no grants or rebates',
          'outstanding':['author adoption of nesting and common lambda','perfect-foresight expectation extension','final production promotion']},
 saving_maximization='exhaustive_piecewise_linear_continuation',
 numerical_gates={'market':2e-4,'mass':2e-10,'population':2e-10,'stationary_measurement':2e-8,'childless_identity':2e-10,
                  'occupied_value_negative_steps':0,'realized_budget_excess_mass':2e-10,'budget_gap_threshold':1e-9})
for key in ('source_sha256','code_bundle_sha256','target_fingerprint','search_domain'):
 if original is not None and original.get(key)!=c[key]:raise RuntimeError(f'Imported smoke original contract has mixed {key}')
if args.profile == 'parallel32':
 if args.parallel_policy_receipt is None:raise ValueError('A measured complete parallel policy smoke receipt is required')
 timing_path=args.parallel_policy_receipt.resolve(); timing=adapter.read_json(timing_path)
 c['parallel_policy_timing']={'receipt':str(timing_path),'receipt_sha256':adapter.digest(timing_path),
  'projected_full_policy_seconds':timing['elapsed_seconds']*44/8,
  'reserved_history_waves_seconds':7200,'buffer_seconds':1200,
  'interpretation':'Linear projection of measured four-process eight-date smoke to forty-four dates; a runtime forecast, not a completion guarantee'}
 import run_e5f_joint_nested_long_search as search
 search.verify_parallel_policy_budget(c)
elif args.parallel_policy_receipt is not None:
 raise ValueError('Parallel timing evidence is only accepted for the parallel32 profile')
adapter.write_json(out/'contract.json',c)
print(json.dumps(dict(bundle=bundle,contract_sha256=adapter.digest(out/'contract.json'),remote=remote,finish_epoch=c['absolute_finish_epoch'])))
