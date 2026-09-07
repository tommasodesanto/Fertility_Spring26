import sys,json,hashlib,time,re,argparse
from pathlib import Path
root=Path(__file__).resolve().parents[3]
parser=argparse.ArgumentParser(description="Build an immutable experimental run contract after scientific source review")
parser.add_argument('--remote-root',required=True)
parser.add_argument('--seed-summary',type=Path,required=True)
parser.add_argument('--outdir',type=Path,required=True)
parser.add_argument('--finish-epoch',type=float,required=True)
parser.add_argument('--expected-history-seconds',type=float,required=True)
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
 code_bundle_sha256=bundle,search_domain=adapter.SEARCH_DOMAIN,max_workers=12,case_timeout_seconds=3600,max_histories=360,
 max_search_seconds=32400,max_total_seconds=43200,population_size=32,max_generations=8,polish_rounds=2,smoke_histories=4,
 random_seed=20260906,absolute_finish_epoch=args.finish_epoch,
 expected_history_seconds=args.expected_history_seconds,expected_history_solve_count_upper=160,
 runtime_estimate_status="provisional until complete exhaustive-saving smoke is measured",
 run_size='120 wealth x6 housing x1 market x17 ages x15 income x4 parity x4 child counts; up to360 complete histories (including4smoke), each five cleared historical dates and normalized old steady state; 22 final Jacobian probes and2 exact repeats are within360',
 estimated_search_wall_hours=360*args.expected_history_seconds/12/3600,policy_path_dates=11,policy_path_cases=4,production_promoted=False,
 closure={'expectations':'current-date prices treated as permanent; temporary equilibrium, not perfect foresight',
          'post2023_population':'maintained closed national: M=0,rho=1; inherited four-slot birth queue',
          'historical_population':'unchanged Census totals and ACS householder-age bridge 2007-2023',
          'fixed_supply_elasticity':.63,'replacement_conversion':1/2.1,'fiscal':'tax revenue discarded; no grants or rebates',
          'outstanding':['author adoption of nesting and common lambda','perfect-foresight expectation extension','final production promotion']},
 saving_maximization='exhaustive_piecewise_linear_continuation',
 numerical_gates={'market':2e-4,'mass':2e-10,'population':2e-10,'stationary_measurement':2e-8,'childless_identity':2e-10,
                  'occupied_value_negative_steps':0,'realized_budget_excess_mass':2e-10,'budget_gap_threshold':1e-9})
adapter.write_json(out/'contract.json',c)
print(json.dumps(dict(bundle=bundle,contract_sha256=adapter.digest(out/'contract.json'),remote=remote,finish_epoch=c['absolute_finish_epoch'])))
