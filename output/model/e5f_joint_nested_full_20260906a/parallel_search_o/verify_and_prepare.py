"""Fail-closed source, original-state and full-smoke checks before search."""
from pathlib import Path
import argparse,json,sys,subprocess,os,time
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/'code/model'),str(root/'code/model/tools')]
import run_e5f_joint_overnight_case as a
import run_e5f_joint_nested_long_search as search
import run_e5f_transition_calibration as cal
from intergen_eqscale_seq_optimized import solver
ap=argparse.ArgumentParser();ap.add_argument('--source-only',action='store_true');ap.add_argument('--run-search-after-verification',action='store_true');args=ap.parse_args()
for rel,sha in a.read_json(root/'source_manifest.json').items():a.verify(root/rel,sha)
assert cal.code_fingerprint_contract(solver)['bundle_sha256']==a.BUNDLE
subprocess.run([sys.executable,str(root/'code/model/tools/test_e5f_joint_nested_long_search.py')],check=True)
if args.source_only:
 print(json.dumps(dict(status='pass_source_only',bundle=a.BUNDLE)));sys.exit(0)
j=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907m/output/model/joint_nested_overnight')
a.verify(j/'contract.json','cb22a9d7ee19da10e558eada48cd5fad5bb1055eb8dc7e01220699481c426c15')
proof_path=j/'support_repair/original_state_verification.json'
proof=a.read_json(proof_path)
assert proof['status']=='pass_original_failed_policy_state' and proof['new_bundle']==a.BUNDLE
assert proof['budget_before']['budget_excess_mass']==3.727996554261621e-10 and proof['budget_after']['budget_excess_mass']==0
assert proof['policy_arrays']['occupied_negative_steps']==0 and proof['projection_mass']<=1e-6
checkpoint=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907j/output/model/joint_nested_overnight/policy_loop_smoke/supply-plus-20/date_2023/dated_state.pkl.gz')
assert proof['checkpoint_sha256']=='ffdd4ee394db635f0278a99306177e01eee9ee6d42312670c22c7e0a4b4e6639'
a.verify(checkpoint,proof['checkpoint_sha256'])
default=a.read_json(j/'default_off_reference/baseline_reference.json');assert default['status']=='exact' and len(default['arrays'])==10
smoke=j/'smoke';history=a.read_json(smoke/'smoke_verification.json');policy=a.read_json(smoke/'policy_loop_verification.json')
assert history['status']==policy['status']=='pass'
policy_receipt_path=Path(policy['receipt']);a.verify(policy_receipt_path,policy['sha256']);pol=a.read_json(policy_receipt_path)
# Check every dated gate and archive hashes of the newly generated artifacts.
policy_hashes={}
for name,branch in pol['cases'].items():
 folder=policy_receipt_path.parent/name
 assert a.read_json(folder/'receipt.json')==branch
 assert [int(float(r['calendar_year'])) for r in a.read_csv(folder/'policy_path.csv')]==[2023,2027]
 for year in (2023,2027):
  date=folder/f'date_{year}';budget=a.read_json(date/'budget_summary.json');arrays=a.read_json(date/'policy_array_summary.json')
  assert budget['budget_excess_mass']<=2e-10 and arrays['occupied_negative_steps']==0
  for p in arrays['probabilities'].values():assert p['nonfinite']==0 and 0<=p['minimum']<=p['maximum']<=1
  graphs=list((date/'standard_diagnostics').glob('*.png'));assert len(graphs)==17
  for f in (date/'dated_state.pkl.gz',date/'budget_summary.json',date/'policy_array_summary.json',*graphs):policy_hashes[str(f)]=a.digest(f)
 for f in (folder/'receipt.json',folder/'policy_path.csv'):policy_hashes[str(f)]=a.digest(f)
for f in (policy_receipt_path,policy_receipt_path.parent/'inherited_state_verification.json'):policy_hashes[str(f)]=a.digest(f)
assert len(policy_hashes)==170
for path in history['anchor_receipts']+history['probe_receipts']:
 r=a.read_json(path)
 assert r['budget_diagnostic']['budget_excess_mass']<=2e-10 and r['policy_array_diagnostic']['occupied_negative_steps']==0
# The actual require_smoke method below additionally checks all historical
# receipts, source/target/domain/closure/gates, exact-repetition proof, policies.
seed=Path(a.read_json(smoke/'best_so_far.json')['best']['summary'])
times=[]
for path in history['anchor_receipts']+history['probe_receipts']:
 r=a.read_json(path);times.append(float(r['elapsed_seconds']))
subprocess.run([sys.executable,str(root/'code/model/tools/build_e5f_joint_nested_long_contract.py'),
 '--remote-root',str(root),'--seed-summary',str(seed),'--outdir',str(root/'output/model/joint_nested_overnight'),
 '--finish-epoch','1788788100','--expected-history-seconds',str(max(times)),
 '--runtime-estimate-status','measured','--profile','parallel32_overlap','--policy-workers','4',
 '--parallel-policy-receipt',str(policy_receipt_path),'--imported-smoke-root',str(smoke),
 '--imported-smoke-contract',str(j/'contract.json'),'--imported-smoke-contract-sha256',a.digest(j/'contract.json'),
 '--imported-smoke-verification-sha256',a.digest(smoke/'smoke_verification.json'),
 '--imported-policy-verification-sha256',a.digest(smoke/'policy_loop_verification.json')],check=True)
contract_path=root/'output/model/joint_nested_overnight/contract.json';contract_sha=a.digest(contract_path)
contract=search.verify_contract(contract_path,contract_sha)
# Use an isolated preflight output directory; never initialize the actual search twice.
preflight=dict(contract,output_root=str(root/'output/model/preflight'))
run=search.Search(preflight,'search');run.require_smoke()
result=dict(status='pass_complete_preflight',epoch=time.time(),scientific_bundle=a.BUNDLE,source_files=len(a.read_json(root/'source_manifest.json')),
 original_state_proof_sha256=a.digest(proof_path),original_checkpoint_sha256=proof['checkpoint_sha256'],
 historical_smoke_sha256=a.digest(smoke/'smoke_verification.json'),policy_smoke_sha256=a.digest(policy_receipt_path),
 policy_artifact_sha256=policy_hashes,measured_history_seconds=times,projected_full_policy_seconds=pol['elapsed_seconds']*44/8,contract_sha256=contract_sha,
 selection_rule=contract['final_selection_rule'],hard_cutoff_epoch=1788788100,production_promoted=False)
a.write_json(root/'preflight_verification.json',result);print(json.dumps(result),flush=True)
if args.run_search_after_verification:
 os.environ.update(E5F_JOINT_MODE='search',E5F_JOINT_CONTRACT=str(contract_path),E5F_JOINT_CONTRACT_SHA256=contract_sha)
 os.execv('/bin/bash',['bash',str(root/'code/cluster/submit_e5f_joint_nested_long.sh')])
