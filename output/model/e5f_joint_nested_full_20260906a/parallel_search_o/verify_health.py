"""Recheck completed verification evidence and reduce active-case health."""
from pathlib import Path
import sys,json,time
root=Path.cwd();sys.path[:0]=[str(root/'code/model'),str(root/'code/model/tools')]
import run_e5f_joint_overnight_case as a
import run_e5f_joint_nested_long_search as s
p=a.read_json(root/'preflight_verification.json');assert p['status']=='pass_complete_preflight'
for rel,sha in a.read_json(root/'source_manifest.json').items():a.verify(root/rel,sha)
for path,sha in p['policy_artifact_sha256'].items():a.verify(path,sha)
c=s.verify_contract(root/'output/model/joint_nested_overnight/contract.json',p['contract_sha256'])
proof=a.read_json(Path(c['imported_smoke']['root'])/'smoke_verification.json')
anchor_hashes={a.digest(Path(path).parent/'summary.json') for path in proof['anchor_receipts']}
pol=a.read_json(c['parallel_policy_timing']['receipt'])
s.validate_policy_receipt(pol,c,smoke=True,selected_hashes=anchor_hashes)
folder=root/'output/model/joint_nested_overnight/search/initial_population';plan=folder/'plan.json'
a.load_plan(plan,a.digest(plan))
health=[]
for path in sorted(folder.glob('task_*/heartbeat.json')):
 r=a.read_json(path);health.append(dict(case=path.parent.name,heartbeat_age_seconds=time.time()-r['epoch'],elapsed_seconds=r['elapsed_seconds']))
assert len(health)==32 and max(r['heartbeat_age_seconds'] for r in health)<300
result=dict(status='healthy_initial_population',epoch=time.time(),job=17106283,contract_sha256=p['contract_sha256'],verified_source_files=70,verified_policy_artifacts=len(p['policy_artifact_sha256']),imported_histories=4,new_completed_histories=len(list(folder.glob('task_*/case_receipt.json'))),active_case_heartbeats=health,policy_smoke_elapsed_seconds=pol['elapsed_seconds'],production_promoted=False)
a.write_json(root/'health_verification.json',result)
print(json.dumps({k:v for k,v in result.items() if k!='active_case_heartbeats'}))
print('Maximum heartbeat age:',max(r['heartbeat_age_seconds'] for r in health))
