"""Read-only full receipt/anchor check; writes one verification receipt."""
from pathlib import Path
import sys,json
root=Path.cwd();sys.path[:0]=[str(root/'code/model'),str(root/'code/model/tools')]
import run_e5f_joint_overnight_case as a
smoke=root/'output/model/joint_nested_overnight/smoke';proof=a.read_json(smoke/'smoke_verification.json')
assert proof['status']=='pass'
rows=[];hashes=0
for path in proof['anchor_receipts']+proof['probe_receipts']:
 receipt=a.read_json(path);folder=Path(path).parent
 planpath=folder.parent/'plan.json';plan=a.load_plan(planpath,receipt['plan_sha256'])
 case=next(c for c in plan['cases'] if c['id']==receipt['case_id'])
 for rel,sha in receipt['artifact_sha256'].items():a.verify(folder/rel,sha);hashes+=1
 a.validate_result(folder,plan,case)
 assert receipt['budget_diagnostic']['budget_excess_mass']<=2e-10
 assert receipt['policy_array_diagnostic']['occupied_negative_steps']==0
 rows.append(dict(case=receipt['case_id'],loss=receipt['loss'],elapsed_seconds=receipt['elapsed_seconds'],artifact_count=len(receipt['artifact_sha256']),budget_excess_mass=receipt['budget_diagnostic']['budget_excess_mass']))
anchors=[Path(p).parent for p in proof['anchor_receipts']]
exact=a.compare_reference(anchors[0],anchors[1]/'summary.json')
for p in (anchors[0]/'standard_diagnostics').glob('*.png'):assert a.digest(p)==a.digest(anchors[1]/'standard_diagnostics'/p.name)
assert hashes==84
result=dict(status='pass_four_completed_histories',scientific_bundle=a.BUNDLE,histories=rows,verified_artifacts=hashes,exact=exact,exact_anchor_pngs=17,smoke_verification_sha256=a.digest(smoke/'smoke_verification.json'),production_promoted=False)
a.write_json(root/'output/model/joint_nested_overnight/completed_history_verification.json',result);print(json.dumps(result))
