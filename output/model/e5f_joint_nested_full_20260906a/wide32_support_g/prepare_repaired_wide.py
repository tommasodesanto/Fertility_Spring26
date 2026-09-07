"""Verify original completed smoke components and prepare the wide contract; no submission."""
from pathlib import Path
import json, subprocess, sys

root=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907g')
old=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907e')
sys.path[:0]=[str(root/'code/model'),str(root/'code/model/tools')]
import run_e5f_joint_nested_long_search as search
a=search.adapter
for name,sha in a.read_json(root/'prepared_wide_source_manifest.json').items():a.verify(root/name,sha)
original_path=old/'output/model/joint_nested_overnight/contract.json'
original_sha='9f512dfcaf7dc85b359e0b65d4f8f1543221f4c59f88d82888b91a1ebf4d3e41'
a.verify(original_path,original_sha); original=a.read_json(original_path)

# The low-scale canary must be adjudicated before enabling the large run.
canary=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907f/output/model/joint_nested_overnight/corner_canary')
canary_assessment=a.read_json(root/'canary_assessment.json')
assert canary_assessment['job_id']==17091265 and canary_assessment['ready_for_wide_search'] is True
assert canary_assessment['source_bundle']==original['code_bundle_sha256']
a.verify(canary/'plan.json',canary_assessment['plan_sha256'])
for rel,sha in canary_assessment['evidence_sha256'].items():a.verify(canary/rel,sha)
assert canary_assessment['evidence_sha256']

hist=Path(original['output_root'])/'smoke'
policy=Path(original['output_root'])/'policy_loop_smoke'
proof=a.read_json(hist/'smoke_verification.json')
assert proof['status']=='pass' and len(proof['anchor_receipts'])==len(proof['probe_receipts'])==2
elapsed=[]; case_hashes={}
for value in proof['anchor_receipts']+proof['probe_receipts']:
    path=Path(value); receipt=a.read_json(path); out=path.parent
    plan=a.load_plan(out.parent/'plan.json',receipt['plan_sha256'])
    case=next(c for c in plan['cases'] if c['id']==receipt['case_id'])
    for rel,sha in receipt['artifact_sha256'].items():a.verify(out/rel,sha)
    a.verify(out.parent/case['center'],case['center_sha256'])
    summary,_,_=a.validate_result(out,plan,case)
    assert receipt['status']=='complete' and receipt['loss']==summary['best_candidate']['transition_loss']
    assert receipt['policy_array_diagnostic']['occupied_negative_steps']==0
    assert receipt['budget_diagnostic']['budget_excess_mass']<=2e-10
    for b in receipt['policy_array_diagnostic']['probabilities'].values():
        assert b['nonfinite']==0 and 0<=b['minimum']<=b['maximum']<=1
    case_hashes[value]=a.digest(path);elapsed.append(receipt['elapsed_seconds'])
anchors=[Path(p).parent for p in proof['anchor_receipts']]
a.compare_reference(anchors[1],anchors[0]/'summary.json')
graphs=list((anchors[0]/'standard_diagnostics').glob('*.png'));assert len(graphs)==17
for graph in graphs:a.verify(anchors[1]/'standard_diagnostics'/graph.name,a.digest(graph))
base=a.read_json(anchors[0]/'summary.json')['panel_design']['unit_vector']
for path in proof['probe_receipts']:
    u=a.read_json(Path(path).parent/'summary.json')['panel_design']['unit_vector']
    assert len(u)==11 and sum(x!=y for x,y in zip(u,base))==11
pr=a.read_json(policy/'equilibrium_receipt.json')
search.validate_policy_receipt(pr,original,smoke=True,selected_hashes={a.digest(p/'summary.json') for p in anchors})
a.verify(old/pr['selected_summary'],pr['selected_summary_sha256'])
a.verify(policy/'inherited_state_verification.json',pr['inherited_state_verification_sha256'])
a.verify(root/'code/model/tools/run_e5f_joint_nested_finalize.py',original['finalizer_sha256'])
policy_artifacts={}
for name,case in pr['cases'].items():
    folder=policy/name;assert a.read_json(folder/'receipt.json')==case
    rows=a.read_csv(folder/'policy_path.csv');assert [int(float(r['calendar_year'])) for r in rows]==[2023,2027]
    for year in (2023,2027):
        date=folder/f'date_{year}';budget=a.read_json(date/'budget_summary.json');arrays=a.read_json(date/'policy_array_summary.json')
        assert budget['budget_excess_mass']<=2e-10 and arrays['occupied_negative_steps']==0
        for b in arrays['probabilities'].values():assert b['nonfinite']==0 and 0<=b['minimum']<=b['maximum']<=1
        pngs=list((date/'standard_diagnostics').glob('*.png'));assert len(pngs)==17
        for f in [date/'dated_state.pkl.gz',date/'budget_summary.json',date/'policy_array_summary.json',*pngs]:
            policy_artifacts[str(f)]=a.digest(f)
    for f in (folder/'receipt.json',folder/'policy_path.csv'):policy_artifacts[str(f)]=a.digest(f)
# A repaired intermediate statistic must lead either to a full valid case
# or to the original controller's declared historical-support rejection.
# This specific replay has the latter outcome; preserve all measurement gates.
replay=old/'output/model/joint_nested_overnight/support_replay'
replay_provenance=a.read_json(replay/'replay_provenance.json')
assert not replay_provenance['complete'] and replay_provenance['unavailable_trial_count']==2
failure=a.read_json(replay/'task_026/adapter_failure.json')
assert failure['error']=='The dated first-birth branch has zero treated mass'
assert search.classify_failure(failure['error'],failure['type'])=='undefined_first_birth_support'
assessment=a.read_json(root/'replay_assessment.json')
assert assessment['status']=='expected_controlled_history_support_rejection_after_repair'
a.verify(replay/'replay_provenance.json',assessment['replay_provenance_sha256'])
a.verify(replay/'task_026/adapter_failure.json',assessment['adapter_failure_sha256'])
assert replay_provenance['new_bundle']==original['code_bundle_sha256']
a.verify(old/'code/model/tools/run_e5f_joint_nested_long_search.py',assessment['controller_sha256'])
assert search.classify_failure(failure['error'],failure['type'])=='undefined_first_birth_support'
a.verify(old/'logs/support_replay_17090361.out',assessment['stdout_sha256'])
rp=a.load_plan(replay/'plan.json',replay_provenance['new_plan_sha256'])
assert len(rp['cases'])==1 and rp['cases'][0]['id']==26
a.verify(replay/rp['cases'][0]['center'],replay_provenance['original_center_sha256'])
latest=replay/'task_026/cases/task_001/latest_completed_period.json'
a.verify(latest,assessment['latest_completed_period_sha256'])
assert a.read_json(latest)['completed_period']==1
# The supported anchor must remain exact across the reporting/status repair.
previous=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907c/output/model/joint_nested_overnight/smoke/smoke_anchor/task_001')
a.verify(previous/'summary.json','066d836f5551861a1e2e563c15a910db1ccdb87e754fec14f4c276b7bea16fcc')
source_comparison=a.compare_reference(anchors[0],previous/'summary.json')
for graph in graphs:a.verify(previous/'standard_diagnostics'/graph.name,a.digest(graph))
out=root/'output/model/joint_nested_overnight'; imported=out/'imported_smoke';imported.mkdir(exist_ok=False)
(imported/'smoke_verification.json').write_bytes((hist/'smoke_verification.json').read_bytes())
a.write_json(imported/'policy_loop_verification.json',dict(status='pass',receipt=str(policy/'equilibrium_receipt.json'),
    sha256=a.digest(policy/'equilibrium_receipt.json'),driver_sha256=original['finalizer_sha256'],
    receipt_working_directory=str(old),verification_scope='Completed full-loop smoke job17090362; original receipt unchanged'))
a.write_json(imported/'combined_verification.json',dict(status='pass',history_component_job=17090362,policy_component_job=17090362,
    scope='All required historical and policy components completed in the original single smoke controller; independent verification below',
    original_contract=str(original_path),original_contract_sha256=original_sha,history_receipt_sha256=case_hashes,
    supported_anchor_cross_source_comparison=source_comparison,cross_source_pngs=17,
    failed_case_replay_provenance_sha256=a.digest(replay/'replay_provenance.json'),failed_case_replay_failure_sha256=a.digest(replay/'task_026/adapter_failure.json'),
    failed_case_replay_assessment_sha256=a.digest(root/'replay_assessment.json'),
    canary_assessment_sha256=a.digest(root/'canary_assessment.json'),
    historical_proof_sha256=a.digest(hist/'smoke_verification.json'),policy_receipt_sha256=a.digest(policy/'equilibrium_receipt.json'),
    policy_artifact_sha256=policy_artifacts,history_seconds=elapsed,policy_smoke_seconds=pr['elapsed_seconds'],
    expected_full_policy_seconds_at_smoke_rate=pr['elapsed_seconds']*44/8,production_promoted=False))
subprocess.run([sys.executable,str(root/'code/model/tools/build_e5f_joint_nested_long_contract.py'),
    '--remote-root',str(root),'--seed-summary',str(anchors[0]/'summary.json'),'--outdir',str(out),
    '--finish-epoch','1788788100','--expected-history-seconds',str(max(elapsed)),
    '--runtime-estimate-status','measured','--profile','wide32',
    '--imported-smoke-root',str(imported),'--imported-smoke-contract',str(original_path),
    '--imported-smoke-contract-sha256',original_sha,
    '--imported-smoke-verification-sha256',a.digest(imported/'smoke_verification.json'),
    '--imported-policy-verification-sha256',a.digest(imported/'policy_loop_verification.json')],check=True)
contract=out/'contract.json';c=search.verify_contract(contract,a.digest(contract))
run=search.Search(c,'preflight');run.require_smoke();run.summary('imported_components_verified')
print(json.dumps(dict(status='preflight_passed',contract=str(contract),contract_sha256=a.digest(contract),
    budget=c['budget_estimate'],history_cases=run.completed,production_promoted=False)),flush=True)
