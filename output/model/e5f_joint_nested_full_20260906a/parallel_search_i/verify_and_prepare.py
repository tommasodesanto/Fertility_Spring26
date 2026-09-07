"""Verify completed repaired smoke and original-case replay, then prepare; never submit."""
from pathlib import Path
import json,math,subprocess,sys
root=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907i')
fresh=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907h')
old=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907e')
sys.path[:0]=[str(root/'code/model'),str(root/'code/model/tools')]
import run_e5f_joint_nested_long_search as search
import run_e5f_transition_calibration as cal
from intergen_eqscale_seq_optimized import solver
a=search.adapter
for name,sha in a.read_json(root/'source_manifest.json').items():a.verify(root/name,sha)
assert cal.code_fingerprint_contract(solver)['bundle_sha256']==a.BUNDLE
if '--source-only' in sys.argv:
    subprocess.run([sys.executable,str(root/'code/model/tools/test_e5f_joint_nested_long_search.py')],check=True)
    report=dict(status='source_and_controller_tests_passed',file_count=70,source_bundle=a.BUNDLE,
        source_manifest_sha256=a.digest(root/'source_manifest.json'),controller_tests=18,
        complete_scientific_smoke_verified=False,contract_created=False)
    a.write_json(root/'source_preflight.json',report);print(json.dumps(report));sys.exit(0)
original_path=fresh/'output/model/joint_nested_overnight/contract.json'
original_sha='117169ecdc8deaca8aa29433a10464e0b32384c4b10dfde0734464914ce021f5'
a.verify(original_path,original_sha);original=a.read_json(original_path)
assert original['policy_workers']==4
assert original['code_bundle_sha256']==a.BUNDLE
hist=Path(original['output_root'])/'smoke';policy=hist.parent/'policy_loop_smoke'
assert a.read_json(hist/'search_state.json')['status']=='smoke_passed'
# The implemented fix must have passed the exact original failed state.
replay=hist.parent/'mass_repair_replay';repair=a.read_json(replay/'saved_state_repair_verification.json')
assert repair['status']=='pass_saved_original_failure_state' and repair['new_bundle']==a.BUNDLE
assert repair['saved_state_sha256']=='4f8e69c794958f45dc0b14a5bfc64415f20a969a869b20b632bfa02626965835'
for gate in repair['gates']:
    assert gate['relative_tolerance']==5e-9 and gate['relative_gap']<=5e-9
provenance=a.read_json(replay/'replay_provenance.json')
assert provenance['new_bundle']==a.BUNDLE
assert provenance['original_plan_sha256']=='3884c6202ecefdc018ae2a53d4f4f506731708524dc79132a7b3877ab6fc9a71'
assert provenance['original_center_sha256']=='40eb3f4df91b5de25925112c12a5aee94464ca2220c9e5042939a3543d75a0c4'
replay_assessment={'original_state_verified':True,'provenance_sha256':a.digest(replay/'replay_provenance.json')}
if provenance['complete']:
    rp=a.load_plan(replay/'plan.json',provenance['new_plan_sha256']);case=rp['cases'][0]
    rr=a.read_json(replay/'task_018/case_receipt.json');assert rr['status']=='complete'
    for rel,sha in rr['artifact_sha256'].items():a.verify(replay/'task_018'/rel,sha)
    summary,_,_=a.validate_result(replay/'task_018',rp,case)
    replay_assessment.update(status='complete_valid_original_case',loss=summary['best_candidate']['transition_loss'])
else:
    failure=a.read_json(replay/'task_018/adapter_failure.json')
    assert search.classify_failure(failure['error'],failure['type'])=='undefined_first_birth_support'
    replay_assessment.update(status='original_controlled_historical_support_rejection_after_repair',failure=failure)
# Preserve and verify every original historical artifact, then compare numeric
# results and the stable graph set with the preceding complete source.
proof=a.read_json(hist/'smoke_verification.json')
assert proof['status']=='pass' and len(proof['anchor_receipts'])==len(proof['probe_receipts'])==2
oldpaths={'anchor_1':'smoke_anchor/task_001','anchor_2':'smoke_anchor/task_002',
          'all_minus':'smoke_all_coordinate_probes/task_001','all_plus':'smoke_all_coordinate_probes/task_002'}
elapsed=[];hashes={};comparisons={};anchors=[]
for value in proof['anchor_receipts']+proof['probe_receipts']:
    path=Path(value);receipt=a.read_json(path);dest=path.parent
    plan=a.load_plan(dest.parent/'plan.json',receipt['plan_sha256'])
    case=next(c for c in plan['cases'] if c['id']==receipt['case_id'])
    a.verify(dest.parent/case['center'],case['center_sha256'])
    for rel,sha in receipt['artifact_sha256'].items():a.verify(dest/rel,sha)
    _,_,_=a.validate_result(dest,plan,case)
    assert receipt['status']=='complete' and receipt['budget_diagnostic']['budget_excess_mass']<=2e-10
    assert receipt['policy_array_diagnostic']['occupied_negative_steps']==0
    ref=old/'output/model/joint_nested_overnight/smoke'/oldpaths[case['label']]
    comparisons[case['label']]=a.compare_reference(dest,ref/'summary.json')
    assert a.read_csv(dest/'parameter_table.csv')==a.read_csv(ref/'parameter_table.csv')
    graphs=list((dest/'standard_diagnostics').glob('*.png'));assert len(graphs)==17
    for graph in graphs:a.verify(ref/'standard_diagnostics'/graph.name,a.digest(graph))
    hashes[value]=a.digest(path);elapsed.append(receipt['elapsed_seconds'])
    if case['label'].startswith('anchor'):anchors.append(dest)
# Policy cases use the same selected state and economics as the serial smoke.
pp=a.read_json(hist/'policy_loop_verification.json');a.verify(Path(pp['receipt']),pp['sha256'])
pr=a.read_json(policy/'equilibrium_receipt.json')
search.validate_policy_receipt(pr,original,smoke=True,selected_hashes={a.digest(p/'summary.json') for p in anchors})
a.verify(fresh/pr['selected_summary'],pr['selected_summary_sha256'])
a.verify(policy/'inherited_state_verification.json',pr['inherited_state_verification_sha256'])
a.verify(root/'code/model/tools/run_e5f_joint_nested_finalize.py',original['finalizer_sha256'])
oldpolicy=old/'output/model/joint_nested_overnight/policy_loop_smoke'
oldinventory=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907g/policy_source_artifact_hashes.json')
a.verify(oldinventory,'9c47c9852774a5cc25ac07cbcff1a312dc6c905793bc69877ea130aab570c058')
for rel,sha in a.read_json(oldinventory)['artifact_sha256'].items():a.verify(oldpolicy/rel,sha)
policyhashes={};numeric=0
for name,receipt in pr['cases'].items():
    folder=policy/name;ref=oldpolicy/name
    assert a.read_json(folder/'receipt.json')==receipt
    rows=a.read_csv(folder/'policy_path.csv');before=a.read_csv(ref/'policy_path.csv')
    assert rows==before  # Exact complete CSV content, including queue strings.
    assert [int(float(r['calendar_year'])) for r in rows]==[2023,2027]
    for row in rows:
        for val in row.values():
            try:float(val);numeric+=1
            except ValueError:pass
    for year in (2023,2027):
        date=folder/f'date_{year}';prior=ref/f'date_{year}'
        budget=a.read_json(date/'budget_summary.json');arrays=a.read_json(date/'policy_array_summary.json')
        assert budget['budget_excess_mass']<=2e-10 and arrays['occupied_negative_steps']==0
        for b in arrays['probabilities'].values():assert b['nonfinite']==0 and 0<=b['minimum']<=b['maximum']<=1
        graphs=list((date/'standard_diagnostics').glob('*.png'));assert len(graphs)==17
        for graph in graphs:a.verify(prior/'standard_diagnostics'/graph.name,a.digest(graph))
        for f in (date/'dated_state.pkl.gz',date/'budget_summary.json',date/'policy_array_summary.json',*graphs):policyhashes[str(f)]=a.digest(f)
    for f in (folder/'receipt.json',folder/'policy_path.csv'):policyhashes[str(f)]=a.digest(f)
for f in (policy/'equilibrium_receipt.json',policy/'inherited_state_verification.json'):policyhashes[str(f)]=a.digest(f)
assert len(policyhashes)==170
out=root/'output/model/joint_nested_overnight';out.mkdir(parents=True,exist_ok=True)
# Directly import immutable h proofs, whose policy paths are absolute.
subprocess.run([sys.executable,str(root/'code/model/tools/build_e5f_joint_nested_long_contract.py'),
 '--remote-root',str(root),'--seed-summary',str(sorted(anchors)[0]/'summary.json'),'--outdir',str(out),
 '--finish-epoch','1788788100','--expected-history-seconds',str(max(elapsed)),
 '--runtime-estimate-status','measured','--profile','parallel32','--policy-workers','4',
 '--parallel-policy-receipt',str(policy/'equilibrium_receipt.json'),
 '--imported-smoke-root',str(hist),'--imported-smoke-contract',str(original_path),
 '--imported-smoke-contract-sha256',original_sha,'--imported-smoke-verification-sha256',a.digest(hist/'smoke_verification.json'),
 '--imported-policy-verification-sha256',a.digest(hist/'policy_loop_verification.json')],check=True)
contract=out/'contract.json';c=search.verify_contract(contract,a.digest(contract));run=search.Search(c,'preflight');run.require_smoke()
assert run.completed==4 and len(run.ledger)==4 and run.best is not None
run.summary('imported_repaired_smoke_verified')
report=dict(status='preflight_passed',contract_sha256=a.digest(contract),source_manifest_sha256=a.digest(root/'source_manifest.json'),
 source_bundle=a.BUNDLE,original_case_replay=replay_assessment,history_receipt_sha256=hashes,history_seconds=elapsed,
 cross_source_comparisons=comparisons,history_exact_pngs=68,policy_artifact_sha256=policyhashes,
 policy_exact_pngs=136,policy_exact_numeric_entries=numeric,policy_smoke_seconds=pr['elapsed_seconds'],
 policy_full_projection_seconds=pr['elapsed_seconds']*44/8,budget_estimate=c['budget_estimate'],production_promoted=False)
a.write_json(out/'preflight_verification.json',report);print(json.dumps(report),flush=True)
