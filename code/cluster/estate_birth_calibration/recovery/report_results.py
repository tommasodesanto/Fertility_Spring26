"""Audit and render the terminal Estate-A recovery collection; no model solves."""
import csv,hashlib,json,math,tarfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
OUT=ROOT/'output/model/experiments/birth_count_choice/estate_a_recovery_20261004_v1/collection'
EXPECTED_TARGET='c7a3d185668122e508a6c322bc5ef0715ebb0ecb23948c8d9b184ee25d1cde70'
EXPECTED_WEIGHT='f762ebb5684ab30487b3b8b64fc10977fda396b520035d91c0c5c803255f88e4'
EXPECTED_STAGE='974a14e243da6a2ad0572bb9825b47ab349828f9144cbbad69e74b40a4408b22'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))
def require(v,m):
    if not v:raise RuntimeError(m)
def fmt(v):
    if v in ('',None):return '—'
    try:return f'{float(v):.8g}'
    except (TypeError,ValueError):return str(v)
def csvwrite(p,rows_):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows_[0]));w.writeheader();w.writerows(rows_)
def main():
    if not (OUT/'receipts').exists():
        (OUT/'receipts').mkdir()
        with tarfile.open(OUT/'all_receipts.tar.gz') as archive: archive.extractall(OUT/'receipts')
    coll=read(OUT/'final_collection.json')
    require(coll['inventory_sha256']==EXPECTED_STAGE and coll['target_fingerprint']==EXPECTED_TARGET and coll['weight_fingerprint']==EXPECTED_WEIGHT,'Stage/target/weight drift')
    require(coll['job_id']=='19141024' and len(coll['chains'])==15,'Collection task drift')
    by={(r['arm'],r['chain']):r for r in coll['chains']}
    require(len(by)==15 and {(('binary' if t<5 else 'count3'),(t if t<5 else t-5)) for t in range(15)}==set(by),'Task map drift')
    verified=[r for r in coll['chains'] if r['status']=='selected_numerically_verified']
    require(len(verified)==11,'Verified count drift')
    failed=[r for r in coll['chains'] if r['status']!='selected_numerically_verified']
    require({(r['arm'],r['chain']) for r in failed}=={('count3',j) for j in (0,1,4,6)},'Failure mapping drift')
    for t in range(15):
        arm='binary' if t<5 else 'count3';chain=t if t<5 else t-5
        folder=OUT/'receipts'/f'task_{t:02d}';run=folder/'run';r=by[arm,chain]
        start=read(folder/'launcher_start.json');terminal=read(folder/'launcher_terminal.json');contract=read(run/'start_contract.json')
        require(start['mode']=='production' and start['arm']==arm and start['chain']==chain and start['stage_inventory_sha256']==EXPECTED_STAGE,'Launcher source/task drift')
        require(start['deadline_epoch']<=1791122400 and start['cpus']==1 and start['memory_GiB']==24,'Resource/deadline drift')
        require(terminal['slurm_job_id']==start['slurm_job_id'] and terminal['arm']==arm and terminal['chain']==chain,'Terminal identity drift')
        require(contract['target_fingerprint']==EXPECTED_TARGET and contract['weight_fingerprint']==EXPECTED_WEIGHT and len(contract['bounds'])==10,'Objective contract drift')
        if r['status']=='selected_numerically_verified':
            require(terminal['exit_code']==0 and read(run/'completed.json')=={k:v for k,v in r.items() if k!='task'},'Verified terminal/completed drift')
            require(read(run/'native_postcheck/completed.json')['search_receipt_sha256']==sha(run/'search_completed.json'),'Fresh child search receipt drift')
            require(math.isclose(r['native_loss'],sum(float(z['loss_contribution'] or 0) for z in r['target_fit']),rel_tol=0,abs_tol=1e-8),'Loss sum drift')
            for z in r['target_fit']:
                if z['weight'] and z['loss_contribution']:
                    require(math.isclose(float(z['weight'])*float(z['gap'])**2,float(z['loss_contribution']),rel_tol=0,abs_tol=1e-7),'Row weighted loss drift')
            require(len(r['target_fit'])==14 and len(r['parameters'])==31,'14/31 drift')
            for key,bounds in contract['bounds'].items():
                match=[p for p in r['parameters'] if p['parameter']==key]
                require(len(match)==1 and [float(match[0]['lower']),float(match[0]['upper'])]==bounds,'Free bound drift: '+key)
                require(float(match[0]['estimate'])==r['selected']['parameters'][key],'Free estimate drift: '+key)
        else:
            require(terminal['exit_code']!=0 and not (run/'completed.json').exists(),'Failed chain presented as completed')
            failure=read(run/'failure.json')
            require(failure['type']=='RuntimeError' and failure['message']=='native GE acceptance failed: uncomputed_bounded_budget','Failure type drift')
            best=read(run/'best_so_far.json')['best']
            require(best['status']=='passed' and math.isfinite(float(best['loss'])),'Missing provisional best')
    winners={arm:min((r for r in verified if r['arm']==arm),key=lambda r:r['native_loss']) for arm in ('binary','count3')}
    require(winners['binary']['chain']==1 and winners['count3']['chain']==3,'Winner chain drift')
    audit={}
    for arm,r in winners.items():
        base=OUT/arm;root=base/'selected_root';repeat=base/'selected_repeat_final'
        fit=rows(root/'target_fit_new_contract.csv');repfit=rows(repeat/'target_fit_new_contract.csv')
        params=rows(root/'parameters_estate_a.csv');repparams=rows(repeat/'parameters.csv')
        require(fit==repfit==r['target_fit'] and params==r['parameters'],'Actual selected/repeat report rows drift')
        for p in repparams:
            if p['parameter']=='beta_annual':
                p['lower'],p['upper']='.93','.99'
                p['near_bound']=str(min(float(p['estimate'])-.93,.99-float(p['estimate']))<=.0006)
        require(repparams==params,'Actual exact-repeat parameters drift')
        rootplots={p.name:sha(p) for p in sorted((root/'standard_diagnostics').glob('*.png'))}
        repeatplots={p.name:sha(p) for p in sorted((repeat/'standard_diagnostics').glob('*.png'))}
        require(len(rootplots)==17 and rootplots==repeatplots==r['repeat']['standard_plot_hashes'],'Actual 17 plot hashes drift')
        require(r['repeat']['status']=='exact_full_ge_repeat_passed' and r['selected_postcheck']['status']=='passed','Fresh native gate drift')
        child=read(base/'provenance/native_postcheck/completed.json');search=read(base/'provenance/search_completed.json')
        require(child['search_receipt_sha256']==sha(base/'provenance/search_completed.json') and child['status']=='full_native_postcheck_passed','Selected child provenance drift')
        require(search['selected']==r['selected'] and search['starts_file_sha256']==r['starts_file_sha256'],'Selected search provenance drift')
        require(all(abs(float(a)-float(b))<=1e-10 for a,b in zip(r['selected']['residual'],r['selected_postcheck']['residual'])),'Fresh residual drift')
        audit[arm]=dict(chain=r['chain'],native_loss=r['native_loss'],search_stop_reason=r['search_stop_reason'],objective_calls=r['objective_calls'],target_rows=len(fit),parameter_rows=len(params),actual_matching_plot_hashes=len(rootplots),actual_target_rows_exact_repeat=True,actual_parameter_rows_exact_repeat=True,fresh_child_search_sha256=child['search_receipt_sha256'])
    fitjoint=[]
    for index,(a,b) in enumerate(zip(winners['binary']['target_fit'],winners['count3']['target_fit'])):
        require(a['moment']==b['moment'] and a['target']==b['target'] and a['weight']==b['weight'] and a['role']==b['role'],'Joint target row drift')
        fitjoint.append(dict(moment=a['moment'],role=a['role'],target=a['target'],weight=a['weight'],binary_model=a['model'],binary_gap=a['gap'],binary_loss=a['loss_contribution'],count3_model=b['model'],count3_gap=b['gap'],count3_loss=b['loss_contribution']))
    csvwrite(OUT/'target_fit_joint.csv',fitjoint)
    paramjoint=[]
    free=set(read(OUT/'receipts/task_01/run/start_contract.json')['bounds'])
    for a,b in zip(winners['binary']['parameters'],winners['count3']['parameters']):
        require(a['parameter']==b['parameter'] and a['status']==b['status'],'Joint parameter drift')
        paramjoint.append(dict(parameter=a['parameter'],status='estimated (free; recovery search bounds)' if a['parameter'] in free else a['status'],binary_estimate=a['estimate'],binary_lower=a['lower'],binary_upper=a['upper'],binary_near_bound=a['near_bound'],count3_estimate=b['estimate'],count3_lower=b['lower'],count3_upper=b['upper'],count3_near_bound=b['near_bound']))
    csvwrite(OUT/'parameters_joint.csv',paramjoint)
    outcomes=[]
    for t in range(15):
        arm='binary' if t<5 else 'count3';chain=t if t<5 else t-5;r=by[arm,chain]
        best=read(OUT/'receipts'/f'task_{t:02d}/run/best_so_far.json')['best']
        outcomes.append(dict(task=t,arm=arm,chain=chain,status=r['status'],verified_native_loss=r['native_loss'] if r['status']=='selected_numerically_verified' else '',provisional_saved_best_loss=best['loss'],objective_calls=r.get('objective_calls',''),failure_message=r.get('failure',{}).get('message','') if isinstance(r.get('failure'),dict) else ''))
    csvwrite(OUT/'chain_outcomes.csv',outcomes)
    auditreceipt=dict(status='actual_native_selected_and_exact_repeat_verified',job_id='19141024',verified_chains=11,failed_chains=4,target_fingerprint=EXPECTED_TARGET,weight_fingerprint=EXPECTED_WEIGHT,inventory_sha256=EXPECTED_STAGE,winners=audit,no_adoption=True)
    (OUT/'audit_receipt.json').write_text(json.dumps(auditreceipt,indent=2)+'\n')
    lines=['# Estate-A recovery result, October 4, 2026','',
        'Fifteen isolated chains ran under the experimental Estate-A/new-wealth contract with unchanged economics, 14 reported moments (ten scored), weights, and ten free-parameter bounds. Eleven passed the full fresh native selected-point and exact-repeat gate; four count-three chains failed with `native GE acceptance failed: uncomputed_bounded_budget`. The latter have only provisional saved best cases. No optimization-convergence certificate, grid certification, rank/identification certification, or scientific adoption is claimed.','',
        f'Best **binary** chain 1: verified loss **{winners["binary"]["native_loss"]:.9f}** after {winners["binary"]["objective_calls"]} calls. Best **count-three** chain 3: verified loss **{winners["count3"]["native_loss"]:.9f}** after {winners["count3"]["objective_calls"]} calls. Both searches stopped at the 800-second prelaunch guard and completed a fresh native child with exact-repeat diagnostics. A larger household choice set does not mathematically nest the old aggregate fit, so these losses alone are not a causal interpretation.','',
        'The table contains every target row. Blank weights and loss contributions are normalization rows; the exact strings are in [target_fit_joint.csv](target_fit_joint.csv).','',
        '| Moment | Role | Target | Weight | Binary model | Binary gap | Binary loss | Count-three model | Count-three gap | Count-three loss |',
        '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in fitjoint:lines.append('| '+' | '.join(str(r[k]) if k in ('moment','role') else fmt(r[k]) for k in r)+' |')
    lines += ['','The next table contains all 31 reported parameters. The ten searched coordinates are exactly those in the recovery start-plan bounds; the other rows are fixed, derived, or retained inputs. `Near` is the native report flag. Full-precision values are in [parameters_joint.csv](parameters_joint.csv).','',
              '| Parameter | Status | Binary estimate | Binary bounds | Binary near | Count-three estimate | Count-three bounds | Count-three near |',
              '|---|---|---:|---:|:---:|---:|---:|:---:|']
    for r in paramjoint:
        bb=f"[{fmt(r['binary_lower'])}, {fmt(r['binary_upper'])}]" if r['binary_lower'] else '—';cb=f"[{fmt(r['count3_lower'])}, {fmt(r['count3_upper'])}]" if r['count3_lower'] else '—'
        lines.append(f"| {r['parameter']} | {r['status']} | {fmt(r['binary_estimate'])} | {bb} | {r['binary_near_bound'] or '—'} | {fmt(r['count3_estimate'])} | {cb} | {r['count3_near_bound'] or '—'} |")
    lines += ['','All chain outcomes are in [chain_outcomes.csv](chain_outcomes.csv). The four failed count-three chains (0, 1, 4, 6) retain provisional best losses 152.924831, 550.605309, 105.082178, and 484.766329, respectively; these values did not pass a final native postcheck and are excluded from winner selection.','',
              'Selected actual 14-row CSVs and 31-row parameter CSVs match their receipts. The selected and exact-repeat target rows match, and the 17 selected PNG hashes per arm match both the repeat PNGs and native receipt. The child search SHA and selected residuals also match. See [audit_receipt.json](audit_receipt.json) and the downloaded [binary](binary/selected_root/standard_diagnostics/fertility_by_age.png) and [count-three](count3/selected_root/standard_diagnostics/fertility_by_age.png) diagnostics. The remaining standard diagnostic PNGs are in adjacent directories.','',
              'Source: Slurm array 19141024; pinned recovery inventory `'+EXPECTED_STAGE+'`; target fingerprint `'+EXPECTED_TARGET+'`; weight fingerprint `'+EXPECTED_WEIGHT+'`. The full native receipts are in `receipts/`, and the strict [final collection](final_collection.json) retains all 15 task outcomes. The original parent optimizer failed; its saved checkpoints were used solely as provisional new search starts. The recovery made a numerical time-budget correction and changed the starts; it did not change earnings, initial wealth, transfers, floors, preferences, targets, weights, or housing/estate equations.']
    (OUT/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(status='report_passed',winners=audit,verified=11,failed=4)))
if __name__=='__main__':main()
