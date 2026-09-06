#!/usr/bin/env python3
"""Summarize completed cases, or regenerate graphs from a verified checkpoint.

Run after the remote collector verifies original checkpoint/graph receipts.
These plots supplement, and do not replace, the seventeen standard diagnostics.
"""
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent
REPORT = OUT/'analysis'


def read_csv(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


def main():
    cases, fits, parameters, inputs = [], [], [], {}
    for stage in ('smoke', 'round1', 'round2', 'repeats'):
        plan_path = OUT/stage/'plan.json'
        if not plan_path.exists():
            continue
        plan = json.loads(plan_path.read_text())
        for case in plan['cases']:
            directory = plan_path.parent/case['output']
            receipt_path = directory/'case_receipt.json'
            if not receipt_path.exists():
                continue
            receipt = json.loads(receipt_path.read_text())
            assert receipt['status']=='complete'
            assert receipt['plan_sha256']==hashlib.sha256(plan_path.read_bytes()).hexdigest()
            for name, expected in receipt['artifact_sha256'].items():
                p = directory/name
                if p.suffix=='.gz':
                    continue  # Large checkpoints are verified by the remote collector.
                assert hashlib.sha256(p.read_bytes()).hexdigest()==expected, p
            f = read_csv(directory/'target_fit_long.csv')
            p = read_csv(directory/'parameter_table.csv')
            assert len(f)==12 and len({r['moment'] for r in f})==12
            assert sum(r['is_free_parameter']=='True' for r in p)==11
            for row in f:
                t,m,g,w,c = (float(row[k]) for k in ('target','model','gap','weight','loss_contribution'))
                assert math.isclose(g,m-t,abs_tol=2e-10)
                assert math.isclose(c,w*g*g,abs_tol=2e-9)
            loss=sum(float(r['loss_contribution']) for r in f)
            assert math.isclose(loss,receipt['loss'],abs_tol=2e-8)
            model={r['moment']:float(r['model']) for r in f}
            theta={r['parameter']:float(r['value']) for r in p}
            source_summary=json.loads((directory/'summary.json').read_text())
            branch_path=directory/'cases'/f[0]['candidate']/'dated_first_birth_housing_did.csv'
            branch=read_csv(branch_path)[-1]
            assert float(branch['calendar_year'])==2023
            response=float(branch['treated_mean_housing'])-float(branch['control_mean_housing'])
            assert math.isclose(response,model['housing_increment_0to1'],abs_tol=2e-11)
            entry={'stage':stage,'case_id':case['id'],'label':case['label'],
                'loss':loss,'first_birth_rooms':response,
                'treated_rooms':float(branch['treated_mean_housing']),
                'control_rooms':float(branch['control_mean_housing']),
                'continuation_births_per_treated_household':float(branch['treated_continuation_births'])/float(branch['destination_mass']),
                'mean_rooms':model['aggregate_mean_occupied_rooms_18_85'],
                'three_plus_rooms_gap':model['prime30_55_parent_3plus_minus_1to2_mean_rooms'],
                'ownership':model['own_rate'],'childlessness':model['childless_rate'],
                'H0':theta['H0'],'first_child_jump':theta['hbar_first_child_jump'],
                'per_child_floor':theta['hbar_child_rooms'],
                'one_child_floor':theta['hbar_first_child_jump']+theta['hbar_child_rooms'],
                'two_child_floor':theta['hbar_first_child_jump']+2*theta['hbar_child_rooms'],
                'three_child_floor':theta['hbar_first_child_jump']+3*theta['hbar_child_rooms'],
                'elapsed_seconds':receipt['elapsed_seconds'],
                'dated_model_solve_count':source_summary['best_candidate']['model_solve_count'],
                'old_normalization_stationary_solves':source_summary['old_fertility_normalization']['stationary_solves'],
                'occupied_value_drops':receipt['policy_array_diagnostic']['occupied_negative_steps'],
                'budget_excess_share':receipt['budget_diagnostic']['budget_excess_share'],
                'directory':str(directory.relative_to(OUT))}
            cases.append(entry)
            for row in f: fits.append({'stage':stage,'case_id':case['id'],'label':case['label'],**row})
            for row in p: parameters.append({'stage':stage,'case_id':case['id'],'label':case['label'],**row})
            for path in (plan_path,receipt_path,branch_path,directory/'summary.json',directory/'target_fit_long.csv',directory/'parameter_table.csv'):
                inputs[str(path.relative_to(OUT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    if not cases:
        raise RuntimeError('No complete collected cases')
    REPORT.mkdir(exist_ok=True)
    anchor=next(r for r in cases if r['stage']=='smoke')
    for row in cases:
        row['change_treated_rooms']=row['treated_rooms']-anchor['treated_rooms']
        row['change_control_rooms']=row['control_rooms']-anchor['control_rooms']
        row['change_first_birth_response']=row['first_birth_rooms']-anchor['first_birth_rooms']
        assert math.isclose(row['change_first_birth_response'],row['change_treated_rooms']-row['change_control_rooms'],abs_tol=2e-11)
    cases.sort(key=lambda r:(r['loss'],r['stage'],r['case_id']))
    write_csv(REPORT/'all_cases.csv',cases)
    write_csv(REPORT/'all_target_fits.csv',fits)
    write_csv(REPORT/'all_parameters.csv',parameters)
    summary={'completed_cases':len(cases),'anchor':anchor,'best':cases[0],
        'largest_housing_response':max(cases,key=lambda r:r['first_birth_rooms']),
        'improving_cases':sum(r['loss']<anchor['loss']-1e-8 for r in cases),
        'input_sha256':inputs,'production_promoted':False,
        'scope':'Discrete observed candidates; no global frontier, causal decomposition, or unreachable-target claim.'}
    coordinate=[r for r in cases if r['stage']=='round1']
    if len(coordinate)==23:
        by_id={r['case_id']:r for r in coordinate}
        s=json.loads((OUT/by_id[1]['directory']/'summary.json').read_text())
        domain=s['panel_design']['domain']
        ordered=read_csv(OUT/by_id[1]['directory']/'target_fit_long.csv')
        names=[r['moment'] for r in ordered]
        weighted={}
        for i,x in by_id.items():
            by_moment={r['moment']:float(r['standardized_gap']) for r in read_csv(OUT/x['directory']/'target_fit_long.csv')}
            assert set(by_moment)==set(names)
            weighted[i]=np.array([by_moment[name] for name in names])
        J=np.empty((12,11));detail=[];columns=[]
        for j,spec in enumerate(domain):
            lower=2+2*j;upper=3+2*j
            for i,direction in ((lower,-1),(upper,1)):
                point=json.loads((OUT/by_id[i]['directory']/'summary.json').read_text())
                step=np.array(point['panel_design']['unit_vector'])-np.array(s['panel_design']['unit_vector'])
                expected=np.zeros(11);expected[j]=direction*.0025
                assert np.allclose(step,expected,rtol=0,atol=2e-15)
            left=(weighted[1]-weighted[lower])/.0025
            right=(weighted[upper]-weighted[1])/.0025
            center=(left+right)/2;J[:,j]=center
            norm=float(np.linalg.norm(center))
            discrepancy=float(np.linalg.norm(right-left)/max(norm,1e-12))
            cosine=float(left@right/max(np.linalg.norm(left)*np.linalg.norm(right),1e-12))
            columns.append(dict(parameter=spec['name'],central_column_norm=norm,
                relative_one_sided_discrepancy=discrepancy,one_sided_cosine=cosine))
            for k,name in enumerate(names):
                detail.append(dict(moment=name,parameter=spec['name'],minus_derivative=left[k],
                    plus_derivative=right[k],central_derivative=center[k],
                    units='weighted residual per transformed unit coordinate'))
        singular=np.linalg.svd(J,compute_uv=False)
        write_csv(REPORT/'coordinate_sensitivity.csv',detail)
        write_csv(REPORT/'coordinate_consistency.csv',columns)
        write_csv(REPORT/'jacobian_singular_values.csv',[dict(index=i+1,singular_value=v,relative=v/singular[0]) for i,v in enumerate(singular)])
        summary['local_sensitivity']=dict(rank_at_relative_1e_minus3=int(np.sum(singular>singular[0]*1e-3)),
            condition_number=float(singular[0]/singular[-1]),
            most_inconsistent_columns=sorted(columns,key=lambda r:r['relative_one_sided_discrepancy'],reverse=True),
            qualification='At this center, step, weighting and parameter transform only; neither global nor causal identification. No Jacobian inversion or new search is performed.')
    (REPORT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for stage,color in (('round1','#4682a9'),('round2','#b96a37')):
        rows=[r for r in cases if r['stage']==stage]
        if rows:
            ax[0].scatter([r['first_birth_rooms'] for r in rows],[r['loss'] for r in rows],label=stage,color=color,s=28)
            ax[1].scatter([r['first_birth_rooms'] for r in rows],[r['mean_rooms'] for r in rows],color=color,s=28)
    for axis in ax:
        axis.axvline(.7202462623815278,color='#777777',ls='--',lw=1,label='Target')
        axis.set_xlabel('First-birth housing response (rooms)')
    ax[0].scatter(anchor['first_birth_rooms'],anchor['loss'],marker='*',s=110,color='black',label='Starting candidate')
    ax[1].scatter(anchor['first_birth_rooms'],anchor['mean_rooms'],marker='*',s=110,color='black')
    ax[1].axhline(5.779970481941968,color='#777777',ls='--',lw=1)
    ax[0].set_ylabel('Unchanged twelve-moment loss');ax[1].set_ylabel('Mean occupied rooms')
    ax[0].legend(fontsize=8)
    fig.suptitle('Supplemental: observed fit trade-offs')
    fig.savefig(REPORT/'observed_tradeoffs.png',dpi=170);plt.close(fig)
    print(json.dumps({k:summary[k] for k in ('completed_cases','best','largest_housing_response','improving_cases')},indent=2))


def verified_packet(relative_case):
    root=OUT.parents[2]
    sys.path[:0]=[str(OUT),str(root/'code/model'),str(root/'code/model/tools')]
    directory=(OUT/relative_case).resolve()
    directory.relative_to(OUT)
    plan_path=directory.parent/'plan.json'
    plan=json.loads(plan_path.read_text())
    if plan['schema']=='e5f_fixed_first_child_profile_v1':
        import fixed_jump_adapter as adapter
    else:
        import run_e5f_bounded_calibration_refinement as adapter
    receipt=adapter.read_json(directory/'case_receipt.json')
    adapter.load_plan(plan_path,receipt['plan_sha256'])
    if receipt['status']!='complete':raise RuntimeError('Incomplete saved case')
    case=next(c for c in plan['cases'] if c['output']==directory.name)
    if case['id']!=receipt['case_id']:raise RuntimeError('Receipt belongs to another case')
    adapter.validate_result(directory,plan,case)
    adapter.verify(directory/'dated_state.pkl.gz',receipt['artifact_sha256']['dated_state.pkl.gz'])
    adapter.verify(root/'code/model/tools/run_e5f_independent_numerical_audit.py',plan['helper_sha256']['run_e5f_independent_numerical_audit.py'])
    import run_e5f_independent_numerical_audit as audit
    import run_e5f_transition_calibration as calibration
    _,configured=audit.transition.configure_sequential_model()
    if calibration.code_fingerprint_contract(configured)['bundle_sha256']!=plan['code_bundle_sha256']:
        raise RuntimeError('Use the frozen scientific snapshot named in the plan')
    return adapter,receipt,audit,audit.load_checkpoint(directory/'dated_state.pkl.gz')


def regenerate_graphs(relative_case, destination):
    adapter,receipt,audit,packet=verified_packet(relative_case)
    if destination.exists():raise FileExistsError(destination)
    destination.mkdir(parents=True)
    audit.standard_diagnostics(packet,destination,validate_production_young=False)
    graphs=list((destination/'standard_diagnostics').glob('*.png'))
    if len(graphs)!=17:raise RuntimeError('Expected the unchanged seventeen graphs')
    exact=all(adapter.digest(p)==receipt['artifact_sha256']['standard_diagnostics/'+p.name] for p in graphs)
    adapter.write_json(destination/'regeneration_receipt.json',dict(model_solves=0,graph_count=17,
        exact_original_png_hashes=exact,checkpoint_sha256=receipt['artifact_sha256']['dated_state.pkl.gz']))
    print(json.dumps(dict(graph_count=17,exact_original_png_hashes=exact,model_solves=0)))


def quantile_check():
    rows=[];hashes={}
    for case_id in (1,2,3,10,11,12,13,14,15):
        relative=f'round1/task_{case_id:03d}'
        adapter,receipt,audit,packet=verified_packet(relative)
        from intergen_eqscale_seq_optimized import solver as model
        # The calibration passes post-fertility, pre-transaction mass as asset_g.
        P=packet['parameters'];e=packet['evaluation'];g=np.asarray(e.g_post_fertility)
        bg=np.asarray(packet['b_grid']);prices=np.asarray(e.policy.price)
        values=[];weights=[]
        # Match add_annual_gross_old_wealth_moments, including its order of cells.
        for j in range(P.J):
            age=float(P.age_start)+j*float(P.da)
            if not 76<=age<=84:continue
            for i in range(P.I):
                for z in range(g.shape[4]):
                    income=model.annual_gross_income_at_state(P,i,j,float(P.z_grid[z]))
                    if not np.isfinite(income) or income<=0:continue
                    for tenure in range(g.shape[1]):
                        house=prices[i]*P.H_own[tenure-1] if tenure else 0.
                        v=(bg+house)/income
                        for n in range(P.n_parity):
                            for c in range(P.n_child_states):
                                mass=g[:,tenure,i,j,z,n,c]
                                positive=np.isfinite(mass)&(mass>0)
                                if positive.any():values.append(v[positive]);weights.append(mass[positive])
        v=np.concatenate(values);w=np.concatenate(weights)
        q50=model.weighted_quantile(v,w,.5);q90=model.weighted_quantile(v,w,.9)
        saved={r['moment']:float(r['model']) for r in read_csv(OUT/relative/'target_fit_long.csv')}
        ratio=q90/max(q50,1e-12)
        assert math.isclose(ratio,saved['old_total_wealth_to_annual_income_p90_p50_7684'],rel_tol=0,abs_tol=2e-11)
        row=dict(case_id=case_id,p50=q50,p90=q90,p90_p50=ratio)
        for q,label in ((q50,'p50'),(q90,'p90')):
            row[label+'_cdf_strictly_below']=float(w[v<q].sum()/w.sum())
            row[label+'_cdf_at_or_below']=float(w[v<=q].sum()/w.sum())
            row[label+'_previous_distinct_value']=float(v[v<q].max()) if np.any(v<q) else None
            row[label+'_next_distinct_value']=float(v[v>q].min()) if np.any(v>q) else None
        rows.append(row);hashes[relative]=receipt['artifact_sha256']['dated_state.pkl.gz']
        print(json.dumps(row),flush=True)
    REPORT.mkdir(exist_ok=True)
    write_csv(REPORT/'quantile_support.csv',rows)
    adapter.write_json(REPORT/'quantile_support_receipt.json',dict(status='complete',model_solves=0,
        checkpoints=len(rows),checkpoint_sha256=hashes,all_saved_tail_ratios_reproduced=True,
        scope='Unchanged inverse-CDF statistic and source; diagnostic reconstruction only'))



def summarize_profile():
    """Keep the conditional fixed-jump diagnostics separate from calibration."""
    rows=[];inputs={}
    references=['profile_smoke/task_001']+[f'profile_grid/task_{i:03d}' for i in (1,2,3)]
    for relative in references:
        directory=OUT/relative;r=json.loads((directory/'case_receipt.json').read_text())
        if r['status']!='complete':raise RuntimeError('Profile point incomplete')
        for name,expected in r['artifact_sha256'].items():
            if not name.endswith('.pkl.gz') and hashlib.sha256((directory/name).read_bytes()).hexdigest()!=expected:
                raise RuntimeError(f'Profile collected hash differs: {relative}/{name}')
        fit=read_csv(directory/'target_fit_long.csv');s=json.loads((directory/'summary.json').read_text())
        theta=s['best_candidate']['theta'];model={x['moment']:float(x['model']) for x in fit}
        branch=read_csv(directory/'cases'/fit[0]['candidate']/'dated_first_birth_housing_did.csv')[-1]
        treated=float(branch['treated_mean_housing']);control=float(branch['control_mean_housing'])
        if abs(treated-control-model['housing_increment_0to1'])>2e-11:raise RuntimeError('Profile branch identity differs')
        j,c=theta['hbar_first_child_jump'],theta['hbar_child_rooms']
        rows.append(dict(directory=relative,jump=j,per_child_floor=c,one_child_floor=j+c,two_child_floor=j+2*c,
            three_child_floor=j+3*c,loss=r['loss'],first_birth_rooms=model['housing_increment_0to1'],
            mean_rooms=model['aggregate_mean_occupied_rooms_18_85'],ownership=model['own_rate'],
            treated_rooms=treated,control_rooms=control,
            continuation_births_per_treated_household=float(branch['treated_continuation_births'])/float(branch['destination_mass']),
            old_psi=s['old_psi_child'],new_psi=s['best_candidate']['new_psi_child'],
            occupied_value_drops=r['policy_array_diagnostic']['occupied_negative_steps'],
            budget_excess_share=r['budget_diagnostic']['budget_excess_share']))
        for p in (directory/'case_receipt.json',directory/'target_fit_long.csv',directory/'summary.json'):
            inputs[str(p.relative_to(OUT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    REPORT.mkdir(exist_ok=True);write_csv(REPORT/'profile_cases.csv',rows)
    (REPORT/'profile_analysis_receipt.json').write_text(json.dumps(dict(input_sha256=inputs,model_solves=0,
        scope='Four conditional points including anchor, not re-estimated calibrations'),indent=2)+'\n')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,3,figsize=(10.3,3.2),layout='constrained')
    for axis,key,label,target in zip(ax,('first_birth_rooms','mean_rooms','loss'),
            ('First-birth response (rooms)','Mean rooms','Twelve-row objective'),(.7202462623815278,5.779970481941968,None)):
        axis.plot([r['jump'] for r in rows],[r[key] for r in rows],'-o',color='#173e55',lw=1.8)
        if target is not None:axis.axhline(target,color='#b96a37',ls='--',lw=1,label='Target');axis.legend(fontsize=8)
        axis.set_xlabel('Fixed first-child floor');axis.set_ylabel(label);axis.set_xticks([.5,.6,.75,.9])
    fig.suptitle('Supplemental: preserve the three-child floor; hold nine coordinates fixed',fontsize=11)
    fig.savefig(REPORT/'conditional_floor_profile.png',dpi=190);plt.close(fig)
    print(json.dumps(rows,indent=2))

def verify_artifacts(local=False):
    """Verify original case receipts remotely, then hashes of collected files."""
    def digest(path):
        h=hashlib.sha256()
        with path.open('rb') as f:
            for b in iter(lambda:f.read(1<<20),b''):h.update(b)
        return h.hexdigest()
    manifest_path=OUT/'artifact_manifest.json'
    if local:
        manifest=json.loads(manifest_path.read_text());checked=0;remote=[]
        for name,expected in manifest['files'].items():
            path=OUT/name
            if not path.exists() and name.endswith('.pkl.gz'):
                remote.append(name);continue
            if digest(path)!=expected:raise RuntimeError(f'Collected artifact differs: {name}')
            checked+=1
        result=dict(status='pass',locally_verified_files=checked,remote_only_checkpoints=remote,
            remote_manifest_sha256=digest(manifest_path),remote_validation=manifest['validation'])
        (OUT/'collection_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(dict(status='pass',locally_verified_files=checked,remote_only_checkpoints=len(remote))))
        return
    root=OUT.parents[2]
    sys.path[:0]=[str(OUT),str(root/'code/model'),str(root/'code/model/tools')]
    import run_e5f_bounded_calibration_refinement as bounded
    import fixed_profile as profile
    import build_e5f_bounded_refinement_plan as planner
    counts={'smoke':2,'round1':23,'round2':12,'repeats':2,
        'profile_smoke':2,'profile_grid':3,'profile_repeats':2}
    stop_path=OUT/'profile_stop.json'
    if stop_path.exists():
        stop=json.loads(stop_path.read_text())
        if stop['status']!='complete_no_improving_extended_point' or stop['production_promoted']:
            raise RuntimeError('Invalid stopped-profile classification')
        if digest(OUT/stop['unused_repeat_plan'])!=stop['unused_repeat_plan_sha256']:
            raise RuntimeError('Unused repeat plan changed')
        if list((OUT/'profile_repeats').glob('task_*/case_receipt.json')):
            raise RuntimeError('Unexpected extra repeat after declared stop')
        best=json.loads((OUT/'profile_best_so_far.json').read_text())
        grid=read_csv(OUT/'profile_grid/report/all_candidates.csv')
        if best['stage']!='conditional_profile_two_exact_smokes' or min(float(r['loss']) for r in grid)<=best['best']['loss']:
            raise RuntimeError('Cannot omit final repeats for an improved profile point')
        del counts['profile_repeats']
    receipts=[]
    for stage,expected in counts.items():
        path=OUT/stage/'plan.json'
        collector=profile.collect if stage.startswith('profile') else planner.collect
        _,status,rows=collector(path,digest(path))
        if status['status']!='complete' or len(rows)!=expected:
            raise RuntimeError(f'Incomplete planned stage {stage}')
        found=sorted((OUT/stage).glob('task_*/case_receipt.json'))
        if len(found)!=expected:raise RuntimeError(f'Unexpected case receipts in {stage}')
        receipts.extend(found)
    followup=OUT/'empirical_room_gap'
    if (followup/'grid/plan.json').exists():
        import collect_empirical_room_gap as gap
        stages=[('grid',4)]
        if not (followup/'stop.json').exists():stages.append(('repeats',2))
        for stage,expected in stages:
            path=followup/stage/'plan.json'
            _,status,rows=gap.collect(path,digest(path))
            if status['status']!='complete' or len(rows)!=expected:raise RuntimeError('Incomplete observed-gap follow-up')
            found=sorted((followup/stage).glob('task_*/case_receipt.json'))
            if len(found)!=expected:raise RuntimeError('Unexpected observed-gap case count')
            receipts.extend(found)
        if (followup/'stop.json').exists():
            stop=json.loads((followup/'stop.json').read_text())
            anchor=json.loads((followup/'grid/plan.json').read_text())['anchor_loss']
            if stop['status']!='complete_no_improvement' or stop['best']['loss']<anchor:
                raise RuntimeError('Invalid follow-up stop')
    original=0
    for path in receipts:
        r=json.loads(path.read_text())
        if len(r['artifact_sha256'])!=21 or r['standard_graph_count']!=17:
            raise RuntimeError('Incomplete case artifact set')
        for name,expected in r['artifact_sha256'].items():
            if digest(path.parent/name)!=expected:raise RuntimeError(f'Original hash differs: {path.parent/name}')
            original+=1
        if r['reference'] and (not r['reference']['exact_twelve_row_fit'] or not r['reference']['exact_parameters'] or r['reference']['exact_numeric_history_entries']!=253):
            raise RuntimeError('Incomplete repeat verification')
    saving=json.loads((OUT/'saving_flag_check/summary.json').read_text())
    if saving['status']!='complete':raise RuntimeError('Fixed-price saving check incomplete')
    folders=tuple(counts)+('saving_flag_check','analysis','graph_replay_smoke','graph_replay_profile','empirical_room_gap','graph_replay_empirical')
    files={}
    for folder in folders:
        for path in sorted((OUT/folder).rglob('*')):
            if path.is_file() and path.suffix in ('.json','.csv','.png','.pdf','.gz'):
                files[str(path.relative_to(OUT))]=digest(path)
    validation=dict(status='pass',historical_cases=len(receipts),profile_stop=json.loads(stop_path.read_text()) if stop_path.exists() else None,original_recorded_hashes_verified=original,
        fixed_price_solve_packets=2,standard_graphs=17*(len(receipts)+2),
        manifest_files=len(files),analysis_driver_sha256=digest(Path(__file__)),
        scope='Original case receipts and complete collected result folders; fixed-price output hashes recorded at final collection')
    manifest_path.write_text(json.dumps(dict(validation=validation,files=files),indent=2,sort_keys=True)+'\n')
    print(json.dumps(validation))


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--graphs',help='Saved case relative to this experiment directory')
    parser.add_argument('--graphs-out',type=Path)
    parser.add_argument('--quantile-check',action='store_true')
    parser.add_argument('--profile-summary',action='store_true')
    parser.add_argument('--verify-artifacts',action='store_true')
    parser.add_argument('--verify-collected',action='store_true')
    args=parser.parse_args()
    if args.profile_summary:summarize_profile()
    elif args.verify_artifacts:verify_artifacts()
    elif args.verify_collected:verify_artifacts(local=True)
    elif args.quantile_check:quantile_check()
    elif args.graphs:
        if args.graphs_out is None:parser.error('--graphs-out is required with --graphs')
        regenerate_graphs(args.graphs,args.graphs_out)
    else:main()
