#!/usr/bin/env python3
"""Torch-only, bounded Jacobian and conditional reoptimization experiment.

No household or equilibrium equations are implemented here. Evaluation, failure
classification and authentication are delegated to the pinned night controller
and its unchanged evening runtime. Every probe normalizes completed fertility.
"""
from __future__ import annotations
import argparse, csv, json, math, os, random, shutil, sys, time
from pathlib import Path
import run_e5f_night_calibration as core

FREE = core.FREE
LANES = ('primary', 'early10', 'early100', 'profile_half', 'profile_double', 'profile_quadruple')
SCHEMA = 'e5f_fertility_identification_v1'
read, write, sha, canon = core.read, core.write, core.sha, core.canon


def verify(path):
    c = read(path)
    assert not sys.flags.optimize and c['schema'] == SCHEMA
    assert os.environ.get('EXPECTED_E5F_IDENTIFICATION_SHA256') == sha(path)
    for pin in list(c['files'].values()) + [c['source_contract'], c['source_manifest']]:
        assert sha(pin['path']) == pin['sha256'], pin['path']
    assert sha(__file__) == c['files']['driver']['sha256']
    assert sha(core.__file__) == c['files']['night_helper']['sha256']
    old = read(c['source_contract']['path'])
    for name in ('source_manifest', 'source_root', 'native_ancestry_contract', 'reference_case',
                 'normalization', 'fixed', 'validation_rows', 'standard_diagnostic_names', 'log_parameters'):
        assert c[name] == old[name], 'Scientific invariant changed: ' + name
    for name, pin in old['files'].items():
        if name != 'driver': assert c['files'][name] == pin
    objs = {lane: read(c['lanes'][lane]['objective']['path']) for lane in LANES}
    assert set(c['lanes']) == set(LANES) and set(c['initial_point']) == set(FREE)
    reference = read(old['lanes']['primary']['objective']['path'])
    assert sha(old['lanes']['primary']['objective']['path']) == old['lanes']['primary']['objective']['sha256']
    assert objs['primary'] == reference
    early = c['early_moment']
    base_rows = {r['restriction_id']: r for r in reference['target_rows']}
    assert early in base_rows and base_rows[early]['actual_weight'] > 0
    assert len(base_rows) == 14
    bounds = {r['parameter']: (r['lower'], r['upper']) for r in reference['parameter_restrictions']}
    assert set(bounds) == set(FREE)
    multipliers = dict(zip(LANES, (1, 10, 100, 10, 10, 10)))
    factors = dict(zip(LANES[3:], (.5, 2., 4.)))
    for lane, obj in objs.items():
        spec = c['lanes'][lane]
        assert sha(spec['objective']['path']) == spec['objective']['sha256']
        assert canon(obj) == spec['canonical_sha256']
        assert canon(obj['target_rows']) == spec['target_weight_fingerprint']
        assert obj['parameter_restrictions'] == reference['parameter_restrictions']
        assert {k:v for k,v in obj.items() if k != 'target_rows'} == {k:v for k,v in reference.items() if k != 'target_rows'}
        assert spec['early_multiplier'] == multipliers[lane]
        assert len(obj['target_rows']) == 14
        for row in obj['target_rows']:
            original = base_rows[row['restriction_id']]
            assert {k:v for k,v in row.items() if k != 'actual_weight'} == {k:v for k,v in original.items() if k != 'actual_weight'}
            expected = original['actual_weight']
            if row['restriction_id'] == early: expected *= multipliers[lane]
            assert row['actual_weight'] == expected
        expected_fixed = {} if lane not in factors else {'kappa_fert_continuation': c['initial_point']['kappa_fert_continuation'] * factors[lane]}
        assert spec['fixed'] == expected_fixed
        for k, value in spec['fixed'].items(): assert bounds[k][0] <= value <= bounds[k][1]
    for k, value in c['initial_point'].items():
        lo, hi = bounds[k]; step = c['jacobian_steps'][k]
        assert math.isfinite(value) and math.isfinite(step) and step > 0 and lo <= value-step < value+step <= hi
    assert set(c['jacobian_steps']) == set(FREE)
    assert set(c['proposal_widths']) == set(FREE) and all(math.isfinite(x) and x > 0 for x in c['proposal_widths'].values())
    b = c['budget']
    assert 1 <= b['workers'] <= 24 and b['objective_cap_seconds'] == 1800
    assert b['population'] == 16 and b['generations'] == 4 and b['max_search_cases'] == 480
    assert b['max_objective_cases'] == 538
    assert (b['total_seconds'], b['search_seconds'], b['repeat_seconds']) == (21600, 18000, 21000)
    core.verify_anchor_pins(c)
    assert read(c['anchor']['receipt']['path'])['point'] == c['initial_point']
    assert len(c['seed_points']) == len(c['seed_artifacts']) >= 2
    assert c['seed_points'][0] == c['initial_point']
    for point, artifact in zip(c['seed_points'], c['seed_artifacts']):
        assert sha(artifact['path']) == artifact['sha256']
        assert read(artifact['path'])['point'] == point and set(point) == set(FREE)
        assert all(bounds[k][0] <= v <= bounds[k][1] for k,v in point.items())
    return c, objs


def probe_points(c):
    result = []
    for name in FREE:
        for scale in (1., .5):
            for sign in (-1, 1):
                point = dict(c['initial_point']); point[name] += sign*scale*c['jacobian_steps'][name]
                result.append(dict(lane='primary', point=point, design=f'jacobian:{name}:{scale}:{sign}', coordinate=name, scale=scale, sign=sign))
    return result


def transformed(point, c):
    return {k: math.log(v) if k in c['log_parameters'] else v for k,v in point.items()}


def from_transformed(point, c, obj, lane):
    bounds = {r['parameter']: (r['lower'],r['upper']) for r in obj['parameter_restrictions']}
    result = {}
    for k, x in point.items():
        lo,hi = bounds[k]
        if k in c['log_parameters']: x = math.exp(min(math.log(hi), max(math.log(lo), x)))
        result[k] = min(hi, max(lo, x))
    result.update(c['lanes'][lane]['fixed'])
    return result


def initial_population(c, obj, lane):
    rng = random.Random(c['seed'] + LANES.index(lane)*100003)
    points = []
    for i in range(c['budget']['population']):
        p = transformed(c['seed_points'][i % len(c['seed_points'])], c)
        if i >= len(c['seed_points']):
            for k in FREE: p[k] += rng.gauss(0, c['proposal_widths'][k]*(.5 if i%2 else 1.))
        points.append(from_transformed(p,c,obj,lane))
    assert len({canon(p) for p in points}) == len(points), 'Initial population collapsed at bounds'
    return points


def de_trials(c, obj, lane, population, generation):
    """DE/rand/1/bin; failed parents retain locations with infinite selection cost."""
    rng = random.Random(c['seed'] + LANES.index(lane)*100003 + generation*1009)
    free = [k for k in FREE if k not in c['lanes'][lane]['fixed']]
    transformed_pop = [transformed(row['point'],c) for row in population]
    trials = []
    for i, old in enumerate(transformed_pop):
        a,b,d = rng.sample([j for j in range(len(population)) if j != i], 3)
        forced = rng.choice(free); p = dict(old)
        for k in free:
            if k == forced or rng.random() < .8:
                p[k] = transformed_pop[a][k] + .7*(transformed_pop[b][k]-transformed_pop[d][k])
        trials.append(from_transformed(p,c,obj,lane))
    return trials


def choose(parent, trial):
    if trial['status'] == 'success' and (parent['status'] != 'success' or trial['loss'] <= parent['loss']): return trial
    return parent


def jacobian_export(c, records, output):
    indexed = {r['design']: r for r in records if r['design'].startswith('jacobian:')}
    fields = ['parameter','moment','step','full_derivative','half_derivative','difference','full_curvature','half_curvature','status']
    rows = []
    anchor_table = core.table(Path(c['anchor']['case_path'])/'target_fit.csv','moment')
    moments = list(anchor_table) + ['psi_child']
    anchor_values = {m:float(r['model']) for m,r in anchor_table.items()}
    anchor_values['psi_child'] = read(c['anchor']['receipt']['path'])['normalization']['psi_child']
    for name in FREE:
        values = {}; curvatures = {}; statuses = []
        for scale in (1., .5):
            pair = [indexed.get(f'jacobian:{name}:{scale}:{sign}') for sign in (-1,1)]
            if all(r and r['status']=='success' for r in pair):
                tables = [core.table(Path(r['case_path'])/'target_fit.csv','moment') for r in pair]
                pair_values = [{m:float(t[m]['model']) for m in anchor_table} for t in tables]
                for r,v in zip(pair,pair_values):v['psi_child']=read(Path(r['case_path'])/'receipt.json')['normalization']['psi_child']
                h=scale*c['jacobian_steps'][name]
                values[scale] = {m:(pair_values[1][m]-pair_values[0][m])/(2*h) for m in moments}
                curvatures[scale] = {m:(pair_values[1][m]-2*anchor_values[m]+pair_values[0][m])/(h*h) for m in moments}
            else: statuses.append(str(scale)+':'+'/'.join(r['status'] if r else 'unrun' for r in pair))
        for moment in moments:
            a = values.get(1.,{}).get(moment,''); b=values.get(.5,{}).get(moment,'')
            rows.append(dict(parameter=name,moment=moment,step=c['jacobian_steps'][name],full_derivative=a,half_derivative=b,difference=b-a if a!='' and b!='' else '',full_curvature=curvatures.get(1.,{}).get(moment,''),half_curvature=curvatures.get(.5,{}).get(moment,''),status='complete' if not statuses else ';'.join(statuses)))
    with (output/'jacobian.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    # Numerical linear algebra executes only in this Torch controller, never Mac.
    import numpy as np
    objective=read(c['lanes']['primary']['objective']['path'])
    scored=[r for r in objective['target_rows'] if r['actual_weight'] is not None and r['actual_weight']>0]
    lookup={(r['parameter'],r['moment']):r for r in rows}
    scales={k:abs(c['initial_point'][k]) for k in FREE}
    report=dict(column_scaling='derivative times absolute anchor parameter; rows multiplied by sqrt(primary weight)',column_scales=scales,moments=[r['restriction_id'] for r in scored],parameters=list(FREE),complete=False)
    if all(lookup[(k,r['restriction_id'])]['status']=='complete' for k in FREE for r in scored):
        reports={}
        for field in ('full_derivative','half_derivative'):
            matrix=np.array([[lookup[(k,r['restriction_id'])][field]*scales[k]*math.sqrt(r['actual_weight']) for k in FREE] for r in scored])
            _, singular, vt=np.linalg.svd(matrix,full_matrices=False)
            tolerance=float(singular[0]*max(matrix.shape)*np.finfo(float).eps)
            reports[field]=dict(matrix=matrix.tolist(),singular_values=singular.tolist(),right_vectors=vt.tolist(),numerical_rank=int((singular>tolerance).sum()),rank_tolerance=tolerance,condition=float(singular[0]/singular[-1]) if singular[-1]>0 else None)
        delta=np.array(reports['half_derivative']['matrix'])-np.array(reports['full_derivative']['matrix'])
        noise=float(np.linalg.norm(delta,2))
        report.update(complete=True,steps=reports,step_difference_spectral_norm=noise,half_singular_values_above_step_difference=int(sum(s>noise for s in reports['half_derivative']['singular_values'])),interpretation='Local numerical sensitivity only; full rank does not establish global identification. Singular values exceeding the step-difference norm are a sensitivity heuristic, not a statistical rank test.')
    write(output/'jacobian_scaled_svd.json',report)


def controller(a,c,objs):
    pin=sha(a.contract); budget=c['budget']; helper=core.module('identification_supervisor',c['files']['recovery_search']['path'])
    if a.stage == 'smoke':
        start=time.time(); timing=dict(start=start,search_cutoff=start+budget['search_seconds'],repeat_cutoff=start+budget['repeat_seconds'],end=start+budget['total_seconds'])
        prior=[]
    else:
        assert sha(a.approval)==a.approval_sha256
        approval=read(a.approval);assert approval['status']=='approved_search' and approval['contract_sha256']==pin
        smoke_pin=approval['smoke_receipt'];assert sha(smoke_pin['path'])==smoke_pin['sha256']
        smoke=read(smoke_pin['path']);assert smoke['status']=='exact_loop_smoke_passed' and smoke['contract_sha256']==pin
        timing=smoke['clock'];prior=smoke['records'];assert len(prior)==6
        for row in prior:
            req=read(row['request_path']); actual=core.validate(Path(row['case_path']).parent,c,objs,req)
            assert all(row[k]==v for k,v in actual.items())
            core.compare_anchor(c['anchor']['case_path'],row['case_path'])
        for lane in LANES[:3]:
            pair=[r for r in prior if r['lane']==lane];assert len(pair)==2
            core.compare_tables(pair[0]['case_path'],pair[1]['case_path'])
        assert timing['end']-timing['start']==budget['total_seconds']
    assert timing['start'] <= time.time() < timing['search_cutoff'], 'No late restart'
    a.output.mkdir(parents=True,exist_ok=False);write(a.output/'clock.json',timing)
    records=[];best={};fatal=False;last=0.;populations={};unrun=[]
    def heartbeat(**extra):
        nonlocal last
        now=time.time()
        if now-last >= 5 or extra.get('force'):
            write(a.output/'heartbeat.json',dict(epoch=now,completed=len(records),stage=a.stage,clock=timing,best={k:v['loss'] for k,v in best.items()},**extra));last=now
    def batch(items,stage,deadline,graphs=False):
        nonlocal fatal
        assert len(prior)+len(records)+len(items) <= budget['max_objective_cases']
        reqs=[]
        for item in items:
            name=f'{stage}_{len(records)+len(reqs):04d}_{item["lane"]}'
            ctx=dict(candidate_id=name,stage=stage,contract_sha256=pin,source_sha256=c['files']['driver']['sha256'],target_sha256=c['lanes'][item['lane']]['objective']['sha256'],point_sha256=canon(item['point']))
            reqs.append(dict(item,id=name,context=ctx))
        def launch(req,batch_deadline):
            end=min(batch_deadline,time.time()+budget['objective_cap_seconds']);folder=a.output/req['id']
            payload=dict(req,contract_sha256=pin,scientific_candidate_id=core.identity(c,req['lane'],req['point']),normalization_inputs=c['normalization'],controller_pid=os.getpid(),deadline_epoch=end,graphs=graphs)
            path=a.output/(req['id']+'.request.json');write(path,payload);req['payload']=payload;req['request_path']=str(path.resolve())
            env=os.environ.copy();env.update({key:'1' for key in core.THREADS})
            return helper.ManagedProcess([sys.executable,__file__,'--stage','evaluate','--contract',str(a.contract),'--output',str(folder),'--request',str(path)],a.output/(req['id']+'.log'),end,env)
        def finish(req,proc,code):
            nonlocal fatal
            try:
                verify(a.contract)
                status,data,error=core.classify(a.output/req['id'],c,objs,req['payload'],proc,code)
            except Exception as exc: status,data,error='fatal',{},str(exc)
            row=dict(case=req['id'],lane=req['lane'],point=req['point'],design=req['design'],request_path=req['request_path'],status=status,error=error,returncode=code,deadline=proc.deadline,**data)
            if status=='success':
                fits=core.table(Path(row['case_path'])/'target_fit.csv','moment');early=fits[c['early_moment']]
                primary=next(r for r in objs['primary']['target_rows'] if r['restriction_id']==c['early_moment'])
                row.update(early_model=float(early['model']),early_gap=float(early['gap']),other_primary_loss=row['primary_rescore']-primary['actual_weight']*float(early['gap'])**2)
                if stage not in ('repeat','jacobian','smoke') and (row['lane'] not in best or row['loss']<best[row['lane']]['loss']):best[row['lane']]=row
            if status=='fatal' or stage in ('smoke','repeat') and status!='success':fatal=True
            row['halt_new_dispatch']=fatal;records.append(row)
            write(a.output/'latest_completed.json',row);write(a.output/'best_so_far.json',best)
            write(a.output/'checkpoint.json',dict(records=records,best=best,populations=populations,clock=timing,contract_sha256=pin));heartbeat(force=True)
            return row
        result=helper.run_batch(reqs,workers=min(budget['workers'],len(reqs)),deadline=deadline,launch=launch,finish=finish,heartbeat=heartbeat,poll_seconds=1,allowed_statuses={'success'} if stage in ('smoke','repeat') else {'success','inadmissible','censored_timeout','censored_late_completion'},guard=lambda:'fatal_stop' if fatal else None)
        finished={r['case'] for r in records};unrun.extend(dict(id=r['id'],lane=r['lane'],point=r['point'],design=r['design']) for r in reqs if r['id'] not in finished)
        write(a.output/'unrun.json',unrun);write(a.output/f'batch_{stage}_{len(records):04d}.json',result)
        if any(r['status']=='failed' for r in result['results']):fatal=True
        return result
    try:
        if a.stage=='smoke':
            batch([dict(lane=lane,point=c['initial_point'],design='anchor_smoke') for lane in LANES[:3] for _ in range(2)],'smoke',timing['search_cutoff'],True)
            assert not fatal and len(records)==6 and all(r['status']=='success' for r in records)
            for row in records:core.compare_anchor(c['anchor']['case_path'],row['case_path'])
            for lane in LANES[:3]:
                pair=[r for r in records if r['lane']==lane]
                core.compare_tables(pair[0]['case_path'],pair[1]['case_path'])
            write(a.output/'complete.json',dict(status='exact_loop_smoke_passed',records=records,clock=timing,contract_sha256=pin));return
        batch(probe_points(c),'jacobian',timing['search_cutoff'])
        jacobian_export(c,records,a.output)
        assert not fatal,'Fatal Jacobian case; no search dispatched'
        # All lanes progress generation by generation, with interleaved slots.
        for generation in range(budget['generations']+1):
            if time.time()>=timing['search_cutoff']:break
            proposals={lane:initial_population(c,objs[lane],lane) if generation==0 else de_trials(c,objs[lane],lane,populations[lane],generation) for lane in LANES}
            items=[dict(lane=lane,point=proposals[lane][slot],design=f'de:{generation}:{slot}') for slot in range(budget['population']) for lane in LANES]
            before=len(records);result=batch(items,'search',timing['search_cutoff'])
            new=records[before:];lookup={(r['lane'],r['design']):r for r in new}
            if fatal or not result['complete']:break
            for lane in LANES:
                trials=[lookup[(lane,f'de:{generation}:{slot}')] for slot in range(budget['population'])]
                populations[lane]=trials if generation==0 else [choose(p,t) for p,t in zip(populations[lane],trials)]
            write(a.output/f'generation_{generation:02d}.json',dict(populations=populations,best=best))
        assert not fatal,'Fatal search case; final selection withheld'
        attempt_counts={lane:sum(r['lane']==lane and r['case'].startswith('search_') for r in records) for lane in LANES}
        write(a.output/'search_budget_accounting.json',dict(planned_per_lane=budget['population']*(budget['generations']+1),attempted_per_lane=attempt_counts,unrun_per_lane={lane:budget['population']*(budget['generations']+1)-n for lane,n in attempt_counts.items()},all_attempts_including_smokes=len(prior)+len(records)))
        selected=dict(best);write(a.output/'selected.json',dict(selected=selected,unavailable_lanes=[k for k in LANES if k not in selected],frozen_before_repeats=True))
        assert selected,'No admissible search candidate'
        before=len(records);batch([dict(lane=lane,point=row['point'],design='repeat:'+lane) for lane,row in selected.items() for _ in range(2)],'repeat',timing['repeat_cutoff'],True)
        repeats=records[before:];assert not fatal and len(repeats)==2*len(selected)
        for lane,row in selected.items():
            pair=[r for r in repeats if r['lane']==lane];assert len(pair)==2
            for repeat in pair:core.compare_tables(row['case_path'],repeat['case_path']);assert repeat['loss']==row['loss']
            for name in c['standard_diagnostic_names']:
                assert sha(Path(pair[0]['case_path'])/'standard_diagnostics'/name)==sha(Path(pair[1]['case_path'])/'standard_diagnostics'/name),'Repeat plots differ'
            source=Path(pair[0]['case_path']);out=a.output/'selected_export'/lane;out.mkdir(parents=True)
            for file in source.iterdir():
                assert time.time()<timing['end']
                if file.is_file() and file.suffix in ('.json','.csv'):shutil.copy2(file,out/file.name)
            shutil.copytree(source/'standard_diagnostics',out/'standard_diagnostics')
            (out/'initial_state.pkl.gz').symlink_to((source/'initial_state.pkl.gz').resolve())
            write(out/'export_receipt.json',dict(selected=row,repeats=pair,contract_sha256=pin,diagnostic_only=True))
            write(out/'parameter_roles.json',dict(searched=[k for k in FREE if k not in c['lanes'][lane]['fixed']],fixed_for_profile=c['lanes'][lane]['fixed'],normalized=['psi_child'],original_31_parameter_table_statuses_preserved=True))
            write(out/'artifact_hashes.json',{str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='initial_state.pkl.gz'})
        assert time.time()<timing['end']
        write(a.output/'complete.json',dict(status='bounded_experiment_complete',records=records,selected=selected,unavailable_lanes=[k for k in LANES if k not in selected],clock=timing,contract_sha256=pin,unrun=unrun))
    except Exception as exc:
        write(a.output/'complete.json',dict(status='incomplete_or_fatal_stop',error_type=type(exc).__name__,error=str(exc),records=records,best=best,clock=timing,contract_sha256=pin,unrun=unrun));raise


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=('prepare','smoke','run','evaluate'),required=True)
    for name in ('contract','output'):parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--request',type=Path);parser.add_argument('--approval',type=Path);parser.add_argument('--approval-sha256')
    a=parser.parse_args();assert os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm only'
    c,objs=verify(a.contract)
    if a.stage=='prepare':
        a.output.mkdir(parents=True,exist_ok=False);write(a.output/'prepared.json',dict(status='verified_zero_solves',contract_sha256=sha(a.contract)));return
    if a.stage=='evaluate':core.evaluate(a,c,objs)
    else:controller(a,c,objs)

if __name__=='__main__':main()
