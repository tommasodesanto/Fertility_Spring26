#!/usr/bin/env python3
"""Authenticated, bounded continuation of the unfinished shares_020 match only."""
import argparse,csv,hashlib,json,os,shutil,time
from pathlib import Path
import numpy as np
import run_e5f_housing_fertility_cost_diagnostic as d
import run_e5f_first_child_loading_probe as loading
import run_e5f_soft_housing_probe as common

MAX_TOTAL=16
MAX_NEW=8
SOLVE_SECONDS=600

def fingerprint(root):
    return {str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in root.rglob('*') if p.is_file() and p.suffix in ('.json','.csv')}

def bracket(rows):
    low=max((r for r in rows if float(r['fertility'])<d.TARGET),key=lambda r:float(r['b']))
    high=min((r for r in rows if float(r['fertility'])>d.TARGET),key=lambda r:float(r['b']))
    assert float(low['b'])<float(high['b'])
    return float(low['b']),float(high['b'])

def copy_case(source,dest):
    dest.mkdir(parents=True,exist_ok=False)
    for p in source.iterdir():
        if p.name in ('state.pkl.gz','standard_diagnostics'): (dest/p.name).symlink_to(p.resolve(),target_is_directory=p.is_dir())
        elif p.is_file() and p.suffix in ('.json','.csv'):shutil.copy2(p,dest/p.name)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['smoke','run'],required=True)
    for name in ('source','output','combined'):ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    assert os.environ.get('SLURM_JOB_ID'),'Torch Slurm execution required'
    assert a.source.resolve()!=a.output.resolve() and a.source.resolve()!=a.combined.resolve()
    design=json.loads((a.source/'design.json').read_text());old=json.loads((a.source/'complete.json').read_text())
    history=list(csv.DictReader((a.source/'solve_history.csv').open()))
    prior=[r for r in history if r['case'].startswith('matchsearch_shares_020_')]
    assert len(prior)==9 and not any(r['phase']=='matched' and r['arm']=='shares_020' for r in old['records'])
    assert len(old['records'])==10 and design['contract']==common.PIN
    assert design['bracket']==list(d.B_BRACKET) and design['target']==d.TARGET and design['tolerance']==d.TOL
    before=fingerprint(a.source)
    rt,tax,obj,selected=d.setup(a.output/('context_'+a.stage+'_'+os.environ['SLURM_JOB_ID']))
    assert selected['pins']==design['source_pins']
    control=selected['packet'];P=control['parameters'];grid=control['b_grid'];price=np.asarray(control['solution'].p_eq)
    assert P.utility_reference_rent==design['frozen_reference_rent'] and P.pension==design['fixed_pension'] and P.tau_pay==design['fixed_payroll_tax']
    np.testing.assert_array_equal(price,design['baseline_price'])
    lo,hi=bracket(prior);mid=(lo+hi)/2
    assert all(abs(mid-float(r['b']))>1e-14 for r in prior)
    # Exercise exact bounded-loop indexing, bracket updates, and checkpoint labels without a solve.
    trial=list(prior);synthetic=[]
    for index in range(len(prior),MAX_TOTAL):
        l,h=bracket(trial);b=(l+h)/2
        assert not any(abs(b-float(r['b']))<1e-14 for r in trial)
        synthetic.append({'case':f'matchsearch_shares_020_{index:02d}','b':b,'fertility':d.TARGET-.01})
        trial.append(synthetic[-1])
    assert len(synthetic)+1<=MAX_NEW and synthetic[0]['b']==mid
    Q=loading.configure(P,'shares_020',.2);Q.psi_child=mid
    assert Q.pension==P.pension and Q.utility_reference_rent==P.utility_reference_rent
    plan=dict(stage=a.stage,source=str(a.source),output=str(a.output),combined=str(a.combined),source_pins=selected['pins'],initial_bracket=[lo,hi],first_midpoint=mid,prior_evaluations=len(prior),max_total_b=MAX_TOTAL,max_new_native=MAX_NEW,solve_seconds=SOLVE_SECONDS,synthetic_loop=synthetic,driver_sha256=hashlib.sha256(Path(d.__file__).read_bytes()).hexdigest(),helper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    common.write(a.output/('smoke.json' if a.stage=='smoke' else 'continuation_manifest.json'),plan)
    if a.stage=='smoke':
        assert before==fingerprint(a.source);print(json.dumps(plan));return
    assert (a.output/'smoke.json').exists() and not (a.output/'solve_history.csv').exists()
    smoke=json.loads((a.output/'smoke.json').read_text())
    for key in ('source_pins','helper_sha256','driver_sha256','initial_bracket'):assert smoke[key]==plan[key]
    start=time.monotonic();new=[];records=[];matched=None;stop='matching_budget_exhausted'
    def evaluate(b,factor,label):
        assert len(new)<MAX_NEW
        if time.monotonic()-start>SOLVE_SECONDS:raise TimeoutError('Numerical reserve reached')
        assert not any(abs(float(r['b'])-b)<1e-14 and float(r['factor'])==factor and r['arm']=='shares_020' for r in history+new)
        Q=loading.configure(P,'shares_020',.2);Q.psi_child=b
        common.write(a.output/'heartbeat.json',dict(case=label,status='solving',elapsed=time.monotonic()-start,job=os.environ.get('SLURM_JOB_ID')))
        packet=d.solve(rt,Q,grid,price*factor)
        folder=a.output/label;folder.mkdir();d.save(folder/'state.pkl.gz',packet)
        tfr=float(rt['chain'].extract_moments(packet['solution'],Q)['tfr'])
        row=dict(case=label,arm='shares_020',b=b,factor=factor,fertility=tfr,elapsed_seconds=time.monotonic()-start)
        new.append(row);common.table(a.output/'solve_history.csv',new);common.write(a.output/'latest_completed.json',row)
        common.write(a.output/'best_so_far.json',min([r for r in new if r['factor']==1],key=lambda r:abs(r['fertility']-d.TARGET)))
        print(json.dumps(row),flush=True);return packet,tfr
    try:
        for index in range(len(prior),MAX_TOTAL):
            lo,hi=bracket(prior+new);b=(lo+hi)/2
            packet,tfr=evaluate(b,1.,f'matchsearch_shares_020_{index:02d}')
            if abs(tfr-d.TARGET)<=d.TOL:matched=(b,packet,tfr);break
        if matched:
            b,packet,tfr=matched
            folder=a.output/'matched_shares_020_1.0';folder.mkdir()
            records.append(d.report(packet,rt,tax,obj,selected['case'],folder,'shares_020',.2,'matched',1.,control))
            common.write(a.output/'comparison_summary.json',records)
            shocked,_=evaluate(b,1.1,'matched_shares_020_1.1')
            records.append(d.report(shocked,rt,tax,obj,selected['case'],a.output/'matched_shares_020_1.1','shares_020',.2,'matched',1.1,control))
            stop='completed'
    except TimeoutError as exc:stop=str(exc)
    assert before==fingerprint(a.source),'Source compact outputs changed'
    a.combined.mkdir(parents=True,exist_ok=False)
    allrecords=old['records']+records
    for rec in allrecords:copy_case((a.source if rec in old['records'] else a.output)/rec['case'],a.combined/rec['case'])
    for name in ('design.json','control_replay.json'):shutil.copy2(a.source/name,a.combined/name)
    common.table(a.combined/'solve_history.csv',history+new)
    matches=old['matching_status'].copy()
    matches['shares_020']=dict(status='matched' if len(records)==2 else stop,search_evaluations=len(prior)+len([r for r in new if r['factor']==1]),b=matched[0] if matched else None,fertility=matched[2] if matched else None)
    complete=dict(status=stop,records=allrecords,matching_status=matches,evaluations=old['evaluations']+len(new),elapsed_seconds=old['elapsed_seconds']+time.monotonic()-start,job=os.environ.get('SLURM_JOB_ID'),source_job=old['job'],max_b_evaluations=MAX_TOTAL,continuation_source=str(a.output),original_source=str(a.source),numerical_budget_extension_only=len(prior)+len([r for r in new if r['factor']==1])>12)
    common.write(a.combined/'complete.json',complete);common.write(a.combined/'comparison_summary.json',allrecords)
    common.write(a.output/'completion.json',dict(status=stop,new_evaluations=len(new),source_unchanged=True,combined=str(a.combined),elapsed=time.monotonic()-start))
    print(json.dumps(dict(status=stop,new_evaluations=len(new),reported_cases=len(allrecords),job=os.environ.get('SLURM_JOB_ID'))),flush=True)

if __name__=='__main__':main()
