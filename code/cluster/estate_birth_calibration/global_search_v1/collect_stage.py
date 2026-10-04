"""Collect stage-1 receipts and propose diverse points; never certify a winner."""
import argparse,csv,hashlib,json,math
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')
def normalized(point,plan):
    row=[]
    for k in sorted(plan['bounds']):
        lo,hi=plan['bounds'][k];v=point[k]
        if plan['transforms'][k]=='log':lo,hi,v=math.log(lo),math.log(hi),math.log(v)
        row.append((v-lo)/(hi-lo))
    return row
def distance(a,b):return math.sqrt(sum((x-y)**2 for x,y in zip(a,b)))
def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    plan=json.loads((a.root/'control/plan.json').read_text());plan_sha=sha(a.root/'control/plan.json')
    rows=[];tasks=[]
    for task in range(16):
        run=a.root/'results'/f'production_task_{task}'
        if not run.exists():tasks.append(dict(task=task,status='not_started'));continue
        term=run/'launcher_terminal.json';done=run/'run/completed.json';cases=run/'run/cases.json'
        t=json.loads(term.read_text()) if term.exists() else None
        d=json.loads(done.read_text()) if done.exists() else None
        cc=json.loads(cases.read_text()) if cases.exists() else []
        if d:assert d['plan_sha256']==plan_sha and d['target_fingerprint']==plan['target_fingerprint'] and d['weight_fingerprint']==plan['weight_fingerprint']
        for row in cc:
            ix=row['index'];assert 4*task<=ix<4*task+4 and row['parameters']==plan['points'][ix]
            if row['status']=='passed':
                assert row['valid_loss'] and math.isfinite(row['loss']) and len(row['case_result']['target_fit'])==14 and len(row['case_result']['parameter_table'])==31
                assert row['case_result']['target_fingerprint']==plan['target_fingerprint'] and row['case_result']['weight_fingerprint']==plan['weight_fingerprint']
            else:assert row['valid_loss'] is False and 'loss' not in row and row['error_type']
            rows.append(dict(task=task,**row))
        tasks.append(dict(task=task,status=d['status'] if d else 'running_or_failed',terminal_exit_code=t['exit_code'] if t else None,cases=len(cc),
                          last_heartbeat=json.loads((run/'run/heartbeat.json').read_text()) if (run/'run/heartbeat.json').exists() else None))
    assert len({r['index'] for r in rows})==len(rows),'Duplicate global point'
    feasible=sorted((r for r in rows if r['status']=='passed'),key=lambda r:r['loss'])
    candidate=[];coords=[]
    for row in feasible:
        unit=normalized(row['parameters'],plan)
        if all(distance(unit,x)>=0.25 for x in coords):
            candidate.append(dict(source='stage1_sobol',index=row['index'],task=row['task'],provisional_loss=row['loss'],parameters=row['parameters']))
            coords.append(unit)
            if len(candidate)==8:break
    all_terminal=all(t['terminal_exit_code'] is not None for t in tasks)
    summary=dict(status=('terminal_with_failures' if any(t['terminal_exit_code']!=0 for t in tasks) else 'terminal') if all_terminal else 'in_progress',
        plan_sha256=plan_sha,completed_points=len(rows),expected_points=64,passed_points=len(feasible),
        rejected_or_timed_out=len(rows)-len(feasible),tasks=tasks,
        best_provisional=dict(index=feasible[0]['index'],task=feasible[0]['task'],loss=feasible[0]['loss']) if feasible else None,
        seed_proposal=candidate,seed_distance_floor_transformed=0.25,
        verified_incumbent=plan['incumbent_control'],
        warning='All stage-1 losses provisional until fresh native selected-point and exact-repeat checks; Stage 2 requires lead seed review.')
    a.out.mkdir(parents=True,exist_ok=True);write(a.out/'stage1_summary.json',summary)
    if feasible:
        best=feasible[0]['case_result']
        for name,key in [('best_provisional_fit.csv','target_fit'),('best_provisional_parameters.csv','parameter_table')]:
            rows2=best[key]
            with (a.out/name).open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(rows2[0]));w.writeheader();w.writerows(rows2)
    print(json.dumps(dict(status=summary['status'],completed_points=len(rows),passed_points=len(feasible),best_provisional=summary['best_provisional'])))
if __name__=='__main__':main()
