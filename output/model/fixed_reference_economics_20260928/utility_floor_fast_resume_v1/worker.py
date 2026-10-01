"""One unchanged native full GE from a controller-authenticated task."""
import argparse, importlib.util, json, os, signal, sys, time
from pathlib import Path

def load_runtime(path):
    sys.path.insert(0,str(path))
    spec=importlib.util.spec_from_file_location('original_utility_runner',path/'runner.py')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def main():
    p=argparse.ArgumentParser();p.add_argument('--runtime',type=Path,required=True);p.add_argument('--task-file',type=Path,required=True);p.add_argument('--index',type=int,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();task=json.loads(a.task_file.read_text())['tasks'][a.index]
    m=load_runtime(a.runtime);m.verify_sources()
    plan=json.loads(Path(__file__).with_name('plan.json').read_text())
    for n,h in plan['runtime_files'].items():m.require(m.inputs.sha(a.runtime/n)==h,'Runtime authentication drift '+n)
    m.require(m.source_fingerprint()==plan['source_fingerprint'],'Source fingerprint drift')
    m.require(m.inputs.canonical(m.PLAN['target_contract'])==plan['target_fingerprint'],'Target fingerprint drift')
    m.require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit() and os.environ.get('SLURM_CPUS_PER_TASK','1')=='1','Native run requires one-core Torch')
    deadline=float(task['deadline_epoch']);m.require(deadline>time.time(),'Global budget expired')
    a.out.mkdir(parents=True,exist_ok=False)
    signal.signal(signal.SIGALRM,m.alarm);signal.setitimer(signal.ITIMER_REAL,deadline-time.time())
    try:
        P,grid=m.inputs.proposal('floor_s2');P,entry=m.inputs.entry(P,grid,'nonnegative_mean')
        Q=m.utility_checks(P,grid,'floor_s2',a.out)
        evaluate=m.native_evaluator(a.out,'floor_s2',Q,grid,deadline,float(task['starting_price']))
        result=evaluate(task['label'],task['parameters'],deadline)
        row=dict(label=task['label'],kind=task['kind'],parameters=task['parameters'],**result)
        if row['status']=='passed':row['loss']=float(sum(x*x for x in row['residual']))
        m.write(a.out/'case.json',row)
        m.require(row['status'] in ('passed','budget_exhausted','inadmissible_numerical'),'Native GE critical failure: '+row['status'])
    except BaseException as e:
        m.write(a.out/'failure.json',dict(type=type(e).__name__,message=str(e),no_auto_retry=True));raise
    finally:signal.setitimer(signal.ITIMER_REAL,0)
if __name__=='__main__':main()
