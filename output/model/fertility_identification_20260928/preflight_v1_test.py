#!/usr/bin/env python3
"""Torch-only synthetic controller tests. No household solves or imports."""
import argparse, copy, os, signal, sys, tempfile, time, types, unittest
from pathlib import Path
from unittest import mock
if not os.environ.get('SLURM_JOB_ID','').isdigit():raise RuntimeError('Torch Slurm only')
import run_e5f_fertility_identification as d
fx=d.core.module('identification_fixtures',Path(__file__).with_name('test_e5f_evening_calibration.py'))


def worker():
    p=argparse.ArgumentParser();p.add_argument('--worker',action='store_true');p.add_argument('--stage');p.add_argument('--contract',type=Path);p.add_argument('--output',type=Path);p.add_argument('--request',type=Path);a=p.parse_args()
    c=d.read(a.contract);objs={k:d.read(v['objective']['path']) for k,v in c['lanes'].items()};req=d.read(a.request)
    assert req['controller_pid']==os.getppid()
    fx.success(a.output,c,objs,req,pid=os.getpid())
    start=d.read(a.output/'startup.json');start['parent_pid']=os.getppid();d.write(a.output/'startup.json',start)


class Tests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
        c,objs,self.path=fx.fixture(self.root);base=objs['primary'];c['lanes']={};self.objs={}
        for lane in d.LANES:
            obj=copy.deepcopy(base);pin=self.root/(lane+'.json');d.write(pin,obj)
            self.objs[lane]=obj;c['lanes'][lane]=dict(objective=dict(path=str(pin),sha256=d.sha(pin)),canonical_sha256=d.canon(obj),target_weight_fingerprint=d.canon(obj['target_rows']),fixed={} if not lane.startswith('profile') else {'kappa_fert_continuation':.25})
        c.update(early_moment='m0',jacobian_steps=dict.fromkeys(d.FREE,.002),seed_points=[dict(c['initial_point']),dict(c['initial_point'],chi=.6)])
        c['budget'].update(population=4,generations=1,max_objective_cases=538,search_seconds=18000,repeat_seconds=21000,objective_cap_seconds=1800)
        c['anchor']={'case_path':str(self.root/'anchor/case'),'receipt':{'path':str(self.root/'anchor/case/receipt.json')}}
        self.c=c;d.write(self.path,c)
        fx.success(self.root/'anchor',c,self.objs,fx.req(c,self.path,self.root/'anchor'))
    def test_forty_symmetric_points(self):
        points=d.probe_points(self.c);self.assertEqual(len(points),40)
        for key in d.FREE:
            rows=[r for r in points if r['coordinate']==key];self.assertEqual(len(rows),4)
            for row in rows:self.assertAlmostEqual(row['point'][key]-self.c['initial_point'][key],row['scale']*row['sign']*.002)
    def test_multistart_and_fixed_profile(self):
        for lane in d.LANES:
            points=d.initial_population(self.c,self.objs[lane],lane)
            self.assertEqual(points[0]['chi'],.5);self.assertEqual(points[1]['chi'],.6)
            pop=[dict(point=p,status='inadmissible') for p in points]
            trial=d.de_trials(self.c,self.objs[lane],lane,pop,1)
            self.assertEqual(trial,d.de_trials(self.c,self.objs[lane],lane,pop,1))
            for point in points+trial:
                for k,v in self.c['lanes'][lane]['fixed'].items():self.assertEqual(point[k],v)
                self.assertTrue(.001 <= point['tenure_choice_kappa'] <= .1)
    def test_selection_replaces_inadmissible_and_never_worsens(self):
        self.assertEqual(d.choose({'status':'inadmissible'},{'status':'success','loss':100})['loss'],100)
        self.assertEqual(d.choose({'status':'success','loss':1},{'status':'success','loss':2})['loss'],1)
    def test_unknown_kill_not_timeout(self):
        folder=self.root/'unknown';folder.mkdir();req=fx.req(self.c,self.path,folder)
        proc=types.SimpleNamespace(process=types.SimpleNamespace(pid=10),observed_running_at_expiry=False,deadline_kill_reaped=False,deadline_expired=True)
        self.assertEqual(d.core.classify(folder,self.c,self.objs,req,proc,-signal.SIGKILL)[0],'fatal')
    def test_owned_real_process_deadline(self):
        def launch(req,deadline):return fx.supervisor.ManagedProcess([sys.executable,'-c','import time; time.sleep(60)'],self.root/'timeout.log',min(deadline,time.time()+.1),os.environ.copy())
        seen=[]
        def finish(req,proc,code):
            seen.append((proc.observed_running_at_expiry,proc.deadline_kill_reaped,code));return {'status':'censored_timeout'}
        fx.supervisor.run_batch([{}],workers=1,deadline=time.time()+3,launch=launch,finish=finish,heartbeat=lambda **kw:None,poll_seconds=.05,allowed_statuses={'censored_timeout'})
        self.assertEqual(seen,[(True,True,-signal.SIGKILL)])
    def test_actual_controller_process_loop(self):
        started=[]
        def managed(command,log,deadline,env):
            self.assertTrue(all(env[k]=='1' for k in d.core.THREADS))
            rewritten=[sys.executable,str(Path(__file__).resolve()),'--worker']+command[2:]
            proc=fx.supervisor.ManagedProcess(rewritten,log,deadline,env);started.append(proc.process.pid);return proc
        proxy=types.SimpleNamespace(ManagedProcess=managed,run_batch=fx.supervisor.run_batch)
        a=argparse.Namespace(stage='smoke',contract=self.path,output=self.root/'smoke')
        with mock.patch.object(d.core,'module',return_value=proxy),mock.patch.object(d,'verify',return_value=(self.c,self.objs)):
            d.controller(a,self.c,self.objs)
            smoke=a.output/'complete.json';self.assertEqual(d.read(smoke)['status'],'exact_loop_smoke_passed')
            approval=self.root/'approval.json';d.write(approval,dict(status='approved_search',contract_sha256=d.sha(self.path),smoke_receipt=dict(path=str(smoke),sha256=d.sha(smoke))))
            b=argparse.Namespace(stage='run',contract=self.path,output=self.root/'run',approval=approval,approval_sha256=d.sha(approval))
            try:d.controller(b,self.c,self.objs)
            except Exception:
                import shutil
                destination=Path(os.environ.get('IDENTIFICATION_TEST_FAILURE_DIR','/tmp'))/('identification_fixture_'+str(os.getpid()))
                shutil.copytree(self.root,destination);print('PRESERVED FAILED FIXTURE',destination);raise
        done=d.read(b.output/'complete.json');self.assertEqual(done['status'],'bounded_experiment_complete')
        self.assertEqual(done['clock'],d.read(smoke)['clock']);self.assertEqual(len(done['selected']),6)
        self.assertEqual(len(started),6+40+48+12)
        self.assertTrue(d.read(b.output/'jacobian_scaled_svd.json')['complete'])
        for lane in d.LANES:
            exported=b.output/'selected_export'/lane
            self.assertEqual(len(list((exported/'standard_diagnostics').glob('*.png'))),17)
            self.assertEqual(len(d.read(exported/'export_receipt.json')['repeats']),2)
    def test_anchor_mutation_rejected(self):
        original=self.root/'anchor/case';other=self.root/'other'
        fx.success(other,self.c,self.objs,fx.req(self.c,self.path,other))
        rows=list(d.core.table(other/'case/target_fit.csv','moment').values());rows[12]['model']='9';fx.csvwrite(other/'case/target_fit.csv',rows)
        with self.assertRaises(AssertionError):d.core.compare_anchor(original,other/'case')

if __name__=='__main__':
    if '--worker' in sys.argv:worker()
    else:unittest.main()
