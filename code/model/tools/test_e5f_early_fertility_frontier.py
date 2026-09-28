#!/usr/bin/env python3
"""Torch-only synthetic frontier tests, including real owned subprocesses."""
import argparse,copy,importlib.util,json,os,sys,tempfile,time,types,unittest
from pathlib import Path
from unittest import mock
if not os.environ.get('SLURM_JOB_ID','').isdigit():raise RuntimeError('Torch Slurm only')
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('frontier_test_subject',HERE/'run_e5f_early_fertility_frontier.py');f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f)
# PYTHONPATH points to the immutable original tools directory on Torch.
import test_e5f_evening_calibration as fx
import run_e5f_night_calibration as night

def fixture(root):
    c,objs,contract=fx.fixture(root);c['initial_point'].update(kappa_fert=.1,kappa_fert_continuation=.1)
    for lane,obj in objs.items():
        obj['target_rows'][0]['restriction_id']='early_fertility';pin=c['lanes'][lane];f.write(pin['objective']['path'],obj);pin['objective']['sha256']=f.sha(pin['objective']['path']);pin['canonical_sha256']=night.canon(obj);pin['target_weight_fingerprint']=night.canon(obj['target_rows'])
    f.write(contract,c);anchors={}
    for name,h0 in (('primary_best',.5),('early_best',.6)):
        folder=root/name;point=dict(c['initial_point'],H0=h0);req=fx.req(c,contract,folder,'primary',point)
        data=fx.success(folder,c,objs,req);request=root/(name+'.request.json');f.write(request,req)
        row=dict(case=name,point=point,lane='primary',status='success',request_path=str(request),**data)
        early=night.table(folder/'case/target_fit.csv','moment')['early_fertility'];anchors[name]=dict(record=row,primary_loss=data['primary_rescore'],early_gap=float(early['gap']),table_pins={key:f.sha(folder/'case'/key) for key in ('target_fit.csv','parameters.csv','receipt.json')})
    return c,objs,contract,anchors

def worker():
    p=argparse.ArgumentParser();p.add_argument('--synthetic-worker',action='store_true');p.add_argument('--stage');p.add_argument('--contract',type=Path);p.add_argument('--output',type=Path);p.add_argument('--request',type=Path);a=p.parse_args()
    c=f.read(a.contract);objs={lane:f.read(c['lanes'][lane]['objective']['path']) for lane in night.LANES};req=f.read(a.request)
    assert req['controller_pid']==os.getppid() and req['graphs'] is True
    fx.success(a.output,c,objs,req,pid=os.getpid());startup=f.read(a.output/'startup.json');startup['parent_pid']=os.getppid();f.write(a.output/'startup.json',startup)
    if os.environ.get('FRONTIER_TEST_BAD_ANCHOR')=='1' and req['id']=='anchor_primary_best':
        case=a.output/'case';rows=list(night.table(case/'target_fit.csv','moment').values());rows[0].update(model=3.,gap=2.,loss_contribution=float(rows[0]['weight'])*4);fx.csvwrite(case/'target_fit.csv',rows)
        receipt=f.read(case/'receipt.json');receipt['loss']=sum(float(r['loss_contribution']) for r in rows if r['loss_contribution']!='');f.write(case/'receipt.json',receipt)
        side=f.read(a.output/'success.json');side.update(loss=receipt['loss'],receipt_sha256=f.sha(case/'receipt.json'));f.write(a.output/'success.json',side)
    f.write(a.output/'synthetic_execution.json',dict(pid=os.getpid(),parent_pid=os.getppid(),model_solves=0))

class FrontierTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name);self.c,self.objs,self.contract,self.anchors=fixture(self.root)
    def test_recipes_preserve_unvaried_coordinates(self):
        rows=f.requests_for(self.anchors,self.c,self.objs['primary']);self.assertEqual(len(rows),10);self.assertFalse(any(r.get('duplicate_of') for r in rows))
        for row in rows[2:]:
            for key in ('H0','beta_annual','chi','theta0','child_benefit_curvature','tenure_choice_kappa'):self.assertEqual(row['point'][key],self.anchors['primary_best']['record']['point'][key])
        self.assertEqual(rows[-1]['point']['first_birth_fixed_cost'],0);self.assertEqual(rows[-1]['point']['delta_alpha_jump'],0)
    def test_clipping_and_duplicate_without_replacement(self):
        self.anchors['primary_best']['record']['point'].update(kappa_fert=0.,kappa_fert_continuation=1.,first_birth_fixed_cost=0.,delta_alpha_jump=0.)
        rows=f.requests_for(self.anchors,self.c,self.objs['primary']);self.assertEqual(len(rows),10)
        self.assertTrue(rows[2]['parameter_changes']['kappa_fert_continuation']['clipped']);self.assertTrue(rows[2]['duplicate_of']);self.assertTrue(all(r.get('duplicate_of') for r in rows[2:]))
    def test_plan_fingerprint_refuses_mutation_before_imports(self):
        path=self.root/'plan.json';f.write(path,{'schema':f.SCHEMA});pin=f.sha(path);f.write(path,{'schema':'changed'})
        with self.assertRaises(AssertionError):f.verify_plan(path,pin)
    def test_fixed_window_and_fresh_ended_proof(self):
        p=dict(not_before=1790593200,hard_end=1790595900,objective_end=1790595600,total_seconds=2400)
        t=f.timing(p,1790593200);self.assertEqual(t['end'],1790595600);self.assertEqual(t['objective_end'],1790595300)
        with self.assertRaises(AssertionError):f.timing(p,1790593199)
    def test_nondominance_uses_other_moment_loss(self):
        rows=[dict(id='a',early_gap=.2,primary_loss=100.,other_moments_primary_loss=1.),dict(id='b',early_gap=.3,primary_loss=2.,other_moments_primary_loss=2.)]
        self.assertEqual(f.nondominated(rows),['a'])
    def test_resource_gate_requires_fresh_at_most8(self):
        now=1790593200.;p=dict(not_before=now,main_heartbeat=str(self.root/'main_heartbeat.json'));clock=dict(latest_dispatch=now+1)
        f.write(p['main_heartbeat'],dict(epoch=now,active=8))
        with mock.patch.object(f.time,'time',return_value=now):self.assertEqual(f.wait_resources(p,clock,self.root)['active'],8)
        f.write(p['main_heartbeat'],dict(epoch=now,active=24))
        with mock.patch.object(f.time,'time',side_effect=[now,now,now+2]),mock.patch.object(f.time,'sleep'):
            with self.assertRaises(TimeoutError):f.wait_resources(p,clock,self.root)
    def test_stale_heartbeat_needs_exact_fresh_completion_proof(self):
        now=1790593200.;complete=self.root/'main_complete.json';proof=self.root/'proof.json'
        p=dict(not_before=now,main_job_id='123',main_complete=str(complete),main_heartbeat=str(self.root/'main_heartbeat.json'));clock=dict(latest_dispatch=now+1)
        f.write(p['main_heartbeat'],dict(epoch=now-31,active=0));f.write(complete,dict(status='bounded_search_complete',contract_sha256=f.CONTRACT_SHA))
        evidence=dict(job_id='123',state='COMPLETED',exit_code='0:0',verified_epoch=now,complete_sha256=f.sha(complete));f.write(proof,evidence)
        with mock.patch.object(f.time,'time',return_value=now):self.assertEqual(f.wait_resources(p,clock,self.root,proof,f.sha(proof))['method'],'wrapper_verified_ended_main')
        with self.assertRaises(AssertionError):f.completion_proof(p,proof,'badpin',now)
        with self.assertRaises(AssertionError):f.completion_proof(p,proof,f.sha(proof),now+61)
    def test_compact_scan_does_not_read_checkpoint(self):
        row=self.anchors['primary_best']['record'];source=Path(row['case_path']);(source/'initial_state.pkl.gz').unlink()
        with mock.patch.object(f,'CONTRACT_SHA',f.sha(self.contract)):
            audit=f.compact_case(row,self.c,self.objs,night)
        self.assertEqual(audit['primary_loss'],row['primary_rescore'])
    def test_compact_loss_allows_only_float_roundoff(self):
        row=dict(self.anchors['primary_best']['record']);row['loss']+=1e-14
        with mock.patch.object(f,'CONTRACT_SHA',f.sha(self.contract)):
            f.compact_case(row,self.c,self.objs,night)
            row['loss']+=.001
            with self.assertRaises(AssertionError):f.compact_case(row,self.c,self.objs,night)
    def real_batch(self,bad_anchor):
        output=self.root/('bad_batch' if bad_anchor else 'good_batch');output.mkdir();p=dict(night_driver=dict(path=night.__file__,sha256=f.sha(night.__file__)),contract=dict(path=str(self.contract),sha256=f.sha(self.contract)))
        launched=[]
        def managed(command,log,deadline,env):
            self.assertTrue(all(env[k]=='1' for k in f.THREADS));rewritten=[sys.executable,str(Path(__file__).resolve()),'--synthetic-worker']+command[2:]
            proc=fx.supervisor.ManagedProcess(rewritten,log,deadline,env);launched.append(proc.process.pid);return proc
        proxy=types.SimpleNamespace(ManagedProcess=managed,run_batch=fx.supervisor.run_batch);now=time.time();clock=dict(start=now,end=now+120,objective_end=now+100,latest_dispatch=now+60)
        with mock.patch.object(f,'CONTRACT_SHA',f.sha(self.contract)),mock.patch.object(night,'verify',return_value=(self.c,self.objs)),mock.patch.dict(os.environ,FRONTIER_TEST_BAD_ANCHOR='1' if bad_anchor else '0'):
            result=f.run(p,self.c,self.objs,night,proxy,output,clock,self.anchors)
        self.assertEqual(len(launched),10);self.assertEqual(len(set(launched)),10);self.assertEqual(len(result['records']),10)
        self.assertTrue(all(r['status']=='success' for r in result['records']),result['records'])
        self.assertEqual(len(list(output.glob('*/case/standard_diagnostics/*.png'))),170)
        self.assertEqual(len(list(__import__('csv').DictReader((output/'full_target_comparison.csv').open()))),14*12)
        self.assertEqual(len(list(__import__('csv').DictReader((output/'full_parameter_comparison.csv').open()))),31*12)
        self.assertEqual(result['accepted'],not bad_anchor)
        self.assertEqual(bool(result['summary']['nondominated_ids']),not bad_anchor)
        print(json.dumps(dict(test='real_frontier_batch',owned_processes=10,anchor_mismatch_injected=bad_anchor,accepted=result['accepted'],model_solves=0)))
    def test_real_owned_batch_all10_and_complete_reports(self):self.real_batch(False)
    def test_real_owned_batch_anchor_failure_rejects_all_claims(self):self.real_batch(True)

if __name__=='__main__':
    if '--synthetic-worker' in sys.argv:worker()
    else:unittest.main()
