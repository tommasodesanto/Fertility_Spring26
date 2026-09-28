#!/usr/bin/env python3
"""Draft synthetic night-controller checks. Torch only; no household solves.

These supplement, rather than replace, the inherited exact-loop/supervisor tests.
Six fresh numerical smokes and end-to-end saved rendering remain launch gates.
"""
import argparse,copy,gzip,importlib.util,json,os,pickle,shutil,sys,tempfile,time,types,unittest
from unittest import mock
from pathlib import Path
if not os.environ.get('SLURM_JOB_ID','').isdigit():raise RuntimeError('Torch Slurm only')
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('night_subject',HERE/'run_e5f_night_calibration.py');d=importlib.util.module_from_spec(s);s.loader.exec_module(d)
fx=d.module('evening_test_fixtures',HERE/'test_e5f_evening_calibration.py')

def synthetic_case(folder,c,objs,req):
    """Deterministic synthetic artifacts; never import the household runtime."""
    fx.success(folder,c,objs,req,pid=os.getpid())
    case=folder/'case';h=req['point']['H0']
    pair=(10.,10.) if h==c['initial_point']['H0'] else (2.,2.) if h in (.11,.12) else (2.,0.) if h==.21 else (0.,2.) if h==.22 else (3.,3.)
    rows=list(d.table(case/'target_fit.csv','moment').values())
    for row in rows:
        gap=pair[0] if row['moment']=='m0' else pair[1] if row['moment']=='m1' else 0.
        row.update(model=float(row['target'])+gap,gap=gap)
        if row['weight']!='':row['loss_contribution']=float(row['weight'])*gap*gap
    fx.csvwrite(case/'target_fit.csv',rows)
    with gzip.open(case/'initial_state.pkl.gz','wb') as f:pickle.dump({'synthetic_packet':True},f)
    receipt=d.read(case/'receipt.json');receipt['loss']=sum(float(r['loss_contribution']) for r in rows if r['loss_contribution']!='');receipt['case_checkpoint_sha256']=d.sha(case/'initial_state.pkl.gz');d.write(case/'receipt.json',receipt)
    side=d.read(folder/'success.json');side.update(loss=receipt['loss'],receipt_sha256=d.sha(case/'receipt.json'),checkpoint_sha256=receipt['case_checkpoint_sha256']);d.write(folder/'success.json',side)
    startup=d.read(folder/'startup.json');startup['parent_pid']=os.getppid();d.write(folder/'startup.json',startup)

def subprocess_worker():
    ap=argparse.ArgumentParser();ap.add_argument('--synthetic-worker',action='store_true');ap.add_argument('--stage');ap.add_argument('--contract',type=Path);ap.add_argument('--output',type=Path);ap.add_argument('--request',type=Path);a=ap.parse_args()
    c=d.read(a.contract);objs={lane:d.read(c['lanes'][lane]['objective']['path']) for lane in d.LANES};req=d.read(a.request)
    assert req['controller_pid']==os.getppid()
    if a.stage=='render':
        calls=[]
        def diagnostics(packet,out,validate_production_young):
            assert packet=={'synthetic_packet':True} and validate_production_young is False
            calls.append('standard_diagnostics');graphs=out/'standard_diagnostics';graphs.mkdir()
            for name in c['standard_diagnostic_names']:(graphs/name).write_bytes(b'synthetic graph')
        runtime=types.SimpleNamespace(setup=lambda *args:types.SimpleNamespace(rt={'audit':types.SimpleNamespace(standard_diagnostics=diagnostics)}))
        with mock.patch.object(d,'module',return_value=runtime):d.render_saved(a,c)
        assert calls==['standard_diagnostics']
    else:synthetic_case(a.output,c,objs,req)
    d.write(a.output/'synthetic_execution.json',dict(pid=os.getpid(),parent_pid=os.getppid(),stage=a.stage,model_solves=0))

class NightTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
        self.c,self.objs,self.contract=fx.fixture(self.root)
        # render_saved evaluates the pinned path before module() is substituted.
        self.c['files']['runtime']={'path':str(Path(__file__).resolve()),'sha256':d.sha(__file__)}
        self.c['budget'].update(absolute_start_epoch=1790550000,absolute_end_epoch=1790596800,total_seconds=46800,search_reserve_seconds=3600,export_reserve_seconds=600,max_search_cases=720,max_objective_cases=746,maximum_diagnostic_objectives=12,objective_cap_seconds=1800)
        self.c['initial_design']=[];self.c['seed']=20260928
    def test_fixed_end_and_two_reserves(self):
        a=d.clock(self.c,1790560000);b=d.clock(self.c,1790565000)
        self.assertEqual(a,b);self.assertEqual(a['search_cutoff'],1790593200);self.assertEqual(a['repeat_cutoff'],1790596200)
        with self.assertRaises(AssertionError):d.clock(self.c,1790596800)
    def test_proposals_never_call_uniform_fullbox(self):
        from unittest import mock
        seen=set()
        with mock.patch.object(d.random.Random,'uniform',side_effect=AssertionError('Full-box proposals forbidden')):
            rows=[d.propose(self.c,self.objs['primary'],'primary',self.c['initial_point'],i,seen) for i in range(240)]
        self.assertEqual(len(seen),240);self.assertTrue(all(.001<=r['tenure_choice_kappa']<=.1 for r in rows))
        self.assertTrue(all('broad' not in d.proposal_label(self.c,i) for i in range(240)))
    def test_seed_and_lane_independence(self):
        a=d.propose(self.c,self.objs['primary'],'primary',self.c['initial_point'],0,set())
        b=d.propose(self.c,self.objs['identity'],'identity',self.c['initial_point'],0,set())
        self.assertNotEqual(a,b)
        other=copy.deepcopy(self.c);other['seed']=20260927
        self.assertNotEqual(a,d.propose(other,self.objs['primary'],'primary',other['initial_point'],0,set()))
    def test_all_ten_coordinates_receive_coordinate_draws(self):
        names={d.proposal_label(self.c,i).removeprefix('coordinate_') for i in range(240) if d.proposal_label(self.c,i).startswith('coordinate_')}
        self.assertEqual(names,set(d.FREE))
    def test_common_primary_winner_gets_two_additional_repeats(self):
        best={lane:dict(lane=lane,point={'x':i},loss=1.) for i,lane in enumerate(d.LANES)}
        common=dict(lane='identity',point={'x':99},primary_rescore=.1)
        selected,key=d.final_selection(best,common)
        self.assertEqual(key,'common_primary');self.assertEqual(len(selected)*2,8);self.assertEqual(selected[key],common)
        selected,key=d.final_selection(best,best['block']);self.assertEqual(len(selected)*2,6);self.assertEqual(key,'block')
    def test_anchor_comparison_ignores_only_lane_weights(self):
        paths=[]
        for lane in ('primary','block'):
            path=self.root/lane;r=fx.req(self.c,self.contract,path,lane);fx.success(path,self.c,self.objs,r);paths.append(path/'case')
        d.compare_anchor(*paths)
        rows=list(d.table(paths[1]/'parameters.csv','parameter').values());rows[0]['estimate']=.49;fx.csvwrite(paths[1]/'parameters.csv',rows)
        with self.assertRaises(AssertionError):d.compare_anchor(*paths)
    def test_anchor_comparison_includes_validation_moments(self):
        paths=[]
        for name in ('a','b'):
            path=self.root/name;r=fx.req(self.c,self.contract,path);fx.success(path,self.c,self.objs,r);paths.append(path/'case')
        rows=list(d.table(paths[1]/'target_fit.csv','moment').values());rows[10]['model']=999;fx.csvwrite(paths[1]/'target_fit.csv',rows)
        with self.assertRaises(AssertionError):d.compare_anchor(*paths)
    def test_anchor_and_repeat_csv_mutation_fails_pins(self):
        cases=[]
        for name in ('anchor','repeat_one','repeat_two'):
            folder=self.root/name;fx.success(folder,self.c,self.objs,fx.req(self.c,self.contract,folder));case=folder/'case'
            cases.append(dict(case_path=str(case),**{key:dict(path=str(case/file),sha256=d.sha(case/file)) for key,file in (('receipt','receipt.json'),('target_fit','target_fit.csv'),('parameters','parameters.csv'))}))
        self.c['anchor']=dict(cases[0],repeats=cases[1:]);d.verify_anchor_pins(self.c)
        changed=Path(cases[2]['parameters']['path']);changed.write_text(changed.read_text()+'\n')
        with self.assertRaises(AssertionError):d.verify_anchor_pins(self.c)
    def test_real_subprocess_loop_four_selections_eight_repeats_and_hourly_render(self):
        """Real owned processes/run_batch; synthetic science and reporter only."""
        self.c['budget'].update(workers=6,max_search_cases=6)
        self.c['anchor']={'case_path':str(self.root/'anchor/case')}
        for lane in d.LANES:
            for row in self.objs[lane]['target_rows']:
                if row['restriction_id'] in ('m0','m1'):
                    row['actual_weight']=9. if (lane=='primary' and row['restriction_id']=='m0') or (lane=='identity' and row['restriction_id']=='m1') else 1.
            spec=self.c['lanes'][lane];d.write(spec['objective']['path'],self.objs[lane]);spec['objective']['sha256']=d.sha(spec['objective']['path']);spec['canonical_sha256']=d.canon(self.objs[lane]);spec['target_weight_fingerprint']=d.canon(self.objs[lane]['target_rows'])
        d.write(self.contract,self.c)
        synthetic_case(self.root/'anchor',self.c,self.objs,fx.req(self.c,self.contract,self.root/'anchor','block'))
        original_time=time.time;offset=[0.];started=[]
        def managed(command,log,deadline,env):
            self.assertTrue(all(env[key]=='1' for key in d.THREADS))
            rewritten=[sys.executable,str(Path(__file__).resolve()),'--synthetic-worker']+command[2:]
            process=fx.supervisor.ManagedProcess(rewritten,log,deadline,env);started.append(process.process.pid);return process
        def run_batch(candidates,**kwargs):
            result=fx.supervisor.run_batch(candidates,**kwargs)
            if candidates and candidates[0]['context']['stage']=='initial':offset[0]=3601.
            return result
        proxy=types.SimpleNamespace(ManagedProcess=managed,run_batch=run_batch)
        fake_time=types.SimpleNamespace(time=lambda:original_time()+offset[0])
        def planned(c,obj,lane,center,index,seen):
            point=dict(c['initial_point']);point['H0']={'primary':(.11,.12),'identity':(.21,.22),'block':(.31,.32)}[lane][index];seen.add(d.identity(c,lane,point));return point
        now=original_time();timing=dict(start=now,end=now+9000,search_cutoff=now+7200,repeat_cutoff=now+8400)
        a=argparse.Namespace(stage='smoke',contract=self.contract,output=self.root/'smoke')
        with mock.patch.object(d,'module',return_value=proxy),mock.patch.object(d,'verify',return_value=(self.c,self.objs)),mock.patch.object(d,'clock',return_value=timing),mock.patch.object(d,'time',fake_time),mock.patch.object(d,'propose',side_effect=planned):
            d.controller(a,self.c,self.objs)
            smoke=a.output/'complete.json';self.assertEqual(d.read(smoke)['status'],'exact_loop_smoke_passed')
            approval=self.root/'approval.json';d.write(approval,dict(status='approved_search',contract_sha256=d.sha(self.contract),smoke_receipt=dict(path=str(smoke),sha256=d.sha(smoke))))
            a=argparse.Namespace(stage='search',contract=self.contract,output=self.root/'search',approval=approval,approval_sha256=d.sha(approval))
            try:d.controller(a,self.c,self.objs)
            except Exception:
                retained=Path(os.environ.get('NIGHT_TEST_FAILURE_DIR',tempfile.gettempdir()))/f"night_fixture_{os.environ['SLURM_JOB_ID']}_{os.getpid()}"
                retained.parent.mkdir(parents=True,exist_ok=True);shutil.copytree(self.root,retained);print('RETAINED SYNTHETIC FIXTURE',retained)
                for path in sorted(a.output.glob('*.log')):
                    if path.stat().st_size:print('SUBPROCESS LOG',path.name,path.read_text()[-6000:])
                for filename in ('latest_completed.json','hourly_reports.json'):
                    path=a.output/filename
                    if path.exists():print('CONTROLLER EVIDENCE',filename,path.read_text()[-10000:])
                raise
        complete=d.read(a.output/'complete.json');self.assertEqual(complete['status'],'bounded_search_complete');self.assertEqual(len(complete['selected']),4);self.assertEqual(complete['common_primary_key'],'common_primary')
        self.assertEqual(complete['selected']['common_primary']['point']['H0'],.22)
        repeats=[r for r in complete['records'] if r['case'].startswith('repeat_')];self.assertEqual(len(repeats),8)
        self.assertEqual(len(complete['hourly_reports']),3);self.assertTrue(all(r['status']=='success' for r in complete['hourly_reports']))
        for key,row in complete['selected'].items():
            exported=a.output/'selected_export'/key;d.compare_tables(row['case_path'],exported)
            self.assertEqual(len(list((exported/'standard_diagnostics').glob('*.png'))),17)
            self.assertEqual(len(d.read(exported/'export_receipt.json')['repeats']),2)
        executions=[d.read(p) for p in self.root.glob('**/synthetic_execution.json')]
        self.assertEqual(len(executions),23);self.assertEqual(len(set(started)),23);self.assertTrue(all(x['parent_pid']==os.getpid() and x['model_solves']==0 for x in executions))
        print(json.dumps(dict(test='real_subprocess_night_controller',owned_processes=23,smokes=6,search_objectives=6,final_repeats=8,hourly_renders=3,selected_cases=4,model_solves=0,status='passed')))
if __name__=='__main__':
    if '--synthetic-worker' in sys.argv:subprocess_worker()
    else:unittest.main()
