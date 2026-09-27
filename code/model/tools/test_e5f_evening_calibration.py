#!/usr/bin/env python3
"""Torch-only synthetic evening-controller tests; no household solves."""
import argparse,csv,importlib.util,os,signal,sys,tempfile,time,types,unittest
from pathlib import Path
from unittest import mock
if not os.environ.get('SLURM_JOB_ID','').isdigit():raise RuntimeError('Torch Slurm only')
p=Path(__file__).with_name('run_e5f_evening_calibration.py');s=importlib.util.spec_from_file_location('evening_controller_test',p);d=importlib.util.module_from_spec(s);s.loader.exec_module(d)
reviewed=Path(os.environ['REVIEWED_RECOVERY_SEARCH']);sys.path.insert(0,str(reviewed.parent));supervisor=d.module('reviewed_evening_test_supervisor',reviewed)
builder=d.module('evening_contract_builder_test',Path(__file__).with_name('prepare_e5f_evening_contract.py'))

def csvwrite(path,rows):
    with Path(path).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def fixture(root):
    bounds=[dict(parameter=k,lower=.001 if k=='tenure_choice_kappa' else 0.,upper=.1 if k=='tenure_choice_kappa' else 1.) for k in d.FREE]
    now=time.time();c=dict(schema=d.SCHEMA,seed=27,initial_point={k:.01 if k=='tenure_choice_kappa' else .5 for k in d.FREE},normalization=dict(initial_psi=.1,initial_step=.03,warm_price=True,maximum_stationary_solves=23),validation_rows=['m10','m11','m12'],fixed=dict(due=True,delta_alpha=0,sigma=2),proposal_widths=dict.fromkeys(d.FREE,.1),log_parameters=['tenure_choice_kappa'],standard_diagnostic_names=[f'g{i:02}.png' for i in range(17)],budget=dict(absolute_start_epoch=now-1,absolute_end_epoch=now-1+21600,total_seconds=21600,workers=6,max_search_cases=3,max_objective_cases=384,maximum_native_prerequisite_cases=2,objective_cap_seconds=900,search_reserve_seconds=2700,export_reserve_seconds=300),files={'driver':dict(path=str(p),sha256=d.sha(p)),'recovery_search':dict(path=str(reviewed),sha256=d.sha(reviewed))},lanes={})
    manifest=root/'manifest.json';d.write(manifest,{'files':{}});c['source_manifest']=dict(path=str(manifest),sha256=d.sha(manifest));objs={}
    for lane,mult in zip(d.LANES,(1.,2.,3.)):
        rows=[dict(restriction_id=f'm{i}',target=1.,actual_weight=mult if i<10 else 0. if i<13 else None) for i in range(14)]
        obj=dict(target_rows=rows,parameter_restrictions=bounds);path=root/f'{lane}.json';d.write(path,obj);objs[lane]=obj
        c['lanes'][lane]=dict(objective=dict(path=str(path),sha256=d.sha(path)),canonical_sha256=d.canon(obj),target_weight_fingerprint=d.canon(rows))
    contract=root/'contract.json';d.write(contract,c);return c,objs,contract

def req(c,contract,folder,lane='primary',point=None,graphs=True):
    point=dict(point or c['initial_point']);ctx=dict(candidate_id=folder.name,stage='smoke',contract_sha256=d.sha(contract),source_sha256=c['files']['driver']['sha256'],target_sha256=c['lanes'][lane]['objective']['sha256'],point_sha256=d.canon(point))
    return dict(point=point,lane=lane,context=ctx,contract_sha256=d.sha(contract),scientific_candidate_id=d.identity(c,lane,point),normalization_inputs=c['normalization'],controller_pid=os.getpid(),graphs=graphs)
def success(folder,c,objs,r,pid=123):
    folder.mkdir();case=folder/'case';case.mkdir();d.write(folder/'startup.json',dict(context=r['context'],pid=pid,parent_pid=os.getpid(),normalization_inputs=c['normalization']))
    (case/'initial_state.pkl.gz').write_bytes(b'synthetic');d.write(case/'stationary_solves.json',[dict(status='completed',psi_child=.1),dict(status='completed',psi_child=.13)])
    rows=[]
    for x in objs[r['lane']]['target_rows']:
        w=x['actual_weight'];rows.append(dict(moment=x['restriction_id'],target=1.,model=2.,gap=1.,weight='' if w is None else w,loss_contribution='' if w is None else w,role='validation' if w==0 else 'fit'))
    csvwrite(case/'target_fit.csv',rows)
    bounds={x['parameter']:x for x in objs[r['lane']]['parameter_restrictions']};params=[dict(parameter=k,estimate=r['point'][k],lower=bounds[k]['lower'],upper=bounds[k]['upper'],near_bound=False) for k in d.FREE]
    params.extend(dict(parameter=f'fixed{i}',estimate=1.,lower='',upper='',near_bound='') for i in range(21));csvwrite(case/'parameters.csv',params)
    if r['graphs']:
        graphs=case/'standard_diagnostics';graphs.mkdir()
        for name in c['standard_diagnostic_names']:(graphs/name).write_bytes(b'synthetic')
    receipt=dict(contract_sha256=r['contract_sha256'],scientific_identity=d.science(c),scientific_candidate_id=d.identity(c,r['lane'],r['point']),lane=r['lane'],point=r['point'],normalization_inputs=c['normalization'],normalization=dict(psi_child=.13,stationary_solves=2),objective_stationary_solves=2,target_system_sha256=c['lanes'][r['lane']]['objective']['sha256'],target_weight_fingerprint=c['lanes'][r['lane']]['target_weight_fingerprint'],source_manifest_sha256=c['source_manifest']['sha256'],case_checkpoint_sha256=d.sha(case/'initial_state.pkl.gz'),loss=sum(x['loss_contribution'] for x in rows if x['loss_contribution']!=''))
    d.write(case/'receipt.json',receipt);d.write(folder/'success.json',dict(context=r['context'],receipt_sha256=d.sha(case/'receipt.json'),checkpoint_sha256=receipt['case_checkpoint_sha256'],loss=receipt['loss']))
    return d.validate(folder,c,objs,r)

class Tests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name);self.c,self.objs,self.contract=fixture(self.root)
    def test_verify_pins_and_three_validation_rows(self):
        with mock.patch.dict(os.environ,EXPECTED_E5F_EVENING_SHA256=d.sha(self.contract)):
            d.verify(self.contract);d.write(self.c['lanes']['primary']['objective']['path'],{})
            with self.assertRaises(AssertionError):d.verify(self.contract)
    def test_clock_never_restarts(self):
        a=d.clock(self.c,time.time());b=d.clock(self.c,time.time()+1000);self.assertEqual(a,b);self.assertEqual(a['search_cutoff'],a['end']-2700)
        with self.assertRaises(AssertionError):d.clock(self.c,a['end'])
    def test_losses_stay_separate_primary_rescore_common(self):
        out=[]
        for lane in d.LANES:
            folder=self.root/lane;r=req(self.c,self.contract,folder,lane);out.append(success(folder,self.c,self.objs,r))
        self.assertEqual([r['loss'] for r in out],[10.,20.,30.]);self.assertEqual([r['primary_rescore'] for r in out],[10.,10.,10.])
    def test_graph_names_not_only_count(self):
        folder=self.root/'one';r=req(self.c,self.contract,folder);success(folder,self.c,self.objs,r)
        graphs=folder/'case/standard_diagnostics';(graphs/'g00.png').rename(graphs/'renamed.png')
        with self.assertRaises(AssertionError):d.validate(folder,self.c,self.objs,r)
    def test_complete_numeric_repeats(self):
        folders=[self.root/'a',self.root/'b']
        for folder in folders:success(folder,self.c,self.objs,req(self.c,self.contract,folder))
        d.compare_tables(folders[0]/'case',folders[1]/'case')
        rows=list(csv.DictReader((folders[1]/'case/target_fit.csv').open()));rows[10]['model']='9';csvwrite(folders[1]/'case/target_fit.csv',rows)
        with self.assertRaises(AssertionError):d.compare_tables(folders[0]/'case',folders[1]/'case')
    def test_owner_timeout_and_unknown_failure(self):
        folder=self.root/'one';folder.mkdir();r=req(self.c,self.contract,folder)
        proc=types.SimpleNamespace(process=types.SimpleNamespace(pid=123),observed_running_at_expiry=False,deadline_kill_reaped=False,deadline_expired=True)
        self.assertEqual(d.classify(folder,self.c,self.objs,r,proc,-signal.SIGKILL)[0],'fatal')
        proc.observed_running_at_expiry=True;proc.deadline_kill_reaped=True
        self.assertEqual(d.classify(folder,self.c,self.objs,r,proc,-signal.SIGKILL)[0],'censored_timeout')
        d.write(folder/'failure.json',dict(context=r['context'],status='fatal'))
        self.assertEqual(d.classify(folder,self.c,self.objs,r,proc,-signal.SIGKILL)[0],'fatal')
    def test_proposals_finite_log_bound_and_lane_independent(self):
        seen=set();rows=[d.propose(self.c,self.objs['primary'],'primary',self.c['initial_point'],i,seen) for i in range(128)]
        self.assertEqual(len(seen),128);self.assertTrue(all(.001<=r['tenure_choice_kappa']<=.1 for r in rows))
        self.assertNotEqual(d.identity(self.c,'primary',rows[0]),d.identity(self.c,'identity',rows[0]))
    def test_exact_lane_weight_formulas(self):
        ids=list(builder.SCALE_FLOORS)+list(builder.VALIDATION)+['initial_normalization']
        original=dict(target_rows=[dict(restriction_id=k,target=.007 if k=='bequest_wealth' else 2.,actual_weight=None if k=='initial_normalization' else 120.) for k in ids],parameter_restrictions=[dict(parameter=k,lower=0.,upper=1.) for k in d.FREE[:-1]])
        lanes=builder.lane_objectives(original);rows={lane:{r['restriction_id']:r for r in obj['target_rows']} for lane,obj in lanes.items()}
        self.assertEqual(rows['identity']['bequest_wealth']['actual_weight'],10000.)
        self.assertEqual(rows['block']['cps_childlessness']['actual_weight'],10.)
        self.assertEqual(rows['block']['wealth_earnings']['actual_weight'],20.)
        self.assertEqual(rows['primary']['mean_rooms']['actual_weight'],120.)
        for lane in d.LANES:
            for k in builder.VALIDATION:self.assertEqual(rows[lane][k]['actual_weight'],0.)
        self.assertEqual(original['target_rows'][0]['actual_weight'],120.)
    def test_declared_profiles_consume_proposal_slots(self):
        bounds={r['parameter']:(r['lower'],r['upper']) for r in self.objs['primary']['parameter_restrictions']}
        design,omitted=builder.initial_design(self.c['initial_point'],bounds,self.c['proposal_widths'])
        self.c['initial_design']=design;self.assertFalse(omitted);seen=set()
        for i,row in enumerate(design):
            self.assertEqual(d.propose(self.c,self.objs['primary'],'primary',dict(self.c['initial_point'],chi=.6),i,seen),row['point'])
            self.assertEqual(d.proposal_label(self.c,i),row['label'])
        self.assertEqual([x['point']['tenure_choice_kappa'] for x in design[:3]],[.05,.004,.006])
    def test_exact_controller_smoke_search_repeat_and_export(self):
        owner=self;started=[];pid=[100]
        class Process:
            def __init__(self,command,log,deadline,env):
                pid[0]+=1;self.process=types.SimpleNamespace(pid=pid[0]);self.deadline=deadline;self.observed_running_at_expiry=False;self.deadline_kill_reaped=False;self.deadline_expired=False
                r=d.read(command[command.index('--request')+1]);folder=Path(command[command.index('--output')+1]);started.append(r)
                owner.assertTrue(all(env[k]=='1' for k in d.THREADS));success(folder,owner.c,owner.objs,r,pid[0])
            def poll(self):return 0
            def close(self):pass
        proxy=types.SimpleNamespace(ManagedProcess=Process,run_batch=supervisor.run_batch)
        a=argparse.Namespace(stage='smoke',contract=self.contract,output=self.root/'smoke')
        with mock.patch.object(d,'module',return_value=proxy),mock.patch.object(d,'verify',return_value=(self.c,self.objs)):
            d.controller(a,self.c,self.objs);receipt=a.output/'complete.json';self.assertEqual(len(started),6)
            approval=self.root/'approval.json';d.write(approval,dict(status='approved_search',contract_sha256=d.sha(self.contract),smoke_receipt=dict(path=str(receipt),sha256=d.sha(receipt))))
            a=argparse.Namespace(stage='search',contract=self.contract,output=self.root/'search',approval=approval,approval_sha256=d.sha(approval));d.controller(a,self.c,self.objs)
        self.assertEqual(len(started),15);complete=d.read(a.output/'complete.json');self.assertEqual(complete['status'],'bounded_search_complete')
        for lane in d.LANES:
            original=complete['selected'][lane];self.assertIn('/smoke/',original['case_path']);self.assertEqual(len(list((a.output/'selected_export'/lane/'standard_diagnostics').glob('*.png'))),17)
    def test_search_rejects_unapproved_file(self):
        approval=self.root/'approval.json';d.write(approval,dict(status='not_approved',contract_sha256=d.sha(self.contract)))
        a=argparse.Namespace(stage='search',contract=self.contract,output=self.root/'search',approval=approval,approval_sha256=d.sha(approval))
        with self.assertRaises(AssertionError):d.controller(a,self.c,self.objs)
        self.assertFalse(a.output.exists())
    def test_reviewed_loop_fatal_stops_dispatch_drains_started(self):
        started=[];closed=[]
        class Process:
            def __init__(self,i):self.i=i
            def poll(self):return 0
            def close(self):closed.append(self.i)
        def launch(r,end):started.append(r['id']);return Process(r['id'])
        result=supervisor.run_batch([{'id':i} for i in range(3)],workers=2,deadline=time.time()+10,launch=launch,finish=lambda r,p,c:dict(status='fatal' if r['id']==0 else 'success',halt_new_dispatch=r['id']==0),heartbeat=lambda **kw:None,poll_seconds=0)
        self.assertEqual(started,[0,1]);self.assertEqual(closed,[0,1]);self.assertEqual(result['unrun_ids'],[2])
if __name__=='__main__':unittest.main()
