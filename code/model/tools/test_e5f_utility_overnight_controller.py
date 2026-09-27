#!/usr/bin/env python3
"""Synthetic controller checks; no household solves.

REVIEWED_RECOVERY_SEARCH must point to the pinned reviewed supervisor. Its real
run_batch loop is exercised with synthetic owned-process objects.
Local checks require ALLOW_LOCAL_CONTROLLER_TESTS=1; never spoof Slurm metadata.
"""
from __future__ import annotations
import argparse,csv,importlib.util,json,os,signal,sys,tempfile,time,types,unittest
from pathlib import Path
from unittest import mock

if not os.environ.get('SLURM_JOB_ID') and os.environ.get('ALLOW_LOCAL_CONTROLLER_TESTS')!='1':
    raise RuntimeError('Run on Slurm or explicitly authorize local synthetic checks')

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('overnight_test_subject',HERE/'run_e5f_utility_overnight_calibration.py')
d=importlib.util.module_from_spec(spec);spec.loader.exec_module(d)
reviewed=Path(os.environ['REVIEWED_RECOVERY_SEARCH']).resolve()
sys.path.insert(0,str(reviewed.parent))
supervision=d.module('overnight_test_reviewed_supervision',reviewed)


def table(path,rows):
    with Path(path).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def fixture(root):
    obj={'target_rows':[{'restriction_id':'fit','target':2.,'actual_weight':1.}, {'restriction_id':'normalization','target':2.1,'actual_weight':None}],
         'parameter_restrictions':[{'parameter':k,'lower':0.,'upper':1.} for k in d.FREE]}
    objective=root/'objective.json';d.write(objective,obj)
    c={'schema':d.SCHEMA,'status':'reviewed_smoke','objective':{'path':str(objective),'sha256':d.sha(objective)},
       'objective_canonical_sha256':d.canon(obj),'target_weight_fingerprint':d.canon(obj['target_rows']),
       'normalization':{'initial_psi':.1,'initial_step':.03,'warm_price':True,'maximum_stationary_solves':23},
       'initial_point':dict.fromkeys(d.FREE,.5),'files':{'driver':{'sha256':'driver'},'recovery_policy':{'path':str(reviewed.parent/'e5f_utility_recovery_policy_v1.py')},'recovery_search':{'path':str(reviewed)}},
       'runtime_tools':str(reviewed.parent),'seed':2026,'proposal_widths':dict.fromkeys(d.FREE,.1),
       'budget':{'total_seconds':28800,'search_seconds':23400,'repeat_seconds':3600,'export_seconds':1800,'smoke_seconds':120,'objective_cap_seconds':60,'rounds':1,'points_per_round':2,'workers':2}}
    contract=root/'contract.json';d.write(contract,c)
    return c,obj,contract


def request(c,contract,folder,point=None,stage='smoke',graphs=True):
    p=dict(point or c['initial_point']);context=dict(candidate_id=folder.name,stage=stage,contract_sha256=d.sha(contract),source_sha256='driver',target_sha256=c['objective']['sha256'],point_sha256=d.canon(p))
    return dict(id=folder.name,point=p,context=context,stage=stage,point_sha256=d.canon(p),scientific_candidate_id=d.candidate_id(c,p),normalization_inputs=d.norm_inputs(c),contract_sha256=d.sha(contract),controller_pid=os.getpid(),deadline_epoch=time.time()+100,graphs=graphs)


def success(folder,c,req,pid=123,loss=1.):
    folder.mkdir();case=folder/'case';case.mkdir()
    d.write(folder/'startup.json',dict(context=req['context'],pid=pid,parent_pid=os.getpid(),normalization_inputs=d.norm_inputs(c)))
    (case/'initial_state.pkl.gz').write_bytes(b'synthetic checkpoint')
    ledger=[dict(status='completed',psi_child=.1),dict(status='completed',psi_child=.13)]
    d.write(case/'stationary_solves.json',ledger)
    table(case/'target_fit.csv',[dict(moment='fit',target=2.,model=3.,gap=1.,weight=1.,loss_contribution=loss),dict(moment='normalization',target=2.1,model=2.1,gap=0.,weight='',loss_contribution='')])
    table(case/'parameters.csv',[dict(parameter=k,estimate=req['point'][k],lower=0.,upper=1.,near_bound='False') for k in d.FREE])
    if req['graphs']:
        graphs=case/'standard_diagnostics';graphs.mkdir()
        for i in range(17):(graphs/f'graph_{i:02}.png').write_bytes(b'synthetic')
    receipt=dict(point=req['point'],loss=loss,normalization=dict(psi_child=.13,stationary_solves=2),normalization_inputs=d.norm_inputs(c),objective_stationary_solves=2,overnight_contract_sha256=req['contract_sha256'],scientific_identity=d.science_id(c),scientific_candidate_id=d.candidate_id(c,req['point']),target_system_sha256=c['objective']['sha256'],target_weight_fingerprint=c['target_weight_fingerprint'],objective_canonical_sha256=c['objective_canonical_sha256'],case_checkpoint_sha256=d.sha(case/'initial_state.pkl.gz'))
    d.write(case/'receipt.json',receipt)
    d.write(folder/'success.json',dict(context=req['context'],loss=loss,receipt_sha256=d.sha(case/'receipt.json'),checkpoint_sha256=receipt['case_checkpoint_sha256']))
    return dict(status='success',point=req['point'],**d.validate_success(folder,c,req))


class ControllerTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
        self.c,self.obj,self.contract=fixture(self.root)
    def test_approval_does_not_change_science_but_normalization_does(self):
        c=json.loads(json.dumps(self.c));c['status']='approved_production';c['approval']={'reviewer':'test'};c['production_blockers']=[]
        self.assertEqual(d.science_id(c),d.science_id(self.c));c['normalization']['initial_step']=.04
        self.assertNotEqual(d.candidate_id(c,c['initial_point']),d.candidate_id(self.c,self.c['initial_point']))
    def test_local_execution_is_explicit_host_pinned_and_not_slurm(self):
        self.c['execution']=dict(kind='local',hostname='approved-host',authorization_id='night-run')
        with mock.patch.dict(os.environ,{'E5F_LOCAL_EXECUTION_AUTHORIZATION':'night-run'},clear=True),mock.patch.object(d.socket,'gethostname',return_value='approved-host'):
            d.verify_execution(self.c)
            with mock.patch.dict(os.environ,{'SLURM_JOB_ID':'123'}):
                with self.assertRaises(AssertionError):d.verify_execution(self.c)
            with mock.patch.object(d.socket,'gethostname',return_value='other-host'):
                with self.assertRaises(AssertionError):d.verify_execution(self.c)
            os.environ['E5F_LOCAL_EXECUTION_AUTHORIZATION']='wrong'
            with self.assertRaises(AssertionError):d.verify_execution(self.c)
    def test_legacy_execution_remains_slurm_only(self):
        with mock.patch.dict(os.environ,{},clear=True):
            with self.assertRaises(AssertionError):d.verify_execution(self.c)
        with mock.patch.dict(os.environ,{'SLURM_JOB_ID':'123'},clear=True):d.verify_execution(self.c)
    def test_absolute_deadline_preserves_repeat_and_export_reserves(self):
        start=100000.;self.c['budget']['absolute_end_epoch']=start+20000
        clock=d.build_clock(self.c,start,'search','contract')
        self.assertEqual(clock['end'],start+20000)
        self.assertEqual(clock['repeat_cutoff'],start+18200)
        self.assertEqual(clock['search_cutoff'],start+14600)
        self.assertEqual(d.build_clock(self.c,start,'smoke','contract')['search_cutoff'],start+120)
        with self.assertRaises(AssertionError):d.build_clock(self.c,start+20000,'search','contract')
    def test_clean_runtime_seam(self):
        expected=tuple(range(8));runtime=types.SimpleNamespace(setup=mock.Mock(return_value=expected))
        self.c['files']['calibration_runtime']={'path':'explicit-runtime.py'}
        with mock.patch.object(d,'module',return_value=runtime):self.assertEqual(d.setup(self.c,self.obj,{},self.root),expected)
        runtime.setup.assert_called_once_with(self.c,self.obj,{},self.root)
    def test_contract_and_target_pins_fail_closed(self):
        self.c['files']['driver']['path']=d.__file__
        for item in self.c['files'].values():item['sha256']=d.sha(item['path'])
        base=self.root/'base.json';d.write(base,{'fixed':'reference'})
        self.c.update(base_contract={'path':str(base),'sha256':d.sha(base)},fixed={'delta_alpha':0,'sigma':2},pending_observer_mismatches=['declared'],economic_changes=['declared'])
        d.write(self.contract,self.c)
        with mock.patch.dict(os.environ,{'EXPECTED_UTILITY_OVERNIGHT_SHA256':d.sha(self.contract)}):
            d.verify(self.contract)
            changed=dict(self.obj);changed['target_rows']=[dict(self.obj['target_rows'][0],target=9.),self.obj['target_rows'][1]]
            d.write(self.c['objective']['path'],changed)
            with self.assertRaises(AssertionError):d.verify(self.contract)
    def test_worker_caps_follow_execution_contract(self):
        self.c['files']['driver']['path']=d.__file__
        for item in self.c['files'].values():item['sha256']=d.sha(item['path'])
        base=self.root/'base.json';d.write(base,{'fixed':'reference'})
        self.c.update(base_contract={'path':str(base),'sha256':d.sha(base)},fixed={'delta_alpha':0,'sigma':2},pending_observer_mismatches=['declared'],economic_changes=['declared'])
        for kind,workers,valid in [('slurm',24,True),('slurm',25,False),('local',10,True),('local',11,False)]:
            self.c['execution']=dict(kind=kind,hostname='host',authorization_id='token')
            self.c['budget'].update(workers=workers,points_per_round=workers)
            d.write(self.contract,self.c)
            with mock.patch.dict(os.environ,{'EXPECTED_UTILITY_OVERNIGHT_SHA256':d.sha(self.contract)}):
                if valid:d.verify(self.contract)
                else:
                    with self.assertRaises(AssertionError):d.verify(self.contract)
    def pair(self):
        cases=[]
        for name in ('original','repeat'):
            folder=self.root/name;req=request(self.c,self.contract,folder);record=success(folder,self.c,req);cases.append((folder,req,record))
        return cases
    def test_every_target_cell_checked_even_equal_loss(self):
        (a,_,_),(b,_,_)=self.pair();d.compare_tables(a/'case',b/'case')
        rows=list(csv.DictReader((b/'case/target_fit.csv').open()));rows[1]['model']='2.2';table(b/'case/target_fit.csv',rows)
        with self.assertRaises(AssertionError):d.compare_tables(a/'case',b/'case')
    def test_parameter_bound_and_duplicate_checked(self):
        (a,_,_),(b,_,_)=self.pair();rows=list(csv.DictReader((b/'case/parameters.csv').open()));rows[0]['upper']='2';table(b/'case/parameters.csv',rows)
        with self.assertRaises(AssertionError):d.compare_tables(a/'case',b/'case')
        table(b/'case/parameters.csv',rows+[rows[0]])
        with self.assertRaises(AssertionError):d.keyed_csv(b/'case/parameters.csv','parameter')
    def test_checkpoint_and_normalization_audited(self):
        folder=self.root/'one';req=request(self.c,self.contract,folder);success(folder,self.c,req)
        ledger=d.read(folder/'case/stationary_solves.json');ledger[0]['psi_child']=.11;d.write(folder/'case/stationary_solves.json',ledger)
        with self.assertRaises(AssertionError):d.validate_success(folder,self.c,req)
        ledger[0]['psi_child']=.1;d.write(folder/'case/stationary_solves.json',ledger);(folder/'case/initial_state.pkl.gz').write_bytes(b'changed')
        with self.assertRaises(AssertionError):d.validate_success(folder,self.c,req)
    def test_solve_cap_checked(self):
        folder=self.root/'one';req=request(self.c,self.contract,folder);success(folder,self.c,req)
        receipt=d.read(folder/'case/receipt.json');receipt['objective_stationary_solves']=24;d.write(folder/'case/receipt.json',receipt)
        side=d.read(folder/'success.json');side['receipt_sha256']=d.sha(folder/'case/receipt.json');d.write(folder/'success.json',side)
        with self.assertRaises(AssertionError):d.validate_success(folder,self.c,req)
    def test_only_beta_roundtrip_has_small_tolerance(self):
        folder=self.root/'one';req=request(self.c,self.contract,folder);success(folder,self.c,req)
        rows=list(csv.DictReader((folder/'case/parameters.csv').open()))
        for row in rows:
            if row['parameter']=='beta_annual':row['estimate']='0.5000000000000001'
        table(folder/'case/parameters.csv',rows);d.validate_success(folder,self.c,req)
        for row in rows:
            if row['parameter']=='chi':row['estimate']='0.5000000000000001'
        table(folder/'case/parameters.csv',rows)
        with self.assertRaises(AssertionError):d.validate_success(folder,self.c,req)
    def test_timeout_requires_positive_owned_evidence(self):
        folder=self.root/'missing';folder.mkdir();req=request(self.c,self.contract,folder);req['request']=dict(req)
        proc=types.SimpleNamespace(process=types.SimpleNamespace(pid=123),observed_running_at_expiry=False,deadline_kill_reaped=False,deadline_expired=True)
        self.assertEqual(d.classify(folder,self.c,req,proc,-signal.SIGKILL)[0],'fatal')
        proc.observed_running_at_expiry=True;proc.deadline_kill_reaped=True
        self.assertEqual(d.classify(folder,self.c,req,proc,-signal.SIGKILL)[0],'censored_timeout')
        d.write(folder/'failure.json',dict(context=req['context'],status='fatal'))
        self.assertEqual(d.classify(folder,self.c,req,proc,-signal.SIGKILL)[0],'fatal')
    def test_zero_exit_cannot_be_inadmissible(self):
        folder=self.root/'one';req=request(self.c,self.contract,folder);success(folder,self.c,req);req['request']=dict(req)
        d.write(folder/'failure.json',dict(context=req['context'],status='inadmissible',classification='existing_native_prefix'))
        proc=types.SimpleNamespace(process=types.SimpleNamespace(pid=123),observed_running_at_expiry=False,deadline_kill_reaped=False,deadline_expired=False)
        self.assertEqual(d.classify(folder,self.c,req,proc,0)[0],'fatal')
    def test_proposals_are_bounded_deterministic_and_distinct(self):
        seen={d.candidate_id(self.c,self.c['initial_point'])}
        a=d.proposals(self.c,self.obj,self.c['initial_point'],10,0,set(seen));b=d.proposals(self.c,self.obj,self.c['initial_point'],10,0,set(seen))
        self.assertEqual(a,b);self.assertEqual(len({d.canon(x) for x in a}),10)
        self.assertTrue(all(0<=v<=1 for row in a for v in row.values()))
    def test_reviewed_batch_drains_siblings_after_fatal(self):
        launched=[];closed=[]
        class Process:
            def __init__(self,name):self.name=name
            def poll(self):return 0
            def close(self):closed.append(self.name)
        def launch(req,deadline):launched.append(req['id']);return Process(req['id'])
        def finish(req,process,code):return dict(id=req['id'],status='fatal' if req['id']==0 else 'success',halt_new_dispatch=req['id']==0)
        result=supervision.run_batch([{'id':i} for i in range(3)],workers=2,deadline=time.time()+10,launch=launch,finish=finish,heartbeat=lambda **kw:None,poll_seconds=0)
        self.assertEqual(launched,[0,1]);self.assertEqual(sorted(closed),[0,1]);self.assertEqual(result['unrun_ids'],[2])
    def test_exact_controller_loop_smoke_reuse_and_frozen_export(self):
        owner=self;starts=[];fake_pid=[1000]
        class Process:
            def __init__(self,command,log_path,deadline,env):
                self.deadline=deadline;self.observed_epoch=time.time();self.observed_running_at_expiry=False;self.deadline_kill_reaped=False;self.deadline_expired=False
                fake_pid[0]+=1;self.process=types.SimpleNamespace(pid=fake_pid[0])
                folder=Path(command[command.index('--output')+1]);req=d.read(command[command.index('--request')+1]);starts.append(req)
                owner.assertTrue(all(env[k]=='1' for k in d.THREADS));success(folder,owner.c,req,pid=self.process.pid)
            def poll(self):return 0
            def close(self):pass
        proxy=types.SimpleNamespace(ManagedProcess=Process,run_batch=supervision.run_batch)
        smoke=argparse.Namespace(stage='smoke',output=self.root/'smoke',contract=self.contract,smoke_receipt=None,smoke_sha256=None)
        with mock.patch.object(d,'module',return_value=proxy),mock.patch.object(d,'verify',return_value=(self.c,self.obj)):
            d.controller(smoke,self.c,self.obj)
            receipt=smoke.output/'complete.json';self.assertEqual(len(starts),2);d.load_smoke(receipt,d.sha(receipt),self.c)
            production=dict(self.c,status='approved_production');production_path=self.root/'production.json';d.write(production_path,production)
            search=argparse.Namespace(stage='search',output=self.root/'search',contract=production_path,smoke_receipt=receipt,smoke_sha256=d.sha(receipt))
            # All synthetic losses tie. The smoke seed remains the canonical selection.
            d.controller(search,production,self.obj)
        self.assertEqual(len(starts),6) # two smoke + two new search + two repeats
        selected=d.read(search.output/'selected.json')['selected'];complete=d.read(search.output/'complete.json')
        self.assertEqual(selected['case_path'],d.read(receipt)['records'][0]['case_path'])
        self.assertEqual(selected,complete['selected']);self.assertEqual(complete['status'],'bounded_search_complete')
        self.assertEqual(len(list((search.output/'selected_export/standard_diagnostics').glob('*.png'))),17)
        self.assertEqual(d.read(search.output/'selected_export/receipt.json'),d.read(Path(selected['case_path'])/'receipt.json'))
    def run_interrupted_search(self,mode):
        owner=self;starts=[];fake_pid=[1000]
        source=self.root/'source';source.mkdir();(source/'model.py').write_text('pinned source')
        manifest=self.root/'source_manifest.json';d.write(manifest,dict(files={'model.py':d.sha(source/'model.py')}))
        self.c.update(source_root=str(source),source_manifest=dict(path=str(manifest),sha256=d.sha(manifest)))
        d.write(self.contract,self.c)
        class Process:
            def __init__(self,command,log_path,deadline,env):
                self.deadline=deadline;self.observed_epoch=time.time();self.observed_running_at_expiry=False;self.deadline_kill_reaped=False;self.deadline_expired=False
                fake_pid[0]+=1;self.process=types.SimpleNamespace(pid=fake_pid[0]);self.code=0
                folder=Path(command[command.index('--output')+1]);req=d.read(command[command.index('--request')+1]);starts.append(req)
                success(folder,owner.c,req,pid=self.process.pid)
                if req['stage']=='initial' and folder.name.endswith('0000'):
                    self.code=1
                    d.write(folder/'failure.json',dict(context=req['context'],status='fatal',classification='unknown_or_integrity_failure',phase='setup' if mode=='setup' else 'objective',error_type='InfeasibleThetaError',error='Unfamiliar native census'))
                    if mode=='source_tamper':(source/'model.py').write_text('changed source')
            def poll(self):return self.code
            def close(self):pass
        proxy=types.SimpleNamespace(ManagedProcess=Process,run_batch=supervision.run_batch)
        smoke=argparse.Namespace(stage='smoke',output=self.root/'smoke',contract=self.contract,smoke_receipt=None,smoke_sha256=None)
        with mock.patch.object(d,'module',return_value=proxy),mock.patch.object(d,'verify',return_value=(self.c,self.obj)):
            d.controller(smoke,self.c,self.obj)
            receipt=smoke.output/'complete.json'
            search=argparse.Namespace(stage='search',output=self.root/'search',contract=self.contract,smoke_receipt=receipt,smoke_sha256=d.sha(receipt))
            if mode=='objective':self.assertEqual(d.controller(search,self.c,self.obj),1)
            else:
                with self.assertRaises((RuntimeError,AssertionError)):d.controller(search,self.c,self.obj)
        return search,starts
    def test_fatal_search_stops_then_verifies_and_labels_saved_best(self):
        search,starts=self.run_interrupted_search('objective')
        complete=d.read(search.output/'complete.json')
        self.assertEqual(complete['status'],'interrupted_search_verified_best')
        self.assertEqual(sum(r['status']=='fatal' for r in complete['records']),1)
        self.assertEqual([r['stage'] for r in starts],['smoke','smoke','initial','initial','repeat','repeat'])
        export=d.read(search.output/'selected_export/export_receipt.json')
        self.assertEqual(export['status'],'verified_selected_export_from_interrupted_search')
        self.assertEqual(len(export['search_interruption']['fatal_records']),1)
    def test_setup_integrity_failure_never_exports_or_dispatches_repeats(self):
        search,starts=self.run_interrupted_search('setup')
        self.assertFalse((search.output/'selected_export').exists())
        self.assertTrue(d.read(search.output/'complete.json')['integrity_compromised'])
        self.assertFalse(any(r['stage']=='repeat' for r in starts))
    def test_changed_source_blocks_post_fatal_verification_and_export(self):
        search,starts=self.run_interrupted_search('source_tamper')
        self.assertFalse((search.output/'selected_export').exists())
        self.assertFalse(any(r['stage']=='repeat' for r in starts))

if __name__=='__main__':unittest.main()
