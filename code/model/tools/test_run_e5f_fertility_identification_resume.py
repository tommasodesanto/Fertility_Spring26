"""Pure scheduler fixtures. Run on Torch; no model evaluation is invoked."""
import copy
import unittest
from types import SimpleNamespace
from run_e5f_fertility_identification_resume import ReplayScheduler, inherit_clock


class ResumeTests(unittest.TestCase):
    def fixture(self, status='success'):
        req = dict(id='search_0040_primary', lane='primary', point={'x':1.},
                   design='de:0:0', context={'stage':'search'})
        payload = dict(req, deadline_epoch=100)
        row = dict(deadline=100, returncode=0 if status=='success' else -9,
                   status=status, case=req['id'])
        replay = {req['id']:dict(row=row,request=payload,request_path='/old/request.json')}
        calls = []; finished = []
        def run_batch(requests, **kwargs):
            calls.extend(copy.deepcopy(requests))
            return dict(complete=True,results=[{'status':'success'} for _ in requests])
        adapter = ReplayScheduler(replay,{req['id']:(status,{},None)},run_batch)
        def finish(request,proc,code):
            got = adapter.classify(None,None,None,request['payload'],proc,code,
                                   lambda *args: self.fail('Replay reached real classifier'))
            result = dict(case=request['id'],status=got[0]); finished.append(result)
            return result
        return req, adapter, calls, finished, finish

    def test_completed_success_reused_new_case_only_dispatched(self):
        req,a,calls,finished,finish=self.fixture()
        new=dict(id='search_0041_early10',lane='early10',point={'x':2.},
                 design='de:0:0',context={'stage':'search'})
        result=a.run_batch([req,new],finish=finish)
        self.assertEqual([r['id'] for r in calls],['search_0041_early10'])
        self.assertEqual(calls[0]['context']['stage'],'de')
        self.assertEqual(len(result['results']),2)
        self.assertEqual(len(a.used),1)
        self.assertEqual(a.replay[req['id']]['request']['context']['stage'],'search')

    def test_timeout_never_retried(self):
        req,a,calls,finished,finish=self.fixture('censored_timeout')
        a.run_batch([req],finish=finish)
        self.assertFalse(calls)
        self.assertEqual(finished[0]['status'],'censored_timeout')

    def test_duplicate_proposal_rejected(self):
        req,a,calls,finished,finish=self.fixture()
        a.run_batch([req],finish=finish)
        with self.assertRaises(AssertionError):a.run_batch([req],finish=finish)

    def test_point_mismatch_rejected(self):
        req,a,calls,finished,finish=self.fixture()
        req=dict(req,point={'x':99})
        with self.assertRaises(AssertionError):a.run_batch([req],finish=finish)
        self.assertFalse(calls)

    def test_failure_overlay_is_rejection_and_original_preserved(self):
        req,a,calls,finished,finish=self.fixture('fatal')
        original=copy.deepcopy(a.replay)
        a.classifications[req['id']]=('inadmissible',{}, {'original_status':'fatal'})
        a.run_batch([req],finish=finish)
        self.assertEqual(finished[0]['status'],'inadmissible')
        self.assertEqual(a.replay,original)
        self.assertFalse(calls)

    def test_alias_jacobian_only_new_context(self):
        req,a,calls,finished,finish=self.fixture()
        a.run_batch([dict(id='new_probe',context={'stage':'jacobian'})],finish=finish)
        self.assertEqual(calls[0]['context']['stage'],'initial')

    def test_new_classifier_uses_unchanged_core(self):
        _,a,_,_,_=self.fixture()
        sentinel=object()
        self.assertIs(a.classify(None,None,None,{},SimpleNamespace(),0,
                                lambda *args:sentinel),sentinel)

    def test_deadline_inherited_no_restart(self):
        clock=dict(start=1790608994.4615145,search_cutoff=1790626994.4615145,
                   repeat_cutoff=1790629994.4615145,end=1790630594.4615145)
        self.assertEqual(inherit_clock(clock,clock,clock['start']+100),clock)
        with self.assertRaises(AssertionError):inherit_clock(clock,clock,clock['search_cutoff'])
        changed=dict(clock,end=clock['end']+1)
        with self.assertRaises(AssertionError):inherit_clock(changed,clock,clock['start']+100)



class AuthenticationTests(unittest.TestCase):
    def fixture(self, root):
        from pathlib import Path
        from run_e5f_fertility_identification_resume import core
        old=root/'old';old.mkdir()
        clock=dict(start=1790608994.4615145,search_cutoff=1790626994.4615145,
                   repeat_cutoff=1790629994.4615145,end=1790630594.4615145)
        smoke=root/'smoke.json';core.write(smoke,dict(clock=clock))
        approval=root/'approval.json';core.write(approval,dict(smoke_receipt=dict(path=str(smoke),sha256=core.sha(smoke))))
        rows=[];pins=[]
        for i in range(41):
            name=f'jacobian_{i:04d}_primary' if i<40 else 'search_0040_primary'
            request=old/(name+'.request.json')
            req=dict(id=name,lane='primary',point={'x':1},design=str(i),contract_sha256='contract',scientific_candidate_id='identity',context={'stage':'jacobian' if i<40 else 'search'})
            core.write(request,req)
            row=dict(case=name,lane='primary',point=req['point'],design=str(i),request_path=str(request),status='success',error=None,case_path=str(old/name/'case'),loss=1)
            if i==40:
                failure=dict(context=req['context'],error_type='NonpositiveNormalizedBenefit',phase='objective',status='fatal',classifier_error='candidate identity/stage required')
                path=old/name/'failure.json';core.write(path,failure)
                pins.append(dict(path=str(path),sha256=core.sha(path)))
                row.update(status='fatal',error=failure)
            rows.append(row);pins.append(dict(path=str(request),sha256=core.sha(request)))
        complete=old/'complete.json';core.write(complete,dict(status='incomplete_or_fatal_stop',contract_sha256='contract',clock=clock,records=rows))
        pins.extend(dict(path=str(p),sha256=core.sha(p)) for p in (complete,approval))
        manifest=dict(status='lead_reviewed_resume',contract_sha256='contract',old_job_terminal=True,pins=pins,old_complete=str(complete),original_approval=str(approval),rejection_overlays={'search_0040_primary':dict(status='inadmissible',reason='verified_nonpositive_benefit_stage_label_only')})
        return manifest,dict(budget={'max_objective_cases':538}),clock

    def call(self,manifest,c,clock):
        from unittest.mock import patch
        from run_e5f_fertility_identification_resume import authenticate,core
        with patch.object(core,'identity',return_value='identity'),patch.object(core,'validate',return_value={'loss':1}):
            return authenticate(c,{},manifest,'contract',clock['start']+100)

    def test_authentication_and_fail_closed(self):
        import tempfile
        from pathlib import Path
        from run_e5f_fertility_identification_resume import core
        with tempfile.TemporaryDirectory() as temp:
            manifest,c,clock=self.fixture(Path(temp))
            rows,statuses=self.call(manifest,c,clock)
            self.assertEqual(len(rows),41)
            self.assertEqual(statuses['search_0040_primary'][0],'inadmissible')
            bad=copy.deepcopy(manifest);bad['pins'][0]['sha256']='wrong'
            with self.assertRaises(AssertionError):self.call(bad,c,clock)
            bad=copy.deepcopy(manifest);bad['rejection_overlays']={}
            with self.assertRaises(AssertionError):self.call(bad,c,clock)
            extra=Path(temp)/'old/omitted.request.json';core.write(extra,{})
            with self.assertRaises(AssertionError):self.call(manifest,c,clock)
            extra.unlink()
            failure_path=Path(temp)/'old/search_0040_primary/failure.json'
            failure=core.read(failure_path);failure['classifier_error']='some other classifier error';core.write(failure_path,failure)
            for pin in manifest['pins']:
                if pin['path']==str(failure_path):pin['sha256']=core.sha(failure_path)
            with self.assertRaises(AssertionError):self.call(manifest,c,clock)


class OriginalLoopTests(unittest.TestCase):
    def test_original_loop_reconstructs_all_generations_without_duplicate_solves(self):
        import tempfile
        from pathlib import Path
        from unittest.mock import patch
        from contextlib import ExitStack
        from run_e5f_fertility_identification_resume import original,core
        with tempfile.TemporaryDirectory() as temp, ExitStack() as stack:
            root=Path(temp); case=root/'case';(case/'standard_diagnostics').mkdir(parents=True)
            clock=dict(start=1790608994.4615145,search_cutoff=1790626994.4615145,repeat_cutoff=1790629994.4615145,end=1790630594.4615145)
            c=dict(budget=dict(workers=24,population=16,generations=4,max_objective_cases=538,total_seconds=21600,objective_cap_seconds=1800),files={'driver':{'sha256':'pin'},'recovery_search':{'path':'helper'}},normalization={},early_moment='early',anchor={'case_path':str(case)},lanes={k:dict(objective={'sha256':'pin'},fixed={}) for k in original.LANES},standard_diagnostic_names=[])
            objs={k:{'target_rows':[{'restriction_id':'early','actual_weight':1}]} for k in original.LANES}
            def point(lane,slot,generation=0):return {'x':1+original.LANES.index(lane)*100+slot+generation*1000}
            probes=[dict(lane='primary',point={'x':i+1},design='jacobian:'+str(i)) for i in range(40)]
            initial=[dict(lane=k,point=point(k,s),design=f'de:0:{s}') for s in range(16) for k in original.LANES]
            replay={};classifications={}
            for i,item in enumerate(probes+initial[:12]):
                stage='jacobian' if i<40 else 'search';name=f'{stage}_{i:04d}_{item["lane"]}'
                ctx=dict(candidate_id=name,stage=stage,contract_sha256='pin',source_sha256='pin',target_sha256='pin',point_sha256=core.canon(item['point']))
                payload=dict(item,id=name,context=ctx)
                status='censored_timeout' if i==41 else 'success'
                replay[name]=dict(row={'deadline':clock['search_cutoff'],'returncode':0},request=payload,request_path='old/'+name)
                classifications[name]=(status,dict(loss=item['point']['x'],primary_rescore=item['point']['x'],case_path=str(case)) if status=='success' else {},None)
            new_calls=[];parent_checks=[]
            def managed(*args):return SimpleNamespace(deadline=clock['search_cutoff'])
            def schedule(reqs,**kw):
                out=[]
                for req in reqs:
                    new_calls.append(req['id']);proc=kw['launch'](req,kw['deadline']);out.append(kw['finish'](req,proc,0))
                return dict(complete=True,results=out)
            adapter=ReplayScheduler(replay,classifications,schedule)
            helper=SimpleNamespace(run_batch=adapter.run_batch,ManagedProcess=managed)
            def classify(folder,c,objs,req,proc,code):
                return adapter.classify(folder,c,objs,req,proc,code,lambda *args:('success',dict(loss=req['point']['x'],primary_rescore=req['point']['x'],case_path=str(case)),None))
            prior=[]
            for lane in original.LANES[:3]:
                for i in range(2):
                    path=root/f'smoke_{lane}_{i}.json';core.write(path,{})
                    prior.append(dict(lane=lane,request_path=str(path),case_path=str(case)))
            smoke=root/'smoke.json';core.write(smoke,dict(status='exact_loop_smoke_passed',contract_sha256='pin',clock=clock,records=prior))
            approval=root/'approval.json';core.write(approval,dict(status='approved_search',contract_sha256='pin',smoke_receipt={'path':str(smoke),'sha256':'pin'}))
            def trials(c,obj,lane,population,generation):
                self.assertEqual(len(population),16)
                self.assertEqual({r['lane'] for r in population},{lane})
                parent_checks.append((lane,generation))
                return [point(lane,s,generation) for s in range(16)]
            patches=[(original,'sha',lambda p:'pin'),(core,'sha',lambda p:'pin'),(core,'module',lambda *args:helper),(core,'validate',lambda *args:{}),(core,'compare_anchor',lambda *args:None),(core,'compare_tables',lambda *args:None),(core,'classify',classify),(core,'table',lambda *args:{'early':{'model':'0.5','gap':'-0.3'}}),(core,'identity',lambda *args:'identity'),(original,'verify',lambda *args:(c,objs)),(original,'probe_points',lambda c:probes),(original,'initial_population',lambda c,obj,lane:[point(lane,s) for s in range(16)]),(original,'de_trials',trials),(original,'jacobian_export',lambda *args:None),(original.time,'time',lambda:clock['start']+100)]
            for obj,name,value in patches:stack.enter_context(patch.object(obj,name,value))
            a=SimpleNamespace(stage='run',contract=root/'contract',approval=approval,approval_sha256='pin',output=root/'new')
            original.controller(a,c,objs)
            complete=core.read(a.output/'complete.json')
            self.assertEqual(complete['status'],'bounded_experiment_complete')
            self.assertEqual(len(complete['records'])+6,538)
            self.assertEqual(len(new_calls),480-12+12)
            self.assertEqual(len(set(new_calls)),len(new_calls))
            self.assertEqual(adapter.used,set(replay))
            self.assertEqual(len(parent_checks),24)
            self.assertEqual(complete['clock'],clock)

if __name__=='__main__':unittest.main()
