#!/usr/bin/env python3
"""Pure controller tests; intentionally do not run the model."""
import math
import json
import copy
import random
import tempfile
import unittest
from unittest import mock
from concurrent.futures import Future
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
import run_e5f_joint_nested_long_search as search

class ControllerTests(unittest.TestCase):
    def test_fixed_final_records_better_probe_without_selecting_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);plan=out/'plan.json';obj=object.__new__(search.Search)
            obj.c={'run_profile':'parallel32_fixed'};obj.best={'loss':10.};obj.ledger=[];obj.completed=0
            receipt={'status':'complete','plan_sha256':'sha','artifact_sha256':{},'case_id':1,'loss':5.,'elapsed_seconds':10.}
            summary={'best_candidate':{'transition_loss':5.},'panel_design':{'unit_vector':[.5]*11}}
            search.write_json(out/'case_receipt.json',receipt)
            with mock.patch.object(search.adapter,'validate_result',return_value=(summary,None,None)), \
                    mock.patch.object(search.adapter,'load_plan',return_value={}):
                obj._record_completed(plan,'sha',{'id':1,'label':'jacobian_0_minus'},out)
                self.assertEqual(obj.best,{'loss':10.});self.assertEqual(obj.ledger[0]['loss'],5.)
                obj._record_completed(plan,'sha',{'id':1,'label':'initial_1'},out)
                self.assertEqual(obj.best['loss'],5.)

    def test_combined_final_wave_requires_exact_repeats_and_discloses_better_probes(self):
        for missing_repeat in (False,True):
            with tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);obj=object.__new__(search.Search);obj.root=root
                obj.c={**search.RUN_PROFILES['parallel32_fixed'],'run_profile':'parallel32_fixed'}
                obj.completed=32;obj.best={'loss':10.,'unit_vector':[.5]*11,'summary':str(root/'anchor/summary.json')}
                for name in ('anchor','repeat1','repeat2'):
                    graphs=root/name/'standard_diagnostics';graphs.mkdir(parents=True)
                    for i in range(17):(graphs/f'{i}.png').write_bytes(bytes([i]))
                def batch(stage,vectors,labels,**kwargs):
                    self.assertEqual(stage,'final_joint_verification');self.assertEqual(len(vectors),24)
                    self.assertEqual(vectors[-2:],[[.5]*11]*2)
                    rows=[{'label':'jacobian_0_minus','loss':9.,'summary':str(root/'probe/summary.json')}]
                    rows += [{'label':f'selected_repeat_{i}','loss':10.,'summary':str(root/f'repeat{i}/summary.json')}
                             for i in range(1,2 if missing_repeat else 3)]
                    return rows
                with mock.patch.object(obj,'can_fit',return_value=True),mock.patch.object(obj,'batch',side_effect=batch) as run, \
                        mock.patch.object(obj,'write_jacobian'),mock.patch.object(search.adapter,'compare_reference') as compare, \
                        mock.patch.object(obj,'reports'),mock.patch.object(obj,'run_finalizer_if_pinned') as finalizer, \
                        mock.patch.object(obj,'summary'):
                    if missing_repeat:
                        with self.assertRaisesRegex(RuntimeError,'Two final repetitions'):obj.final_assessment()
                        finalizer.assert_not_called();self.assertFalse((root/'final_verification.json').exists())
                    else:
                        obj.final_assessment();self.assertEqual(compare.call_count,2);finalizer.assert_called_once()
                        proof=search.adapter.read_json(root/'final_verification.json')
                        self.assertEqual(proof['selected']['loss'],10.);self.assertEqual(proof['exact_repeats'],2)
                        diag=search.adapter.read_json(root/'fixed_selection_diagnostics.json')
                        self.assertEqual(diag['unselected_lower_loss_probes'][0]['loss'],9.)
                    run.assert_called_once()

    def test_smoke_only_contract_rejects_search_before_writing(self):
        with self.assertRaisesRegex(RuntimeError, 'verification only'):
            search.Search({'authorized_mode':'smoke'}, 'search')

    def test_smoke_probes_move_every_coordinate_at_bounds(self):
        center=[0.,1.,.5,.0001,.9999,.2,.3,.4,.6,.7,.8]
        lo,hi=search.smoke_probes(center)
        for j,x in enumerate(center):
            self.assertTrue(0 <= lo[j] <= 1 and 0 <= hi[j] <= 1)
            self.assertNotEqual(lo[j],x);self.assertNotEqual(hi[j],x)
            self.assertNotEqual(lo[j],hi[j])

    def test_parallel_population_keeps_v1_proposals_and_fits_one_wave(self):
        center=[.5]*11; domain=search.adapter.SEARCH_DOMAIN
        original=search.initial_population(center,domain,random.Random(20260906),'v1')
        parallel=search.initial_population(center,domain,random.Random(20260906),'parallel32')
        self.assertEqual(parallel,original)
        self.assertEqual(search.initial_population(center,domain,random.Random(20260906),'parallel32_fixed'),original)
        self.assertEqual(search.RUN_PROFILES['parallel32_fixed']['final_reserve_seconds'],3600+4200+1200)
        self.assertEqual(len(parallel),32)
        self.assertTrue(all(any(u[j]!=center[j] for u in parallel) for j in range(11)))
        obj=object.__new__(search.Search);obj.completed=4;obj.c=search.RUN_PROFILES['parallel32']
        obj.finish=20000.;obj.search_finish=10000.
        with mock.patch.object(search.time,'time',return_value=6400.):self.assertTrue(obj.can_fit(32))
        with mock.patch.object(search.time,'time',return_value=6400.1):self.assertFalse(obj.can_fit(32))
        self.assertEqual(obj.c['final_reserve_seconds'],7200+4200+1200)

    def test_parallel_reserve_requires_matching_complete_measured_policy_smoke(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);path=folder/'equilibrium_receipt.json';proof_path=folder/'policy_loop_verification.json'
            names=('baseline','supply-plus-20','dependent-child-ltv95','property-tax-2pct-no-rebate')
            branch={'status':'complete','dates':2,'source_summary_sha256':'selected',
                'gates':{'maximum_market_residual':1e-5,'maximum_mass_residual':1e-15}}
            receipt={'status':'complete','smoke':True,'failures':{},'scientific_bundle':'bundle',
                'target_fingerprint':'target','selected_summary_sha256':'selected','policy_workers':4,
                'elapsed_seconds':400.,'cases':{name:copy.deepcopy(branch) for name in names}}
            contract={'code_bundle_sha256':'bundle','target_fingerprint':'target','policy_workers':4}
            def write_evidence(value):
                path.write_text(json.dumps(value))  # Includes deliberately malformed nonfinite timing.
                search.write_json(proof_path,{'receipt':str(path),'sha256':search.digest(path)})
                contract['imported_smoke']={'root':str(folder),'policy_loop_verification_sha256':search.digest(proof_path)}
                contract['parallel_policy_timing']={'receipt':str(path),'receipt_sha256':search.digest(path),
                    'projected_full_policy_seconds':value['elapsed_seconds']*44/8}
            write_evidence(receipt)
            self.assertEqual(search.verify_parallel_policy_budget(contract),2200.)
            for mutation in ('slow','partial','mixed_bundle','serial','nan','nonpositive'):
                bad=copy.deepcopy(receipt)
                if mutation=='slow':bad['elapsed_seconds']=800.
                elif mutation=='partial':bad['status']='partial_policy_failures'
                elif mutation=='mixed_bundle':bad['scientific_bundle']='other'
                elif mutation=='serial':bad['policy_workers']=1
                elif mutation=='nan':bad['elapsed_seconds']=float('nan')
                else:bad['elapsed_seconds']=0.
                write_evidence(bad)
                with self.subTest(mutation=mutation),self.assertRaises(RuntimeError):search.verify_parallel_policy_budget(contract)
            write_evidence(receipt)
            search.write_json(proof_path,{'receipt':str(folder/'other.json'),'sha256':search.digest(path)})
            contract['imported_smoke']['policy_loop_verification_sha256']=search.digest(proof_path)
            with self.assertRaisesRegex(RuntimeError,'imported complete smoke'):search.verify_parallel_policy_budget(contract)

    def test_smoke_runs_four_required_histories_in_one_bounded_batch(self):
        with tempfile.TemporaryDirectory() as tmp:
            obj=object.__new__(search.Search);obj.root=Path(tmp);obj.c={}
            u=[.5]*11;obj.seed={'panel_design':{'unit_vector':u}}
            def batch(stage,vectors,labels,**kwargs):
                self.assertEqual(stage,'smoke_histories');self.assertTrue(kwargs['smoke'])
                self.assertEqual(labels,['anchor_1','anchor_2','all_minus','all_plus'])
                self.assertEqual(vectors[:2],[u,u])
                self.assertTrue(all(sum(a!=b for a,b in zip(v,u))==11 for v in vectors[2:]))
                rows=[]
                for i,(vector,label) in enumerate(zip(vectors,labels),1):
                    folder=obj.root/f'task_{i:03d}';graphs=folder/'standard_diagnostics';graphs.mkdir(parents=True)
                    for k in range(17):(graphs/f'graph_{k:02d}.png').write_bytes(bytes([k]))
                    rows.append({'id':i,'label':label,'unit_vector':vector,'summary':str(folder/'summary.json')})
                return list(reversed(rows))
            with mock.patch.object(obj,'batch',side_effect=batch) as launch, \
                    mock.patch.object(search.adapter,'compare_reference',return_value={'exact':True}) as exact, \
                    mock.patch.object(obj,'summary') as summary:
                obj.smoke()
            launch.assert_called_once();exact.assert_called_once();summary.assert_called_once_with('smoke_passed')
            proof=search.adapter.read_json(obj.root/'smoke_verification.json')
            self.assertEqual(len(proof['anchor_receipts']),2);self.assertEqual(len(proof['probe_receipts']),2)
            self.assertEqual(proof['exact_standard_pngs'],17)

    def test_smoke_rejects_an_incomplete_four_case_batch(self):
        obj=object.__new__(search.Search);obj.seed={'panel_design':{'unit_vector':[.5]*11}}
        with mock.patch.object(obj,'batch',return_value=[{'label':'anchor_1'}]):
            with self.assertRaisesRegex(RuntimeError,'Four completed smoke histories required'):obj.smoke()

    def test_three_timeouts_stop_new_search_but_drain_an_active_valid_case(self):
        with tempfile.TemporaryDirectory() as tmp:
            obj=object.__new__(search.Search);obj.root=Path(tmp)
            obj.c={**search.RUN_PROFILES['wide32'],'max_workers':4,
                'controller_sha256':'controller','adapter_sha256':'adapter','planner_sha256':'planner',
                'search_domain':search.adapter.SEARCH_DOMAIN,'source_sha256':search.adapter.SOURCE,
                'target_fingerprint':search.adapter.TARGET,
                'base_plan':{'schema':'e5f_joint_nested_overnight_v1','source_sha256':search.adapter.SOURCE,
                    'target_fingerprint':search.adapter.TARGET,'code_bundle_sha256':search.adapter.BUNDLE,
                    'search_domain':search.adapter.SEARCH_DOMAIN,'first_child_jump_upper':2.}}
            obj.seed={'old_psi_child':.3,'best_candidate':{'theta':{},'new_psi_child':.1}}
            obj.finish=obj.search_finish=1e20;obj.started=search.time.monotonic()
            obj.ledger=[];obj.rejects=[];obj.completed=4;obj.consecutive_timeouts=0;obj.search_stop_reason=None
            launched=[];pending=[]
            class Pool:
                def submit(self,fn,plan,sha,case):
                    launched.append(case['id']);future=Future();out=plan.parent/case['output']
                    if case['id']==4:pending.append((future,(case,out,'complete','')))
                    else:future.set_result((case,out,'timeout','case timeout'))
                    return future
                def shutdown(self,**kwargs):pass
            original_reject=obj._reject
            def reject(*args):
                original_reject(*args)
                if obj.search_stop_reason and pending:
                    future,result=pending.pop();future.set_result(result)
            def record(plan,sha,case,out):
                obj.ledger.append({'id':case['id'],'label':case['label'],'plan_sha256':sha,'loss':1.})
                obj.completed+=1;obj.consecutive_timeouts=0
            with mock.patch.object(search,'ThreadPoolExecutor',return_value=Pool()), \
                    mock.patch.object(obj,'_reject',side_effect=reject), \
                    mock.patch.object(obj,'_record_completed',side_effect=record), \
                    mock.patch.object(obj,'reports'),mock.patch.object(obj,'state'), \
                    mock.patch.object(search.planner,'collect') as collect:
                rows=obj.batch('initial_population',[[.5]*11 for _ in range(8)],[f'case_{i}' for i in range(8)])
            self.assertEqual(launched,[1,2,3,4])
            self.assertEqual([r['id'] for r in rows],[4])
            self.assertEqual(obj.completed,8)  # Four inherited cases plus four actually attempted.
            self.assertEqual(len(obj.rejects),3)
            self.assertTrue(all('loss' not in r for r in obj.rejects))
            unstarted=search.adapter.read_json(obj.root/'initial_population_unstarted.json')
            self.assertEqual([r['id'] for r in unstarted['cases']],[5,6,7,8])
            self.assertFalse(unstarted['counted_as_completed_or_rejected'])
            self.assertFalse(collect.call_args.kwargs['require_complete'])
            self.assertEqual(obj.consecutive_timeouts,0)
            self.assertEqual(obj.search_stop_reason['reason'],'three_consecutive_timeouts')
            self.assertFalse(obj.can_fit(1))  # A later success does not silently resume search.
            self.assertTrue(obj.can_fit(22,final=True))
            self.assertTrue(obj.can_fit(2,final=True))

    def test_timeout_stop_does_not_change_final_repeat_rejection(self):
        for stage,status,detail,smoke in (
                ('final_repeats','timeout','case timeout',True),
                ('initial_population','failed',{'type':'RuntimeError','error':'mass accounting gate failed'},False)):
            with self.subTest(stage=stage),tempfile.TemporaryDirectory() as tmp:
                obj=object.__new__(search.Search);obj.root=Path(tmp);obj.rejects=[]
                obj.c={'max_workers':1};obj.finish=1e20
                obj.search_stop_reason={'reason':'three_consecutive_timeouts'} if smoke else None
                case={'id':1,'output':'task_001'};plan=obj.root/'plan.json'
                with mock.patch.object(obj,'can_fit',return_value=True), \
                        mock.patch.object(obj,'new_plan',return_value=(plan,'sha')), \
                        mock.patch.object(search.adapter,'load_plan',return_value={'cases':[case]}), \
                        mock.patch.object(obj,'_child',return_value=(case,obj.root/'task_001',status,detail)), \
                        mock.patch.object(obj,'_reject') as reject,mock.patch.object(obj,'stop_active') as stop:
                    with self.assertRaisesRegex(RuntimeError,'Required verification case rejected|Fatal candidate failure'):
                        obj.batch(stage,[[.5]*11],['case'],smoke=smoke)
                stop.assert_called_once()
                self.assertEqual(reject.call_count,1 if smoke else 0)

    def test_wide_reserve_reduces_search_time_before_hard_cutoff(self):
        with tempfile.TemporaryDirectory() as tmp:
            seed=Path(tmp)/'seed.json';search.write_json(seed,{})
            contract={**search.RUN_PROFILES['wide32'],'output_root':tmp,'seed_center':str(seed),
                      'absolute_finish_epoch':31000.}
            with mock.patch.object(search.time,'time',return_value=1000.):
                obj=search.Search(contract,'test')
            self.assertEqual(obj.finish,31000.)
            self.assertEqual(obj.search_finish,14800.)
            self.assertEqual(obj.finish-obj.search_finish,16200.)

    def test_final_repeats_do_not_hide_skipped_policy_finalization(self):
        with tempfile.TemporaryDirectory() as tmp:
            obj=object.__new__(search.Search);obj.root=Path(tmp);obj.c={'finalizer_driver':'pinned.py'}
            self.assertEqual(obj.completion_status(),'best_valid_without_final_repeats')
            search.write_json(obj.root/'final_verification.json',{'status':'pass'})
            self.assertEqual(obj.completion_status(),'calibration_verified_policy_incomplete')
            search.write_json(obj.root/'finalizer_status.json',{'status':'skipped_insufficient_time'})
            self.assertEqual(obj.completion_status(),'calibration_verified_policy_incomplete')
            search.write_json(obj.root/'finalizer_status.json',{'status':'complete'})
            self.assertEqual(obj.completion_status(),'complete_verified')

    def test_wide_batch_uses_unchanged_adapter_and_preserves_population_slots(self):
        with tempfile.TemporaryDirectory() as tmp:
            obj=object.__new__(search.Search); obj.root=Path(tmp)
            obj.c={**search.RUN_PROFILES['wide32'], 'controller_sha256':'controller',
                'adapter_sha256':'adapter','planner_sha256':'planner', 'search_domain':search.adapter.SEARCH_DOMAIN,
                'base_plan':{'schema':'e5f_joint_nested_overnight_v1','source_sha256':search.adapter.SOURCE,
                    'target_fingerprint':search.adapter.TARGET,'code_bundle_sha256':search.adapter.BUNDLE,
                    'search_domain':search.adapter.SEARCH_DOMAIN,'first_child_jump_upper':2.}}
            obj.seed={'old_psi_child':.3,'best_candidate':{'theta':{},'new_psi_child':.1}}
            obj.finish=obj.search_finish=1e20; obj.ledger=[]; obj.rejects=[]
            def child(plan,sha,case):
                # Real adapter validation rejects the original oversized-plan bug.
                self.assertLessEqual(len(search.adapter.load_plan(plan,sha)['cases']),40)
                return case,plan.parent/case['output'],'complete',''
            def record(plan,sha,case,out):
                obj.ledger.append({'id':case['id'],'label':case['label'],'plan_sha256':sha})
            with mock.patch.object(obj,'can_fit',return_value=True), mock.patch.object(obj,'_child',side_effect=child), \
                    mock.patch.object(obj,'_record_completed',side_effect=record), mock.patch.object(obj,'reports'), \
                    mock.patch.object(obj,'state'), mock.patch.object(search.planner,'collect') as collect:
                rows=obj.batch('initial_population',[[.5]*11 for _ in range(64)],[f'slot_{i}' for i in range(64)])
            self.assertEqual(sorted(r['id'] for r in rows),list(range(1,65)))
            self.assertTrue(all(r['label']==f"slot_{r['id']-1}" for r in rows))
            self.assertEqual(collect.call_count,2)
            self.assertEqual(len({r['plan_sha256'] for r in rows}),2)

    def test_policy_receipt_rejects_partial_mixed_or_incomplete_paths(self):
        names=('baseline','supply-plus-20','dependent-child-ltv95','property-tax-2pct-no-rebate')
        branch={'status':'complete','dates':2,'source_summary_sha256':'selected',
            'gates':{'maximum_market_residual':1e-5,'maximum_mass_residual':1e-15}}
        receipt={'status':'complete','smoke':True,'failures':{},'scientific_bundle':'bundle',
            'target_fingerprint':'target','selected_summary_sha256':'selected',
            'cases':{name:copy.deepcopy(branch) for name in names}}
        contract={'code_bundle_sha256':'bundle','target_fingerprint':'target'}
        search.validate_policy_receipt(receipt,contract,smoke=True,selected_hashes={'selected'})
        for mutation in ('partial','wrong_source','short','nan','concurrency'):
            bad=copy.deepcopy(receipt)
            if mutation=='partial':bad['status']='partial_policy_failures'
            elif mutation=='wrong_source':bad['selected_summary_sha256']='other'
            elif mutation=='short':bad['cases']['baseline']['dates']=1
            elif mutation=='concurrency':bad['policy_workers']=4
            else:bad['cases']['baseline']['gates']['maximum_market_residual']=float('nan')
            with self.assertRaises(RuntimeError):
                search.validate_policy_receipt(bad,contract,smoke=True,selected_hashes={'selected'})

    def test_failure_classification(self):
        self.assertEqual(search.classify_failure("Housing market did not clear"), "market_nonconvergence")
        self.assertIsNone(search.classify_failure("mass accounting gate failed"))
        self.assertIsNone(search.classify_failure("fertility accounting tolerance failed"))
        self.assertEqual(search.classify_failure("Old-steady-state fertility normalization missed tolerance: x"), "fertility_normalization")
        self.assertEqual(search.classify_failure("x", "InfeasibleThetaError"), "infeasible_theta")

    def test_population_changes_all_coordinates_and_is_bounded(self):
        domain=[{"name":f"x{i}","lower":.1,"upper":10.,"transform":"log"} for i in range(11)]
        domain[1]["name"]="tenure_choice_kappa"; domain[2]["name"]="joint_nest_lambda"; domain[9]["name"]="hbar_first_child_jump"
        center=[.5]*11; pop=search.initial_population(center,domain,random.Random(20260906))
        self.assertEqual(len(pop),32)
        self.assertTrue(all(len(u)==11 and all(0<=x<=1 for x in u) for u in pop))
        self.assertTrue(all(any(x != .5 for x in u) for u in pop[1:]))

    def test_wide32_population_is_deterministic_and_covers_unique_grid(self):
        domain=[{"name":f"x{i}","lower":.001,"upper":10.,"transform":"log"} for i in range(11)]
        domain[1]["name"]="tenure_choice_kappa"; domain[2]["name"]="joint_nest_lambda"; domain[9]["name"]="hbar_first_child_jump"
        domain[10]={"name":"psi_child_change_2023","lower":-1.5,"upper":.2,"transform":"asinh"}
        center=[.5]*11
        first=search.initial_population(center,domain,random.Random(20260906),"wide32")
        second=search.initial_population(center,domain,random.Random(20260906),"wide32")
        self.assertEqual(first,second)
        self.assertEqual(len(first),64); self.assertEqual(first[0],center)
        self.assertTrue(all(len(u)==11 and all(0<=x<=1 for x in u) for u in first))
        grid={(round(search.transform(u[1],domain[1]),8),round(search.transform(u[2],domain[2]),8)) for u in first[1:48]}
        expected={(k,l) for k in (.01,.03,.1,.3,1.,2.,4.,8.) for l in (.02,.05,.2,.5,.8,1.) if (k,l)!=(2.,.8)}
        self.assertEqual(grid,{(round(k,8),round(l,8)) for k,l in expected})
        self.assertNotIn((2.,.8),grid)
        self.assertTrue(all(abs(u[j]-center[j])<=.08 for u in first[1:48] for j in range(11) if j not in (1,2)))
        paired=[u for u in first[1:48] if any(math.isclose(search.transform(u[2],domain[2]),l) for l in (.02,.2))]
        self.assertEqual(len(paired),16)
        self.assertEqual(len({tuple(u) for u in first}),64)
        for original,proposal in zip(paired,first[48:]):
            self.assertEqual(original[:10],proposal[:10])
            old_delta=search.transform(original[10],domain[10]); new_delta=search.transform(proposal[10],domain[10])
            self.assertLess(abs(new_delta),abs(old_delta))
            self.assertEqual(math.copysign(1,new_delta),math.copysign(1,old_delta))

    def test_wide_proposals_cover_small_changes_without_removing_original_grid(self):
        domain=search.adapter.SEARCH_DOMAIN
        center=[.5]*11
        center[1]=math.log(2./.005)/math.log(10./.005)
        center[2]=math.log(.8/.02)/math.log(1./.02)
        center[10]=(math.asinh(-.32871390689556357)-math.asinh(-1.5))/(math.asinh(.2)-math.asinh(-1.5))
        population=search.initial_population(center,domain,random.Random(20260906),"wide32")
        original_changes=[search.transform(u[10],domain[10]) for u in population[1:48]]
        paired_changes=[search.transform(u[10],domain[10]) for u in population[48:]]
        self.assertTrue(all(delta < -.2 for delta in original_changes))
        self.assertTrue(any(-.0001 < delta < 0 for delta in paired_changes))
        self.assertTrue(all(-.23 < delta < 0 for delta in paired_changes))
        # A subsequent DE proposal still varies the preference coordinate;
        # it is not constrained to the starting-point scaling relation.
        trials=search.de_trials(population,list(range(64)),random.Random(17))
        self.assertTrue(any(u[10] != p[10] for u,p in zip(trials,population)))

    def test_profiles_and_wide_budget_reserve_final_verification(self):
        self.assertEqual(search.RUN_PROFILES["wide32"]["max_workers"],32)
        self.assertEqual(search.RUN_PROFILES["wide32"]["population_size"],64)
        self.assertEqual(search.RUN_PROFILES["wide32"]["max_histories"],640)
        obj=object.__new__(search.Search); obj.completed=580
        obj.c=search.RUN_PROFILES["wide32"]; obj.finish=1e20; obj.search_finish=1e20
        self.assertFalse(obj.can_fit(64))  # 580 + 64 would crowd out 24 required final histories.
        obj.completed=606
        self.assertFalse(obj.can_fit(2*search.N))
        obj.completed=616
        self.assertTrue(obj.can_fit(2*search.N,final=True))
        self.assertTrue(obj.can_fit(2,final=True))

    def test_stale_imported_smoke_contract_is_rejected_before_receipts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); original=root/"original_contract.json"
            search.write_json(original,{"source_sha256":"stale","code_bundle_sha256":"bundle","target_fingerprint":"target","search_domain":[]})
            obj=object.__new__(search.Search); obj.c={"source_sha256":"current","code_bundle_sha256":"bundle",
                "target_fingerprint":"target","search_domain":[],"imported_smoke":{"root":str(root),
                "original_contract":str(original),"original_contract_sha256":search.digest(original),
                "smoke_verification_sha256":"unused"}}
            with self.assertRaisesRegex(RuntimeError,"original contract has mixed source_sha256"):
                obj.require_smoke()

    def test_final_repeat_plan_copies_original_center_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); source=root/"source"; source.mkdir(); original_center=source/"center.json"
            original_center.write_bytes(b'{"original":"exact bytes"}\n')
            origin={"plan":str(source/"plan.json"),"plan_sha256":"plan-sha","id":7}
            source_plan={"cases":[{"id":7,"center":"center.json","panel_task_id":3,"panel_size":5,
                "panel_design":"mixed","panel_seed":17,"radius":.05}]}
            domain=[{"name":f"x{i}","lower":.1,"upper":10.,"transform":"log"} for i in range(11)]
            obj=object.__new__(search.Search); obj.root=root/"run"; obj.root.mkdir()
            obj.c={"base_plan":{},"controller_sha256":"controller","adapter_sha256":"adapter","planner_sha256":"planner",
                   "search_domain":domain,"max_search_seconds":1e12}; obj.finish=1e12; obj.search_finish=1e12
            obj.seed={"old_psi_child":1.,"best_candidate":{"theta":{},"new_psi_child":2.}}
            obj.repeat_origin=origin
            with mock.patch.object(search.adapter,"load_plan",return_value=source_plan):
                plan_path,_=obj.new_plan("final_repeats",[[.1]*11,[.1]*11],["selected_repeat_1","selected_repeat_2"])
            plan=search.adapter.read_json(plan_path)
            self.assertEqual((plan_path.parent/plan["cases"][0]["center"]).read_bytes(),original_center.read_bytes())
            self.assertEqual(plan["cases"][0]["panel_task_id"],3)
            self.assertEqual(plan["cases"][1]["panel_seed"],17)
            with mock.patch.object(search.adapter,'load_plan',return_value=source_plan):
                combined,_=obj.new_plan('final_joint_verification',[[.2]*11,[.1]*11,[.1]*11],
                    ['jacobian_0_minus','selected_repeat_1','selected_repeat_2'])
            cp=search.adapter.read_json(combined)
            self.assertNotEqual((combined.parent/cp['cases'][0]['center']).read_bytes(),original_center.read_bytes())
            for case in cp['cases'][1:]:
                self.assertEqual((combined.parent/case['center']).read_bytes(),original_center.read_bytes())
                self.assertEqual(case['panel_task_id'],3)

    def test_budget_and_valid_incumbent_logic(self):
        obj=object.__new__(search.Search); obj.completed=358; obj.c={"max_histories":360,"max_workers":12,"case_timeout_seconds":3600}; obj.finish=1e20; obj.search_finish=1e20
        self.assertTrue(obj.can_fit(2, final=True)); self.assertFalse(obj.can_fit(3, final=True)); self.assertFalse(obj.can_fit(2))
        valid=[{"loss":3.}, {"loss":2.}]
        rejected={"rejection_type":"market_nonconvergence", "error":"no calibrated loss"}
        self.assertEqual(search.best_completed(valid)["loss"],2.)
        self.assertNotIn("loss", rejected)  # A rejected proposal cannot become an incumbent.

if __name__ == "__main__": unittest.main()
