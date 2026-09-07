#!/usr/bin/env python3
"""Pure controller tests; intentionally do not run the model."""
import math
import copy
import random
import tempfile
import unittest
from unittest import mock
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
import run_e5f_joint_nested_long_search as search

class ControllerTests(unittest.TestCase):
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
        for mutation in ('partial','wrong_source','short','nan'):
            bad=copy.deepcopy(receipt)
            if mutation=='partial':bad['status']='partial_policy_failures'
            elif mutation=='wrong_source':bad['selected_summary_sha256']='other'
            elif mutation=='short':bad['cases']['baseline']['dates']=1
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
        self.assertTrue(all(abs(u[j]-center[j])<=.12 for u in first[48:] for j in range(11)))

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

    def test_budget_and_valid_incumbent_logic(self):
        obj=object.__new__(search.Search); obj.completed=358; obj.c={"max_histories":360,"max_workers":12,"case_timeout_seconds":3600}; obj.finish=1e20; obj.search_finish=1e20
        self.assertTrue(obj.can_fit(2, final=True)); self.assertFalse(obj.can_fit(3, final=True)); self.assertFalse(obj.can_fit(2))
        valid=[{"loss":3.}, {"loss":2.}]
        rejected={"rejection_type":"market_nonconvergence", "error":"no calibrated loss"}
        self.assertEqual(search.best_completed(valid)["loss"],2.)
        self.assertNotIn("loss", rejected)  # A rejected proposal cannot become an incumbent.

if __name__ == "__main__": unittest.main()
