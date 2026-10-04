"""Mock exact-loop checks for distinct fixed shock workers and strict collection."""
import copy
import importlib
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
sys.path.insert(0,str(Path(__file__).resolve().parent))
import test_transition_driver as fixtures
r=importlib.import_module('experiments.birth_count_choice.transition_panel')


class FakeController:
    seen=[]
    def __init__(self,plan,adapter,output):
        self.plan=plan;self.runtime=adapter;self.out=Path(output);self.policy_calls=0
        self.last=None;self.count=0;self.deadline=0
    def prepare(self):self.policy_calls=62
    def run(self):
        self.prepare()
        result=self.evaluate(self.plan['fit_start_psi'])
        self.seed={'policy_calls':60}
        self.last['fertility']=[dict(period_tfr_topcode_adjusted=v) for v in result['payload']['models']]
        return dict(status='matched',fitted_parameter=dict(estimate=self.plan['fit_start_psi']))
    def evaluate(self,psi):
        self.seen.append(psi);self.count=1;self.policy_calls+=8
        target=self.plan['target_contract']['rows'][3]['target']
        models=[2.,1.9,1.8,target+psi-.2]
        self.last=dict(root_pass=True,replay_pass=True,accounting_valid=True,stationary_pass=True,
            market_maximum_residual=0.,fiscal_maximum_residual=0.,replay_maximum_gap=0.,
            stationary_renewal_gap=0.,terminal_passes=[False,False],
            horizon_comparison=dict(passed=True),state_horizon_comparison=dict(passed=True))
        return dict(certified=True,psi=psi,model=models[3],gap=models[3]-target,
            loss_contribution=(models[3]-target)**2,payload=dict(models=models))


class PanelTests(unittest.TestCase):
    def fixture(self,folder):
        plan,_=fixtures.DriverTests().fixture(folder)
        plan_path=Path(folder)/'plan.json';r.control.write(plan_path,plan)
        config=dict(schema='estate_a_transition_panel_v1',identity=plan['identity'],
            plan=r.driver.pin(plan_path),panel_source=r.driver.pin(r.__file__),
            guesses=[dict(index=i,psi=.1+.02*i) for i in range(12)])
        config_path=Path(folder)/'panel.json';r.control.write(config_path,config)
        return plan,r.driver.pin(config_path)
    def worker(self,plan,config,index,out):
        psi=.1+.02*index
        with patch.object(r.driver,'runtime_module',return_value=SimpleNamespace(
            CurrentEstateARuntime=SimpleNamespace(from_handoff=lambda *_:SimpleNamespace(total_native_calls=0)))),\
             patch.object(r.control,'NativeAdapter',return_value=SimpleNamespace()),\
             patch.object(r.control,'Controller',FakeController):
            return r.evaluate_candidate(plan,psi,index,out,panel_config=config)
    def test_distinct_starts_call_controller_run_and_preserve_plan(self):
        with tempfile.TemporaryDirectory() as folder:
            plan,config=self.fixture(folder);old=copy.deepcopy(plan);FakeController.seen=[]
            a=self.worker(plan,config,0,Path(folder)/'a');b=self.worker(plan,config,5,Path(folder)/'b')
            self.assertEqual(FakeController.seen,[.1,.2]);self.assertEqual(plan,old)
            self.assertEqual(len(a['fit_rows']),4);self.assertEqual([x['weight'] for x in a['fit_rows']],[0.,0.,0.,1.])
            self.assertEqual(a['actual_policy_calls'],70)
            effective=json.loads((Path(folder)/'a/effective_plan.json').read_text())
            self.assertEqual(effective,{**plan,'fit_start_psi':.1})
            self.assertTrue(a['scalar_optimizer_used']);self.assertFalse(b['production_ready'])
            collected=r.collect_candidates([Path(folder)/'a',Path(folder)/'b'],Path(folder)/'collected')
            self.assertEqual(collected['ranking'][0]['index'],5)
    def test_mixed_contract_and_failed_gate_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            plan,config=self.fixture(folder)
            self.worker(plan,config,0,Path(folder)/'a');self.worker(plan,config,1,Path(folder)/'b')
            path=Path(folder)/'b/candidate_receipt.json';rec=json.loads(path.read_text())
            rec['contract']['identity']['engine_sha256']='b'*64
            rec['contract_sha256']=r.fingerprint(rec['contract']);r.control.write(path,rec)
            with self.assertRaisesRegex(ValueError,'Mixed'):r.collect_candidates([Path(folder)/'a',Path(folder)/'b'],Path(folder)/'c')
            path=Path(folder)/'a/candidate_receipt.json';rec=json.loads(path.read_text());rec['gates']['root_pass']=False;r.control.write(path,rec)
            with self.assertRaisesRegex(ValueError,'gates failed'):r.collect_candidates([path],Path(folder)/'d')
    def test_config_mismatch_stops_before_native_callbacks(self):
        with tempfile.TemporaryDirectory() as folder:
            plan,config=self.fixture(folder)
            with patch.object(r.driver,'runtime_module') as module:
                with self.assertRaisesRegex(ValueError,'index/psi'):
                    r.evaluate_candidate(plan,.3,0,Path(folder)/'a',panel_config=config)
                module.assert_not_called()
            wrong=copy.deepcopy(config);wrong['sha256']='a'*64
            with self.assertRaisesRegex(ValueError,'changed pin'):r.validate_candidate(plan,.1,0,wrong)
    def test_rejected_proposals_never_ranked(self):
        with tempfile.TemporaryDirectory() as folder:
            plan,config=self.fixture(folder)
            self.worker(plan,config,0,Path(folder)/'a')
            path=Path(folder)/'a/candidate_receipt.json';rec=json.loads(path.read_text())
            rec.update(accepted=False,status='rejected');r.control.write(path,rec)
            result=r.collect_candidates([path],Path(folder)/'c')
            self.assertEqual(result['ranking'],[]);self.assertEqual(len(result['rejected_candidates']),1)


if __name__=='__main__':unittest.main()
