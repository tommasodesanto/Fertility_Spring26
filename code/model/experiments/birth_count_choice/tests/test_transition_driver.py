"""Zero-solve plan, identity and unchanged-controller delegation checks."""
import copy
import importlib
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
r = importlib.import_module('experiments.birth_count_choice.transition')
CASE = r.ROOT/'output/model/experiments/birth_count_choice/estate_a_v1/single/cases/20261003T212605039706Z_a739edc3'


class DriverTests(unittest.TestCase):
    def fixture(self, folder):
        record = r.build_handoff(CASE,[Path(r.__file__),r.HERE/'transition_runtime.py',
            r.HERE/'model/engine/household.py'])
        path = Path(folder)/'handoff.json'; r.controller.write(path,record)
        identity = {k:'a'*64 for k in r.controller.IDENTITY_KEYS-{'source_pins'}}
        identity.update(reference_sha256=r.controller.sha(path),source_pins=record['source_pins'])
        rt = SimpleNamespace(identity=lambda:identity,P=SimpleNamespace(psi_child=.2))
        budget = {k:30 for k in r.BUDGET_KEYS}; budget['maximum_policy_calls']=500
        plan = r.build_plan(r.pin(path),rt,mode='diagnostic',horizons=[24,32],budget=budget,
            fit_max_evaluations=5,path_max_evaluations=6,endpoint_max_evaluations=2,
            fiscal_relaxation_authorized=True,execution_enabled=True)
        return plan,rt

    def test_full_retained_targets_gates_and_explicit_identity(self):
        with tempfile.TemporaryDirectory() as folder:
            p,_=self.fixture(folder)
            self.assertEqual(r.preflight(p)['native_calls'],0)
            self.assertEqual(p['target_contract']['rows'][3]['target'],1.64575)
            self.assertEqual(len(p['target_contract']['rows']),4)
            self.assertEqual(p['gates'],r.controller.GATES)
            self.assertEqual(p['source_files']['controller'],r.pin(r.SHARED/'one_shock_floor.py'))
            self.assertEqual(p['psi_bound_ratios'],[.01,2.])
            self.assertFalse(p['policy_contract_closed'])
            self.assertEqual(len(p['standard_plot_names']),17)

    def test_historical_preparation_and_changed_source_are_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            p,_=self.fixture(folder)
            for key in ('prepared_native_inputs','prepared_consumer_compatibility','diagnostic_measurement_reuse'):
                changed=copy.deepcopy(p); changed[key]={'path':'old','sha256':'a'*64}
                with self.assertRaisesRegex(ValueError,'historical reuse'):r.preflight(changed)
            changed=copy.deepcopy(p);changed['source_files']['driver']['sha256']='a'*64
            with self.assertRaisesRegex(ValueError,'driver changed'):r.preflight(changed)
            changed=copy.deepcopy(p);changed['identity']['source_pins']={'old':'a'*64}
            with self.assertRaisesRegex(ValueError,'source pins'):r.preflight(changed)

    def test_production_horizons_and_fiscal_flag_remain_required(self):
        with tempfile.TemporaryDirectory() as folder:
            p,_=self.fixture(folder)
            changed=copy.deepcopy(p);changed['mode']='production'
            with self.assertRaisesRegex(ValueError,'104/128'):r.preflight(changed)
            changed['horizons']=[104,128]
            self.assertEqual(r.preflight(changed)['mode'],'production')
            changed['fiscal_relaxation_authorized']=False
            with self.assertRaisesRegex(ValueError,'fiscal'):r.preflight(changed)

    def test_execute_calls_literal_retained_adapter_and_controller(self):
        with tempfile.TemporaryDirectory() as folder:
            p,rt=self.fixture(folder)
            module=SimpleNamespace(CurrentEstateARuntime=SimpleNamespace(from_handoff=lambda *_:rt))
            with patch.object(r,'runtime_module',return_value=module), \
                 patch.object(r.controller,'NativeAdapter') as adapter, \
                 patch.object(r.controller,'Controller') as runner:
                runner.return_value.run.return_value={'status':'fake'}
                self.assertEqual(r.execute(p,folder),{'status':'fake'})
                adapter.assert_called_once_with(rt,p)
                runner.assert_called_once_with(p,adapter.return_value,folder)
                runner.return_value.run.assert_called_once_with()

    def test_missing_budget_and_execution_flag_fail_before_callbacks(self):
        with tempfile.TemporaryDirectory() as folder:
            p,rt=self.fixture(folder)
            p['execution_enabled']=False
            with patch.object(r,'runtime_module') as module:
                with self.assertRaisesRegex(ValueError,'enable execution'):r.execute(p,folder)
                module.assert_not_called()
            kwargs=dict(mode='diagnostic',horizons=[24,32],budget={},fit_max_evaluations=5,
                path_max_evaluations=6,endpoint_max_evaluations=2,fiscal_relaxation_authorized=True)
            with self.assertRaisesRegex(ValueError,'budget'):r.build_plan(p['handoff'],rt,**kwargs)


if __name__=='__main__':unittest.main()
