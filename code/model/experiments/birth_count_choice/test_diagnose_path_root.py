"""No native imports, constructors, or solves."""
import importlib.util
import json
from pathlib import Path
import tempfile
import time
from types import SimpleNamespace
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('diagnose_path_root',HERE/'diagnose_path_root.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


class RootDiagnosticTests(unittest.TestCase):
    def test_original_controls_with_twelve_evaluations(self):
        package = Path(__file__).resolve().parents[4] / 'output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/execution_smoke_v5'
        plan = json.loads((package/'inputs/fit_manifest.json').read_text())
        controls = m.root_controls(plan)
        self.assertEqual(controls['max_evaluations'],12)
        self.assertEqual(controls['market_tolerance'],2e-4)
        self.assertEqual(controls['fiscal_tolerance'],2e-5)
        self.assertEqual(controls['final_reproduction_tolerance'],1e-10)
        self.assertEqual(controls['price_bound_ratios'],plan['path']['price_bound_ratios'])
        self.assertEqual(controls['pension_bound_ratios'],plan['path']['pension_bound_ratios'])
        plan['gates']['market_tolerance'] *= 2
        with self.assertRaisesRegex(ValueError,'gates differ'):
            m.root_controls(plan)

    def test_root_callback_cap_includes_replay(self):
        for completed in range(12):
            m.require_root_callback_slot(completed)
        with self.assertRaisesRegex(ValueError,'Twelve-call'):
            m.require_root_callback_slot(12)
        m.require_native_slot(719,time.monotonic()+10)
        with self.assertRaisesRegex(ValueError,'720-native-call'):
            m.require_native_slot(720,time.monotonic()+10)
        with self.assertRaisesRegex(ValueError,'deadline'):
            m.require_native_slot(0,time.monotonic()-1)

    def test_accelerator_pin_uses_actual_identity_schema(self):
        package = Path(__file__).resolve().parents[4] / 'output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/execution_smoke_v5'
        plan = json.loads((package/'inputs/fit_manifest.json').read_text())
        expected = package/'frozen/source'/m.ACCEL_RELATIVE
        self.assertEqual(m.require_accelerator_source(SimpleNamespace(__file__=str(expected)),package,plan)['sha256'],
                         plan['identity']['source_pins'][m.ACCEL_RELATIVE])
        plan['identity']['source_pins'][m.ACCEL_RELATIVE]='0'*64
        with self.assertRaisesRegex(ValueError,'implementation differs'):
            m.require_accelerator_source(SimpleNamespace(__file__=str(expected)),package,plan)

    def test_v2_config_and_self_source_pins_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)/'initial_path_diagnostic_v2';directory.mkdir()
            v2config=directory/'config.json';v2config.write_text(json.dumps({'driver':{'path':'wrong','sha256':'0'*64}}))
            helper=directory/'diagnose_initial_path.py';helper.write_text('')
            config=dict(v2_config={'path':str(v2config),'sha256':m.sha(v2config)},
                        v2_helper={'path':str(helper),'sha256':m.sha(helper)},
                        driver={'path':str(Path(m.__file__).resolve()),'sha256':m.sha(m.__file__)})
            with self.assertRaisesRegex(ValueError,'Approved v2 helper pin differs'):
                m.load_contract(config)
            config['driver']['sha256']='0'*64
            with self.assertRaisesRegex(ValueError,'Pinned file differs'):
                m.load_contract(config)
            output=Path(tmp)/'failed_run'
            with self.assertRaisesRegex(ValueError,'Pinned file differs'):
                m.execute(config,output)
            self.assertEqual(json.loads((output/'failure.json').read_text())['phase'],'config')

    def test_five_fresh_native_seed_records_required(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp).resolve()
            (folder/'baseline_checks.json').write_text(json.dumps({'valid':True,'terminal':{'all_checks_pass':True}}))
            pins=[]
            for i in range(1,6):
                record=folder/f'map_{i:03d}'/'native_record.json'
                record.parent.mkdir();record.write_text(json.dumps({'accounting_valid':True,
                    'gates':{key:True for key in m.NATIVE_GATE_KEYS}}))
                pins.append({'path':str(record),'sha256':m.sha(record)})
            d=SimpleNamespace(pinned=m.pinned)
            runtime=SimpleNamespace(rt=SimpleNamespace(total_native_calls=10,
                                                       identity=lambda:{'source':'frozen'}))
            seed=dict(accounting_valid=True,policy_calls=10,mapping_count=5,
                      horizon=12,perturbed_date=5,source_evidence=pins,
                      identity=dict(source='frozen',stage_start_year=2007,
                                    inherited_state_sha256='state'))
            self.assertEqual(len(m.validate_seed(seed,folder,runtime,d,'state',0)),5)
            seed['source_evidence']=pins[:-1]
            with self.assertRaisesRegex(ValueError,'Five pinned'):
                m.validate_seed(seed,folder,runtime,d,'state',0)
            seed['source_evidence']=pins
            first=folder/'map_001/native_record.json'
            first.write_text(json.dumps({'accounting_valid':True,'gates':{'mass':True}}))
            seed['source_evidence'][0]['sha256']=m.sha(first)
            with self.assertRaisesRegex(ValueError,'native gates failed'):
                m.validate_seed(seed,folder,runtime,d,'state',0)

    def test_lossless_numeric_receipts_and_unsupported_type(self):
        matrix=np.arange(48*48,dtype=np.float64).reshape(48,48)/7.
        payload=dict(seed_matrix=matrix,root=dict(history=[dict(step=np.int64(2),
                    accepted=np.bool_(True),residual=np.array([.1,.2,-.3],dtype=np.float64))]),
                    scalar=np.float64(1.25))
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'receipt.json'
            m.write_output(path,payload)
            restored=json.loads(path.read_text())
            self.assertEqual(restored['seed_matrix'],matrix.tolist())
            self.assertEqual(len(restored['seed_matrix']),48)
            self.assertEqual(len(restored['seed_matrix'][0]),48)
            self.assertEqual(restored['root']['history'][0]['residual'],[.1,.2,-.3])
            self.assertEqual(restored['root']['history'][0]['step'],2)
            self.assertIs(restored['root']['history'][0]['accepted'],True)
            self.assertEqual(restored['scalar'],1.25)
            self.assertNotIn('...',path.read_text())
            original=path.read_bytes()
            with self.assertRaisesRegex(TypeError,'Unsupported receipt type'):
                m.write_output(path,{'object':object()})
            self.assertEqual(path.read_bytes(),original)
            with self.assertRaisesRegex(ValueError,'Only finite numeric arrays'):
                m.write_output(path,{'array':np.array([1.,np.nan])})
            self.assertEqual(path.read_bytes(),original)
            with self.assertRaisesRegex(ValueError,'Only finite numeric arrays'):
                m.write_output(path,{'array':np.array(['unsupported'])})
            self.assertEqual(path.read_bytes(),original)


if __name__=='__main__':
    unittest.main()
