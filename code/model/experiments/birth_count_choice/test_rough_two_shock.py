"""Focused zero-solve fake-loop checks for the isolated rough two-shock route."""
from __future__ import annotations
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import test_two_shock as old_tests
import rough_two_shock as rough
import rough_two_shock_runtime as runtime


def plan(*,smoke=False):
    p=old_tests.plan()
    p.update(schema=rough.SCHEMA,kind='two_unanticipated_permanent_rough',mode='rough_diagnostic',
        baseline_psi=rough.BASELINE_EXPECTED,absolute_psi_bounds=rough.BOUNDS,
        stages=copy.deepcopy(rough.STAGES),horizons=[6] if smoke else [16],smoke=smoke,
        smoke_seed_endpoint_padding=smoke,gates=dict(rough.ROUGH_GATES),
        stage_starts=[dict(initial=x,bounds=rough.BOUNDS) for x in
                      ((rough.BASELINE_EXPECTED,rough.BASELINE_EXPECTED) if smoke else (.14736308634876963,.12))],
        budget=dict(rough.SMOKE_BUDGET if smoke else rough.ROUGH_BUDGET))
    p['seed']=dict(horizon=4 if smoke else 12,perturbed_date=1 if smoke else 5,log_step=1e-5)
    p['fit']['max_evaluations']=12;p['fit']['fertility_tolerance']=.005
    p['path']['max_evaluations']=3 if smoke else 12
    return p


class Tests(unittest.TestCase):
    def execute(self,smoke=False):
        with tempfile.TemporaryDirectory() as d:
            p=plan(smoke=smoke);fake=old_tests.FakeRuntime(p)
            controller=rough.Controller(p,fake,Path(d)/'result')
            with patch.object(rough,'preflight',return_value=dict(status='PASS')):
                result=controller.run()
            return result,fake,controller

    def test_h16_sequential_fit_and_exact_handoff(self):
        result,fake,controller=self.execute()
        self.assertEqual(result['status'],'rough_matched')
        self.assertEqual([x[1] for x in fake.calls if x[0]=='seed'],[2007,2015])
        self.assertEqual([x[2] for x in fake.calls if x[0]=='export'],[4,2])
        self.assertEqual(len(result['fit_table']),4)
        self.assertEqual([x['weight'] for x in result['fit_table']],[0,1,0,1])
        self.assertFalse(result['original_horizon_gates_tested'])
        self.assertEqual(len(fake.initialization[0]),12)
        self.assertEqual(fake.replayed['start_year'],2007)
        np.testing.assert_array_equal(fake.replayed['boundary_value'],np.array([2015.]))
        self.assertGreater(result['actual_policy_calls'],0)

    def test_failed_accounting_is_never_scored(self):
        p=plan();fake=old_tests.FakeRuntime(p)
        old=fake.evaluate_stage
        def bad(**kw):
            reply=old(**kw);reply['accounting_valid']=False;return reply
        fake.evaluate_stage=bad
        with tempfile.TemporaryDirectory() as d:
            controller=rough.Controller(p,fake,Path(d))
            controller.current_state=fake.initial_state()
            with self.assertRaisesRegex(ValueError,'accounting'):
                controller.evaluate(0,.14736308634876963)

    def test_absolute_bounds_and_h16_slice(self):
        self.assertEqual(rough.BOUNDS,[0.001789207206604163,0.35784144132083257])
        self.assertEqual(rough.ROUGH_BUDGET['maximum_policy_calls'],19898)
        self.assertEqual(len(np.arange(16)[2:14]),12)
        self.assertEqual(len(np.arange(12)[2:14]),10)
        self.assertEqual(len(runtime.no_arbitrage_prices(old_tests.NS(R_gross=1.02,delta=.01,tau_H=.01,
            user_cost_rate=.04),.8,1.,16)),16)

    def test_tiny_smoke_exercises_both_clocks_without_fitting(self):
        result,fake,_=self.execute(smoke=True)
        self.assertEqual(result['status'],'execution_passed')
        self.assertFalse(result['empirical_fitted'])
        self.assertFalse(result['empirical_h16_seed_tested'])
        self.assertEqual([x[1] for x in fake.calls if x[0]=='seed'],[2007,2015])

    def test_owned_superseded_packet_pruning_keeps_numeric_record(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)/'map_001';packet=folder/'date_000'/'diagnostic_packet.pkl.gz'
            packet.parent.mkdir(parents=True);packet.write_bytes(b'large-packet')
            record=dict(diagnostic_packets=[dict(path=str(packet),sha256=old_tests.m.sha(packet))],
                        rows=[dict(calendar_year=2007)],market_residual=[0.],fiscal_residual=[0.])
            old_tests.m.write(folder/'native_record.json',record)
            self.assertEqual(runtime.prune_map_packets(folder,'superseded'),1)
            self.assertFalse(packet.exists())
            retained=json.loads((folder/'native_record.json').read_text())
            self.assertEqual(retained['rows'],record['rows'])
            self.assertEqual(retained['diagnostic_packets'],[])

    def test_flat_bridge_target_fingerprints_and_new_wealth_row(self):
        with tempfile.TemporaryDirectory() as d:
            case=Path(d)
            rows=['moment,target','wealth_earnings,4.45838713455674']
            rows.extend(f'other_{i},1' for i in range(13))
            (case/'target_fit_new_contract.csv').write_text('\n'.join(rows)+'\n')
            contract=dict(target_fingerprint='c7a3d185668122e508a6c322bc5ef0715ebb0ecb23948c8d9b184ee25d1cde70',
                weight_fingerprint='f762ebb5684ab30487b3b8b64fc10977fda396b520035d91c0c5c803255f88e4')
            (case/'input_contract.json').write_text(json.dumps(contract))
            self.assertEqual(rough.check_new_wealth_target(case)['sha256'],
                             old_tests.m.sha(case/'target_fit_new_contract.csv'))
            contract['weight_fingerprint']='wrong'
            (case/'input_contract.json').write_text(json.dumps(contract))
            with self.assertRaisesRegex(ValueError,'fingerprints'):
                rough.check_new_wealth_target(case)

if __name__=='__main__':unittest.main()
