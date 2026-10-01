"""Bounded orchestration tests using the actual unchanged two-block root."""
import gzip
import copy
import importlib
import json
from pathlib import Path
import pickle
import sys
import tempfile
import time
from types import SimpleNamespace as NS
import unittest

import numpy as np
import joint_six_date as joint


class JointTests(unittest.TestCase):
    def test_actual_seed_source_compatibility_and_unknown_change_rejected(self):
        base=joint.frozen.ROOT/'output/model/transition_readiness_v1/joint_preparation'
        config=json.loads((base/'config.json').read_text())
        seed=json.loads((base/'measured_seed/receipt.json').read_text())
        proof=joint.seed_compatibility(config,seed)
        self.assertEqual(proof['status'],'PASS')
        changed=copy.deepcopy(seed);changed['source_pins']['run_e5f_preference_transition.py']='0'*64
        with self.assertRaisesRegex(ValueError,'Unexpected'):
            joint.seed_compatibility(config,changed)

    def exercise(self, market, fiscal):
        sys.path.insert(0,str(joint.frozen.HERE/'pinned_tools'))
        acceleration=importlib.import_module('e5f_four_shock_acceleration')
        accounting=importlib.import_module('run_e5f_preference_budget_diagnostic')
        calls=[];P=NS(pension=1.);packet=dict(parameters=P,solution=NS(p_eq=[1.]));terminal=dict(parameters=P)
        def checkpoint(path,value):
            with gzip.open(path,'wb') as stream:pickle.dump(value,stream)
        def mapping(*args,**kwargs):
            q,b=args[4:6];calls.append((q.copy(),b.copy()))
            record=dict(market_residual=(market-np.log(q)).tolist(),fiscal_residual=(fiscal-np.log(b)).tolist(),
                gates=dict(budget=True,transaction=True,funded=True,mass=True),rows=[],fertility=[],diagnostic_packets=[])
            return NS(terminal_state={},rows=[]),record
        inner=NS(mapping=mapping,write=accounting.write,dump_checkpoint=checkpoint,plain=lambda x:x,
            terminal_checks=lambda *args:dict(all_checks_pass=False,raw_queue_pass=False))
        plan=dict(path=dict(price_bound_ratios=[.05,20.],pension_bound_ratios=[.05,20.],market_tolerance=2e-4,
            fiscal_tolerance=1e-6,market_slope=1.,fiscal_slope=1.,max_log_step=.15,damping=.7,
            final_reproduction_tolerance=1e-10),fit=dict(max_condition_number=1e8,worsening_factor=1.5))
        with tempfile.TemporaryDirectory() as td:
            output=Path(td)
            receipt,latest=joint.solve_paths(inner=inner,acceleration=acceleration,packet=packet,evaluator=None,
                terminal=terminal,endpoint=dict(price=1.),plan=plan,jacobian=-np.eye(12),output=output,deadline=time.monotonic()+1950,native=False)
            assert (output/'latest_completed.json').is_file() and (output/'best_so_far.json').is_file()
            assert len(list(output.glob('map_*_checkpoint.pkl.gz')))==len(calls)
            assert (output/'root.json').is_file() and (output/'joint_receipt.json').is_file()
            return receipt,calls

    def test_actual_root_one_step_and_fresh_repeat(self):
        receipt,calls=self.exercise(.0003,.000002)
        self.assertEqual(len(calls),3)
        self.assertTrue(receipt['root_certified'])
        np.testing.assert_array_equal(calls[1][0],calls[2][0])
        np.testing.assert_array_equal(calls[1][1],calls[2][1])
        self.assertFalse(receipt['terminal_pass'])
        self.assertTrue(receipt['terminal_failure_interpretable'])
        self.assertFalse(receipt['production_ready'])

    def test_unqualified_candidate_stops_before_replay(self):
        receipt,calls=self.exercise(.1,.001)
        self.assertEqual(len(calls),2)
        self.assertFalse(receipt['root_certified'])
        self.assertFalse(receipt['terminal_failure_interpretable'])

    def test_deadline_replay_and_render_reserve(self):
        budget=joint.MapBudget(time.monotonic()+1200)
        with self.assertRaisesRegex(TimeoutError,'reserving'):
            budget.reserve(False)
        self.assertEqual(budget.started,0)

    def test_three_map_cap(self):
        budget=joint.MapBudget(time.monotonic()+1950)
        for _ in range(3):budget.reserve(True)
        with self.assertRaisesRegex(TimeoutError,'Three'):
            budget.reserve(True)


if __name__=='__main__':unittest.main()
