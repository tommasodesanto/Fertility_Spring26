"""Adaptive price-loop tests with synthetic accounting; zero native solves."""
import csv
import importlib.util
import json
import sys
import tempfile
import time
import types
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

# Import only the owned driver, with no household/runtime dependency loaded.
_SPEC = importlib.util.spec_from_file_location('utility_price_test_driver', Path(__file__).with_name('phase_b_pilot.py'))
phase = importlib.util.module_from_spec(_SPEC)
with patch.dict(sys.modules, {'single_price': types.SimpleNamespace(solve_fixed_price=None)}):
    _SPEC.loader.exec_module(phase)


class PriceChecks(unittest.TestCase):
    def run_mock(self, residual, *, start=None, remaining=32, fail_at=None,
                 repeat_array_drift=False, repeat_field_drift=False,
                 repeat_residual_drift=False, remaining_seconds=10000,
                 strict_writer=False):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        out = Path(self.tmp.name)
        writes = {}
        def write(path, value):
            if not strict_writer:
                path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_suffix(path.suffix+'.tmp')
            temporary.write_text(json.dumps(value))
            temporary.replace(path)
            writes[path.name] = value
        context = {name: None for name in ('prepared', 'manifest', 'objective', 'runtime', 'reference')}
        context.update(fp=types.SimpleNamespace(write=write), P=types.SimpleNamespace(psi_child=.1),
                       q_ref=1., out=out)
        if start is not None:
            context['price_start'] = start
        budget = types.SimpleNamespace(remaining_lifecycle=remaining,
                                       deadline_epoch=time.time()+remaining_seconds,
                                       stage_deadline_seconds=1)
        calls = []
        def solve(ctx, credit, q, b, label, stage_dir):
            if strict_writer:
                # Real fp.write requires the parent directory to exist. Its
                # initial atomic receipt precedes stage-directory creation.
                self.assertTrue((out/'phase_b_ge'/'price_search.json').is_file())
                self.assertEqual(writes['price_search.json']['attempts'][-1]['status'],'started')
            calls.append((q, label))
            b.remaining_lifecycle -= 1
            if fail_at == len(calls):
                raise RuntimeError('Native feasibility projection fails')
            fields = {'a': np.asarray([q])}
            if repeat_array_drift and label == 'selected_repeat':
                fields['a'] = np.asarray([q+1e-10])
            if repeat_field_drift and label == 'selected_repeat':
                fields['extra'] = 1.
            return dict(price=np.asarray([q]), sol=types.SimpleNamespace(**fields),
                        sd=types.SimpleNamespace(b=np.asarray([1.])), case_deadline_epoch=time.time()+300)
        def observe(ctx, live, label, *, final=False):
            q = float(live['price'][0])
            r = float(residual(q))
            if repeat_residual_drift and label.startswith('selected_repeat'):
                r += 1e-12
            if final:
                phase._require(abs(r) <= phase.RENEWAL_TOL, 'Birth renewal root fails')
                dest = out/'phase_b_ge'/label
                dest.mkdir(parents=True, exist_ok=True)
                for filename, n in [('target_fit.csv', 14), ('parameters.csv', 31)]:
                    with (dest/filename).open('w', newline='') as f:
                        writer=csv.DictWriter(f,fieldnames=['name','value'])
                        writer.writeheader()
                        writer.writerows(dict(name=str(i),value='1') for i in range(n))
            return dict(price=q,renewal_residual=r,population_scale=1.)
        self.calls, self.writes, self.budget = calls, writes, budget
        with patch.object(phase, 'solve_fixed_price', solve), patch.object(phase, 'observe_price', observe):
            return phase.run_phase_b(context, dict(selected_d_bar=0., selected_live={'bad':'stale'}), budget)

    def assert_passed(self, result, root):
        self.assertEqual(result['status'], 'passed')
        self.assertLessEqual(abs(result['selected_price']-root), 1e-6)
        self.assertEqual(self.calls[-1], (result['selected_price'], 'selected_repeat'))
        self.assertEqual(len(self.calls), 32-self.budget.remaining_lifecycle)
        self.assertEqual(len(result['price_search']['attempts']), len(self.calls))
        self.assertTrue(all(p['status']=='observed' for p in result['price_search']['attempts']))
        self.assertLessEqual(abs(result['selected']['renewal_residual']), phase.RENEWAL_TOL)

    def test_root_below_old_bracket(self):
        result=self.run_mock(lambda q: .4-q)
        self.assert_passed(result,.4)
        self.assertLess(min(q for q,_ in self.calls),.85)

    def test_root_above_old_bracket(self):
        result=self.run_mock(lambda q: 2.5-q)
        self.assert_passed(result,2.5)
        self.assertGreater(max(q for q,_ in self.calls),1.15)

    def test_existing_inside_root(self):
        self.assert_passed(self.run_mock(lambda q: 1.02-q),1.02)

    def test_native_strict_writer_first_receipt_precedes_first_solve(self):
        self.assert_passed(self.run_mock(lambda q: 1.02-q,strict_writer=True),1.02)

    def test_optional_start_fresh_and_global_bounds_origin(self):
        result=self.run_mock(lambda q: 2.6-q,start=2.5)
        self.assert_passed(result,2.6)
        self.assertEqual(self.calls[0],(2.5,'price_start'))
        self.assertEqual([result['price_search'][k] for k in ['lower','upper']], [.125,8.])

    def test_invalid_optional_start_falls_back(self):
        for start in [99.,float('nan'),'invalid']:
            result=self.run_mock(lambda q: 1.-q,start=start)
            self.assert_passed(result,1.)
            self.assertEqual(self.calls[0][0],1.)

    def test_absent_root_hits_finite_caps_and_records_both_sides(self):
        result=self.run_mock(lambda q: 1.)
        self.assertEqual(result['status'],'uncomputed_price_unbracketed')
        self.assertEqual(result['price_search']['termination_reason'],'both_diagnostic_price_caps')
        prices=[q for q,_ in self.calls]
        self.assertIn(.125,prices); self.assertIn(8.,prices)
        self.assertLessEqual(len(prices),21)
        self.assertNotIn('selected_repeat',[label for _,label in self.calls])

    def test_opposite_direction_can_find_nonmonotone_root(self):
        result=self.run_mock(lambda q: (q-.4)*(q+1))
        self.assert_passed(result,.4)
        self.assertGreater(self.calls[1][0],1.)
        self.assertLess(self.calls[2][0],1.)
        self.assertTrue(any('slope' in e['reason'] for e in result['price_search']['events']))

    def test_feasibility_failure_fatal_recorded(self):
        with self.assertRaisesRegex(RuntimeError,'Native feasibility'):
            self.run_mock(lambda q: 2.-q,fail_at=2)
        attempt=self.writes['price_search.json']['attempts'][-1]
        self.assertEqual(attempt['status'],'failed_fatal')
        self.assertEqual(attempt['price'],1.6)

    def test_lifecycle_repeat_reserve(self):
        result=self.run_mock(lambda q: 10.-q,remaining=3)
        self.assertEqual(result['status'],'uncomputed_price_unbracketed')
        self.assertEqual(self.budget.remaining_lifecycle,1)
        self.assertEqual(result['price_search']['termination_reason'],'budget_or_repeat_reserve')

    def test_time_reserve_unchanged(self):
        with self.assertRaisesRegex(RuntimeError,'No time reserve'):
            self.run_mock(lambda q: 1.-q,remaining_seconds=400)
        self.assertEqual(self.calls,[])

    def test_exact_repeat_array_shape_and_tolerance(self):
        for kw,message in [(dict(repeat_array_drift=True),'solution array'),
                           (dict(repeat_field_drift=True),'field set'),
                           (dict(repeat_residual_drift=True),'repeat differs in renewal')]:
            with self.assertRaisesRegex(RuntimeError,message):
                self.run_mock(lambda q: 1.-q,**kw)

    def test_tolerance_at_fresh_seed_still_repeats(self):
        result=self.run_mock(lambda q: .999e-6)
        self.assertEqual(result['status'],'passed')
        self.assertEqual(len(self.calls),2)
        self.assertEqual(phase.RENEWAL_TOL,1e-6)

    def test_outside_tolerance_is_not_accepted(self):
        result=self.run_mock(lambda q: 1.001e-6)
        self.assertEqual(result['status'],'uncomputed_price_unbracketed')

    def test_integrated_zero_lifecycle_smoke(self):
        with tempfile.TemporaryDirectory() as temp:
            def write(path, value):
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(json.dumps(value))
            context={name: None for name in ('prepared','manifest','objective','runtime','reference')}
            context.update(fp=types.SimpleNamespace(write=write),P=types.SimpleNamespace(psi_child=.1),
                           selected_d_bar=0.,q_ref=.8,out=temp)
            budget=types.SimpleNamespace(remaining_lifecycle=32,deadline_epoch=time.time()+10000,
                                         stage_deadline_seconds=1)
            result=phase.smoke_phase_b(context,budget)
            self.assertEqual(result['status'],'passed_mock_zero_lifecycle')
            self.assertEqual(result['lifecycle_solves'],0)
            self.assertEqual(result['selected_price_factor'],1.02)

if __name__=='__main__':
    unittest.main()
