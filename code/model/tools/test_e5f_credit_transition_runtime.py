import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
import numpy as np
import e5f_credit_transition_runtime as runtime
import e5f_exact_policy_cache as cache

# Load the pure native queue file directly; do not import model package/Numba.
path = Path(__file__).resolve().parents[1] / 'intergen_eqscale_seq_optimized/adult_entry.py'
spec = importlib.util.spec_from_file_location('credit_transition_test_adult_entry', path)
entry = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = entry
spec.loader.exec_module(entry)

class RuntimeTests(unittest.TestCase):
    def test_off_is_noop(self):
        obj = SimpleNamespace()
        self.assertTrue(runtime.install_credit(obj, None, '/not-created')['baseline_noop'])
        self.assertTrue(runtime.install_split_queue(obj, '/not-created', enabled=False)['baseline_noop'])
        self.assertEqual(vars(obj), {})

    def test_exact_impulse(self):
        slots = [0.] * 7
        results = []
        for period in range(7):
            due, slots = runtime.split_queue_step(slots, 2.1 if period == 0 else 0, 1/2.1, entry.SplitBirthEntryQueue)
            results.append(due)
        self.assertEqual(results, [0.,0.,0.,.5,.5,0.,0.])
        self.assertEqual(sum(results), 1.)

    def test_constant_actual_births(self):
        slots = runtime.flatten_queue(entry.SplitBirthEntryQueue.constant_prehistory(4.2))
        for _ in range(8):
            due, future = runtime.split_queue_step(slots, 4.2, 1/2.1, entry.SplitBirthEntryQueue)
            self.assertEqual(due,2.)
            self.assertEqual(slots,future)
            slots = future

    def test_bad_conversion_or_slots(self):
        for q,c in (([0.]*7,1.),([0.]*4,1/2.1)):
            with self.assertRaises(ValueError):
                runtime.split_queue_step(q, 2.1,c,entry.SplitBirthEntryQueue)

    def test_initial_state_preserves_distribution_and_raw_units(self):
        pf = SimpleNamespace(PFInitialState=SimpleNamespace)
        g = np.ones((1,1,1,2,1,1,1))
        state = runtime.initialize_state(pf,entry.SplitBirthEntryQueue,g,4.2,2.1)
        np.testing.assert_array_equal(state.g_pre,g)
        self.assertIsNot(state.g_pre,g)
        self.assertEqual(state.scheduled_entries,[1.]*7)
        self.assertEqual(state.scheduled_raw_entries,[.5]*7)

    def test_continuation_changes_only_guard(self):
        guard='if continuation_V is not None or not bool(getattr(P, "exhaustive_saving_control", False)):'
        text='    '+guard+'\n        raise ValueError("guard")\n'
        result=runtime.continuation_source(text)
        self.assertEqual(result,text.replace('continuation_V is not None or ',''))
        with self.assertRaises(ValueError):runtime.continuation_source(result)

    def test_queue_source_only_two_calls_and_guards(self):
        text='    started = time.perf_counter()\n'+('    x=transition.advance_birth_vintage_queue(q,b,c)\n'*2)
        patched=runtime.queue_source(text)
        self.assertEqual(patched.count('_credit_split_queue_step('),2)
        self.assertIn('historical_conditioning is not None',patched)
        with self.assertRaises(ValueError):runtime.queue_source(text+'transition.advance_birth_vintage_queue(q,b,c)')

    def test_existing_cache_exact_argument_and_result(self):
        calls=[]
        def solve_date_policy(*,continuation_V):
            calls.append(1)
            return continuation_V.copy()+1
        pf=SimpleNamespace(solve_date_policy=solve_date_policy)
        with runtime.exact_cache_context(pf,cache,enabled=True,max_bytes=100000) as stats:
            a=pf.solve_date_policy(continuation_V=np.array([1.,2.]))
            b=pf.solve_date_policy(continuation_V=np.array([1.,2.]))
            np.testing.assert_array_equal(a,b)
            self.assertEqual(len(calls),1)
            self.assertEqual(stats.hits,1)
        self.assertIs(pf.solve_date_policy,solve_date_policy)

if __name__=='__main__':unittest.main()
