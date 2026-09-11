import json
import unittest
from run_observed_bracket import collect_observer, jsonable


class ObserverTests(unittest.TestCase):
    def test_actual_numpy_diagnostics_are_json_serializable(self):
        import numpy as np
        value = {'mappings': [[{'rate': np.float64(1.5),
            'ages': np.array([18., 22.]), 'flows': np.array([[.1, .2]])}]]}
        saved = json.loads(json.dumps(jsonable(value)))
        self.assertEqual(saved['mappings'][0][0]['flows'], [[.1, .2]])
        self.assertEqual(saved['mappings'][0][0]['rate'], 1.5)

    def test_preserves_existing_observer_and_dated_mapping_sequence(self):
        original_calls, measured, saved = [], [], []
        e, p, g, s = object(), object(), object(), object()
        def original(*args):
            original_calls.append(args)
        def measure(a, b):
            self.assertIs(a, e); self.assertIs(b, p)
            self.assertEqual(len(original_calls), len(measured)+1)
            measured.append(1)
            return {'rate': 2.0}
        observer = collect_observer(original, measure, lambda x: saved.append(str(x)))
        for index in [0, 1, 0, 1]:
            observer(index, e, p, g, s)
        self.assertEqual(len(saved), 4)
        self.assertEqual(original_calls[1], (1, e, p, g, s))
        self.assertIn('2011', saved[-1])

    def test_rejects_missing_date(self):
        observer = collect_observer(None, lambda *a: {}, lambda x: None)
        with self.assertRaises(RuntimeError):
            observer(1, None, None, None, None)

    def test_measurement_failure_is_not_silently_ignored(self):
        def fail(*args):
            raise ValueError('Accounting mismatch')
        observer = collect_observer(None, fail, lambda x: None)
        with self.assertRaisesRegex(ValueError, 'Accounting mismatch'):
            observer(0, None, None, None, None)


if __name__ == '__main__':
    unittest.main()
