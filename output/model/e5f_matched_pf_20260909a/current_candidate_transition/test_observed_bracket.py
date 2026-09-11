import unittest
from run_observed_bracket import collect_observer


class ObserverTests(unittest.TestCase):
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
