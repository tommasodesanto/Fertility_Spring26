"""Bounded source preservation checks; no model solve or plotting."""
import ast
import importlib.util
import json
import unittest
import signal
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("fast_objective", HERE / "fast_objective.py")
fast = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fast)


def function(path, name):
    source = path.read_text()
    node = next(n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef) and n.name == name)
    return "\n".join(source.splitlines()[node.lineno-1:node.end_lineno]) + "\n"


class SourcePreservation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original = (function(fast.NATIVE / "phase_b_pilot.py", "observe_price"),
                        function(fast.NATIVE / "phase_b_pilot.py", "run_phase_b"),
                        function(fast.NATIVE / "runner.py", "native_evaluator"))
        cls.modified = fast.variants(*cls.original)

    def test_variants_compile(self):
        for source in self.modified:
            compile(source, "generated", "exec")

    def test_root_numerical_prefix_exact(self):
        boundary = '    repeat, repeat_live = trial(root["price"]'
        prefix = self.original[1].split(boundary)[0]
        self.assertTrue(self.modified[1].startswith(prefix))
        self.assertIn('final=True', prefix)
        self.assertNotIn('repeat_live', self.modified[1])

    def test_observer_gates_and_all_measurements_retained(self):
        self.assertEqual(self.original[0].split('        # Native plotting expects')[0],
                         self.modified[0].split('        result["target_fit_rows"]')[0])
        for token in ('native_population_step(', 'validate_parameter_estimates(', 'fp.gates(',
                      'runtime.score_targets(', 'len(fits) == 14 and len(params) == 31'):
            self.assertIn(token, self.modified[0])
        self.assertNotIn('standard_diagnostics(', self.modified[0])

    def test_native_initializer_preserved(self):
        expected = self.original[2].split('    def evaluate(')[0].replace(
            '    import phase_b_pilot as ge\n', '    ge = _get_fast_ge(out)\n')
        self.assertEqual(expected, self.modified[2].split('    def evaluate(')[0])
        self.assertIn('ge.validate_parameter_estimates(candidate', self.modified[2])

    def test_ambiguous_boundary_rejected(self):
        with self.assertRaises(RuntimeError):
            fast.variants(self.original[0] * 2, *self.original[1:])

    def test_deadline_wrapper_resolves_exploratory_observer(self):
        source = function(fast.NATIVE / 'phase_b_pilot.py', '_observe_with_deadline')
        original = dict(time=time, signal=signal, _alarm=lambda *_: None,
                        _require=lambda ok, message: self.assertTrue(ok, message),
                        observe_price=lambda *a, **kw: 'original_with_plots')
        exec(compile(source, 'native_wrapper', 'exec'), original)
        exploratory = dict(original)
        exploratory['observe_price'] = lambda *a, **kw: ('exploratory_without_plots', kw['final'])
        exec(compile(source, 'exploratory_wrapper', 'exec'), exploratory)
        result = exploratory['_observe_with_deadline'](
            dict(deadline_epoch=time.time()+10), {}, 'test', final=True)
        self.assertEqual(result, ('exploratory_without_plots', True))

    def test_saved_comparison_all_rows(self):
        full = fast.NATIVE / 'local_run/floor/smoke/000_baseline/phase_b_ge/selected_root'
        if not full.exists():
            self.skipTest('Local full baseline unavailable')
        receipt = fast.compare_saved_baseline({'report': str(full)}, full)
        self.assertEqual(receipt['checks']['target_fit.csv']['rows'], 14)
        self.assertEqual(receipt['checks']['parameters.csv']['rows'], 31)


if __name__ == '__main__':
    unittest.main()
