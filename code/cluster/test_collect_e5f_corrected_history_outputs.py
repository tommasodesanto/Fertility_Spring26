"""Collection-contract tests with synthetic files; no economic solve."""
import importlib.util
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

spec = importlib.util.spec_from_file_location(
    'collector', Path(__file__).with_name('collect_e5f_corrected_history_outputs.py'))
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


class CollectionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.b = Path(self.tmp.name).resolve()
        b = self.b
        self.case = b / 'case'
        self.case.mkdir()
        self.args = SimpleNamespace(scientific_manifest=b / 'scientific.json',
            helper_root=b / 'helper', reader=b / 'reader.py', profile=b / 'profile.py',
            once=True, command_timeout=30)
        self.source = b / 'model'
        (self.source / 'code/model/tools').mkdir(parents=True)
        (b / 'helper').mkdir()
        self.kernel = self.source / 'code/model/tools/observer.py'
        for p in [b / 'helper/runtime.py', b / 'reader.py', b / 'profile.py', self.kernel]:
            p.write_text('# pinned\n')
        prior = b / 'plan.json'
        c.atomic_json(prior, {'source_root': str(self.source),
                             'file_sha256': {str(self.kernel): c.sha(self.kernel)}})
        c.atomic_json(self.args.scientific_manifest, {'prior_plan': str(prior),
            'file_sha256': {str(p): c.sha(p) for p in [prior, b / 'helper/runtime.py', self.kernel]}})
        c.atomic_json(self.case / 'contract_receipt.json', {'source_root': str(self.source),
            'manifest_sha256': c.sha(self.args.scientific_manifest), 'case': 'A0', 'count': 6})
        self.fits = [dict(year=y, folder=str(self.case / f'window_{y}/trial_00'),
                         target=1.8, model=1.801, gap=.001, psi=.1) for y in c.WINDOW]

    def tearDown(self):
        self.tmp.cleanup()

    def completed_fixture(self):
        c.atomic_json(self.case / 'realized_fit.json', self.fits)
        c.atomic_json(self.case / 'finite_history_complete.json', {'realized': self.fits})
        folder = Path(self.fits[-1]['folder']) / 'alternative'
        folder.mkdir(parents=True, exist_ok=True)
        for name in ['accepted_forecast.pkl.gz', 'first_period_diagnostics.pkl.gz',
                     'native_2023_snapshot.pkl.gz']:
            (folder / name).write_bytes(b'fixture-only; not a native checkpoint')
        c.atomic_json(folder / 'root_receipt.json', dict(start_year=2019, count=6,
            case='A0', converged=True, finite_horizon_market_fiscal_converged=True,
            final={'mapping_valid': True}, final_reproduction_max_abs=0, psi=.1))
        c.atomic_json(folder / 'fertility.json', [dict(calendar_year=2019,
                                                     period_tfr_topcode_adjusted=1.801)])

    def inspect(self):
        return c.inspect_case(self.case, self.args, self.b / 'out')

    def test_pending_is_not_a_terminal_error(self):
        self.assertEqual(self.inspect()['status'], 'PENDING')

    def test_alternative_and_pinned_imports(self):
        self.completed_fixture()
        state = self.inspect()
        self.assertEqual(state['status'], 'READY', state)
        manifest = c.load(Path(state['manifest']))
        self.assertIn('/alternative/', manifest['accepted_forecast'])
        self.assertEqual(manifest['runtime_paths'][0], str(self.b / 'helper'))
        self.assertIn(str(self.kernel), manifest['kernel_files'])

    def test_earlier_window_gap_is_checked(self):
        self.fits[0]['gap'] = .006
        self.completed_fixture()
        self.assertEqual(self.inspect()['status'], 'ERROR')

    def test_changed_scientific_kernel_is_rejected(self):
        self.completed_fixture()
        (self.b / 'helper/runtime.py').write_text('# changed\n')
        self.assertEqual(self.inspect()['status'], 'ERROR')

    def test_completed_case_is_not_rerun(self):
        self.completed_fixture()
        c.atomic_json(self.b / 'out/case/status.json', dict(status='COMPLETE', proof='saved'))
        self.assertEqual(self.inspect()['proof'], 'saved')


if __name__ == '__main__':
    unittest.main()
