"""Pure shape, fixed-input provenance and supplied-path plumbing checks."""
import copy
from types import SimpleNamespace
import unittest
from unittest import mock

import numpy as np
import run_e5f_matched_pf_baseline as baseline
import run_e5f_matched_pf_preference_pilot as pilot


class PreferenceShapeTest(unittest.TestCase):
    def setUp(self):
        self.old = np.linspace(.29, -.031, 5)

    def test_zero_is_exact_original_array(self):
        np.testing.assert_array_equal(pilot.preference_shape(self.old, 100, 0.),
            np.r_[self.old, np.full(95, self.old[-1])])

    def test_three_shapes_have_exact_endpoints_and_flat_tail(self):
        for a in pilot.SHAPES:
            p = pilot.preference_shape(self.old, 100, a)
            self.assertEqual(p[0], self.old[0])
            np.testing.assert_array_equal(p[4:], np.full(96, self.old[-1]))
            self.assertTrue(np.all(np.diff(p) <= 0))
            x = np.arange(5)/4
            np.testing.assert_allclose(p[:5], self.old[0]+(self.old[-1]-self.old[0])*(x+a*x*(1-x)), atol=1e-16)

    def test_increasing_and_constant_endpoints(self):
        for old in (self.old[::-1], np.full(5, .2)):
            for a in pilot.SHAPES:
                p = pilot.preference_shape(old, 6, a)
                self.assertEqual(p[0], old[0])
                self.assertEqual(p[-1], old[-1])

    def test_uncontracted_shape_rejected(self):
        for a in (float('nan'), float('inf'), .1, 1.):
            with self.assertRaises(ValueError):
                pilot.preference_shape(self.old, 6, a)

    def test_invalid_override_rejected(self):
        valid = pilot.preference_shape(self.old, 6, .5)
        variants = [valid[:-1], np.r_[np.nan, valid[1:]], valid.copy(), valid.copy(), valid.copy()]
        variants[2][0] += 1e-9
        variants[3][-1] += 1e-9
        variants[4][2] = valid[1]+.01
        for values in variants:
            with self.assertRaises(ValueError):
                baseline.validated_preference_path(self.old, 6, values)

    def test_override_is_copied(self):
        p = pilot.preference_shape(self.old, 6, .5)
        out = baseline.validated_preference_path(self.old, 6, p)
        out[1] = 0
        self.assertNotEqual(out[1], p[1])

    def test_probe_passes_exact_override_to_joined_operator(self):
        old = SimpleNamespace(psi_path=self.old)
        expected = pilot.preference_shape(self.old, 6, .5)
        # Stop after the first consumer: proves the hook is invoked by the real
        # probe before any household computation, without mocking a model solve.
        c = dict(path_date_count=6, probe_coordinate=-1, probe_log_step=.01,
            terminal_preference_rule='hold_normalized_2023_intercept',
            initial_price_rule='log_old_to_selected_2023_then_terminal')
        with mock.patch.object(baseline, 'load_normalized', return_value=(old, None)), \
             mock.patch.object(baseline, 'load_terminal', return_value={}), \
             mock.patch.object(baseline, 'validated_preference_path', side_effect=RuntimeError('hook reached')) as hook:
            with self.assertRaisesRegex(RuntimeError, 'hook reached'):
                baseline.run_history_probe(None, c, SimpleNamespace(arm='sequential'),
                    None, {}, None, 0., psi_path_override=expected)
            hook.assert_called_once_with(self.old, 6, expected)


class ParentContractTest(unittest.TestCase):
    def setUp(self):
        keys = ('checkpoint_sha256', 'selected_summary_sha256', 'normalized_checkpoint_sha256',
            'normalized_summary_sha256', 'normalized_contract_sha256', 'terminal_checkpoint_sha256',
            'terminal_summary_sha256', 'terminal_contract_sha256', 'target_fingerprint')
        source = {pilot.DRIVER: 'a'*64, 'code/model/model.py': 'b'*64,
            'code/model/tools/e5f_matched_pf_path_root.py': 'c'*64,
            'code/model/tools/e5f_matched_pf_moments.py': 'd'*64}
        self.parent = dict({k:'1'*64 for k in keys}, arm='sequential', path_date_count=100,
            demographic_sources={'five':'pinned'}, initial_price_rule='inherited',
            terminal_preference_rule='hold_normalized_2023_intercept', probe_log_step=.01,
            source_sha256=source)
        self.c = copy.deepcopy(self.parent)
        self.c['preconditioner'] = 'verified_parent_broyden'
        self.c['source_sha256'][pilot.DRIVER] = 'e'*64
        for name in ('run_e5f_matched_pf_preference_pilot.py', 'test_run_e5f_matched_pf_preference_pilot.py'):
            self.c['source_sha256']['code/model/tools/'+name] = 'f'*64
        self.c['reviewed_preference_driver_change'] = dict(path=pilot.DRIVER,
            from_sha256='a'*64, to_sha256='e'*64, scope=pilot.CHANGE_SCOPE)
        final = dict(prices=[1.]*100, residual=[1e-5]*100, score=1e-5,
            mapping_valid=True, payload={'directory':'/parent/evaluation_003'})
        self.history = dict(final=final, converged=True, final_reproduction_max_abs=0.,
            final_jacobian=(-np.eye(100)).tolist())
        self.summary = dict(best=final, arm='sequential', status='converged',
            finite_horizon_market_converged=True, final_reproduction_max_abs=0.)

    def validate(self):
        return pilot.validate_parent(self.parent, self.history, self.summary, self.c, 'sequential')

    def test_valid_parent_is_seed_only(self):
        p,J,F = self.validate()
        self.assertEqual(p.shape,(100,))
        np.testing.assert_array_equal(J,-np.eye(100))
        np.testing.assert_array_equal(F,np.full(100,1e-5))

    def test_changed_inputs_rejected(self):
        for key in ('normalized_checkpoint_sha256','terminal_checkpoint_sha256','target_fingerprint'):
            old=self.c[key];self.c[key]='changed'
            with self.assertRaisesRegex(ValueError,'input'):
                self.validate()
            self.c[key]=old

    def test_changed_scientific_source_rejected(self):
        self.c['source_sha256']['code/model/model.py']='changed'
        with self.assertRaisesRegex(ValueError,'source changed'):
            self.validate()

    def test_missing_new_source_or_review_rejected(self):
        self.c['reviewed_preference_driver_change']['scope']='arbitrary'
        with self.assertRaisesRegex(ValueError,'reviewed'):
            self.validate()
        self.c['reviewed_preference_driver_change']['scope']=pilot.CHANGE_SCOPE
        del self.c['source_sha256']['code/model/tools/run_e5f_matched_pf_preference_pilot.py']
        with self.assertRaisesRegex(ValueError,'source pins'):
            self.validate()

    def test_unreproduced_or_uncleared_parent_rejected(self):
        self.summary['final_reproduction_max_abs']=1e-9
        with self.assertRaises(ValueError):self.validate()
        self.summary['final_reproduction_max_abs']=0.
        self.history['final']['residual'][0]=.001
        with self.assertRaisesRegex(ValueError,'market gate'):self.validate()

    def test_six_date_smoke_is_explicit_and_diagonal(self):
        self.c['path_date_count']=6
        with self.assertRaisesRegex(ValueError,'six-date'):self.validate()
        self.c['scope']='six_date_plumbing_smoke_only'
        with self.assertRaisesRegex(ValueError,'complete matching horizon'):self.validate()
        self.c['preconditioner']='diagonal'
        p,J,F=self.validate()
        self.assertEqual(len(p),6);self.assertIsNone(J)

    def test_nonfinite_parent_jacobian_rejected(self):
        self.history['final_jacobian'][0][0]=float('nan')
        with self.assertRaisesRegex(ValueError,'preconditioner'):self.validate()


if __name__ == '__main__':
    unittest.main()
