"""No-solve contract tests for simple fertility nests and retained coordinates."""
import contextlib
import copy
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code/model'), str(ROOT / 'code/model/tools')]
import run_e5f_transition_calibration as cal
import run_e5f_bounded_calibration_refinement as adapter


class FertilityNestContractTests(unittest.TestCase):
    def args(self):
        return SimpleNamespace(fertility_nest_choice=True, fixed_tenure_choice_kappa=.005,
            housing_supply_elasticity=.63, model_profile=cal.REPAIRED_MODEL_PROFILE,
            estimate_first_child_room_jump=True, target_profile=cal.E5_BASELINE_TARGET_PROFILE,
            first_child_room_jump_upper=.5)

    def test_original_coordinates_and_menu(self):
        theta = {'kappa_fert': 2.1681730392479377, 'kappa_fert_continuation': 1.7364706586958831}
        domain, overrides, profile = cal.activate_model_profile(cal.REPAIRED_MODEL_PROFILE, theta)
        domain, overrides, profile = cal.configure_first_child_room_jump(
            model_profile_name=cal.REPAIRED_MODEL_PROFILE, theta=theta, active_domain=domain,
            profile_overrides=overrides, model_profile=profile, fixed_jump=None, estimate_jump=True)
        before = copy.deepcopy((domain, theta))
        cal.configure_fertility_nest_calibration(self.args(), overrides, profile)
        self.assertEqual((domain, theta), before)
        self.assertEqual(len(domain), 11)
        self.assertEqual([r for r in domain if r[0] == 'hbar_first_child_jump'][0][2], .5)
        names = {r[0] for r in domain}
        self.assertTrue({'kappa_fert', 'kappa_fert_continuation'} <= names)
        self.assertTrue({'tenure_choice_kappa', 'joint_nest_lambda'}.isdisjoint(names))
        self.assertTrue(overrides['joint_nested_choice'] and overrides['fertility_nest_choice'])
        self.assertFalse(overrides['two_shock_choice'])
        self.assertEqual(profile['fertility_nest']['housing_products'], 6)
        self.assertFalse(profile['fertility_nest']['tenure_committed_across_conception'])
        self.assertNotIn('joint_nested', profile)

    def test_rejects_changed_external_restrictions_targets_and_bounds(self):
        for key, value in [('fixed_tenure_choice_kappa', .01), ('housing_supply_elasticity', 1.),
                           ('estimate_first_child_room_jump', False), ('target_profile', 'young_ownership'),
                           ('first_child_room_jump_upper', 2.)]:
            args = self.args()
            setattr(args, key, value)
            with self.subTest(key=key), self.assertRaises(ValueError):
                cal.configure_fertility_nest_calibration(args, {}, {})

    def test_scale_order_rejected_without_modifying_coordinates(self):
        for bad in (.004, float('nan'), float('inf')):
            theta = {'tenure_choice_kappa': .005, 'kappa_fert': bad, 'kappa_fert_continuation': 1.7}
            with self.subTest(scale=bad), self.assertRaises(ValueError):
                cal.validate_fertility_nest_scales(theta)
            self.assertIs(theta['kappa_fert'], bad)
        cal.validate_fertility_nest_scales(dict(tenure_choice_kappa=.005, kappa_fert=.005,
                                              kappa_fert_continuation=.005))

    def test_choice_flags_are_mutually_exclusive(self):
        with patch.object(sys, 'argv', ['driver', '--fertility-nest-choice']):
            args = cal.parse_args()
            self.assertTrue(args.fertility_nest_choice)
            self.assertFalse(args.joint_nested_choice or args.two_shock_choice)
        for incompatible in ('--two-shock-choice', '--joint-nested-choice', '--exhaustive-saving-control'):
            with patch.object(sys, 'argv', ['driver', '--fertility-nest-choice', incompatible]), \
                    contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                cal.parse_args()

    def test_preflight_keeps_full_target_domain_and_checks_pins(self):
        theta = dict(beta=.9806088988890026, kappa_fert=2.1681730392479377,
            kappa_fert_continuation=1.7364706586958831, chi=1.0434717373613915,
            H0=14.562959141565095, theta0=.5284284711333161, theta1=.10724930821495539,
            hbar_child_rooms=.2822101230841891, first_birth_fixed_cost=4.5591384147193255,
            hbar_first_child_jump=.3649311350610432, psi_child=0., tenure_choice_kappa=.005)
        model = SimpleNamespace(__file__=str(ROOT / 'code/model/intergen_eqscale_seq_optimized/solver.py'))
        target = cal.e5_target_system_for_profile(cal.E5_BASELINE_TARGET_PROFILE)
        self.assertEqual(target.fingerprint, adapter.TARGET)
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            source = folder / 'source.json'
            source.write_text('{}')
            center = folder / 'center.json'
            center.write_text(json.dumps(dict(best_candidate=dict(theta=theta,
                old_psi_child=.1, new_psi_child=-.040697054175712566))))
            argv = ['driver', '--source', str(source), '--outdir', str(folder / 'out'),
                '--fertility-nest-choice', '--no-plots', '--validate-only',
                '--model-profile', cal.REPAIRED_MODEL_PROFILE, '--estimate-first-child-room-jump',
                '--first-child-room-jump-upper', '.5', '--fixed-tenure-choice-kappa', '.005',
                '--housing-supply-elasticity', '.63', '--panel-center-json', str(center),
                '--panel-task-id', '1', '--panel-size', '1', '--expected-source-sha256', cal.source_sha256(source),
                '--expected-target-set', target.name, '--expected-target-fingerprint', target.fingerprint,
                '--expected-code-bundle-sha256', cal.code_fingerprint_contract(model)['bundle_sha256']]
            with patch.object(cal.transition, 'configure_sequential_model', return_value=(None, model)), \
                 patch.object(cal.closure, 'load_winner', return_value=({}, {})), \
                 patch.object(cal.closure, 'theta_from_winner', side_effect=lambda _: dict(theta)), \
                 patch.object(cal.closure, 'make_overrides', return_value={}), \
                 patch.object(cal, 'solve_old_steady_state', side_effect=AssertionError('preflight must not solve')), \
                 contextlib.redirect_stdout(io.StringIO()):
                with patch.object(sys, 'argv', argv):
                    cal.main()
                receipt = json.loads((folder / 'out/validated_inputs.json').read_text())
                self.assertEqual(receipt['target_count'], 12)
                self.assertEqual(len(receipt['search_domain']), 11)
                self.assertEqual(receipt['base_choice_flags'], dict(joint_nested_choice=True,
                    fertility_nest_choice=True, two_shock_choice=False, exhaustive_saving_control=False))
                expected_domain = [r for r in cal.E5F_INCOME_ENTRY_DOMAIN if r[0] != 'psi_child']
                expected_domain += [('hbar_first_child_jump', 0., .5, 'softzero'),
                                    ('psi_child_change_2023', -1.5, .2, 'asinh')]
                self.assertEqual(receipt['search_domain'], [list(r) for r in expected_domain])
                self.assertEqual(receipt['theta']['kappa_fert'], theta['kappa_fert'])
                self.assertEqual(receipt['theta']['kappa_fert_continuation'], theta['kappa_fert_continuation'])
                self.assertEqual(receipt['targets'], [dict(moment=n, target=t, weight=w) for n, t, w in
                    zip(target.moment_names, target.target_values, target.weights)])
                for flag in ('--expected-source-sha256', '--expected-target-fingerprint', '--expected-code-bundle-sha256'):
                    changed = argv.copy()
                    changed[changed.index(flag) + 1] = 'invalid'
                    with self.subTest(pin=flag), patch.object(sys, 'argv', changed), self.assertRaises(RuntimeError):
                        cal.main()

    def test_invalid_center_scale_cannot_be_clipped_into_admissibility(self):
        with tempfile.TemporaryDirectory() as temporary:
            center = Path(temporary) / 'center.json'
            center.write_text(json.dumps(dict(theta=dict(kappa_fert=.004,
                kappa_fert_continuation=1.7), old_psi_child=.1, new_psi_child=0.)))
            args = SimpleNamespace(panel_task_id=1, panel_size=1, panel_local_radius=.02,
                panel_center_json=center, fertility_nest_choice=True, fixed_tenure_choice_kappa=.005)
            with self.assertRaisesRegex(ValueError, 'must be >= housing'):
                cal.panel_candidate({}, args)

    def test_sequential_exhaustive_control_is_separate(self):
        args = self.args()
        args.exhaustive_saving_control = True
        overrides, profile = {}, {}
        cal.configure_exhaustive_saving_control(args, overrides, profile)
        self.assertTrue(overrides['exhaustive_saving_control'])
        self.assertFalse(any(overrides[name] for name in
            ('fertility_nest_choice', 'two_shock_choice', 'joint_nested_choice')))
        self.assertIn('sequential_exhaustive', profile)
        self.assertNotIn('fertility_nest', profile)
        args.first_child_room_jump_upper = 2.
        with self.assertRaises(ValueError):
            cal.configure_exhaustive_saving_control(args, {}, {})
        with patch.object(sys, 'argv', ['driver', '--exhaustive-saving-control']):
            parsed = cal.parse_args()
            self.assertTrue(parsed.exhaustive_saving_control)
            self.assertFalse(parsed.fertility_nest_choice or parsed.joint_nested_choice or parsed.two_shock_choice)

    def test_adapter_requires_pinned_source_and_explicit_classification(self):
        plan = dict(schema='e5f_bounded_refinement_v1', source_sha256=adapter.SOURCE,
            target_fingerprint=adapter.TARGET, code_bundle_sha256='f' * 64,
            choice_model='fertility_nest', suppress_plots=True, cases=[dict(id=1)])
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'plan.json'
            def load(payload):
                adapter.write_json(path, payload)
                return adapter.load_plan(path, adapter.digest(path))
            with patch.object(adapter, 'FERTILITY_NEST_BUNDLE', None), self.assertRaises(RuntimeError):
                load(plan)
            with patch.object(adapter, 'FERTILITY_NEST_BUNDLE', 'f' * 64):
                self.assertEqual(load(plan), plan)
                control = {**plan, "choice_model": "sequential_exhaustive"}
                self.assertEqual(load(control), control)
                for change in (dict(suppress_plots=False), dict(choice_model='two_shock'),
                               dict(choice_model=None), dict(source_sha256='changed'),
                               dict(target_fingerprint='changed')):
                    with self.subTest(change=change), self.assertRaises(RuntimeError):
                        load({**plan, **change})
                with self.assertRaises(RuntimeError):
                    adapter.load_plan(path, 'wrong-plan-hash')


if __name__ == '__main__':
    unittest.main()
