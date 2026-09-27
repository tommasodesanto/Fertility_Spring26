"""Synthetic age-25 observation checks; no model import or solve is required.

Run on Torch with the tools directory on PYTHONPATH. The lead also runs the
existing initial-fertility tests against real accounting primitives there.
"""
import copy
import json
import sys
import unittest
from types import ModuleType, SimpleNamespace
from unittest.mock import patch

import numpy as np

import e5f_initial_fertility_observer as observer


def fixture(*, age25_pre=(1., 0., 0., 0.), age25_flows=(1., 0., 0.), masses=None):
    """Hand-constructed ever-born transitions; current children are all zero."""
    masses = np.ones(7) if masses is None else np.asarray(masses, dtype=float)
    stock = np.tile([.4, .3, .2, .1], (7, 1))
    stock[1] = age25_pre
    before = masses[:, None] * stock
    flows = np.zeros((7, 3))
    flows[0, 0] = .1 * masses[0]  # Keep timing defined outside the age-25 cell.
    flows[1] = masses[1] * np.asarray(age25_flows)
    after = before.copy()
    after[:, :3] -= flows
    after[:, 1:] += flows
    pre = np.zeros((1, 1, 1, 7, 1, 4, 4))
    post = pre.copy()
    pre[0, 0, 0, :, 0, :, 0] = before
    post[0, 0, 0, :, 0, :, 0] = after
    P = SimpleNamespace(J=7, age_start=18., da=4., period_years=4., n_parity=4,
        sequential_births=True, fertility_units="literal_topcode",
        child_state_mode="independent_count", A_f_start=1., A_f_end=7.,
        tfr_top_bin_weight=3.602359422009)
    ev = SimpleNamespace(g_pre=pre, g_post_fertility=post, births=float(flows.sum()))
    first = flows[:, 0]
    risk = before[:, 0]
    hazard = np.divide(first, risk, out=np.zeros(7), where=risk > 0.)
    ages = np.arange(18., 46., 4.)
    helper = ModuleType("run_e5f_transition_calibration")
    helper.first_birth_accounting_by_age = lambda evaluation, parameters: {
        "flow": first, "at_risk": risk, "hazard": hazard,
    }
    helper.period_fertility_diagnostics = lambda evaluation, parameters: {
        "birth_flow_first": flows[:, 0],
        "birth_flow_second": flows[:, 1],
        "birth_flow_third_bin_entry": flows[:, 2],
        "age_cell_start": ages,
        "age_cell_midpoint": ages + 2.,
        "age_mass": masses,
        "period_first_birth_mean_age": float(np.dot(ages + 2., first) / first.sum()),
        "period_first_birth_share_age30plus": float(first[ages >= 30.].sum() / first.sum()),
    }
    return ev, P, helper


def observe(ev, P, helper, projection="uniform_birth_time"):
    # Only substitute the two supplied accounting functions; never import the
    # calendar runtime, household solver, or any checkpoint.
    with patch.dict(sys.modules, {"run_e5f_transition_calibration": helper}):
        return observer.observe_initial_fertility(ev, P, age_projection=projection)


class EarlyFertilityObserverTests(unittest.TestCase):
    def test_exact_age25_uniform_interpolation_and_overlap_mass(self):
        packet = observe(*fixture())
        self.assertEqual(packet["moments"][observer.EARLY_FERTILITY_MOMENT], .875)
        self.assertEqual(packet["ever_born_shares_age25"],
                         {"0": .125, "1": .875, "2": 0., "3plus": 0.})
        self.assertEqual(packet["accounting"]["age25_population_mass"], .25)
        self.assertEqual(packet["accounting"]["age25_overlap_weights"],
                         [0., .25, 0., 0., 0., 0., 0.])
        self.assertEqual(packet["accounting"]["age25_post_fertility_interpolation_share"],
                         [None, .875, None, None, None, None, None])

    def test_top_state_is_literal_three_not_completed_fertility_weight(self):
        for top_weight in (3., 3.602359422009, 20.):
            with self.subTest(top_weight=top_weight):
                ev, P, helper = fixture(age25_pre=(0., 0., 1., 0.),
                                        age25_flows=(0., 0., 1.))
                P.tfr_top_bin_weight = top_weight
                packet = observe(ev, P, helper)
                self.assertEqual(packet["moments"][observer.EARLY_FERTILITY_MOMENT], 2.875)
                self.assertEqual(packet["ever_born_shares_age25"]["3plus"], .875)

    def test_ever_born_stock_does_not_count_children_currently_at_home(self):
        ev, P, helper = fixture(age25_pre=(0., 0., 0., 1.), age25_flows=(0., 0., 0.))
        # Both distributions have all their mass at the zero-at-home index.
        packet = observe(ev, P, helper)
        self.assertEqual(packet["moments"][observer.EARLY_FERTILITY_MOMENT], 3.)
        self.assertFalse(packet["metadata"]["early_fertility_age25"]["uses_children_at_home"])

    def test_age25_is_independent_of_legacy_constant_post_cell_selector(self):
        args = fixture()
        uniform = observe(*args)
        constant = observe(*args, projection="constant_post_cell")
        key = observer.EARLY_FERTILITY_MOMENT
        self.assertEqual(uniform["moments"][key], .875)
        self.assertEqual(constant["moments"][key], .875)
        self.assertEqual(constant["metadata"]["age_projection"], "constant_post_cell")
        self.assertEqual(constant["metadata"]["early_fertility_age25"]["age_projection"],
                         "uniform_birth_time")

    def test_age_mass_is_preserved_and_other_age_cells_do_not_enter(self):
        a = observe(*fixture())
        b = observe(*fixture(masses=[2., 8., 5., 10., 20., 30., 40.]))
        key = observer.EARLY_FERTILITY_MOMENT
        self.assertEqual(a["moments"][key], b["moments"][key])
        self.assertEqual(b["accounting"]["age25_population_mass"], 2.)
        self.assertEqual(sum(b["accounting"]["age25_ever_born_mass"]), 2.)

    def test_old_moment_definitions_metadata_and_inputs_are_preserved(self):
        ev, P, helper = fixture()
        pre, post, parameters = ev.g_pre.copy(), ev.g_post_fertility.copy(), copy.deepcopy(vars(P))
        packet = observe(ev, P, helper)
        self.assertAlmostEqual(packet["moments"]["childless_rate_40_44"], .4)
        self.assertAlmostEqual(packet["moments"]["exactly_one_among_mothers_40_44"], .5)
        self.assertAlmostEqual(packet["moments"]["period_mean_age_first_birth"], 260. / 11.)
        self.assertEqual(packet["moments"]["period_share_first_births_age30plus"], 0.)
        self.assertEqual(packet["metadata"]["observer_id"], "initial_cps_nchs_fertility_diagnostic_v1")
        self.assertFalse(packet["metadata"]["production_smm_eligible"])
        self.assertFalse(packet["metadata"]["weights_or_standard_errors_adopted"])
        meta = packet["metadata"]["early_fertility_age25"]
        self.assertEqual(meta["children_ever_born_weights"], [0., 1., 2., 3.])
        self.assertFalse(meta["uses_tfr_top_bin_weight"])
        self.assertIn("early_fertility_target_20260926", meta["empirical_source"])
        np.testing.assert_array_equal(ev.g_pre, pre)
        np.testing.assert_array_equal(ev.g_post_fertility, post)
        self.assertEqual(vars(P), parameters)
        json.dumps(packet, allow_nan=False)

    def test_zero_age25_population_marks_only_new_moment_unavailable(self):
        packet = observe(*fixture(masses=[1., 0., 1., 1., 1., 1., 1.]))
        self.assertIsNone(packet["moments"][observer.EARLY_FERTILITY_MOMENT])
        self.assertIsNone(packet["ever_born_shares_age25"])
        self.assertEqual(packet["accounting"]["age25_population_mass"], 0.)
        self.assertEqual(packet["accounting"]["age25_ever_born_mass"], [0., 0., 0., 0.])
        meta = packet["metadata"]["early_fertility_age25"]
        self.assertFalse(meta["available"])
        self.assertIn("zero projected population mass", meta["unavailable_reason"])
        self.assertIn("registry must reject", meta["required_scored_moment_rule"])
        self.assertAlmostEqual(packet["moments"]["childless_rate_40_44"], .4)
        self.assertAlmostEqual(packet["moments"]["exactly_one_among_mothers_40_44"], .5)
        self.assertEqual(packet["moments"]["period_mean_age_first_birth"], 20.)
        self.assertEqual(packet["moments"]["period_share_first_births_age30plus"], 0.)
        json.dumps(packet, allow_nan=False)

    def test_invalid_geometry_is_rejected_before_using_flow_helpers(self):
        for name, value in (("age_start", 20.), ("da", 2.), ("period_years", 1.),
                            ("J", 6), ("n_parity", 5), ("age_start", float("nan"))):
            with self.subTest(name=name, value=value):
                ev, P, helper = fixture()
                setattr(P, name, value)
                with self.assertRaises(ValueError):
                    observe(ev, P, helper)

    def test_inconsistent_pre_post_age_mass_fails(self):
        ev, P, helper = fixture()
        ev.g_post_fertility[0, 0, 0, 1, 0, 1, 0] += .01
        with self.assertRaisesRegex(RuntimeError, "changes age mass"):
            observe(ev, P, helper)


if __name__ == "__main__":
    unittest.main()
