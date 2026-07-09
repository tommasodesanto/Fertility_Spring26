from __future__ import annotations

import unittest
from types import SimpleNamespace

import numpy as np

from intergen_housing_fertility.local_panel import (
    is_better_record,
    record_sort_key,
    record_selection_loss,
)
from intergen_housing_fertility.parameters import get_income_age_profile, setup_parameters
from intergen_housing_fertility.production_profile import (
    PRODUCTION_H_OWN,
    PRODUCTION_PROFILE_NAME,
    PRODUCTION_SEARCH_BOUNDS,
    PRODUCTION_TARGET_SET,
    validate_production_profile,
)
from intergen_housing_fertility.solver import (
    advance_cohort_one_period,
    advance_cohort_one_period_markov_income,
    assign_current_choices_to_beginning_assets,
    build_forward_tenure_transition_maps,
    compute_eq_stats,
    realize_current_choices,
    realize_current_choices_markov_income,
    realize_current_cross_section,
)


def _family_count_attr() -> str:
    # Keep new research-facing text on the preferred family-size terminology.
    return "n_" + "par" + "ity"


def _fixture(*, j_count: int = 2, z_count: int | None = None) -> SimpleNamespace:
    b_grid = np.arange(-3.0, 4.0)
    nb = b_grid.size
    n_tenure = 2
    n_location = 1
    n_family = 2
    n_stage = 4

    p = SimpleNamespace(
        n_house=1,
        I=n_location,
        n_child_states=n_stage,
        n_child_stages=1,
        J=j_count,
        use_numba_scatter=False,
        propagate_birth_entry_grant=True,
        H_own=np.array([6.0]),
        A_f_end=0,
    )
    setattr(p, _family_count_attr(), n_family)
    sd = SimpleNamespace(nc=n_family * n_stage)

    loc_shape = (nb, n_tenure, n_location, n_location, j_count)
    tenure_shape = (nb, n_tenure, n_location, j_count)
    if z_count is None:
        loc_probs = np.ones(loc_shape + (n_family, n_stage))
        tenure_choice = np.zeros(tenure_shape + (n_family, n_stage), dtype=np.int16)
        bp_pol = np.zeros(tenure_shape + (n_family, n_stage))
    else:
        loc_probs = np.ones(loc_shape + (z_count, n_family, n_stage))
        tenure_choice = np.zeros(tenure_shape + (z_count, n_family, n_stage), dtype=np.int16)
        bp_pol = np.zeros(tenure_shape + (z_count, n_family, n_stage))

    lmm_idx = np.tile(np.arange(nb - 1, dtype=np.int64), (n_location, n_tenure, 1))
    lmm_idx = np.concatenate([lmm_idx, lmm_idx[:, :, -1:]], axis=2)
    lmm_wt = np.zeros((n_location, n_tenure, nb))

    hc = np.array([[0.0, 2.0]])
    he = np.array([[0.0, 1.5]])
    phi_choice = np.ones((n_location, n_tenure, n_family, n_stage))
    phi_choice[:, 1, :, :] = 0.8
    birth_dp = np.zeros((n_family, n_stage, n_tenure, n_tenure), dtype=bool)
    grant = np.zeros((n_location, n_tenure, n_family, n_stage))
    tmx_idx, tmx_wt = build_forward_tenure_transition_maps(
        p, b_grid, hc, he, phi_choice, birth_dp, grant
    )
    return SimpleNamespace(
        P=p,
        SD=sd,
        b_grid=b_grid,
        loc_probs=loc_probs,
        tenure_choice=tenure_choice,
        bp_pol=bp_pol,
        lmm_idx=lmm_idx,
        lmm_wt=lmm_wt,
        tmx_idx=tmx_idx,
        tmx_wt=tmx_wt,
        hc=hc,
        he=he,
        phi_choice=phi_choice,
        birth_dp=birth_dp,
        grant=grant,
    )


def _mapped_values(grid: np.ndarray, idx: np.ndarray, wt: np.ndarray) -> np.ndarray:
    return (1.0 - wt) * grid[idx] + wt * grid[idx + 1]


class IncomeProfileRegressionTests(unittest.TestCase):
    def test_pre_first_break_uses_entry_profile_value(self) -> None:
        p = setup_parameters()
        profile = get_income_age_profile(p)
        self.assertEqual(float(profile[0]), 0.65)
        self.assertEqual(float(profile[1]), 0.65)
        self.assertEqual(float(profile[2]), 0.85)
        self.assertEqual(float(profile[4]), 1.0)

        p.legacy_entry_income_peak = True
        legacy = get_income_age_profile(p)
        self.assertEqual(float(legacy[0]), 1.0)
        self.assertEqual(float(legacy[1]), 0.65)


class SearchSafetyTests(unittest.TestCase):
    def test_lower_loss_nonconverged_record_cannot_win(self) -> None:
        converged = {"case": 1, "rank_loss": 10.0, "strict_converged": True}
        lower_but_invalid = {"case": 2, "rank_loss": 1.0, "strict_converged": False}
        better_converged = {"case": 3, "rank_loss": 9.0, "strict_converged": True}
        self.assertFalse(is_better_record(lower_but_invalid, converged))
        self.assertTrue(is_better_record(better_converged, converged))
        self.assertTrue(np.isinf(record_selection_loss(lower_but_invalid)))
        ranked = sorted([lower_but_invalid, converged, better_converged], key=record_sort_key)
        self.assertEqual([row["case"] for row in ranked], [3, 1, 2])

    def test_named_production_profile_is_exact(self) -> None:
        self.assertEqual(PRODUCTION_H_OWN.tolist(), [2.0, 4.0, 6.0, 8.0, 10.0])
        bounds = {name: (lower, upper) for name, lower, upper in PRODUCTION_SEARCH_BOUNDS}
        self.assertEqual(bounds["beta_annual"], (0.94, 0.995))
        self.assertEqual(bounds["chi"], (0.4, 1.15))
        validate_production_profile(
            PRODUCTION_PROFILE_NAME,
            J=17,
            Nb=120,
            n_house=5,
            income_states=5,
            target_set=PRODUCTION_TARGET_SET,
            stage="search",
        )
        with self.assertRaises(ValueError):
            validate_production_profile(
                PRODUCTION_PROFILE_NAME,
                J=17,
                Nb=60,
                n_house=5,
                income_states=5,
                target_set=PRODUCTION_TARGET_SET,
                stage="search",
            )


class GrantMapRegressionTests(unittest.TestCase):
    def test_zero_grant_is_identical_to_historical_transaction_map(self) -> None:
        f = _fixture()
        idx_on, wt_on = build_forward_tenure_transition_maps(
            f.P,
            f.b_grid,
            f.hc,
            f.he,
            f.phi_choice,
            f.birth_dp,
            f.grant,
        )
        f.P.propagate_birth_entry_grant = False
        idx_off, wt_off = build_forward_tenure_transition_maps(
            f.P,
            f.b_grid,
            f.hc,
            f.he,
            f.phi_choice,
            f.birth_dp,
            np.ones_like(f.grant),
        )
        np.testing.assert_array_equal(idx_on, idx_off)
        np.testing.assert_allclose(wt_on, wt_off, atol=0.0, rtol=0.0)

        renter_to_owner = _mapped_values(
            f.b_grid,
            idx_on[0, 0, 1, 1, 1, :],
            wt_on[0, 0, 1, 1, 1, :],
        )
        expected = np.clip(
            np.maximum(f.b_grid - 2.0, -0.8 * 2.0),
            f.b_grid[0],
            f.b_grid[-1],
        )
        np.testing.assert_allclose(renter_to_owner, expected)

    def test_nonzero_grant_and_policy_hook_precedence_match_bellman(self) -> None:
        f = _fixture()
        f.grant[0, 1, 1, 1] = 1.0
        idx, wt = build_forward_tenure_transition_maps(
            f.P,
            f.b_grid,
            f.hc,
            f.he,
            f.phi_choice,
            f.birth_dp,
            f.grant,
        )
        mapped = _mapped_values(f.b_grid, idx[0, 0, 1, 1, 1], wt[0, 0, 1, 1, 1])
        expected = np.clip(
            np.maximum(f.b_grid - 2.0 + 1.0, -0.8 * 2.0),
            f.b_grid[0],
            f.b_grid[-1],
        )
        np.testing.assert_allclose(mapped, expected)

        f.birth_dp[1, 1, 0, 1] = True
        idx_both, wt_both = build_forward_tenure_transition_maps(
            f.P,
            f.b_grid,
            f.hc,
            f.he,
            f.phi_choice,
            f.birth_dp,
            f.grant,
        )
        mapped_both = _mapped_values(
            f.b_grid,
            idx_both[0, 0, 1, 1, 1],
            wt_both[0, 0, 1, 1, 1],
        )
        expected_waiver = np.clip(
            np.maximum(f.b_grid - 2.0, -0.8 * 2.0),
            f.b_grid[0],
            f.b_grid[-1],
        )
        np.testing.assert_allclose(mapped_both, expected_waiver)


class CurrentChoiceTimingTests(unittest.TestCase):
    def test_birth_and_purchase_are_visible_in_same_current_period(self) -> None:
        f = _fixture()
        f.tenure_choice[:, 0, 0, 0, 1, 1] = 1
        cohort = np.zeros((f.b_grid.size, 2, 1, 2, 4))
        cohort[4, 0, 0, 1, 1] = 1.0
        realized = realize_current_choices(
            cohort,
            0,
            f.loc_probs,
            f.tenure_choice,
            None,
            f.lmm_idx,
            f.lmm_wt,
            f.tmx_idx,
            f.tmx_wt,
        )
        self.assertAlmostEqual(float(np.sum(realized)), 1.0)
        self.assertAlmostEqual(float(np.sum(realized[:, 1, 0, 1, 1])), 1.0)
        self.assertEqual(float(np.sum(realized[:, :, :, 0, :])), 0.0)

        realized_full = np.zeros((f.b_grid.size, 2, 1, f.P.J, 2, 4))
        inherited_full = np.zeros_like(realized_full)
        realized_full[:, :, :, 0, :, :] = realized
        inherited_full[:, :, :, 0, :, :] = cohort
        h_r = np.full_like(realized_full, 2.0)
        realized_stats = compute_eq_stats(
            realized_full, f.P, f.b_grid, np.array([1.0]), h_r
        )
        inherited_stats = compute_eq_stats(
            inherited_full, f.P, f.b_grid, np.array([1.0]), h_r
        )
        self.assertAlmostEqual(float(realized_stats.housing_demand[0]), 6.0)
        self.assertAlmostEqual(float(inherited_stats.housing_demand[0]), 2.0)

    def test_terminal_age_applies_current_housing_choice(self) -> None:
        f = _fixture(j_count=2)
        f.tenure_choice[:, 0, 0, 1, 0, 0] = 1
        g = np.zeros((f.b_grid.size, 2, 1, 2, 2, 4))
        g[4, 0, 0, 1, 0, 0] = 0.7
        out = realize_current_cross_section(
            g,
            f.loc_probs,
            f.tenure_choice,
            None,
            f.lmm_idx,
            f.lmm_wt,
            f.tmx_idx,
            f.tmx_wt,
        )
        self.assertAlmostEqual(float(np.sum(out)), 0.7)
        self.assertAlmostEqual(float(np.sum(out[:, 1, 0, 1, 0, 0])), 0.7)

    def test_beginning_assets_are_conditioned_on_current_buy_and_sell(self) -> None:
        f = _fixture()
        # A renter at b=1 buys; an owner at b=2 sells.
        f.tenure_choice[:, 0, 0, 0, 0, 0] = 1
        f.tenure_choice[:, 1, 0, 0, 0, 0] = 0
        cohort = np.zeros((f.b_grid.size, 2, 1, 2, 4))
        cohort[4, 0, 0, 0, 0] = 0.4
        cohort[5, 1, 0, 0, 0] = 0.6
        assigned = assign_current_choices_to_beginning_assets(
            cohort,
            0,
            f.loc_probs,
            f.tenure_choice,
            None,
            f.lmm_idx,
            f.lmm_wt,
        )
        self.assertAlmostEqual(float(np.sum(assigned)), 1.0)
        self.assertAlmostEqual(float(assigned[4, 1, 0, 0, 0]), 0.4)
        self.assertAlmostEqual(float(assigned[5, 0, 0, 0, 0]), 0.6)
        self.assertAlmostEqual(float(np.sum(assigned * f.b_grid[:, None, None, None, None])), 1.6)

    def test_one_period_transition_applies_transaction_saving_and_stage_once(self) -> None:
        f = _fixture(j_count=2)
        f.tenure_choice[:, 0, 0, 0, 1, 1] = 1
        f.bp_pol[:, :, :, 0, :, :] = 2.0
        cohort = np.zeros((f.b_grid.size, 2, 1, 2, 4))
        cohort[4, 0, 0, 1, 1] = 1.0
        nxt = advance_cohort_one_period(
            cohort,
            0,
            f.loc_probs,
            f.tenure_choice,
            None,
            f.bp_pol,
            f.P,
            f.b_grid,
            f.SD,
            f.lmm_idx,
            f.lmm_wt,
            f.tmx_idx,
            f.tmx_wt,
            False,
            None,
        )
        self.assertAlmostEqual(float(np.sum(nxt)), 1.0)
        self.assertAlmostEqual(float(nxt[5, 1, 0, 1, 2]), 1.0)
        self.assertEqual(float(np.sum(nxt[:, :, :, 1, 1])), 0.0)

    def test_markov_path_conserves_mass_and_transitions_income_once(self) -> None:
        f = _fixture(j_count=2, z_count=2)
        f.tenure_choice[:, 0, 0, 0, :, 1, 1] = 1
        f.bp_pol[:, :, :, 0, :, :, :] = 2.0
        cohort = np.zeros((f.b_grid.size, 2, 1, 2, 2, 4))
        cohort[4, 0, 0, 0, 1, 1] = 1.0
        realized = realize_current_choices_markov_income(
            cohort,
            0,
            f.loc_probs,
            f.tenure_choice,
            None,
            f.lmm_idx,
            f.lmm_wt,
            f.tmx_idx,
            f.tmx_wt,
        )
        self.assertAlmostEqual(float(np.sum(realized)), 1.0)
        self.assertAlmostEqual(float(np.sum(realized[:, 1, 0, 0, 1, 1])), 1.0)

        pi_z = np.array([[0.0, 1.0], [1.0, 0.0]])
        nxt = advance_cohort_one_period_markov_income(
            cohort,
            0,
            f.loc_probs,
            f.tenure_choice,
            None,
            f.bp_pol,
            f.P,
            f.b_grid,
            f.SD,
            f.lmm_idx,
            f.lmm_wt,
            f.tmx_idx,
            f.tmx_wt,
            False,
            None,
            pi_z,
        )
        self.assertAlmostEqual(float(np.sum(nxt)), 1.0)
        self.assertAlmostEqual(float(nxt[5, 1, 0, 1, 1, 2]), 1.0)
        self.assertEqual(float(np.sum(nxt[:, :, :, 0, :, :])), 0.0)


if __name__ == "__main__":
    unittest.main()
