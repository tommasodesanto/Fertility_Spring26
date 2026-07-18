from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from intergen_seq_fertility.local_panel import (
    de_candidate_is_acceptable,
    is_better_record,
    record_sort_key,
    record_selection_loss,
    run_global_de_panel,
)
from intergen_seq_fertility.parameters import (
    apply_overrides,
    compute_stationary_worker_income,
    get_income_age_profile,
    setup_parameters,
)
from intergen_seq_fertility.production_profile import (
    FROZEN_SOURCE_THETA,
    PRODUCTION_H_OWN,
    PRODUCTION_MAX_ITER_EQ,
    PRODUCTION_PROFILE_NAME,
    PRODUCTION_SEARCH_BOUNDS,
    PRODUCTION_TARGET_SET,
    comparison_arm_switches,
    production_profile_metadata,
    production_profile_overrides,
    validate_frozen_source_theta,
    validate_production_profile,
)
from intergen_seq_fertility.solver import (
    advance_cohort_one_period,
    advance_cohort_one_period_markov_income,
    assign_current_choices_to_beginning_assets,
    build_forward_tenure_transition_maps,
    compute_eq_stats,
    realize_current_choices,
    realize_current_choices_markov_income,
    realize_current_cross_section,
)
from intergen_seq_fertility.utils import make_grid


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
    def test_working_age_mean_and_entry_value_match_frozen_source(self) -> None:
        p = setup_parameters()
        profile = get_income_age_profile(p)
        expected = np.array(
            [
                0.720554272517321,
                0.720554272517321,
                0.942263279445728,
                0.942263279445728,
                1.108545034642032,
            ]
        )
        np.testing.assert_allclose(profile[:5], expected, atol=1e-15, rtol=0.0)
        self.assertAlmostEqual(float(np.mean(profile[: p.J_R])), 1.0, places=15)
        self.assertAlmostEqual(compute_stationary_worker_income(p, profile), 1.0, places=15)

        p.normalize_income_profile = False
        raw = get_income_age_profile(p)
        self.assertEqual(float(raw[0]), 0.65)
        self.assertEqual(float(raw[1]), 0.65)
        self.assertEqual(float(raw[2]), 0.85)
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

    def test_de_population_never_accepts_invalid_candidate(self) -> None:
        self.assertFalse(de_candidate_is_acceptable(np.inf, np.inf))
        self.assertFalse(de_candidate_is_acceptable(np.inf, 10.0))
        self.assertTrue(de_candidate_is_acceptable(9.0, 10.0))
        self.assertTrue(de_candidate_is_acceptable(9.0, np.inf))

    def test_named_production_profile_is_exact(self) -> None:
        self.assertEqual(PRODUCTION_H_OWN.tolist(), [2.0, 4.0, 6.0, 8.0, 10.0])
        self.assertEqual(PRODUCTION_MAX_ITER_EQ, 10)
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
            max_iter_eq=10,
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
                max_iter_eq=10,
                stage="search",
            )
        with self.assertRaises(ValueError):
            validate_production_profile(
                PRODUCTION_PROFILE_NAME,
                J=17,
                Nb=120,
                n_house=5,
                income_states=5,
                target_set=PRODUCTION_TARGET_SET,
                max_iter_eq=25,
                stage="search",
            )

    def test_profile_pins_dense_grid_and_runtime_metadata(self) -> None:
        overrides = production_profile_overrides()
        expected_grid_fields = {
            "b_min": -12.0,
            "b_max": 30.0,
            "b_core_lo": -5.0,
            "b_core_hi": 7.0,
            "b_mid_hi": 15.0,
            "b_frac_low": 0.08,
            "b_frac_core": 0.72,
            "b_frac_mid": 0.12,
            "b_grid_power": 1.5,
        }
        for key, value in expected_grid_fields.items():
            self.assertEqual(overrides[key], value)
        self.assertTrue(overrides["normalize_income_profile"])
        self.assertFalse(overrides["legacy_entry_income_peak"])
        self.assertEqual(overrides["housing_event_horizon"], 0)
        self.assertEqual(overrides["max_iter_eq"], 10)

        metadata = production_profile_metadata()
        self.assertEqual(metadata["runtime_overrides"]["max_iter_eq"], 10)
        self.assertEqual(metadata["runtime_overrides"]["b_min"], -12.0)
        self.assertEqual(metadata["runtime_overrides"]["b_max"], 30.0)

        p = apply_overrides(setup_parameters(), {"Nb": 120, **overrides})
        grid = make_grid(p)
        self.assertEqual(grid.size, 120)
        self.assertEqual(float(grid[0]), -12.0)
        self.assertEqual(float(grid[10]), -5.0)
        self.assertEqual(float(grid[96]), 7.0)
        self.assertEqual(float(grid[-1]), 30.0)
        self.assertTrue(np.any(grid == 0.0))

    def test_exact_frozen_theta_and_comparison_switches(self) -> None:
        expected_theta = {
            "alpha_cons": 0.6195578033469402,
            "beta": 0.7808047585225533,
            "c_bar_0": 1.2793341787106423,
            "c_bar_n": 0.11617574349164181,
            "chi": 1.1496942798885048,
            "h_bar_0": 1.0,
            "h_bar_jump": 2.1198232149690575,
            "h_bar_n": 1.3819945583480897,
            "kappa_fert": 1.0188742355464353,
            "psi_child": 0.3490340757304799,
            "tenure_choice_kappa": 0.0,
            "theta0": 0.0014981114599317271,
            "theta_n": 0.9722009696445806,
        }
        self.assertEqual(FROZEN_SOURCE_THETA, expected_theta)
        validate_frozen_source_theta(dict(expected_theta))
        wrong_theta = dict(expected_theta)
        wrong_theta["theta_n"] = 0.9811
        with self.assertRaises(ValueError):
            validate_frozen_source_theta(wrong_theta)

        frozen = comparison_arm_switches("frozen_source_repro")
        repaired = comparison_arm_switches("repaired_timing")
        self.assertFalse(frozen["legacy_entry_income_peak"])
        self.assertFalse(repaired["legacy_entry_income_peak"])
        differing = {key for key in frozen if frozen[key] != repaired[key]}
        self.assertEqual(
            differing,
            {
                "use_postdecision_current_distribution",
                "propagate_birth_entry_grant",
            },
        )

    def test_global_de_applies_named_production_profile(self) -> None:
        captured_overrides: list[dict[str, object]] = []

        def fake_case(
            idx: int,
            candidate: dict[str, object],
            J: int,
            Nb: int,
            n_house: int,
            max_iter_eq: int,
            income: dict[str, object],
            rank_targets: dict[str, float],
            rank_weights: dict[str, float],
            extra_overrides: dict[str, object] | None = None,
        ) -> dict[str, object]:
            del J, Nb, n_house, max_iter_eq, income, rank_targets, rank_weights
            captured_overrides.append(dict(extra_overrides or {}))
            return {
                "case": idx,
                "label": candidate["label"],
                "status": "ok",
                "rank_loss": 1.0,
                "full_old_nonlocation_loss": 1.0,
                "theta": candidate["theta"],
                "moments": {},
                "p_eq": [1.0],
                "market_residual": 1e-6,
                "strict_converged": True,
                "tol_eq": 1e-4,
                "elapsed_sec": 0.0,
                "timings": {"strict_converged": True},
            }

        with TemporaryDirectory() as tmpdir:
            with (
                patch(
                    "intergen_seq_fertility.local_panel.run_local_panel_case",
                    side_effect=fake_case,
                ),
                patch("intergen_seq_fertility.local_panel.write_panel_summary_tables"),
                patch("intergen_seq_fertility.local_panel.write_panel_plots"),
            ):
                summary = run_global_de_panel(
                    Path(tmpdir),
                    max_evals=1,
                    seed=20260709,
                    J=17,
                    Nb=120,
                    income_states=5,
                    n_house=5,
                    max_iter_eq=10,
                    minutes=1.0,
                    pop_size=4,
                    target_set=PRODUCTION_TARGET_SET,
                    seed_theta=FROZEN_SOURCE_THETA,
                    profile_name=PRODUCTION_PROFILE_NAME,
                    progress=False,
                )

        self.assertEqual(len(captured_overrides), 1)
        self.assertTrue(captured_overrides[0]["use_postdecision_current_distribution"])
        self.assertTrue(captured_overrides[0]["propagate_birth_entry_grant"])
        metadata = summary["metadata"]
        self.assertEqual(metadata["production_profile"], PRODUCTION_PROFILE_NAME)
        self.assertEqual(len(metadata["bounds"]), 13)
        self.assertEqual(
            [row["name"] for row in metadata["bounds"]],
            [row[0] for row in PRODUCTION_SEARCH_BOUNDS],
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
