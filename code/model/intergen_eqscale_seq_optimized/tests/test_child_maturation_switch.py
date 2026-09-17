"""Tests for the default-off parent-age child-maturation switch (2026-09-17).

Exemption implementation: m-d draws (option B in the task).  A child born
in the current period (birth indicator ``d = 1``) is safe; the maturation
draw applies to the remaining ``m - d`` children at home.  No one-bit
``born_this_period`` flag is carried: it would double the child state and
every downstream policy array.  Instead the standard and exempt binomial
rows are blended by the newborn share of each post-birth cell (exact for
cell totals and entrant flows), and Bellman birth values carry a uniform
envelope shift (first-order; exact re-optimization is second-order).
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from intergen_eqscale_seq_optimized.parameters import (
    apply_overrides,
    child_exit_prob_by_age,
    setup_parameters,
)
from intergen_eqscale_seq_optimized.solver import (
    advance_cohort_one_period_markov_income,
    run_model_cp_dt,
)
from intergen_eqscale_seq_optimized.tests.test_eqscale_seq import _tiny_markov


def _tiny_independent(**extra: object) -> dict[str, object]:
    base = {
        **_tiny_markov(),
        "J": 17,
        "J_R": 12,
        "A_f_start": 1,
        "A_f_end": 7,
        "sequential_births": True,
        "child_state_mode": "independent_count",
        "n_parity": 4,
        "preference_spec": "eqscale",
        "delta_alpha": 0.05,
        "delta_alpha_jump": 0.03,
        "gamma_e": 0.3,
        "psi_child": -0.2,
        "kappa_fert": 1.5,
    }
    base.update(extra)
    return base


def test_constant_mode_arrays_bitwise_identical() -> None:
    """Task test 1: constant mode reproduces current arrays bit for bit."""
    implicit = apply_overrides(setup_parameters(), _tiny_independent(Nb=20))
    explicit = apply_overrides(
        setup_parameters(),
        _tiny_independent(
            Nb=20,
            child_maturation_mode="constant",
            mu_young=0.05,
            a_rise=34.0,
            a_full=62.0,
        ),
    )
    assert not hasattr(implicit, "Pi_child_by_age")
    assert not hasattr(explicit, "Pi_child_by_age")
    assert np.array_equal(implicit.Pi_child, explicit.Pi_child)
    sol_implicit, _, _ = run_model_cp_dt(_tiny_independent(Nb=20), verbose=False)
    sol_explicit, _, _ = run_model_cp_dt(
        _tiny_independent(Nb=20, child_maturation_mode="constant"), verbose=False
    )
    assert np.array_equal(sol_implicit.V, sol_explicit.V)
    assert np.array_equal(sol_implicit.g, sol_explicit.g)


def test_default_off_bitwise() -> None:
    """Default (no maturation keys) matches explicit constant on the solve."""
    sol_default, _, _ = run_model_cp_dt(_tiny_independent(Nb=20), verbose=False)
    sol_off, _, _ = run_model_cp_dt(
        _tiny_independent(Nb=20, child_maturation_mode="constant"), verbose=False
    )
    assert np.array_equal(sol_default.V, sol_off.V)
    assert np.array_equal(sol_default.g, sol_off.g)
    assert float(np.sum(sol_default.g)) == pytest.approx(
        float(np.sum(sol_off.g)), rel=0.0, abs=0.0
    )


def test_parent_age_matrices_hazard_and_exemption() -> None:
    """Task test 2: row-stochastic, per-child hazard mu(a), newborn identity."""
    P = apply_overrides(
        setup_parameters(),
        {
            "J": 17,
            "n_parity": 4,
            "child_state_mode": "independent_count",
            "use_stochastic_aging": True,
            "child_maturation_mode": "parent_age",
        },
    )
    npar = int(P.n_parity)
    by_age = np.asarray(P.Pi_child_by_age)
    ex_age = np.asarray(P.Pi_child_exempt_by_age)
    assert by_age.shape == (17, npar, npar, npar)
    assert ex_age.shape == (17, npar, npar, npar)
    for j in range(int(P.J)):
        mu = child_exit_prob_by_age(P, j)
        std = by_age[j]
        exm = ex_age[j]
        np.testing.assert_allclose(std.sum(axis=1), 1.0, rtol=0.0, atol=1e-14)
        np.testing.assert_allclose(exm.sum(axis=1), 1.0, rtol=0.0, atol=1e-14)
        for nn in range(npar):
            for m in range(npar):
                if m > nn:
                    continue
                stay = sum(nxt * std[m, nxt, nn] for nxt in range(npar))
                assert stay == pytest.approx(m * (1.0 - mu), rel=1e-12, abs=1e-14)
                stay_ex = sum(nxt * exm[m, nxt, nn] for nxt in range(npar))
                expect_ex = 0.0 if m == 0 else 1.0 + (m - 1) * (1.0 - mu)
                assert stay_ex == pytest.approx(expect_ex, rel=1e-12, abs=1e-14)
                if m >= 1:
                    # Newborn component: no mass can fall to zero children.
                    assert float(exm[m, 0, nn]) == 0.0
        # Single-newborn cell is the identity: the birth-period child stays.
        for nn in range(1, npar):
            row = exm[1, :, nn]
            assert float(row[1]) == pytest.approx(1.0, abs=1e-15)
            assert float(np.sum(row) - row[1]) == pytest.approx(0.0, abs=1e-15)
    # Hazard profile: flat young, rising, absorbing at the top.
    assert child_exit_prob_by_age(P, 0) == pytest.approx(0.05)
    assert child_exit_prob_by_age(P, 4) == pytest.approx(0.05)
    mid = child_exit_prob_by_age(P, 8)
    assert 0.05 < mid < 1.0
    assert child_exit_prob_by_age(P, 11) == pytest.approx(1.0)
    assert child_exit_prob_by_age(P, 16) == pytest.approx(1.0)


def _advance_fixture(mode: str):
    P = apply_overrides(
        setup_parameters(),
        {
            "J": 6,
            "I": 1,
            "Nb": 20,
            "n_parity": 4,
            "child_state_mode": "independent_count",
            "use_stochastic_aging": True,
            "use_numba_scatter": False,
            "H_own": np.array([2.0, 4.0]),
            "child_maturation_mode": mode,
        },
    )
    Nb, nt, I, J, Nz, npar, ncs = 20, 3, 1, 6, 1, 4, 4
    b_grid = np.linspace(0.0, 3.0, Nb)
    SD = SimpleNamespace(nc=npar * ncs)
    rng = np.random.default_rng(7)
    gj = rng.random((Nb, nt, I, Nz, npar, ncs)) + 0.01
    gj = gj / float(np.sum(gj))
    loc_probs = np.full((Nb, nt, I, I, J, Nz, npar, ncs), 1.0 / I)
    tenure_choice = np.zeros((Nb, nt, I, J, Nz, npar, ncs), dtype=np.int64)
    bp_pol = np.zeros((Nb, nt, I, J, Nz, npar, ncs))
    for b in range(Nb):
        bp_pol[b] = b_grid[b]
    lmm_idx = np.zeros((I, nt, Nb), dtype=np.int64)
    lmm_wt = np.zeros((I, nt, Nb))
    for b in range(Nb):
        lmm_idx[:, :, b] = b
    tmx_idx = np.zeros((I, nt, nt, npar, ncs, Nb), dtype=np.int64)
    tmx_wt = np.zeros((I, nt, nt, npar, ncs, Nb))
    for b in range(Nb):
        tmx_idx[..., b] = b
    Pi_z = np.eye(Nz)
    Pia = P.Pi_child
    return P, SD, b_grid, gj, loc_probs, tenure_choice, bp_pol, lmm_idx, lmm_wt, tmx_idx, tmx_wt, Pia, Pi_z


def test_population_operator_conserves_mass() -> None:
    """Task test 3: one-period advance conserves mass in both modes."""
    for mode in ("constant", "parent_age"):
        (
            P, SD, b_grid, gj, loc_probs, tenure_choice, bp_pol,
            lmm_idx, lmm_wt, tmx_idx, tmx_wt, Pia, Pi_z,
        ) = _advance_fixture(mode)
        out = advance_cohort_one_period_markov_income(
            gj, 2, loc_probs, tenure_choice, None, bp_pol, P, b_grid, SD,
            lmm_idx, lmm_wt, tmx_idx, tmx_wt, True, Pia, Pi_z,
        )
        np.testing.assert_allclose(float(np.sum(out)), float(np.sum(gj)), rtol=0.0, atol=1e-12)
        assert float(np.min(out)) >= -1e-14
    # Blended (newborn-share) rows also conserve mass exactly.
    (
        P, SD, b_grid, gj, loc_probs, tenure_choice, bp_pol,
        lmm_idx, lmm_wt, tmx_idx, tmx_wt, Pia, Pi_z,
    ) = _advance_fixture("parent_age")
    rng = np.random.default_rng(11)
    frac = rng.random((int(P.n_parity), int(P.n_child_states)))
    out = advance_cohort_one_period_markov_income(
        gj, 2, loc_probs, tenure_choice, None, bp_pol, P, b_grid, SD,
        lmm_idx, lmm_wt, tmx_idx, tmx_wt, True, Pia, Pi_z, newborn_frac=frac,
    )
    np.testing.assert_allclose(float(np.sum(out)), float(np.sum(gj)), rtol=0.0, atol=1e-12)


def test_smoke_solves_tiny_config_both_modes() -> None:
    """Task test 4: tiny (Nb=20, J=17, n_parity=4) solves in both modes."""
    sol_c, _, _ = run_model_cp_dt(_tiny_independent(Nb=20), verbose=False)
    sol_p, P_p, _ = run_model_cp_dt(
        _tiny_independent(Nb=20, child_maturation_mode="parent_age"), verbose=False
    )
    assert np.all(np.isfinite(sol_p.V))
    assert float(np.sum(sol_p.g)) == pytest.approx(1.0, abs=1e-9)
    assert hasattr(P_p, "Pi_child_by_age")
    assert sol_p.g.shape == sol_c.g.shape
    # Parent-age hazard empties old-age cells: no child at home at the top age.
    old_home = float(np.sum(sol_p.g[:, :, :, -1, :, :, 1:]))
    assert old_home == pytest.approx(0.0, abs=1e-9)


def test_maturation_switch_rejects_bad_combos() -> None:
    with pytest.raises(ValueError):
        apply_overrides(setup_parameters(), {"child_maturation_mode": "bogus"})
    with pytest.raises(ValueError):
        apply_overrides(
            setup_parameters(),
            {"child_state_mode": "shared_clock", "child_maturation_mode": "parent_age"},
        )
