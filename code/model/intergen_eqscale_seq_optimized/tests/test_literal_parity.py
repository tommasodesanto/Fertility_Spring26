"""Tests for the gated literal-parity (L4) extension.

The extension adds four default-off conventions to the sequential/eqscale
architecture: ``n_parity=4`` literal birth states 0/1/2/3+, a generalized
per-parity upward-attempt margin, an ``entrant_conversion_factor`` on the
matured-children entrant flow, a ``fertility_units`` convention for the tfr
moment, and a ``child_bin_high_cutoff`` for the family-size room bins.  With
all defaults the solver must reproduce the pre-extension solution bitwise
(checked against a saved golden by ``golden_l4_check.py``; here we check the
gates are inert when passed explicitly at their defaults).
"""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from intergen_eqscale_seq_optimized.calibration import extract_moments, literal_topcode_tfr
from intergen_eqscale_seq_optimized.solver import (
    current_child_bin_dt,
    get_completed_fertility,
    run_model_cp_dt as run_fork,
)
from intergen_eqscale_seq_optimized.tests.test_eqscale_seq import _tiny_markov


def _seq_eqscale(**extra: object) -> dict[str, object]:
    base = {
        **_tiny_markov(),
        "sequential_births": True,
        "preference_spec": "eqscale",
        "delta_alpha": 0.05,
        "delta_alpha_jump": 0.03,
        "gamma_e": 0.3,
        "psi_child": -0.2,
        "kappa_fert": 1.5,
    }
    base.update(extra)
    return base


L4 = {
    "n_parity": 4,
    "fertility_units": "literal_topcode",
    "tfr_top_bin_weight": 3.4,
    "entrant_conversion_factor": 0.5,
    "child_bin_high_cutoff": 3,
}


def test_gates_passed_at_defaults_are_inert() -> None:
    a, _, _ = run_fork(_seq_eqscale(), verbose=False)
    b, _, _ = run_fork(
        _seq_eqscale(
            n_parity=3,
            fertility_units="parity2x",
            tfr_top_bin_weight=3.0,
            entrant_conversion_factor=1.0,
            child_bin_high_cutoff=2,
        ),
        verbose=False,
    )
    assert np.array_equal(a.V, b.V)
    assert np.array_equal(a.g, b.g)
    assert float(a.childless_rate) == float(b.childless_rate)


def test_completed_fertility_reads_parity() -> None:
    P = SimpleNamespace(n_child_stages=1)
    assert get_completed_fertility(0, 0, P) == 0
    assert get_completed_fertility(1, 1, P) == 1
    assert get_completed_fertility(1, 2, P) == 1  # K+1 matured-one
    assert get_completed_fertility(2, 3, P) == 2  # K+2 matured-2plus
    assert get_completed_fertility(3, 3, P) == 3  # literal third child preserved
    assert get_completed_fertility(3, 1, P) == 3  # raising, parity 3
    # Legacy zero-mass cell mappings preserved verbatim for bitwise nesting
    # (they feed the precomputed bequest table even at zero mass).
    assert get_completed_fertility(0, 2, P) == 1
    assert get_completed_fertility(0, 3, P) == 2
    assert get_completed_fertility(2, 2, P) == 1  # unreachable (nn=2, K+1)


def test_child_bin_high_cutoff() -> None:
    dep_last = 1
    assert current_child_bin_dt(1, 1, dep_last) == 3
    assert current_child_bin_dt(2, 1, dep_last) == 4  # default cutoff 2
    assert current_child_bin_dt(2, 1, dep_last, high_cutoff=3) == 3  # 1-2 bin
    assert current_child_bin_dt(3, 1, dep_last, high_cutoff=3) == 4  # 3+ bin
    assert current_child_bin_dt(0, 0, dep_last, high_cutoff=3) == 2


def test_l4_solves_and_preserves_age_mass() -> None:
    l3, _, _ = run_fork(_seq_eqscale(), verbose=False)
    l4, _, _ = run_fork(_seq_eqscale(**L4), verbose=False)
    assert l4.g.shape[5] == 4
    m3 = np.sum(l3.g, axis=(0, 1, 2, 4, 5, 6))
    m4 = np.sum(l4.g, axis=(0, 1, 2, 4, 5, 6))
    np.testing.assert_allclose(m4, m3, rtol=0.0, atol=1e-12)
    parity_mass = np.sum(l4.g, axis=(0, 1, 2, 3, 4, 6))
    assert np.all(parity_mass >= -1e-15)
    assert parity_mass[3] > 0.0  # third births actually happen


def test_l4_no_same_period_chaining() -> None:
    _, P, _ = run_fork(_seq_eqscale(**L4), verbose=False)
    assert float(P._second_births_by_age[0]) == 0.0
    assert float(P._third_births_by_age[0]) == 0.0
    assert float(P._third_births_by_age[1]) == 0.0
    assert float(np.sum(P._third_births_by_age)) > 0.0


def test_l4_third_birth_hazard_respects_fecundity() -> None:
    sol, _, _ = run_fork(
        _seq_eqscale(**L4, fecundity_omega1=0.5, fecundity_omega2=0.0),
        verbose=False,
    )
    np.testing.assert_allclose(
        sol.third_birth_hazard_by_age,
        0.5 * sol.third_attempt_hazard_by_age,
        atol=1e-12,
        rtol=0.0,
    )
    assert 0.0 <= float(sol.parity_progression_2to3_flow) <= 1.0


def test_entrant_conversion_factor_scales_accounting() -> None:
    a, _, _ = run_fork(_seq_eqscale(), verbose=False)
    b, _, _ = run_fork(_seq_eqscale(entrant_conversion_factor=0.5), verbose=False)
    # Fixed-price PE with fixed entry shares: the factor is accounting-only.
    assert np.array_equal(a.g, b.g)
    np.testing.assert_allclose(
        float(b.entrants_mature_total),
        0.5 * float(a.entrants_mature_total),
        rtol=1e-12,
        atol=0.0,
    )


def test_literal_topcode_tfr_arithmetic() -> None:
    parity = np.array([0.2, 0.2, 0.3, 0.3])
    expected = 0.2 * 1.0 + 0.3 * 2.0 + 0.3 * 3.4
    np.testing.assert_allclose(literal_topcode_tfr(parity, 3.4), expected, rtol=0.0, atol=1e-15)
    # Degenerate two-state input still works.
    np.testing.assert_allclose(literal_topcode_tfr(np.array([0.5, 0.5]), 1.0), 0.5, rtol=0.0, atol=1e-15)


def test_extract_moments_tfr_convention() -> None:
    l4, P4, _ = run_fork(_seq_eqscale(**L4), verbose=False)
    m_lit = extract_moments(l4, P4)
    parity = np.asarray(l4.parity_dist, dtype=float)
    np.testing.assert_allclose(
        m_lit["tfr"], literal_topcode_tfr(parity, 3.4), rtol=0.0, atol=1e-12
    )
    l3, P3, _ = run_fork(_seq_eqscale(), verbose=False)
    m_leg = extract_moments(l3, P3)
    np.testing.assert_allclose(m_leg["tfr"], 2.0 * float(l3.mean_parity), rtol=0.0, atol=1e-15)
