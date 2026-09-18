"""Tests for the default-off size-dependent rental wedge switch.

Rent per room is r + w0 + w1 * max(0, h - hk): a real operating cost paid to
the outside landlord. It enters the renter budget only (not the
property-tax base, not the rebate, not wealth) and leaves the user cost
itself unchanged.
"""

from __future__ import annotations

import numpy as np
import pytest

from intergen_eqscale_seq_optimized.kernels import renter_wedge_flow
from intergen_eqscale_seq_optimized.parameters import (
    apply_overrides,
    rental_wedge_total_cost,
    setup_parameters,
)
from intergen_eqscale_seq_optimized.solver import (
    renter_wedge_flow_py,
    run_model_cp_dt,
)
from intergen_eqscale_seq_optimized.tests.test_eqscale_seq import _tiny_markov


WEDGE = {"rental_wedge_intercept": 0.02, "rental_wedge_slope": 0.05, "rental_wedge_knee": 6.0}


def test_off_bitwise_identical() -> None:
    """Zeros leave every solved array bit for bit unchanged."""
    base, _, _ = run_model_cp_dt(_tiny_markov(), verbose=False)
    explicit, _, _ = run_model_cp_dt(
        {**_tiny_markov(), **WEDGE, "rental_wedge_intercept": 0.0, "rental_wedge_slope": 0.0},
        verbose=False,
    )
    assert np.array_equal(base.V, explicit.V)
    assert np.array_equal(base.g, explicit.g)
    assert np.array_equal(base.c_pol, explicit.c_pol)
    assert np.array_equal(base.hR_pol, explicit.hR_pol)
    assert np.array_equal(base.bp_pol, explicit.bp_pol)


def test_cost_of_eight_rooms() -> None:
    """With w0 = 0.02, w1 = 0.05, hk = 6, eight rooms cost 8 * (r + 0.12)."""
    P = apply_overrides(setup_parameters(), {**_tiny_markov(), **WEDGE})
    r = 0.2
    assert rental_wedge_total_cost(8.0, r, P) == pytest.approx(8.0 * (r + 0.02 + 0.10))
    assert rental_wedge_total_cost(4.0, r, P) == pytest.approx(4.0 * (r + 0.02))


def test_numba_python_mirror_agree() -> None:
    """The numba intratemporal helper matches its Python mirror."""
    numba_flow = getattr(renter_wedge_flow, "py_func", renter_wedge_flow)
    rng = np.random.default_rng(3)
    S_raw = rng.uniform(0.5, 6.0, size=25)
    hbc = rng.uniform(0.1, 1.0, size=25)
    args = dict(ri=0.2, w0=0.02, w1=0.05, hk=6.0, hRmax=11.0, al=0.7, oms=-1.0, es=1.0)
    u_py, ct_py, ht_py = renter_wedge_flow_py(S_raw, 0.4, hbc, args["ri"], 0.02, 0.05, 6.0, 11.0, 0.7, -1.0, 1.0)
    for k in range(25):
        u_nb, ct_nb, ht_nb = numba_flow(
            float(S_raw[k]), 0.4, float(hbc[k]), args["ri"], 0.02, 0.05, 6.0, 11.0, 0.7, -1.0, 1.0
        )
        assert u_nb == pytest.approx(u_py[k], rel=1e-12)
        assert ct_nb == pytest.approx(ct_py[k], rel=1e-12)
        assert ht_nb == pytest.approx(ht_py[k], rel=1e-12)


def test_optimal_renter_housing_weakly_lower() -> None:
    """The wedge (weakly) lowers rented rooms: intratemporally at every
    (surplus, committed-rooms) point, and in equilibrium aggregate renter
    demand. (A literal policy-function comparison can show higher h at a few
    young-saver states: the wedge destroys the motive to save for a down
    payment, so those households save much less and spend more today. The
    substitution margin itself is monotone everywhere, as checked here.)"""
    ri = 0.2875
    for surplus in np.linspace(0.3, 10.0, 40):
        for hbc in (0.1, 0.25, 0.5, 1.0, 2.0):
            _, _, ht_base = renter_wedge_flow_py(
                surplus, 0.05, hbc, ri, 0.0, 0.0, 6.0, 6.0, 0.7, -1.0, 1.0
            )
            _, _, ht_wedged = renter_wedge_flow_py(
                surplus, 0.05, hbc, ri, 0.02, 0.05, 6.0, 6.0, 0.7, -1.0, 1.0
            )
            assert float(ht_wedged) <= float(ht_base) + 1e-9
    base, _, _ = run_model_cp_dt(_tiny_markov(), verbose=False)
    wedged, _, _ = run_model_cp_dt({**_tiny_markov(), **WEDGE}, verbose=False)
    demand_base = float(np.sum(base.g[:, 0] * base.hR_pol[:, 0]))
    demand_wedged = float(np.sum(wedged.g[:, 0] * wedged.hR_pol[:, 0]))
    assert demand_wedged <= demand_base + 1e-9
    assert np.all(np.isfinite(wedged.V))


def test_python_branch_matches() -> None:
    """The non-compiled renter path honors the wedge and conserves mass."""
    sol, _, _ = run_model_cp_dt(
        {**_tiny_markov(), **WEDGE, "use_full_kernel": False}, verbose=False
    )
    assert np.all(np.isfinite(sol.V))
    assert float(np.sum(sol.g)) == pytest.approx(1.0, abs=1e-9)
    base, _, _ = run_model_cp_dt(
        {**_tiny_markov(), "use_full_kernel": False}, verbose=False
    )
    demand_base = float(np.sum(base.g[:, 0] * base.hR_pol[:, 0]))
    demand_wedged = float(np.sum(sol.g[:, 0] * sol.hR_pol[:, 0]))
    assert demand_wedged <= demand_base + 1e-9


def test_tiny_smoke_both_modes() -> None:
    """The tiny smoke solves with the wedge on; mass is conserved."""
    sol, _, _ = run_model_cp_dt({**_tiny_markov(), **WEDGE}, verbose=False)
    assert np.all(np.isfinite(sol.V))
    assert float(np.sum(sol.g)) == pytest.approx(1.0, abs=1e-9)


def test_wedge_rejects_negative() -> None:
    with pytest.raises(ValueError):
        apply_overrides(setup_parameters(), {"rental_wedge_intercept": -0.01})
    with pytest.raises(ValueError):
        apply_overrides(setup_parameters(), {"rental_wedge_slope": -0.01})
