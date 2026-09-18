"""Tests for the default-off mortgage block switches.

Two independent switches, both off by default. Origination-only applies the
collateral floor solely on transactions (where the down-payment test fires);
a stayer faces a no-cash-out rule instead (floor b when in debt, the
collateral floor when not, with the taper/line rollover bypassed). The
amortization rate forces a stayer in debt to retire at least that share per
period; combined with origination-only the stayer floor is the max of the
two, and with amortization alone it is the max of the amortization floor
and the standard floor.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from intergen_eqscale_seq_optimized.kernels import full_owner_block_kernel
from intergen_eqscale_seq_optimized.parameters import apply_overrides, setup_parameters
from intergen_eqscale_seq_optimized.solver import owner_borrowing_floor, run_model_cp_dt
from intergen_eqscale_seq_optimized.tests.test_eqscale_seq import _tiny_markov


GS_ALPHA1 = (3.0 - math.sqrt(5.0)) / 2.0
GS_ALPHA2 = (math.sqrt(5.0) - 1.0) / 2.0
DEBT_GRID = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])


def _tiny_debt(**extra: object) -> dict[str, object]:
    base = {**_tiny_markov(), "b_min": -2.0, "b_core_lo": -2.0}
    base.update(extra)
    return base


def _owner_kernel_call(
    grid: np.ndarray,
    income: float,
    collateral_floor: float,
    *,
    stay_on: int = 0,
    stay_orig: int = 0,
    amort: float = 0.0,
):
    kernel = getattr(full_owner_block_kernel, "py_func", full_owner_block_kernel)
    nb = grid.size
    resources = grid + income
    return kernel(
        np.ascontiguousarray(resources),
        np.ascontiguousarray(resources),
        np.zeros((nb, 1)),
        np.zeros((nb, 1)),
        0,
        np.ascontiguousarray(grid),
        np.array([0.4]),
        np.array([0.5]),
        np.array([0.0]),
        np.array([0.0]),
        np.array([0.7]),
        np.array([1.0]),
        np.array([collateral_floor]),
        0.1,
        2.0,
        1.0,
        1.0,
        0.04,
        0.7,
        -1.0,
        0.9,
        1.0,
        0.0,
        GS_ALPHA1,
        GS_ALPHA2,
        1e-8,
        0,
        0,
        np.zeros(1),
        0,
        stay_on,
        stay_orig,
        amort,
    )


def test_off_bitwise_identical() -> None:
    """Both switches off leave every solved array bit for bit unchanged."""
    base, _, _ = run_model_cp_dt(_tiny_debt(), verbose=False)
    explicit, _, _ = run_model_cp_dt(
        {**_tiny_debt(), "mortgage_origination_only": False, "mortgage_amortization": 0.0},
        verbose=False,
    )
    assert np.array_equal(base.V, explicit.V)
    assert np.array_equal(base.g, explicit.g)
    assert np.array_equal(base.c_pol, explicit.c_pol)
    assert np.array_equal(base.hR_pol, explicit.hR_pol)
    assert np.array_equal(base.bp_pol, explicit.bp_pol)
    assert base.bp_pol_stay is None
    assert explicit.bp_pol_stay is None


def test_origination_only_stayer_cannot_raise_debt() -> None:
    """With (a) on, a stayer in debt cannot choose b' below its balance."""
    # Near-zero income forces both regimes onto their floors, isolating them.
    _, bp_off, _ = _owner_kernel_call(DEBT_GRID, 0.01, -8.0)
    _, bp_on, _ = _owner_kernel_call(DEBT_GRID, 0.01, -8.0, stay_on=1, stay_orig=1)
    indebted = DEBT_GRID < 0.0
    assert np.all(bp_on[indebted, 0] >= DEBT_GRID[indebted] - 1e-12)
    # Without the switch the same stayer can borrow against collateral.
    assert bool(np.any(bp_off[indebted, 0] < DEBT_GRID[indebted] - 1e-9))
    # With no debt a first mortgage on the owned house is an origination.
    nodebt = DEBT_GRID >= 0.0
    assert np.all(bp_on[nodebt, 0] >= -8.0 - 1e-12)


def test_amortization_floor() -> None:
    """With (b) at 0.11, a stayer with b = -1 cannot choose b' below -0.89."""
    _, bp, _ = _owner_kernel_call(
        DEBT_GRID, 0.01, -8.0, stay_on=1, stay_orig=0, amort=0.11
    )
    assert float(bp[0, 0]) >= -0.89 - 1e-12
    assert float(bp[0, 0]) == pytest.approx(-0.89)  # the floor binds above the balance
    floor = owner_borrowing_floor(
        apply_overrides(setup_parameters(), _tiny_debt()),
        np.array([-1.0]),
        np.array([-8.0]),
        0,
        stay_on=True,
        stay_orig=True,
        amort=0.11,
    )
    assert float(floor[0]) == pytest.approx(-0.89)


def test_buyer_feasible_set_unchanged() -> None:
    """A buyer's feasible set never reads stayer values: with Vd_stay far
    from Vd, origin-renter tenure rows match bit for bit."""
    from intergen_eqscale_seq_optimized.kernels import tenure_choice_kernel

    kernel = getattr(tenure_choice_kernel, "py_func", tenure_choice_kernel)
    nb = 6
    b_grid = np.linspace(-1.0, 3.0, nb)
    nt, npar, ncs = 3, 3, 4
    rng = np.random.default_rng(0)
    Vd = rng.random((nb, nt, 1, npar, ncs))
    Vd_stay = Vd - 5.0
    heq = np.zeros((1, nt))
    heq[:, 1:] = 0.5
    hcost = np.zeros((1, nt))
    hcost[:, 1:] = 2.0
    dp = np.full((1, nt, npar, ncs), 0.5)
    bmo = np.full((1, nt, npar, ncs), -8.0)
    birth_dp = np.zeros((npar, ncs, nt, nt), dtype=bool)
    grant = np.zeros((1, nt, npar, ncs))
    VH1, tc1 = kernel(Vd, b_grid, heq, hcost, dp, bmo, birth_dp, grant, Vd)
    VH2, tc2 = kernel(Vd, b_grid, heq, hcost, dp, bmo, birth_dp, grant, Vd_stay)
    assert np.array_equal(tc1[:, 0], tc2[:, 0])
    assert np.array_equal(VH1[:, 0], VH2[:, 0])


def test_tiny_smoke_both_modes() -> None:
    """The tiny smoke solves with each switch on; mass is conserved."""
    for overrides in (
        {"mortgage_origination_only": True},
        {"mortgage_amortization": 0.11},
        {"mortgage_origination_only": True, "mortgage_amortization": 0.11},
    ):
        sol, _, _ = run_model_cp_dt({**_tiny_debt(), **overrides}, verbose=False)
        assert np.all(np.isfinite(sol.V))
        assert float(np.sum(sol.g)) == pytest.approx(1.0, abs=1e-9)
        assert sol.bp_pol_stay is not None
        assert float(np.sum(sol.housing_demand)) > 0.0
    off, _, _ = run_model_cp_dt(_tiny_debt(), verbose=False)
    on, _, _ = run_model_cp_dt(
        {**_tiny_debt(), "mortgage_origination_only": True}, verbose=False
    )
    assert not np.array_equal(off.V, on.V)


def test_amortization_rejects_bad_values() -> None:
    with pytest.raises(ValueError):
        apply_overrides(setup_parameters(), {"mortgage_amortization": 1.0})
    with pytest.raises(ValueError):
        apply_overrides(setup_parameters(), {"mortgage_amortization": -0.1})
