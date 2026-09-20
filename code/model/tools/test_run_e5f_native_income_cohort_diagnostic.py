from __future__ import annotations
import numpy as np
from types import SimpleNamespace
from run_e5f_native_income_cohort_diagnostic import synthetic_mass_flow, housing_stats


def test_synthetic_cohort_mass_flow():
    pre = np.array([1.0, .5]); post = np.array([.8, .4]); deaths = np.array([.1, .05]); entry = np.array([.1, .05])
    assert abs(synthetic_mass_flow(pre, post, deaths, entry) - .3) < 1e-12


def test_zero_flow_is_exact():
    x = np.ones((2, 3)); assert synthetic_mass_flow(x, x, np.zeros_like(x), np.zeros_like(x)) == 0.0


def test_invalid_negative_mass_rejected():
    try:
        synthetic_mass_flow(np.array([1.0]), np.array([-0.1]), np.zeros(1), np.zeros(1))
    except ValueError:
        return
    raise AssertionError("negative mass must be rejected")


def test_native_cohort_shape_and_housing_arithmetic():
    # Current tenure axis contains realized renter/owner mass.
    shape = (1, 3, 1, 2, 1, 1, 1)
    mass = np.zeros(shape); mass[:, 0, :, 0, ...] = .25
    mass[:, 1, :, 0, ...] = .5; mass[:, 2, :, 0, ...] = .25
    policy = SimpleNamespace(hR_pol=np.full(shape, 2.0))
    rooms, owner = housing_stats(SimpleNamespace(g_current=mass, policy=policy), SimpleNamespace(H_own=np.array([4.0, 8.0])))
    assert abs(rooms - 4.5) < 1e-12
    assert abs(owner - .75) < 1e-12
