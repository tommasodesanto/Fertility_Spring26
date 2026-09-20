import numpy as np

from types import SimpleNamespace
import pytest
from analyze_e5f_native_financing_groups import _first, _housing, _weighted_quantiles, group_definition


def test_weighted_quantiles_and_duplicate_cuts_are_explicit():
    cuts = _weighted_quantiles(np.array([1., 2., 3.]), np.array([1., 0., 3.]))
    assert cuts.tolist() == [1.0, 3.0, 3.0]


def test_groups_are_disjoint_and_exhaustive_with_owner_products():
    # (wealth, tenure, location, age, m, n, shock)
    g = np.zeros((4, 3, 1, 3, 1, 1, 1))
    g[:, :, :, :, :, :, :] = 1.0
    groups, meta = group_definition(g, np.arange(4.), np.array([0, 1, 2]), np.array([1]))
    eligible = np.zeros_like(g, dtype=bool); eligible[:, :, :, :, :, 0, :] = True
    union = np.zeros_like(eligible)
    for x in groups:
        assert not np.any(union & x["mask"])
        union |= x["mask"]
    assert np.array_equal(union, eligible)
    assert meta["eligible_mass"] == 36.0


def test_first_birth_loss_is_exact_and_nonnegative():
    pre = np.zeros((2, 2, 1, 2, 1, 1, 1)); post = pre.copy()
    pre[..., 0, 0, :] = 0.4; post[..., 0, 0, :] = 0.1
    assert np.isclose(_first(pre, post, np.array([0, 1])), 2.4)


def test_first_birth_increase_is_an_error():
    pre = np.ones((1, 1, 1, 1, 1, 1, 1)); post = pre + 1e-4
    with pytest.raises(ValueError, match="mass increases"):
        _first(pre, post, np.array([0]))


def test_housing_hand_mapping_reconciles_rooms_and_six_room_masses():
    gc = np.zeros((1, 2, 1, 1, 1, 1, 1)); gc[:, 0] = 2.; gc[:, 1] = 3.
    h = np.zeros_like(gc); h[:, 0] = 6.
    out = _housing(gc, h, SimpleNamespace(H_own=np.array([7.])))
    assert out["renter_mass"] == 2 and out["owner_mass"] == 3
    assert out["rooms_total"] == 33 and out["mass_renting_ge6_rooms"] == 2 and out["mass_owning_ge6_rooms"] == 3
