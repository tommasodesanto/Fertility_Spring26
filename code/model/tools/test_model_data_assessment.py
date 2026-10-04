"""Focused measurement checks for the saved-case comparison."""
import numpy as np
import pytest

from model_data_assessment import (financial_positions, invert_independent_birth_map,
                                   require_sha256, select_grouped_national_acs,
                                   thirds, weighted_ecdf, within_cell_stock)


def test_ecdf_ties_and_negative_values():
    x, f = weighted_ecdf([-2, 1, -2], [1, 3, 2])
    np.testing.assert_allclose(x, [-2, 1])
    np.testing.assert_allclose(f, [.5, 1])


def test_within_cell_projection_uses_age_exposure():
    # Both ages fall in the 20–23 cell; the post stock of the next cell is irrelevant.
    got = within_cell_stock([0, 100], [4, 200], [20, 24], [20, 23], [1, 3], 4)
    assert got == pytest.approx(2.75)


def test_first_birth_flow_map():
    post = np.zeros((1, 1, 1, 1, 1, 4, 4))
    post[..., 0, 0], post[..., 1, 1] = .8, .2
    first = np.zeros((1, 1, 1, 1, 1, 4)); first[..., 1] = .2
    later = np.zeros((1, 1, 1, 1, 1, 2, 2, 4))
    pre, rate = invert_independent_birth_map(post, first, later, [1.], range(1))
    assert pre[..., 0, 0].item() == pytest.approx(1.)
    assert (pre * rate)[..., 0, 0].sum() == pytest.approx(.2)


def test_mortgage_included_in_financial_position():
    position, total = financial_positions([100], [80])
    assert position.item() == 20
    assert total.item() == 100


def test_fractional_tied_terciles():
    allocation = thirds([2, 1])
    np.testing.assert_allclose((allocation * [2, 1]).sum(axis=1), [1, 1, 1])
    np.testing.assert_allclose(allocation[:, 0], [.5, .5, 0])


def test_bad_source_hash_rejected(tmp_path):
    p = tmp_path / 'source'; p.write_text('source')
    with pytest.raises(RuntimeError, match='SHA-256 mismatch'):
        require_sha256(p, '0' * 64, 'source')


def test_annual_national_acs_only():
    base = {'geography_scope': 'national', 'sample': 'all_structures', 'age_kind': 'annual'}
    rows = [base, dict(base, age_kind='four_year'), dict(base, sample='DUE')]
    assert select_grouped_national_acs(rows) == [base]
