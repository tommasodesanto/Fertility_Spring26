"""Bitwise off-value tests for the sandbox's mechanism switches.

Run with: PYTHONPATH=.:sandbox code/model/.venv/bin/python -m pytest sandbox/tests -q
(from code/model/), or via `make sandbox-test`.

Each switch must reproduce the untouched package's precompute_shared() output
bitwise when left at its default ("linear" / "deflate"). This is checked
against a genuine (tiny-grid) live parameter object built the same way
run_ss.py builds one, not a hand-rolled stub, so the test exercises the real
predicate functions (independent_child_maturation_active, etc.).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

SANDBOX_ROOT = Path(__file__).resolve().parents[1]
MODEL_ROOT = SANDBOX_ROOT.parent
TOOLS_ROOT = MODEL_ROOT / "tools"
_PACKAGE_ROOT = os.environ.get("SANDBOX_PACKAGE_ROOT")
sys.path[:0] = [str(SANDBOX_ROOT), str(MODEL_ROOT), str(TOOLS_ROOT)]
if _PACKAGE_ROOT:
    # Test the sandbox's monkeypatch targets against an alternate package
    # snapshot (e.g. the September 13 corrected_initial source), inserted
    # ahead of MODEL_ROOT above so it wins import resolution. See
    # sandbox/README.md's package-version notice.
    sys.path.insert(0, str(Path(_PACKAGE_ROOT).resolve()))

import mechanisms  # noqa: E402
from intergen_eqscale_seq_optimized import solver as _solver  # noqa: E402


def _tiny_parameters():
    """A small live parameter object built the same way run_ss.py builds one."""
    import argparse
    import audit_closed_reproductive_closure as closure

    chain = closure.load_chain(profile="e5f-floor")
    overrides = closure.make_overrides(
        chain, {"psi_child": 0.15}, nb=5, profile="e5f-floor",
    )
    overrides["hbar_child_rooms"] = 0.0
    from intergen_eqscale_seq_optimized.parameters import apply_overrides, setup_parameters
    P = setup_parameters()
    P = apply_overrides(P, overrides)
    P.beta = 0.99 ** 4
    P.rho = 1 / P.beta - 1
    P.rho_hat = P.rho
    P.user_cost_rate = P.q + P.delta + P.tau_H
    P.phi = np.asarray(P.phi, dtype=float).reshape(-1)
    if P.phi.size == 1:
        P.phi = P.phi.item() * np.ones(P.n_parity)
    b_grid = _solver.make_grid(P)
    return P, b_grid


def test_child_benefit_form_off_is_bitwise_identical():
    P, b_grid = _tiny_parameters()
    P.child_benefit_form = "linear"
    P.scale_weighting = "deflate"
    expected = _solver.precompute_shared(P, b_grid)
    actual = mechanisms.sandboxed_precompute_shared(P, b_grid)
    assert actual is expected or np.array_equal(actual.psi_v, expected.psi_v)
    np.testing.assert_array_equal(actual.psi_v, expected.psi_v)
    np.testing.assert_array_equal(actual.c_bar, expected.c_bar)
    np.testing.assert_array_equal(actual.h_bar, expected.h_bar)
    np.testing.assert_array_equal(actual.escale_flat, expected.escale_flat)
    np.testing.assert_array_equal(actual.type_map, expected.type_map)


def test_scale_weighting_off_is_bitwise_identical():
    P, b_grid = _tiny_parameters()
    P.child_benefit_form = "linear"
    P.scale_weighting = "deflate"
    P.eqscale_form = "power"
    expected = _solver.precompute_shared(P, b_grid)
    actual = mechanisms.sandboxed_precompute_shared(P, b_grid)
    np.testing.assert_array_equal(actual.escale_flat, expected.escale_flat)


def test_context_manager_restores_original_function():
    original = _solver.precompute_shared
    with mechanisms.sandbox_context():
        assert _solver.precompute_shared is mechanisms.sandboxed_precompute_shared
    assert _solver.precompute_shared is original


def test_child_benefit_form_log_changes_psi_v_where_children_present():
    P, b_grid = _tiny_parameters()
    P.child_benefit_form = "log"
    P.scale_weighting = "deflate"
    baseline = _solver.precompute_shared(P, b_grid)
    with mechanisms.sandbox_context():
        modified = _solver.precompute_shared(P, b_grid)
    present = baseline.psi_v != 0.0
    assert present.any(), "test fixture has no child-having cells to compare"
    assert not np.array_equal(modified.psi_v[present], baseline.psi_v[present])
    # Childless cells are untouched.
    np.testing.assert_array_equal(modified.psi_v[~present], baseline.psi_v[~present])


def test_scale_weighting_multiply_changes_escale_under_power_form():
    P, b_grid = _tiny_parameters()
    P.child_benefit_form = "linear"
    P.scale_weighting = "multiply"
    P.eqscale_form = "power"
    baseline = _solver.precompute_shared(P, b_grid)
    with mechanisms.sandbox_context():
        modified = _solver.precompute_shared(P, b_grid)
    assert not np.array_equal(modified.escale_flat, baseline.escale_flat)


def test_child_earnings_penalty_is_rejected():
    import pytest
    with pytest.raises(NotImplementedError):
        mechanisms.check_switches_supported({"child_earnings_penalty": 0.10})
    # A zero/absent penalty is a no-op and must not raise.
    mechanisms.check_switches_supported({"child_earnings_penalty": 0.0})
    mechanisms.check_switches_supported({})
