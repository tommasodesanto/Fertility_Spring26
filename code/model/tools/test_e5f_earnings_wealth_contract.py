"""Small contract tests for the frozen-source earnings/purchase adapter."""
from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[3]
FROZEN_MODEL = ROOT / "tmp/income_refinement_local_20260921/source/code/model"
if str(FROZEN_MODEL) not in sys.path:
    sys.path.insert(0, str(FROZEN_MODEL))

TOOLS = ROOT / "code/model/tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))
contract = importlib.import_module("e5f_earnings_wealth_contract")
solver = importlib.import_module("intergen_eqscale_seq_optimized.solver")
kernels = importlib.import_module("intergen_eqscale_seq_optimized.kernels")


def _reference_model():
    def entry(b_grid, P, *, i=0, j=0, z_value=1.0):
        if z_value < 1.0:
            return np.array([0, 1]), np.array([0.7, 0.3])
        return np.array([0, 1]), np.array([0.2, 0.8])

    def income_transition_values(P):
        return np.array([0.5, 1.5]), np.array([0.4, 0.6]), np.eye(2)

    return SimpleNamespace(
        I=1,
        z_grid=np.array([0.5, 1.5]),
        entry_wealth_grid_weights=entry,
        income_transition_values=income_transition_values,
    )


def test_rank_coupled_entry_preserves_wealth_and_iid_transitory_independence():
    model = _reference_model()
    grid = np.array([0.0, 1.0])
    conditional, receipt = contract.rank_coupled_entry(
        model, model, grid, np.array([0.25, 0.75]), np.array([0.3, 0.7])
    )
    assert conditional.shape == (2, 4)
    np.testing.assert_allclose(
        receipt["candidate_wealth_marginal"], receipt["reference_wealth_marginal"], atol=2e-14
    )
    np.testing.assert_allclose(conditional[:, 0], conditional[:, 1])
    np.testing.assert_allclose(conditional[:, 2], conditional[:, 3])


def test_fixed_entry_requires_exact_reference_grid_and_income_node():
    model = _reference_model()
    contract.install_fixed_entry(model)
    P = SimpleNamespace(
        fixed_reference_entry_conditional=np.array([[0.6, 0.2], [0.4, 0.8]]),
        fixed_reference_entry_grid=np.array([0.0, 1.0]),
        z_grid=np.array([0.5, 1.5]),
    )
    idx, wt = model.entry_wealth_grid_weights(P.fixed_reference_entry_grid, P, z_value=1.5)
    np.testing.assert_array_equal(idx, [0, 1])
    np.testing.assert_allclose(wt, [0.2, 0.8])
    with pytest.raises(AssertionError):
        model.entry_wealth_grid_weights(np.array([0.0, 2.0]), P, z_value=1.5)
    with pytest.raises(ValueError):
        model.entry_wealth_grid_weights(P.fixed_reference_entry_grid, P, z_value=1.25)


def test_install_purchase_income_compiles_frozen_markov_solver(tmp_path):
    # This exercises only the source-anchor rewrite; it does not call the model.
    original = solver.solve_bellman_full_markov_income
    maps = solver.build_forward_tenure_transition_maps
    tenure = solver.tenure_choice_kernel, solver.tenure_logit_kernel
    try:
        diff = contract.install_purchase_income(solver, tmp_path / "purchase_income.diff")
    finally:
        solver.solve_bellman_full_markov_income = original
        solver.build_forward_tenure_transition_maps = maps
        solver.tenure_choice_kernel, solver.tenure_logit_kernel = tenure
    assert "diagnostic/solve_bellman_full_markov_income" in diff
    assert "dp_choice = dp_arr - income_for_purchase" in diff
    assert (tmp_path / "purchase_income.diff").exists()
    assert (tmp_path / "purchase_income.tenure.py").exists()


@pytest.mark.parametrize("boundary", ["upper_sale", "lower_purchase", "exact_endpoint"])
def test_transaction_support_applies_to_argmax_and_logit(tmp_path, boundary):
    grid = np.array([-1., 0., 1.])
    values = np.zeros((3, 2, 1, 1, 1))
    cost = np.array([[0., .3]])
    sale = cost.copy()
    # Nonbinding economic thresholds isolate the numerical support condition.
    dp = np.full((1, 2, 1, 1), -10.)
    floor = dp.copy()
    birth_dp = np.zeros((1, 1, 2, 2), dtype=bool)
    grants = np.zeros((1, 2, 1, 1))
    if boundary == "upper_sale":
        values[:, 0] = 10.
        b, old, desired = 2, 1, 0
    else:
        values[:, 1] = 10.
        b, old, desired = 0, 0, 1
        if boundary == "exact_endpoint":
            cost[:, 1] = 1.
            b = 1  # b=0 minus price=1 reaches lower endpoint exactly.
    args = (values, grid, sale, cost, dp, floor, birth_dp, grants)
    original = solver.tenure_choice_kernel, solver.tenure_logit_kernel
    reference_choice = original[0](*args)[1]
    assert reference_choice[b, old, 0, 0, 0] == desired
    try:
        contract._install_transaction_support(solver, tmp_path / "support.diff")
        deterministic = solver.tenure_choice_kernel(*args)
        probabilistic = solver.tenure_logit_kernel(*args, .1)
    finally:
        solver.tenure_choice_kernel, solver.tenure_logit_kernel = original
    if boundary == "exact_endpoint":
        assert deterministic[1][b, old, 0, 0, 0] == desired
        assert probabilistic[2][b, old, 0, 0, 0, desired] > .999
    else:
        assert deterministic[1][b, old, 0, 0, 0] != desired
        assert probabilistic[2][b, old, 0, 0, 0, desired] == 0.
    # Every branch at b=0 stays on-grid in the .3-price examples, so the
    # correction must preserve both conditional values and probabilities there.
    if boundary != "exact_endpoint":
        native_logit = original[1](*args, .1)
        for got, expected in zip(probabilistic, native_logit):
            np.testing.assert_array_equal(got[1], expected[1])
        native_deterministic = original[0](*args)
        for got, expected in zip(deterministic, native_deterministic):
            np.testing.assert_array_equal(got[1], expected[1])


def test_income_shift_changes_tenure_feasibility_threshold():
    b_grid = np.array([-1.0, -0.8, -0.6, 0.0, 0.8, 1.0])
    # Owner value dominates when feasible; hcost=.8, phi=.75 gives floor=-.6.
    Vd = np.zeros((len(b_grid), 2, 1, 1, 1))
    Vd[:, 1, 0, 0, 0] = 10.0
    heq = np.zeros((1, 2))
    hcost = np.array([[0.0, 0.8]])
    dp = np.array([[[[0.0]], [[0.2]]]])
    bmo = np.array([[[[0.0]], [[-0.6]]]])
    birth_dp = np.zeros((1, 1, 2, 2), dtype=bool)
    grants = np.zeros((1, 2, 1, 1))
    _, choice0 = kernels.tenure_choice_kernel(
        Vd, b_grid, heq, hcost, dp, bmo, birth_dp, grants
    )
    # At b=0, x=-.8: without income both tests fail.
    assert choice0[3, 0, 0, 0, 0] == 0
    dp_income = dp - 0.25
    bmo_income = np.maximum(bmo - 0.25, b_grid[0])
    _, choice1 = kernels.tenure_choice_kernel(
        Vd, b_grid, heq, hcost, dp_income, bmo_income, birth_dp, grants
    )
    assert choice1[3, 0, 0, 0, 0] == 1
    _, insufficient = kernels.tenure_choice_kernel(
        Vd, b_grid, heq, hcost, dp - 0.1, np.maximum(bmo - 0.1, b_grid[0]), birth_dp, grants)
    assert insufficient[3, 0, 0, 0, 0] == 0
    _, zero_shift = kernels.tenure_choice_kernel(
        Vd, b_grid, heq, hcost, dp - 0., np.maximum(bmo - 0., b_grid[0]), birth_dp, grants)
    np.testing.assert_array_equal(choice0, zero_shift)


def test_forward_map_preserves_supported_transaction_wealth_below_collateral_limit():
    P = SimpleNamespace()
    grid = np.array([-1.0, -0.2, 0.0, 1.0])
    hc = np.array([[0.0, 0.3]])
    he = np.zeros_like(hc)
    phi = np.full((1, 2, 1, 1), 0.8)
    birth_dp = np.zeros((1, 1, 2, 2), dtype=bool)
    grants = np.zeros((1, 2, 1, 1))
    idx, wt = contract._ordinary_forward_maps(P, grid, hc, he, phi, birth_dp, grants)
    i, w = idx[0, 0, 1, 0, 0, 2], wt[0, 0, 1, 0, 0, 2]
    implied = float((1-w) * grid[i] + w * grid[i+1])
    assert implied == pytest.approx(-0.3)


def _audit_fixture():
    grid = np.array([-1., -.2, 0., 1.])
    P = SimpleNamespace(I=1, J=1, z_grid=[1.], H_own=[.3], psi=0., R_gross=1.,
                        n_parity=1, n_child_states=1)
    shared = SimpleNamespace(phi_choice=np.full((1, 2, 1, 1), .8))
    idx, wt = contract._ordinary_forward_maps(P, grid, np.array([[0., .3]]),
        np.array([[0., .3]]), shared.phi_choice, np.zeros(1), np.zeros(1))
    shape = (4, 2, 1, 1, 1, 1, 1)
    post = np.zeros(shape)
    post[2, 0, 0, 0, 0, 0, 0] = 1.
    probs = np.zeros(shape + (2,))
    probs[..., 1] = 1.
    current = np.zeros(shape)
    i, w = idx[0, 0, 1, 0, 0, 2], wt[0, 0, 1, 0, 0, 2]
    current[i, 1, 0, 0, 0, 0, 0] = 1-w
    current[i+1, 1, 0, 0, 0, 0, 0] = w
    policy = SimpleNamespace(price=[1.], tenure_probs=probs, loc_probs=np.ones((4, 2, 1, 1, 1, 1, 1, 1)),
        bp_pol=np.zeros(shape), maps=SimpleNamespace(tmx_idx=idx, tmx_wt=wt))
    evaluation = SimpleNamespace(policy=policy, g_post_fertility=post, g_current=current)
    model = SimpleNamespace(income_at_state=lambda *args: .2)
    return evaluation, P, shared, grid, model


def test_accounting_accepts_income_purchase_with_unchanged_final_debt_limit():
    result = contract.audit_purchase_accounting(*_audit_fixture())
    assert result['purchase_threshold_violation_mass'] == 0
    assert result['end_mortgage_floor_violation_mass'] == 0


def test_native_owner_kernel_does_not_roll_over_purchase_shortfall():
    grid = np.array([-1., -.3, -.24, 0., 1.])
    zero = np.zeros(1)
    one = np.ones(1)
    continuation = np.zeros((len(grid), 1))
    # At transaction wealth -.3, income .5 leaves .2 resources. The final
    # mortgage floor is -.24, even though transaction wealth starts below it.
    resources = grid + .5
    value, bp, consumption = kernels.full_owner_block_kernel(
        resources, resources, continuation, continuation, 0, grid,
        zero, zero, zero, zero, .7*one, one, -.24*one,
        0., 1., 1., 1., 0., .7, -1., .96, 0., 0.,
        .3819660112501051, .6180339887498949, 1e-8, 0, 1)
    assert value[1, 0] > -1e9
    assert bp[1, 0] == pytest.approx(-.24)
    assert consumption[1, 0] + bp[1, 0] == pytest.approx(resources[1])


@pytest.mark.parametrize('defect', ['end_debt', 'threshold', 'wealth_creation', 'unsupported'])
def test_accounting_rejects_bad_occupied_branches(defect):
    evaluation, P, shared, grid, model = _audit_fixture()
    if defect == 'end_debt':
        evaluation.policy.bp_pol[:, 1] = -.4
    elif defect == 'threshold':
        model.income_at_state = lambda *args: 0.
    elif defect == 'wealth_creation':
        evaluation.policy.maps.tmx_wt[0, 0, 1, 0, 0, 2] += .1
    else:
        evaluation.g_post_fertility[:] = 0.
        evaluation.g_post_fertility[0, 0, 0, 0, 0, 0, 0] = 1.
    with pytest.raises(RuntimeError, match='purchase accounting gate failed'):
        contract.audit_purchase_accounting(evaluation, P, shared, grid, model)
