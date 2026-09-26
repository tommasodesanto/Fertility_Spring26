"""Sparse, adjoint wealth-grid transport for a one-period estate receipt.

The operator applies a liquid receipt at the *start* of a receiving period:
``b -> b + X``.  It does not alter income.  For a branch with probability
``p``, the shifted value is linearly interpolated between its neighbouring
wealth-grid nodes.  Thus ``backward_expectation`` is ``Q @ values`` and
``forward_transport`` is exactly ``Q.T @ mass`` for the same implicit,
row-stochastic matrix ``Q``.

At a grid boundary the numerical operator clips shifted wealth to the endpoint.
This is an interpolation approximation, not a value bound.  The plan retains
the true-versus-clipped difference for every source node and branch; forward
transport rejects occupied clipping by default rather than concealing it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class EstateReceiptJumpPlan:
    """Interpolation indices and audit data for an IID receipt lottery.

    Arrays with shape ``(n_branches, n_wealth)`` index source wealth nodes by
    their second axis. ``true_minus_clipped`` is positive at upper clipping and
    negative at lower clipping (the latter cannot arise with nonnegative
    receipts from a source grid node, but is recorded explicitly).
    """

    wealth_grid: np.ndarray
    amounts: np.ndarray
    probabilities: np.ndarray
    lower_index: np.ndarray
    upper_index: np.ndarray
    upper_weight: np.ndarray
    shifted_wealth: np.ndarray
    clipped_wealth: np.ndarray
    lower_clipped: np.ndarray
    upper_clipped: np.ndarray
    true_minus_clipped: np.ndarray
    is_identity: bool


class OccupiedClippingError(ValueError):
    """Raised when occupied transport exceeds the caller's clipping allowance."""

    def __init__(self, message: str, account: dict[str, float]) -> None:
        super().__init__(message)
        self.account = account


def build_estate_receipt_jump_plan(
    wealth_grid: Any, amounts: Any, probabilities: Any
) -> EstateReceiptJumpPlan:
    """Build a plan from a finite increasing grid and a finite receipt lottery."""
    grid = np.asarray(wealth_grid, dtype=float)
    receipt_amounts = np.asarray(amounts, dtype=float)
    receipt_probabilities = np.asarray(probabilities, dtype=float)
    if grid.ndim != 1 or grid.size < 2 or not np.all(np.isfinite(grid)):
        raise ValueError("wealth_grid must be a finite one-dimensional grid with at least two nodes")
    if not np.all(np.diff(grid) > 0.0):
        raise ValueError("wealth_grid must be strictly increasing")
    if receipt_amounts.ndim != 1 or receipt_probabilities.ndim != 1 or receipt_amounts.size == 0:
        raise ValueError("amounts and probabilities must be nonempty one-dimensional arrays")
    if receipt_amounts.shape != receipt_probabilities.shape:
        raise ValueError("amounts and probabilities must have the same shape")
    if not np.all(np.isfinite(receipt_amounts)) or np.any(receipt_amounts < 0.0):
        raise ValueError("receipt amounts must be finite and nonnegative")
    if not np.all(np.isfinite(receipt_probabilities)) or np.any(receipt_probabilities < 0.0):
        raise ValueError("receipt probabilities must be finite and nonnegative")
    if not np.isclose(receipt_probabilities.sum(), 1.0, rtol=1e-12, atol=1e-12):
        raise ValueError("receipt probabilities must sum to one")

    shifted = grid[None, :] + receipt_amounts[:, None]
    clipped = np.clip(shifted, grid[0], grid[-1])
    lower_clipped = shifted < grid[0]
    upper_clipped = shifted > grid[-1]
    upper_index = np.searchsorted(grid, clipped, side="left")
    upper_index = np.clip(upper_index, 0, grid.size - 1)
    lower_index = np.maximum(upper_index - 1, 0)
    exact_node = grid[upper_index] == clipped
    lower_index = np.where(exact_node, upper_index, lower_index)
    spacing = grid[upper_index] - grid[lower_index]
    upper_weight = np.zeros_like(clipped)
    non_node = spacing > 0.0
    upper_weight[non_node] = (
        (clipped[non_node] - grid[lower_index[non_node]]) / spacing[non_node]
    )

    return EstateReceiptJumpPlan(
        wealth_grid=grid,
        amounts=receipt_amounts,
        probabilities=receipt_probabilities,
        lower_index=lower_index,
        upper_index=upper_index,
        upper_weight=upper_weight,
        shifted_wealth=shifted,
        clipped_wealth=clipped,
        lower_clipped=lower_clipped,
        upper_clipped=upper_clipped,
        true_minus_clipped=shifted - clipped,
        is_identity=bool(np.all(receipt_amounts == 0.0)),
    )


def _validate_leading_dimension(array: Any, n_wealth: int, name: str) -> np.ndarray:
    result = np.asarray(array)
    if result.ndim < 1 or result.shape[0] != n_wealth:
        raise ValueError(f"{name} must have wealth as leading axis of length {n_wealth}")
    if not np.issubdtype(result.dtype, np.number):
        raise TypeError(f"{name} must have a numeric dtype")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be finite")
    return result


def backward_expectation(plan: EstateReceiptJumpPlan, values: Any) -> np.ndarray:
    """Apply the receipt lottery to values, with endpoint clipping if needed."""
    values_array = _validate_leading_dimension(values, plan.wealth_grid.size, "values")
    if plan.is_identity:
        return values_array.copy()
    output = np.zeros_like(values_array, dtype=np.result_type(values_array.dtype, float))
    tail = (None,) * (values_array.ndim - 1)
    for branch, probability in enumerate(plan.probabilities):
        weight = plan.upper_weight[branch][(slice(None),) + tail]
        interpolated = (
            (1.0 - weight) * values_array[plan.lower_index[branch]]
            + weight * values_array[plan.upper_index[branch]]
        )
        output += probability * interpolated
    return output


def forward_transport(
    plan: EstateReceiptJumpPlan, mass: Any, max_clipped_wealth: float = 0.0
) -> tuple[np.ndarray, dict[str, float]]:
    """Transport mass by the adjoint interpolation operator.

    ``max_clipped_wealth`` is the maximum allowed *occupied expected wealth
    loss* from endpoint clipping.  It defaults to zero, so any occupied
    overflow fails.  The exception has an ``account`` attribute containing the
    transport audit.  Unoccupied boundary nodes do not trigger rejection.
    """
    if np.isnan(max_clipped_wealth) or max_clipped_wealth < 0.0:
        raise ValueError("max_clipped_wealth must be nonnegative")
    mass_array = _validate_leading_dimension(mass, plan.wealth_grid.size, "mass")
    if np.any(mass_array < 0.0):
        raise ValueError("mass must be nonnegative")
    if plan.is_identity:
        total_mass = float(mass_array.sum())
        wealth_before = float(np.tensordot(plan.wealth_grid, mass_array.sum(axis=tuple(range(1, mass_array.ndim)), dtype=float) if mass_array.ndim > 1 else mass_array, axes=1))
        return mass_array.copy(), {
            "transported_mass": total_mass,
            "expected_receipt_flow": 0.0,
            "clipped_wealth_loss": 0.0,
            "mass_encountering_clipping": 0.0,
            "wealth_before": wealth_before,
            "wealth_after": wealth_before,
        }

    result = np.zeros_like(mass_array, dtype=np.result_type(mass_array.dtype, float))
    axes_after_wealth = tuple(range(1, mass_array.ndim))
    mass_by_wealth = mass_array.sum(axis=axes_after_wealth, dtype=float) if axes_after_wealth else mass_array.astype(float)
    for branch, probability in enumerate(plan.probabilities):
        np.add.at(result, plan.lower_index[branch], probability * (1.0 - plan.upper_weight[branch])[(slice(None),) + (None,) * len(axes_after_wealth)] * mass_array)
        np.add.at(result, plan.upper_index[branch], probability * plan.upper_weight[branch][(slice(None),) + (None,) * len(axes_after_wealth)] * mass_array)

    expected_receipt_flow = float(mass_by_wealth.sum() * np.dot(plan.probabilities, plan.amounts))
    clipped_wealth_loss = float(np.sum(
        plan.probabilities[:, None] * mass_by_wealth[None, :] * plan.true_minus_clipped
    ))
    clipping_mass = float(np.sum(
        plan.probabilities[:, None]
        * mass_by_wealth[None, :]
        * (plan.lower_clipped | plan.upper_clipped)
    ))
    wealth_before = float(np.dot(plan.wealth_grid, mass_by_wealth))
    mass_after_by_wealth = result.sum(axis=axes_after_wealth, dtype=float) if axes_after_wealth else result
    wealth_after = float(np.dot(plan.wealth_grid, mass_after_by_wealth))
    account = {
        "transported_mass": float(result.sum()),
        "expected_receipt_flow": expected_receipt_flow,
        "clipped_wealth_loss": clipped_wealth_loss,
        "mass_encountering_clipping": clipping_mass,
        "wealth_before": wealth_before,
        "wealth_after": wealth_after,
    }
    if clipped_wealth_loss > max_clipped_wealth:
        raise OccupiedClippingError(
            "occupied endpoint clipping exceeds max_clipped_wealth", account
        )
    return result, account


def adjointness_residual(plan: EstateReceiptJumpPlan, values: Any, mass: Any) -> float:
    """Return ``<mass, backward(values)> - <forward(mass), values>`` for audit.

    This utility allows clipping because it audits the numerical endpoint
    convention itself; callers concerned with economic truncation should use
    ``forward_transport`` with its default zero clipping allowance.
    """
    transported, _ = forward_transport(plan, mass, max_clipped_wealth=np.inf)
    return float(np.vdot(np.asarray(mass), backward_expectation(plan, values)) - np.vdot(transported, np.asarray(values)))
