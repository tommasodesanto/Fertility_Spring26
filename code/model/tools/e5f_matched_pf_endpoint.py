"""Explicit-input stationary endpoint evaluation for matched PF experiments.

No history reconstruction, policy solve, supply anchoring, or fiscal choice is
performed here. The supplied stationary policy must have been solved at the
declared price/transfer/preferences. A price-root caller must supply a NEW policy
at each trial; changing its price label is not a policy solve.

This reuses the existing finite-level household/person stationary mapping. It
does not solve a positive-growth balanced-growth path or select terminal data.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np


@dataclass(frozen=True)
class EndpointControls:
    """All numerical budgets and acceptance gates are caller-owned."""

    maximum_inner_iterations: int
    inner_damping: float
    distribution_tolerance: float
    birth_rate_tolerance: float
    one_step_tolerance: float
    market_tolerance: float
    fiscal_absolute_tolerance: float
    accounting_absolute_tolerance: float


@dataclass
class EndpointEvaluation:
    fixed_point: Any
    residuals: dict[str, float]
    root_coordinates: tuple[str, ...]
    root_residuals: np.ndarray
    gates: dict[str, bool]
    contract: dict[str, Any]

    @property
    def accepted(self) -> bool:
        return all(self.gates.values())

    @property
    def mapping_valid(self) -> bool:
        """Whether the inner evaluation is usable by an outer market root."""
        return all(value for name, value in self.gates.items()
                   if name not in {"housing_market", "fiscal_regime"})


def _positive(name: str, value: float) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _finite(name: str, value: float) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _solve_inner(**kwargs: Any) -> Any:
    # Lazy import keeps contract tests independent of model initialization.
    from run_e5f_perfect_foresight_person_demography import (
        solve_terminal_household_person_fixed_point,
    )
    return solve_terminal_household_person_fixed_point(**kwargs)


def evaluate_endpoint(
    *,
    parameters: Any,
    b_grid: np.ndarray,
    policy: Any,
    asset_price: float,
    transfer: float,
    psi_child: float,
    demographic_primitives: Any,
    initial_g_pre: np.ndarray,
    supply_rule: Any,
    fiscal_regime: str,
    controls: EndpointControls,
) -> EndpointEvaluation:
    """Evaluate one supplied policy under frozen terminal demographics.

    equal_rebate: root in (log price, transfer), residuals
      ((housing demand - supply)/supply, tax revenue - transfer * heads).
    fixed_transfer: root in log price only; tax revenue less transfers is
      reported as government surplus, without imposing an equal-rebate gate.

    The returned fiscal root residual is in model period-resource units, with
    no hidden rescaling. Inputs are neither reanchored nor recalibrated.
    """
    if fiscal_regime not in {"equal_rebate", "fixed_transfer"}:
        raise ValueError("An explicit equal_rebate or fixed_transfer regime is required")
    price = _positive("asset_price", asset_price)
    transfer = _finite("transfer", transfer)
    psi = _finite("psi_child", psi_child)
    if transfer < 0.0:
        raise ValueError("transfer must be nonnegative")
    if (isinstance(controls.maximum_inner_iterations, bool)
            or int(controls.maximum_inner_iterations) != controls.maximum_inner_iterations
            or controls.maximum_inner_iterations < 1):
        raise ValueError("maximum_inner_iterations must be a positive integer")
    if not 0.0 < _positive("inner_damping", controls.inner_damping) <= 1.0:
        raise ValueError("inner_damping must not exceed one")
    for name in ("distribution_tolerance", "birth_rate_tolerance",
                 "one_step_tolerance", "market_tolerance",
                 "fiscal_absolute_tolerance", "accounting_absolute_tolerance"):
        _positive(name, getattr(controls, name))
    grid = np.asarray(b_grid, dtype=float)
    initial = np.asarray(initial_g_pre, dtype=float)
    if (grid.ndim != 1 or grid.size < 2 or not np.all(np.isfinite(grid))
            or np.any(np.diff(grid) <= 0.0)):
        raise ValueError("b_grid must be finite and strictly increasing")
    if (initial.ndim != 7 or not np.all(np.isfinite(initial))
            or np.any(initial < 0.0) or float(initial.sum()) <= 0.0):
        raise ValueError("initial_g_pre must be a finite nonnegative seven-axis distribution")
    if int(parameters.I) != 1 or initial.shape[3] != int(parameters.J):
        raise ValueError("Endpoint requires one market and the declared age dimension")
    if grid.size != int(parameters.Nb) or initial.shape[0] != grid.size:
        raise ValueError("The supplied grid and distribution do not match parameters.Nb")
    policy_price = np.asarray(policy.price, dtype=float).reshape(-1)
    if policy_price.size != 1 or not np.isfinite(policy_price[0]):
        raise ValueError("Policy must contain one finite asset price")
    for name, actual, expected in (
        ("policy asset price", policy_price[0], price),
        ("parameter transfer", parameters.property_tax_lump_sum_transfer, transfer),
        ("parameter psi_child", parameters.psi_child, psi),
    ):
        if not math.isclose(_finite(name, actual), expected, rel_tol=0.0, abs_tol=1e-13):
            raise ValueError(f"Supplied {name} differs from its declared endpoint value")
    _positive("stationary user cost", parameters.user_cost_rate)
    if _finite("tau_H", parameters.tau_H) < 0.0:
        raise ValueError("tau_H must be nonnegative")
    if supply_rule.mode != "static-elastic":
        raise ValueError("Pass the inherited static-elastic supply rule explicitly")
    for name in ("initial_price", "initial_stock", "elasticity"):
        _positive(f"supply {name}", getattr(supply_rule, name))
    if int(demographic_primitives.start_year) != 2023:
        raise ValueError("The frozen person/head mapping must be anchored in 2023")
    if int(demographic_primitives.last_empirical_year) < 2023:
        raise ValueError("Terminal demographic year must not precede the 2023 anchor")

    inner = _solve_inner(
        policy=policy, parameters=parameters, b_grid=grid,
        initial_g_pre=initial, demographic_primitives=demographic_primitives,
        supply_rule=supply_rule,
        distribution_tolerance=controls.distribution_tolerance,
        birth_rate_tolerance=controls.birth_rate_tolerance,
        one_step_tolerance=controls.one_step_tolerance,
        maximum_iterations=controls.maximum_inner_iterations,
        damping=controls.inner_damping,
    )
    heads = float(np.sum(inner.g_pre))
    persons = float(np.sum(inner.persons.persons))
    demand, supply = float(inner.housing_demand), float(inner.housing_supply)
    market = (demand - supply) / supply if supply > 0.0 else math.nan
    fiscal = float(inner.government_budget_residual)
    outlays = transfer * heads
    revenue = fiscal + outlays
    residuals = {
        "housing_relative": market,
        "fiscal_absolute": fiscal,
        "transfer_per_head_gap": float(inner.equal_transfer_gap),
        "household_heads": heads,
        "resident_persons": persons,
        "housing_demand": demand,
        "housing_supply": supply,
        "tax_revenue": revenue,
        "transfer_outlays": outlays,
        "annual_births_per_head": float(inner.annual_births_per_head),
        "renewal_ratio": float(inner.renewal_ratio),
        "distribution_relative_l1": float(inner.distribution_mapping_relative_l1),
        "birth_rate_relative_gap": float(inner.annual_birth_rate_relative_gap),
        "person_one_step_relative_l1": float(inner.person_one_step_relative_l1),
        "head_one_step_relative_l1": float(inner.head_one_step_relative_l1),
        "age_head_gap": float(inner.age_head_one_step_max_abs),
        "household_person_head_gap": float(inner.household_person_head_gap),
        "returned_head_mass_gap": heads - float(np.sum(inner.persons.heads)),
        "inherited_supply_gap": supply - float(supply_rule.initial_stock) * (
            price / float(supply_rule.initial_price)) ** float(supply_rule.elasticity),
    }
    finite = all(math.isfinite(v) for v in residuals.values())
    nonnegative_arrays = all(
        np.all(np.isfinite(a)) and np.all(np.asarray(a) >= 0.0)
        for a in (inner.g_pre, inner.persons.persons, inner.persons.heads)
    )
    scale = max(1.0, heads)
    gates = {
        "finite_positive_state": finite and nonnegative_arrays and heads > 0 and persons > 0 and supply > 0,
        "nonnegative_flows": demand >= 0.0 and revenue >= -controls.accounting_absolute_tolerance
        and residuals["annual_births_per_head"] >= 0.0 and residuals["renewal_ratio"] >= 0.0,
        "inner_converged": bool(inner.converged),
        "distribution": abs(residuals["distribution_relative_l1"]) <= controls.distribution_tolerance,
        "birth_rate": abs(residuals["birth_rate_relative_gap"]) <= controls.birth_rate_tolerance,
        "person_one_step": abs(residuals["person_one_step_relative_l1"]) <= controls.one_step_tolerance,
        "head_one_step": abs(residuals["head_one_step_relative_l1"]) <= controls.one_step_tolerance,
        "age_head_identity": abs(residuals["age_head_gap"]) <= controls.one_step_tolerance * scale,
        "household_person_identity": abs(residuals["household_person_head_gap"]) <= controls.accounting_absolute_tolerance,
        "returned_head_identity": abs(residuals["returned_head_mass_gap"]) <= controls.accounting_absolute_tolerance,
        "inherited_supply": abs(residuals["inherited_supply_gap"]) <= controls.accounting_absolute_tolerance * max(1.0, abs(supply)),
        "housing_market": math.isfinite(market) and abs(market) <= controls.market_tolerance,
        "fiscal_regime": fiscal_regime == "fixed_transfer" or (
            math.isfinite(fiscal) and abs(fiscal) <= controls.fiscal_absolute_tolerance),
    }
    return EndpointEvaluation(
        fixed_point=inner, residuals=residuals,
        root_coordinates=("log_asset_price", "transfer") if fiscal_regime == "equal_rebate" else ("log_asset_price",),
        root_residuals=np.asarray([market, fiscal] if fiscal_regime == "equal_rebate" else [market]),
        gates=gates,
        contract={
            "schema": "e5f_matched_pf_endpoint_v1",
            "scope": "supplied-policy finite-level stationary endpoint evaluation",
            "population_growth_factor": 1.0,
            "fiscal_regime": fiscal_regime,
            "asset_price": price, "transfer": transfer, "psi_child": psi,
            "renter_price": price * float(parameters.user_cost_rate),
            "tau_H_period": float(parameters.tau_H),
            "supply": {name: getattr(supply_rule, name) for name in (
                "mode", "initial_price", "initial_stock", "elasticity")},
            "demographic_start_year": int(demographic_primitives.start_year),
            "demographic_terminal_year": int(demographic_primitives.last_empirical_year),
            "choice_flags": {name: bool(getattr(parameters, name, False)) for name in (
                "joint_nested_choice", "two_shock_choice", "fertility_nest_choice",
                "exhaustive_saving_control")},
            "controls": vars(controls).copy(),
        },
    )
