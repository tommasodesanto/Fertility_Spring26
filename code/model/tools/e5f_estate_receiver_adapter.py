"""Experimental stationary estate accounting on the frozen B15 runtime.

Install after the reviewed current-income purchase adapter. Frozen files and
the empirical target observers are unchanged. These are equal, anticipated
period transfers to the age cells starting in [45, 65], not genealogical or
stochastic inheritance. No transition law or population closure is supplied.
All imports, source generation, tests and numerical work belong on Torch.
"""
from __future__ import annotations

import copy
import difflib
import functools
import hashlib
import inspect
import json
import math
from pathlib import Path
import time

import numpy as np

CASES = ("control", "net_valuation", "net_valuation_transfer")
TRANSFER_CASE = CASES[2]
FLOW_TOLERANCE = 1e-10


def _write(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def recipient_ages(P):
    ages = float(P.age_start) + np.arange(int(P.J)) * float(P.da)
    return ages, (ages >= 45.) & (ages <= 65.)


def configure(parameters, case, transfer=0.0):
    """Copy actual serialized primitives, adding only diagnostic settings."""
    if case not in CASES or not math.isfinite(transfer) or transfer < 0:
        raise ValueError("Unknown estate case or invalid transfer")
    if case != TRANSFER_CASE and transfer != 0:
        raise ValueError("Only the receiver case can pay estates")
    if (str(getattr(parameters, "estate_receiver", "none")) != "none"
            or bool(getattr(parameters, "bequest_net_of_selling_cost", False))
            or float(getattr(parameters, "estate_lump_sum_transfer", 0.)) != 0.
            or float(getattr(parameters, "estate_tax_rate", 0.)) != 0.):
        raise ValueError("Reference already changes estate accounting")
    if (int(parameters.I) != 1
            or not bool(getattr(parameters, "use_postdecision_current_distribution", True))
            or np.any(np.asarray(getattr(parameters, "child_earnings_penalty", 0.)) != 0.)
            or not 0. <= float(parameters.psi) < 1.):
        raise ValueError("Unsupported reference or death-distribution convention")
    out = copy.deepcopy(parameters)
    out.estate_probe_case = case
    out.estate_probe_transfer = float(transfer)
    return out


def transfer_at_age(P, age_index):
    if getattr(P, "estate_probe_case", "control") != TRANSFER_CASE:
        return 0.
    age = float(P.age_start) + int(age_index) * float(P.da)
    return float(P.estate_probe_transfer) if 45. <= age <= 65. else 0.


def install(model, output_dir):
    """Hook budgets and one reviewed Bellman estate-table construction.

    Current income already enters purchase eligibility as Y/R. Estate receipts
    follow that existing timing, enter resources once, and are not labor income.
    The control calls the original Bellman and income functions directly.
    """
    if getattr(model, "_estate_probe_installed", False):
        raise RuntimeError("Estate adapter must be installed exactly once")
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    original_bellman = model.solve_bellman_full_markov_income
    original_income = model.income_at_state
    original_annual = model.annual_gross_income_at_state
    source = inspect.getsource(original_bellman)
    # Require the reviewed purchase-income adapter, not a bare package solve.
    if "income_for_purchase" not in source or "bmo_purchase" not in source:
        raise RuntimeError("Reviewed current-income purchase runtime is missing")
    anchor = "            hv = p_hat[i] * P.H_own[ten - 1] if ten > 0 else 0.0\n"
    if source.count(anchor) != 1:
        raise RuntimeError("Frozen Bellman estate-table anchor changed")
    name = "def solve_bellman_full_markov_income("
    if source.count(name) != 1:
        raise RuntimeError("Unexpected Bellman function definition")
    revised = source.replace(name, "def _estate_probe_bellman_net(", 1).replace(
        anchor, anchor + "            hv *= 1.0 - float(P.psi)\n", 1)
    generated = output / "estate_bellman.generated.py"
    generated.write_text(revised)
    diff = "".join(difflib.unified_diff(source.splitlines(True), revised.splitlines(True),
                                     fromfile="reviewed_purchase_bellman",
                                     tofile="experimental_net_estate_bellman"))
    (output / "estate_bellman.diff").write_text(diff)
    exec(compile(revised, str(generated), "exec"), model.__dict__)
    net_bellman = model._estate_probe_bellman_net

    @functools.wraps(original_bellman)
    def bellman(*args, **kwargs):
        P = kwargs.get("P", args[2] if len(args) > 2 else None)
        if P is None:
            raise TypeError("Bellman parameter object is missing")
        case = getattr(P, "estate_probe_case", "control")
        if case not in CASES:
            raise ValueError("Unknown estate case")
        return (original_bellman if case == "control" else net_bellman)(*args, **kwargs)

    def income(P, i, j, z_value):
        base = original_income(P, i, j, z_value)
        if getattr(P, "estate_probe_case", "control") != TRANSFER_CASE:
            return base
        return base + transfer_at_age(P, j)

    def annual_income(P, i, j, z_value):
        if getattr(P, "estate_probe_case", "control") != TRANSFER_CASE:
            return original_annual(P, i, j, z_value)
        # Freeze the original measurement definition: exclude estate receipts
        # rather than silently grossing them up as taxable labor earnings.
        period_years = float(getattr(P, "period_years", getattr(P, "da", 1.0)))
        annual = original_income(P, i, j, z_value) / max(period_years, 1e-12)
        return annual / max(1. - float(getattr(P, "tau_pay", 0.)), 1e-12) if j < int(P.J_R) else annual

    model.solve_bellman_full_markov_income = bellman
    model.income_at_state = income
    model.annual_gross_income_at_state = annual_income
    model._estate_probe_installed = True
    receipt = dict(
        adapter="experimental_stationary_estate_receiver_v1",
        original_bellman_sha256=hashlib.sha256(source.encode()).hexdigest(),
        net_bellman_sha256=hashlib.sha256(revised.encode()).hexdigest(),
        control_calls_original_bellman=True,
        empirical_gross_estate_observer_unchanged=True,
        estate_receipts_excluded_from_annual_income_observer=True,
        transfer_timing="period resources and existing Y/R purchase eligibility; credited once",
        recipient_rule="equal per-period transfer to age cells starting in [45,65]",
        dynamics_implemented=False)
    _write(output / "estate_adapter.json", receipt)
    return receipt


def estate_accounts(sol, P, prices):
    """Post-saving estates from the same current measure as the native observer.

    D = sum_j (1-s_j) sum_x g_current(x,j) max(b'(x,j)+(1-psi)q h,0).
    At terminal age death is certain. Both D and transfer*M are PERIOD flows.
    Native reported bequests stay gross for comparable empirical target rows.
    """
    g, bp = np.asarray(sol.g), np.asarray(sol.bp_pol)
    if (g.ndim != 7 or bp.shape != g.shape or g.shape[2] != 1
            or g.shape[3] != int(P.J) or g.shape[1] != 1 + len(P.H_own)
            or not np.isfinite(g).all() or np.any(g < 0.) or not np.isfinite(bp).all()):
        raise ValueError("Invalid income-resolved estate distribution/policy")
    if not bool(getattr(P, "use_postdecision_current_distribution", True)):
        raise ValueError("Estates require the post-transaction current measure")
    price = np.asarray(prices, dtype=float).reshape(-1)
    if price.shape != (1,) or not np.isfinite(price).all() or price[0] <= 0.:
        raise ValueError("A positive one-market estate price is required")
    ages, recipients = recipient_ages(P)
    if bool(getattr(P, "use_age_survival", False)):
        survival = np.asarray(P.survival_probs, dtype=float)
        if (survival.shape != (int(P.J) - 1,) or not np.isfinite(survival).all()
                or np.any(survival < 0.) or np.any(survival > 1.)):
            raise ValueError("Invalid death schedule")
    else:
        survival = np.ones(int(P.J) - 1)
    death = np.r_[1. - survival, 1.]
    gross = net = 0.
    for j, probability in enumerate(death):
        if probability == 0.:
            continue
        for ten in range(g.shape[1]):
            hv = price[0] * float(P.H_own[ten - 1]) if ten else 0.
            cell = g[:, ten, 0, j]
            choice = bp[:, ten, 0, j]
            gross += float(probability) * float(np.sum(cell * np.maximum(choice + hv, 0.)))
            net += float(probability) * float(np.sum(cell * np.maximum(choice + (1. - float(P.psi)) * hv, 0.)))
    years = float(getattr(P, "period_years", P.da))
    native_gross = float(sol.annual_bequest_flow) * years
    if not math.isfinite(native_gross) or abs(gross - native_gross) > 1e-10 * max(1., abs(gross)):
        raise RuntimeError("Estate death measure disagrees with native gross-flow observer")
    mass = float(g[:, :, :, recipients].sum())
    transfer = float(getattr(P, "estate_probe_transfer", 0.))
    if not math.isfinite(transfer) or transfer < 0.:
        raise ValueError("Invalid estate payment")
    paid = transfer * mass
    return dict(generated_net_period=net, generated_gross_period=gross,
                recipient_mass=mass, paid_period=paid, residual=paid-net,
                transfer=transfer, recipient_ages=ages[recipients].tolist(),
                period_years=years, native_gross_flow_gap=gross-native_gross,
                net_external_outflow_period=net-paid,
                funded_transfer_required=getattr(P, "estate_probe_case", "control") == TRANSFER_CASE)


def solve_case(*, model, native_solver, parameters, b_grid, initial_prices,
               case, deadline_epoch, output_dir, max_solves=40):
    """Bounded stationary transfer fixed point; preferences and entry stay fixed."""
    if not getattr(model, "_estate_probe_installed", False):
        raise RuntimeError("Install and verify the estate adapter first")
    if not 1 <= int(max_solves) <= 40:
        raise ValueError("The estate loop permits at most forty native solves")
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    if (output / "estate_iterations.json").exists():
        raise FileExistsError("Do not overwrite an earlier estate loop")
    records = []
    transfer = 0.
    initial = np.asarray(initial_prices, dtype=float)
    limit = int(max_solves) if case == TRANSFER_CASE else 1
    for index in range(limit):
        if time.time() >= deadline_epoch:
            raise TimeoutError("Estate stage exhausted its wall-clock budget")
        P = configure(parameters, case, transfer)
        record = dict(index=index+1, status="started", transfer=transfer,
                      started_epoch=time.time())
        records.append(record)
        _write(output / "estate_iterations.json", records)
        tick = time.monotonic()
        try:
            sol, P, price, fiscal = native_solver(
                model=model, parameters=P, b_grid=b_grid, initial_prices=initial,
                payroll_tax=float(parameters.tau_pay), marginal_tolerance=1e-9,
                fiscal_tolerance=1e-6)
            if float(P.psi_child) != float(parameters.psi_child) or float(P.tau_pay) != float(parameters.tau_pay):
                raise RuntimeError("Estate diagnostic changed preferences or fiscal setting")
            accounts = estate_accounts(sol, P, price)
            if case == TRANSFER_CASE and accounts["recipient_mass"] <= 0.:
                raise RuntimeError("Positive recipient exposure is required")
            desired = accounts["generated_net_period"] / accounts["recipient_mass"] if case == TRANSFER_CASE else 0.
            relative = abs(desired-transfer) / max(abs(transfer), 1e-12)
            converged = case != TRANSFER_CASE or (
                abs(accounts["residual"]) <= FLOW_TOLERANCE and relative <= 1e-6)
            record.update(status="completed", seconds=time.monotonic()-tick,
                          price=np.asarray(price).tolist(), accounts=accounts,
                          next_undamped_transfer=desired, converged=converged)
            _write(output / "estate_iterations.json", records)
            _write(output / "latest_completed_iteration.json", record)
            if time.time() >= deadline_epoch:
                raise TimeoutError("Estate stage exhausted its wall-clock budget after solve")
            if converged:
                receipt = dict(**accounts, case=case, native_solves=len(records),
                               converged=True, flow_tolerance=FLOW_TOLERANCE,
                               preferences_fixed=True, population_closure="normalized-entry diagnostic; replacement must be reported")
                _write(output / "estate_receipt.json", receipt)
                return sol, P, price, fiscal, receipt
            transfer = 0.5 * (transfer + desired)
            # Reuse only prices; no unverified value-function warm start.
            initial = np.asarray(price, dtype=float)
        except Exception as exc:
            record.update(status="failed", seconds=time.monotonic()-tick,
                          error_type=type(exc).__name__, error=str(exc))
            _write(output / "estate_iterations.json", records)
            raise
    raise RuntimeError("Estate transfer did not balance within the contracted solve count")
