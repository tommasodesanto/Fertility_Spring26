"""Experimental IID inheritance receipts as pre-choice liquid wealth.

Install only after the source-pinned purchase and net-estate adapters. No
frozen file is changed. The Bellman expectation and both forward operators
use the same positive interpolation operator. This supplies stationary and
cohort measurement support, not a funded general-equilibrium closure.
"""
from __future__ import annotations

import copy
import csv
import difflib
import functools
import hashlib
import inspect
from pathlib import Path

import numpy as np

from e5f_estate_receipt_jump import (
    backward_expectation, build_estate_receipt_jump_plan, forward_transport,
)

CASES = ("no_receipt", "conditional_mean", "receipt_lottery")
GROUP_WEIGHTS = {"bottom50": .5, "middle40": .4, "top10": .1}


def pooled_profile(path, ages):
    """Pool published usual-income groups; explicitly zero unsupported ages.

    The tail rule is an experimental support restriction, not an estimate.
    Each positive receipt is represented by its age-conditional mean amount.
    """
    rows = list(csv.DictReader(Path(path).open()))
    indexed = {}
    for row in rows:
        key = (int(row["age"]), row["income_group"])
        if key in indexed:
            raise ValueError("Duplicate published age/group")
        p = float(row["probability_4y"])
        mu = float(row["mean_amount_4y"])
        amount = float(row["conditional_amount_4y"])
        if not (np.isfinite([p, mu, amount]).all() and 0 <= p <= 1 and mu >= 0
                and amount >= 0 and np.isclose(p * amount, mu, rtol=1e-12, atol=1e-12)):
            raise ValueError("Invalid mapped receipt profile")
        indexed[key] = row
    if set(indexed) != {(a, g) for a in range(25, 81) for g in GROUP_WEIGHTS}:
        raise ValueError("Published age/group support differs")
    result = []
    for age in ages:
        if not float(age).is_integer():
            raise ValueError("This diagnostic uses exact integer age nodes")
        supported = 25 <= age <= 80
        probability = mean = 0.
        if supported:
            for group, weight in GROUP_WEIGHTS.items():
                row = indexed[(int(age), group)]
                p = float(row["probability_4y"])
                mu = float(row["mean_amount_4y"])
                probability += weight * p
                mean += weight * mu
        result.append(dict(age=float(age), published_age_supported=bool(supported),
                           probability=probability, relative_mean=mean,
                           relative_positive_amount=mean / probability if probability else 0.))
    return result


def configure(parameters, case, profile, scale, grid):
    if case not in CASES or not np.isfinite(scale) or scale < 0:
        raise ValueError("Unknown case or invalid funding scale")
    if (getattr(parameters, "estate_probe_case", None) != "net_valuation"
            or getattr(parameters, "estate_probe_transfer", None) != 0.):
        raise ValueError("A net-valuation reference without income transfers is required")
    if len(profile) != int(parameters.J):
        raise ValueError("Receipt age support differs from model")
    P = copy.deepcopy(parameters)
    plans = []
    for j, row in enumerate(profile):
        if row["age"] != float(P.age_start) + j * float(P.da):
            raise ValueError("Receipt and household ages differ")
        p = float(row["probability"])
        mu = scale * float(row["relative_mean"])
        amount = scale * float(row["relative_positive_amount"])
        if not np.isclose(p * amount, mu, rtol=1e-12, atol=1e-12):
            raise ValueError("The lottery and deterministic conditional means differ")
        support, weights = ([0.], [1.]) if case == "no_receipt" else (
            ([mu], [1.]) if case == "conditional_mean" else ([0., amount], [1. - p, p]))
        plans.append(build_estate_receipt_jump_plan(grid, support, weights))
    if not plans[0].is_identity:
        raise ValueError("This diagnostic does not change the entrant wealth distribution")
    P.estate_receipt_risk_case = case
    P.estate_receipt_risk_scale = float(scale)
    P.estate_receipt_risk_profile = copy.deepcopy(profile)
    P._estate_receipt_jump_plans = plans
    P._estate_receipt_forward_ledger = []
    return P


def _backward(P, receiving_age, values):
    return backward_expectation(P._estate_receipt_jump_plans[receiving_age], values)


def _forward(P, receiving_age, mass, *, record=False):
    moved, account = forward_transport(P._estate_receipt_jump_plans[receiving_age], mass)
    if abs(float(moved.sum()) - float(mass.sum())) > 1e-11 * max(1., float(mass.sum())):
        raise RuntimeError("Receipt operator loses population mass")
    residual = account["wealth_after"] - account["wealth_before"] - account["expected_receipt_flow"]
    if abs(residual) > 1e-11 * max(1., abs(account["wealth_after"])):
        raise RuntimeError("Receipt operator fails the financial-wealth identity")
    if record:
        P._estate_receipt_forward_ledger.append(dict(age_index=int(receiving_age), **account))
    return moved


def install(model, output_dir):
    """Generate minimal, fail-closed hooks in the reviewed frozen functions."""
    if not getattr(model, "_estate_probe_installed", False) or getattr(model, "_estate_receipt_risk_installed", False):
        raise RuntimeError("Install once, after the reviewed net-estate adapter")
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=False)
    model._estate_receipt_backward = _backward
    model._estate_receipt_forward = _forward
    net = inspect.getsource(model._estate_probe_bellman_net)
    anchor = """                        Vnr += transition_weight * next_values[
                            :, :, :, j + 1, znext, :, :
                        ]"""
    if net.count(anchor) != 1:
        raise RuntimeError("Net Bellman continuation anchor differs")
    revised_net = net.replace("def _estate_probe_bellman_net(", "def _estate_risk_bellman(", 1).replace(
        anchor, """                        Vnr += transition_weight * _estate_receipt_backward(
                            P, j + 1, next_values[:, :, :, j + 1, znext, :, :]
                        )""", 1)
    native_forward = inspect.getsource(model.forward_distribution_markov_income)
    anchor_forward = "    for j in range(J - 1):\n        _gate_dead_mass_at_age("
    anchor_terminal = "    _gate_dead_mass_at_age(\n        g[:, :, :, J - 1, :, :, :],"
    if native_forward.count(anchor_forward) != 1 or native_forward.count(anchor_terminal) != 1:
        raise RuntimeError("Native forward receipt anchors differ")
    revised_forward = native_forward.replace(
        "def forward_distribution_markov_income(", "def _estate_risk_forward_distribution(", 1).replace(
        anchor_forward, "    P._estate_receipt_forward_ledger = []\n    for j in range(J - 1):\n"
        "        g[:, :, :, j, :, :, :] = _estate_receipt_forward(P, j, g[:, :, :, j, :, :, :], record=True)\n"
        "        _gate_dead_mass_at_age(", 1).replace(
        anchor_terminal, "    g[:, :, :, J - 1, :, :, :] = _estate_receipt_forward(\n"
        "        P, J - 1, g[:, :, :, J - 1, :, :, :], record=True)\n" + anchor_terminal, 1)
    original_cohort = model.advance_cohort_one_period_markov_income
    cohort_source = inspect.getsource(original_cohort)
    pins = {}
    for name, before, after in (("bellman", net, revised_net), ("forward", native_forward, revised_forward)):
        generated = output / (name + ".generated.py")
        generated.write_text(after)
        (output / (name + ".diff")).write_text("".join(difflib.unified_diff(
            before.splitlines(True), after.splitlines(True), fromfile="frozen_" + name,
            tofile="experimental_receipt_" + name)))
        exec(compile(after, str(generated), "exec"), model.__dict__)
        pins[name] = {"original": hashlib.sha256(before.encode()).hexdigest(),
                      "generated": hashlib.sha256(after.encode()).hexdigest()}

    @functools.wraps(original_cohort)
    def cohort(*args, **kwargs):
        P = kwargs.get("P", args[6] if len(args) > 6 else None)
        age = kwargs.get("j", args[1] if len(args) > 1 else None)
        result = original_cohort(*args, **kwargs)
        return _forward(P, int(age) + 1, result)

    model.solve_bellman_full_markov_income = model._estate_risk_bellman
    model.forward_distribution_markov_income = model._estate_risk_forward_distribution
    model.advance_cohort_one_period_markov_income = cohort
    model._estate_receipt_risk_installed = True
    return dict(function_sha256=pins, cohort_source_sha256=hashlib.sha256(cohort_source.encode()).hexdigest(),
                timing="start of receiving period, before interest and choices; paid once as liquid wealth",
                income_observer_unchanged=True, additional_persistent_states=0,
                clipping_tolerance=0., risk="IID across receiving periods; no family linkage or lifetime receipt cap")


def normalized_ledger(P, solution):
    """Normalize raw native-forward receipt flows exactly as its population."""
    rows = P._estate_receipt_forward_ledger
    if len(rows) != int(P.J) or [r["age_index"] for r in rows] != list(range(int(P.J))):
        raise RuntimeError("Missing or repeated native receipt transitions")
    raw_mass = sum(row["transported_mass"] for row in rows)
    normalized_mass = float(np.asarray(solution.g).sum())
    scale = normalized_mass / raw_mass
    age_mass = np.asarray(solution.g).sum(axis=(0, 1, 2, 4, 5, 6))
    for j, row in enumerate(rows):
        if abs(row["transported_mass"] * scale - age_mass[j]) > 1e-10:
            raise RuntimeError("Receipt ledger and current age masses differ")
    return dict(normalization_scale=scale, by_age=[{
        key: value if key == "age_index" else value * scale for key, value in row.items()
    } for row in rows], paid_period=sum(row["expected_receipt_flow"] for row in rows) * scale,
        clipped_wealth_loss=sum(row["clipped_wealth_loss"] for row in rows) * scale)
