#!/usr/bin/env python3
"""Pure micro-tests for the prepared fixed-credit overlay; no Bellman or GE run."""
from __future__ import annotations

import ast
import importlib.util
import inspect
import os
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
OVERLAY = HERE / "overlay"


def repo_root() -> Path:
    return next(parent for parent in HERE.parents if (parent / "code/model").is_dir())


def isolated_solver_credit_functions():
    tree = ast.parse((OVERLAY / "solver.py").read_text())
    wanted = {"debt_rule_at_age", "fixed_unsecured_credit_active", "renter_borrowing_floor"}
    nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in wanted]
    module = ast.Module(body=nodes, type_ignores=[])
    ast.fix_missing_locations(module)
    ns = {"np": np, "SimpleNamespace": SimpleNamespace, "Any": Any,
          "unsecured_debt_floor": lambda u, s, d: np.minimum(s * np.minimum(np.asarray(u, float), 0.0), -d)}
    exec(compile(module, "isolated_credit_functions", "exec"), ns)
    return ns


def load_kernels():
    spec = importlib.util.spec_from_file_location("overlay_kernels", OVERLAY / "kernels.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_credit_floor() -> None:
    fn = isolated_solver_credit_functions()["renter_borrowing_floor"]
    P = SimpleNamespace(J=17, debt_taper_weights=np.linspace(1.0, 0.0, 18), debt_caps=np.linspace(5.0, 0.0, 18), unsecured_credit_limit=None, use_age_survival=False)
    current = np.array([-4.0, -1.0, 2.0])
    for j in (0, 6, 11):
        expected = np.minimum(P.debt_taper_weights[j + 1] * np.minimum(current, 0.0), -P.debt_caps[j + 1])
        np.testing.assert_array_equal(fn(P, current, j), expected)  # absent/None legacy identity
    for D in (0.0, 2.5):
        P.unsecured_credit_limit = D
        for j in (0, 6, 11):
            np.testing.assert_array_equal(fn(P, current, j), np.full_like(current, -D))
    P.unsecured_credit_limit = 2.5
    assert np.all(fn(P, current, 16) == 0.0)  # terminal estate restriction
    P.use_age_survival = True
    P.survival_probs = np.ones(17)
    P.survival_probs[6] = 0.999
    assert np.all(fn(P, current, 6) == 0.0)  # positive current death probability


def test_parameter_contract() -> None:
    spec = importlib.util.spec_from_file_location("overlay_parameters", OVERLAY / "parameters.py")
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for bad in (-0.01, float("inf"), float("nan"), [1.0]):
        P = mod.setup_parameters()
        try:
            mod.apply_overrides(P, {"unsecured_credit_limit": bad})
        except ValueError:
            pass
        else:
            raise AssertionError(f"invalid credit limit accepted: {bad}")
    P = mod.setup_parameters()
    assert mod.apply_overrides(P, {"unsecured_credit_limit": 0.0}).unsecured_credit_limit == 0.0
    delattr(P, "unsecured_credit_limit")
    assert mod.apply_overrides(P, {"unsecured_credit_limit": 0.0}).unsecured_credit_limit == 0.0
    P.unsecured_credit_limit = [1.0]
    try:
        mod.build_debt_caps(P)
    except ValueError:
        pass
    else:
        raise AssertionError("direct debt-cap rebuild accepted a list")
    P = mod.setup_parameters()
    try:
        mod.apply_overrides(P, {"unsecured_credit_limit": 1.0, "native_solvency_credit": True})
    except ValueError:
        pass
    else:
        raise AssertionError("natural-credit conflict did not fail fast")


def test_raw_sale_gate_and_structure() -> None:
    k = load_kernels()
    b = np.array([-1.0, 0.0, 1.0])
    V = np.zeros((3, 2, 1, 1, 1))
    heq = np.array([[0.0, 0.5]])  # raw owner-sale balances: -0.5, 0.5, 1.5
    h = np.zeros((1, 2)); dp = np.zeros((1, 2, 1, 1)); bm = np.full((1, 2, 1, 1), -99.0)
    birth = np.zeros((1, 1, 2, 2), dtype=np.bool_); grant = np.zeros((1, 2, 1, 1))
    _, choice_off = k.tenure_choice_kernel(V, b, heq, h, dp, bm, birth, grant, V, False, False, False)
    _, choice_on = k.tenure_choice_kernel(V, b, heq, h, dp, bm, birth, grant, V, False, False, True)
    assert choice_off[0, 1, 0, 0, 0] == 0
    assert choice_on[0, 1, 0, 0, 0] != 0
    _, _, prob_on = k.tenure_logit_kernel(V, b, heq, h, dp, bm, birth, grant, 0.1, V, False, True)
    assert prob_on[0, 1, 0, 0, 0, 0] == 0.0
    # Raw balances at -1, 0, +1: the equality boundary remains feasible.
    heq_boundary = np.array([[0.0, 1.0]])
    _, choice_boundary = k.tenure_choice_kernel(V, b, heq_boundary, h, dp, bm, birth, grant, V, False, False, True)
    _, _, prob_boundary = k.tenure_logit_kernel(V, b, heq_boundary, h, dp, bm, birth, grant, 0.1, V, False, True)
    assert choice_boundary[0, 1, 0, 0, 0] == 0 and prob_boundary[0, 1, 0, 0, 0, 0] > 0.0
    if os.environ.get("NUMBA_DISABLE_JIT", "0") != "1":
        assert k.NUMBA_AVAILABLE and k.tenure_choice_kernel.signatures and k.tenure_logit_kernel.signatures
    original = ast.parse((repo_root() / "code/model/intergen_eqscale_seq_optimized/solver.py").read_text())
    prepared = ast.parse((OVERLAY / "solver.py").read_text())
    original_kernels = ast.parse((repo_root() / "code/model/intergen_eqscale_seq_optimized/kernels.py").read_text())
    prepared_kernels = ast.parse((OVERLAY / "kernels.py").read_text())
    def nodes(tree, name):
        return [ast.dump(n, include_attributes=False) for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == name]
    assert nodes(original_kernels, "full_owner_block_kernel") == nodes(prepared_kernels, "full_owner_block_kernel")
    def calls(tree):
        return [ast.dump(n, include_attributes=False) for n in ast.walk(tree) if isinstance(n, ast.Call) and getattr(n.func, "id", None) == "full_owner_block_kernel"]
    assert calls(original) == calls(prepared)


def test_native_renter_floor_replaces_legacy_line() -> None:
    k = load_kernels()
    grid = np.array([-2.0, 0.0, 2.0])
    args = (np.array([2.0, 2.0, 2.0]), np.array([2.0, 2.0, 2.0]), np.array([[100.0], [0.0], [-100.0]]), np.zeros((3, 1)), 0, grid,
            np.array([0.0]), np.array([0.0]), np.array([0.0]), np.array([0.0]), np.array([0.5]), np.array([1.0]),
            1.0, 10.0, 1e-6, 0.0, 0.0, 0.5, 0.5, 0.95, 0.0, 0.0, 0.381966, 0.618034, 1e-5)
    legacy = k.full_renter_block_kernel(*args)[1]
    positive = k.full_renter_block_kernel(*args, fixed_renter_floor=-1.0)[1]
    zero = k.full_renter_block_kernel(*args, fixed_renter_floor=0.0)[1]
    assert legacy[1, 0] >= 0.0 and positive[1, 0] < -0.5 and zero[1, 0] >= 0.0
    if os.environ.get("NUMBA_DISABLE_JIT", "0") != "1":
        assert k.NUMBA_AVAILABLE and k.full_renter_block_kernel.signatures


if __name__ == "__main__":
    test_credit_floor(); test_parameter_contract(); test_raw_sale_gate_and_structure(); test_native_renter_floor_replaces_legacy_line()
    print("fixed-credit pure contract tests passed")
