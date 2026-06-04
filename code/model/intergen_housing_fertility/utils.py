"""Numerical utilities for the intergenerational housing/fertility model."""

from __future__ import annotations

import json
from types import SimpleNamespace
from typing import Any

import numpy as np


def make_asset_grid(P: SimpleNamespace) -> np.ndarray:
    """Build a liquid-asset grid with extra resolution near zero."""

    n = int(P.Nb)
    if n < 8:
        raise ValueError("Nb must be at least 8.")
    n1 = max(2, round(0.20 * n))
    n2 = max(3, round(0.45 * n))
    n3 = n - n1 - n2
    if n3 < 3:
        n3 = 3
        n2 = n - n1 - n3
    low = np.linspace(P.b_min, 0.0, n1 + 1)[:-1]
    mid = np.linspace(0.0, P.b_mid, n2 + 1)[:-1]
    u = np.linspace(0.0, 1.0, n3 + 1)[1:]
    high = P.b_mid + (P.b_max - P.b_mid) * (u ** P.b_grid_power)
    grid = np.concatenate([low, mid, high]).astype(float)
    grid[np.argmin(np.abs(grid))] = 0.0
    return np.unique(grid)


def interpolate_age_profile(P: SimpleNamespace) -> np.ndarray:
    ages = P.age_start + np.arange(P.J)
    return np.interp(ages, P.income_age_breaks, P.income_age_values)


def housing_need(P: SimpleNamespace, children: int) -> float:
    return float(P.h_need_base + P.h_need_child * children)


def child_goods_cost(P: SimpleNamespace, children: int, income: float) -> float:
    return float(P.child_cost_0 * children + P.child_cost_income * income * children)


def jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {k: jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, SimpleNamespace):
        return {k: jsonable(v) for k, v in vars(value).items()}
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def dumps_json(value: Any) -> str:
    return json.dumps(jsonable(value), indent=2, sort_keys=True)
