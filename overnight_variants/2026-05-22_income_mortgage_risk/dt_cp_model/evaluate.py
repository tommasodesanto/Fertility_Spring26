"""Direct single-parameter-vector model evaluation.

This module intentionally avoids SMM loss construction and geography inversion.
It is for the narrow question: given a structural parameter vector, does the
Python model solve and how long does it take?
"""

from __future__ import annotations

from types import SimpleNamespace
from typing import Sequence

import numpy as np

from .parameters import build_calibration_setup
from .solver import run_model_cp_dt
from .theta import apply_theta


def solve_theta(
    theta: Sequence[float],
    setup_mode: str = "fast",
    max_iter_eq: int | None = None,
    verbose: bool = False,
) -> tuple[SimpleNamespace, SimpleNamespace, np.ndarray]:
    setup = build_calibration_setup(setup_mode)
    P = apply_theta(setup.P_base, theta, setup.names)
    if max_iter_eq is not None:
        P.max_iter_eq = int(max_iter_eq)
    return run_model_cp_dt(P, verbose=verbose)
