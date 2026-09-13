"""Minimal sequence-space bridge for the original four-vintage queue experiment.

This is deliberately an *adapter*, not a replacement household block.  Its
``evaluator`` must call the frozen native mapping, including the exact queue
law, household feasibility checks, population levels, and fiscal/market gates.
The bridge only packs dated unknowns and the already-scaled residual vector in
the ordering used by ``sequence_jacobian``:

    unknowns = (log_house_price, pension, rebate), shape (3, T)
    targets  = (housing_relative_imbalance, PAYGO_relative_imbalance,
                rebate_relative_imbalance), shape (3, T)

The native mapping owns the economics.  In particular, ``state`` is carried as
an opaque object and is never normalized to unit mass here.  It must contain
the level-valued household distribution and both four-entry queues.  Prices
are exponentiated only at the native-call boundary; pensions and rebates remain
in levels because that is the current root's coordinate convention.

If the official ``sequence_jacobian`` package is available, ``factor_ssj``
constructs its ``JacobianDict`` / ``FactoredJacobianDict`` and
``nonlinear_update`` applies its signed Newton update.  The package is not a
fast-news implementation by itself: obtaining that requires household-policy,
distribution/queue, and aggregate derivative objects listed in the companion
note.  A dense finite-difference Jacobian may be supplied for a *small*
directional smoke only, but it is not represented as a performance result.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping, Sequence

import numpy as np


UNKNOWN_NAMES = ("log_house_price", "pension", "rebate")
TARGET_NAMES = (
    "housing_relative_imbalance",
    "paygo_relative_imbalance_scaled_200",
    "rebate_relative_imbalance_scaled_200",
)


@dataclass(frozen=True)
class NativePathEvaluation:
    """Exact native result needed by the bridge.

    ``residual`` is the native solver's scaled vector, ordered as ``(3, T)``;
    the two fiscal rows already include the multiplier 200.  ``state`` is an
    opaque native terminal state with ``g_pre``, ``scheduled_entries``, and
    ``scheduled_raw_entries``.  ``receipt`` must report unchanged native gates.
    """

    residual: np.ndarray
    state: Any
    receipt: Mapping[str, Any]
    rows: Sequence[Mapping[str, Any]] = ()


NativeEvaluator = Callable[[np.ndarray, np.ndarray, np.ndarray], NativePathEvaluation]


def pack_unknowns(log_prices: np.ndarray, pensions: np.ndarray, rebates: np.ndarray) -> np.ndarray:
    """Validate and stack the dated root coordinates as a ``(3, T)`` array."""
    arrays = tuple(np.asarray(x, dtype=float) for x in (log_prices, pensions, rebates))
    if any(x.ndim != 1 for x in arrays) or len({x.size for x in arrays}) != 1 or arrays[0].size < 1:
        raise ValueError("log prices, pensions, and rebates must be equally sized one-dimensional paths")
    if not all(np.isfinite(x).all() for x in arrays):
        raise ValueError("dated root coordinates must be finite")
    if np.any(arrays[1] < 0.0) or np.any(arrays[2] < 0.0):
        raise ValueError("pensions and rebates must be nonnegative levels")
    return np.stack(arrays)


def unpack_unknowns(unknowns: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    unknowns = np.asarray(unknowns, dtype=float)
    if unknowns.ndim != 2 or unknowns.shape[0] != len(UNKNOWN_NAMES) or unknowns.shape[1] < 1:
        raise ValueError("unknowns must have shape (3, T)")
    if not np.isfinite(unknowns).all() or np.any(unknowns[1:] < 0.0):
        raise ValueError("invalid dated root coordinates")
    return tuple(unknowns[i].copy() for i in range(3))  # type: ignore[return-value]


class OriginalQueueSequenceSpaceBridge:
    """Call the exact native mapping without changing its state normalization."""

    def __init__(self, evaluator: NativeEvaluator):
        self.evaluator = evaluator

    def evaluate(self, unknowns: np.ndarray) -> NativePathEvaluation:
        log_prices, pensions, rebates = unpack_unknowns(unknowns)
        value = self.evaluator(np.exp(log_prices), pensions, rebates)
        residual = np.asarray(value.residual, dtype=float)
        if residual.shape != unknowns.shape or not np.isfinite(residual).all():
            raise ValueError("native residual must be finite and have shape (3, T)")
        for field in ("g_pre", "scheduled_entries", "scheduled_raw_entries"):
            if not hasattr(value.state, field):
                raise ValueError("native state loses original-queue object: " + field)
        return NativePathEvaluation(residual=residual, state=value.state,
                                    receipt=value.receipt, rows=value.rows)


def dense_directional_jacobian(bridge: OriginalQueueSequenceSpaceBridge, unknowns: np.ndarray,
                               direction: np.ndarray, step: float) -> np.ndarray:
    """Central directional derivative of the exact native residual.

    This deliberately costs two full native mappings and is only admissible for
    the authorized horizon-two smoke.  It does not claim fake-news efficiency.
    """
    unknowns = np.asarray(unknowns, dtype=float)
    direction = np.asarray(direction, dtype=float)
    if direction.shape != unknowns.shape or not np.isfinite(direction).all() or step <= 0.0:
        raise ValueError("direction must match unknowns and step must be positive")
    return (bridge.evaluate(unknowns + step * direction).residual
            - bridge.evaluate(unknowns - step * direction).residual) / (2.0 * step)


def factor_ssj(jacobian: np.ndarray, horizon: int):
    """Return official toolkit's factored Jacobian for the dated native system.

    ``jacobian`` has shape ``(3T, 3T)`` with target-major / unknown-major
    packing.  Importing is intentionally local and optional: no package is
    installed or vendored by this prototype.
    """
    matrix = np.asarray(jacobian, dtype=float)
    expected = len(TARGET_NAMES) * horizon
    if matrix.shape != (expected, expected) or not np.isfinite(matrix).all():
        raise ValueError("Jacobian must have shape (3T, 3T)")
    try:
        from sequence_jacobian.classes import FactoredJacobianDict, JacobianDict
    except ImportError as exc:
        raise RuntimeError("Install sequence-jacobian in an isolated environment before SSJ factoring") from exc
    nested = {}
    for oi, target in enumerate(TARGET_NAMES):
        nested[target] = {}
        for ui, unknown in enumerate(UNKNOWN_NAMES):
            nested[target][unknown] = matrix[oi*horizon:(oi+1)*horizon, ui*horizon:(ui+1)*horizon]
    return FactoredJacobianDict(JacobianDict(nested, TARGET_NAMES, UNKNOWN_NAMES, T=horizon), horizon)


def nonlinear_update(factored_jacobian: Any, residual: np.ndarray) -> np.ndarray:
    """Apply the official signed Newton step ``-H_U^{-1} H`` to a residual."""
    residual = np.asarray(residual, dtype=float)
    if residual.ndim != 2 or residual.shape[0] != len(TARGET_NAMES):
        raise ValueError("residual must have shape (3, T)")
    try:
        from sequence_jacobian.classes import ImpulseDict
    except ImportError as exc:
        raise RuntimeError("sequence-jacobian is required for this update") from exc
    result = factored_jacobian.apply(ImpulseDict({name: residual[i]
                                                  for i, name in enumerate(TARGET_NAMES)}))
    return np.stack([np.asarray(result[name], dtype=float) for name in UNKNOWN_NAMES])


def native_smoke_command() -> str:
    """Exact bounded job design; deliberately not submitted by this prototype."""
    return (
        "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 "
        "/share/apps/anaconda3/2025.06/bin/python "
        "code/model/tools/e5f_sequence_space_prototype.py --native-smoke "
        "--spec /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/"
        "batches/afternoon_original_queue_20260913a/spec.json --horizon 2 "
        "--max-evaluations 6 --wall-minutes 10 --cpus 1 --mem-gib 24 "
        "--account torch_pr_570_general"
    )


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Print the bounded native SSJ smoke contract.")
    parser.add_argument("--native-smoke", action="store_true")
    parser.add_argument("--spec")
    parser.add_argument("--horizon", type=int)
    parser.add_argument("--max-evaluations", type=int)
    parser.add_argument("--wall-minutes", type=int)
    parser.add_argument("--cpus", type=int)
    parser.add_argument("--mem-gib", type=int)
    parser.add_argument("--account")
    args = parser.parse_args()
    if args.native_smoke:
        if (args.horizon, args.max_evaluations, args.wall_minutes, args.cpus, args.mem_gib) != (2, 6, 10, 1, 24):
            raise ValueError("authorized native smoke is exactly horizon 2, <=6 evaluations, 10m, 1 CPU, 24 GiB")
        print(native_smoke_command())
    else:
        parser.print_help()
