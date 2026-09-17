"""Bounded scalar price root over fresh, explicitly supplied endpoint solves.

Only fixed-transfer endpoints are accepted. No model imports, fiscal closure,
preference selection, price reanchoring, or tolerance relaxation occurs here.
Cache entries are scalar records: large endpoint arrays are retained only for
the best point. A successful search must pass a fresh final replay.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import time
from typing import Any, Callable


@dataclass
class PriceRootResult:
    status: str
    evaluations: int
    cache_hits: int
    records: list[dict[str, Any]]
    best_record: dict[str, Any] | None
    best_endpoint: Any | None
    final_endpoint: Any | None = None
    replay_residual_gap: float | None = None

    @property
    def converged(self) -> bool:
        return self.status == "complete_reproduced_root"


class _Stop(Exception):
    pass


def solve_price_root(
    *,
    fresh_evaluate: Callable[[float], Any],
    start_price: float,
    bound_ratios: tuple[float, float],
    maximum_evaluations: int,
    deadline_monotonic: float,
    market_tolerance: float,
    replay_tolerance: float,
    progress_callback: Callable[[dict[str, Any], dict[str, Any] | None], None] | None = None,
) -> PriceRootResult:
    """Root (D-S)/S in log price, reserving one uncached final evaluation.

    deadline_monotonic is an absolute time.monotonic() deadline. All endpoint
    calls, including failures and replay, count toward maximum_evaluations.
    Negative excess demand searches below the start price; positive excess
    demand searches above it. Search visits the log midpoint toward that bound,
    then the bound itself if needed, before safeguarded root refinement.

    A narrow bracket is NEVER an acceptance criterion. Endpoint mapping_valid,
    endpoint.accepted, the explicit market gate, and residual reproduction must
    all hold for a completed root. A caller callback should persist the supplied
    latest and best scalar records; callback errors are propagated.
    """
    start = float(start_price)
    low_ratio, high_ratio = map(float, bound_ratios)
    if not all(math.isfinite(v) and v > 0 for v in (start, low_ratio, high_ratio)):
        raise ValueError("Price and bound ratios must be finite and positive")
    if not low_ratio < 1.0 < high_ratio:
        raise ValueError("Bound ratios must strictly surround one")
    if (isinstance(maximum_evaluations, bool)
            or int(maximum_evaluations) != maximum_evaluations or maximum_evaluations < 2):
        raise ValueError("At least two integer evaluations are needed, including final replay")
    for name, value in (("market_tolerance", market_tolerance), ("replay_tolerance", replay_tolerance)):
        if not math.isfinite(float(value)) or float(value) <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if market_tolerance > 2e-4:
        raise ValueError("The maintained 2e-4 market gate cannot be relaxed")
    if not math.isfinite(float(deadline_monotonic)):
        raise ValueError("An explicit finite monotonic deadline is required")
    lower, upper = math.log(start) + math.log(low_ratio), math.log(start) + math.log(high_ratio)
    try:
        lower_price, upper_price = math.exp(lower), math.exp(upper)
    except OverflowError as exc:
        raise ValueError("Price bounds overflow") from exc
    if not (0 < lower_price < start < upper_price < math.inf):
        raise ValueError("Price bounds must be representable and surround the start")
    result = PriceRootResult("running", 0, 0, [], None, None)
    cache: dict[float, dict[str, Any]] = {}

    def publish(record):
        result.records.append(record)
        if progress_callback is not None:
            progress_callback(dict(record), None if result.best_record is None else dict(result.best_record))

    def stop(status):
        result.status = status
        raise _Stop

    def evaluate(log_price, *, replay=False):
        if time.monotonic() >= deadline_monotonic:
            stop("incomplete_deadline")
        price = math.exp(log_price)
        if not replay and price in cache:
            result.cache_hits += 1
            # Distinct log floats can map to the same executable price.
            return dict(cache[price], log_price=log_price)
        limit = maximum_evaluations if replay else maximum_evaluations - 1
        if result.evaluations >= limit:
            stop("incomplete_budget")
        result.evaluations += 1
        began = time.monotonic()
        record = dict(evaluation=result.evaluations, price=price, log_price=log_price,
                      phase="final_replay" if replay else "search", status="running")
        try:
            endpoint = fresh_evaluate(price)
            if endpoint.contract.get("fiscal_regime") != "fixed_transfer":
                raise ValueError("Endpoint must explicitly use fixed_transfer")
            if float(endpoint.contract["asset_price"]) != price:
                raise ValueError("Endpoint was evaluated at a different asset price")
            if tuple(endpoint.root_coordinates) != ("log_asset_price",):
                raise ValueError("Endpoint must expose only the log-price root coordinate")
            residual = float(endpoint.residuals["housing_relative"])
            root_vector = endpoint.root_residuals
            if len(root_vector) != 1 or float(root_vector[0]) != residual:
                raise ValueError("Endpoint root residual differs from its housing residual")
            valid = bool(endpoint.mapping_valid) and math.isfinite(residual)
            accepted = bool(endpoint.accepted)
            record.update(housing_residual=residual, mapping_valid=valid,
                          endpoint_accepted=accepted, endpoint_gates=dict(endpoint.gates),
                          status="evaluated" if valid else "invalid_mapping")
        except Exception as exc:
            record.update(status="evaluation_failed", error=repr(exc),
                          elapsed_seconds=time.monotonic() - began)
            publish(record)
            stop("incomplete_evaluation_failure")
        record["elapsed_seconds"] = time.monotonic() - began
        if valid and (result.best_record is None or
                      abs(residual) < abs(result.best_record["housing_residual"])):
            result.best_record = dict(record)
            result.best_endpoint = endpoint
        publish(record)
        if not valid:
            stop("incomplete_invalid_mapping")
        if time.monotonic() >= deadline_monotonic:
            stop("incomplete_deadline")
        if replay:
            result.final_endpoint = endpoint
        else:
            cache[price] = record
        return record

    def meets_market(record):
        return abs(record["housing_residual"]) <= market_tolerance

    def certify(candidate):
        if not candidate["endpoint_accepted"]:
            stop("incomplete_endpoint_acceptance")
        reference_residual = candidate["housing_residual"]
        replay = evaluate(candidate["log_price"], replay=True)
        result.replay_residual_gap = abs(replay["housing_residual"] - reference_residual)
        if (not replay["endpoint_accepted"] or not meets_market(replay)
                or result.replay_residual_gap > replay_tolerance):
            stop("incomplete_replay_failure")
        result.status = "complete_reproduced_root"

    try:
        center = evaluate(math.log(start))
        if meets_market(center):
            certify(center)
            return result
        bound = lower if center["housing_residual"] < 0 else upper
        previous = center
        bracket = None
        for trial in ((center["log_price"] + bound) / 2, bound):
            current = evaluate(trial)
            if meets_market(current):
                certify(current)
                return result
            if previous["housing_residual"] * current["housing_residual"] < 0:
                bracket = sorted((previous, current), key=lambda row: row["log_price"])
                break
            previous = current
        if bracket is None:
            stop("incomplete_no_bracket")
        left, right = bracket
        while True:
            xl, xr = left["log_price"], right["log_price"]
            fl, fr = left["housing_residual"], right["housing_residual"]
            width = xr - xl
            trial = xl - fl * width / (fr - fl)
            # Keep every trial away from a bracket endpoint to guarantee shrinkage.
            if not math.isfinite(trial) or not xl + .1 * width <= trial <= xr - .1 * width:
                trial = (xl + xr) / 2
            if not xl < trial < xr:
                stop("incomplete_price_resolution")
            current = evaluate(trial)
            if meets_market(current):
                certify(current)
                return result
            if fl * current["housing_residual"] < 0:
                right = current
            else:
                left = current
    except _Stop:
        return result
