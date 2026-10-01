"""Bounded scalar price proposals; the caller retains equilibrium certification.

This module does not solve households, change tolerances, or cache solutions.
It only proposes prices for an explicit one-market warm start. A failed search
returns its best point for the caller's ordinary bracket/refinement routine.
"""
from __future__ import annotations

import math
from typing import Any, Callable


def search_warm_price(
    evaluate: Callable[[float], tuple[float, float, Any]],
    *,
    initial_price: float,
    initial_slope: float | None,
    lower_bound: float,
    upper_bound: float,
    tolerance: float,
) -> tuple[float, Any, float, dict[str, Any]]:
    """Try the supplied price and at most eight safeguarded secant proposals.

    ``evaluate`` returns signed excess demand, relative market error, and a
    solution payload. Only the best payload is retained. Global price bounds
    always apply; without a bracket each proposal stays within ten percent
    of its predecessor. No sign or monotonicity assumption certifies success:
    only the caller's unchanged relative-error tolerance does.
    """
    values = (initial_price, lower_bound, upper_bound, tolerance)
    if (not all(math.isfinite(v) for v in values)
            or not 0 < lower_bound <= initial_price <= upper_bound
            or lower_bound >= upper_bound or tolerance <= 0):
        raise ValueError("Warm price search requires finite positive bounds, price and tolerance")
    if initial_slope is not None and not math.isfinite(initial_slope):
        raise ValueError("Warm price slope must be finite or None")

    trace: list[dict[str, float]] = []
    seen: set[float] = set()
    best_price, best_payload, best_metric = initial_price, None, math.inf

    def visit(price: float) -> tuple[float, float]:
        nonlocal best_price, best_payload, best_metric
        excess, metric, payload = evaluate(price)
        if not math.isfinite(excess) or not math.isfinite(metric) or metric < 0:
            raise RuntimeError("Warm price search received a nonfinite or negative market metric")
        seen.add(price)
        trace.append(dict(price=price, excess=excess, metric=metric))
        if metric < best_metric:
            best_price, best_payload, best_metric = price, payload, metric
        return excess, metric

    def clipped_step(price: float, proposal: float) -> float:
        return min(max(proposal, lower_bound, price / 1.10), upper_bound, price * 1.10)

    price = float(initial_price)
    excess, metric = visit(price)
    if metric >= tolerance:
        proposal = (price - excess / initial_slope
                    if initial_slope not in (None, 0.0)
                    else price * (1.02 if excess > 0 else 1.0 / 1.02))
        if not math.isfinite(proposal):
            proposal = price * (1.02 if excess > 0 else 1.0 / 1.02)
        proposal = clipped_step(price, proposal)
        for _ in range(8):
            if proposal in seen:
                break
            next_excess, next_metric = visit(proposal)
            if next_metric < tolerance:
                break
            candidate = (proposal - next_excess * (proposal - price) / (next_excess - excess)
                         if next_excess != excess else 0.5 * (proposal + price))
            ordered = sorted(trace, key=lambda row: row["price"])
            brackets = [(a["price"], b["price"]) for a, b in zip(ordered, ordered[1:])
                        if (a["excess"] < 0 < b["excess"]
                            or b["excess"] < 0 < a["excess"])]
            if brackets:
                lo, hi = min(brackets, key=lambda pair: pair[1] - pair[0])
                if not math.isfinite(candidate) or not lo < candidate < hi or candidate in seen:
                    candidate = 0.5 * (lo + hi)
            else:
                if not math.isfinite(candidate):
                    candidate = proposal * (1.02 if next_excess > 0 else 1.0 / 1.02)
                candidate = clipped_step(proposal, candidate)
            price, excess, proposal = proposal, next_excess, float(candidate)

    return best_price, best_payload, best_metric, {
        "used": True,
        "initial_price": float(initial_price),
        "initial_slope": initial_slope,
        "price_evaluations": len(trace),
        "trace": trace,
        "fallback_required": best_metric >= tolerance,
    }
