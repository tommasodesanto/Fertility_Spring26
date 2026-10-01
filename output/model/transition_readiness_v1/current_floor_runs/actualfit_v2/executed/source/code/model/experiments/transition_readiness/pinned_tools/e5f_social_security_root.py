"""Bounded joint housing/Social Security roots with caller-owned economics.

This module neither chooses a fiscal closure nor constructs household policies.
The evaluator must solve households with the entire proposed price and pension
or payroll-tax path anticipated, then return residuals from that population.
Property-tax rebates are a separate budget and are not a residual here.
"""
from __future__ import annotations

import copy
from collections.abc import Mapping

import numpy as np

from e5f_matched_pf_path_root import solve_price_path


class CandidateDomainError(ValueError):
    """An explicitly inadmissible economic trial, not a programming failure.

    An evaluator may raise this to return an invalid mapping to the bounded
    root. Other exceptions propagate; no implicit retry hides a broken solve.
    """


def solve_social_security_path(
    *, closure, initial_prices, initial_fiscal_values, evaluate, project_prices,
    fiscal_bounds, market_tolerance, fiscal_tolerance, market_slope, fiscal_slope,
    max_log_step, damping, max_evaluations, deadline_monotonic,
    max_condition_number, worsening_factor, final_reproduction_tolerance,
    callback=None, initial_jacobian=None,
):
    """Solve equal-length housing and Social Security residual paths.

    ``closure`` is explicitly ``fixed_tax`` (positive period pension is the
    fiscal unknown) or ``fixed_pension`` (payroll tax is the fiscal unknown).
    The caller owns the fixed object. ``fiscal_bounds=(lower, upper)`` are
    finite, positive, scalar bounds; tax bounds must lie strictly below one.
    Initial fiscal values must already satisfy those bounds.

    ``evaluate(prices, fiscal_values)`` returns a mapping containing vectors
    ``market_residual`` and ``fiscal_residual``, Boolean ``mapping_valid``, and
    optional lightweight ``payload``. Housing residuals are (demand-supply)/
    supply. Fiscal residuals are revenue minus pension outlays, with any scale
    explicitly chosen by the caller and identical across all evaluations.
    ``fiscal_tolerance`` applies to that supplied fiscal residual directly.
    The evaluator must return an invalid mapping or raise CandidateDomainError
    for an inadmissible trial; other exceptions propagate.

    Positive slopes describe approximate magnitudes of the own derivatives
    with respect to log price and log fiscal value. The tax residual is sign
    reversed internally. Fiscal residuals are scaled by market_tolerance /
    fiscal_tolerance. A block-diagonal default Jacobian accounts for this
    scaling after every reset. Both blocks retain physical log coordinates,
    and ``max_log_step`` bounds both price and fiscal proportional updates.

    Optional ``initial_jacobian`` and returned ``final_jacobian`` use the
    original, unsigned residuals [housing, revenue-minus-outlays] and physical
    log coordinates [log prices, log fiscal values]. Every evaluation, including
    rejected trials and the reserved fresh replay, counts toward the budget.
    Final replay must reproduce each original residual block within the same
    explicit absolute ``final_reproduction_tolerance``. Callers provide a
    watchdog for the individual evaluations and persist callback records.
    """
    if closure not in {"fixed_tax", "fixed_pension"}:
        raise ValueError("Explicit fixed_tax or fixed_pension closure required")
    prices = np.asarray(initial_prices, dtype=float)
    fiscal = np.asarray(initial_fiscal_values, dtype=float)
    if (prices.ndim != 1 or not prices.size or fiscal.shape != prices.shape
            or not np.isfinite(prices).all() or np.any(prices <= 0)
            or not np.isfinite(fiscal).all() or np.any(fiscal <= 0)):
        raise ValueError("Price and fiscal paths must be equal positive finite vectors")
    prices, fiscal = prices.copy(), fiscal.copy()
    try:
        lower, upper = map(float, fiscal_bounds)
    except (TypeError, ValueError) as exc:
        raise ValueError("Explicit scalar lower and upper fiscal bounds required") from exc
    if (not np.isfinite([lower, upper]).all() or not 0 < lower < upper
            or (closure == "fixed_pension" and upper >= 1)):
        raise ValueError("Fiscal bounds must be positive; payroll-tax upper bound must be below one")
    if np.any(fiscal < lower) or np.any(fiscal > upper):
        raise ValueError("Initial fiscal values are outside the explicit fiscal bounds")
    scales = np.asarray([market_tolerance, fiscal_tolerance, market_slope, fiscal_slope], dtype=float)
    if not np.isfinite(scales).all() or np.any(scales <= 0):
        raise ValueError("Separate market/fiscal tolerances and slopes must be finite and positive")
    if not callable(evaluate) or not callable(project_prices) or (callback is not None and not callable(callback)):
        raise ValueError("Evaluation, price projection and optional callback must be callable")
    if not np.isfinite(final_reproduction_tolerance) or final_reproduction_tolerance < 0:
        raise ValueError("Replay tolerance must be finite and nonnegative")
    n = prices.size
    residual_scale = float(market_tolerance / fiscal_tolerance)
    fiscal_default_slope = float(residual_scale * fiscal_slope)
    if not np.isfinite([residual_scale, fiscal_default_slope]).all() or min(residual_scale, fiscal_default_slope) <= 0:
        raise ValueError("Fiscal residual scale and default slope are not representable")
    sign = 1.0 if closure == "fixed_tax" else -1.0
    row_scale = np.r_[np.ones(n), np.full(n, sign * residual_scale)]
    default_jacobian = -np.diag(np.r_[np.full(n, market_slope), np.full(n, fiscal_default_slope)])
    evaluation_count = 0
    domain_rejections = {}

    def decode(values):
        values = np.asarray(values, dtype=float)
        if values.shape != (2*n,) or not np.isfinite(values).all() or np.any(values <= 0):
            raise ValueError("Joint root coordinates must be positive and finite")
        if np.any(values[n:] < lower) or np.any(values[n:] > upper):
            raise ValueError("Fiscal candidate lies outside the explicit fiscal bounds")
        return values[:n].copy(), values[n:].copy()

    def project(values):
        values = np.asarray(values, dtype=float)
        if values.shape != (2*n,) or not np.isfinite(values).all() or np.any(values <= 0):
            raise ValueError("Projection received invalid joint root coordinates")
        projected_prices = np.asarray(project_prices(values[:n].copy()), dtype=float)
        if (projected_prices.shape != (n,) or not np.isfinite(projected_prices).all()
                or np.any(projected_prices <= 0)):
            raise ValueError("Price projection must return a positive finite vector of unchanged shape")
        return np.r_[projected_prices, np.clip(values[n:], lower, upper)]

    def evaluate_joint(values):
        nonlocal evaluation_count
        evaluation_count += 1
        trial_prices, trial_fiscal = decode(values)
        try:
            reply = evaluate(trial_prices.copy(), trial_fiscal.copy())
        except CandidateDomainError as exc:
            domain_rejections[evaluation_count] = str(exc)
            return dict(residual=np.full(2*n, np.nan), mapping_valid=False,
                        payload={"candidate_domain_rejection": str(exc)})
        if not isinstance(reply, Mapping):
            raise ValueError("Economic evaluator must return a mapping")
        market = np.asarray(reply["market_residual"], dtype=float)
        budget = np.asarray(reply["fiscal_residual"], dtype=float)
        if (market.shape != (n,) or budget.shape != (n,)
                or not isinstance(reply["mapping_valid"], (bool, np.bool_))):
            raise ValueError("Evaluator returned invalid market/fiscal shapes or non-Boolean mapping gate")
        return dict(residual=np.r_[market, budget] * row_scale,
                    mapping_valid=bool(reply["mapping_valid"]), payload=reply.get("payload"))

    def physical_record(record):
        if record is None:
            return None
        result = copy.deepcopy(record)
        if record.get("evaluation") in domain_rejections:
            result["candidate_domain_rejection"] = domain_rejections[record["evaluation"]]
        if "prices" in record:
            result["prices"], result["fiscal_values"] = decode(record["prices"])
            result["fiscal_unknown"] = "period_pension" if closure == "fixed_tax" else "payroll_tax"
        if "x" in record:
            result["x"] = np.log(np.r_[result["prices"], result["fiscal_values"]])
        if "residual" in record:
            residual = np.asarray(record["residual"], dtype=float) / row_scale
            result["residual"] = residual
            result["market_residual"] = residual[:n].copy()
            result["fiscal_residual"] = residual[n:].copy()
            result["market_gate"] = bool(np.isfinite(residual[:n]).all()
                and np.max(np.abs(residual[:n])) <= market_tolerance)
            result["fiscal_gate"] = bool(np.isfinite(residual[n:]).all()
                and np.max(np.abs(residual[n:])) <= fiscal_tolerance)
        for name in ("score", "best_score"):
            if name in record and record[name] is not None:
                result[name] = float(record[name]) / market_tolerance
        result["closure"] = closure
        return result

    jacobian = None
    if initial_jacobian is not None:
        jacobian = np.asarray(initial_jacobian, dtype=float)
        if jacobian.shape != (2*n, 2*n) or not np.isfinite(jacobian).all():
            raise ValueError("Initial Jacobian must be a finite 2N by 2N physical-log-coordinate matrix")
        jacobian = row_scale[:, None] * jacobian

    def progress(record):
        # Completion is published below, after checking both unscaled replay
        # blocks. The generic core has only one tolerance for its scaled vector.
        if callback is not None and record.get("event") != "complete":
            callback(physical_record(record))

    result = solve_price_path(
        initial_prices=np.r_[prices, fiscal], evaluate=evaluate_joint,
        project=project, slope=float(market_slope), market_tolerance=float(market_tolerance),
        max_log_step=max_log_step, damping=damping, max_evaluations=max_evaluations,
        deadline_monotonic=deadline_monotonic, max_condition_number=max_condition_number,
        worsening_factor=worsening_factor,
        final_reproduction_tolerance=float(final_reproduction_tolerance) * max(1.0, residual_scale),
        callback=progress,
        initial_jacobian=jacobian, default_jacobian=default_jacobian,
    )
    best, final = physical_record(result["best"]), physical_record(result["final"])
    gaps = None if best is None or final is None else np.abs(final["residual"] - best["residual"])
    market_replay = None if gaps is None else float(np.max(gaps[:n]))
    fiscal_replay = None if gaps is None else float(np.max(gaps[n:]))
    gates = dict(
        mapping=bool(final is not None and final["mapping_valid"]),
        housing=bool(final is not None and final["market_gate"]),
        social_security=bool(final is not None and final["fiscal_gate"]),
        market_replay=bool(market_replay is not None and np.isfinite(market_replay)
                          and market_replay <= final_reproduction_tolerance),
        fiscal_replay=bool(fiscal_replay is not None and np.isfinite(fiscal_replay)
                          and fiscal_replay <= final_reproduction_tolerance),
    )
    final_jacobian = np.asarray(result["final_jacobian"]) / row_scale[:, None]
    result.update(best=best, final=final,
        history=[physical_record(row) for row in result["history"]],
        final_jacobian=final_jacobian, closure=closure, gates=gates,
        domain_rejections=dict(domain_rejections),
        converged=bool(result["converged"] and all(gates.values())),
        market_reproduction_max_abs=market_replay, fiscal_reproduction_max_abs=fiscal_replay,
        final_reproduction_max_abs=None if gaps is None else float(np.max(gaps)),
        contract=dict(closure=closure, fiscal_bounds=(lower, upper),
            market_tolerance=float(market_tolerance), fiscal_tolerance=float(fiscal_tolerance),
            market_slope=float(market_slope), fiscal_slope=float(fiscal_slope),
            fiscal_residual_scale=residual_scale, max_log_step=float(max_log_step),
            final_reproduction_tolerance=float(final_reproduction_tolerance),
            score_definition="maximum absolute residual divided by its own tolerance"))
    if result["status"] == "converged" and not result["converged"]:
        result["status"] = "final_joint_gate_failed"
    if callback is not None:
        callback(dict(event="complete", converged=result["converged"],
            status=result["status"], evaluations=result["evaluations"], closure=closure,
            gates=copy.deepcopy(gates)))
    return result
