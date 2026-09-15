"""Isolated diagnostic copy of ``e5f_matched_pf_path_root.solve_price_path``.

The ONLY change is the step-length rule: instead of clipping every coordinate
of the damped Newton step to ``[-max_log_step, max_log_step]`` independently,
the whole step is scaled by a common factor so that its largest coordinate
equals ``max_log_step``.  The direction is preserved, which matters when the
supplied Jacobian couples dates and blocks.  Budget, projection, Broyden
update, safeguards, and final certification are byte-identical to the retained
solver.  Not a production file.
"""
from __future__ import annotations

import copy
import time
import numpy as np


def solve_price_path_scaled(*, initial_prices, evaluate, project, slope,
                     market_tolerance, max_log_step, damping, max_evaluations,
                     deadline_monotonic, max_condition_number, worsening_factor,
                     final_reproduction_tolerance, callback=None,
                     initial_jacobian=None, default_jacobian=None,
                     tolerance_vector=None):
    """Return a dict; only ``converged=True`` certifies the fresh final mapping.

    ``evaluate`` returns ``residual``, Boolean ``mapping_valid`` and optionally
    ``payload``. The final evaluation is reserved within max_evaluations. The
    deadline is checked before every call; a caller watchdog bounds each call.
    ``callback`` receives a copied ledger record after every evaluation, with
    ``new_best`` and the current best score, plus the final completion record.
    ``initial_jacobian`` optionally supplies d(residual)/d(log price), a finite
    NxN matrix. ``default_jacobian`` optionally specifies a nonsingular reset
    matrix for differently scaled equations; the default remains ``-slope * I``.
    ``tolerance_vector`` (diagnostic copy only) gives a per-coordinate gate; when
    supplied the ledger ``score`` is ``max_i |r_i| / tol_i`` and convergence
    means ``score <= 1``, while ``raw_max_abs`` keeps the retained max-abs
    residual for comparability.  ``market_tolerance`` is then unused.
    """
    p0 = np.asarray(initial_prices, dtype=float)
    if p0.ndim != 1 or not p0.size or not np.isfinite(p0).all() or np.any(p0 <= 0):
        raise ValueError('Initial prices must be a positive finite vector')
    limits = [slope, market_tolerance, max_log_step, damping, max_condition_number,
              worsening_factor, deadline_monotonic, final_reproduction_tolerance]
    if (not np.isfinite(limits).all() or min(limits[:5]) <= 0 or damping > 1
            or max_condition_number <= 1 or worsening_factor <= 1
            or final_reproduction_tolerance < 0 or type(max_evaluations) is not int
            or max_evaluations < 2 or deadline_monotonic <= time.monotonic()):
        raise ValueError('Invalid explicit numerical/time/evaluation limits')
    if not callable(evaluate) or not callable(project) or (callback is not None and not callable(callback)):
        raise ValueError('Evaluation, projection and optional callback must be callable')
    if tolerance_vector is not None:
        tolerance_vector = np.asarray(tolerance_vector, dtype=float)
        if tolerance_vector.shape != p0.shape or not np.isfinite(tolerance_vector).all() or np.any(tolerance_vector <= 0):
            raise ValueError('tolerance_vector must be a positive finite vector matching the coordinates')
        gate = 1.0
    else:
        gate = float(market_tolerance)
    started = time.monotonic()
    ledger, best, evaluation_count = [], None, 0
    J0 = -float(slope) * np.eye(p0.size)
    if default_jacobian is not None:
        J0 = np.asarray(default_jacobian, dtype=float).copy()
        if J0.shape != (p0.size, p0.size) or not np.isfinite(J0).all():
            raise ValueError('Default Jacobian must be a finite nonsingular NxN matrix')
        try:
            if not np.isfinite(np.linalg.solve(J0, np.ones(p0.size))).all():
                raise np.linalg.LinAlgError('nonfinite reset step')
        except np.linalg.LinAlgError as exc:
            raise ValueError('Default Jacobian must be a finite nonsingular NxN matrix') from exc
    J = J0.copy() if initial_jacobian is None else np.asarray(initial_jacobian, dtype=float).copy()
    if J.shape != J0.shape or not np.isfinite(J).all():
        raise ValueError('Initial Jacobian must be a finite NxN matrix for log-price residual derivatives')
    active_damping = float(damping)

    def projected(prices):
        values = np.asarray(project(np.asarray(prices, dtype=float).copy()), dtype=float)
        if values.shape != p0.shape or not np.isfinite(values).all() or np.any(values <= 0):
            raise ValueError('Projection must return a positive finite price vector of unchanged shape')
        return values.copy()

    def emit(record):
        if callback is not None:
            callback(copy.deepcopy(record))

    def sample(prices, *, phase, **diagnostics):
        nonlocal best, evaluation_count
        if evaluation_count >= max_evaluations or time.monotonic() >= deadline_monotonic:
            raise TimeoutError('PF root evaluation/time budget exhausted')
        began = time.monotonic()
        evaluation_count += 1
        reply = evaluate(prices.copy())
        residual = np.asarray(reply['residual'], dtype=float)
        if residual.shape != p0.shape or not isinstance(reply['mapping_valid'], (bool, np.bool_)):
            raise ValueError('Evaluator returned invalid residual shape or non-Boolean mapping gate')
        valid = bool(reply['mapping_valid']) and bool(np.isfinite(residual).all())
        raw_max_abs = float(np.max(np.abs(residual))) if valid else float('inf')
        score = (float(np.max(np.abs(residual) / tolerance_vector)) if (valid and tolerance_vector is not None)
                 else raw_max_abs)
        point = dict(prices=prices.copy(), x=np.log(prices), residual=residual.copy(),
                     score=score, raw_max_abs=raw_max_abs, mapping_valid=valid,
                     payload=copy.deepcopy(reply.get('payload')))
        improved = phase != 'final' and valid and (best is None or score < best['score'])
        if improved:
            best = copy.deepcopy(point)
        record = dict(evaluation=evaluation_count, phase=phase, prices=prices.copy(),
            residual=residual.copy(), score=score, raw_max_abs=raw_max_abs, mapping_valid=valid,
            new_best=improved, best_score=None if best is None else best['score'],
            elapsed_seconds=time.monotonic() - started,
            evaluation_seconds=time.monotonic() - began, **diagnostics)
        ledger.append(record)
        emit(record)
        return point

    reason, final = 'evaluation_budget', None
    try:
        current = sample(projected(p0), phase='initial')
        if not current['mapping_valid']:
            reason = 'invalid_initial_mapping'
        while best is not None and best['score'] > gate and evaluation_count < max_evaluations - 1:
            reset_reason = None
            try:
                if not np.isfinite(J).all() or np.linalg.cond(J) > max_condition_number:
                    raise np.linalg.LinAlgError('ill-conditioned Jacobian')
                raw_step = -np.linalg.solve(J, current['residual'])
                if not np.isfinite(raw_step).all():
                    raise np.linalg.LinAlgError('nonfinite Newton step')
            except np.linalg.LinAlgError:
                J = J0.copy()
                raw_step = (current['residual'] / float(slope) if default_jacobian is None
                            else -np.linalg.solve(J0, current['residual']))
                reset_reason = 'jacobian_reset_positive_residual_price_increase'
            step = active_damping * raw_step
            largest = float(np.max(np.abs(step)))
            if largest > max_log_step:
                step = step * (max_log_step / largest)
            trial_prices = projected(np.exp(current['x'] + step))
            actual_step = np.log(trial_prices) - current['x']
            if float(actual_step @ actual_step) <= np.finfo(float).tiny:
                reason = 'projection_stalled_without_market_gate'
                break
            previous_best_score = best['score']
            trial = sample(trial_prices, phase='iterate', damping=active_damping,
                requested_log_step=step.copy(), actual_log_step=actual_step.copy(), reset_reason=reset_reason)
            if not trial['mapping_valid'] or trial['score'] > worsening_factor * previous_best_score:
                current = copy.deepcopy(best)
                active_damping *= .5
                J = J0.copy()
                ledger[-1]['safeguard'] = 'restored_best_halved_damping_reset_jacobian'
                emit(dict(ledger[-1], event='safeguard'))
                continue
            y = trial['residual'] - current['residual']
            J += np.outer(y - J @ actual_step, actual_step) / float(actual_step @ actual_step)
            current = trial
        if best is not None:
            if best['score'] <= gate:
                reason = 'candidate_market_gate'
            # Deliberately do not project again: certify these identical prices.
            final = sample(best['prices'], phase='final')
    except TimeoutError:
        reason = 'time_or_evaluation_budget'
    reproduction = None if final is None or best is None else float(np.max(np.abs(final['residual'] - best['residual'])))
    converged = bool(final is not None and best is not None and best['score'] <= gate
        and final['mapping_valid'] and final['score'] <= gate
        and np.isfinite(reproduction) and reproduction <= final_reproduction_tolerance)
    if final is not None and best is not None and best['score'] <= gate and not converged:
        reason = ('final_mapping_failed' if not final['mapping_valid'] else
                  'final_market_gate_failed' if final['score'] > gate else
                  'final_reproduction_failed')
    result = dict(converged=converged, status='converged' if converged else reason,
        best=best, final=final, final_reproduction_max_abs=reproduction,
        evaluations=evaluation_count, elapsed_seconds=time.monotonic() - started,
        history=ledger, final_jacobian=J.copy(), final_damping=active_damping,
        tolerance_vector=None if tolerance_vector is None else tolerance_vector.copy(), gate=gate)
    emit(dict(event='complete', converged=converged, status=result['status'], evaluations=evaluation_count))
    return result
