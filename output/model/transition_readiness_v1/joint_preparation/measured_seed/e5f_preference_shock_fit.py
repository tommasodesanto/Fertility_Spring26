"""Scalar shock fitting and sequential surprises; no model or economic defaults."""
from __future__ import annotations

import time
import numpy as np
from e5f_ssj_scaled_step_root import solve_price_path_scaled


def fit_rows(targets, models, kind):
    if kind not in ('one_permanent', 'four_successive') or len(targets) != 4 or len(models) != 4:
        raise ValueError('One/four shocks require the complete four-window readout')
    result=[]
    for i,(target,model) in enumerate(zip(targets,models)):
        weight=1. if kind=='four_successive' or i==3 else 0.
        gap=None if model is None else float(model)-float(target['target'])
        if gap is not None and not np.isfinite(gap):raise ValueError('Nonfinite fertility measurement')
        result.append(dict(target,model=model,gap=gap,weight=weight,
            loss_contribution=None if gap is None else weight*gap**2,
            role='fitted' if weight else 'validation'))
    return result


def fit_one(*, evaluate, target, initial_level, bounds, controls, callback=None):
    """Fit one permanent level to one moment, conditional on the inherited state.

    The evaluator accepts ONLY a scalar current preference, returns ``certified``
    and ``model``, and solves its own PF equilibrium. Invalid equilibria never
    supply residuals. A central log derivative initializes the bounded root;
    convergence requires a fresh reproduction of the final equilibrium.
    """
    x=float(initial_level);limits=np.asarray(bounds,float)
    if limits.shape!=(2,) or not np.isfinite(np.r_[x,limits,float(target)]).all() or not 0<limits[0]<x<limits[1]:
        raise ValueError('Positive finite ordered bounds with an interior initial level required')
    maximum=controls['max_evaluations'];step=float(controls['log_difference_step'])
    if type(maximum) is not int or maximum<5 or not 0<step<.25:
        raise ValueError('Budget must include derivative samples and a fresh final replay')
    seconds=float(controls['total_seconds'])
    if not np.isfinite(seconds) or seconds<=0:raise ValueError('Finite positive fit time budget required')
    deadline=time.monotonic()+seconds;count=0
    def sample(levels):
        nonlocal count
        if count>=maximum or time.monotonic()>=deadline:raise TimeoutError('Shock-fit budget exhausted')
        count+=1;reply=evaluate(float(levels[0]))
        if not reply['certified']:
            return dict(mapping_valid=False,residual=np.array([np.nan]),payload=reply.get('payload'))
        model=float(reply['model'])
        if not np.isfinite(model):raise ValueError('Finite fertility observation required')
        return dict(mapping_valid=True,residual=np.array([float(target)-model]),
            payload=dict(reply.get('payload',{}),model=model,outer_evaluation=count))
    base=sample([x])
    if not base['mapping_valid']:raise RuntimeError('Initial shock proposal has no certified equilibrium')
    minus=max(limits[0],x*np.exp(-step));plus=min(limits[1],x*np.exp(step))
    lo,hi=sample([minus]),sample([plus])
    if not lo['mapping_valid'] or not hi['mapping_valid']:
        raise RuntimeError('Uncertified derivative proposal; stop without scoring its fertility')
    derivative=float((hi['residual'][0]-lo['residual'][0])/np.log(plus/minus))
    if abs(derivative)<=1e-10:raise RuntimeError('Shock is locally underidentified at the initial proposal')
    root=solve_price_path_scaled(initial_prices=np.array([x]),evaluate=sample,
        project=lambda values:np.clip(values,limits[0],limits[1]),slope=1.,
        market_tolerance=controls['fertility_tolerance'],max_log_step=controls['max_log_step'],
        damping=controls['damping'],max_evaluations=maximum-count,deadline_monotonic=deadline,
        max_condition_number=controls['max_condition_number'],worsening_factor=controls['worsening_factor'],
        final_reproduction_tolerance=controls['reproduction_tolerance'],
        initial_jacobian=np.array([[derivative]]),callback=callback)
    final=root['final'];parameter=None
    if final is not None and final['mapping_valid']:
        value=float(final['prices'][0]);lower,upper=map(float,limits)
        parameter=dict(estimate=value,lower=lower,upper=upper,
            near_bound=min(value-lower,upper-value)<=.01*(upper-lower))
    return dict(status='matched' if root['converged'] else 'not_matched',converged=root['converged'],
        root=root,parameter=parameter,equilibrium_evaluations=count,
        initial_fertility_derivative_log_psi=-derivative)


def fit_sequence(*, evaluate_factory, advance, targets, initial_level, bounds, controls, callback=None):
    """Estimate four surprises in order; advance only after each certified fit.

    Later shocks are absent from the evaluator API. ``advance`` implements the
    first period under the accepted forecast and preserves its complete state.
    Each stage has its own evaluation cap and shares one overall deadline.
    """
    if len(targets)!=4:raise ValueError('Exactly four successive fertility targets required')
    deadline=time.monotonic()+controls['total_seconds'];stages=[];level=float(initial_level)
    for stage,target in enumerate(targets):
        result=fit_one(evaluate=evaluate_factory(stage),target=target['target'],initial_level=level,
            bounds=bounds,controls=dict(controls,total_seconds=deadline-time.monotonic()),
            callback=None if callback is None else lambda row:callback(stage,row))
        stages.append(result)
        if not result['converged']:break
        advance(stage,result)
        level=result['parameter']['estimate']
    return dict(status='matched' if len(stages)==4 and all(s['converged'] for s in stages) else 'not_matched',
        converged=len(stages)==4 and all(s['converged'] for s in stages),stages=stages,
        equilibrium_evaluations=sum(s['equilibrium_evaluations'] for s in stages))
