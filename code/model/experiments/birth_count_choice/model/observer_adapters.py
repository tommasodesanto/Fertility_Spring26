"""Recorded, in-memory observer clones installed only after frozen authentication.

No frozen source file is edited. All probability inputs come from the observed
policy's owned count snapshot. Uniform within-cell interpolation is retained.
"""
from __future__ import annotations
import hashlib
import inspect
import json
import sys
from pathlib import Path
import numpy as np
from .engine.birth_count import birth_count_transition


def apply_count_fertility(pre, first, P, continuation=None, *, birth_count_realized_probs=None):
    probabilities = birth_count_realized_probs
    if probabilities is None:
        probabilities = getattr(P, 'birth_count_realized_probs', None)
    if probabilities is None:
        raise RuntimeError('Count observer requires owned realized-count probabilities')
    result = birth_count_transition(pre, probabilities,
                                    cap=getattr(P, 'birth_count_choice_cap', 3))
    locations = np.sum(result['expected_births_by_cell'], axis=(0, 1, 3, 4))
    return result['post'], result['expected_births'], locations


def first_birth_accounting_by_age(evaluation, P):
    pre = evaluation.g_pre
    probabilities = evaluation.policy.birth_count_realized_probs
    flows = np.zeros(int(P.J)); risk = np.zeros_like(flows)
    for j in range(int(P.J)):
        result = birth_count_transition(pre[:, :, :, j], probabilities[:, :, :, j],
                                        cap=getattr(P, 'birth_count_choice_cap', 3))
        flows[j] = result['births_by_order'][0]
        risk[j] = np.sum(pre[:, :, :, j, :, 0, :])
    hazard = np.divide(flows, risk, out=np.zeros_like(flows), where=risk > 1e-15)
    if np.any(hazard < -1e-14) or np.any(hazard > 1 + 1e-12):
        raise RuntimeError('Invalid count-menu first-birth household hazard')
    return {'flow': flows, 'at_risk': risk, 'hazard': np.clip(hazard, 0, 1)}


def count_flows_and_risks(evaluation, P):
    flows = np.zeros((int(P.J), 3)); risk = np.zeros_like(flows)
    for j in range(int(P.J)):
        result = birth_count_transition(evaluation.g_pre[:, :, :, j],
                                      evaluation.policy.birth_count_realized_probs[:, :, :, j],
                                      cap=getattr(P, 'birth_count_choice_cap', 3))
        if np.max(np.abs(result['post'] - evaluation.g_post_fertility[:, :, :, j])) > 2e-10:
            raise RuntimeError('Count observer does not replay its evaluation')
        flows[j] = result['births_by_order']; risk[j] = result['at_risk_by_order']
    return flows, risk


def _once(source, old, new, name):
    if source.count(old) != 1:
        raise RuntimeError(f'Observer source anchor drift ({name}): {old[:80]!r}')
    return source.replace(old, new, 1)


def _clone(function, out, changes, injected, receipt):
    """Keep wrapper-held function identities, replacing only copied code objects."""
    original_path = inspect.getsourcefile(function)
    before = inspect.getsource(function)
    after = before
    for old, new in changes:
        after = _once(after, old, new, function.__name__)
    path = Path(out) / ('birth_count_' + function.__name__ + '.py')
    path.write_text(after)
    namespace = function.__globals__
    namespace.update(injected)
    exec(compile(after, str(path), 'exec'), namespace)
    replacement = namespace[function.__name__]
    function.__code__ = replacement.__code__
    namespace[function.__name__] = function
    receipt.append(dict(function=function.__name__, original_source=original_path,
                        clone_source=str(path), original_sha256=hashlib.sha256(before.encode()).hexdigest(),
                        clone_sha256=hashlib.sha256(after.encode()).hexdigest(),
                        exact_replacements=[dict(before=a, after=b) for a, b in changes]))
    return function


def install_birth_count_observers(context, facade, out):
    rt = context['prepared'].rt
    if getattr(rt['primitive'].pf.calendar, '_birth_count_observers_installed', False):
        return
    cal = rt['primitive'].pf.calendar
    transition = rt['primitive'].pf.transition
    calibration = rt['calibration']
    receipt = []
    original_factory = cal.policy_from_solution
    def policy_from_solution(solution, price, P, grid, shared):
        policy = original_factory(solution, price, P, grid, shared)
        for key in ('birth_count_action_probs', 'birth_count_realized_probs'):
            setattr(policy, key, np.asarray(getattr(solution, key)).copy())
        return policy
    cal.policy_from_solution = policy_from_solution
    cal.apply_fertility = transition.apply_sequential_fertility = apply_count_fertility
    # Both active calendar replay paths receive this policy's owned snapshot.
    _clone(cal.evaluate_period, out, [(
        '            gated, policy.fert_probs, P, continuation\n',
        '            gated, policy.fert_probs, P, continuation,\n'
        '            birth_count_realized_probs=policy.birth_count_realized_probs,\n')], {}, receipt)
    _clone(cal.reconstruct_stationary_pre_fertility, out, [(
        '            g_pre, policy.fert_probs, P, policy_continuation_birth_probs(policy, P)\n',
        '            g_pre, policy.fert_probs, P, policy_continuation_birth_probs(policy, P),\n'
        '            birth_count_realized_probs=policy.birth_count_realized_probs,\n')], {}, receipt)
    # Preserve all inherited rate, age, top-code and target arithmetic, replacing
    # only the one-birth flow loop and its binary risk bound.
    source = inspect.getsource(calibration.period_fertility_diagnostics)
    start = source.index('    for j in range(int(P.J)):\n')
    end = source.index('    explicit_by_age =', start)
    _clone(calibration.period_fertility_diagnostics, out, [(
        source[start:end], '    flows, count_risk = _count_flows_and_risks(evaluation, P)\n'),
        ('        "birth_flow_first": flows[:, 0],\n',
         '        "birth_order_risk": count_risk,\n        "birth_flow_first": flows[:, 0],\n')],
        {'_count_flows_and_risks': count_flows_and_risks}, receipt)
    calibration.first_birth_accounting_by_age = first_birth_accounting_by_age
    # The frozen initial-fertility function imports helpers at call time.
    initial = next((m for m in sys.modules.values() if str(getattr(m, '__file__', '')).endswith('e5f_initial_fertility_observer.py')), None)
    if initial is None:
        raise RuntimeError('Authenticated initial fertility module absent')
    _clone(initial.observe_initial_fertility, out, [(
        '            or np.any(flows > pre_parity[:, :3] + FLOW_ATOL)):',
        '            or np.any(flows > np.asarray(period["birth_order_risk"]) + FLOW_ATOL)):'),
        ('"one parity transition per cell, uniformly distributed within the four-year age interval"',
         '"count-menu pre/post stocks interpolated uniformly within the four-year age interval; ordered within-cell dates unavailable"')], {}, receipt)
    # First-birth event starts from origin n=0; destinations retain n=m=X.
    fn = calibration.first_birth_housing_response
    _clone(fn, out, [(
        '            attempt = policy.fert_probs[:, :, :, j, zz, 1]\n'
        '            realized = float(fecundity[j]) * childless * attempt\n',
        '            origin = np.zeros(shape)\n'
        '            origin[:, :, :, zz, 0, settled] = childless\n'
        '            tagged = _count_transition(origin, policy.birth_count_realized_probs[:, :, :, j], cap=P.birth_count_choice_cap)["first_birth_tagged_post"]\n'
        '            realized = np.sum(tagged, axis=(-2, -1))[:, :, :, zz]\n'),
        ('            birth_cohort[:, :, :, zz, 1, 1] = realized\n',
         '            birth_cohort = tagged\n'),
        ('    second birth inside the four-year window.\n',
         '    further fertility choice inside the four-year window; origin births retain their actual count mixture.\n')],
        {'_count_transition': birth_count_transition}, receipt)
    _clone(calibration.finish_dated_first_birth_housing_branch, out, [(
        '            calendar.policy_continuation_birth_probs(evaluation.policy, P),\n        )\n        control_post = control_pre',
        '            calendar.policy_continuation_birth_probs(evaluation.policy, P),\n'
        '            birth_count_realized_probs=evaluation.policy.birth_count_realized_probs,\n'
        '        )\n        control_post = control_pre')], {}, receipt)
    from .engine.diagnostics import write_diagnostics
    rt['audit'].model = facade
    rt['audit'].write_diagnostics = write_diagnostics
    _clone(rt['audit'].standard_diagnostics, out, [(
        '    stats.fert2_probs = calendar.policy_continuation_birth_probs(p, P)\n',
        '    stats.fert2_probs = calendar.policy_continuation_birth_probs(p, P)\n'
        '    stats.birth_count_action_probs = p.birth_count_action_probs.copy()\n'
        '    stats.birth_count_realized_probs = p.birth_count_realized_probs.copy()\n'
        '    stats.birth_count_pre_distribution = e.g_pre.copy()\n')], {}, receipt)
    # Kept for the inherited dated branch implementation, though the maintained
    # independent-count SS observer takes the same-policy branch above.
    fn = calibration.begin_dated_first_birth_housing_branch
    _clone(fn, out, [(
        '                attempt = policy.fert_probs[:, :, :, j, zz, 1]\n'
        '                realized = float(fecundity[j]) * childless * attempt\n'
        '                treated[:, :, :, j, zz, 1, 1] = realized\n',
        '                origin = np.zeros_like(evaluation.g_pre[:, :, :, j])\n'
        '                origin[:, :, :, zz, 0, settled] = childless\n'
        '                tagged = _count_transition(origin, policy.birth_count_realized_probs[:, :, :, j], cap=P.birth_count_choice_cap)["first_birth_tagged_post"]\n'
        '                realized = np.sum(tagged, axis=(-2, -1))[:, :, :, zz]\n'
        '                treated[:, :, :, j] += tagged\n')], {'_count_transition': birth_count_transition}, receipt)
    wrapped = rt['observe_recent_parent_flow']
    recent = next((cell.cell_contents for cell in (wrapped.__closure__ or ())
                   if hasattr(cell.cell_contents, 'observe_recent_parent_flow')), None)
    if recent is None:
        raise RuntimeError('Recent observer wrapper shape drift')
    _clone(recent.observe_recent_parent_flow, out, [(
        '        return transition.apply_sequential_fertility(mass, first, P, continuation)',
        '        return transition.apply_sequential_fertility(mass, first, P, continuation,\n'
        '            birth_count_realized_probs=policy.birth_count_realized_probs)'),
        ('    birth_post = empty_post * ((parity > 0) & (child_state == 1))\n',
         '    birth_post = empty_post * ((parity > 0) & (child_state > 0))\n'
         '    selected_households = float(np.sum(empty_pre * np.sum(policy.birth_count_realized_probs[..., 1:], axis=-1)))\n'),
        ('selected_birth_flow_error=abs(float(birth_post.sum()) - selected_births)',
         'selected_birth_flow_error=abs(float(birth_post.sum()) - selected_households)'),
        ('    weights = uniform_age_cell_overlap(P, 30., 56.)\n',
         '    first_post = _count_transition(pre * never, policy.birth_count_realized_probs, cap=P.birth_count_choice_cap)["first_birth_tagged_post"]\n'
         '    former_post, _, _ = fertility(pre * former)\n'
         '    first_current = transport(first_post)\n'
         '    continuation_current = transport(former_post * ((parity > 0) & (child_state > 0)))\n'
         '    _check(_error(first_current + continuation_current, birth_current), "origin-tag partition")\n'
         '    weights = uniform_age_cell_overlap(P, 30., 56.)\n'),
        ('        first_birth=_group(birth_current, weights, parity == 1),\n'
         '        continuation_birth=_group(birth_current, weights, parity >= 2),\n',
         '        first_birth=_group(first_current, weights),\n'
         '        continuation_birth=_group(continuation_current, weights),\n'),
        ('"supplied policy.fert_probs and owned policy.fert2_probs"',
         '"owned policy.birth_count_realized_probs; first/continuation groups tagged by origin"'),
        ('"At most one explicit birth per four-year period; top-code representative does not scale household mass"',
         '"Up to the configured birth-count cap per period; successful households counted once; first-birth families retain actual count mixtures"')],
        {'_count_transition': birth_count_transition}, receipt)
    cal._birth_count_observers_installed = True
    context['birth_count_observer_adapters'] = receipt
    Path(out, 'birth_count_observer_receipt.json').write_text(json.dumps(dict(
        frozen_sources_edited=False, authentication_before_adaptation=True,
        acceptance_tolerances_unchanged=True, adapters=receipt,
        maintained_conventions=['post-interest transactions', 'soft credit', 'top-count correction',
                                'linear within-age-cell stock interpolation', 'unchanged empirical targets and weights'],
        revised_conventions=['owned Binomial count kernel', 'crossed-order flows',
                             'first-birth origin tags and actual destination mixture',
                             'birth-household weights counted once']), indent=2) + '\n')
