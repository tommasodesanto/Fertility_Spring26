"""Opt-in numerical initialization adapter; no household or solver equation edits."""
import copy
import math


def apply(evaluator, original_contract, initial_psi):
    """Call only after unmodified runtime.setup authenticated original sources."""
    if evaluator.c != original_contract:
        raise RuntimeError('Runtime contract changed before numerical adapter')
    if not math.isfinite(initial_psi) or initial_psi <= 0:
        raise ValueError('Positive finite diagnostic initial psi required')
    changed = copy.deepcopy(original_contract)
    changed['normalization']['initial_psi'] = initial_psi
    restored = copy.deepcopy(changed)
    restored['normalization']['initial_psi'] = original_contract['normalization']['initial_psi']
    if restored != original_contract:
        raise RuntimeError('Adapter changed more than normalization initial guess')
    evaluator.c = changed
    return changed
