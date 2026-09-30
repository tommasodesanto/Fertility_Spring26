"""Narrow bridge to the sibling two-price debt-rule packet.

This module deliberately owns no preference, grid, or credit primitive.  It
only applies the two overrides already defined in ``credit_rule_quick_v1``.
``prepare`` returns the companion packet's exact 160-node grid and
``apply_rule`` mutates only its documented three debt objects and, for
``author``, installs the companion's compiled sale gate.
"""
from pathlib import Path
import importlib.util
import sys
import numpy as np

SIBLING = Path(__file__).resolve().parents[1] / 'credit_rule_quick_v1'

def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec); sys.modules[name] = mod
    spec.loader.exec_module(mod); return mod

def prepare(reference):
    return np.asarray(reference['b_grid']).copy()

def apply_rule(model, P, case):
    runner = _load(SIBLING / 'run_credit_rule_quick.py', '_credit_rule_quick_source')
    P.lambda_d = 0.; P.debt_taper_weights = np.zeros(P.J + 1) if case == 'author' else runner.weights(P)
    P.debt_caps = np.zeros(P.J + 1)
    if case == 'author':
        strict = _load(SIBLING / 'strict_tenure.py', '_credit_rule_quick_strict')
        model.tenure_logit_kernel = strict.tenure_logit_kernel
        model.tenure_choice_kernel = strict.tenure_choice_kernel
    elif case != 'ours':
        raise RuntimeError('Unknown debt rule: ' + str(case))
    return {'lambda_d': 0., 'debt_taper_weights': P.debt_taper_weights.tolist(),
            'debt_caps': P.debt_caps.tolist(), 'strict_sale_gate': case == 'author'}
