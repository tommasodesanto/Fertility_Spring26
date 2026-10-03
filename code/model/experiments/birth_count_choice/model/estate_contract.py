"""Estate-A flags and wealth-only rescoring; native old-target reports are retained."""
from __future__ import annotations
import copy
import csv
import json
import math
from pathlib import Path
import numpy as np
from . import calibration

OLD_WEALTH_TARGET = '6.92658379107299'
NEW_WEALTH_TARGET = '4.45838713455674'
OBSERVER_CONTRACT = 'estate_a_postsaving_net_selling_cost_v1'
FLAG_NAMES = {'birth_count_choice_enabled', 'birth_count_choice_cap',
              'bequest_net_of_selling_cost', 'estate_flow_net_of_selling_cost'}


def experiment_flags(birth_cap):
    if isinstance(birth_cap, bool) or birth_cap not in (1, 3):
        raise ValueError('Estate-A birth cap must be 1 or 3')
    return dict(birth_count_choice_enabled=True, birth_count_choice_cap=birth_cap,
                bequest_net_of_selling_cost=True, estate_flow_net_of_selling_cost=True)


def apply_experiment_flags(P, flags):
    if not isinstance(flags, dict) or set(flags) != FLAG_NAMES:
        raise ValueError('Estate-A requires exactly the four whitelisted flags')
    if flags != experiment_flags(flags['birth_count_choice_cap']):
        raise ValueError('Estate-A flags must enable both net-estate margins')
    if any(type(flags[k]) is not bool for k in FLAG_NAMES - {'birth_count_choice_cap'}):
        raise ValueError('Estate-A switches require Boolean values')
    if type(flags['birth_count_choice_cap']) is not int:
        raise ValueError('Estate-A birth cap must be an integer')
    if P is not None:
        for key, value in flags.items(): setattr(P, key, value)
        if str(getattr(P, 'estate_receiver', 'none')).lower() != 'none':
            raise ValueError('Estate-A does not authorize a recipient mapping')
    return P


def contract():
    base, _, _, bounds = calibration._contract()
    rows = copy.deepcopy(base)
    wealth = [row for row in rows if row['moment'] == 'wealth_earnings']
    if len(rows) != 14 or len(wealth) != 1 or wealth[0]['target'] != OLD_WEALTH_TARGET:
        raise RuntimeError('Authenticated old wealth target drift')
    if wealth[0]['weight'] != '7.595098472533724':
        raise RuntimeError('Wealth weight drift')
    wealth[0]['target'] = NEW_WEALTH_TARGET
    return (rows, calibration._canonical(rows),
            calibration._canonical(dict(base_contract=rows, multipliers={})),
            dict(bounds, beta_annual=(.93, .99)))


def _write_rows(path, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def rescore_rows(rows):
    base = calibration._contract()[0]
    identity = [{k: r[k] for k in ('moment','target','weight','role')} for r in rows]
    if identity != base: raise RuntimeError('Native complete old target contract drift')
    revised = copy.deepcopy(rows)
    for r in revised:
        if r['role'] == 'scored':
            gap = float(r['model']) - float(r['target'])
            if not math.isclose(gap, float(r['gap']), abs_tol=1e-10): raise RuntimeError('Native gap arithmetic drift')
            if not math.isclose(float(r['weight'])*gap*gap, float(r['loss_contribution']), abs_tol=1e-8): raise RuntimeError('Native loss arithmetic drift')
        if r['moment'] == 'wealth_earnings':
            r['target'] = NEW_WEALTH_TARGET
            r['gap'] = str(float(r['model']) - float(NEW_WEALTH_TARGET))
            r['loss_contribution'] = str(float(r['weight']) * float(r['gap'])**2)
    expected, target_pin, weight_pin, _ = contract()
    if [{k:r[k] for k in ('moment','target','weight','role')} for r in revised] != expected:
        raise RuntimeError('New complete target contract drift')
    residual = np.asarray([math.sqrt(float(r['weight']))*float(r['gap']) for r in revised if r['role']=='scored'])
    if residual.shape != (10,) or not np.isfinite(residual).all(): raise RuntimeError('New residual invalid')
    return revised, residual


def rescore_report(report_directory, destination=None):
    report = Path(report_directory); out = Path(destination or report)
    old_residual, old_rows, parameters = calibration.residual_from_report(report)
    rows, residual = rescore_rows(old_rows)
    base, old_target, old_weight, _ = calibration._contract()
    _, target_pin, weight_pin, _ = contract()
    _write_rows(out / 'target_fit_new_contract.csv', rows)
    receipt = dict(schema='estate_a_target_and_observer_identity_v1',
        observer_contract=OBSERVER_CONTRACT, target_fingerprint=target_pin, weight_fingerprint=weight_pin,
        native_original_target_fingerprint=old_target, native_original_weight_fingerprint=old_weight,
        native_original_loss=float(old_residual@old_residual), loss=float(residual@residual),
        original_wealth_target=OLD_WEALTH_TARGET, wealth_target=NEW_WEALTH_TARGET,
        all_other_targets_and_weights_unchanged=True,
        native_original_target_fit_preserved=True, empirical_bequest_target_unchanged=True,
        estate_flow='positive post-saving bp+(1-psi)*price*chosen_house; death weighted and annualized once',
        aggregate_wealth_stock='beginning living household b+price*house; gross housing stock unchanged',
        receiver_mapping='none', extra_interest_on_bp=False)
    (out / 'estate_a_rescore_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    return rows, residual, parameters, receipt


def make_evaluator(out, lane, P, grid, deadline, price_start=None, *, birth_cap,
                   native_runner=None, exploratory=False, solver=None,
                   target_fingerprint=None, weight_fingerprint=None):
    """Use the copied direct population-one GE, then score its saved moments."""
    _, target_pin, weight_pin, bounds = contract()
    if target_fingerprint is not None and target_fingerprint != target_pin: raise RuntimeError('Estate-A target fingerprint drift')
    if weight_fingerprint is not None and weight_fingerprint != weight_pin: raise RuntimeError('Estate-A weight fingerprint drift')
    if P is None:
        if grid is not None: raise ValueError('grid requires P')
        from .inputs import load_inputs
        P, grid = load_inputs()
    P = copy.deepcopy(P); apply_experiment_flags(P, experiment_flags(birth_cap))
    native = calibration.make_evaluator(out,lane,P,grid,deadline,price_start,
        native_runner=native_runner,exploratory=exploratory,solver=solver,bounds_override=bounds)
    def evaluate(label, point, end):
        result = native(label, point, end)
        if result.get('status') != 'passed': return result
        rows, residual, parameters, receipt = rescore_report(result['report'])
        # Search restrictions are experiment-owned; retain native parameter rows separately.
        parameters = copy.deepcopy(parameters)
        for row in parameters:
            if row['parameter'] == 'beta_annual':
                row['lower'], row['upper'] = '.93', '.99'
                row['near_bound'] = str(min(float(row['estimate'])-.93,.99-float(row['estimate']))<=.0006)
        _write_rows(Path(result['report'])/'parameters_estate_a.csv', parameters)
        result.update(native_original_target_fit=result['target_fit'],native_original_loss=result['loss'],
            target_fit=rows,parameter_table=parameters,residual=residual.tolist(),loss=float(residual@residual),
            objective=float(residual@residual),target_fingerprint=target_pin,weight_fingerprint=weight_pin,
            observer_contract=OBSERVER_CONTRACT,experiment_flags=experiment_flags(birth_cap),
            search_bound_changes={'beta_annual':[.93,.99]})
        return result
    return evaluate
