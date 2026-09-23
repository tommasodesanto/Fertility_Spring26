"""Isolated matched preference adapter; frozen solver and numeric gates stay pinned.

prepare_plan generates auditable contract/schema extensions before any run. The
alternative changes housing needs to child-dependent spending shares, retaining
current power equivalence scales. No solve occurs on import or preparation.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import run_e5f_earnings_wealth_candidate as earnings

EARNINGS_PATH = Path(earnings.__file__).resolve()
MODES = ('floor_control', 'child_dependent_shares')
DELTA_DOMAIN = (('delta_alpha_jump', 0., .25, 'softzero'), ('delta_alpha', 0., .25, 'softzero'))


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def replace_once(text, before, after):
    if text.count(before) != 1:
        raise ValueError('Source replacement is not unique: ' + before[:100])
    return text.replace(before, after, 1)


def preference_mode(plan):
    spec = plan['preference_specification']
    if spec.get('author_decision') != 'approved_diagnostic' or spec.get('mapping') not in MODES:
        raise ValueError('Explicit approved diagnostic preference mapping required')
    if plan.get('entry_specification', {}).get('rule') != 'zero_assets':
        raise ValueError('This matched experiment authorizes only explicit zero-asset entry')
    if spec['mapping'] == 'child_dependent_shares' and not plan.get('wealth_grid_specification'):
        raise ValueError('Share probe requires the declared V5 extended-grid generator')
    return spec['mapping']


def share_domain(parent):
    return tuple(row for row in parent.PARENTHOOD_SEARCH_DOMAIN if row[0] != 'h_P') + DELTA_DOMAIN


def install_share_utility(parent):
    """Replace only utility migration/binding/schema before native probe import."""
    domain = share_domain(parent)
    names = tuple(row[0] for row in domain)
    fixed = dict(parent._FIXED_UTILITY, child_room_floor=False)
    fixed.pop('delta_alpha'); fixed.pop('delta_alpha_jump')
    fixed.update(hbar_first_child_jump=0., hbar_child_rooms=0.)
    scalar = parent._finite_scalar

    def validate(P):
        parent._validate_lifecycle(P)
        for name, expected in fixed.items():
            value = getattr(P, name, None)
            valid = value == expected if isinstance(expected, str) else (value is expected if isinstance(expected, bool) else scalar(value, name) == expected)
            if not valid:
                raise ValueError('Share utility requires ' + name + '=' + repr(expected))
        for name, low, high, _ in DELTA_DOMAIN:
            if not low <= scalar(getattr(P, name), name) <= high:
                raise ValueError('Share coordinate outside bounds: ' + name)
        if not 0. < scalar(P.beta, 'beta') < 1.:
            raise ValueError('Period beta outside (0,1)')
        scalar(P.psi_child, 'psi_child')

    def initialize(P):
        parent._validate_lifecycle(P)
        result = copy.deepcopy(P)
        for name, value in fixed.items(): setattr(result, name, value)
        for name, _, _, _ in DELTA_DOMAIN: setattr(result, name, 0.)
        validate(result)
        return result

    def candidate(values=None, *, require_complete=False):
        values = {} if values is None else dict(values)
        if set(values) - set(names) or (require_complete and set(values) != set(names)):
            raise ValueError('Share candidate must use exactly the declared ten structural coordinates')
        out = {}
        for name, low, high, _ in domain:
            if name in values:
                value = scalar(values[name], name)
                if not low <= value <= high: raise ValueError('Coordinate outside bounds: ' + name)
                out[name] = value
        return out

    def bind(P, structural=None, *, copy_parameters=True):
        validate(P)
        checked = candidate(structural)
        result = copy.deepcopy(P)
        for name, value in checked.items():
            if name == 'beta_annual':
                result.beta = value ** 4
                result.rho = 1. / result.beta - 1.
                result.rho_hat = result.rho
            elif name == 'H0':
                if np.asarray(P.H0).size != 1: raise ValueError('One-market scalar H0 required')
                result.H0 = np.full(np.asarray(P.H0).shape, value) if np.asarray(P.H0).shape else value
            else: setattr(result, name, value)
        if 'kappa_fert' in checked: result.eps_fert = result.kappa_fert
        validate(result)
        if copy_parameters: return result
        for name, value in vars(result).items(): setattr(P, name, value)
        return P

    parent.initialize_parenthood_utility = initialize
    parent.validate_parenthood_utility = validate
    parent.validate_parenthood_candidate = candidate
    parent.bind_parenthood_utility = bind
    parent.PARENTHOOD_SEARCH_DOMAIN = domain
    parent.PARENTHOOD_SEARCH_NAMES = names
    parent.UTILITY_CONTRACT = 'e5f_power_scale_child_spending_shares_diagnostic_v1'
    parent.parenthood_utility_metadata = lambda: dict(contract=parent.UTILITY_CONTRACT,
        search_domain=domain, free_parameter_count=10, fixed_utility=fixed,
        alpha_formula='alpha(0)=0.733; alpha(m)=clip(0.733-delta_alpha_jump-delta_alpha*m,0.05,0.95) for m>0',
        scale_formula='((2+0.7*m)/2)**0.7', adoption='experimental only')
    return domain


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': earnings.digest(path)}


def prepare_plan(plan, directory):
    """Return a pinned plan; generated files live outside the immutable source."""
    plan = copy.deepcopy(plan)
    mode = preference_mode(plan)
    directory = Path(directory).resolve(); directory.mkdir(parents=True, exist_ok=False)
    for name,item in plan['files'].items():
        if earnings.digest(item['path']) != item['sha256']:
            raise ValueError('Changed parent plan file: '+name)
    old_contract = earnings.read(plan['files']['run_contract']['path'])
    for name in ('working_objective','scorer','validator'):
        item=old_contract[name]
        if earnings.digest(item['path']) != item['sha256']:
            raise ValueError('Changed parent scoring file: '+name)
    objective = earnings.read(old_contract['working_objective']['path'])
    if fingerprint(objective) != earnings.OBJECTIVE:
        raise ValueError('Preparation requires the untouched V5 objective')
    target_payload = {k:v for k,v in objective.items() if k != 'parameter_restrictions'}
    plan['target_system_sha256'] = fingerprint(target_payload)
    plan['parent_objective_canonical_sha256'] = earnings.OBJECTIVE
    plan['files']['earnings_adapter'] = pin(EARNINGS_PATH)
    plan['files']['adapter'] = pin(__file__)
    plan['adapter_path'] = str(Path(__file__).resolve())
    plan['adapter_sha256'] = earnings.digest(__file__)
    if mode == 'child_dependent_shares':
        restrictions = [row for row in objective['parameter_restrictions'] if row['parameter'] != 'h_P']
        if len(restrictions) != 8: raise ValueError('Unexpected parent restriction set')
        for name, low, high, transform in DELTA_DOMAIN:
            restrictions.append(dict(parameter=name, lower=low, upper=high, transform=transform))
        objective['parameter_restrictions'] = restrictions
        new_hash = fingerprint(objective)
        objective_path = directory / 'objective.json'; earnings.write(objective_path, objective)
        scorer_path = directory / 'scorer.generated.py'
        scorer = Path(old_contract['scorer']['path']).read_text()
        scorer = replace_once(scorer, "'first_birth_fixed_cost', 'h_P')", "'first_birth_fixed_cost', 'delta_alpha_jump', 'delta_alpha')")
        scorer = replace_once(scorer, "'free_parameter_count': 9", "'free_parameter_count': len(PARAMETERS)")
        scorer_path.write_text(scorer)
        validator_path = directory / 'validator.generated.py'
        validator = Path(old_contract['validator']['path']).read_text()
        validator = replace_once(validator, "{'hbar_child_rooms':0.,'payroll_tax'", "{'hbar_first_child_jump':0.,'hbar_child_rooms':0.,'payroll_tax'")
        validator += "\n# Explicit ten-coordinate alternative; all numeric gate functions unchanged.\nNAMES=tuple(n for n in NAMES if n!='h_P')+('delta_alpha_jump','delta_alpha')\nBOUNDS={k:v for k,v in BOUNDS.items() if k!='h_P'}\nBOUNDS.update(delta_alpha_jump=(0.,.25),delta_alpha=(0.,.25))\n"
        validator_path.write_text(validator)
        wrapper = Path(plan['files']['wrapper']['path']).read_text()
        wrapper = replace_once(wrapper, "APPROVED_OBJECTIVE = '" + earnings.OBJECTIVE + "'", "APPROVED_OBJECTIVE = '" + new_hash + "'")
        wrapper = replace_once(wrapper, "len(parameters)==17", "len(parameters)==19")
        wrapper = replace_once(wrapper, "len(scored['parameters'])==17", "len(scored['parameters'])==19")
        wrapper_path = directory / 'wrapper.generated.py'; wrapper_path.write_text(wrapper)
        plan['files']['wrapper'] = pin(wrapper_path)
        old_contract['wrapper_sha256'] = earnings.digest(wrapper_path)
        old_contract['working_objective'] = dict(pin(objective_path), canonical_sha256=new_hash)
        old_contract['scorer'] = pin(scorer_path); old_contract['validator'] = pin(validator_path)
        for name,path in [('preference_objective',objective_path),('preference_scorer',scorer_path),('preference_validator',validator_path)]:
            plan['files'][name] = pin(path)
        plan['objective_canonical_sha256'] = new_hash
    else:
        plan['objective_canonical_sha256'] = earnings.OBJECTIVE
    contract_path = directory / 'run_contract.json'; earnings.write(contract_path, old_contract)
    plan['files']['run_contract'] = pin(contract_path)
    earnings.write(directory / 'preference_contract.json', dict(mapping=mode,
        parent_objective_sha256=earnings.OBJECTIVE, objective_sha256=plan['objective_canonical_sha256'],
        target_system_sha256=plan['target_system_sha256'], source_manifest_unchanged=True,
        runtime_extension=plan['files']['adapter'], generated_files={k:v for k,v in plan['files'].items() if k.startswith('preference_') or k=='wrapper'},
        parameter_rows=19 if mode=='child_dependent_shares' else 17,
        preference_changes=['remove housing floor','child-dependent housing expenditure shares'] if mode=='child_dependent_shares' else [],
        unchanged=['power equivalence scale','earnings','zero entry assets','current-income purchase eligibility','target rows and weights','numerical gates']))
    plan['files']['preference_contract'] = pin(directory / 'preference_contract.json')
    return plan


def verify_plan(plan):
    preference_mode(plan)
    check = copy.deepcopy(plan)
    check['files']['adapter'] = plan['files']['earnings_adapter']
    check['objective_canonical_sha256'] = earnings.OBJECTIVE
    earnings.verify_plan(check)
    if Path(plan['files']['adapter']['path']).resolve() != Path(__file__).resolve():
        raise ValueError('Wrong preference adapter path')
    if earnings.digest(__file__) != plan['files']['adapter']['sha256']:
        raise ValueError('Preference adapter hash mismatch')
    run = earnings.read(plan['files']['run_contract']['path'])
    objective = earnings.read(run['working_objective']['path'])
    if fingerprint(objective) != plan['objective_canonical_sha256']:
        raise ValueError('Derived objective mismatch')
    if fingerprint({k:v for k,v in objective.items() if k!='parameter_restrictions'}) != plan['target_system_sha256']:
        raise ValueError('Target/weight/measurement contract changed')


def run_probe(plan, initial_path, output):
    source = Path(plan['source_root'])
    sys.path[:0] = [str(source/'code/model/tools'), str(source/'code/model')]
    if preference_mode(plan) == 'child_dependent_shares':
        import e5f_parenthood_utility as parent
        install_share_utility(parent)
        original = earnings.accounting.rewrite_probe_grid
        def rewrite(text, nodes):
            text = original(text, nodes)
            return replace_once(text, "('hbar_child_rooms',P.hbar_child_rooms,", "('hbar_first_child_jump',P.hbar_first_child_jump,'zero restriction'),\n                ('hbar_child_rooms',P.hbar_child_rooms,")
        earnings.accounting.rewrite_probe_grid = rewrite
    earnings.run_probe(plan, initial_path, output)


def run_case(plan_path, plan, arm, output, repetitions, preflight=False):
    # The parent's subprocess command must re-enter this explicit adapter.
    previous = earnings.__file__
    earnings.__file__ = __file__
    try:
        result = earnings.run_case(plan_path, plan, arm, output, repetitions, preflight)
    finally:
        earnings.__file__ = previous
    if not preflight:
        path = Path(output)/'runtime_contract.json'; receipt = earnings.read(path)
        receipt['preference_specification'] = plan['preference_specification']
        receipt['target_system_sha256'] = plan['target_system_sha256']
        receipt['additional_runtime_files'].update({k:v for k,v in plan['files'].items() if k.startswith('preference_') or k=='earnings_adapter'})
        if preference_mode(plan)=='child_dependent_shares':
            receipt['changed_economic_objects'] += ['no housing needs floor','child-dependent expenditure shares']
        earnings.write(path, receipt)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True); parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mode', choices=('pilot','probe-child'), default='pilot')
    parser.add_argument('--arm', choices=earnings.ARMS); parser.add_argument('--repetitions',type=int)
    parser.add_argument('--contract',type=Path); parser.add_argument('--preflight',action='store_true')
    args = parser.parse_args(); plan = earnings.read(args.plan); verify_plan(plan)
    if args.mode=='probe-child': run_probe(plan,args.contract,args.output)
    elif args.preflight:
        results=[run_case(args.plan,plan,c['arm'],args.output/c['id'],c['repetitions'],True) for c in plan['cases']]
        args.output.mkdir(exist_ok=True,parents=True); earnings.write(args.output/'receipt.json',{'cases':results,'solves':0})
    else: print(json.dumps(run_case(args.plan,plan,args.arm,args.output,args.repetitions)))

if __name__=='__main__': main()
