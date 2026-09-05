#!/usr/bin/env python3
"""Save diagnostic copies for the three tax equilibria and eight Shapley cells."""
import argparse
import copy
import gzip
import math
import pickle
from pathlib import Path

import run_e5f_candidate_policy_comparison as comparison
import run_e5f_independent_numerical_audit as audit
import run_e5f_post2023_rebated_property_tax_smoke as tax
import run_e5f_post2023_coven_property_tax_smoke as fixed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--plan-sha256', required=True)
    args = parser.parse_args()
    plan = comparison.verify_plan(args.plan, args.plan_sha256)
    root = Path(plan['output_root']) / 'rebate_standard_diagnostics'
    if root.exists():
        raise FileExistsError(root)
    root.mkdir()
    active = False
    cell_number = 0
    original_case = tax.run_impact_case
    original_advance = tax.baseline.advance_from_evaluation
    original_joint = tax.solve_joint_rebated_period
    original_fixed = fixed.evaluate_fixed_price

    def save(name, packet, scope):
        out = root / name
        out.mkdir()
        packet = copy.deepcopy(packet)
        with gzip.open(out / 'dated_state.pkl.gz', 'wb', compresslevel=1) as stream:
            pickle.dump(packet, stream, protocol=5)
        audit.standard_diagnostics(packet, out, validate_production_young=False)
        budget = audit.budget_audit(packet, out)
        arrays = audit.policy_array_audit(packet, out)
        graphs = sorted((out / 'standard_diagnostics').glob('*.png'))
        if len(graphs) != 17:
            raise RuntimeError('Incomplete tax diagnostic graph set')
        comparison.write(out / 'receipt.json', dict(status='complete', scope=scope,
            calendar_year=2023, budget=budget, policy_arrays=arrays,
            hashes={str(p.relative_to(out)):comparison.sha(p) for p in [out/'dated_state.pkl.gz', *graphs]}))

    def observed_joint(**kwargs):
        solved = original_joint(**kwargs)
        actual = float(kwargs['P'].property_tax_lump_sum_transfer)
        expected = float(solved.transfer)
        record = dict(case=kwargs['case'].name, selected_transfer=expected,
            parameter_transfer=actual, gap=actual-expected)
        comparison.write(root / f"{kwargs['case'].name}_selected_fiscal_state.json", record)
        if not math.isclose(actual, expected, rel_tol=0, abs_tol=1e-12):
            raise RuntimeError(f'Selected fiscal state differs from parameter state: {record}')
        return solved

    def observed_case(**kwargs):
        nonlocal active
        active = True
        try:
            return original_case(**kwargs)
        finally:
            active = False

    def observed_advance(**kwargs):
        result = original_advance(**kwargs)
        if active:
            save(kwargs['label'], dict(parameters=kwargs['P'], b_grid=kwargs['b_grid'],
                evaluation=kwargs['evaluation'], shared=kwargs['shared'],
                supply_rule=kwargs['supply_rule'], state=kwargs['state']),
                '2023 equilibrium at the declared tax/rebate convention')
        return result

    def observed_fixed(**kwargs):
        nonlocal cell_number
        solved = original_fixed(**kwargs)
        save(f'component_{cell_number:03d}', dict(parameters=kwargs['P'], b_grid=kwargs['b_grid'],
            evaluation=solved.evaluation, shared=solved.shared, state=kwargs['state']),
            'Fixed tax/price/rebate component cell; market clearing is not imposed')
        cell_number += 1
        return solved

    tax.run_impact_case = observed_case
    tax.solve_joint_rebated_period = observed_joint
    tax.baseline.advance_from_evaluation = observed_advance
    fixed.evaluate_fixed_price = observed_fixed
    try:
        args.stage, args.case = 'rebate', 'baseline'
        comparison.run(args)
        if cell_number != 8:
            raise RuntimeError('The eight-cell diagnostic loop is incomplete')
    finally:
        tax.run_impact_case = original_case
        tax.solve_joint_rebated_period = original_joint
        tax.baseline.advance_from_evaluation = original_advance
        fixed.evaluate_fixed_price = original_fixed


if __name__ == '__main__':
    main()
