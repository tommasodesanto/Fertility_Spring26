"""One six-path diagnostic using the active exhaustive-saving fixture.

Preserve legacy failure 17357939. Run setup exactly once (24 dated Bellman
calls); run individual test instances afterward without unittest class setup.
This records failure evidence and candidate states. It does not alter tests or
select a replacement acceptance probe. Invoke only on the cluster.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import traceback
import unittest
from unittest.mock import patch

import numpy as np


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def native(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def value_and_control_row(module, control_probe, changed_probe, state):
    before_inputs, after_inputs = control_probe['inputs'], changed_probe['inputs']
    wealth, tenure, _, _, _, parity, children = state
    if children != 0:
        raise ValueError('Diagnostic scalar observer is declared for zero dependents')
    before = float(control_probe['output'][1][wealth, parity])
    after = float(changed_probe['output'][1][wealth, parity])
    evaluate = module.conditional_probe_objective
    b = evaluate(before_inputs, state, before)
    a = evaluate(after_inputs, state, after)
    old_under_new = evaluate(after_inputs, state, before)
    new_under_old = evaluate(before_inputs, state, after)
    key = 'Vc_flat' if tenure == 0 else 'Vco_flat'
    current_inputs_equal = all(np.array_equal(before_inputs[k], after_inputs[k])
                               for k in before_inputs if k != key)
    before_interior = b['lower'] + 1e-6 < before < min(b['upper'], before_inputs['b_grid'][-1])-1e-6
    after_interior = a['lower'] + 1e-6 < after < min(a['upper'], after_inputs['b_grid'][-1])-1e-6
    before_gain = b['value'] - new_under_old['value']
    after_gain = a['value'] - old_under_new['value']
    return dict(
        state=list(state), saving_control=before, saving_changed=after,
        saving_change=after-before, control_objective=b, changed_objective=a,
        old_control_under_new_continuation=old_under_new,
        new_control_under_old_continuation=new_under_old,
        control_revealed_preference_gain=before_gain,
        changed_revealed_preference_gain=after_gain,
        control_kernel_value=float(control_probe['output'][0][wealth, parity]),
        changed_kernel_value=float(changed_probe['output'][0][wealth, parity]),
        control_kernel_value_error=float(control_probe['output'][0][wealth, parity])-b['value'],
        changed_kernel_value_error=float(changed_probe['output'][0][wealth, parity])-a['value'],
        control_interior=bool(before_interior), changed_interior=bool(after_interior),
        current_inputs_equal=bool(current_inputs_equal),
        continuation_max_change=float(np.max(np.abs(after_inputs[key][:,parity]-before_inputs[key][:,parity]))),
        control_on_wealth_knot=bool(np.min(np.abs(before_inputs['b_grid']-before)) < 1e-10),
        changed_on_wealth_knot=bool(np.min(np.abs(after_inputs['b_grid']-after)) < 1e-10),
        informative_both_interior=bool(before_interior and after_interior
            and abs(after-before)>1e-8 and before_gain>1e-10 and after_gain>1e-10),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, default=Path.cwd())
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    sys.path[:0] = [str(root/'code/model/tools'), str(root/'code/model')]
    import test_e5f_social_security_compiled as module

    T = module.CompiledSocialSecurityTests
    started = time.monotonic()
    packet = dict(status='started', original_failed_job='17357939',
        scope='same six two-date conditional paths, active saving method only',
        expected_cases=6, expected_dated_bellman_calls=24,
        source_sha256={str(path.relative_to(root)): file_hash(path) for path in (
            root/'code/model/tools/test_e5f_social_security_compiled.py',
            root/'code/model/tools/test_run_e5f_perfect_foresight_transition.py',
            root/'code/model/intergen_eqscale_seq_optimized/parameters.py',
            root/'code/model/intergen_eqscale_seq_optimized/solver.py',
            root/'code/model/intergen_eqscale_seq_optimized/kernels.py')},
        diagnostic_sha256=file_hash(Path(__file__)),
        selected_new_probe=None,
        interpretation='Diagnostic only; failure lists are not waived acceptance tests')
    def save():
        packet['elapsed_seconds'] = time.monotonic()-started
        (out/'active_optimizer_diagnostic.json').write_text(
            json.dumps(packet, indent=2, allow_nan=False, default=native)+'\n')
    save()
    original_overrides = module.fixture._tiny_overrides
    seen_overrides = []
    def active_overrides():
        original = original_overrides()
        active = dict(original, exhaustive_saving_control=True)
        changed_keys = [key for key in set(original)|set(active)
                        if key not in original or key not in active
                        or not np.array_equal(original[key], active[key])]
        if changed_keys != ['exhaustive_saving_control']:
            raise RuntimeError(f'Unexpected fixture override differences: {changed_keys}')
        seen_overrides.append(dict(original=original, active=active,
                                   changed_keys=changed_keys))
        return active
    try:
        print('Starting exactly one active-optimizer six-path setup.', flush=True)
        # Set the flag BEFORE the initial fixture solution and terminal policy
        # are created. No after-the-fact parameter mutation is permitted.
        with patch.object(module.fixture, '_tiny_overrides', side_effect=active_overrides):
            T.setUpClass()
        if len(seen_overrides) != 1 or not T.parameters.exhaustive_saving_control:
            raise RuntimeError('Active optimizer was not bound exactly once at fixture construction')
        packet['fixture_override_receipt'] = seen_overrides[0]
        packet['active_exhaustive_saving_control'] = bool(T.parameters.exhaustive_saving_control)
        if len(T.results) != 6:
            raise RuntimeError('Expected six and only six conditional paths')
        packet['dated_bellman_calls'] = sum(result['path'].bellman_solves for result in T.results.values())
        if packet['dated_bellman_calls'] != 24:
            raise RuntimeError('Expected exactly 24 dated Bellman calls')
        flags = [int(probe['inputs']['exhaustive_saving'])
                 for result in T.results.values() for call in result['bellman_calls']
                 for probe in call['probes'].values()]
        if not flags or set(flags) != {1}:
            raise RuntimeError(f'Actual compiled kernel flags are not all active: {set(flags)}')
        packet['observed_kernel_flag_count'] = len(flags)
        packet['setup_completed'] = True
        save()
        print('Setup complete; checking all seven methods without another setup.', flush=True)
        result = unittest.TestResult()
        for method in unittest.defaultTestLoader.getTestCaseNames(T):
            T(method).run(result)  # TestCase.run calls no setUpClass.
        packet['tests'] = dict(methods_run=result.testsRun,
            failures=[dict(test=str(test), traceback=detail) for test, detail in result.failures],
            errors=[dict(test=str(test), traceback=detail) for test, detail in result.errors],
            skipped=[dict(test=str(test), reason=reason) for test, reason in result.skipped],
            all_pass=result.wasSuccessful())
        (out/'test_failures.txt').write_text('\n\n'.join(
            str(test)+'\n'+detail for test,detail in result.failures+result.errors)+'\n')
        packet['paths'] = {name:dict(
            bellman_calls=entry['path'].bellman_solves,
            maximum_policy_reproduction_error=entry['path'].maximum_policy_reproduction_error,
            maximum_mass_accounting_error=entry['path'].maximum_mass_accounting_error,
            maximum_feasibility_projection_mass=entry['path'].maximum_feasibility_projection_mass,
            budgets=[date['budget'] for date in entry['dated']])
            for name,entry in T.results.items()}
        packet['cases'] = []
        for changed_name, control_name in (('future_pension','fixed_tax'),('future_tax','fixed_pension')):
            changed,control = T.results[changed_name],T.results[control_name]
            sign = 1.0 if changed_name=='future_pension' else -1.0
            case = dict(case=changed_name, control=control_name, value_extrema=[])
            for date in (0,1):
                delta = changed['path'].values[date]-control['path'].values[date]
                signed = sign*delta
                index = tuple(int(x) for x in np.unravel_index(np.argmin(signed),signed.shape))
                case['value_extrema'].append(dict(date=date, minimum_signed_change=float(signed[index]),
                    worst_state=list(index), control_value=float(control['path'].values[date][index]),
                    changed_value=float(changed['path'].values[date][index]),
                    wrong_sign_count=int(np.count_nonzero(signed < -2e-10)),
                    maximum_signed_change=float(np.max(signed)),
                    minimum_by_age=[float(np.min(signed[:,:,:,age])) for age in range(T.parameters.J)]))
            left,right = changed['dated'][0]['evaluation'],control['dated'][0]['evaluation']
            occupied=right.g_current>1e-12
            delta_bp = left.policy.bp_pol-right.policy.bp_pol
            max_index=tuple(int(x) for x in np.unravel_index(np.abs(delta_bp).argmax(),delta_bp.shape))
            case['all_policy_response'] = dict(maximum_abs_saving_change=float(np.max(np.abs(delta_bp))),
                maximum_abs_occupied_saving_change=float(np.max(np.abs(delta_bp[occupied]))),
                maximum_state=list(max_index), control=float(right.policy.bp_pol[max_index]),
                changed=float(left.policy.bp_pol[max_index]), control_mass=float(right.g_current[max_index]))
            probes = [entry['bellman_calls'][1]['probes'][changed_name] for entry in (control,changed)]
            state = module.ANTICIPATION_PROBES[changed_name]
            try:
                case['original_declared_probe'] = value_and_control_row(module,*probes,state)
            except Exception:
                case['original_declared_probe_error'] = traceback.format_exc()
            # Diagnostic scan of already-observed blocks, no new model calls.
            # Keep every zero-dependent state, including failed/uninformative
            # ones, and never promote a new acceptance probe in this script.
            case['same_block_zero_dependent_scan'] = []
            for wealth in range(len(T.b_grid)):
                for parity in range(T.parameters.n_parity):
                    candidate=(wealth,*state[1:5],parity,0)
                    try:
                        row=value_and_control_row(module,*probes,candidate)
                        row['control_mass']=float(right.g_current[candidate])
                    except Exception as exc:
                        row=dict(state=list(candidate), observer_error=str(exc))
                    case['same_block_zero_dependent_scan'].append(row)
            case['informative_both_interior_states']=[row['state'] for row in case['same_block_zero_dependent_scan']
                if row.get('informative_both_interior')]
            packet['cases'].append(case)
            save()
        packet['status'] = 'diagnostic_complete_tests_passed' if result.wasSuccessful() else 'diagnostic_complete_tests_failed'
        save()
        print(json.dumps(dict(status=packet['status'],methods_run=result.testsRun,
            failures=len(result.failures),errors=len(result.errors),dated_bellman_calls=packet['dated_bellman_calls'])),flush=True)
        # A complete diagnostic returns zero even if acceptance tests fail;
        # its explicit status/individual traces distinguish that outcome.
        return 0
    except Exception:
        packet['status']='diagnostic_failed'
        packet['error']=traceback.format_exc()
        save()
        print(packet['error'],file=sys.stderr,flush=True)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
