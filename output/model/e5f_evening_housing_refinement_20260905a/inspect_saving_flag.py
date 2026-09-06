#!/usr/bin/env python3
"""One flagged calibration state, two saving methods, fixed prices and inherited mass.

The first pair is the exact-loop smoke. Stop on reproduction/dominance failure.
No parameter search or new market equilibrium is performed.
"""
import copy
import gzip
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT/'code/model'), str(ROOT/'code/model/tools')]
import run_e5f_candidate_policy_comparison as common
import run_e5f_global_saving_quantification as global_check
import run_e5f_bounded_calibration_refinement as base
import run_e5f_transition_calibration as calibration
from intergen_eqscale_seq_optimized import solver


if __name__ == '__main__':
    manifest_path = Path(__file__).with_name('saving_flag_plan.json')
    if len(sys.argv)!=2 or common.sha(manifest_path)!=sys.argv[1]:
        raise RuntimeError('Pinned value-check plan required')
    manifest = common.read(manifest_path)
    if common.sha(__file__) != manifest['driver_sha256']:
        raise RuntimeError('Diagnostic source changed')
    plan = base.load_plan(Path(manifest['parent_plan']), manifest['parent_plan_sha256'])
    if calibration.code_fingerprint_contract(solver)['bundle_sha256'] != plan['code_bundle_sha256']:
        raise RuntimeError('Scientific bundle changed')
    if len(manifest['cases']) != 1:
        raise RuntimeError('This diagnostic permits only one pair of fixed-price solves')
    for path, digest in manifest['helper_hashes'].items():
        if common.sha(path) != digest:
            raise RuntimeError(f'Diagnostic helper changed: {path}')
    out = Path(__file__).with_name('saving_flag_check')
    if out.exists():
        raise FileExistsError(out)
    out.mkdir()
    audit = global_check.audit
    global_check.transition.configure_sequential_model()
    global_check.install_sequential_operators()
    global_check.operator_contract()
    results = []
    started = time.time()
    try:
        for index, item in enumerate(manifest['cases']):
            path = Path(item['checkpoint'])
            if common.sha(path) != item['sha256']:
                raise RuntimeError('Checkpoint changed')
            case = next(c for c in plan['cases'] if c['id']==item['case_id'])
            receipt = base.read_json(path.parent/'case_receipt.json')
            if receipt['plan_sha256'] != manifest['parent_plan_sha256'] or receipt['artifact_sha256']['dated_state.pkl.gz'] != item['sha256']:
                raise RuntimeError('Checkpoint is not the declared completed case')
            base.validate_result(path.parent, plan, case)
            packet = audit.load_checkpoint(path)
            original = packet['evaluation']
            reference = global_check.quantities(original, packet['parameters'])
            for method in ('local', 'global'):
                common.write(out/'heartbeat.json', dict(case=item['label'], method=method, phase='solving', elapsed_seconds=time.time()-started))
                global_check.set_method(method)
                P = copy.deepcopy(packet['parameters'])
                shared = global_check.model.precompute_shared(P, packet['b_grid'])
                e = global_check.calendar.evaluate_period(original.policy.price.copy(),
                    original.g_pre.copy(), P, packet['b_grid'], shared,
                    global_check.calendar.SolveCounter(), packet['supply_rule'])
                if method == 'local':
                    for field in ('V', 'bp_pol', 'hR_pol', 'fert_probs', 'tenure_probs', 'fert2_probs'):
                        if not np.array_equal(getattr(e.policy, field), getattr(original.policy, field), equal_nan=True):
                            raise RuntimeError(f'Fresh local policy failed exact reproduction: {item["label"]} {field}')
                    for field in ('g_pre', 'g_current', 'g_post_fertility'):
                        if not np.array_equal(getattr(e, field), getattr(original, field)):
                            raise RuntimeError(f'Fresh local distribution failed exact reproduction: {field}')
                occupied = original.g_pre > 1e-12
                change = e.policy.V - original.policy.V
                if method == 'global' and float(change[occupied].min()) < -1e-7:
                    raise RuntimeError('Exhaustive-saving value dominance failed')
                case_out = out / f'{index:02d}_{method}'
                case_out.mkdir()
                case_packet = dict(packet, parameters=P, evaluation=e, shared=shared)
                arrays = audit.policy_array_audit(case_packet, case_out)
                audit.standard_diagnostics(case_packet, case_out, validate_production_young=False)
                quantities = global_check.quantities(e, P)
                record = dict(case=item['label'], method=method, checkpoint_sha256=item['sha256'],
                    quantities=quantities, policy_arrays=arrays,
                    births_difference_percent=100*(quantities['births_per_household']/reference['births_per_household']-1),
                    rooms_difference_percent=100*(quantities['rooms_per_household']/reference['rooms_per_household']-1),
                    ownership_difference_pp=100*(quantities['ownership']-reference['ownership']),
                    minimum_occupied_value_change=float(change[occupied].min()), elapsed_seconds=time.time()-started)
                with gzip.open(case_out/'dated_state.pkl.gz', 'wb', compresslevel=1) as stream:
                    pickle.dump(case_packet, stream, protocol=5)
                common.write(case_out/'receipt.json', record)
                results.append(record)
                common.write(out/'latest_completed_case.json', record)
                common.write(out/'best_so_far.json', {'scope':'diagnostic, no search', 'largest_abs_birth_change_percent':max(abs(r['births_difference_percent']) for r in results)})
            if index == 0:
                common.write(out/'exact_loop_smoke.json', {'status':'pass','local_reproduction':'exact','global_occupied_value_dominance':True})
        common.write(out/'summary.json', dict(status='complete', cases=results,
            scope='Fixed prices, parameters and inherited distributions; no new equilibrium or calibration',
            elapsed_seconds=time.time()-started, model_solves=len(results)))
    except BaseException as error:
        common.write(out/'failure.json', {'error':str(error),'elapsed_seconds':time.time()-started})
        raise
    finally:
        global_check.set_method('local')
