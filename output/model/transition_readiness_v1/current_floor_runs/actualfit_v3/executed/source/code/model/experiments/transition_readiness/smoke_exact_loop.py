#!/usr/bin/env python3
"""One bounded control-flow smoke of the real run() loop; all native calls stubbed.

This checks orchestration artifacts only and cannot certify economics/numerics.
"""
import gzip
import importlib.util
import json
import os
from pathlib import Path
import pickle
from types import SimpleNamespace as NS
from unittest.mock import patch

import numpy as np
import legacy_changed_psi as driver


def smoke(output):
    spec = importlib.util.spec_from_file_location('smoke_accounting', driver.HERE / 'pinned_tools/run_e5f_preference_budget_diagnostic.py')
    accounting = importlib.util.module_from_spec(spec); spec.loader.exec_module(accounting)
    calls = []; parameters = NS(psi_child=driver.PSI, pension=1.)
    packet = dict(parameters=parameters, solution=NS(p_eq=[1.]))
    state = NS(g_pre=np.zeros(2), scheduled_entries=np.zeros(2), scheduled_raw_entries=np.zeros(2))
    def checkpoint(path, value):
        with gzip.open(path, 'wb') as stream:
            pickle.dump(value, stream)
    def mapping(*args, **kwargs):
        calls.append(str(args[8].name))
        rows = [dict(asset_price=1., payroll_tax_revenue=1., pension_outlays=1.,
            pension_period_units=1., housing_demand=1., housing_supply=1., scaled_pension_budget_residual=0.) for _ in range(6)]
        record = dict(rows=rows, fertility=[dict(period_tfr_topcode_adjusted=2.1) for _ in range(6)],
            market_residual=[0.]*6, fiscal_residual=[0.]*6,
            gates=dict(a=True, b=True, c=True, d=True), cache={}, diagnostic_packets=[])
        return NS(terminal_state=state, rows=rows), record
    inner = NS(write=accounting.write, plain=lambda v:v, dump_checkpoint=checkpoint,
        load_reference=lambda out: (dict(standard_diagnostic_names=[]), packet, NS(rt=dict(audit=None,
            primitive=NS(pf=NS(birth_queue_values=lambda v:v))))),
        shock_path=lambda spec,horizon: spec['levels']*horizon, mapping=mapping,
        terminal_checks=lambda *args: dict(all_checks_pass=True),
        render_diagnostics=lambda *args: dict(stubbed_renderer=True))
    class Native:
        def __init__(self, *args):
            self.deadline=0
        def endpoint(self, psi):
            return packet, dict(price=1., psi_child=psi)
        def guarded(self, seconds, call):
            return call()
    estimator = NS(NativeEstimator=Native, draft_plan=lambda kind: dict(housing='fixed_stock', budget={}, path={}))
    modules = {'run_e5f_preference_transition.py': inner, 'run_e5f_preference_estimation.py': estimator,
               'run_e5f_preference_budget_diagnostic.py': accounting}
    env = {key:'1' for key in ('NUMBA_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')}
    env['SLURM_JOB_ID']='123'
    with patch.object(driver, 'import_pinned', return_value=modules), patch.object(driver.sys, 'platform', 'linux'), patch.dict(os.environ, env):
        driver.run(output)
    complete=json.loads((output/'complete.json').read_text())
    assert calls == ['changed_psi_input','pension_trial','fresh_replay'], calls
    for name in calls:
        assert (output/(name+'_checkpoint.pkl.gz')).is_file()
    for name in ('latest_completed.json','best_so_far.json','heartbeat.json','complete.json','preflight.json'):
        assert (output/name).is_file()
    assert complete['numerical_certified'] is True
    assert complete['joint_price_pension_root_verified'] is False
    assert complete['production_ready'] is False
    proof=dict(status='PASS', native_model_calls=0, real_run_loop_map_invocations=calls,
               checkpoints=3, heartbeat=True, latest_completed=True, best_so_far=True,
               stubbed_endpoint_and_renderer=True, scientific_validation=False)
    accounting.write(output/'control_flow_smoke_receipt.json',proof)
    return proof


if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    print(json.dumps(smoke(parser.parse_args().output),indent=2))
