"""Torch-only saved-object inspection. No household or equilibrium solve."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

import numpy as np

BASE = Path('/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a')
BUNDLE = BASE / 'utility_four_arm_preparation_20260925_v2'
CONTRACT_HASH = 'c3fa4d5b7925a54a747c72511538b8182029e1493e4d02e5ef703a30fcecc28a'
CHECKPOINT_HASH = 'c2863ae2fe153a92043df45657d5dc7982bb69753afda0557efa85d1c8871f18'
CASE = BUNDLE / 'results/run_001/floor_linear/floor_linear_initial_0001/case'
os.environ['EXPECTED_UTILITY_COMPARISON_CONTRACT_SHA256'] = CONTRACT_HASH
sys.path.insert(0, str(BUNDLE / 'tools'))
import collect_e5f_utility_comparison as collector


def load():
    contract, runner = collector._load_contract(BUNDLE / 'launch_v1/contract.json')
    assert collector.sha256(CASE / 'initial_state.pkl.gz') == CHECKPOINT_HASH
    scratch = Path(tempfile.mkdtemp(prefix='incentives_runtime_', dir=Path.cwd()))
    _, _, _, tax, _, _, runtime, _ = runner.setup(contract, 'floor_linear', scratch / 'setup')
    saved = collector.scientific_checkpoint(CASE, contract, 'floor_linear', runtime, tax)
    return contract, runtime, saved


def describe(obj):
    result = {}
    for name, value in vars(obj).items():
        if isinstance(value, np.ndarray):
            result[name] = dict(shape=list(value.shape), dtype=str(value.dtype))
        elif isinstance(value, (str, int, float, bool)) or value is None:
            result[name] = value
        else:
            result[name] = dict(type=type(value).__name__)
    return result


if __name__ == '__main__':
    contract, runtime, saved = load()
    packet = saved['packet']
    P, solution, evaluation = packet['parameters'], packet['solution'], packet['evaluation']
    result = dict(native_solves=0, checkpoint_sha256=CHECKPOINT_HASH,
                  packet_fields=list(packet), parameters=describe(P),
                  solution=describe(solution), evaluation=describe(evaluation),
                  policy=describe(evaluation.policy),
                  fecundity=runtime['model'].get_fecundity_by_age(P).tolist())
    for name in ('g_pre', 'g_current'):
        array = getattr(evaluation, name, None)
        if array is not None:
            result[name + '_summary'] = dict(sum=float(np.sum(array)), minimum=float(np.min(array)), finite=bool(np.isfinite(array).all()))
    result['source_functions'] = {}
    import inspect
    for name in ('readiness_gate_active', 'readiness_settled_state', 'independent_child_maturation_active', 'get_fecundity_by_age'):
        fun = getattr(runtime['model'], name)
        result['source_functions'][name] = dict(file=inspect.getsourcefile(fun), source=inspect.getsource(fun))
    Path('inspection.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(dict(status='saved_objects_inspected', native_solves=0)))
