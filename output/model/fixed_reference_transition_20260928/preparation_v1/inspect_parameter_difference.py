"""Zero-solve diagnosis of the exact serialized-parameter gate, Torch only."""
import gzip
import importlib.util
import os
import pickle
import sys
from pathlib import Path
import run_preparation as prep

prep.require(sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Torch only')
out = prep.HERE / ('parameter_diagnosis_' + os.environ['SLURM_JOB_ID'])
out.mkdir(exist_ok=False)
prep.require(prep.sha(prep.MANIFEST) == prep.MANIFEST_SHA, 'Manifest pin')
m = prep.read(prep.MANIFEST)
export = Path(m['local_export'])
prep.require(prep.sha(export/'initial_state.pkl.gz') == m['checkpoint']['sha256'], 'Checkpoint pin')
c = prep.read(m['contract']['path'])
sys.path[:0] = [str(prep.ROOT/'code/model/tools'), str(prep.ROOT/'tmp/e5f_overnight_local_20260927/portable/tools_v4')]
import e5f_evening_calibration_runtime as runtime
ev = runtime.setup(dict(c, objective=m['objective']), prep.read(m['objective']['path']), out/'runtime')
with gzip.open(export/'initial_state.pkl.gz', 'rb') as stream:
    packet = pickle.load(stream)
spec = importlib.util.spec_from_file_location('auth_helpers', prep.REFERENCE/'fixed_reference_authenticate.py')
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
raw = helper.json_value(vars(packet['parameters']))
expected = m['actual_serialized_parameters']
diff = {key: dict(saved=raw.get(key), manifest=expected.get(key)) for key in sorted(set(raw)|set(expected)) if raw.get(key) != expected.get(key)}
before = set(vars(packet['parameters']))
ev.tax.actual_parameters(packet['parameters'])
after = helper.json_value(vars(packet['parameters']))
prep.write(out/'result.json', dict(status='diagnostic_only_zero_solves', raw_fields=len(raw), manifest_fields=len(expected),
    differences=diff, after_actual_parameters_matches=(after == expected), added_fields=sorted(set(after)-before),
    reference_checkpoint_sha256=m['checkpoint']['sha256'], manifest_sha256=prep.MANIFEST_SHA))
print('DIAGNOSIS WRITTEN', out/'result.json', flush=True)
