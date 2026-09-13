"""Freeze only the population-free backward-boundary adapter; no model changes.

Run on Torch after copying the five revised files into the new helper directory.
Submission requires a separate passing native comparison receipt.
"""
import copy
import hashlib
import json
import math
from pathlib import Path
import sys
import time

B = Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913')
old = B/'history_source_extended_seed_v1'
new = B/'history_source_policy_seed_v1'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
read = lambda p: json.loads(Path(p).read_text())
def save(p, value):
    if Path(p).exists():
        assert read(p) == value, 'Refusing to replace a different frozen artifact'
        return
    with Path(p).open('x') as f:
        json.dump(value, f, indent=2); f.write('\n')

changed = sorted(p.name for p in new.glob('*.py') if sha(p) != sha(old/p.name))
assert changed == sorted(['e5f_closed_finite_boundary.py',
    'run_e5f_final_rebated_history.py', 'test_e5f_closed_finite_boundary.py',
    'test_run_e5f_final_rebated_history.py', 'test_e5f_rebated_surprises.py'])
base_path = B/'history_manifest_corrected_extended_seed_v1.json'
base = read(base_path); manifest = copy.deepcopy(base)
pins = {str(new/Path(p).name) if Path(p).parent == old else p: h
        for p, h in base['file_sha256'].items()}
pins.update({str(p): sha(p) for p in new.glob('*.py')})
manifest['file_sha256'] = dict(sorted(pins.items()))
same = copy.deepcopy(manifest); same['file_sha256'] = base['file_sha256']
assert same == base
destination = B/'history_manifest_corrected_policy_seed_v1.json'
save(destination, manifest)
sys.path.insert(0, str(new))
from run_e5f_final_rebated_history import verify_pins
verify_pins(pins)

remaining = math.floor(1789322400 - time.time())
assert 10860 < remaining <= 43200
parent = B/'history_array_manifest_corrected_extended_100.json'
array = read(parent)
runroot = B/'histories_corrected_policy_seed_100'
array['runroot'] = str(runroot); array['horizon_hours'] = remaining/3600
array['environment']['PYTHONPATH'] = array['environment']['PYTHONPATH'].replace(str(old), str(new))
array_pins = {p: h for p, h in array['source_pins'].items()
              if Path(p).parent != old and p != str(base_path)}
array_pins.update(pins); array_pins[str(destination)] = sha(destination)
array['source_pins'] = dict(sorted(array_pins.items()))
for stage in array['stages']:
    argv = stage['command']; stage['seconds'] = remaining
    for flag, value in [('--manifest', str(destination)),
            ('--output', str(runroot/stage['name'])), ('--seconds', str(remaining-60))]:
        assert argv.count(flag) == 1
        argv[argv.index(flag)+1] = value
array_path = B/'history_array_manifest_corrected_policy_seed_100.json'
save(array_path, array)
receipt = dict(status='prepared_not_submitted_native_validation_required',
    old_helper=str(old), new_helper=str(new), changed_files=changed,
    source_sha256={str(new/name):sha(new/name) for name in changed},
    manifest=str(destination), manifest_sha256=sha(destination),
    array_manifest=str(array_path), array_sha256=sha(array_path),
    model_kernels_changed=False, target_parameters_or_gates_changed=False,
    boundary_change='Construct lifetime policy without substituting initial population; audit every actual dated and carried endpoint population.',
    sizing='Two independent100-period histories; four windows, maximum6 preference trials/window and24 mappings/trial plus retained bounded alternative. Observed nonconstant100-period mapping about30minutes. Full completion by18UTC uncertain.',
    safety='Native valid-point and exact6-period-loop equivalence required before launch. Same checkpoints, per-case summaries, five-minute heartbeat,30-minute health threshold,18UTC deadline and2hour policy reserve.')
save(B/'policy_seed_preparation.json', receipt)
print(json.dumps(receipt))
