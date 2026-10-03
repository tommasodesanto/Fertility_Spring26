"""Explicit zero-lifecycle authenticated export; never imported by the loader."""
from pathlib import Path
import hashlib, importlib.util, json, sys, tempfile
import numpy as np
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'code/model/tools'))
import model_playground
model_playground._install_read_only_overlay()
spec = importlib.util.spec_from_file_location('export_timing_calibrate', ROOT / 'code/model/experiments/purchase_timing_sandbox/calibrate.py')
calibrate = importlib.util.module_from_spec(spec); spec.loader.exec_module(calibrate)
v2, timing, manifest, _ = calibrate.checked_inputs('alternative')
source = ROOT / 'output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/native_postcheck/completed.json'
selected = json.loads(source.read_text())
order = model_playground.PARAMETER_ORDER
records = {row['parameter']: float(row['estimate']) for row in selected['parameters']}
point = {key: records[key] for key in order}
_, bounds, _ = v2.inputs.seed_and_bounds('floor_s0')
bounds = {key: tuple(value) for key, value in bounds.items()}; bounds.update(h_P=(.1, 2.6), psi_child=tuple(v2.CONFIG['psi_bounds']))
v2.inputs.LANES['floor_s0'].update(seed=point, bounds=bounds, free_coordinates=list(point))
P, grid = v2.inputs.proposal('floor_s0'); P, entry = v2.inputs.entry(P, grid, 'nonnegative_mean')
with tempfile.TemporaryDirectory() as directory:
    calibrate.install_timing_observer(v2, timing, Path(directory), P)
    P = v2.native.utility_checks(P, grid, 'floor_s0', Path(directory))
from small_credit_lab import credit
credit.bind_engine_credit(P, 'corrected', 0.)
P.H0 = np.array([records['H0']])
# The historical constructor's evidence folder is ephemeral, not an economic input.
P.native_inherited_distribution_evidence_dir = ''
from refactor_lab.inputs import encode
arrays = {'b_grid': grid}
fields = {key: encode(value, key, arrays) for key, value in vars(P).items() if not key.startswith('_')}
folder = Path(__file__).resolve().parent
np.savez_compressed(folder/'arrays.npz', **arrays)
receipt = dict(parameters=fields, source=str(source.relative_to(ROOT)), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    selected_parameters=point, H0=records['H0'], price=0.7760569760205563,
    lifecycle_solves=0, entry=entry, target_fingerprint=manifest['target_fingerprint'], weight_fingerprint=manifest['weight_fingerprint'])
receipt['arrays_sha256'] = hashlib.sha256((folder/'arrays.npz').read_bytes()).hexdigest()
(folder/'bundle.json').write_text(json.dumps(receipt, indent=2, sort_keys=True)+'\n')
print(json.dumps(dict(fields=len(fields), dimensions=[P.Nb,P.Nz], mean_entry_wealth=entry['mean_wealth'], negative_share=entry['negative_wealth_share'], lifecycle_solves=0)))
