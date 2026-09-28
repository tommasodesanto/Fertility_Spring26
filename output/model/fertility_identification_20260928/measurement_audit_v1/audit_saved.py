"""Read-only Torch inspection of the September 28 primary checkpoint; zero solves.

Run in the pinned September 28 container, writing only to this audit directory.
The runtime authenticates source ancestry before deserializing model classes.
No constructor is used as evidence of the saved parameter configuration.
"""
import os
import sys
from pathlib import Path
import csv
import copy
import gzip
import hashlib
import json
import math
import pickle

assert sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit()
ROOT = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
B = ROOT / 'output/model/fertility_identification_20260928'
OUT = B / 'measurement_audit_v1'
SOURCE = B / 'resume_v1/selected_export/primary'
sys.path.insert(0, str(ROOT / 'code/model/tools'))
import numpy as np
import run_e5f_fertility_identification as driver
import e5f_evening_calibration_runtime as runtime


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


print('Authenticating contract, all export files, and checkpoint', flush=True)
c, objectives = driver.verify(B / 'contract_v1/contract.json')
hashes = read(SOURCE / 'artifact_hashes.json')
for name, digest in hashes.items():
    assert sha(SOURCE / name) == digest, name
receipt = read(SOURCE / 'receipt.json')
checkpoint = SOURCE / 'initial_state.pkl.gz'
assert sha(checkpoint) == receipt['case_checkpoint_sha256']
assert hashes['parameters.csv'] == c['anchor']['parameters']['sha256']
evaluator = runtime.setup(dict(c, objective=c['lanes']['primary']['objective']),
                          objectives['primary'], OUT / ('runtime_' + os.environ['SLURM_JOB_ID']))
with gzip.open(checkpoint, 'rb') as stream:
    packet = pickle.load(stream)
P = packet['parameters']
matrix = np.asarray(P.Pi_z).copy()
z = np.asarray(P.z_grid).copy()
w = np.asarray(P.z_weights).copy()
print('Saved parameter object loaded; independently reconstructing B15', flush=True)
external_path = OUT / 'inputs/single_process_external_estimate.json'
external = read(external_path)
approved = external['point_fit']['pure_ar1_method_of_moments']
rho = approved['rho_four_year']
variance = approved['stationary_log_variance']
n = 15
assert matrix.shape == (n, n) and z.shape == w.shape == (n,)
# Independent Rouwenhorst recursion with p=q=(1+rho)/2.
p = (1 + rho) / 2
expected = np.array([[p, 1-p], [1-p, p]])
for size in range(3, n + 1):
    old = expected
    expected = np.zeros((size, size))
    expected[:-1, :-1] += p * old
    expected[:-1, 1:] += (1-p) * old
    expected[1:, :-1] += (1-p) * old
    expected[1:, 1:] += p * old
    expected[1:-1] /= 2
expected_w = np.array([math.comb(n-1, i) / 2**(n-1) for i in range(n)])
expected_z = np.exp(np.linspace(-math.sqrt((n-1)*variance), math.sqrt((n-1)*variance), n))
expected_z /= expected_w @ expected_z
effective_z, effective_w, effective_matrix = evaluator.rt['model'].income_transition_values(copy.deepcopy(P))
redraw = rho * np.eye(n) + (1-rho) * np.tile(w, (n, 1))
logz = np.log(z); centered = logz - w @ logz
covariances = [float((w * centered) @ np.linalg.matrix_power(matrix, lag) @ centered) for lag in range(5)]
level = z - w @ z
level_covariances = [float((w * level) @ np.linalg.matrix_power(matrix, lag) @ level) for lag in range(5)]
matrix_copies = {}
for label in ('solution', 'shared', 'evaluation'):
    value = packet.get(label)
    mapping = value if isinstance(value, dict) else vars(value) if hasattr(value, '__dict__') else {}
    for key, item in mapping.items():
        if isinstance(item, np.ndarray) and item.shape == matrix.shape:
            matrix_copies[label + '.' + key] = {'max_abs_gap_saved_Pi_z': float(np.max(np.abs(item-matrix))),
                                               'sha256_raw_array': hashlib.sha256(item.tobytes()).hexdigest()}
fields = {}
for name in ('period_years', 'age_start', 'da', 'J', 'J_R', 'use_income_types',
             'income_type_transition', 'income_shock_persistence', 'retirement_income_z_scale',
             'permanent_income_levels_enabled', 'permanent_income_log_variance',
             'psi_child', 'n_max', 'n_top', 'R', 'tau_pay'):
    value = getattr(P, name, None)
    fields[name] = value.item() if isinstance(value, np.generic) else value
for name in ('income', 'income_gross', 'age_efficiency', 'z_grid', 'z_weights'):
    value = getattr(P, name, None)
    if value is not None:
        fields[name] = np.asarray(value).tolist()
income = dict(saved_fields=fields, approved_external_estimate=approved,
    external_source=str(external_path), external_source_sha256=sha(external_path),
    matrix_shape=list(matrix.shape), row_sum_max_gap=float(np.max(np.abs(matrix.sum(axis=1)-1))),
    minimum_probability=float(matrix.min()), stationary_distribution_max_gap=float(np.max(np.abs(w@matrix-w))),
    mean_income_multiplier=float(w@z), saved_vs_effective_matrix_max_gap=float(np.max(np.abs(matrix-effective_matrix))),
    saved_vs_rouwenhorst_matrix_max_gap=float(np.max(np.abs(matrix-expected))),
    saved_vs_redraw_matrix_max_gap=float(np.max(np.abs(matrix-redraw))),
    saved_vs_approved_grid_max_gap=float(np.max(np.abs(z-expected_z))),
    saved_vs_binomial_weights_max_gap=float(np.max(np.abs(w-expected_w))),
    log_covariances_lag0_to4=covariances, level_covariances_lag0_to4=level_covariances,
    implied_log_ar1_rho=covariances[1]/covariances[0],
    implied_log_innovation_sd=math.sqrt(covariances[0]*(1-(covariances[1]/covariances[0])**2)),
    saved_matrix_copies=matrix_copies)
write('income_audit.json', income)
np.savetxt(OUT/'saved_income_transition.csv', matrix, delimiter=',', fmt='%.17g')
np.savetxt(OUT/'approved_rouwenhorst_transition.csv', expected, delimiter=',', fmt='%.17g')

# Reconstruct saved fertility count distributions without rerunning the KFE.
evaluation = packet['evaluation']
observer = read(SOURCE/'observers.json')['fertility']['uniform_birth_time']
accounting = observer['accounting']
pre = np.asarray(accounting['pre_parity_mass_by_age'])
post = np.asarray(accounting['post_parity_mass_by_age'])
for j in range(P.J):
    # Arrays at a fixed age: wealth, tenure, location, income, children-ever-born, birth-age.
    for actual, saved in ((evaluation.g_pre[:,:,:,j], pre[j]), (evaluation.g_current[:,:,:,j], post[j])):
        counts = actual.sum(axis=(0, 1, 2, 3, 5))
        assert np.max(np.abs(counts-saved)) < 1e-12
age25 = .125*pre[1] + .875*post[1]
shares = age25 / age25.sum()
observer_shares = np.array([observer['ever_born_shares_age25'][key] for key in ('0','1','2','3plus')])
assert np.max(np.abs(shares-observer_shares)) < 1e-12
fertility = dict(model_mother_share_age25=float(1-shares[0]),
    model_children_capped3_age25=float(shares@np.arange(len(shares))),
    model_children_capped3_given_mother_age25=float(shares@np.arange(len(shares))/(1-shares[0])),
    model_age25_count_shares=shares.tolist(), checkpoint_observer_count_replay=True,
    adult_entry_gate=receipt['adult_entry_gate'], normalization=receipt['normalization'])
write('fertility_checkpoint_replay.json', fertility)
fit = list(csv.DictReader((SOURCE/'target_fit.csv').open()))
loss = sum(float(row['weight'])*float(row['gap'])**2 for row in fit if row['role']=='scored')
assert abs(loss-receipt['loss']) < 1e-10
identity = dict(label='2007 stationary reference — block0506, September 28 verified export',
    inspected_checkpoint=str(checkpoint), resolved_checkpoint=str(checkpoint.resolve()),
    checkpoint_sha256=sha(checkpoint), contract_sha256=sha(B/'contract_v1/contract.json'),
    source_manifest_sha256=sha(c['source_manifest']['path']), source_root=c['source_root'],
    scientific_identity=receipt['scientific_identity'], scientific_candidate_id=receipt['scientific_candidate_id'],
    parameter_table_sha256=sha(SOURCE/'parameters.csv'), parameter_identity_matches_block0506=True,
    objective_sha256=sha(c['lanes']['primary']['objective']['path']),
    target_weight_fingerprint=receipt['target_weight_fingerprint'], export_hashes_verified=len(hashes),
    standard_plots_verified=sum(name.endswith('.png') for name in hashes),
    recomputed_primary_loss=loss, source_changes=[], economic_changes=[], model_solves=0,
    script_sha256=sha(__file__), slurm_job=os.environ['SLURM_JOB_ID'])
write('reference_identity.json', identity)
print(json.dumps(dict(identity=identity, income={k:v for k,v in income.items() if k!='saved_fields'}, fertility=fertility)), flush=True)
