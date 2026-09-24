#!/usr/bin/env python3
"""Collect aggregate outputs of the window / -2-only room designs and plot.

Usage:
  collect_rooms_window_designs.py                 # collect from Torch, verify, plot
  collect_rooms_window_designs.py --local DIR     # verify/plot a local results folder
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'output/first_birth_correction_review/window_designs'
REMOTE = '/scratch/td2248/projects/Fertility_Spring26_rooms_window_designs_20260923'
DATA_SHA256 = '3466e1b734d43f3991d96b3360e261fad6c0df770e00aec864f1dcace93b340b'
ARMS = ['window_original', 'window_aligned', 'm2only_original', 'm2only_aligned']
FILES = ['run_receipt.json', 'coefficients.csv', 'covariance.csv', 'input_support.csv',
         'fitted_support.csv', 'fit_receipt.csv', 'sample_key_hashes.json',
         'completion.txt', 'estimation.log']
WINDOW_NAMES = ['Dleft', 'Dm7', 'Dm5', None, 'Dm1', 'Dp1', 'Dp3', 'Dp5', 'Dp7', 'Dp9', 'Dright']
WINDOW_LABELS = ['≤−8', '−7/−6', '−5/−4', '−3/−2\nref.', '−1/0', '+1/+2', '+3/+4', '+5/+6',
                 '+7/+8', '+9/+10', '≥+11']
ANNUAL = [('F7event', -7), ('F6event', -6), ('F5event', -5), ('F4event', -4), ('F3event', -3),
          (None, -2), ('F1event', -1)] + [(f'L{k}event', k) for k in range(12)]


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def collect(out):
    code = f'''import json,pathlib
root=pathlib.Path({REMOTE!r})/'results'
payload={{}}
for arm in {ARMS!r}:
 p=root/arm
 if (p/'run_receipt.json').exists():
  payload[arm]={{n:(p/n).read_text() for n in {FILES!r} if (p/n).exists()}}
print(json.dumps(payload))
'''
    response = subprocess.run(
        ['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=10', 'torch',
         'python3 -c ' + shlex.quote(code)],
        check=True, capture_output=True, text=True, timeout=90)
    payload = json.loads(response.stdout)
    expected = hashlib.sha256((HERE/'audit_rooms_window_designs.do').read_bytes()).hexdigest()
    states = {arm: 'pending' for arm in ARMS}
    for arm, files in payload.items():
        assert arm in ARMS and set(files) <= set(FILES)
        receipt = json.loads(files['run_receipt.json'])
        assert receipt['estimator_do_sha256'] == expected, 'estimator hash mismatch'
        assert receipt['data_sha256'] == DATA_SHA256, 'data hash mismatch'
        directory = out/arm
        directory.mkdir(parents=True, exist_ok=True)
        for name, content in files.items():
            (directory/name).write_text(content)
        states[arm] = receipt['status']
    print(json.dumps(states))
    return states


def load_arm(directory):
    points = rows(directory/'coefficients.csv')
    index = {p['coefficient']: i for i, p in enumerate(points)}
    matrix = np.full((len(points), len(points)), np.nan)
    for row in rows(directory/'covariance.csv'):
        matrix[index[row['coefficient_i']], index[row['coefficient_j']]] = float(row['covariance'])
    assert np.isfinite(matrix).all() and np.allclose(matrix, matrix.T, atol=1e-10)
    assert np.linalg.eigvalsh(matrix).min() >= -1e-8
    estimates = {p['coefficient']: float(p['estimate']) for p in points}
    errors = {p['coefficient']: float(np.sqrt(matrix[index[p['coefficient']], index[p['coefficient']]]))
              for p in points}
    fit = rows(directory/'fit_receipt.csv')[0]
    assert np.isclose(float(fit['headline_effect']), estimates[fit['headline']])
    assert np.isclose(float(fit['headline_se']), errors[fit['headline']])
    assert int(fit['fitted_unsupported_rows']) == 0, 'a treated cohort lost its reference in e(sample)'
    return estimates, errors, fit


def render(out, label='Torch full fits', check_pairs=True):
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/psid_correction_review/matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    summary = []
    colors = {'original': '#777777', 'aligned': '#b5541c'}
    for design, ax in zip(['window', 'm2only'], axes):
        hashes = []
        for assignment in ['original', 'aligned']:
            arm = f'{design}_{assignment}'
            directory = out/arm
            estimates, errors, fit = load_arm(directory)
            if check_pairs:
                receipt = json.loads((directory/'run_receipt.json').read_text())
                assert receipt['status'] == 'pass'
                hashes.append(receipt['sample_keys_sha256'])
            offset = -0.04 if assignment == 'original' else 0.04
            if design == 'window':
                names, labels = WINDOW_NAMES, WINDOW_LABELS
                xs = np.arange(len(names))
            else:
                names = [n for n, _ in ANNUAL]
                labels = [str(k) for _, k in ANNUAL]
                xs = np.array([k for _, k in ANNUAL], dtype=float)
            means = [0 if n is None else estimates[n] for n in names]
            ses = [0 if n is None else errors[n] for n in names]
            text = 'Original year assignment' if assignment == 'original' else 'Verified survey-year assignment'
            ax.errorbar(xs+offset, means, yerr=1.96*np.array(ses), color=colors[assignment],
                        marker='o', capsize=3, linewidth=1.4, label=text)
            record = dict(arm=arm, observations=fit['observations'], clusters=fit['clusters'],
                          headline=fit['headline'], headline_effect=fit['headline_effect'],
                          headline_se=fit['headline_se'], reference_mean_rooms=fit['reference_mean_rooms'],
                          dropped_unsupported_rows=fit['dropped_unsupported_rows'],
                          dropped_cohorts=fit['dropped_cohorts'], runtime_seconds=fit['runtime_seconds'])
            for n in (['Dp1', 'Dp3', 'Dp5'] if design == 'window' else ['L0event', 'L2event', 'L3event', 'L4event', 'L6event', 'L8event']):
                record[f'{n}_estimate'] = estimates[n]
                record[f'{n}_se'] = errors[n]
            summary.append(record)
        if check_pairs:
            assert len(set(hashes)) == 1, f'{design}: paired samples differ'
        ax.axhline(0, color='#bbbbbb', linewidth=0.8)
        if design == 'window':
            ax.set(xticks=np.arange(len(WINDOW_LABELS)), xticklabels=WINDOW_LABELS,
                   ylabel='Rooms relative to the −3/−2 window',
                   xlabel='Calendar years relative to first birth',
                   title='Two-year interview windows, baseline −3/−2')
        else:
            ax.set(xticks=[k for _, k in ANNUAL], ylabel='Rooms relative to year −2',
                   xlabel='Years relative to first birth',
                   title='Annual event time, cohorts observed at −2 only')
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(frameon=False, fontsize=9)
    fig.text(.5, .01, f'{label}. Author\'s original Sun–Abraham specification: all adults, unweighted, '
             'ID/year FE, age and education covariates, last cohort as control, original room codes. '
             'Only the room date assignment differs within a panel. Diagnostic, not a calibration target.',
             ha='center', fontsize=7.5)
    fig.tight_layout(rect=(0, .05, 1, 1))
    for suffix in ['png', 'pdf']:
        fig.savefig(out/f'window_designs_comparison.{suffix}', dpi=180)
    plt.close(fig)
    keys = sorted({k for r in summary for k in r}, key=lambda k: (k not in summary[0], k))
    with (out/'summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]) + [k for k in keys if k not in summary[0]],
                                lineterminator='\n')
        writer.writeheader()
        writer.writerows(summary)
    print((out/'summary.csv').read_text())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--local', type=Path, help='render a local results folder (e.g. toy smoke) without cluster collection')
    args = parser.parse_args()
    if args.local:
        render(args.local.resolve(), label='Local synthetic smoke', check_pairs=False)
    else:
        OUT.mkdir(parents=True, exist_ok=True)
        states = collect(OUT)
        if all(value == 'pass' for value in states.values()):
            render(OUT)
            (OUT/'verification.json').write_text(json.dumps(dict(
                status='pass', data_sha256=DATA_SHA256,
                estimator_do_sha256=hashlib.sha256((HERE/'audit_rooms_window_designs.do').read_bytes()).hexdigest(),
                covariance_symmetric_psd=True, receipt_arithmetic_verified=True,
                paired_samples_identical_within_design=True), indent=2)+'\n')
