#!/usr/bin/env python3
"""Collect the one-change-at-a-time sequence fits from Torch, verify, and plot.

Usage: collect_rooms_sequence.py [--local DIR]
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
OUT = HERE / 'output/first_birth_correction_review/sequence_designs'
REMOTE = '/scratch/td2248/projects/Fertility_Spring26_rooms_sequence_20260923'
DATA_SHA256 = '37cf65e3e9c94f11232383d0f0200d18ac0c30442864757cb48297b6991cfa58'
ARMS = ['S0', 'B21', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S0c', 'B21c', 'S4a', 'S1c', 'S2c', 'S3c', 'S4ac', 'S4c']
ACCEPTED_ESTIMATOR_SHA256 = {'c22816ab3cf2a01336b530d83f26856d595bc08739c92e0a50135e3fbed20a12', '2ec590ea8b13e4e001281b62e4a7f8b8e1a1f41c142578d3c288451a613530ab'}  # first eight arms ran under this revision
LABELS = {'S0': 'S0 original spec, corrected dates (=window_aligned)',
          'B21': 'B21 same, baseline −2/−1',
          'S1': 'S1 + confirmed-childless controls only',
          'S2': 'S2 + women, current ref/spouse',
          'S3': 'S3 + one woman per single-FU household-year',
          'S4': 'S4 + biological first-birth history',
          'S5': 'S5 + non-room codes to missing',
          'S6': 'S6 + PSID weights [pw=IW]',
          'S0c': 'S0c original spec, corrected dates, codes cleaned',
          'B21c': 'B21c same, baseline −2/−1, codes cleaned',
          'S4a': 'S4a S3 + biological history, no entry rule',
          'S1c': 'S1c confirmed-childless controls, codes cleaned',
          'S2c': 'S2c + women, current ref/spouse, codes cleaned',
          'S3c': 'S3c + one woman per single-FU household-year, codes cleaned',
          'S4ac': 'S4ac + biological history (no entry rule), codes cleaned',
          'S4c': 'S4c + first-birth-after-entry rule, codes cleaned'}
FILES = ['run_receipt.json', 'coefficients.csv', 'covariance.csv', 'input_support.csv',
         'fitted_support.csv', 'fit_receipt.csv', 'sample_key_hashes.json', 'completion.txt', 'estimation.log']
NAMES = ['Wleft', 'Wm2', 'Wm1', None, 'Wp1', 'Wp2', 'Wp3', 'Wp4', 'Wp5', 'Wp6', 'Wright']


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
    response = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=10', 'torch',
                               'python3 -c ' + shlex.quote(code)], check=True, capture_output=True, text=True, timeout=120)
    payload = json.loads(response.stdout)
    expected = hashlib.sha256((HERE/'audit_rooms_sequence.do').read_bytes()).hexdigest()
    states = {arm: 'pending' for arm in ARMS}
    for arm, files in payload.items():
        assert arm in ARMS and set(files) <= set(FILES)
        receipt = json.loads(files['run_receipt.json'])
        assert receipt['estimator_do_sha256'] in ACCEPTED_ESTIMATOR_SHA256 | {expected}, 'estimator hash mismatch'
        assert receipt['data_sha256'] == DATA_SHA256, 'data hash mismatch'
        d = out/arm
        d.mkdir(parents=True, exist_ok=True)
        for name, content in files.items():
            (d/name).write_text(content)
        states[arm] = receipt['status']
    print(json.dumps(states))
    return states


def load_arm(d):
    points = rows(d/'coefficients.csv')
    index = {p['coefficient']: i for i, p in enumerate(points)}
    matrix = np.full((len(points), len(points)), np.nan)
    for row in rows(d/'covariance.csv'):
        matrix[index[row['coefficient_i']], index[row['coefficient_j']]] = float(row['covariance'])
    assert np.isfinite(matrix).all() and np.allclose(matrix, matrix.T, atol=1e-10)
    assert np.linalg.eigvalsh(matrix).min() >= -1e-8
    est = {p['coefficient']: float(p['estimate']) for p in points}
    se = {p['coefficient']: float(np.sqrt(matrix[index[p['coefficient']], index[p['coefficient']]])) for p in points}
    fit = rows(d/'fit_receipt.csv')[0]
    assert np.isclose(float(fit['headline_effect']), est[fit['headline']])
    assert np.isclose(float(fit['headline_se']), se[fit['headline']])
    assert int(fit['fitted_unsupported_rows']) == 0
    return est, se, fit


def render(out, label='Torch full fits', check_receipts=True, arms=None, stem='sequence_comparison', title='One change at a time: corrected dates, two-year windows'):
    arms = arms or ARMS
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/psid_correction_review/matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(11, 6))
    summary = []
    cmap = plt.get_cmap('viridis')
    for i, arm in enumerate(arms):
        d = out/arm
        if not (d/'fit_receipt.csv').exists():
            continue
        if check_receipts:
            assert json.loads((d/'run_receipt.json').read_text())['status'] == 'pass'
        est, se, fit = load_arm(d)
        lo, hi = int(fit['baseline_lo']), int(fit['baseline_hi'])
        labels = []
        for n in NAMES:
            if n is None:
                labels.append(f'{lo}/{hi}')
            elif n == 'Wleft':
                labels.append(f'≤{hi-6}')
            elif n == 'Wright':
                labels.append(f'≥{hi+13}')
            else:
                k = int(n[2:]); off = {'Wm2': (-5, -4), 'Wm1': (-3, -2)}.get(n, (2*k-1, 2*k))
                labels.append(f'{hi+off[0]:+d}/{hi+off[1]:+d}')
        means = [0 if n is None else est[n] for n in NAMES]
        errs = [0 if n is None else se[n] for n in NAMES]
        style = dict(color=cmap(i/max(len(arms)-1, 1)), marker='o', linewidth=1.3, capsize=2)
        if arm.startswith('B21'):
            style.update(linestyle='--')
        if arm.endswith('c') or arm == 'S4a':
            style.update(marker='s', linestyle=':')
        ax.errorbar(np.arange(len(NAMES))+(i-3.5)*0.03, means, yerr=1.96*np.array(errs),
                    label=f"{LABELS[arm]} (N={int(float(fit['observations'])):,})", **style)
        rec = dict(arm=arm, label=LABELS[arm], baseline=f'{lo}/{hi}', observations=fit['observations'],
                   clusters=fit['clusters'], never_treated_rows=fit['never_treated_rows'], weights=fit['weights'],
                   headline=fit['headline'], headline_window=labels[NAMES.index(fit['headline'])],
                   headline_effect=fit['headline_effect'], headline_se=fit['headline_se'],
                   reference_mean_rooms=fit['reference_mean_rooms'], dropped_cohorts=fit['dropped_cohorts'],
                   runtime_seconds=fit['runtime_seconds'])
        for n in ['Wm1', 'Wp1', 'Wp2', 'Wp3', 'Wp4']:
            rec[f'{n}_estimate'] = est[n]; rec[f'{n}_se'] = se[n]
        summary.append(rec)
    ax.axhline(0, color='#bbbbbb', linewidth=0.8)
    ax.set(xticks=np.arange(len(NAMES)), xticklabels=['≤−8', '−7/−6', '−5/−4', 'base', '+1/+2', '+3/+4', '+5/+6', '+7/+8', '+9/+10', '+11/+12', '≥+13'],
           ylabel='Rooms relative to the baseline window', xlabel='Windows relative to first birth (labels for the −3/−2 baseline; B21 is shifted one year later)',
           title=title)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False, fontsize=7.5, loc='upper left')
    fig.text(.5, .01, f'{label}. Diagnostic, not a calibration target.', ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .04, 1, 1))
    for suffix in ['png', 'pdf']:
        fig.savefig(out/f'{stem}.{suffix}', dpi=180)
    plt.close(fig)
    with (out/('summary.csv' if stem == 'sequence_comparison' else f'{stem}_summary.csv')).open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=list(summary[0]), lineterminator='\n')
        w.writeheader(); w.writerows(summary)
    for r in summary:
        print(f"{r['arm']:>4} {r['headline_window']:>6} {float(r['headline_effect']):.3f} ({float(r['headline_se']):.3f})  N={int(float(r['observations'])):>7,}  {r['label']}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--local', type=Path)
    parser.add_argument('--subset', help='comma-separated arms to render from already collected outputs')
    args = parser.parse_args()
    if args.subset:
        render(OUT, arms=args.subset.split(','), stem='sequence_cleaned_chain', title='One change at a time, PSID non-room codes cleaned: corrected dates, two-year windows')
    elif args.local:
        render(args.local.resolve(), label='Local synthetic smoke', check_receipts=False)
    else:
        OUT.mkdir(parents=True, exist_ok=True)
        states = collect(OUT)
        if all(v == 'pass' for v in states.values()):
            render(OUT)
            (OUT/'verification.json').write_text(json.dumps(dict(
                status='pass', data_sha256=DATA_SHA256,
                estimator_do_sha256=hashlib.sha256((HERE/'audit_rooms_sequence.do').read_bytes()).hexdigest(),
                covariance_symmetric_psd=True, receipt_arithmetic_verified=True), indent=2)+'\n')
