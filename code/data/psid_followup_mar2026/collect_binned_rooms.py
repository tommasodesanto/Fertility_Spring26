#!/usr/bin/env python3
"""Collect aggregate binned-room diagnostics and verify before plotting."""
import csv
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).parent / 'output/first_birth_correction_review/binned_rooms'
REMOTE = '/scratch/td2248/projects/Fertility_Spring26_binned_rooms_20260912b'
ARMS = ['original_binned', 'aligned_binned']
FILES = ['run_receipt.json', 'coefficients.csv', 'covariance.csv',
         'input_support.csv', 'fitted_support.csv', 'fit_receipt.csv',
         'sample_key_hashes.json', 'completion.txt', 'estimation.log']


def collect():
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
        check=True, capture_output=True, text=True, timeout=60)
    payload = json.loads(response.stdout)
    expected = hashlib.sha256((Path(__file__).parent/'audit_binned_rooms.do').read_bytes()).hexdigest()
    states = {arm: 'pending' for arm in ARMS}
    for arm, files in payload.items():
        assert arm in ARMS and set(files) <= set(FILES)
        receipt = json.loads(files['run_receipt.json'])
        assert receipt['estimator_do_sha256'] == expected
        assert receipt['data_sha256'] == 'ede59bdb1ee42028dc97b6f35ee9d281e913ae0745bed777ce39bd6bdd9ba8eb'
        directory = OUT/arm
        directory.mkdir(parents=True, exist_ok=True)
        for name, content in files.items():
            (directory/name).write_text(content)
        states[arm] = receipt['status']
    print(json.dumps(states))
    if all(value == 'pass' for value in states.values()):
        render()


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def render():
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/psid_correction_review/matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    names = ['Dm7', 'Dm5', None, 'Dm1', 'Dp1', 'Dp3']
    labels = ['−7 / −6', '−5 / −4', '−3 / −2\nreference', '−1 / 0', '+1 / +2', '+3 / +4']
    hashes, summary = [], []
    fig, ax = plt.subplots(figsize=(9, 5.3))
    for arm, color, offset in zip(ARMS, ['#777777', '#176a99'], [-0.035, 0.035]):
        directory = OUT/arm
        receipt = json.loads((directory/'run_receipt.json').read_text())
        assert receipt['status'] == 'pass'
        hashes.append(receipt['sample_keys_sha256'])
        points = rows(directory/'coefficients.csv')
        index = {p['coefficient']: i for i, p in enumerate(points)}
        matrix = np.full((len(points), len(points)), np.nan)
        for row in rows(directory/'covariance.csv'):
            matrix[index[row['coefficient_i']], index[row['coefficient_j']]] = float(row['covariance'])
        assert np.isfinite(matrix).all() and np.allclose(matrix, matrix.T, atol=1e-10)
        assert np.linalg.eigvalsh(matrix).min() >= -1e-8
        estimates = {p['coefficient']: float(p['estimate']) for p in points}
        means = [0 if name is None else estimates[name] for name in names]
        errors = [0 if name is None else np.sqrt(matrix[index[name], index[name]]) for name in names]
        fit = rows(directory/'fit_receipt.csv')[0]
        assert np.isclose(float(fit['post_3_4_vs_pre_3_2']), means[-1])
        assert np.isclose(float(fit['standard_error']), errors[-1])
        summary.append(fit)
        label = 'Original year assignment' if arm == ARMS[0] else 'Verified survey-year assignment'
        ax.errorbar(np.arange(6)+offset, means, yerr=1.96*np.array(errors), color=color,
                    marker='o', capsize=3, linewidth=1.5, label=label)
    assert len(set(hashes)) == 1, 'Paired samples differ'
    ax.axhline(0, color='#bbbbbb', linewidth=0.8)
    ax.set(xticks=np.arange(6), xticklabels=labels, ylabel='Rooms relative to the −3/−2 window',
           xlabel='Calendar years relative to first birth', title='Two-year windows: controlled timing diagnostic')
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False)
    fig.text(.5, .015, 'Same sample and supported birth cohorts; individual/year FE; unweighted.\n'
             '2019 comparison cohort used only before 2019. Original room codes retained. Diagnostic, not a new calibration target.',
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .085, 1, 1))
    for suffix in ['png', 'pdf']:
        fig.savefig(OUT/f'binned_rooms_comparison.{suffix}', dpi=180)
    plt.close(fig)
    with (OUT/'summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(summary)
    (OUT/'verification.json').write_text(json.dumps(dict(status='pass', paired_sample_sha256=hashes[0],
        covariance_symmetric_psd=True, receipt_arithmetic_verified=True), indent=2)+'\n')


if __name__ == '__main__':
    collect()
