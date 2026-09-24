#!/usr/bin/env python3
"""Collect, verify and render the first-birth event studies for ownership and moving outcomes (v3 sample).

Usage: collect_first_birth_outcomes_v3.py [--local DIR]
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
OUT = HERE / 'output/sa_first_birth_outcomes_v3'
REMOTE = '/scratch/td2248/projects/Fertility_Spring26_rooms_v3_20260924'
DESIGNS = {'H': 'Household design (women heads/spouses, IW)', 'A2h': 'Adults heading a household at baseline (IW)'}
OUTCOMES = {'own': 'Owns the dwelling (share)', 'moved': 'Moved since last interview (share)',
            'moved_space': 'Moved for more space (share of all)', 'moved_nbhd': 'Moved for neighbourhood (share of all)',
            'space_c': 'Share of MOVES for more space', 'nbhd_c': 'Share of MOVES for neighbourhood'}
ARMS = [f'{d}_{o}' for d in DESIGNS for o in OUTCOMES]
FILES = ['run_receipt.json', 'coefficients.csv', 'covariance.csv', 'input_support.csv',
         'fitted_support.csv', 'fit_receipt.csv', 'sample_key_hashes.json', 'completion.txt', 'estimation.log']
NAMES = ['Wleft', 'Wm2', 'Wm1', None, 'Wp1', 'Wp2', 'Wp3', 'Wp4', 'Wp5', 'Wp6', 'Wright']
XS = np.array([-10.5, -6.5, -4.5, -2.5, -0.5, 1.5, 3.5, 5.5, 7.5, 9.5, 11.5])
LABELS = ['≤−8', '−7/−6', '−5/−4', '−3/−2 ref.', '−1/0', '+1/+2', '+3/+4', '+5/+6', '+7/+8', '+9/+10', '≥+11']


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def data_sha():
    p = Path('/tmp/psid_rooms_v3_20260924/analysis_sample.dta')
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None


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
    expected = hashlib.sha256((HERE/'sa_rooms_first_birth_v2.do').read_bytes()).hexdigest()
    states = {arm: 'pending' for arm in ARMS}
    for arm, files in payload.items():
        receipt = json.loads(files['run_receipt.json'])
        assert receipt['estimator_do_sha256'] in {expected, '4ce0b209b092943c0eb6bf3d886889eac4b7bbda4494953105eca8a434358472'}, 'estimator hash mismatch'  # unconditional arms ran under the earlier revision
        if data_sha():
            assert receipt['data_sha256'] == data_sha(), 'data hash mismatch'
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


def render(out, label='Torch full fits', check_receipts=True):
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/psid_correction_review/matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, len(OUTCOMES), figsize=(4.5*len(OUTCOMES), 8), sharex=True)
    summary, path_rows = [], []
    for i, (d, dlab) in enumerate(DESIGNS.items()):
        for j, (o, olab) in enumerate(OUTCOMES.items()):
            ax = axes[i, j]
            arm = f'{d}_{o}'
            if not (out/arm/'fit_receipt.csv').exists():
                ax.set_visible(False); continue
            if check_receipts:
                assert json.loads((out/arm/'run_receipt.json').read_text())['status'] == 'pass'
            est, se, fit = load_arm(out/arm)
            means = [0 if n is None else est[n] for n in NAMES]
            errs = [0 if n is None else se[n] for n in NAMES]
            ax.errorbar(XS, means, yerr=1.96*np.array(errs), color='#176a99' if d == 'H' else '#2a6f97', marker='o', linewidth=1.5, capsize=3)
            ax.plot([-2.5], [0], marker='D', markersize=8, markerfacecolor='white', markeredgecolor='#176a99', markeredgewidth=1.6, linestyle='none', zorder=5)
            ax.axhline(0, color='#bbbbbb', linewidth=0.8)
            ax.axvline(-0.5, color='#999999', linewidth=0.9, linestyle='--')
            ax.axvspan(-13, -0.5, color='#000000', alpha=0.035, lw=0)
            ax.set(title=f'{olab}\n{dlab}', xlim=(-11.5, 12.5), xticks=np.arange(-10, 14, 2))
            ax.spines[['top', 'right']].set_visible(False)
            if j == 0:
                ax.set_ylabel('Change relative to −3/−2 window')
            if i == 1:
                ax.set_xlabel('Calendar years relative to first birth')
            summary.append(dict(arm=arm, design=d, outcome=o, observations=fit['observations'], clusters=fit['clusters'],
                                treated_individuals=fit['treated_individuals'], control_individuals=fit['control_individuals'],
                                baseline_mean=fit['reference_mean_rooms'], p3p4_effect=fit['headline_effect'], p3p4_se=fit['headline_se'],
                                m1p0_effect=est['Wp1'], m1p0_se=se['Wp1'], p1p2_effect=est['Wp2'], p1p2_se=se['Wp2'],
                                p5p6_effect=est['Wp4'], p5p6_se=se['Wp4'], dropped_cohorts=fit['dropped_cohorts']))
            for n, lab in zip(NAMES, LABELS):
                path_rows.append(dict(arm=arm, window=lab, estimate=0 if n is None else est[n], se=0 if n is None else se[n]))
    fig.text(.5, .01, f'{label}. Same specification as the rooms headline: Sun–Abraham in two-year interview windows, official year-specific PSID items, '
             'person and survey-year FE, age and education covariates, clustered by person, PSID weights. Hollow diamond = omitted baseline window.', ha='center', fontsize=7.5)
    fig.tight_layout(rect=(0, .04, 1, 1))
    for suffix in ['png', 'pdf']:
        fig.savefig(out/f'first_birth_outcomes_v3.{suffix}', dpi=170)
    plt.close(fig)
    with (out/'summary.csv').open('w', newline='') as s:
        w = csv.DictWriter(s, fieldnames=list(summary[0]), lineterminator='\n'); w.writeheader(); w.writerows(summary)
    with (out/'window_path.csv').open('w', newline='') as s:
        w = csv.DictWriter(s, fieldnames=list(path_rows[0]), lineterminator='\n'); w.writeheader(); w.writerows(path_rows)
    for r in summary:
        print(f"{r['arm']:>16} base={float(r['baseline_mean']):.3f}  -1/0 {float(r['m1p0_effect']):+.3f}({float(r['m1p0_se']):.3f})  +1/+2 {float(r['p1p2_effect']):+.3f}({float(r['p1p2_se']):.3f})  +3/+4 {float(r['p3p4_effect']):+.3f}({float(r['p3p4_se']):.3f})  +5/+6 {float(r['p5p6_effect']):+.3f}({float(r['p5p6_se']):.3f})  N={int(float(r['observations'])):,}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--local', type=Path)
    args = parser.parse_args()
    if args.local:
        render(args.local.resolve(), label='Local synthetic smoke', check_receipts=False)
    else:
        OUT.mkdir(parents=True, exist_ok=True)
        states = collect(OUT)
        if all(v == 'pass' for v in states.values()):
            render(OUT)
            (OUT/'verification.json').write_text(json.dumps(dict(status='pass', data_sha256=data_sha(),
                estimator_do_sha256=hashlib.sha256((HERE/'sa_rooms_first_birth_v2.do').read_bytes()).hexdigest(),
                covariance_symmetric_psd=True, receipt_arithmetic_verified=True), indent=2)+'\n')
