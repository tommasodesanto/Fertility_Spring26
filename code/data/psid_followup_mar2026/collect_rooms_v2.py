#!/usr/bin/env python3
"""Collect, verify and render the version-2 first-birth rooms event study.

Usage: collect_rooms_v2.py            (collect from Torch, verify, render)
       collect_rooms_v2.py --local DIR (render a local results folder, e.g. toy)
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
OUT = HERE / 'output/sa_rooms_first_birth_v2'
REMOTE = '/scratch/td2248/projects/Fertility_Spring26_rooms_v2_20260924'
ARMS = ['H', 'H21', 'A', 'Hentry', 'Hctrl', 'Hshift', 'Ashift', 'Hb4', 'Hb6']
PLOT_ARMS = ['H', 'H21', 'A', 'Hb4', 'Hb6']
ACCEPTED_ESTIMATOR_SHA256 = {'20ca592a4928bcbcc6e0b41c01507939d8cfe822a3e16f98078ebf974b36f2ed', '46010a9b799d6750407e18740913b318be8a9c99f389d55f026fa97ccd9da292'}  # H, H21, A ran under this revision
FILES = ['run_receipt.json', 'coefficients.csv', 'covariance.csv', 'input_support.csv',
         'fitted_support.csv', 'fit_receipt.csv', 'sample_key_hashes.json', 'completion.txt', 'estimation.log']
NAMES = ['Wleft', 'Wm2', 'Wm1', None, 'Wp1', 'Wp2', 'Wp3', 'Wp4', 'Wp5', 'Wp6', 'Wright']


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def data_sha():
    p = Path('/tmp/psid_rooms_v2_20260924/analysis_sample.dta')
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
        assert receipt['estimator_do_sha256'] in ACCEPTED_ESTIMATOR_SHA256 | {expected}, 'estimator hash mismatch'
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


def labels_for(hi):
    out = []
    for n in NAMES:
        if n is None:
            out.append(f'{hi-1}/{hi}\nref.')
        elif n == 'Wleft':
            out.append(f'≤{hi-6}')
        elif n == 'Wright':
            out.append(f'≥{hi+13}')
        else:
            k = int(n[2:]); off = {'Wm2': (-5, -4), 'Wm1': (-3, -2)}.get(n, (2*k-1, 2*k))
            out.append(f'{hi+off[0]:+d}/{hi+off[1]:+d}')
    return out


def render(out, label='Torch full fits', check_receipts=True):
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/psid_correction_review/matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.4), sharey=True)
    path_rows, summary = [], []
    panels = [(axes[0], ['H', 'H21'], 'Household design: one woman per household-year, IW weights'),
              (axes[1], [a for a in ['H', 'Hb4', 'Hb6'] if (out/a/'fit_receipt.csv').exists()], 'Household design: earlier baseline windows'),
              (axes[2], ['A'], 'All-adult design (robustness): every adult, unweighted')]
    seen = set()
    for ax, arms, title in panels:
        for arm in arms:
            d = out/arm
            if check_receipts:
                assert json.loads((d/'run_receipt.json').read_text())['status'] == 'pass'
            est, se, fit = load_arm(d)
            hi = int(fit['baseline_hi'])
            means = [0 if n is None else est[n] for n in NAMES]
            errs = [0 if n is None else se[n] for n in NAMES]
            # x axis in calendar years relative to birth: window midpoint
            xs = np.array([hi-8.5, hi-4.5, hi-2.5, hi-0.5, hi+1.5, hi+3.5, hi+5.5, hi+7.5, hi+9.5, hi+11.5, hi+13.5])
            style = dict(color='#176a99', marker='o', linewidth=1.6, capsize=3, label=f'baseline {hi-1}/{hi} (headline)')
            if arm == 'H21':
                style = dict(color='#b5541c', marker='s', linewidth=1.2, linestyle='--', capsize=3, label=f'baseline {hi-1}/{hi} (sensitivity), x shifted one year')
            if arm == 'A':
                style = dict(color='#555555', marker='o', linewidth=1.6, capsize=3, label=f'baseline {hi-1}/{hi}')
            if arm == 'Hb4':
                style = dict(color='#2a9d8f', marker='^', linewidth=1.3, linestyle='-.', capsize=3, label='baseline −5/−4')
            if arm == 'Hb6':
                style = dict(color='#6a4c93', marker='v', linewidth=1.3, linestyle=':', capsize=3, label='baseline −7/−6')
            ax.errorbar(xs, means, yerr=1.96*np.array(errs), **style)
            ax.plot([hi-0.5], [0], marker='D', markersize=9, markerfacecolor='white', markeredgecolor=style['color'], markeredgewidth=1.8, linestyle='none', zorder=5)
            if arm not in seen:
                seen.add(arm)
                summary.append(dict(arm=arm, baseline=f'{hi-1}/{hi}', observations=fit['observations'], clusters=fit['clusters'],
                                    treated_individuals=fit['treated_individuals'], control_individuals=fit['control_individuals'],
                                    weights=fit['weights'], headline_window=labels_for(hi)[NAMES.index(fit['headline'])].replace('\n', ' '),
                                    headline_effect=fit['headline_effect'], headline_se=fit['headline_se'],
                                    reference_mean_rooms=fit['reference_mean_rooms'], dropped_cohorts=fit['dropped_cohorts'],
                                    runtime_seconds=fit['runtime_seconds']))
                for n, lab in zip(NAMES, labels_for(hi)):
                    path_rows.append(dict(arm=arm, window=lab.replace('\n', ' '), estimate=0 if n is None else est[n],
                                          se=0 if n is None else se[n]))
        ax.axhline(0, color='#bbbbbb', linewidth=0.8)
        ax.axvline(-0.5, color='#999999', linewidth=0.9, linestyle='--')
        ax.axvspan(-13, -0.5, color='#000000', alpha=0.035, lw=0)
        ax.text(-0.7, ax.get_ylim()[1]*0.97 if ax.get_ylim()[1] > 0 else 1.4, 'birth', ha='right', va='top', fontsize=8, color='#666666')
        ax.set(xticks=np.arange(-10, 14, 2), xlim=(-11.5, 12.5), title=title,
               xlabel='Calendar years relative to first birth')
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(frameon=False, fontsize=8, loc='upper left')
    axes[0].set_ylabel('Rooms relative to the baseline window')
    fig.text(.5, .01, f'{label}. Points at two-year window midpoints; hollow diamond = omitted baseline window (zero by construction). Sun–Abraham in two-year interview windows; '
             'rooms from the official year-specific PSID variables; non-room codes to missing; person and survey-year FE; age and education covariates; clustered by person.', ha='center', fontsize=7.5)
    fig.tight_layout(rect=(0, .05, 1, 1))
    for suffix in ['png', 'pdf']:
        fig.savefig(out/f'first_birth_rooms_v2.{suffix}', dpi=180)
    plt.close(fig)
    with (out/'summary.csv').open('w', newline='') as s:
        w = csv.DictWriter(s, fieldnames=list(summary[0]), lineterminator='\n'); w.writeheader(); w.writerows(summary)
    with (out/'window_path.csv').open('w', newline='') as s:
        w = csv.DictWriter(s, fieldnames=list(path_rows[0]), lineterminator='\n'); w.writeheader(); w.writerows(path_rows)
    for arm in [a for a in ARMS if a not in PLOT_ARMS and (out/a/'fit_receipt.csv').exists()]:
        est, se, fit = load_arm(out/arm)
        hi = int(fit['baseline_hi'])
        summary.append(dict(arm=arm, baseline=f'{hi-1}/{hi}', observations=fit['observations'], clusters=fit['clusters'],
                            treated_individuals=fit['treated_individuals'], control_individuals=fit['control_individuals'],
                            weights=fit['weights'], headline_window=labels_for(hi)[NAMES.index(fit['headline'])].replace('\n', ' '),
                            headline_effect=fit['headline_effect'], headline_se=fit['headline_se'],
                            reference_mean_rooms=fit['reference_mean_rooms'], dropped_cohorts=fit['dropped_cohorts'],
                            runtime_seconds=fit['runtime_seconds']))
        for n, lab in zip(NAMES, labels_for(hi)):
            path_rows.append(dict(arm=arm, window=lab.replace('\n', ' '), estimate=0 if n is None else est[n], se=0 if n is None else se[n]))
    with (out/'summary.csv').open('w', newline='') as s:
        w = csv.DictWriter(s, fieldnames=list(summary[0]), lineterminator='\n'); w.writeheader(); w.writerows(summary)
    with (out/'window_path.csv').open('w', newline='') as s:
        w = csv.DictWriter(s, fieldnames=list(path_rows[0]), lineterminator='\n'); w.writeheader(); w.writerows(path_rows)
    for r in summary:
        print(f"{r['arm']:>4} {r['headline_window']:>8} {float(r['headline_effect']):.3f} ({float(r['headline_se']):.3f})  N={int(float(r['observations'])):>7,}  treated={int(float(r['treated_individuals'])):,} controls={int(float(r['control_individuals'])):,}")


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
            (OUT/'verification.json').write_text(json.dumps(dict(
                status='pass', data_sha256=data_sha(),
                estimator_do_sha256=hashlib.sha256((HERE/'sa_rooms_first_birth_v2.do').read_bytes()).hexdigest(),
                covariance_symmetric_psd=True, receipt_arithmetic_verified=True), indent=2)+'\n')
