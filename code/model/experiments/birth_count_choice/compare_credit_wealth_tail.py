"""No-solve recreation of the negative-financial-position age panel from saved credit runs."""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_cap2_at_binary_winner import inputs_and_receipt
OUT = ROOT / 'output/model/experiments/birth_count_choice/credit_at_binary_winner_v2/wealth_tail_comparison'
REFERENCE = ROOT / 'output/model/local_solution/cases/20261003T175652812716Z_b1c72f13/aggregate_plots/model_data_assessment/plotted_long.csv'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    _, _, _, _, grid, _, _ = inputs_and_receipt()
    out = []
    pins = {}
    for version in ('v1', 'v2'):
        for arm in ('phi_080', 'phi_100'):
            folder = ROOT / f'output/model/experiments/birth_count_choice/credit_at_binary_winner_{version}/{arm}'
            arrays = folder / 'solution_arrays.npz'
            receipt = json.loads((folder / 'solve_completed.json').read_text())
            digest = sha(arrays)
            if digest != receipt['arrays_sha256']:
                raise RuntimeError(f'Saved array identity failed: {version}/{arm}')
            pins[f'{version}/{arm}'] = dict(arrays_sha256=digest, solve_receipt_sha256=sha(folder / 'solve_completed.json'))
            with np.load(arrays, allow_pickle=False) as saved:
                mass = np.asarray(saved['g_beginning_distribution'], float)
            if mass.shape[0] != len(grid) or mass.shape[3] != 17:
                raise RuntimeError('Saved distribution and wealth grid differ')
            denom = mass.sum(axis=(0, 1, 2, 4, 5, 6))
            negative = mass[grid < 0].sum(axis=(0, 1, 2, 4, 5, 6))
            if np.any(denom <= 0):
                raise RuntimeError('Empty age cell')
            for j, (m, n) in enumerate(zip(denom, negative)):
                out.append(dict(series=f'{version}_{arm}', age=20 + 4*j, share_negative_financial_position=float(n/m), mass=float(m)))
    reference = list(csv.DictReader(REFERENCE.open()))
    benchmark = {}
    for series in ('Model', 'Data'):
        benchmark[series] = [(float(r['x']), float(r['value'])) for r in reference
                             if r['panel'] == 'resource_negative_age' and r['series'] == series]
        if len(benchmark[series]) != 17:
            raise RuntimeError(f'Screenshot reference rows missing: {series}')
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / 'model_shares_by_age.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(out[0])); writer.writeheader(); writer.writerows(out)
    with (OUT / 'screenshot_reference_by_age.csv').open('w', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(('series', 'age', 'share_negative_financial_position'))
        for series, pairs in benchmark.items():
            for age, value in pairs: writer.writerow((series, age, value))
    fig, axs = plt.subplots(1, 2, figsize=(11.2, 4.5), sharey=True)
    for ax, arm, title in zip(axs, ('phi_080', 'phi_100'), ('Financed share 0.8', 'Financed share 1.0')):
        for version, color, style in (('v1', '#bc6442', '-'), ('v2', '#2864a4', '--')):
            rows = [r for r in out if r['series'] == f'{version}_{arm}']
            ax.plot([r['age'] for r in rows], [r['share_negative_financial_position'] for r in rows],
                    style, color=color, lw=2.2, label=f'Estate-A {version}')
        ax.plot(*zip(*benchmark['Model']), ':', color='#777777', lw=1.8,
                label='Screenshot model (different point)')
        ax.plot(*zip(*benchmark['Data']), 'o', color='#444444', ms=3.2,
                label='Screenshot PSID data')
        ax.set(title=title, xlabel='Age (four-year cell midpoint)', xlim=(20, 84), ylim=(0, 1))
        ax.grid(alpha=.2)
    axs[0].set_ylabel('Share with beginning net financial position $b<0$')
    axs[1].legend(loc='upper left', fontsize=8, frameon=False)
    fig.suptitle('Negative financial position by age: estate fix leaves old-age debt share unchanged')
    fig.tight_layout()
    fig.savefig(OUT / 'negative_financial_position_by_age.png', dpi=170)
    plt.close(fig)
    pin = dict(status='saved_arrays_only_no_model_solves', definition='g_beginning_distribution mass at b<0 / total beginning age mass',
               timing='beginning of period before current saving and tenure choice', ages='20,24,...,84 midpoints',
               reducer_sha256=sha(__file__),
               source_array_receipts=pins, wealth_grid_sha256=hashlib.sha256(np.asarray(grid).tobytes()).hexdigest(),
               screenshot_source=str(REFERENCE.relative_to(ROOT)), screenshot_source_sha256=sha(REFERENCE),
               screenshot_model_is_different_calibration_point=True)
    (OUT / 'source_receipt.json').write_text(json.dumps(pin, indent=2, sort_keys=True)+'\n')
    for arm in ('phi_080','phi_100'):
        b = next(r for r in out if r['series']==f'v1_{arm}' and r['age']==84)
        a = next(r for r in out if r['series']==f'v2_{arm}' and r['age']==84)
        print(arm,'age84',b['share_negative_financial_position'],a['share_negative_financial_position'])

if __name__ == '__main__': main()
