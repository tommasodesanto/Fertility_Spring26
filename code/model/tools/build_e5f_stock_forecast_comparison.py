"""Plot saved fertility-stock and housing observers; no model solve or fitting."""
from pathlib import Path
import csv
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/recovered_sequence'


def read(path):
    return json.loads(path.read_text())


def completed(record):
    """Use the same top-bin representative and age observer as the 2023 table."""
    flow = record['fertility']
    third = np.asarray(flow['birth_flow_third_bin_entry'])
    explicit = np.asarray(flow['birth_flow_explicit'])
    adjusted = np.asarray(flow['birth_flow_topcode_adjusted'])
    top = 3 + (adjusted - explicit)[third > 1e-12] / third[third > 1e-12]
    assert len(top) and np.allclose(top, top[0], rtol=0, atol=1e-12)
    stock = record['fertility_stock_timing']
    assert stock['metadata']['age_projection'] == 'uniform_birth_time'
    shares = stock['parity_shares_40_44']
    assert abs(sum(shares.values()) - 1) < 1e-12
    value = shares['1'] + 2 * shares['2'] + float(top[0]) * shares['3plus']
    mass = np.asarray(stock['accounting']['window_parity_mass'])
    direct = float(mass @ np.array([0, 1, 2, top[0]]) / mass.sum())
    assert abs(value - direct) < 1e-12
    return value


def main():
    source = BASE / 'source/stock_forecast'
    packet = read(source / 'observed_dates.json')
    check = read(source / 'verification.json')
    prior = read(BASE / 'source/readout_2023/model_2023.json')
    assert check['status'] == 'PASS' and check['replay_maximum_abs'] <= 2e-10
    assert packet['forecast_receipt_sha256'] == prior['forecast_receipt_sha256']
    records = packet['dates']
    years = [r['calendar_year'] for r in records]
    assert years == [2019, 2023, 2027, 2031, 2035, 2039]
    observed2023 = records[1]
    for key in ('profile', 'fertility', 'fertility_stock_timing'):
        assert observed2023[key] == prior[key], f'Original 2023 {key} differs'
    rows = []
    for r in records:
        totals = r['profile']['totals']
        rows.append(dict(year=r['calendar_year'], completed_fertility_40_44=completed(r),
                         capped_rooms=totals['capped_rooms'], households=totals['households'],
                         mean_capped_rooms=totals['capped_rooms'] / totals['households']))
    validation = {r['moment_key']: r for r in csv.DictReader((BASE / 'figures/validation_2023.csv').open())}
    assert abs(rows[1]['completed_fertility_40_44'] - float(validation['completed_fertility']['model'])) < 1e-12
    assert abs(rows[1]['mean_capped_rooms'] - float(validation['mean_rooms']['model'])) < 1e-12
    cps = {int(k): v for k, v in read(source / 'cps_history_40_44.json').items() if int(k) >= 2000}
    hp = ROOT / 'output/model/e5f_matched_pf_20260909a/design_research/housing/early_housing_target_candidates.csv'
    housing = {int(r['window']): float(r['point']) for r in csv.DictReader(hp.open())
               if r['window'] in ('2007', '2012', '2023') and r['moment'] == 'aggregate_mean_occupied_rooms_capped9_18_85'}
    period = read(BASE / 'figures/figure_verification.json')
    blue, red = '#1f5fa6', '#c73e3a'
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(15, 5.6))
    a.plot(period['years'], period['model'], 'o-', color=blue, markersize=4, label='Model history')
    a.plot(period['future_years'], period['future_model'], 'o:', color=blue, markersize=4, label='Model continuation')
    a.plot(period['years'][1:], period['data'], 's--', color=red, markersize=4, label='Data')
    a.set(title='Period fertility', xlabel='End of four-year birth window', ylabel='Births per woman', ylim=(1.45, 2.2))
    a.legend(frameon=False, fontsize=8.5, loc='upper right')
    for ax, key in ((b, 'completed_fertility_40_44'), (c, 'mean_capped_rooms')):
        values = [r[key] for r in rows]
        ax.plot(years[:2], values[:2], 'o-', color=blue, markersize=4, label='Model history')
        ax.plot(years[1:], values[1:], 'o:', color=blue, markersize=4, label='Model continuation')
        ax.annotate(f'{values[-1]:.2f}', (years[-1], values[-1]), xytext=(4, 9), textcoords='offset points', color=blue)
    x = sorted(cps)
    b.plot(x, [cps[y] for y in x], 's--', color=red, markersize=4, label='CPS')
    b.set(title='Completed fertility measure', xlabel='Survey / model year', ylabel='Children ever born, ages 40–44', ylim=(1.4, 2.2))
    x = sorted(housing)
    c.plot(x, [housing[y] for y in x], 's--', color=red, markersize=4, label='ACS')
    c.set(title='Housing per household', xlabel='Calendar year', ylabel='Occupied rooms, capped at nine', ylim=(5.2, 6.4))
    for ax in (a, b, c):
        ax.axvline(2023, color='.65', lw=.8, ls='--')
        ax.axvspan(2023, 2044, color=blue, alpha=.035)
        ax.set_xlim(1999 if ax is b else 2005, 2044)
        ax.set_xticks([2000, 2011, 2023, 2031, 2039] if ax is b else [2007, 2015, 2023, 2031, 2039])
        ax.grid(alpha=.16)
    fig.suptitle('Fertility and housing: historical comparison and continuation', fontsize=17, x=.06, ha='left', y=.96)
    fig.text(.06, .035,
             'Sources: Census CPS Historical Table 2; ACS, 42 metros; saved model path. CPS 2022/2024 counts are capped at five.\n'
             'Completed fertility uses ages 40–44 throughout; CPS 2024 is the nearest survey to the model’s 2023 observation. Model stock/housing series start in 2019.\n'
             'Continuation holds the last fertility preference fixed. No property-tax rebate; terminal-horizon sensitivity remains unresolved.',
             fontsize=8.3, color='#555555', linespacing=1.5)
    fig.subplots_adjust(left=.06, right=.97, bottom=.25, top=.81, wspace=.32)
    out = BASE / 'figures'
    for ext in ('png', 'pdf'):
        fig.savefig(out / f'historical_fit_stock_forecast.{ext}', dpi=170, facecolor='white')
    plt.close(fig)
    with (out / 'stock_forecast.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    (out / 'stock_forecast_verification.json').write_text(json.dumps(dict(
        status='PASS', original_2023_exact_match=True, records=rows,
        source_sha256=hashlib.sha256((source / 'observed_dates.json').read_bytes()).hexdigest(),
        replay_maximum_abs=check['replay_maximum_abs'], horizon_verified=False), indent=2)+'\n')
    print(json.dumps(rows))


if __name__ == '__main__':
    main()
