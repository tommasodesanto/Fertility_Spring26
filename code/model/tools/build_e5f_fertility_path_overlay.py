"""Overlay saved historical fertility and the permanent-shock diagnostic.

Display only: no model solve, fitted value change, extrapolation, or deck edit.
Both curves retain the start-of-window dating of the historical slide.
"""
from pathlib import Path
import csv
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
HISTORY = ROOT / 'output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/recovered_sequence'
CURRENT = ROOT / 'output/model/e5f_original_queue_20260913a/terminal_restart_v1/fertility_replay_iter3/output'


def read(path):
    return json.loads(path.read_text())


def main():
    inputs = [HISTORY / 'figures/stock_forecast.csv', HISTORY / 'historical_fit.csv',
              HISTORY / 'figures/figure_verification.json',
              CURRENT / 'fertility.json', CURRENT / 'rows.json',
              CURRENT / 'stationary_reference.json', CURRENT / 'lead_verification.json']
    with inputs[0].open() as stream:
        history = list(csv.DictReader(stream))
    with inputs[1].open() as stream:
        fit = list(csv.DictReader(stream))
    old_check = read(inputs[2])
    fertility, native, reference, check = [read(p) for p in inputs[3:]]
    assert check['status'] == 'passed' and check['iteration'] == 3
    assert reference['stationary_endpoint_verified']
    assert len(fertility) == len(native) == 100
    years = np.array([int(r['year']) for r in history])
    previous = np.array([float(r['period_fertility']) for r in history])
    assert years.tolist() == list(range(2007, 2040, 4))
    original = dict(zip([y - 4 for y in old_check['years'][1:]], old_check['model'][1:]))
    original.update(zip([y - 4 for y in old_check['future_years'][1:]], old_check['future_model'][1:]))
    np.testing.assert_allclose(previous, [original[int(y)] for y in years], rtol=0, atol=2e-10)
    current_by_year = {}
    for f, row in zip(fertility, native):
        assert f['calendar_year'] == row['calendar_year']
        mass = np.asarray(f['age_mass'])
        flow = np.asarray(f['birth_flow_topcode_adjusted'])
        value = float(np.divide(flow, mass, out=np.zeros_like(flow), where=mass > 0).sum())
        assert abs(value - f['period_tfr_topcode_adjusted']) < 2e-10
        assert abs(float(flow.sum()) - row['birth_children_topcode_adjusted']) < 2e-10
        current_by_year[f['calendar_year']] = value
    permanent = np.array([current_by_year[int(y)] for y in years])
    data_years = [int(r['year']) for r in fit]
    data = [float(r['target']) for r in fit]
    assert data_years == years[:4].tolist()
    np.testing.assert_allclose([float(r['model']) for r in fit], previous[:4], rtol=0, atol=2e-10)
    np.testing.assert_allclose(data, old_check['data'], rtol=0, atol=2e-10)
    initial = reference['initial']['fertility']['period_tfr_topcode_adjusted']
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(figsize=(10.4, 6.4))
    blue, red, orange = '#1f5fa6', '#c73e3a', '#d97815'
    model_line, = ax.plot(years[:4], previous[:4], 'o-', color=blue, lw=2,
                         ms=5, label='Prior fitted shock sequence')
    future_line, = ax.plot(years[3:], previous[3:], 'o:', color=blue, lw=2,
                          ms=5, label='Prior constant-preference continuation')
    data_line, = ax.plot(data_years, data, 's--', color=red, lw=1.5,
                        ms=5, label='Data: four-year averages', zorder=4)
    new_x = np.r_[2005, 2007, years]
    new_y = np.r_[initial, initial, permanent]
    new_line, = ax.plot(new_x, new_y, 'o-', color=orange, lw=2.3,
                       ms=5, label='One permanent shock: current diagnostic')
    ax.annotate('Initial steady state: 2.10', (2007, initial), xytext=(2010, 2.11),
                fontsize=10, color=orange)
    ax.annotate(f'{permanent[0]:.3f}', (2007, permanent[0]), xytext=(7, -16),
                textcoords='offset points', color=orange, fontsize=10)
    ax.axvline(2023, color='.65', lw=.8, ls='--')
    ax.axvspan(2023, 2041, color=blue, alpha=.025)
    ax.set(xlim=(2005, 2041), ylim=(1.53, 2.17),
           xlabel='Start of four-year birth window', ylabel='Births per woman')
    ax.set_xticks(years)
    ax.grid(alpha=.16)
    ax.legend(handles=[data_line, model_line, future_line, new_line],
              loc='upper right', bbox_to_anchor=(1, .89), frameon=False, fontsize=9.5)
    fig.suptitle('Fertility: fitted shock sequence versus one permanent shock',
                 fontsize=15, x=.085, ha='left', y=.96)
    fig.text(.085, .13,
             'Same date convention as the existing slide: 2007 labels the first four-year window (data: 2008–2011).',
             fontsize=8.8, color='#555555')
    fig.text(.085, .095,
             'Prior run: no rebate, demographic conditioning, six-period forecasts. Current: equal rebates, closed population, 100 periods.',
             fontsize=8.5, color='#555555')
    fig.text(.085, .06,
             'Current path is unconverged iteration 3; prior horizon is unverified. Differences do not isolate shock timing alone.',
             fontsize=8.8, color='#555555')
    fig.subplots_adjust(left=.085, right=.97, bottom=.23, top=.85)
    for artist, x, y in [(model_line, years[:4], previous[:4]),
                          (future_line, years[3:], previous[3:]),
                          (data_line, data_years, data), (new_line, new_x, new_y)]:
        np.testing.assert_array_equal(artist.get_xdata(), x)
        np.testing.assert_array_equal(artist.get_ydata(), y)
    for ext in ['png', 'pdf']:
        fig.savefig(CURRENT / f'fertility_overlay.{ext}', dpi=180, facecolor='white')
    plt.close(fig)
    with (CURRENT / 'fertility_overlay.csv').open('w') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(['window_start_year', 'prior_sequence_or_continuation', 'permanent_shock_iteration3', 'data'])
        for i, year in enumerate(years):
            writer.writerow([int(year), previous[i], permanent[i], data[i] if i < 4 else ''])
    verification = dict(status='passed', display_only=True, model_solves=0,
                        same_start_of_window_dating=True, identical_source_values=True,
                        artist_data_verified=True, current_iteration=3,
                        controlled_timing_comparison=False,
                        source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs})
    (CURRENT / 'fertility_overlay_verification.json').write_text(json.dumps(verification, indent=2) + '\n')
    print(CURRENT / 'fertility_overlay.png')


if __name__ == '__main__':
    main()
