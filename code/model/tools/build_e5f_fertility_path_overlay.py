"""Overlay saved historical fertility and the permanent-shock diagnostic.

Display only: no model solve, fitted value change, extrapolation, or deck edit.
Both curves retain the start-of-window dating of the historical slide.
"""
from pathlib import Path
import argparse
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


def plot_macro_transition(fertility, rows, reference, fit, inputs):
    """Saved-path mock figures; cohort stocks use the original queue identity.

    Survival depends only on age, entry is childless at age 18, and all other
    transitions preserve children ever born. Therefore the cohort conditional
    mean obeys C[t,j] = C[t-1,j-1] + births[t,j]/mass[t,j]. The initial stationary
    age rates supply pre-2007 history. These are reconstructed model stocks,
    with the same fixed 3+ weight as the native birth-flow diagnostics.
    """
    years=np.array([r['calendar_year'] for r in rows],dtype=int)
    rates=np.array([r['age_specific_birth_rate_topcode_adjusted'] for r in fertility])
    old_rates=np.array(reference['initial']['fertility']['age_specific_birth_rate_topcode_adjusted'])
    ages=np.array(fertility[0]['age_cell_start'])
    np.testing.assert_array_equal(years,2007+4*np.arange(100))
    np.testing.assert_array_equal(ages,18+4*np.arange(len(ages)))
    end_age=int(np.flatnonzero(ages==42)[0])
    assert np.all(rates[:,end_age+1:]==0) and np.all(old_rates[end_age+1:]==0)
    # Check rates against native flows and independently accumulate the stock
    # by recurrence and by cohort-diagonal summation.
    for d,row,rate in zip(fertility,rows,rates):
        masses=np.array(d['age_mass']);flows=np.array(d['birth_flow_topcode_adjusted'])
        np.testing.assert_allclose(rate,np.divide(flows,masses,out=np.zeros_like(flows),where=masses>1e-15),atol=2e-12,rtol=0)
        assert abs(flows.sum()-row['birth_children_topcode_adjusted'])<2e-10
    stocks=[]; previous=np.cumsum(old_rates)
    for rate in rates:
        current=np.r_[0.,previous[:-1]]+rate
        stocks.append(current); previous=current
    completed=np.array(stocks)[:,end_age]
    diagonal=np.array([sum((old_rates[j] if t-end_age+j<0 else rates[t-end_age+j,j])
        for j in range(end_age+1)) for t in range(len(rows))])
    np.testing.assert_allclose(completed,diagonal,atol=2e-12,rtol=0)
    period=rates.sum(axis=1); initial_tfr=float(old_rates.sum())
    initial=reference['initial']['quantities']; terminal=reference['terminal']['quantities']
    housing=np.array([r['housing_demand'] for r in rows])/initial['housing_demand']*100
    supply=np.array([r['housing_supply'] for r in rows])/initial['housing_demand']*100
    households=np.array([r['adult_population'] for r in rows])/initial['adult_population']*100
    per_household=housing/households*100
    terminal_housing=100*terminal['housing_demand']/initial['housing_demand']
    terminal_households=100*terminal['adult_population']/initial['adult_population']
    data_years=np.array([int(r['year']) for r in fit]); data=np.array([float(r['target']) for r in fit])
    plt.rcParams.update({'font.size':13,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.titlesize':15,'axes.labelsize':12})
    orange,blue,grey='#d97815','#245f99','#777777'
    plotted={}
    for last,stem in ((2103,'macro_transition_mock'),(2403,'macro_transition_full_horizon')):
        keep=years<=last; x=np.r_[2003,2007,years[keep]]
        fig,axes=plt.subplots(2,2,figsize=(13.4,7.8))
        def line(ax,values,pre,label,color,style='-'):
            y=np.r_[pre,pre,values[keep]]
            artist,=ax.plot(x,y,color=color,lw=2.3,ls=style,label=label)
            np.testing.assert_array_equal(artist.get_xdata(),x)
            np.testing.assert_array_equal(artist.get_ydata(),y)
            return artist
        ax=axes[0,0]
        line(ax,period,initial_tfr,'Model',orange)
        artist,=ax.plot(data_years,data,'s--',color=blue,lw=1.6,ms=5,label='US data')
        np.testing.assert_array_equal(artist.get_ydata(),data)
        ax.axhline(initial_tfr,color=grey,lw=.9,ls=':')
        ax.set(title='Period fertility',ylabel='Births per woman',ylim=(1.55,2.16));ax.legend(frameon=False,fontsize=11)
        ax=axes[0,1]
        line(ax,completed,initial_tfr,'Model cohort',orange)
        ax.axhline(initial_tfr,color=grey,lw=.9,ls=':')
        ax.set(title='Completed fertility',ylabel='Children by the end of fertile ages',ylim=(1.55,2.16))
        ax=axes[1,0]
        line(ax,housing,100,'Housing services used',orange)
        line(ax,supply,100,'Housing supplied',grey,':')
        line(ax,households,100,'Households',blue)
        ax.set(title='Aggregate housing and households',ylabel='Initial steady state = 100')
        ax.legend(frameon=False,fontsize=10.5,loc='lower left')
        ax=axes[1,1]
        line(ax,per_household,100,'Housing per household',orange)
        ax.set(title='Housing services per household',ylabel='Initial steady state = 100')
        for ax in axes.flat:
            ax.axvline(2023,color='.75',lw=.8,ls='--')
            ax.set_xlim(2003,last+2);ax.grid(axis='y',alpha=.16)
            ax.set_xlabel('Start of four-year period')
            ax.set_xticks([2007,2023,2043,2063,2083,2103] if last==2103 else [2007,2103,2203,2303,2403])
            ax.tick_params(axis='x',labelsize=10.5)
        fig.suptitle('Demographic adjustment after a permanent preference decline',fontsize=17,y=.97)
        fig.text(.07,.028,'Preference falls 38.1% in 2007 and stays fixed. Numerical diagnostic: transition not converged.',fontsize=10,color='.3')
        fig.subplots_adjust(left=.075,right=.98,top=.87,bottom=.13,wspace=.25,hspace=.48)
        for ext in ('png','pdf'):fig.savefig(CURRENT/f'{stem}.{ext}',dpi=180,facecolor='white')
        plt.close(fig)
        plotted[stem]=dict(first_year=2003,last_year=last)
    table=[]
    for i,y in enumerate(years):
        table.append(dict(window_start_year=int(y),period_fertility=float(period[i]),
            completed_fertility_reconstructed=float(completed[i]),housing_index=float(housing[i]),
            supply_index=float(supply[i]),households_index=float(households[i]),
            housing_per_household_index=float(per_household[i]),relative_market_residual=rows[i]['relative_market_residual']))
    with (CURRENT/'macro_transition_mock.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(table[0]),lineterminator='\n');writer.writeheader();writer.writerows(table)
    verification=dict(status='passed',model_solves=0,current_iteration=3,finite_equilibrium_converged=False,
        artist_data_verified=True,cohort_stock_identity_verified=True,
        completed_fertility_definition='Reconstructed top-code-adjusted children ever born after births in model age cell42–45; stationary prehistory, age-only survival and childless entry. Not a direct saved-distribution measurement or an empirical ages40–44 match.',
        data_definition='Four-year fertility averages;2007 denotes2008–2011,2019 denotes2020–2023.',
        housing_definition='Aggregate services demanded and supplied shown separately because markets are not yet cleared.',
        terminal_reference=dict(housing_index=terminal_housing,households_index=terminal_households,
             housing_per_household_index=terminal_housing/terminal_households*100),
        completed_minimum=dict(year=int(years[completed.argmin()]),value=float(completed.min())),
        selected_rows=[table[i] for i in (0,4,14,24,99)],plots=plotted,
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs})
    (CURRENT/'macro_transition_verification.json').write_text(json.dumps(verification,indent=2)+'\n')
    print(json.dumps(verification['selected_rows'],indent=2))


def plot_data_model(years, permanent, data_years, data, initial, inputs):
    """Requested two-series view, preserving all underlying observations."""
    plt.rcParams.update({'font.size': 12, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(figsize=(9, 5.5))
    x, y = np.r_[2005, 2007, years], np.r_[initial, initial, permanent]
    model_line, = ax.plot(x, y, 'o-', color='#d97815', lw=2.3, ms=5, label='Model')
    data_line, = ax.plot(data_years, data, 's--', color='#c73e3a', lw=2, ms=5, label='Data')
    ax.set(title='Fertility: data and model', xlabel='Start of four-year period',
           ylabel='Births per woman', xlim=(2005, 2041), ylim=(1.55, 2.15))
    ax.set_xticks(years)
    ax.grid(alpha=.16)
    ax.legend(handles=[data_line, model_line], frameon=False, loc='upper right')
    fig.text(.11, .035, 'Data: four-year averages. Model: permanent shock, unconverged.',
             fontsize=9, color='#555555')
    fig.subplots_adjust(left=.11, right=.97, top=.88, bottom=.19)
    for artist, expected_x, expected_y in [(model_line, x, y), (data_line, data_years, data)]:
        np.testing.assert_array_equal(artist.get_xdata(), expected_x)
        np.testing.assert_array_equal(artist.get_ydata(), expected_y)
    for ext in ['png', 'pdf']:
        fig.savefig(CURRENT / f'fertility_data_model.{ext}', dpi=180, facecolor='white')
    plt.close(fig)
    with (CURRENT / 'fertility_data_model.csv').open('w') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(['window_start_year', 'model', 'data'])
        for i, year in enumerate(years):
            writer.writerow([int(year), permanent[i], data[i] if i < len(data) else ''])
    (CURRENT / 'fertility_data_model_verification.json').write_text(json.dumps(dict(
        status='passed', artist_data_verified=True, current_iteration=3,
        model_series='permanent_shock', displayed_series=['Data', 'Model'],
        model_solves=0, same_start_of_window_dating=True,
        source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    ), indent=2) + '\n')
    print(CURRENT / 'fertility_data_model.png')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-model-only', action='store_true')
    parser.add_argument('--macro-transition', action='store_true')
    args = parser.parse_args()
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
    if args.macro_transition:
        plot_macro_transition(fertility,native,reference,fit,inputs)
        return
    if args.data_model_only:
        plot_data_model(years, permanent, data_years, data, initial, inputs)
        return
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
