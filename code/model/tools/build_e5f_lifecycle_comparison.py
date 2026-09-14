"""Compare saved 2007 stationary and 2023 transition age profiles; no solves.

Model extraction runs beside the frozen runtime on Torch. Data extraction and
rendering run locally. The existing 2023 observer and ACS builder are reused.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[key] = '1'

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = ROOT / 'output/model/e5f_original_queue_20260913a/lifecycle_comparison'
PROFILE2023 = ROOT / 'output/model/e5f_original_queue_20260913a/terminal_restart_v1/fertility_replay_iter3/profile_2023'
DATA2023 = ROOT / 'output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/patch_readout/data/actual2023_age_housing_levels.csv'


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2,
        default=lambda x: x.tolist() if hasattr(x, 'tolist') else str(x)) + '\n')


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def extract_model(spec_path, helper_path, out):
    import numpy as np
    observer = module(helper_path, 'saved_profile_observer')
    spec = read(spec_path)
    sys.path.insert(0, str(Path(spec['batch']) / 'source'))
    import run_e5f_original_queue_experiments as runner
    c = runner.load_context(str(spec_path))
    e, P = c.packet['evaluation'], c.old.parameters
    np.testing.assert_allclose(e.policy.price, c.old.policy.price, rtol=0, atol=2e-10)
    np.testing.assert_allclose(e.g_pre, c.old.initial_state.g_pre, rtol=0, atol=2e-10)
    profile = observer.profile(e, P, c.old.b_grid)
    assert len(profile['rows']) == P.J == 17
    assert all(r['households'] > 0 for r in profile['rows'])
    save(out / 'model_2007.json', dict(calendar_year=2007, interpretation='Initial stationary economy before preference shocks',
        profile=profile, psi_child=float(P.psi_child), model_solves=0,
        source_spec=str(spec_path), source_spec_sha256=sha(spec_path),
        profile_observer=str(helper_path), profile_observer_sha256=sha(helper_path),
        source_checkpoint=c.base.get('initial_checkpoint', c.plan.get('initial_checkpoint')),
        verified_same_initial_distribution=True, verified_same_prices=True,
        verified_age_housing_aggregates=True))
    print(json.dumps({'model_profiles': len(profile['rows']), 'model_solves': 0}))


def extract_data(out):
    import numpy as np
    builder_path = ROOT / 'code/model/tools/extract_e5f_patch_validation.py'
    builder = module(builder_path, 'acs_profile_builder')
    # Same data pipeline as the 2023 figure; only the source year changes.
    builder.YEAR = 2007
    frame, pooled, _, meta = builder.extract()
    assert meta['year'] == 2007 and meta['year_row_range'][1] > meta['year_row_range'][0]
    assert len(frame) == 17 and np.all(frame.hhwt > 0)
    assert meta['room_code_audit']['records_ge99'] == 0
    assert np.all((frame.ownership_rate >= 0) & (frame.ownership_rate <= 1))
    assert np.all((frame.with_minor_rate >= 0) & (frame.with_minor_rate <= 1))
    # The legacy metadata contains literal 2023 labels unrelated to selection.
    # Do not carry those labels into the new household-only artifact.
    meta.pop('fertility_2023_source')
    meta['sample'] = 200701
    meta['builder'] = str(builder_path)
    meta['builder_sha256'] = sha(builder_path)
    meta['comparison_year'] = 2023
    out.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out / 'actual2007_age_housing_levels.csv', index=False)
    pooled.to_csv(out / 'actual2007_aggregate_validation.csv', index=False)
    save(out / 'actual2007_metadata.json', meta)
    print(json.dumps({'data_year': 2007, 'age_cells': len(frame),
        'households': int(frame.households.sum()), 'year_rows_read': meta['year_row_range']}))


def render(out):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    models = {2007: read(out / 'model_2007.json'), 2023: read(PROFILE2023 / 'model_2023.json')}
    check = read(PROFILE2023 / 'verification.json')
    assert check['status'] == 'PASS' and check['row_reproduction']['passed']
    assert check['model_source_sha256'] == sha(PROFILE2023 / 'model_2023.json')
    data_paths = {2007: out / 'actual2007_age_housing_levels.csv', 2023: DATA2023}
    panels = [('Homeownership', 'Percent of households', 'owners', 'ownership_rate', 100, (0, 100)),
              ('Housing size', 'Rooms, capped at 9', 'capped_rooms', 'mean_capped_rooms', 1, (3, 8)),
              ('Children at home', 'Percent of households', 'with_children', 'with_minor_rate', 100, (0, 80))]
    arrays = {}; table = []
    for year in (2007, 2023):
        age = models[year]['profile']['rows']
        with data_paths[year].open() as stream:
            data = {int(r['age_lower']): r for r in csv.DictReader(stream)}
        assert models[year]['calendar_year'] == year
        x = np.array([r['age'] + (r['age_width'] - 1) / 2 for r in age])
        arrays[year] = {'ages': x.tolist(), 'series': {}}
        for title, ylabel, mk, dk, scale, limits in panels:
            mv = np.array([scale * r[mk] / r['households'] for r in age])
            dv = np.array([scale * float(data[int(r['age'])][dk]) for r in age])
            assert np.isfinite(mv).all() and np.isfinite(dv).all()
            assert min(mv.min(), dv.min()) >= limits[0] and max(mv.max(), dv.max()) <= limits[1]
            arrays[year]['series'][mk] = {'model': mv.tolist(), 'data': dv.tolist()}
            for a, m, d in zip(age, mv, dv):
                table.append(dict(year=year, age_lower=a['age'], age_upper=a['age']+a['age_width']-1,
                    measure=mk, model=m, data=d, gap=m-d))
    prior = read(PROFILE2023 / 'figures/lifecycle_2023_verification.json')['plotted_arrays']
    assert arrays[2023] == prior, '2023 arrays must reproduce the existing lifecycle slide exactly'
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.titlesize': 13, 'legend.fontsize': 10})
    for years, stem in [((2007,), 'lifecycle_2007'), ((2007, 2023), 'lifecycle_2007_2023')]:
        fig, axes = plt.subplots(len(years), 3, figsize=(13, 4 * len(years)), squeeze=False)
        for row, year in enumerate(years):
            for ax, (title, ylabel, mk, _, _, limits) in zip(axes[row], panels):
                series = arrays[year]['series'][mk]; x = arrays[year]['ages']
                for key, color, style, label in [('model', '#1f5fa6', '-', 'Model'), ('data', '#c73e3a', '--', f'ACS {year}')]:
                    line, = ax.plot(x, series[key], color=color, ls=style, lw=2.2, label=label)
                    np.testing.assert_array_equal(line.get_ydata(), series[key])
                ax.set(title=f'{year}: {title}', ylabel=ylabel, xlabel='Age of household head',
                       xticks=[20, 35, 50, 65, 80], ylim=limits)
                ax.grid(alpha=.16)
            axes[row, 0].legend(frameon=False)
        fig.tight_layout(rect=(0, .065 if len(years)==1 else .04, 1, 1))
        fig.text(.5, .015, 'Children at home: model dependents; ACS resident own children under 18.', ha='center', fontsize=9)
        for ext in ('pdf', 'png'):
            fig.savefig(out / f'{stem}.{ext}', dpi=170)
        plt.close(fig)
    with (out / 'lifecycle_comparison.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(table[0]), lineterminator='\n'); writer.writeheader(); writer.writerows(table)
    inputs = [out/'model_2007.json', out/'actual2007_age_housing_levels.csv', out/'actual2007_metadata.json',
              PROFILE2023/'model_2023.json', PROFILE2023/'verification.json', DATA2023]
    save(out/'verification.json', dict(status='PASS', model_solves=0, exact_2023_figure_reproduction=True,
        same_axes=True, plotted_arrays=arrays, source_sha256={str(p):sha(p) for p in inputs},
        interpretation='2007 stationary calibration benchmark; 2023 one-permanent-shock transition cross-section, not converged. Full age curves are descriptive model/data comparisons, not individually targeted moments.'))
    print(json.dumps({'status':'PASS', 'figures':['lifecycle_2007','lifecycle_2007_2023'], 'model_solves':0}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['model', 'data', 'render'], required=True)
    parser.add_argument('--outdir', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--spec', type=Path)
    parser.add_argument('--profile-helper', type=Path)
    args = parser.parse_args()
    if args.mode == 'model':
        if args.spec is None or args.profile_helper is None:
            parser.error('--spec and --profile-helper are required for model extraction')
        extract_model(args.spec, args.profile_helper, args.outdir)
    elif args.mode == 'data': extract_data(args.outdir)
    else: render(args.outdir)


if __name__ == '__main__':
    main()
