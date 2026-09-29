"""Torch-only saved-table comparison of the verified two-birth experiment.

No model imports, checkpoint reads, equilibrium solves, or standard-plot writes.
Run only after run_v1/completion.json confirms successful final export.
The optional --plot writes one supplemental lifecycle PNG in saved_analysis/.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys


COORDINATES = (
    'H0', 'beta_annual', 'chi', 'first_birth_fixed_cost', 'kappa_fert',
    'kappa_fert_continuation', 'theta0', 'delta_alpha_jump',
    'child_benefit_curvature', 'tenure_choice_kappa',
)
HELD_STATUS = 'held at reference estimate during two-birth diagnostic'
CASE_STATUS = 'verified_experimental_two_birth_equilibrium_not_adopted'
LABEL = '2007 stationary reference — block0506, September 28 verified export'


def read(path):
    return json.loads(Path(path).read_text())


def rows(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def close(left, right, label, atol=1e-12):
    assert math.isfinite(float(left)) and math.isfinite(float(right)), label
    assert abs(float(left)-float(right)) <= atol, (label, left, right)


def main():
    assert sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Torch Slurm only'
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
        os.environ[key] = '1'
    import numpy as np

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plot', action='store_true', help='Write one supplemental three-panel PNG')
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    base = here.parent
    reference = base/'resume_v1/selected_export/primary'
    run = here/'run_v1'
    worker = run/'worker'
    case = worker/'case'
    out = run/'saved_analysis'
    audit = base/'measurement_audit_v1'
    sources = {Path(__file__).resolve()}

    def load(path):
        sources.add(path)
        return read(path)

    completion = load(run/'completion.json')
    assert completion['status'] == 'completed' and completion['exit_code'] == 0
    worker_receipt = load(worker/'receipt.json')
    assert worker_receipt['status'] == 'passed' and worker_receipt['scientific_promotion'] is False
    assert Path(worker_receipt['case']).resolve() == case.resolve()
    receipts = {}
    standard_plot_names = {}
    for name, path in (('reference', reference), ('two_birth', case)):
        hashes = load(path/'artifact_hashes.json')
        for filename in ('target_fit.csv', 'parameters.csv', 'observers.json', 'receipt.json'):
            source = path/filename
            sources.add(source)
            assert sha(source) == hashes[filename], (name, filename)
        standard_plot_names[name] = sorted(filename for filename in hashes
                                           if filename.startswith('standard_diagnostics/') and filename.endswith('.png'))
        assert len(standard_plot_names[name]) == 17, (name, 'standard PNG count')
        for filename in standard_plot_names[name]:
            source = path/filename
            sources.add(source)
            assert sha(source) == hashes[filename], (name, filename)
        receipts[name] = read(path/'receipt.json')
    assert standard_plot_names['reference'] == standard_plot_names['two_birth']
    receipt = receipts['two_birth']
    assert receipt['status'] == CASE_STATUS and receipt['scientific_promotion'] is False
    assert receipt['reference_checkpoint_sha256'] == receipts['reference']['case_checkpoint_sha256']
    assert receipt['held_reference_coordinate_count'] == receipt['reference_calibration_coordinate_count'] == 10
    assert receipt['estimated_coordinates_this_experiment'] == receipt['free_count'] == 0
    assert receipt['normalized_coordinates'] == 1
    assert receipt['extra_cache_audit']['status'] == 'passed'
    close(worker_receipt['loss'], receipt['loss'], 'worker loss', 1e-10)
    manifest_path = worker/'effective_source_manifest.json'
    manifest = load(manifest_path)
    assert sha(manifest_path) == receipt['effective_source_manifest_sha256']
    for entry in manifest.values():
        path = Path(entry['generated_source'])
        assert path.resolve().is_relative_to((worker/'effective_sources').resolve())
        sources.add(path)
        assert sha(path) == entry['effective_sha256'], str(path)
    for filename, digest in worker_receipt['source_pins'].items():
        path = here/filename
        sources.add(path)
        assert sha(path) == digest, filename

    raw_fits = {name: rows(path/'target_fit.csv')
                for name, path in (('reference', reference), ('two_birth', case))}
    assert all(len(value) == 14 for value in raw_fits.values())
    fits = {name: {r['moment']: r for r in value} for name, value in raw_fits.items()}
    assert len(fits['reference']) == len(fits['two_birth']) == 14
    assert set(fits['reference']) == set(fits['two_birth'])
    fit_comparison = []
    for moment, old in fits['reference'].items():
        new = fits['two_birth'][moment]
        for field in ('target', 'weight', 'role'):
            assert old[field] == new[field], (moment, field)
        entry = dict(moment=moment, role=old['role'], target=old['target'], weight=old['weight'])
        for name, row in (('reference', old), ('two_birth', new)):
            close(float(row['model'])-float(row['target']), row['gap'], (moment, name, 'gap'))
            if row['loss_contribution']:
                close(float(row['weight'])*float(row['gap'])**2, row['loss_contribution'],
                      (moment, name, 'loss'), 1e-9)
            entry.update({name+'_'+key: row[key] for key in ('model', 'gap', 'loss_contribution')})
        entry['two_birth_minus_reference'] = float(new['model'])-float(old['model'])
        fit_comparison.append(entry)
    assert [sum(r['role'] == role for r in fits['reference'].values())
            for role in ('scored', 'validation', 'normalization')] == [10, 3, 1]
    for name in fits:
        close(sum(float(r['loss_contribution'] or 0) for r in fits[name].values()),
              receipts[name]['loss'], (name, 'complete loss'), 1e-9)

    raw_parameters = {name: rows(path/'parameters.csv')
                      for name, path in (('reference', reference), ('two_birth', case))}
    assert all(len(value) == 31 for value in raw_parameters.values())
    parameters = {name: {r['parameter']: r for r in value} for name, value in raw_parameters.items()}
    assert len(parameters['reference']) == len(parameters['two_birth']) == 31
    assert set(parameters['reference']) == set(parameters['two_birth'])
    assert {p for p, r in parameters['reference'].items() if r['lower'] and r['upper']} == set(COORDINATES)
    parameter_comparison = []
    for parameter, old in parameters['reference'].items():
        new = parameters['two_birth'][parameter]
        for field in ('lower', 'upper', 'near_bound'):
            assert old[field] == new[field], (parameter, field)
        if parameter in COORDINATES:
            assert float(old['estimate']) == float(new['estimate']), parameter
            assert new['status'] == HELD_STATUS, parameter
            assert float(new['lower']) <= float(new['estimate']) <= float(new['upper']), parameter
        else:
            assert new['status'] == old['status'], (parameter, 'restriction status')
        parameter_comparison.append(dict(parameter=parameter,
            reference_estimate=old['estimate'], two_birth_estimate=new['estimate'],
            two_birth_minus_reference=float(new['estimate'])-float(old['estimate']),
            lower=new['lower'], upper=new['upper'], near_bound=new['near_bound'],
            reference_status=old['status'], two_birth_status=new['status'],
            held_reference_coordinate=parameter in COORDINATES))

    observers = {name: read(path/'observers.json')['fertility']['uniform_birth_time']
                 for name, path in (('reference', reference), ('two_birth', case))}
    projection_note = observers['two_birth']['metadata']['early_fertility_age25']['two_birth_projection_convention']
    assert projection_note and 'common-event pre/post interpolation' in projection_note

    def project(observer, lo, hi):
        a = observer['accounting']
        ages = np.asarray(a['age_cell_start'], dtype=float)
        pre = np.asarray(a['pre_parity_mass_by_age'], dtype=float)
        post = np.asarray(a['post_parity_mass_by_age'], dtype=float)
        assert pre.shape == post.shape == (len(ages), 4)
        assert np.all(np.isfinite(pre)) and np.all(np.isfinite(post))
        assert np.all(pre >= 0) and np.all(post >= 0)
        assert np.array_equal(ages, 18+4*np.arange(len(ages)))
        assert np.max(np.abs(pre.sum(axis=1)-post.sum(axis=1))) < 2e-10
        left = np.maximum(ages, lo); right = np.minimum(ages+4, hi+1)
        overlap = np.maximum(right-left, 0)/4
        close(overlap.sum()*4, hi+1-lo, 'age-window coverage')
        alpha = ((left+right)/2-ages)/4
        mass = (overlap[:, None]*((1-alpha[:, None])*pre+alpha[:, None]*post)).sum(axis=0)
        assert mass.sum() > 0 and np.all(mass >= 0)
        shares = mass/mass.sum()
        capped = float(shares@np.arange(4)); mother = float(1-shares[0])
        assert mother > 0
        return dict(capped3=capped, mother_share=mother, given_mother=capped/mother,
                    share_0=float(shares[0]), share_1=float(shares[1]),
                    share_2=float(shares[2]), share_3plus=float(shares[3]))

    lifecycle_source = audit/'fertility_lifecycle_matched_windows.csv'
    sources.add(lifecycle_source)
    saved_profiles = rows(lifecycle_source)
    assert [(int(r['age_lower']), int(r['age_upper'])) for r in saved_profiles] == [(20,24),(25,29),(30,34),(35,39),(40,44)]
    lifecycle = []
    for saved in saved_profiles:
        lo, hi = int(saved['age_lower']), int(saved['age_upper'])
        entry = dict(age_lower=lo, age_upper=hi, role='untargeted cross-sectional lifecycle diagnostic')
        for field in ('capped3', 'mother_share', 'given_mother'):
            entry['data_'+field] = float(saved['data_'+field])
        entry['data_share_3plus'] = float(saved['data_3plus'])
        for name, observer in observers.items():
            projected = project(observer, lo, hi)
            if name == 'reference':
                for field in ('capped3', 'mother_share', 'given_mother'):
                    close(projected[field], saved['model_'+field], ('reference saved lifecycle', lo, field))
                close(projected['share_3plus'], saved['model_3plus'], ('reference 3plus', lo))
            entry.update({name+'_'+field: value for field, value in projected.items()})
        lifecycle.append(entry)

    saved25 = load(audit/'early_fertility_decomposition.json')
    age25 = {'data': dict(capped3=saved25['data_early_fertility'],
                         mother_share=saved25['data_mother_share'],
                         given_mother=saved25['data_capped_children_given_mother'])}
    for name, observer in observers.items():
        age25[name] = project(observer, 25, 25)
        close(age25[name]['capped3'], fits[name]['early_fertility']['model'], (name, 'age25 target replay'))
        for i, key in enumerate(('0', '1', '2', '3plus')):
            close(age25[name]['share_'+key], observer['ever_born_shares_age25'][key], (name, 'age25 count shares', i))
        stock = project(observer, 40, 44)
        close(1-stock['mother_share'], fits[name]['cps_childlessness']['model'], (name, 'childlessness replay'))
        close(stock['share_1']/stock['mother_share'], fits[name]['cps_exactly_one']['model'], (name, 'exactly one replay'))
    for field, key in (('capped3', 'model_early_fertility'), ('mother_share', 'model_mother_share'),
                       ('given_mother', 'model_capped_children_given_mother')):
        close(age25['reference'][field], saved25[key], ('reference age25 saved audit', field))

    decompositions = []
    for left, right in (('reference', 'data'), ('two_birth', 'data'), ('two_birth', 'reference')):
        x, y = age25[left], age25[right]
        extensive = (x['mother_share']-y['mother_share'])*(x['given_mother']+y['given_mother'])/2
        intensive = (x['given_mother']-y['given_mother'])*(x['mother_share']+y['mother_share'])/2
        difference = x['capped3']-y['capped3']
        close(extensive+intensive, difference, (left, right, 'symmetric identity'))
        decompositions.append(dict(comparison=left+'_minus_'+right,
            capped_children_difference=difference, symmetric_motherhood_component=extensive,
            symmetric_conditional_children_component=intensive,
            conditional_children_fraction=intensive/difference if abs(difference)>1e-14 else None,
            interpretation='Arithmetic symmetric decomposition; not causal attribution or a target-reachability result'))
    close(decompositions[0]['symmetric_motherhood_component'], saved25['symmetric_mother_share_component'], 'reference extensive replay')
    close(decompositions[0]['symmetric_conditional_children_component'], saved25['symmetric_children_given_mother_component'], 'reference intensive replay')

    first_source = audit/'first_birth_age_cells.csv'
    sources.add(first_source)
    saved_first = rows(first_source)
    assert len(saved_first) == 7
    first_cells = []
    for i, saved in enumerate(saved_first):
        lo, hi = int(saved['age_lower']), int(saved['age_upper'])
        assert (lo, hi) == (18+4*i, 21+4*i)
        entry = dict(age_lower=lo, age_upper=hi, data_share=float(saved['empirical_share']),
                     empirical_tail_mapping=saved['empirical_tail_mapping'], role='untargeted individual first-birth share')
        for name, observer in observers.items():
            a = observer['accounting']
            assert a['age_cell_start'][i] == lo and a['first_birth_flow'] > 0
            entry[name+'_share'] = float(a['parity_birth_flows_by_age'][i][0]/a['first_birth_flow'])
        close(entry['reference_share'], saved['model_share'], ('reference first-birth share', i))
        first_cells.append(entry)
    timing = {}
    for name in ('data', 'reference', 'two_birth'):
        values = [r[name+'_share'] for r in first_cells]
        assert all(0 <= value <= 1 for value in values)
        close(sum(values), 1, (name, 'first-birth share sum'))
        mean = sum((r['age_lower']+2)*r[name+'_share'] for r in first_cells)
        late = sum(r[name+'_share'] for r in first_cells if r['age_lower'] >= 30)
        fit = fits['reference' if name == 'data' else name]
        field = 'target' if name == 'data' else 'model'
        close(mean, fit['nchs_mean_age'][field], (name, 'mapped first-birth mean'), 1e-10)
        close(late, fit['nchs_share30'][field], (name, 'first-birth share30'), 1e-10)
        timing[name] = dict(mapped_mean_first_birth_age=mean, share_first_births_age30plus=late)

    out.mkdir(parents=True, exist_ok=True)
    outputs = []

    def table(filename, records):
        path = out/filename
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]), lineterminator='\n')
            writer.writeheader(); writer.writerows(records)
        outputs.append(path)

    table('target_fit_comparison.csv', fit_comparison)
    table('parameters_comparison.csv', parameter_comparison)
    table('fertility_lifecycle_comparison.csv', lifecycle)
    table('first_birth_age_cells_comparison.csv', first_cells)
    table('age25_comparison.csv', [dict(series=name, **{key: row[key] for key in
                                  ('capped3', 'mother_share', 'given_mother')}) for name, row in age25.items()])
    table('age25_symmetric_decomposition.csv', decompositions)
    if args.plot:
        os.environ['MPLBACKEND'] = 'Agg'
        import matplotlib.pyplot as plt
        from matplotlib.ticker import PercentFormatter
        plt.rcParams.update({'font.size':12, 'font.family':'DejaVu Sans',
                             'axes.spines.top':False, 'axes.spines.right':False})
        styles = [('data', 'Data: CPS 2004/2006', '#b65e27', 's--'),
                  ('reference', 'Frozen 2007 reference', '#126c8c', 'o-'),
                  ('two_birth', 'Two-birth experiment (not adopted)', '#31916a', '^:')]
        fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.5))
        x = np.arange(len(lifecycle))
        labels = [str(r['age_lower'])+'–'+str(r['age_upper']) for r in lifecycle]
        for ax, field, title, minimum_top in zip(axes, ('capped3','mother_share','given_mother'),
                ('Children per woman','Women who are mothers','Children among mothers'), (2.,1.,2.5)):
            maximum = minimum_top
            for name, label, color, style in styles:
                values = np.array([r[name+'_'+field] for r in lifecycle])
                line, = ax.plot(x, values, style, color=color, lw=2.3, ms=6, label=label)
                np.testing.assert_array_equal(line.get_ydata(), values)
                maximum = max(maximum, float(values.max())*1.05)
            ax.set(title=title, ylim=(0,maximum), xlabel='Age group', xticks=x, xticklabels=labels)
            ax.grid(axis='y', alpha=.2)
            if field == 'mother_share': ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
            else: ax.set_ylabel('Children ever born, capped at 3')
        fig.suptitle('Fertility over the lifecycle: model versus data', x=.075, ha='left', y=.985, fontsize=20)
        fig.text(.075,.914,'2007 stationary approximation; two-birth experiment at the ten reference calibration coordinates', fontsize=11,color='#555555')
        fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=3,bbox_to_anchor=(.5,.09),frameon=False,fontsize=10)
        fig.text(.5,.032,'Identical five-year windows; retained linear pre/post interpolation with a common event-time proxy.\n'
                 'Cross-sectional profiles, not a cohort trajectory or a 2007–2023 transition.',ha='center',fontsize=10,color='#555555')
        fig.subplots_adjust(left=.075,right=.985,top=.81,bottom=.25,wspace=.31)
        path = out/'supplemental_lifecycle.png'
        fig.savefig(path,dpi=160,facecolor='white');plt.close(fig);outputs.append(path)

    summary = dict(status='verified_saved_table_comparison_not_adopted', model_solves=0,
        model_imports=0, checkpoint_reads=0, slurm_job=os.environ['SLURM_JOB_ID'],
        reference_label=LABEL, case_status=receipt['status'], reference_checkpoint_sha256=receipts['reference']['case_checkpoint_sha256'],
        experimental_checkpoint_sha256=receipt['case_checkpoint_sha256'],
        effective_source_manifest_sha256=receipt['effective_source_manifest_sha256'],
        all_ten_reference_coordinates_unchanged=True, held_coordinate_statuses_verified=True,
        all_14_targets_weights_roles_unchanged=True, complete_parameter_rows=31,
        reference_lifecycle_first_birth_and_age25_replay=True, two_birth_projection=projection_note,
        age25=age25, age25_symmetric_decompositions=decompositions, first_birth_timing=timing,
        primary_losses={name: receipts[name]['loss'] for name in receipts},
        economic_changes=receipt['economic_changes'], standard_17_plots_changed=False,
        standard_plot_names=standard_plot_names['reference'], standard_plot_hashes_verified=True,
        standard_plot_scope=receipt['fertility_plot_scope'], transition_run=False, scientific_promotion=False,
        caveats=['CPS observations are pooled cross-sectional 2004/2006 women; model household reproductive-member exposure is approximate.',
                 'First-birth timing uses NCHS period counts 2003–2006 with boundary-tail collapse.',
                 'Two births share the retained within-cell event-time proxy; separate ordered birth dates are unavailable.',
                 'A fixed-coordinate economic experiment is not a recalibrated specification or a causal decomposition.'],
        input_hashes={str(path):sha(path) for path in sorted(sources)},
        output_hashes={path.name:sha(path) for path in outputs})
    (out/'analysis_receipt.json').write_text(json.dumps(summary,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps(dict(status=summary['status'], output=str(out), model_solves=0,
                         target_rows=14, parameter_rows=31, held_coordinates=10, plot=args.plot)))


if __name__ == '__main__':
    main()
