"""Torch-only saved-table comparison of a fixed-benefit mechanism control.

Reads no checkpoints and imports no model code. The fixed-benefit control is
not a demographic steady state or a calibration candidate. It deliberately
reports replacement-fertility and renewal misses. All 17 standard plots remain
untouched. Run after both normalized v2 and fixed_benefit_v1 have completed.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys


COORDINATES = ('H0', 'beta_annual', 'chi', 'first_birth_fixed_cost', 'kappa_fert',
               'kappa_fert_continuation', 'theta0', 'delta_alpha_jump',
               'child_benefit_curvature', 'tenure_choice_kappa')
NORMALIZED_STATUS = 'verified_experimental_two_birth_equilibrium_not_adopted'
CONTROL_STATUS = 'audited_fixed_benefit_control_not_demographic_steady_state'
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
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
        os.environ[key] = '1'
    import numpy as np

    here = Path(__file__).resolve().parent
    base = here.parent
    audit = base/'measurement_audit_v1'
    runs = {'normalized_two_birth': here/'run_v1', 'fixed_benefit_control': here/'fixed_benefit_v1'}
    cases = {'reference': base/'resume_v1/selected_export/primary',
             **{name: path/'worker/case' for name, path in runs.items()}}
    out = runs['fixed_benefit_control']/'saved_comparison'
    sources = {Path(__file__).resolve()}

    def load(path):
        sources.add(path)
        return read(path)

    workers = {}
    for name, run in runs.items():
        completion = load(run/'completion.json')
        assert completion['status'] == 'completed' and completion['exit_code'] == 0, name
        workers[name] = load(run/'worker/receipt.json')
        assert Path(workers[name]['case']).resolve() == cases[name].resolve(), name
        assert workers[name]['scientific_promotion'] is False, name
    assert workers['normalized_two_birth']['status'] == 'passed'
    control_worker = workers['fixed_benefit_control']
    assert control_worker['status'] == 'passed_fixed_benefit_diagnostic_checks'
    assert control_worker['calibration_candidate'] is False
    assert control_worker['demographic_renewal_enforced'] is False
    assert control_worker['model_solves'] == 1
    assert control_worker['source_pins'] == workers['normalized_two_birth']['source_pins']
    for filename, digest in control_worker['source_pins'].items():
        path = here/filename
        sources.add(path)
        assert sha(path) == digest, filename
    control_source = here/'fixed_benefit_control.py'
    sources.add(control_source)
    assert sha(control_source) == control_worker['source_sha256']

    receipts, observers, fits, parameters, graph_names = {}, {}, {}, {}, {}
    for name, case in cases.items():
        hashes = load(case/'artifact_hashes.json')
        graph_names[name] = sorted(p for p in hashes if p.startswith('standard_diagnostics/') and p.endswith('.png'))
        assert len(graph_names[name]) == 17, name
        for filename in ['receipt.json','observers.json','target_fit.csv','parameters.csv']+graph_names[name]:
            path = case/filename
            sources.add(path)
            assert sha(path) == hashes[filename], (name, filename)
        receipts[name] = read(case/'receipt.json')
        observers[name] = read(case/'observers.json')['fertility']['uniform_birth_time']
        fit_rows = rows(case/'target_fit.csv'); parameter_rows = rows(case/'parameters.csv')
        fits[name] = {r['moment']: r for r in fit_rows}
        parameters[name] = {r['parameter']: r for r in parameter_rows}
        assert len(fit_rows) == len(fits[name]) == 14, name
        assert len(parameter_rows) == len(parameters[name]) == 31, name
        assert [sum(r['role'] == role for r in fit_rows) for role in ('scored','validation','normalization')] == [10,3,1], name
    assert graph_names['reference'] == graph_names['normalized_two_birth'] == graph_names['fixed_benefit_control']
    assert receipts['normalized_two_birth']['status'] == NORMALIZED_STATUS
    control = receipts['fixed_benefit_control']
    assert control['status'] == CONTROL_STATUS and control['held_reference_child_benefit'] is True
    assert control['control_source_sha256'] == control_worker['source_sha256']
    assert control['source_pins'] == control_worker['source_pins']
    assert control['normalized_count'] == control['normalized_coordinates'] == 0
    assert receipts['normalized_two_birth']['normalized_coordinates'] == 1
    close(control_worker['actual_fertility'], control['normalization']['completed_fertility'], 'control worker fertility')
    close(control_worker['target_fertility'], 2.1, 'control target')
    manifest_values = {}
    for name, run in runs.items():
        receipt = receipts[name]
        assert receipt['scientific_promotion'] is False, name
        assert receipt['estimated_coordinates_this_experiment'] == receipt['free_count'] == 0, name
        assert receipt['held_reference_coordinate_count'] == 10, name
        assert receipt['extra_cache_audit']['status'] == 'passed', name
        assert receipt['reference_checkpoint_sha256'] == receipts['reference']['case_checkpoint_sha256'], name
        manifest_path = run/'worker/effective_source_manifest.json'
        manifest = load(manifest_path)
        assert sha(manifest_path) == receipt['effective_source_manifest_sha256'], name
        manifest_values[name] = {}
        for source_name, entry in manifest.items():
            path = Path(entry['generated_source'])
            assert path.resolve().is_relative_to((run/'worker/effective_sources').resolve())
            sources.add(path)
            assert sha(path) == entry['effective_sha256'], (name, source_name)
            manifest_values[name][source_name] = (entry['original_sha256'],entry['effective_sha256'])
    assert manifest_values['normalized_two_birth'] == manifest_values['fixed_benefit_control'], 'effective household/reporting source differs'

    fit_comparison, parameter_comparison = [], []
    for name in cases:
        assert set(fits[name]) == set(fits['reference']), name
        assert set(parameters[name]) == set(parameters['reference']), name
        close(sum(float(r['loss_contribution'] or 0) for r in fits[name].values()), receipts[name]['loss'], (name,'loss sum'), 1e-9)
    for moment, old in fits['reference'].items():
        entry = dict(moment=moment, role=old['role'], target=old['target'], weight=old['weight'])
        for name in cases:
            row = fits[name][moment]
            for field in ('target','weight','role'):
                assert row[field] == old[field], (name,moment,field)
            close(float(row['model'])-float(row['target']), row['gap'], (name,moment,'gap'))
            if row['loss_contribution']:
                close(float(row['weight'])*float(row['gap'])**2, row['loss_contribution'], (name,moment,'contribution'), 1e-9)
            entry.update({name+'_'+field: row[field] for field in ('model','gap','loss_contribution')})
        fit_comparison.append(entry)
    assert {p for p,r in parameters['reference'].items() if r['lower'] and r['upper']} == set(COORDINATES)
    for parameter, old in parameters['reference'].items():
        entry = dict(parameter=parameter, lower=old['lower'], upper=old['upper'], near_bound=old['near_bound'],
                     reference_calibration_coordinate=parameter in COORDINATES)
        for name in cases:
            row = parameters[name][parameter]
            assert all(row[field] == old[field] for field in ('lower','upper','near_bound')), (name,parameter,'bounds')
            if parameter in COORDINATES:
                assert float(row['estimate']) == float(old['estimate']), (name,parameter,'held coordinate')
                assert float(row['lower']) <= float(row['estimate']) <= float(row['upper']), (name,parameter)
                if name == 'normalized_two_birth':
                    assert row['status'] == 'held at reference estimate during two-birth diagnostic'
                elif name == 'fixed_benefit_control':
                    assert row['status'] == 'held at reference estimate during fixed-benefit control'
            elif name == 'fixed_benefit_control' and parameter == 'psi_child':
                assert row['status'] == 'held at reference child benefit; diagnostic only'
                assert float(row['estimate']) == float(old['estimate'])
            elif name == 'fixed_benefit_control' and parameter == 'child_benefit_CRRA_coefficient':
                assert row['status'] == 'derived from held reference child benefit'
                close(row['estimate'], old['estimate'], 'held-benefit derived coefficient')
            else:
                assert row['status'] == old['status'], (name,parameter,'status')
            entry[name+'_estimate'] = row['estimate']; entry[name+'_status'] = row['status']
        parameter_comparison.append(entry)
    close(control['normalization']['psi_child'], parameters['reference']['psi_child']['estimate'], 'control normalization psi')
    assert control['normalization']['status'] == 'fixed_reference_benefit_diagnostic_not_normalized'
    assert control['adult_entry_gate']['enforced_for_this_diagnostic'] is False

    def project(observer, lo, hi):
        a = observer['accounting']; ages = np.asarray(a['age_cell_start'],dtype=float)
        pre = np.asarray(a['pre_parity_mass_by_age'],dtype=float)
        post = np.asarray(a['post_parity_mass_by_age'],dtype=float)
        assert pre.shape == post.shape == (len(ages),4)
        assert np.array_equal(ages,18+4*np.arange(len(ages)))
        assert np.isfinite(pre).all() and np.isfinite(post).all() and np.all(pre >= 0) and np.all(post >= 0)
        assert np.max(np.abs(pre.sum(axis=1)-post.sum(axis=1))) < 2e-10
        left = np.maximum(ages,lo); right = np.minimum(ages+4,hi+1)
        overlap = np.maximum(right-left,0)/4
        close(overlap.sum()*4,hi+1-lo,'window coverage')
        alpha = ((left+right)/2-ages)/4
        mass = (overlap[:,None]*((1-alpha[:,None])*pre+alpha[:,None]*post)).sum(axis=0)
        assert mass.sum() > 0 and np.all(mass >= 0)
        shares = mass/mass.sum(); capped = float(shares@np.arange(4)); mothers = float(1-shares[0])
        assert mothers > 0
        return dict(capped3=capped,mother_share=mothers,given_mother=capped/mothers,
                    share_0=float(shares[0]),share_1=float(shares[1]),share_2=float(shares[2]),share_3plus=float(shares[3]))

    projection_notes = {}
    for name in runs:
        projection_notes[name] = observers[name]['metadata']['early_fertility_age25']['two_birth_projection_convention']
        assert 'common-event pre/post interpolation' in projection_notes[name]
    assert projection_notes['normalized_two_birth'] == projection_notes['fixed_benefit_control']
    lifecycle_source = audit/'fertility_lifecycle_matched_windows.csv'
    sources.add(lifecycle_source)
    saved_profiles = rows(lifecycle_source)
    assert [(int(r['age_lower']),int(r['age_upper'])) for r in saved_profiles] == [(20,24),(25,29),(30,34),(35,39),(40,44)]
    lifecycle = []
    for saved in saved_profiles:
        lo,hi = int(saved['age_lower']),int(saved['age_upper'])
        entry = dict(age_lower=lo,age_upper=hi,role='untargeted cross-sectional diagnostic')
        for field in ('capped3','mother_share','given_mother'):
            entry['data_'+field] = float(saved['data_'+field])
        for name,observer in observers.items():
            projected = project(observer,lo,hi)
            if name == 'reference':
                for field in ('capped3','mother_share','given_mother'):
                    close(projected[field],saved['model_'+field],('reference lifecycle replay',lo,field))
                close(projected['share_3plus'],saved['model_3plus'],('reference 3plus replay',lo))
            entry.update({name+'_'+field:value for field,value in projected.items()})
        lifecycle.append(entry)
    saved25 = load(audit/'early_fertility_decomposition.json')
    age25 = {'data':dict(capped3=saved25['data_early_fertility'],mother_share=saved25['data_mother_share'],
                        given_mother=saved25['data_capped_children_given_mother'])}
    for name,observer in observers.items():
        age25[name] = project(observer,25,25)
        close(age25[name]['capped3'],fits[name]['early_fertility']['model'],(name,'age25 fit replay'))
        for key in ('0','1','2','3plus'):
            close(age25[name]['share_'+key],observer['ever_born_shares_age25'][key],(name,'age25 shares',key))
        stock = project(observer,40,44)
        close(1-stock['mother_share'],fits[name]['cps_childlessness']['model'],(name,'childlessness replay'))
        close(stock['share_1']/stock['mother_share'],fits[name]['cps_exactly_one']['model'],(name,'one-child replay'))
    for field,key in (('capped3','model_early_fertility'),('mother_share','model_mother_share'),('given_mother','model_capped_children_given_mother')):
        close(age25['reference'][field],saved25[key],('reference age25 audit replay',field))
    decompositions = []
    for left,right in [('reference','data'),('normalized_two_birth','data'),('fixed_benefit_control','data'),
                       ('fixed_benefit_control','reference'),('normalized_two_birth','reference'),
                       ('normalized_two_birth','fixed_benefit_control')]:
        x,y = age25[left],age25[right]
        extensive = (x['mother_share']-y['mother_share'])*(x['given_mother']+y['given_mother'])/2
        intensive = (x['given_mother']-y['given_mother'])*(x['mother_share']+y['mother_share'])/2
        difference = x['capped3']-y['capped3']
        close(extensive+intensive,difference,(left,right,'symmetric identity'))
        decompositions.append(dict(comparison=left+'_minus_'+right,capped_children_difference=difference,
            symmetric_motherhood_component=extensive,symmetric_conditional_children_component=intensive,
            interpretation='Arithmetic difference, not causal attribution; fixed-benefit population need not satisfy demographic renewal'))

    first_source = audit/'first_birth_age_cells.csv'; sources.add(first_source)
    saved_first = rows(first_source); assert len(saved_first) == 7
    first_cells = []
    for i,saved in enumerate(saved_first):
        lo,hi = int(saved['age_lower']),int(saved['age_upper'])
        assert (lo,hi) == (18+4*i,21+4*i)
        entry = dict(age_lower=lo,age_upper=hi,data_share=float(saved['empirical_share']),
                     empirical_tail_mapping=saved['empirical_tail_mapping'],role='untargeted individual first-birth share')
        for name,observer in observers.items():
            a = observer['accounting']; assert a['age_cell_start'][i] == lo and a['first_birth_flow'] > 0
            entry[name+'_share'] = a['parity_birth_flows_by_age'][i][0]/a['first_birth_flow']
        close(entry['reference_share'],saved['model_share'],('reference first-birth replay',i))
        first_cells.append(entry)
    timing = []
    for name in ['data']+list(cases):
        shares = [r[name+'_share'] for r in first_cells]
        assert all(0 <= s <= 1 for s in shares)
        close(sum(shares),1,(name,'first-birth shares sum'))
        mean = sum((r['age_lower']+2)*r[name+'_share'] for r in first_cells)
        late = sum(r[name+'_share'] for r in first_cells if r['age_lower'] >= 30)
        fit,field = (fits['reference'],'target') if name == 'data' else (fits[name],'model')
        close(mean,fit['nchs_mean_age'][field],(name,'mean age replay'),1e-10)
        close(late,fit['nchs_share30'][field],(name,'share30 replay'),1e-10)
        timing.append(dict(series=name,mapped_mean_first_birth_age=mean,share_first_births_age30plus=late))

    demographic = []
    for name,receipt in receipts.items():
        normalization = receipt['normalization']; renewal = receipt['adult_entry_gate']
        fertility = float(normalization['completed_fertility'])
        target = float(fits[name]['initial_normalization']['target'])
        E,B = float(renewal['entry_E']),float(renewal['potential_B'])
        assert E > 0
        close(target,2.1,(name,'replacement target'))
        close(fertility,fits[name]['initial_normalization']['model'],(name,'normalization fit replay'),1e-10)
        close(fertility-target,fits[name]['initial_normalization']['gap'],(name,'fertility gap replay'),1e-10)
        close(E-B,renewal['entry_residual'],(name,'entry residual replay'),1e-10)
        close(float(renewal['births_per_entry'])-2.1,renewal['fertility_gap'],(name,'renewal fertility gap'),1e-10)
        close(float(renewal['births_per_entry']),fertility,(name,'fertility versus renewal'),5e-9)
        demographic.append(dict(series=name,psi_child=float(parameters[name]['psi_child']['estimate']),
            completed_fertility=fertility,target=target,completed_fertility_gap=fertility-target,
            entry_E=E,potential_birth_generated_entry_B=B,renewal_residual_E_minus_B=E-B,
            renewal_residual_relative_to_E=(E-B)/E,births_per_entry=float(renewal['births_per_entry']),
            renewal_enforced=name!='fixed_benefit_control',
            case_status=receipt.get('status','reference'),
            interpretation='fixed-benefit mechanism control; not demographic steady state or calibration candidate'
                           if name=='fixed_benefit_control' else 'replacement-normalized stationary reference or experiment; no automatic adoption'))

    out.mkdir(parents=True,exist_ok=True)
    outputs = []
    def table(filename,records):
        path = out/filename
        with path.open('w',newline='') as stream:
            writer = csv.DictWriter(stream,fieldnames=list(records[0]),lineterminator='\n')
            writer.writeheader();writer.writerows(records)
        outputs.append(path)
    table('target_fit_comparison.csv',fit_comparison)
    table('parameters_comparison.csv',parameter_comparison)
    table('fertility_lifecycle_comparison.csv',lifecycle)
    table('age25_comparison.csv',[dict(series=name,**{key:value[key] for key in ('capped3','mother_share','given_mother')}) for name,value in age25.items()])
    table('age25_symmetric_decomposition.csv',decompositions)
    table('first_birth_age_cells_comparison.csv',first_cells)
    table('first_birth_timing_comparison.csv',timing)
    table('fertility_and_renewal_comparison.csv',demographic)
    summary = dict(status='verified_saved_comparison_including_nonstationary_demographic_control',
        model_solves=0,model_imports=0,checkpoint_reads=0,slurm_job=os.environ['SLURM_JOB_ID'],
        reference_label=LABEL,control_status=control['status'],control_is_calibration_candidate=False,
        control_demographic_renewal_enforced=False,scientific_promotion=False,transition_run=False,
        all_ten_reference_coordinates_unchanged=True,fixed_benefit_equals_reference=True,
        complete_target_rows=14,complete_parameter_rows=31,all_targets_weights_roles_unchanged=True,
        reference_lifecycle_age25_first_birth_replay=True,source_pins_identical_to_normalized_v2=True,
        household_and_reporting_overlay_hashes_identical=True,projection_notes=projection_notes,
        checkpoint_identities={name:r['case_checkpoint_sha256'] for name,r in receipts.items()},
        effective_source_manifest_identities={name:receipts[name]['effective_source_manifest_sha256'] for name in runs},
        age25=age25,demographic=demographic,first_birth_timing=timing,age25_symmetric_decompositions=decompositions,
        scored_discrepancies={name:r['loss'] for name,r in receipts.items()},
        scored_discrepancy_caution='The fixed-benefit score does not enforce replacement fertility or renewal and cannot establish an admissible calibrated improvement.',
        economic_changes={name:receipts[name]['economic_changes'] for name in runs},
        standard_17_plot_names=graph_names['reference'],standard_17_hashes_authenticated=True,standard_17_plots_changed=False,
        caveats=['Model household reproductive-member exposure approximates CPS female exposure.',
                 'Identical five-year interview-age windows; retained model age weights and linear pre/post interpolation.',
                 'Two within-cell births share the inherited event-time proxy; ordered within-cell dates are unavailable.',
                 'NCHS first-birth shares are pooled period 2003–2006 counts with boundary tails; CPS stocks pool 2004/2006 women.',
                 'The fixed-benefit control uses the maintained age population and reports its renewal residual without enforcing demographic closure; it is not a demographic steady state.'],
        input_hashes={str(path):sha(path) for path in sorted(sources)},output_hashes={path.name:sha(path) for path in outputs})
    (out/'analysis_receipt.json').write_text(json.dumps(summary,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps(dict(status=summary['status'],output=str(out),model_solves=0,target_rows=14,parameter_rows=31)))


if __name__ == '__main__':
    main()
