#!/usr/bin/env python3
"""Render the stable 17 diagnostics with verified reporting-only corrections.

This driver reads solved states. It never solves household or market problems,
changes numerical kernels, or writes into the numerical transition output tree.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import math
from pathlib import Path
import time
import traceback
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import run_e5f_simple_fertility_tax_transition as transition
from intergen_eqscale_seq_optimized import diagnostics, joint_nested

common, adapter, audit = transition.common, transition.adapter, transition.audit
BASE_GRAPHS = {
    'ownership_by_age.png', 'fertility_by_age.png', 'housing_market.png',
    'market_clearing_by_market.png', 'market_clearing_residuals.png',
    'tenure_services.png', 'owner_rungs.png', 'housing_prices.png',
    'income_state_outcomes.png', 'ownership_by_age_income_state.png',
    'liquid_wealth_by_age_income_state.png', 'housing_by_age_income_state.png',
    'fertility_policy_by_age_income_state.png',
}
REQUIRED_ARTIFACTS = (
    'dated_state.pkl.gz', 'endpoint_receipt.json', 'advance.json', 'reproduction.json',
    'quantities.json', 'fiscal_ledger.json', 'budget_summary.json',
    'policy_array_summary.json', 'age_measurement.json', 'family_group_quantities.csv', 'lifecycle.csv',
)


def finite_gate(value, maximum, label):
    number = float(value)
    if not math.isfinite(number) or abs(number) > maximum:
        raise RuntimeError(f'{label} gate failed: {value}')


def expected_graphs(P):
    result = set(BASE_GRAPHS)
    for age in np.asarray(getattr(P, 'diagnostic_policy_ages', [30., 42.]), dtype=float).reshape(-1):
        j = diagnostics.age_to_index(P, float(age))
        node = int(round(P.age_start + j*P.da))
        result.add(f'wealth_dist_childless_renter_age{node}.png')
        result.add(f'policy_childless_renter_age{node}.png')
    if len(result) != 17:
        raise RuntimeError('Saved model does not imply the established seventeen-graph set')
    return result


def array_hash(values):
    value = np.ascontiguousarray(values)
    header = f'{value.dtype.str}:{value.shape}:'.encode()
    return hashlib.sha256(header + value.tobytes()).hexdigest()


def owner_probability(joint, P, origin):
    """Ex-ante ownership from wait, failed-attempt and successful-attempt plans."""
    j = int(origin[3])
    pi = float(common.solver.get_fecundity_by_age(P)[j])
    pr = joint.probabilities[origin + (slice(None), slice(None))]
    failure = joint.failure_probabilities[origin + (slice(None),)]
    return (pr[..., 1:, 0] + (1-pi)*failure[..., 1:] + pi*pr[..., 1:, 1]).sum(axis=-1)


class OwnerPlotProbabilities:
    """Narrow rendering accessor for the legacy plot's sole probability lookup.

    Only the obsolete product-one lookup is adapted. The numerical policy and
    original JointChoice object are never replaced in the saved model state.
    """
    def __init__(self, joint, P):
        self.joint, self.P = joint, P

    def __getitem__(self, key):
        if len(key) != 9 or key[-2] != 1 or key[-1] != slice(None):
            raise RuntimeError('Unexpected legacy owner-plot probability access')
        owner = owner_probability(self.joint, self.P, key[:7])
        return np.stack((owner, np.zeros_like(owner)), axis=-1)


def reporting_data(packet, folder):
    P, e = packet['parameters'], packet['evaluation']
    rows, hazards, birth_total, post_gap = [], np.zeros(P.J), 0., 0.
    for j in range(P.J):
        post, effective, births, attempts, risk = joint_nested.factor_age(
            e.g_pre[:, :, :, j], e.policy.joint_choice, P, j)
        birth_total += float(births.sum())
        post_gap += float(np.abs(post-e.g_post_fertility[:, :, :, j]).sum())
        hazards[j] = float(births[0]/risk[0]) if risk[0] > 0 else 0.
        row = dict(age_node=float(P.age_start+j*P.da), first_births=float(births[0]),
                   first_birth_risk_mass=float(risk[0]), first_birth_attempts=float(attempts[0]),
                   first_birth_probability=hazards[j], all_explicit_births=float(births.sum()))
        for n in range(P.n_parity):
            row[f'raw_births_origin_parity_{n}'] = float(births[n])
        rows.append(row)
        del post, effective
    if abs(birth_total-float(e.births)) > 2e-10 or post_gap > 2e-10:
        raise RuntimeError('Reporting birth kernel does not reproduce saved births/population')
    common.policy.baseline.write_csv(folder/'births_by_age_reporting.csv', rows)
    z, _, _ = common.solver.income_transition_values(P)
    income_rows = []
    for j in range(P.J):
        for zz, value in enumerate(z):
            g = e.g_current[:, :, :, j, zz]
            mass, rent_mass, own_mass = float(g.sum()), float(g[:, 0].sum()), float(g[:, 1:].sum())
            rent_h = float(np.sum(g[:, 0]*e.policy.hR_pol[:, 0, :, j, zz]))
            own_h = sum(float(g[:, t].sum())*float(P.H_own[t-1]) for t in range(1,g.shape[1]))
            income_rows.append(dict(age_node=float(P.age_start+j*P.da), z_index=zz,
                z=float(value), permanent_group=int(P.permanent_income_group_index[zz]),
                persistent_state=int(P.permanent_income_base_state_index[zz]),
                working_income=float(P.income[0,j])*float(value),
                actual_income_including_rebate=common.solver.income_at_state(P,0,j,float(value)),
                is_working_age=bool(j < P.J_R), household_mass=mass,
                owner_share=own_mass/mass if mass else None,
                renter_mass=rent_mass, owner_mass=own_mass,
                renter_mean_rooms=rent_h/rent_mass if rent_mass else None,
                owner_mean_rooms=own_h/own_mass if own_mass else None,
                mean_rooms=(rent_h+own_h)/mass if mass else None))
    common.policy.baseline.write_csv(folder/'age_income_housing_audit.csv', income_rows)
    total_h = sum(r['household_mass']*(r['mean_rooms'] or 0.) for r in income_rows)
    if abs(total_h-float(e.demand_by_loc.sum())) > 2e-10:
        raise RuntimeError('Age-income room decomposition does not reproduce aggregate demand')
    return hazards, dict(explicit_birth_gap=abs(birth_total-float(e.births)),
                         post_fertility_l1=post_gap, housing_aggregation_gap=abs(total_h-float(e.demand_by_loc.sum())))


def render_corrected_diagnostics(packet, folder):
    """Adapt rendering only; retain original axes structure and graph filenames."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from cycler import cycler
    P, e = packet['parameters'], packet['evaluation']
    hazards, validation = reporting_data(packet, folder)
    z = np.asarray(P.z_grid)
    colors = list(plt.get_cmap('tab20')(np.linspace(0,1,len(z))))
    if len({tuple(c) for c in colors}) != len(z):
        raise RuntimeError('Income palette does not provide one distinct color per state')
    state_index = {f'z={value:g}': index for index,value in enumerate(z)}
    if len(state_index) != len(z):
        raise RuntimeError('Existing income legend labels do not uniquely identify states')
    original_write = audit.write_diagnostics
    original_policy = diagnostics.plot_policy_childless_renter
    original_save = Figure.savefig
    owner_rows, figure_checks = [], []
    maximum_operator_gap = 0.
    hazard_change = {}

    def corrected_write(sol, params, out):
        view = copy.copy(sol)
        hazard_change.update(old_array_sha256=array_hash(sol.fert_by_age), new_array_sha256=array_hash(hazards))
        view.fert_by_age = hazards.copy()
        return original_write(view, params, out)

    def corrected_policy(sol, params, out, pyplot, b_grid, z_grid, j):
        nonlocal maximum_operator_gap
        joint = sol.joint_choice
        valid = np.asarray(sol.V)[:,0,0,j,:,0,0] > -1e9
        probe = np.zeros_like(e.g_pre[:,:,:,j])
        probe[:,0,0,:,0,0] = valid.astype(float)
        post, effective, _, _, _ = joint_nested.factor_age(probe,joint,params,j)
        native_owner = (post*effective[...,1:].sum(axis=-1)).sum(axis=(1,2,4,5))
        for zz, value in enumerate(z_grid):
            origin = (slice(None),0,0,j,zz,0,0)
            new = owner_probability(joint,params,origin)
            old = joint.probabilities[origin+(1,slice(None))].sum(axis=-1)
            gap = float(np.max(np.abs(np.where(valid[:,zz],new,0.)-native_owner[:,zz])))
            if gap > 2e-11:
                raise RuntimeError('Owner reporting formula differs from native product factorization')
            maximum_operator_gap = max(maximum_operator_gap,gap)
            for b, wealth in enumerate(b_grid):
                owner_rows.append(dict(age_node=float(params.age_start+j*params.da),wealth=float(wealth),
                    z_index=zz,z=float(value),valid_origin=bool(valid[b,zz]),
                    old_product_one_probability=float(old[b]),
                    corrected_exante_owner_probability=float(new[b]),
                    native_factorization_owner_probability=float(native_owner[b,zz]) if valid[b,zz] else None))
        view = copy.copy(sol)
        view.joint_choice = SimpleNamespace(probabilities=OwnerPlotProbabilities(joint,params))
        return original_policy(view,params,out,pyplot,b_grid,z_grid,j)

    def line_hashes(fig):
        return [[(array_hash(line.get_xdata()),array_hash(line.get_ydata())) for line in ax.lines] for ax in fig.axes]

    def styled_save(fig, filename, *args, **kwargs):
        name = Path(filename).name
        before = line_hashes(fig)
        if name == 'fertility_by_age.png':
            fig.axes[0].set_ylabel('First-birth probability among at-risk households\n(4-year period)')
            fig.axes[0].set_title('First-birth probability by age')
        if name == 'income_state_outcomes.png':
            fig.axes[1].set_ylabel('Completed parity (3+ top-coded)')
            fig.axes[0].tick_params(axis='x',labelrotation=55,labelsize=7)
            legend = fig.axes[0].get_legend()
            if legend is not None:
                for text in legend.get_texts():
                    if text.get_text() == 'children':
                        text.set_text('completed parity')
        income_plot = (name.endswith('_by_age_income_state.png')
                       or name.startswith('wealth_dist_childless_renter_')
                       or name.startswith('policy_childless_renter_'))
        if income_plot:
            for ax in fig.axes:
                for line in ax.lines:
                    if line.get_label() in state_index:
                        line.set_color(colors[state_index[line.get_label()]])
            handles, labels = fig.axes[0].get_legend_handles_labels()
            for ax in fig.axes:
                legend = ax.get_legend()
                if legend is not None:
                    legend.remove()
            labels = [f's{state_index[label]+1}: {label}' if label in state_index else label for label in labels]
            fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,0.),ncol=3,fontsize=7)
            fig.tight_layout(rect=(0,.24 if len(fig.axes)==1 else .17,1,.97))
        else:
            fig.tight_layout()
        after = line_hashes(fig)
        if before != after:
            raise RuntimeError('Legend/label formatting changed numerical line arrays')
        figure_checks.append(dict(filename=name,line_array_sha256=after,
                                  formatting_preserved_line_arrays=True))
        return original_save(fig,filename,*args,**kwargs)

    with mpl.rc_context({'axes.prop_cycle':cycler(color=colors)}), \
            patch.object(audit,'write_diagnostics',corrected_write), \
            patch.object(diagnostics,'plot_policy_childless_renter',corrected_policy), \
            patch.object(Figure,'savefig',styled_save):
        audit.standard_diagnostics(packet,folder,validate_production_young=False)
    common.policy.baseline.write_csv(folder/'owner_probability_reporting_correction.csv',owner_rows)
    receipt = dict(status='complete',reporting_only=True,model_kernels_unchanged=True,
        corrections=['all-owner-product conception-weighted probability','pre-choice first-birth hazard',
                     'topcoded parity label','distinct income colors and external legends'],
        maximum_native_owner_operator_gap=maximum_operator_gap,first_birth_hazard_change=hazard_change,
        validation=validation,figure_checks=figure_checks,
        unchanged_other_arrays_scope='Original plotting functions/data remain in use; only owner and first-birth y-arrays are adapted. Formatting line-array equality checked for all figures.',
        income_labels='z is permanent-level times persistent-state working-age multiplier; state order retained',
        income_palette_rgba=[list(c) for c in colors])
    audit.save_json(folder/'plot_correction_receipt.json',receipt)
    return receipt


def verify_date(c, contract_hash, case, year, stage):
    """Verify the completed stage receipt, then only artifacts used for this date."""
    folder = Path(c['output_root'])/stage/case
    path = folder/'receipt.json'
    receipt_hash = (folder/'receipt.sha256').read_text().strip()
    adapter.verify(path, receipt_hash)
    receipt = adapter.read_json(path)
    years = transition.YEARS[:2] if stage == 'smoke' else transition.YEARS
    if (receipt['status'] != 'complete' or receipt['stage'] != stage or receipt['case'] != case
            or receipt['years'] != years or year not in years
            or receipt['contract_sha256'] != contract_hash):
        raise RuntimeError('Completed matching transition stage/date required')
    if (receipt['selected_summary_sha256'] != c['selected_summary_sha256']
            or receipt['selected_checkpoint_sha256'] != c['checkpoint_sha256']
            or receipt['scientific_bundle_sha256'] != c['code_bundle_sha256']
            or receipt['closure'] != c['closure']):
        raise RuntimeError('Transition stage provenance differs from its frozen contract')
    adapter.verify(c['selected_summary'], c['selected_summary_sha256'])
    dated = folder/f'date_{year}'
    verified = {}
    for name in REQUIRED_ARTIFACTS:
        relative = f'date_{year}/{name}'
        if relative not in receipt['artifact_sha256']:
            raise RuntimeError(f'Date artifact absent from completed receipt: {relative}')
        expected = receipt['artifact_sha256'][relative]
        adapter.verify(dated/name, expected)
        verified[name] = expected
    endpoint = adapter.read_json(dated/'endpoint_receipt.json')
    advance = adapter.read_json(dated/'advance.json')
    reproduction = adapter.read_json(dated/'reproduction.json')
    if (endpoint['status'] != 'complete' or endpoint['policy'] != case
            or endpoint['calendar_year'] != year or advance['calendar_year'] != year
            or advance['next_calendar_year'] != year+4
            or advance['row']['calendar_year'] != year):
        raise RuntimeError('Dated validity or calendar metadata missing')
    if endpoint['checkpoint_sha256'] != verified['dated_state.pkl.gz']:
        raise RuntimeError('Endpoint checkpoint digest differs from completed-stage digest')
    if reproduction != receipt['reproduction'][str(year)]:
        raise RuntimeError('Date reproduction receipt differs from completed stage')
    if year == 2023 and 'verified_impact' not in reproduction:
        raise RuntimeError('Initial date lacks verified impact reproduction')
    if stage == 'full' and year in transition.YEARS[:2] and 'smoke' not in reproduction:
        raise RuntimeError('Full path lacks required smoke reproduction')
    finite_gate(advance['row']['mass_accounting_residual'], 2e-10, 'cohort mass')
    finite_gate(endpoint['budget']['budget_excess_mass'], 2e-10, 'household budget')
    if endpoint['policy_arrays']['occupied_negative_steps'] != 0:
        raise RuntimeError('Occupied value audit failed')
    for value in endpoint['normalized_root_residuals']:
        finite_gate(value, 1e-4, 'coupled root')
    finite_gate(endpoint['ledger']['government_budget_residual'], 2.5e-5, 'fiscal balance')
    finite_gate(endpoint['quantities']['market_residual'], 2e-4, 'housing market')
    finite_gate(endpoint['quantities']['feasibility_projection'], 1e-6, 'feasibility projection')
    for value in endpoint['gates']['mass_gaps'].values():
        finite_gate(value, 2e-10, 'dated mass')
    for bounds in endpoint['policy_arrays']['probabilities'].values():
        if bounds['nonfinite'] or bounds['minimum'] < 0 or bounds['maximum'] > 1:
            raise RuntimeError('Saved marginal probability audit failed')
    return dated, endpoint, dict(stage_receipt_sha256=receipt_hash, verified_date_artifact_sha256=verified)


def execute(args):
    c = transition.load_contract(args.transition_contract, args.transition_contract_sha256, launching=False)
    output_root = args.output_root.resolve()
    if not args.output_root.is_absolute():
        raise RuntimeError('Graph output root must be absolute')
    numerical_root = Path(c['output_root']).resolve()
    if output_root == numerical_root or output_root.is_relative_to(numerical_root):
        raise RuntimeError('Graphs must be written outside the numerical output tree')
    dated, endpoint, provenance = verify_date(c,args.transition_contract_sha256,args.case,args.year,args.stage)
    folder = output_root/args.case/f'date_{args.year}'
    folder.mkdir(parents=True,exist_ok=False)
    started = time.monotonic()
    source_paths = (Path(__file__).resolve(),Path(audit.__file__).resolve(),Path(diagnostics.__file__).resolve())
    graphics_hashes = {str(path):adapter.digest(path) for path in source_paths}
    try:
        packet = transition.load_packet(dated/'dated_state.pkl.gz')
        if packet['calendar_year'] != args.year or packet['policy_case'] != args.case:
            raise RuntimeError('Saved checkpoint date or case differs')
        common.assert_model(packet['parameters'],packet['supply_rule'])
        names = expected_graphs(packet['parameters'])
        corrections = render_corrected_diagnostics(packet,folder)
        adapter.verify(dated/'dated_state.pkl.gz',endpoint['checkpoint_sha256'])
        corrections['input_checkpoint_sha256_unchanged'] = endpoint['checkpoint_sha256']
        audit.save_json(folder/'plot_correction_receipt.json',corrections)
        graph_folder = folder/'standard_diagnostics'
        graphs = sorted(graph_folder.glob('*.png'))
        if {p.name for p in graphs} != names:
            raise RuntimeError('Stable seventeen-graph filename set is incomplete or changed')
        units = adapter.read_json(folder/'market_quantity_units.json')
        transition.verify_sources(c)
        for path,digest in graphics_hashes.items():
            adapter.verify(path,digest)
        files = sorted(p for p in folder.rglob('*') if p.is_file())
        manifest = dict(status='complete',case=args.case,calendar_year=args.year,stage=args.stage,
                        transition_contract_sha256=args.transition_contract_sha256,
                        selected_summary_sha256=c['selected_summary_sha256'],
                        selected_checkpoint_sha256=c['checkpoint_sha256'],
                        scientific_bundle_sha256=c['code_bundle_sha256'],
                        input_checkpoint=str(dated/'dated_state.pkl.gz'),
                        input_checkpoint_sha256=endpoint['checkpoint_sha256'],
                        input_provenance=provenance,graphics_source_sha256=graphics_hashes,
                        standard_png_count=len(graphs),standard_filenames=sorted(names),
                        graph_sha256={p.name:adapter.digest(p) for p in graphs},
                        artifact_sha256={str(p.relative_to(folder)):adapter.digest(p) for p in files},
                        units=units,age_measurement=endpoint['age_measurement'],
                        reporting_corrections=corrections,
                        elapsed_seconds=time.monotonic()-started,
                        scope='Read-only rendering of verified temporary-equilibrium household-unit diagnostic; not a resident-population forecast',
                        numerical_solves=0,kernel_changes=False,production_promoted=False)
        audit.save_json(folder/'graph_manifest.json',manifest)
    except BaseException as error:
        audit.save_json(folder/'failure.json',dict(error=str(error),type=type(error).__name__,traceback=traceback.format_exc()))
        raise


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--transition-contract',type=Path,required=True)
    parser.add_argument('--transition-contract-sha256',required=True)
    parser.add_argument('--case',choices=transition.CASES,required=True)
    parser.add_argument('--year',type=int,choices=transition.YEARS,required=True)
    parser.add_argument('--stage',choices=('smoke','full'),default='full')
    parser.add_argument('--output-root',type=Path,required=True)
    execute(parser.parse_args())
