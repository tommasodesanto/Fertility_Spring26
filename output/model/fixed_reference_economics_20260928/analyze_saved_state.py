"""Supplemental occupied-state anatomy; one authenticated saved checkpoint, zero solves.

Torch/Slurm only, inside the original pinned container mount. The lead supplies
--output (a new directory) and the launch/time/memory budget. This script calls
no Bellman, equilibrium, KFE, normalization, or transition solver. Standard
plots and reference artifacts are read-only. Source accounting references:
run_e5f_open_population_transition.apply_sequential_fertility;
solver.realize_current_choices / realize_stayer_cross_section;
solver.add_aggregate_wealth_bequest_flow_moments.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import sys
import time

ROOT = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
BASE = ROOT / 'output/model/fertility_identification_20260928'
SOURCE = BASE / 'resume_v1/selected_export/primary'
LABEL = '2007 stationary reference — block0506, September 28 verified export'
CONTRACT_SHA = '68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf'
CHECKPOINT_SHA = 'b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def main(out):
    if sys.platform != 'linux' or not os.environ.get('SLURM_JOB_ID', '').isdigit():
        raise RuntimeError('Run only in a bounded Torch Slurm allocation, never on the Mac')
    out = out.resolve()
    if out == SOURCE.resolve() or SOURCE.resolve() in out.parents:
        raise ValueError('Reference export is read-only')
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    sys.path.insert(0, str(ROOT / 'code/model/tools'))
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import run_e5f_fertility_identification as driver
    import e5f_evening_calibration_runtime as runtime

    def write(name, value):
        (out / name).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')

    def table(name, rows):
        with (out / name).open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    def ratio(num, den):
        return float(num / den) if den > 0 else None

    def check_close(actual, expected, label, tol=2e-10):
        gap = float(np.max(np.abs(np.asarray(actual) - np.asarray(expected))))
        if not np.isfinite(gap) or gap > tol:
            raise ValueError(f'{label}: {gap:g} exceeds {tol:g}')
        checks[label] = gap

    def quantiles(weights, qs=(.005, .05, .5, .95, .995)):
        if weights.sum() <= 0:
            return [None] * len(qs)
        cdf = np.cumsum(weights) / weights.sum()
        return [float(bg[min(np.searchsorted(cdf, q), len(bg)-1)]) for q in qs]

    def savefig(fig, name):
        fig.suptitle('Supplemental: ' + LABEL, fontsize=9)
        fig.tight_layout(rect=(0, 0, 1, .965))
        fig.savefig(out / name, dpi=150)
        plt.close(fig)

    print('Authenticating frozen sources and checkpoint', flush=True)
    contract_path = BASE / 'contract_v1/contract.json'
    assert sha(contract_path) == CONTRACT_SHA
    c, objectives = driver.verify(contract_path)
    manifest_path = BASE / 'fixed_reference_manifest.json'
    manifest = read(manifest_path)
    assert manifest['label'] == LABEL and manifest['checkpoint']['sha256'] == CHECKPOINT_SHA
    assert manifest['contract']['sha256'] == CONTRACT_SHA
    hashes = read(SOURCE / 'artifact_hashes.json')
    for name, digest in hashes.items():
        assert sha(SOURCE / name) == digest, name
    assert sum(k.endswith('.png') for k in hashes) == 17
    receipt = read(SOURCE / 'receipt.json')
    checkpoint = SOURCE / 'initial_state.pkl.gz'
    assert sha(checkpoint) == CHECKPOINT_SHA == receipt['case_checkpoint_sha256']
    evaluator = runtime.setup(dict(c, objective=c['lanes']['primary']['objective']),
                              objectives['primary'], out / 'runtime_preparation')
    with gzip.open(checkpoint, 'rb') as stream:
        packet = pickle.load(stream)
    P, e, bg = packet['parameters'], packet['evaluation'], np.asarray(packet['b_grid'])
    p, model = e.policy, evaluator.rt['model']
    gp, post, current = map(np.asarray, (e.g_pre, e.g_post_fertility, e.g_current))
    shape = (P.Nb, 1+P.n_house, P.I, P.J, P.Nz, P.n_parity, P.n_child_states)
    assert gp.shape == post.shape == current.shape == shape
    assert P.I == 1 and P.sequential_births and not P.joint_nested_choice
    assert P.child_state_mode == 'independent_count' and P.n_parity == 4
    assert model.readiness_settled_state(P) == 0
    assert float(P.psi_child) == manifest['economic_contract']['frozen_psi_child']
    checks = {}
    for name, g in (('pre', gp), ('post_birth', post), ('current', current)):
        assert np.isfinite(g).all() and g.min() >= 0
        check_close(g.sum(), 1., name + '_mass')
    check_close(gp.sum(axis=(0,1,2,4,5,6)), current.sum(axis=(0,1,2,4,5,6)), 'age_mass')
    observer = read(SOURCE / 'observers.json')
    fertility_saved = observer['fertility']['uniform_birth_time']['accounting']
    for g, key in ((gp, 'pre_parity_mass_by_age'), (current, 'post_parity_mass_by_age')):
        check_close(g.sum(axis=(0,1,2,4,6)), fertility_saved[key], key)
    fit = list(csv.DictReader((SOURCE / 'target_fit.csv').open()))
    loss = sum(float(r['weight'])*float(r['gap'])**2 for r in fit if r['role']=='scored')
    check_close(loss, manifest['common_primary_loss'], 'primary_loss')
    # Small tables accompany every saved-solution readout; no new target definitions.
    for name in ('target_fit.csv', 'parameters.csv'):
        (out / name).write_bytes((SOURCE / name).read_bytes())

    ages = P.age_start + P.da*np.arange(P.J)
    fec = np.asarray(model.get_fecundity_by_age(P))
    birth_by_age = np.zeros((P.J, 3))
    birth_rows, birth_cells = [], []
    print('Extracting births and occupied distributions', flush=True)
    # Origin wealth/income/tenure are BEFORE this period's birth and transaction.
    # No same-period chain births: each flow uses its original gp pool.
    for j, age in enumerate(ages):
        total_b = gp[:,:,:,j].sum(axis=(1,2,3,4,5))
        midpoint_cdf = (np.cumsum(total_b)-.5*total_b)/total_b.sum()
        quintile = np.minimum((5*midpoint_cdf).astype(int), 4)
        fertile = P.A_f_start <= j+1 <= P.A_f_end
        for n in range(3):
            risk = np.zeros((P.Nb, 1+P.n_house, P.I, P.Nz))
            attempts = np.zeros_like(risk)
            if fertile:
                for m in range(n+1):
                    pool = gp[:,:,:,j,:,n,m]
                    pr = (p.fert_probs[:,:,:,j,:,1] if n == 0 else
                          p.fert2_probs[:,:,:,j,:,1,n-1,m])
                    risk += pool
                    attempts += pool*pr
            born = attempts * fec[j] if fertile else attempts
            birth_by_age[j,n] = born.sum()
            birth_rows.append(dict(age_left=float(age), age_right=float(age+P.da), birth_order=n+1,
                at_risk_mass=float(risk.sum()), attempt_mass=float(attempts.sum()),
                birth_mass=float(born.sum()), probability_per_at_risk_cell=ratio(born.sum(),risk.sum()),
                births_per_living_age_household=ratio(born.sum(),total_b.sum()),
                fertile_cell=bool(fertile)))
            if not fertile:
                continue
            for zz, z in enumerate(P.z_grid):
                for tenure, ts in (('renter', slice(0,1)), ('owner', slice(1,None))):
                    for q in range(5):
                        select = quintile == q
                        rr = risk[select,ts,:,zz].sum()
                        bb = born[select,ts,:,zz].sum()
                        aa = attempts[select,ts,:,zz].sum()
                        birth_cells.append(dict(age_left=float(age), birth_order=n+1,
                            origin_wealth_quintile=q+1, income_index=zz, income_multiplier=float(z),
                            origin_tenure=tenure, at_risk_mass=float(rr), attempt_mass=float(aa),
                            birth_mass=float(bb), birth_probability=ratio(bb,rr)))
    for n, word in enumerate(('first','second','third')):
        check_close(birth_by_age[:,n], getattr(P, '_'+word+'_births_by_age'), word+'_births_saved')
    check_close(birth_by_age.sum(), e.births, 'explicit_births')
    table('supplemental_births_by_age.csv', birth_rows)
    table('supplemental_births_by_origin.csv', birth_cells)

    from e5f_overnight_estate_audit import policy_mass_branches
    branches = policy_mass_branches(e, P)
    stay = np.asarray(e.g_stay_distribution)
    assert stay.shape == shape and np.isfinite(stay).all() and stay.min() >= 0
    assert np.max(stay-current) < 2e-10 and np.max(stay[:,0]) == 0
    check_close(stay.sum(), manifest['inherited_gates']['purchase_accounting']['stayer_mass'], 'stayer_mass')
    price = float(p.price[0])
    houses = np.r_[0., P.H_own]
    hvalue = price*houses
    lifecycle, wealth_rows, income_rows = [], [], []
    for j, age in enumerate(ages):
        age_mass = current[:,:,:,j].sum()
        cur_bt = current[:,:,:,j].sum(axis=(2,3,4,5))
        begin_bt = post[:,:,:,j].sum(axis=(2,3,4,5))
        stay_b = stay[:,:,:,j].sum(axis=(1,2,3,4,5))
        current_b = cur_bt.sum(axis=1)
        room_b = (cur_bt*houses).sum(axis=1)
        room_b += np.sum(np.where(current[:,0,:,j]>0,
            current[:,0,:,j]*p.hR_pol[:,0,:,j], 0.), axis=(1,2,3,4))
        qb = quantiles(current_b)
        row = dict(age_left=float(age), age_right=float(age+P.da), household_mass=float(age_mass),
            retired=bool(j>=P.J_R), ownership=float(cur_bt[:,1:].sum()/age_mass),
            owner_stayer_mass=float(stay_b.sum()), buyer_mass=float(cur_bt[:,1:].sum()-stay_b.sum()),
            realized_rooms=float(room_b.sum()/age_mass),
            beginning_net_financial_wealth=float((begin_bt*bg[:,None]).sum()/age_mass),
            beginning_gross_housing_value=float((begin_bt*hvalue).sum()/age_mass),
            beginning_net_worth=float((begin_bt*(bg[:,None]+hvalue)).sum()/age_mass),
            post_transaction_net_financial_wealth=float(current_b@bg/age_mass),
            post_transaction_gross_housing_value=float((cur_bt*hvalue).sum()/age_mass),
            post_transaction_net_worth=float((cur_bt*(bg[:,None]+hvalue)).sum()/age_mass),
            realized_consumption=float(sum(np.sum(mass[:,:,:,j]*cons[:,:,:,j]) for mass,saving,cons in branches)/age_mass),
            realized_end_financial_wealth=float(sum(np.sum(mass[:,:,:,j]*saving[:,:,:,j]) for mass,saving,cons in branches)/age_mass),
            net_financial_q005=qb[0], net_financial_q05=qb[1], net_financial_q50=qb[2],
            net_financial_q95=qb[3], net_financial_q995=qb[4])
        lifecycle.append(row)
        for b, wealth in enumerate(bg):
            mm = current_b[b]
            wealth_rows.append(dict(age_left=float(age), net_financial_wealth=float(wealth),
                household_mass=float(mm), share_within_age=float(mm/age_mass),
                owner_mass=float(cur_bt[b,1:].sum()), owner_stayer_mass=float(stay_b[b]),
                buyer_mass=float(cur_bt[b,1:].sum()-stay_b[b]), ownership=ratio(cur_bt[b,1:].sum(),mm),
                realized_rooms=ratio(room_b[b],mm), inside_q005_q995=bool(qb[0]<=wealth<=qb[4])))
        for zz, z in enumerate(P.z_grid):
            mass = current[:,:,:,j,zz].sum()
            owner = current[:,1:,:,j,zz].sum()
            rooms = sum(current[:,t,:,j,zz].sum()*houses[t] for t in range(1,len(houses)))
            rooms += np.sum(np.where(current[:,0,:,j,zz]>0,
                current[:,0,:,j,zz]*p.hR_pol[:,0,:,j,zz], 0.))
            income_rows.append(dict(age_left=float(age), income_index=zz, income_multiplier=float(z),
                household_mass=float(mass), ownership=ratio(owner,mass), realized_rooms=ratio(rooms,mass)))
    table('supplemental_lifecycle.csv', lifecycle)
    table('supplemental_realized_by_wealth.csv', wealth_rows)
    table('supplemental_realized_by_income.csv', income_rows)
    moments = observer['housing_wealth']['moments']
    mass = sum(r['household_mass'] for r in lifecycle)
    check_close(sum(r['realized_rooms']*r['household_mass'] for r in lifecycle)/mass,
                moments['aggregate_mean_occupied_rooms_ahs_uncapped_18_85'], 'mean_rooms')
    stats = type('Stats', (), {})()
    model.add_aggregate_wealth_bequest_flow_moments(stats, post, current, p.bp_pol, P, bg, p.price)
    check_close(sum(r['beginning_net_worth']*r['household_mass'] for r in lifecycle),
                stats.aggregate_wealth, 'beginning_wealth_stock')
    check_close(stats.aggregate_wealth_to_annual_gross_labor_earnings,
                moments['aggregate_wealth_to_annual_gross_labor_earnings'], 'wealth_earnings')

    # Exact one-market origin conditional expectation: interpolate the actual
    # transaction map to destination policies, then average menu probabilities.
    # Renter origins cannot be owner stayers. Birth/no-birth branches are shown
    # separately; they are not unconditional before-birth policy curves.
    conditional = []
    for age in (30.,42.):
        j = int(np.flatnonzero(ages == age)[0])
        for zz, z in enumerate(P.z_grid):
            origin_mass = gp[:,0,0,j,zz,0,0]
            birth_pr = fec[j]*p.fert_probs[:,0,0,j,zz,1]
            for family, n, m, branch_pr in (('no_birth',0,0,1-birth_pr), ('first_birth',1,1,birth_pr)):
                # Native realization casts stored menus to float64 BEFORE
                # summation/normalization; float32 normalization changes choices.
                menu = np.asarray(p.tenure_probs[:,0,0,j,zz,n,m,:], dtype=float)
                menu_sum = menu.sum(axis=1)
                weights = np.divide(menu,menu_sum[:,None],out=np.zeros_like(menu),where=menu_sum[:,None]>0)
                values = {k:np.zeros(len(bg)) for k in ('rooms','consumption','saving','buyer_rooms','buyer_consumption')}
                for tn in range(len(houses)):
                    idx = p.maps.tmx_idx[0,0,tn,n,m,:]
                    wt = p.maps.tmx_wt[0,0,tn,n,m,:]
                    prob = weights[:,tn]
                    def mapped(array):
                        v = np.asarray(array[:,tn,0,j,zz,n,m])
                        return (1-wt)*v[idx]+wt*v[idx+1]
                    hh = np.full(len(bg),houses[tn]) if tn else mapped(p.hR_pol)
                    cc, bb = mapped(p.c_pol), mapped(p.bp_pol)
                    for key, val in (('rooms',hh),('consumption',cc),('saving',bb)):
                        values[key] += np.where(prob>0,prob*val,0.)
                    if tn:
                        values['buyer_rooms'] += np.where(prob>0,prob*hh,0.)
                        values['buyer_consumption'] += np.where(prob>0,prob*cc,0.)
                owner_pr = weights[:,1:].sum(axis=1)
                for b, wealth in enumerate(bg):
                    valid = bool(np.isfinite(p.V[b,0,0,j,zz,n,m]) and p.V[b,0,0,j,zz,n,m]>-1e9 and menu_sum[b]>0)
                    conditional.append(dict(age_left=age, income_index=zz,income_multiplier=float(z),
                        net_financial_wealth=float(wealth), origin_pre_birth_mass=float(origin_mass[b]),
                        family_branch=family, family_branch_probability=float(branch_pr[b]),
                        branch_mass=float(origin_mass[b]*branch_pr[b]), feasible=valid,
                        owner_probability=float(owner_pr[b]) if valid else None,
                        expected_rooms=float(values['rooms'][b]) if valid else None,
                        expected_consumption=float(values['consumption'][b]) if valid else None,
                        expected_saving=float(values['saving'][b]) if valid else None,
                        rooms_conditional_on_buying=ratio(values['buyer_rooms'][b],owner_pr[b]) if valid else None,
                        consumption_conditional_on_buying=ratio(values['buyer_consumption'][b],owner_pr[b]) if valid else None,
                        owner_stayer_probability=0.))
    table('supplemental_childless_renter_choices.csv', conditional)
    # Independent linear replay of just the initially childless renter cohort:
    # no saving, aging, KFE or solve. It tests interpolation and choice weighting.
    for age in (30.,42.):
        j = int(np.flatnonzero(ages==age)[0])
        cohort = np.zeros_like(post[:,:,:,j])
        for zz in range(P.Nz):
            mass0 = gp[:,0,0,j,zz,0,0]
            birth = fec[j]*p.fert_probs[:,0,0,j,zz,1]
            cohort[:,0,0,zz,0,0] = mass0*(1-birth)
            cohort[:,0,0,zz,1,1] = mass0*birth
        realized = model.realize_current_choices_markov_income(cohort,j,p.loc_probs,
            p.tenure_choice,p.tenure_probs,p.maps.lmm_idx,p.maps.lmm_wt,
            p.maps.tmx_idx,p.maps.tmx_wt,use_compiled_scatter=False)
        rows = [r for r in conditional if r['age_left']==age and r['feasible']]
        check_close(sum(r['branch_mass'] for r in rows),realized.sum(),f'childless_renter_{age:g}_mass')
        rooms = sum(realized[:,t].sum()*houses[t] for t in range(1,len(houses)))
        rooms += np.sum(np.where(realized[:,0]>0,realized[:,0]*p.hR_pol[:,0,:,j],0.))
        for key,actual in (('owner_probability',realized[:,1:].sum()),('expected_rooms',rooms),
                ('expected_consumption',np.sum(realized*p.c_pol[:,:,:,j])),
                ('expected_saving',np.sum(realized*p.bp_pol[:,:,:,j]))):
            check_close(sum(r['branch_mass']*r[key] for r in rows),actual,
                        f'childless_renter_{age:g}_{key}')

    fig, axes = plt.subplots(1,2,figsize=(11,4))
    for n in range(3):
        rows = [r for r in birth_rows if r['birth_order']==n+1]
        axes[0].plot(ages,[r['births_per_living_age_household'] for r in rows],label=f'Birth {n+1}')
        axes[1].plot(ages,[r['probability_per_at_risk_cell'] for r in rows],label=f'Birth {n+1}')
    axes[0].set_ylabel('Births per living household in age cell')
    axes[1].set_ylabel('Birth probability per at-risk household / 4-year cell')
    for ax in axes:
        ax.set_xlim(18,46); ax.set_xlabel('Age at left of model cell'); ax.legend(); ax.grid(alpha=.2)
    savefig(fig,'supplemental_birth_flows.png')
    fig, axes = plt.subplots(1,3,figsize=(13,4))
    for key,label in (('beginning_net_financial_wealth','Net financial wealth'),
                      ('beginning_gross_housing_value','Gross house value'),('beginning_net_worth','Net worth')):
        axes[0].plot(ages,[r[key] for r in lifecycle],label=label)
    axes[0].set_ylabel('Beginning balance sheet: mean per household');axes[0].legend(fontsize=8)
    axes[1].plot(ages,[r['ownership'] for r in lifecycle],label='All owners')
    axes[1].plot(ages,[r['owner_stayer_mass']/r['household_mass'] for r in lifecycle],label='Owner stayers')
    axes[1].plot(ages,[r['buyer_mass']/r['household_mass'] for r in lifecycle],label='Buyers')
    axes[1].legend(fontsize=8);axes[1].set_ylabel('Share of current age-cell households')
    axes[2].plot(ages,[r['realized_rooms'] for r in lifecycle]);axes[2].set_ylabel('Realized rooms / household')
    for ax in axes:
        ax.axvline(ages[P.J_R],color='gray',ls=':'); ax.set_xlabel('Age at left of model cell'); ax.grid(alpha=.2)
    savefig(fig,'supplemental_lifecycle_balance_sheets.png')
    fig, axes = plt.subplots(2,2,figsize=(11,7))
    for col,age in enumerate((30.,42.)):
        rows = [r for r in wealth_rows if r['age_left']==age and r['inside_q005_q995']]
        x=[r['net_financial_wealth'] for r in rows]
        axes[0,col].plot(x,[r['ownership'] for r in rows]); axes[0,col].set_ylabel('Realized owner share')
        axes[1,col].plot(x,[r['realized_rooms'] for r in rows]); axes[1,col].set_ylabel('Realized rooms / household')
        for ax in axes[:,col]:
            ax.set_title(f'Age {age:g}; 0.5–99.5% occupied wealth support');ax.set_xlabel('Post-transaction net financial wealth');ax.grid(alpha=.2)
    savefig(fig,'supplemental_occupied_wealth_shapes.png')
    fig, axes = plt.subplots(2,2,figsize=(11,7))
    for col,age in enumerate((30.,42.)):
        selected=[r for r in conditional if r['age_left']==age]
        j=int(np.flatnonzero(ages==age)[0]);mass_b=gp[:,0,0,j,:,0,0].sum(axis=1)
        lo,hi=quantiles(mass_b,(.005,.995))
        for family,style in (('no_birth','-'),('first_birth','--')):
            xx, own, hh=[],[],[]
            for b,wealth in enumerate(bg):
                if not lo<=wealth<=hi:continue
                rows=[r for r in selected if r['net_financial_wealth']==wealth and r['family_branch']==family and r['feasible']]
                denom=sum(r['branch_mass'] for r in rows)
                xx.append(wealth);own.append(ratio(sum(r['branch_mass']*r['owner_probability'] for r in rows),denom))
                hh.append(ratio(sum(r['branch_mass']*r['expected_rooms'] for r in rows),denom))
            axes[0,col].plot(xx,own,style,label=family.replace('_',' '));axes[1,col].plot(xx,hh,style,label=family.replace('_',' '))
        axes[0,col].set_ylabel('Owner probability conditional on birth outcome')
        axes[1,col].set_ylabel('Expected rooms over realized tenure menu')
        for ax in axes[:,col]:
            ax.set_title(f'Initially childless renters, age {age:g}');ax.set_xlabel('Pre-birth net financial wealth');ax.grid(alpha=.2);ax.legend(fontsize=8)
    savefig(fig,'supplemental_childless_renter_occupied_choices.png')
    # Recheck original small files/17 plots after writing supplements.
    for name,digest in hashes.items():
        assert sha(SOURCE/name)==digest,name
    write('supplemental_receipt.json',dict(label=LABEL,status='PASS',model_solves=0,
        economic_changes=[],source_changes=[],slurm_job=os.environ['SLURM_JOB_ID'],
        elapsed_seconds=time.monotonic()-started,script_sha256=sha(__file__),
        checkpoint=dict(path=str(checkpoint),sha256=CHECKPOINT_SHA),
        reference_manifest=dict(path=str(manifest_path),sha256=sha(manifest_path)),
        contract_sha256=CONTRACT_SHA,source_manifest=c['source_manifest'],
        scientific_identity=receipt['scientific_identity'],target_weight_fingerprint=receipt['target_weight_fingerprint'],
        all_parameters_pinned_by_checkpoint=True,common_primary_loss=loss,checks=checks,
        standard_plots_unchanged=17,source_fit_table=str(SOURCE/'target_fit.csv'),
        source_parameter_table=str(SOURCE/'parameters.csv'),
        definitions=dict(state_axes=['net financial wealth','tenure/product','location','age','income','children ever born','children at home'],
            births='Actual explicit birth flows during the four-year cell; origin characteristics are pre-birth. Third births do not receive top-bin demographic weights.',
            birth_quintiles='Whole wealth nodes assigned by midpoint weighted CDF of all pre-birth households within age; ties mean bins need not contain exactly 20%.',
            wealth_quantiles='Left inverse of occupied current-wealth CDF; no interpolation.',
            buyer='Current owners minus exact owner-stayer mass; includes owner-to-different-product purchases.',
            realized_consumption_saving='Integrates disjoint date-owned buyer/renter and owner-stayer masses using policy_mass_branches, never buyer policies for owner stayers.',
            conditional='Pre-birth childless renters; no-birth/first-birth branches and actual stochastic tenure menu. Buyer-conditional columns use positive purchase probabilities only.',
            balance_sheet='Beginning net worth b + pH uses post-fertility pre-transaction mass, matching authoritative wealth observer; current balance sheets shown separately.',
            housing_equity='Unavailable separately: checkpoint stores net financial assets b and gross housing value pH, not independently identified gross financial assets and mortgage debt. No fabricated equity series.',
            age='Model age-cell left endpoints, not interpolated annual ages or observed cohort paths.'),
        inherited_measurement_audit=str(BASE/'measurement_audit_v1/README.md'),
        limitations=['Descriptive saved-state accounting does not establish causal mechanisms.',
            'No housing-equity decomposition without an additional balance-sheet convention.',
            'No new transition, equilibrium counterfactual, parameter change or fertility renormalization.',
            'Standard plots retained; all new views are supplemental.']))
    print(json.dumps(dict(status='PASS',output=str(out),checks=checks,elapsed_seconds=time.monotonic()-started)),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True,help='New supplemental result directory; refuses overwrite')
    main(parser.parse_args().output)
