"""Saved 2007 household policies: extract on Torch; render from CSV locally."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
import numpy as np


def extract(spec_path, out):
    spec = json.loads(Path(spec_path).read_text())
    sys.path.insert(0, str(Path(spec['batch']) / 'source'))
    import run_e5f_original_queue_experiments as runner
    c = runner.load_context(spec_path)  # Frozen inputs, no Bellman or GE solve.
    e, P = c.packet['evaluation'], c.old.parameters
    p, model = e.policy, c.primitive.calendar.model
    grid = np.asarray(c.old.b_grid, float)
    assert P.I == 1 and P.sequential_births and P.child_state_mode == 'independent_count'
    assert p.joint_choice is None
    j = int(np.argmin(abs(P.age_start + P.da * np.arange(P.J) - 30)))
    settled = int(model.readiness_settled_state(P))
    fec = float(model.get_fecundity_by_age(P)[j])
    z = np.asarray(P.z_grid, float)
    zw = np.asarray(P.z_weights, float)
    groups = np.asarray(P.permanent_income_group_index)
    group_values = np.unique(groups)
    group = group_values[np.argmax([zw[groups == k].sum() for k in group_values])]
    indices = np.flatnonzero(groups == group)
    indices = indices[np.argsort(z[indices])]
    cum = np.cumsum(zw[indices]) / zw[indices].sum()
    low, med, high = [int(indices[np.searchsorted(cum, v)]) for v in (.25, .5, .75)]
    assert len({low, med, high}) == 3
    births = {key: fec * np.asarray(p.fert_probs[:, 0, 0, j, zz, 1], float)
              for key, zz in [('low', low), ('middle', med), ('high', high)]}
    def housing(nn, cs):
        probs = np.asarray(p.tenure_probs[:, 0, 0, j, med, nn, cs, :], float)
        h = np.asarray(p.hR_pol[:, 0, 0, j, med, nn, cs], float)
        rooms = probs[:, 0] * h + probs[:, 1:] @ np.asarray(P.H_own)
        # Independent explicit tenure-by-tenure sum for the plotted mixture.
        explicit = probs[:, 0] * h
        for t, hh in enumerate(P.H_own, 1):
            explicit = explicit + probs[:, t] * hh
        np.testing.assert_allclose(rooms, explicit, rtol=0, atol=1e-12)
        return rooms, probs
    h0, tp0 = housing(0, settled)
    h1, tp1 = housing(1, 1)
    valid = grid >= 0
    for zz in (low, med, high):
        valid &= np.asarray(p.V[:, 0, 0, j, zz, 0, settled]) > -1e9
    valid &= (tp0.sum(1) > .999999) & (tp1.sum(1) > .999999)
    prob_error = float(max(abs(tp0[valid].sum(1)-1).max(), abs(tp1[valid].sum(1)-1).max()))
    assert prob_error < 2e-7  # Saved tenure probabilities are float32.
    assert all(np.all((v[valid] >= 0) & (v[valid] <= 1)) for v in births.values())
    assert not any(bool(getattr(P, k, False)) for k in
                   ['birth_dp_grant', 'birth_entry_grant', 'parent_dp_waiver', 'use_pti_constraint'])
    phi = np.asarray(P.phi, float)
    assert np.all(phi == phi[0])
    q = float(p.price[0])
    thresholds = (1-phi[0]) * q * np.asarray(P.H_own, float)
    for probs in (tp0, tp1):
        for t, threshold in enumerate(thresholds, 1):
            assert np.max(probs[(grid < threshold) & valid, t], initial=0) < 1e-12
    mass = np.asarray(e.g_pre[:, 0, 0, j, :, 0, settled])[:, indices].sum(1)
    mass = np.where(valid, mass, 0)
    assert mass.sum() > 0
    p99 = float(grid[np.searchsorted(np.cumsum(mass) / mass.sum(), .99)])
    xmax = max(p99, float(thresholds[2] * 1.1))
    xmax = min(xmax, float(grid[valid].max()))
    rows = [dict(liquid_wealth=float(b), valid=int(valid[k]),
                 first_birth_low=float(births['low'][k]), first_birth_middle=float(births['middle'][k]),
                 first_birth_high=float(births['high'][k]), housing_no_birth=float(h0[k]),
                 housing_birth=float(h1[k]), ownership_no_birth=float(tp0[k, 1:].sum()),
                 ownership_birth=float(tp1[k, 1:].sum()), slice_mass=float(mass[k]))
            for k, b in enumerate(grid)]
    out.mkdir(parents=True, exist_ok=True)
    with (out/'household_mechanism_2007.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    summary = json.loads(Path(c.manifest['initial_summary']).read_text())['checkpoint']
    meta = dict(status='PASS', no_new_solve=True, age=float(P.age_start + j*P.da),
                age_index=j, settled_state=settled, permanent_income_group=int(group),
                income_indices=[low, med, high], income_multipliers=z[[low, med, high]].tolist(),
                income_relative_to_middle=(z[[low, med, high]]/z[med]).tolist(), fecundity=fec,
                owner_sizes=np.asarray(P.H_own).tolist(), asset_price=q, financed_share=float(phi[0]),
                downpayment_thresholds=thresholds.tolist(), saved_probability_sum_error=prob_error,
                xlim=[0, xmax], slice_wealth_p99=p99, nonnegative_occupied_mass_retained=float(mass[grid<=xmax].sum()/mass.sum()),
                snapshot=summary, spec_sha256=hashlib.sha256(Path(spec_path).read_bytes()).hexdigest(),
                model_source=str(model.__file__),
                birth_formula='fecundity[j] * fert_probs[b,0,0,j,z,1] among settled childless renters',
                housing_formula='tenure_probs[...,0]*hR_pol[...,n,cs] + sum(owner_probabilities*owner_rooms)',
                branch_states={'no_birth':[0,settled], 'birth':[1,1]},
                financing_marker='Purchase eligibility: (1-phi)*price*owner_rooms, not proof of a binding desired borrowing constraint',
                selection='Age nearest30, initial renter, sole location, modal permanent-income group; 25/50/75 percentiles of its income-state weights',
                full_grid=[float(grid.min()),float(grid.max())])
    (out/'verification.json').write_text(json.dumps(meta, indent=2)+'\n')
    print(json.dumps({k:meta[k] for k in ['status','age','income_indices','income_relative_to_middle','xlim','downpayment_thresholds']}), flush=True)


def render(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    m = json.loads((out/'verification.json').read_text())
    assert m['status']=='PASS'
    data = np.genfromtxt(out/'household_mechanism_2007.csv', delimiter=',', names=True)
    keep = data['valid'].astype(bool) & (data['liquid_wealth'] <= m['xlim'][1])
    x = data['liquid_wealth'][keep]
    plt.rcParams.update({'font.size':12, 'axes.spines.top':False, 'axes.spines.right':False})
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.0))
    fig.subplots_adjust(left=.07, right=.985, top=.78, bottom=.18, wspace=.25)
    colors = ['#2968a2','#a7a7a7','#d27328']
    for i, key in enumerate(['low','middle','high']):
        ratio = m['income_relative_to_middle'][i]
        label = f'{key.capitalize()} income ({ratio:.2f}'+r'$\times$'+')'
        ax[0].plot(x, 100*data['first_birth_'+key][keep], label=label, color=colors[i], lw=2.3,
                   ls='--' if key=='middle' else '-')
    ax[0].set(title='First-birth probability', ylabel='Probability over four years (%)')
    ax[0].legend(frameon=False, fontsize=10, loc='best')
    ax[1].plot(x, data['housing_no_birth'][keep], color=colors[0], lw=2.3, label='Without a birth')
    ax[1].plot(x, data['housing_birth'][keep], color=colors[2], lw=2.3, label='With a birth')
    ax[1].set(title='Housing chosen at middle income', ylabel='Housing services (rooms)')
    ax[1].legend(frameon=False, fontsize=10, loc='best')
    for a in ax:
        a.set_xlim(m['xlim']); a.set_xlabel('Liquid wealth (model units)')
        a.set_title(a.get_title(), pad=23)
        a.axvspan(0, m['downpayment_thresholds'][0], color='#ddd6c9', alpha=.4, zorder=-10)
        a.grid(axis='y', color='.9', lw=.6)
        for size, threshold in zip(m['owner_sizes'],m['downpayment_thresholds']):
            if threshold < m['xlim'][1]:
                a.axvline(threshold, color='.6', alpha=.6, ls=':', lw=.8)
                a.text(threshold, 1.02, f'{size:g}', ha='center', transform=a.get_xaxis_transform(), color='.45', fontsize=8)
    fig.suptitle('Fertility, housing, and down payments', fontsize=18, y=.98)
    fig.text(.5,.90,f'2007 steady state · Age {m["age"]:g} · Initially childless renters · Permanent income held fixed',ha='center',fontsize=11)
    fig.text(.5,.035,'Dotted lines: cash needed to buy homes with 2, 4, 6, 8, or 10 rooms. Shading: no owner product affordable.',ha='center',fontsize=9,color='.35')
    fig.savefig(out/'household_mechanism_2007.png',dpi=200)
    fig.savefig(out/'household_mechanism_2007.pdf')
    plt.close(fig)
    (out/'README.md').write_text('''# 2007 household mechanism\n\nSupplemental saved-policy figure; no model solve or calibration change. The CSV retains the full wealth grid; the figure displays valid nonnegative renter wealth through at least the slice\'s 99th percentile, extending to the six-room down payment if needed. Quantile and source details are in verification.json.\n\nBirth probability is fecundity times the saved first-birth attempt probability. Housing averages the actual tenure choice conditional on a realized birth or no birth from the same initial childless state. Current income varies within the same permanent-income group. Dotted lines indicate down-payment eligibility, not proof of a binding desired borrowing constraint. This is a state-conditioned policy comparison, not an income-fertility cross-sectional regression. Curves join the existing grid points without smoothing.\n\nRe-render with:\n```\npython code/model/tools/build_e5f_ss2007_household_mechanism.py --render-only --outdir output/model/e5f_original_queue_20260913a/household_mechanism_2007\n```\n''')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--spec');ap.add_argument('--outdir',type=Path,required=True);ap.add_argument('--render-only',action='store_true');ap.add_argument('--extract-only',action='store_true')
    a=ap.parse_args()
    if not a.render_only:
        if not a.spec: ap.error('--spec required for extraction')
        extract(a.spec,a.outdir)
    if not a.extract_only: render(a.outdir)

if __name__=='__main__': main()
