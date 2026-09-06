#!/usr/bin/env python3
"""Verify and plot a completed diagnostic panel without importing model code."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def collect(root):
    contract = json.loads((root/'final_contract.json').read_text())
    expected = digest(root/'final_contract.json')
    cluster = root/'final_run'
    smoke, panel = cluster/'smoke', cluster/'panel'
    checks = {}
    for label, folder, count in [('smoke', smoke, 4), ('panel', panel, 26)]:
        receipt = json.loads((folder/'summary.json').read_text())
        assert receipt['status'] == 'complete' and receipt['completed_cases'] == count
        assert receipt['contract_sha256'] == expected
        checks[f'{label}_seconds'] = receipt['elapsed_seconds']
    checks['capture_sha256'] = digest(smoke/'captured_values.npz')
    assert checks['capture_sha256'] == receipt['capture_sha256']
    replays = [json.loads((smoke/f'capture_{r}/reference_replay.json').read_text()) for r in range(2)]
    assert replays[0] == replays[1] and replays[0]['status'] == 'exact'
    assert replays[0]['maximum_deterministic_reconstruction_gap'] == 0
    checks['exact_reference_arrays_per_repeat'] = len(replays[0]['arrays'])
    checks['reference_graphs_per_repeat'] = len(list((smoke/'capture_0/reference/standard_diagnostics').glob('*.png')))
    assert checks['reference_graphs_per_repeat'] == 17
    cases = [json.loads(p.read_text()) for p in sorted(panel.glob('case_*/summary.json'))]
    assert len(cases) == 26
    assert {(r['outer_scale'], r['dissimilarity']) for r in cases} == {tuple(p) for p in contract['cases']}
    manifests = {}
    for directory in sorted(panel.glob('case_*')):
        row = json.loads((directory/'summary.json').read_text())
        assert digest(directory/'choice_state.npz') == row['choice_state_sha256']
        assert row['infeasible_mass'] <= 1e-12 and row['mass_error'] <= 2e-10
        assert row['max_probability_error'] <= 1e-12 and row['birth_parity_identity_error'] <= 2e-10
        assert row['added_budget_excess_mass'] <= 2e-10
        assert (directory/'supplemental_choices.png').is_file()
        manifests[str(directory.relative_to(root))] = row['choice_state_sha256']
    checks['choice_state_hashes'] = manifests
    matches = {key: 0.0 for key in ['g_current', 'g_post', 'attempt', 'owner_choice', 'value']}
    for i, (outer, lam) in enumerate(contract['cases']):
        if lam != 1:
            continue
        with np.load(panel/f'case_{i:02d}_joint/choice_state.npz') as joint, np.load(panel/f'case_{i:02d}_sequential/choice_state.npz') as seq:
            for key in matches:
                a, b = joint[key], seq[key]
                assert np.array_equal(np.isfinite(a), np.isfinite(b))
                finite = np.isfinite(a)
                gap = float(np.max(np.abs(a[finite]-b[finite])))
                assert gap < 1e-12, (i, key, gap)
                matches[key] = max(matches[key], gap)
    checks['flat_logit_full_array_maximum_gaps'] = matches
    # Require the same scale/rule cases to exactly reproduce smoke outputs.
    ignored = {'elapsed_seconds'}
    for p in smoke.glob('case_*/summary.json'):
        a = json.loads(p.read_text())
        b = next(r for r in cases if all(r[k] == a[k] for k in ['rule','outer_scale','dissimilarity']))
        assert {k:v for k,v in a.items() if k not in ignored} == {k:v for k,v in b.items() if k not in ignored}
    checks['smoke_panel_exact_case_matches'] = 4
    reference = json.loads((panel/'reference_aggregates.json').read_text())
    comparisons = []
    for outer, lam in contract['cases']:
        pair = {r['rule']:r for r in cases if r['outer_scale']==outer and r['dissimilarity']==lam}
        a, b = pair['joint'], pair['sequential']
        comparisons.append(dict(outer_scale=outer, dissimilarity=lam,
            joint_births_per_100_households=100*a['births_per_household'],
            sequential_births_per_100_households=100*b['births_per_household'],
            difference_births_per_100_households=100*(a['births_per_household']-b['births_per_household']),
            joint_ownership_percent=100*a['ownership'],
            difference_ownership_pp=100*(a['ownership']-b['ownership']),
            joint_rooms=a['rooms'], difference_rooms=a['rooms']-b['rooms'],
            joint_fixed_price_market_gap_percent=100*a['signed_fixed_price_market_gap']))
    with (root/'comparison.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(comparisons[0]));writer.writeheader();writer.writerows(comparisons)
    checks['maximum_absolute_rule_difference_births_per_100_households'] = max(abs(r['difference_births_per_100_households']) for r in comparisons)
    checks['maximum_absolute_rule_difference_ownership_pp'] = max(abs(r['difference_ownership_pp']) for r in comparisons)
    checks['maximum_occupied_value_decreases'] = max(r['occupied_value_decreases'] for r in cases)
    checks['maximum_added_budget_excess_mass'] = max(r['added_budget_excess_mass'] for r in cases)
    checks['reference'] = reference
    (root/'collection_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    fig,axes=plt.subplots(1,3,figsize=(11,3.5),constrained_layout=True)
    colors={.05:'#1b6a83',.5:'#b2642e',1.:'#667143'}
    for ax,(key,title,factor) in zip(axes,[('births_per_household','Birth units per 100 households',100),('ownership','Ownership (%)',100),('rooms','Rooms per household',1)]):
        for lam,color in colors.items():
            for rule,style in [('joint','-'),('sequential','--')]:
                rows=sorted([r for r in cases if r['dissimilarity']==lam and r['rule']==rule],key=lambda r:r['outer_scale'])
                ax.plot([r['outer_scale'] for r in rows],[factor*r[key] for r in rows],style,color=color,marker='o' if rule=='joint' else None,markersize=3,label=f'lambda={lam:g}' if rule=='joint' else None)
        ax.axhline(factor*reference[key],color='#555555',linestyle=':',label='Retained reference')
        ax.set(xscale='log',xlabel='Outer tenure scale (kappa)',title=title)
        ax.grid(alpha=.15)
    axes[0].legend(fontsize=8)
    fig.suptitle('One-date diagnostic: joint choice (solid), sequential control (dashed)\nFixed prices, population and baseline future values; four-year birth units include the 3+ adjustment',fontsize=10)
    fig.savefig(root/'panel_comparison.png',dpi=180);plt.close(fig)
    print(json.dumps({k:v for k,v in checks.items() if k!='choice_state_hashes'},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True)
    collect(parser.parse_args().root)
