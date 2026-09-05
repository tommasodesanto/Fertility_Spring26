#!/usr/bin/env python3
"""Build a concise comparison from completed, collected policy results."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
DEFAULT = ROOT / 'output/model/e5f_candidate_policy_comparison_20260905a'
OLD = ROOT / 'output/model/e5f_post2023_policy_mechanisms_eta063_kappa005_task010_production_20260904a'
OLD_TAX = ROOT / 'output/model/e5f_rebated_tax_shapley_diagnosis_eta063_kappa005_task010_20260904a'
LABELS = {'supply-plus-20':'Supply +20%', 'dependent-child-ltv95':'Family LTV 95%',
          'combined':'Supply and family credit', 'property-tax-2pct-no-rebate':'Tax 1% to 2%, no rebate'}


def read(path):
    return json.loads(Path(path).read_text())


def rows(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, data):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data[0]))
        writer.writeheader(); writer.writerows(data)


def effects(folder, label):
    summary = read(folder / 'report/summary.json')
    assert summary['status'] == 'complete_post2023_policy_mechanism_comparison'
    assert summary['post_2023_periods'] == 10
    paths = rows(folder / 'report/policy_paths.csv')
    base = {int(r['calendar_year']):r for r in paths if r['policy_case']=='baseline'}
    assert sorted(base) == list(range(2023,2064,4))
    results = []
    for row in paths:
        year = int(row['calendar_year']); b = base[year]
        pct = lambda field: 100*(float(row[field])/float(b[field])-1)
        item = dict(calibration=label, policy=row['policy_case'], year=year,
            births_per_household_percent=pct('topcode_adjusted_births_per_adult'),
            total_births_percent=pct('birth_children_topcode_adjusted'),
            household_mass_percent=pct('adult_population'), rooms_percent=pct('housing_demand_per_adult'),
            ownership_pp=100*(float(row['owner_rate'])-float(b['owner_rate'])),
            dependent_child_ownership_pp=100*(float(row['dependent_child_owner_rate'])-float(b['dependent_child_owner_rate'])),
            house_price_percent=pct('asset_price'), housing_user_cost_percent=pct('housing_user_cost'))
        identity=100*((1+item['births_per_household_percent']/100)*(1+item['household_mass_percent']/100)-1)
        assert abs(identity-item['total_births_percent']) < 2e-10
        if year < 2043:
            assert abs(item['household_mass_percent']) < 2e-10
        results.append(item)
    assert len(results)==55
    return results, paths


def table(header, data):
    return '\n'.join(['| '+' | '.join(header)+' |', '| '+' | '.join(['---']*len(header))+' |',
        *['| '+' | '.join(map(str,r))+' |' for r in data]])


def build(out, source):
    contract_fields = ('policy', 'policy_start_year', 'population_closure',
                       'post_2023_periods', 'renewal_lag', 'scope', 'terminal_year')
    for case in ('baseline', *LABELS):
        before, after = read(OLD/case/'summary.json'), read(out/'full'/case/'summary.json')
        for field in contract_fields:
            assert before[field] == after[field], (case, field)
        for field in ('source_sha256', 'target_fingerprint', 'renewal_contract_sha256'):
            assert before['input_hashes'][field] == after['input_hashes'][field], (case, field)
        # The only policy-driver diff adds three reported population columns.
        # The checked Git diff is saved with the historical reconciliation.
        assert before['driver_sha256'] == '628edeada00a711b5fabe54a9f2f5d6b0f857a5566b45b87431e0514a3383be6'
        assert after['driver_sha256'] == '345f2b54718bde99d11ee490e649bf005cc6b20125b6a7dbdb0ab698abf6b2f8'
    assert read(out/'collection_receipt.json')['status'] == 'pass'
    old, old_paths = effects(OLD, 'Current benchmark')
    new, new_paths = effects(out/'full', 'Verified candidate')
    all_rows = old + new
    write_csv(out/'policy_comparison.csv', all_rows)
    index = {(r['calibration'],r['policy'],r['year']):r for r in all_rows}
    getter = lambda label,case,year: index[label,case,year]
    comparisons=[]
    for case in LABELS:
        a,b=getter('Current benchmark',case,2023),getter('Verified candidate',case,2023)
        c,d=getter('Current benchmark',case,2063),getter('Verified candidate',case,2063)
        comparisons.append([LABELS[case], f"{a['births_per_household_percent']:+.3f}", f"{b['births_per_household_percent']:+.3f}",
            f"{c['births_per_household_percent']:+.3f}", f"{d['births_per_household_percent']:+.3f}"])
    impact=[]
    for case in LABELS:
        r=getter('Verified candidate',case,2023)
        impact.append([LABELS[case],f"{r['housing_user_cost_percent']:+.3f}%",f"{r['rooms_percent']:+.3f}%",f"{r['ownership_pp']:+.3f} pp",f"{r['births_per_household_percent']:+.3f}%"])
    total=[]
    for case in LABELS:
        r=getter('Verified candidate',case,2063)
        total.append([LABELS[case], f"{r['births_per_household_percent']:+.3f}%", f"{r['household_mass_percent']:+.3f}%", f"{r['total_births_percent']:+.3f}%"])
    fig,axes=plt.subplots(2,2,figsize=(9.4,5.8),constrained_layout=True)
    for ax,case in zip(axes.flat,LABELS):
        for label,color,style in [('Current benchmark','#888888','--'),('Verified candidate','#1f587a','-')]:
            rr=[r for r in all_rows if r['calibration']==label and r['policy']==case]
            ax.plot([r['year'] for r in rr],[r['births_per_household_percent'] for r in rr],color=color,linestyle=style,label=label)
        ax.set_title(LABELS[case],fontsize=11);ax.set_ylabel('Births per household: change (%)',fontsize=9)
        ax.set_xticks([2023,2043,2063]);ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    fig.savefig(out/'policy_fertility_comparison.png',dpi=180);fig.savefig(out/'policy_fertility_comparison.pdf');plt.close(fig)
    tax_rows=[]
    tax_values={}
    for label,folder in [('Current benchmark',OLD_TAX),('Verified candidate',out/'shapley')]:
        assert read(folder/'summary.json')['status']=='complete'
        rr=[r for r in rows(folder/'shapley_decomposition.csv') if r['metric']=='births_per_adult']
        assert len(rr)==3
        tax_values[label]={r['factor']:float(r['reported_contribution']) for r in rr}
    for name,key in [('Direct tax','tax_rate'),('House-price capitalization','asset_price'),('Equal rebate','equal_rebate')]:
        tax_rows.append([name,*[f"{tax_values[label][key]:+.3f}%" for label in tax_values]])
    tax_rows.append(['Net',*[f"{sum(tax_values[label].values()):+.3f}%" for label in tax_values]])
    diagnostic_receipts=[read(p) for p in (out/'full_diagnostics').glob('*/*/receipt.json')]
    assert len(diagnostic_receipts)==55
    assert all(r['status']=='complete' for r in diagnostic_receipts)
    max_exposure=max(r['policy_arrays']['share_pre_choice_mass_at_negative_steps'] for r in diagnostic_receipts)
    flagged=sum(r['policy_arrays']['occupied_negative_steps']>0 for r in diagnostic_receipts)
    max_budget=max(r['budget']['budget_excess_share'] for r in diagnostic_receipts)
    tax_receipts=[read(p) for p in (out/'rebate_standard_diagnostics').glob('*/receipt.json')]
    assert len(tax_receipts)==11 and all(r['status']=='complete' for r in tax_receipts)
    tax_flagged=sum(r['policy_arrays']['occupied_negative_steps']>0 for r in tax_receipts)
    max_tax_exposure=max(r['policy_arrays']['share_pre_choice_mass_at_negative_steps'] for r in tax_receipts)
    numerical=read(out/'saving_flag_check/summary.json')
    assert numerical['status']=='complete'
    max_fixed_birth=max(abs(r['births_difference_percent']) for r in numerical['cases'])
    original_fit=read(ROOT/'output/model/e5f_transition_calibration_eta063_kappa005_rooms_repair_gauss_newton_coord_20260904a/task_010/summary.json')['best_candidate']['transition_loss']
    candidate_fit=read(out/'plan.json')['candidate_loss']
    max_change=max(abs(getter('Verified candidate',case,year)['births_per_household_percent']-getter('Current benchmark',case,year)['births_per_household_percent']) for case in LABELS for year in [2023,2063])
    sign_same=all(getter('Verified candidate',case,year)['births_per_household_percent']*getter('Current benchmark',case,year)['births_per_household_percent']>0 for case in LABELS for year in range(2023,2064,4))
    for label in ('Current benchmark', 'Verified candidate'):
        for year in range(2023,2064,4):
            ranking = ('property-tax-2pct-no-rebate', 'dependent-child-ltv95', 'supply-plus-20', 'combined')
            values = [getter(label,case,year)['births_per_household_percent'] for case in ranking]
            assert values == sorted(values), (label, year, values)
    if not sign_same:
        raise RuntimeError('Policy signs changed; rewrite the interpretation before producing a report')
    summary=dict(status='complete', old_loss=original_fit, candidate_loss=candidate_fit,
        loss_reduction_percent=100*(1-candidate_fit/original_fit), same_policy_signs=sign_same,
        maximum_birth_response_change_pp=max_change, flagged_policy_dates=flagged,
        maximum_flagged_mass_share=max_exposure, maximum_reporting_budget_excess_share=max_budget,
        flagged_tax_equilibria_or_cells=tax_flagged, maximum_tax_flagged_mass_share=max_tax_exposure,
        maximum_fixed_price_optimizer_birth_change_percent=max_fixed_birth,
        production_promoted=False, policy_rows=len(all_rows), standard_policy_graphs=55*17)
    (out/'comparison_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    fit_note=(ROOT/'docs/model/e5f_first_birth_measurement_review.md').read_text()
    appendix=fit_note.split('# Complete target fit: verified candidate',1)[1].split('# Evidence and reproduction',1)[0]
    base2023=next(r for r in new_paths if r['policy_case']=='baseline' and r['calendar_year']=='2023')
    example_births=10000*float(base2023['topcode_adjusted_births_per_adult'])
    supply=getter('Verified candidate','supply-plus-20',2023)['births_per_household_percent']
    credit=getter('Verified candidate','dependent-child-ltv95',2023)['births_per_household_percent']
    source.write_text(f'''---
title: "Calibration and policy: a clear comparison"
subtitle: "Decision readout for the September 14 presentation"
author: "Prepared for Tommaso De Santo"
date: "5 September 2026"
---

# The result to take into our discussion

**The verified candidate improves calibration fit and preserves the policy ordering.**
The unchanged twelve-target loss falls from **{original_fit:.3f} to {candidate_fit:.3f}**,
a **{100*(1-candidate_fit/original_fit):.2f}%** improvement. The largest change in
the reported impact or 2063 birth response is **{max_change:.3f} percentage points**.
The candidate remains a review point; the retained production calibration has
not been overwritten.

Supply expansion increases housing use and births. Family mortgage credit
increases ownership with little change in housing use or births. An unrebated
property-tax increase reduces births. Equal tax rebates alter that conclusion
because they change household resources. The fiscal convention must accompany
every tax claim.

**The main calibration weakness remains visible:** first-birth rooms are
**0.466 versus 0.720** in the data; mean rooms are **6.329 versus 5.780**.
The completed empirical check supports retaining the 0.720 target. A better
objective value does not resolve these housing discrepancies or establish that
the policy fertility responses are directly identified by causal cost variation.

My recommendation is to use this candidate as the working presentation point
after our review, with these limitations explicit. Keep the current benchmark
as the documented comparison. Another broad calibration search is not the
immediate priority.

# How much the policy numbers change

Entries below are percent changes in **births per household** against each
calibration's own no-policy path. 2063 is a finite transition date.

{table(['Policy','Current 2023','Candidate 2023','Current 2063','Candidate 2063'],comparisons)}

\\clearpage

# What moves on impact

{table(['Candidate policy','Housing user cost','Rooms / HH','Ownership','Births / HH'],impact)}

The supply experiment shifts the supply schedule by 20% at every given price.
Equilibrium prices adjust, so the resulting increase in occupied rooms is the
smaller amount reported above. Family credit applies to households with
dependent children, not only first-time buyers.

For a concrete scale, 10,000 model households have about **{example_births:.1f} births**
in the baseline four-year period. Supply expansion adds about
**{example_births*supply/100:.1f} births**; family credit adds about
**{example_births*credit/100:.2f}**. This is an illustration of model units,
not an annual national forecast.

![Supplemental comparison of the complete policy paths. Dashed grey is the retained benchmark; solid blue is the candidate. The standard seventeen model diagnostics are preserved separately for every candidate policy date.]({(out/'policy_fertility_comparison.png').relative_to(ROOT)}){{width=100%}}

\\clearpage

# Tax results depend on the rebate comparison

The following three-factor Shapley decomposition averages each channel's
contribution across all six orders in which tax, house price and rebate can
change. Entries are percent of births per household in the **rebated 1%**
baseline; the reform is **rebated 2%**. Both equilibria use the same inherited
2023 distribution. The eight component cells reproduce the endpoints and add up.

{table(['Birth contribution','Current benchmark','Verified candidate'],tax_rows)}

This net effect is not the effect relative to the unrebated status quo. The
five-policy table instead compares unrebated 2% with unrebated 1%. The earlier
September 4 chat explicitly drew this distinction; the present review preserves it.

# Birth rates and total births in 2063

{table(['Candidate policy','Births / HH','Household mass','Total births'],total)}

Total births combine the birth rate and the number of households. The exact
identity is checked for every row. Household mass is identical across policies
before 2043, as required by the maintained twenty-year entry lag. The inherited
household queue is a conditional mechanism closure; this table is not a
validated resident-population forecast or a terminal steady-state comparison.

\\clearpage

# What was checked, and what remains open

All five short transition tests passed before the full five paths were launched.
Every full path exactly reproduces its short test's first two dates. The historical
reconstruction, shared starting state, market clearing, population, feasibility,
probability and accounting checks pass at the existing tolerances. Each of the
55 policy dates has a saved state and the unchanged seventeen-graph packet.

The source contains the already verified policy-reuse correction. Its fixed-
parameter history reproduces all twelve target rows, every parameter and bound,
and 253 numeric historical entries of the original candidate exactly. No model
primitive, target, weight, parameter or numerical tolerance was changed here.

The numerical screen asks whether household value ever falls when wealth rises.
It flags **{flagged} of 55 dates**, but the largest affected share is only
**{1e6*max_exposure:.2f} per million** of the inherited household distribution.
For three flagged policy-date distributions, a six-solve check first reproduces
the saved policies exactly, then searches the full saving grid at the same prices.
This removes the flagged value decreases in all three tests. The largest absolute
change in births per household is **{max_fixed_birth:.6f}%**. The eleven tax
equilibria/component cells contain {tax_flagged} flagged packets, with maximum
affected share **{1e6*max_tax_exposure:.2f} per million**. These checks do not bound
the error of a fully re-solved equilibrium path or the tax decomposition.
The largest share affected by the reporting consumption floor exceeding the
budget is **{1e6*max_budget:.3f} per million**. Both rebated equilibria also pass
the check that the selected rebate matches the parameter state used for accounting.

The old chats explain several apparent inconsistencies. Supply elasticity 1.75
initializes the old equilibrium; 0.63 applies to dated markets. This distinction
was already recorded, though its economic justification remains open. A separate,
later demographic branch tracks people more consistently; it does not change
this retained household-entry experiment. The rebated-tax result uses its own
fiscal baseline. Direct author statements support studying policy during the
transition, without settling every inherited population assumption.

For the presentation, the defensible claim is that these mechanisms survive a
meaningful improvement in fit under the same target system. The remaining work
is to decide how prominently to qualify the housing-fit shortfall, population
closure and fiscal scope. This readout does not certify a new perfect-foresight
equilibrium, welfare ranking or globally optimal calibration.

\\clearpage

# Complete target fit: verified candidate
{appendix}
# Reproduction and evidence

The evidence folder is `{out.relative_to(ROOT)}`. It contains the pinned plan,
historical reconciliation, source and input hashes, all policy rows, exact
smoke/full comparison receipts, tax decomposition, numerical-flag check and
complete fit/parameter tables. Its README states the run budgets and commands.
The new comparison CSV contains all 110 policy/year/calibration rows, including
the no-policy paths. The standard graphs are in `full_diagnostics/`, grouped by
policy and date. All paths are relative to:
`{ROOT}`.

Rebuild this readout from the collected outputs:

```sh
python3 code/model/tools/build_e5f_candidate_policy_comparison_report.py
```
''')
    print(source)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=DEFAULT)
    parser.add_argument('--source',type=Path,default=ROOT/'docs/model/e5f_candidate_policy_comparison_review.md')
    args=parser.parse_args();build(args.output_dir.resolve(),args.source.resolve())
