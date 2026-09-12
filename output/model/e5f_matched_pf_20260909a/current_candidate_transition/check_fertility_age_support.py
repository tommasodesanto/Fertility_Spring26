"""Read-only fertility measurement audit; no solver or target mutations.
Run: python check_fertility_age_support.py. Uses cached official NCHS tables.
"""
from pathlib import Path
from decimal import Decimal as D
import csv
import hashlib
import json
import math
import re

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'path_pilot_20260910/fertility_data/sources'
OUT = HERE / 'age_support_check'
OUT.mkdir(exist_ok=True)
FILES = ['nvsr66-1_2015_final_p20.txt', 'nvsr74-1_2023_final_p14.txt']

def parse(name):
    rows = {}
    for line in (SOURCE / name).read_text().splitlines():
        m = re.match(r'^(20\d\d)[\s.]+([\d,]+\.\d)\s+(.*)', line)
        if not m or int(m[1]) in rows:
            continue
        rates = [D(x) for x in m[3].split()]
        assert len(rates) == 10, (name, line)
        total = D(m[2].replace(',', ''))
        assert 5 * sum(rates[i] for i in [0, 1, 4, 5, 6, 7, 8, 9]) == total
        rows[int(m[1])] = (total, rates)
        if int(m[1]) == (2000 if '66-1' in name else 2010):
            break  # National rows only; do not parse race-specific blocks.
    return rows

old, new = [parse(f) for f in FILES]
for year in range(2010, 2016):
    assert old[year] == new[year]
annual = []
for year in range(2007, 2024):
    total, r = (old if year < 2010 else new)[year]
    partial = (2*r[3] + 5*sum(r[4:9]))/1000
    young = (5*r[0] + 3*r[2])/1000
    grouping = (5*r[1] - 3*r[2] - 2*r[3])/1000
    older = 5*r[9]/1000
    assert total/1000 - partial == young + grouping + older
    annual.append(dict(year=year, published_tfr=float(total/1000), partial_18_44=float(partial),
        approximate_18_45=float(partial+r[9]/1000), young_component=float(young),
        teen_grouping_difference=float(grouping), published_older_component=float(older)))

models = {}
for case in ['delta_m0025', 'delta_m005', 'delta_m010']:
    p = HERE / 'return_home_20260911/paths' / case
    root = json.loads((p/'root_receipt.json').read_text())
    mappings = json.loads((p/'dated_household_fertility.json').read_text())['mappings']
    final = mappings[root['final']['payload']['trial']-1]
    assert final == mappings[root['best']['payload']['trial']-1]
    models[case] = {}
    for row in final:
        d = row['diagnostics']
        flows, masses = d['birth_flow_topcode_adjusted'], d['age_mass']
        rates = [f/m if m > 0 else 0 for f,m in zip(flows, masses)]
        assert all(math.isclose(a,b,abs_tol=1e-12) for a,b in zip(rates,d['age_specific_birth_rate_topcode_adjusted']))
        assert all(f == 0 for a,f in zip(d['age_cell_start'],flows) if a >= 46)
        value = sum(rates)
        assert math.isclose(value, d['period_tfr_topcode_adjusted'], abs_tol=1e-12)
        models[case][row['calendar_year']] = value
blocks = []
for decision in [2007,2011,2015,2019]:
    a = [r for r in annual if decision < r['year'] <= decision+4]
    assert len(a) == 4
    row = dict(decision_year=decision, empirical_window=f'{decision+1}–{decision+4}')
    for key in annual[0]:
        if key != 'year':
            row[key] = sum(r[key] for r in a)/4
    for case in models:
        row[case] = models[case][decision]
        row[case+'_gap_published'] = row[case]-row['published_tfr']
        row[case+'_gap_partial'] = row[case]-row['partial_18_44']
    blocks.append(row)
for name, rows in [('annual.csv',annual),('blocks.csv',blocks)]:
    with (OUT/name).open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def decline(rows,key):
    return 100*(1-rows[-1][key]/rows[0][key])
result = {'annual_decline_percent':{k:decline(annual,k) for k in ['published_tfr','partial_18_44','approximate_18_45']},
    'block_decline_percent':{k:decline(blocks,k) for k in ['published_tfr','partial_18_44','approximate_18_45',*models]},
    'checks':'17 annual published TFR sums; 6 overlapping source years; difference decomposition; 3 saved model flow/rate replays passed',
    'source_sha256':{name:hashlib.sha256((SOURCE/name).read_bytes()).hexdigest() for name in FILES}}
(OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
table = '\n'.join(f"| {r['empirical_window']} | {r['published_tfr']:.4f} | {r['partial_18_44']:.4f} | {r['approximate_18_45']:.4f} |" for r in blocks)
(OUT/'README.md').write_text(f'''# Fertility age-support check

The age-support sensitivity reduces the observed decline across the matched four-year windows from **{decline(blocks,'published_tfr'):.2f}% to {decline(blocks,'partial_18_44'):.2f}%**. This matters modestly; it does not eliminate the historical decline or settle the target convention. No calibration target, weight, model equation or production input was changed. No equilibrium was solved; the full continuation remains paused.

| Birth window | Published TFR | Ages 18–44 | Approximate ages 18–45 |
|---|---:|---:|---:|
{table}

Annual 2007-to-2023 declines are {decline(annual,'published_tfr'):.2f}% for published TFR and {decline(annual,'partial_18_44'):.2f}% for ages 18–44. These differ from the four-year-window comparisons above.

## Construction and limitations

The partial index is (2 × rate18–19 + 5 × sum of rates20–24 through40–44)/1000, in births per woman. Four-year observations average the annual indices equally, preserving the existing empirical block convention. The model decision in 2007 maps to births2008–2011; decision2019 maps to births2020–2023.

The model covers ages18–45. The cached tables do not isolate age45. The approximate column adds one year at the published oldest-group rate; this is an explicit approximation, not an exact age45 rate. NCHS computes its oldest rate using births to women45+ divided by women45–49. The full older-group contribution is saved for scale, not used as a mathematical bound on the age45 contribution.

The published15–19 rate need not equal a duration-weighted average of the15–17 and18–19 rates because female exposures differ. Annual CSV separates the young-age contribution, the teen regrouping difference, and the older contribution; their sum reproduces the published-minus-partial difference exactly. No female exposure counts were inferred from rounded rates.

The existing initial first-birth timing contract collapses boundary ages into model endpoint cells and preserves all first births. It is a different moment and was not changed. Treating one model household as one representative potential mother remains an approximation, not a newly established arithmetic error. Moving to a support-restricted production target would require a consistent initial-level and transition observation convention; this audit alone does not authorize that change.

## Saved model check

All three short-path selected mappings exactly match their saved final mappings. Recomputed age-specific rates from top-code-adjusted birth flows divided by model age mass reproduce the saved period fertility indices. No extra factor four is needed: annualization cancels the four-year age-cell width. No births occur outside the modeled fertility cells.

The three saved short paths imply declines of {decline(blocks,'delta_m0025'):.2f}%, {decline(blocks,'delta_m005'):.2f}%, and {decline(blocks,'delta_m010'):.2f}% across these windows. These are shock diagnostics under the earlier unrestricted calibration, not a fitted shock path or results for the new beta cap. Their complete prior calibration fit and parameter tables remain in ../return_home_20260911/READOUT.md. `blocks.csv` reports all four observations for all three paths and gaps against both empirical definitions. Long-horizon acceptance is not established by these short paths.

## Reproduce and sources

Run `python output/model/e5f_matched_pf_20260909a/current_candidate_transition/check_fertility_age_support.py` from the project root. This reads existing cached source text and saved model outputs only. See verification.json for source hashes and checked identities.

- NCHS, Births: Final Data for2015, Table4 (PDF page20): https://www.cdc.gov/nchs/data/nvsr/nvsr66/NVSR66_01.pdf —2007–2009 inputs.
- NCHS, Births: Final Data for2023, Table2 (PDF page14; oldest-age footnote page15): https://www.cdc.gov/nchs/data/nvsr/nvsr74/nvsr74-1.pdf —2010–2023 inputs.
''')
print(json.dumps(result,indent=2))
print(table)
