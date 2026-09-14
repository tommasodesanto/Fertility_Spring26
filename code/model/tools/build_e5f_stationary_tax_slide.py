"""Build a one-slide steady-state tax comparison from verified receipts."""
from pathlib import Path
import csv
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parents[3]
PACKET = ROOT / 'output/model/e5f_original_queue_20260913a'
OUT = PACKET / 'inherited_2023_tax/steady_state_readout'


def read(p):
    return json.loads(p.read_text())


def main():
    paths = [PACKET/'terminal_restart_v1/verified_endpoint_receipt.json',
             PACKET/'inherited_2023_tax/tax2_endpoint_receipt.json']
    receipts = [read(p) for p in paths]
    fresh_path = PACKET/'inherited_2023_tax/tax1_endpoint_verification.json'
    fresh = read(fresh_path)
    assert fresh['passed'] and fresh['fresh_audit']['status'] == 'passed'
    assert all(fresh['fresh_audit']['checks'].values())
    for r in receipts:
        assert r['verified'] and r['stationary_endpoint_verified'] and r['converged']
        assert r['one_step_audit']['status'] == 'passed' and all(r['one_step_audit']['checks'].values())
        assert r['fresh_final_mapping_matches_endpoint']
        q = r['endpoint_reference']
        assert abs(q['renewal_ratio']-1) <= 5e-8
        assert math.isclose(q['housing_demand'], q['housing_supply'], rel_tol=5e-8)
        assert q['payroll_tax_rate'] == .179
    base, policy = [r['endpoint_reference'] for r in receipts]
    def values(q):
        return [('Households', q['population_households'], 'index'),
                ('Total housing services', q['housing_demand'], 'index'),
                ('Housing services per household', q['housing_demand']/q['population_households'], 'index'),
                ('House price', q['asset_price'], 'index'),
                ('Rent per unit of housing', q['renter_price'], 'index'),
                ('Homeownership (percent)', 100*q['owner_rate'], 'percent')]
    rows = []
    for (label,b,unit),(_,p,_) in zip(values(base),values(policy)):
        rows.append(dict(outcome=label,units=unit,baseline_raw=b,policy_raw=p,
            baseline_display=100 if unit=='index' else b,
            policy_display=100*p/b if unit=='index' else p,
            effect=100*(p/b-1) if unit=='index' else p-b,
            effect_units='percent' if unit=='index' else 'percentage points'))
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT/'comparison.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    lines='\n'.join(f"{r['outcome']} & {r['baseline_display']:.1f} & {r['policy_display']:.1f} \\\\" for r in rows)
    tex=r'''\documentclass[11pt,aspectratio=169]{beamer}
\usepackage[T1]{fontenc}
\usepackage{lmodern,booktabs,tabularx,array}
\setbeamertemplate{navigation symbols}{}
\setbeamertemplate{footline}[frame number]
\linespread{1.15}
\begin{document}
\begin{frame}{Property Tax: Steady-State Comparison}
All property-tax revenue is rebated equally to households.
\par\medskip
\renewcommand{\arraystretch}{1.15}
\begin{tabularx}{\textwidth}{@{}Xrr@{}}
\toprule
 & Annual tax: 1\% & Annual tax: 2\% \\
\midrule
__ROWS__
\bottomrule
\end{tabularx}
\par\medskip
More households, but less housing per household.
\vfill
{\footnotesize Steady states at the same post-decline fertility preferences; balanced pensions.\\
Quantity and price indices: 1\% tax = 100. Homeownership is in percent.}
\end{frame}
\end{document}
'''.replace('__ROWS__',lines)
    tex_path=ROOT/'latex/appendix_property_tax_steady_states.tex'
    tex_path.write_text(tex)
    inputs=paths+[fresh_path]
    receipt=dict(status='PASS',model_solves=0,interpretation='Comparison of two verified stationary equilibria; not the impact from inherited2023 or a converged policy transition.',
        annual_taxes=[.01,.02],equal_rebates=True,fixed_payroll_tax=.179,
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        effects=rows,tex_sha256=hashlib.sha256(tex_path.read_bytes()).hexdigest())
    (OUT/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({r['outcome']:round(r['effect'],4) for r in rows}))


if __name__=='__main__':
    main()
