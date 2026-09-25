#!/usr/bin/env python3
"""Torch-only, report-only table layout fixture; no model import or solve."""
import json,csv
from pathlib import Path
import pymupdf as fitz
from render_pair import table_page,fmt
WORK=Path('/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/nightpair_20260925_v1')
obj=json.loads((WORK/'inputs/objective.json').read_text())
fits=[]
for r in obj['target_rows']:
    fits.append([r['restriction_id'],fmt(r['target']),fmt(r['target']),fmt(0),
                 fmt(r['actual_weight']),fmt(0) if r['actual_weight'] is not None else ''])
params=[]
for r in obj['parameter_restrictions']:
    params.append([r['parameter'],fmt((r['lower']+r['upper'])/2),fmt(r['lower']),fmt(r['upper']),
                   'False','experimental free coordinate'])
for name,status in [('theta1','externally fixed B15'),('psi_child','normalized to 2.1'),
    ('payroll_tax','experimental PAYGO rate'),('pension_period','endogenous balanced PAYGO'),
    ('housing_supply_elasticity','retained external setting'),('tenure_choice_kappa','retained external setting'),
    ('alpha_cons','retained external setting'),('sigma','retained external setting'),
    ('selling_cost','retained external setting'),('financed_share','retained external setting'),
    ('annual_depreciation','author adopted input'),('period_depreciation','four-year compounded'),
    ('annual_property_tax','author adopted input'),('period_property_tax','four-year linear source convention'),
    ('income_process','retained B15 persistent-state count'),
    ('entrant_conversion_factor','legacy child-departure diagnostic; inactive in split-birth entry'),
    ('adult_entry_birth_to_household_conversion','effective closed stationary birth conversion')]:
    params.append([name,fmt(.5),'','','',status])
doc=fitz.open()
table_page(doc,'Complete target fit','Fixture, no model results',
    ['Moment','Target','Model','Gap','Weight','Loss'],[164,105,105,105,100,137],fits,31)
for index,start in enumerate((0,13),1):
    table_page(doc,f'Parameters and restrictions ({index}/2)','Fixture, no model results',
        ['Parameter','Estimate','Lower','Upper','Near','Status'],[178,120,80,80,75,183],
        params[start:start+13],32)
assert len(doc)==3 and all(name[0] in doc[0].get_text() for name in fits)
assert all(p[0] in doc[1].get_text()+doc[2].get_text() for p in params)
out=WORK/'results/report_layout_fixture.pdf';out.parent.mkdir(exist_ok=True)
doc.save(out)
print(json.dumps(dict(status='report_only_fixture_passed',pages=3,target_rows=len(fits),
                      parameter_rows=len(params),pdf=str(out)),sort_keys=True))
