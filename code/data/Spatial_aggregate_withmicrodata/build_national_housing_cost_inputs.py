"""Reproduce the national 2007 depreciation candidate from retained official data.

No calibration inputs or model parameters are modified. Run with Python/openpyxl.
"""
import csv
import hashlib
import json
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/national_housing_inputs/depreciation'
SRC = OUT / 'sources'


def series(name):
    with (SRC / (name + '.csv')).open() as f:
        return {int(r['observation_date'][:4]): float(r[name])
                for r in csv.DictReader(f) if r[name] not in ('', '.')}


def main():
    book = openpyxl.load_workbook(SRC / 'BEA_Section5.xlsx', data_only=True)
    def bea(sheet, code, year):
        rows = list(book[sheet].values)
        col = list(rows[7]).index(str(year))
        return float(next(r[col] for r in rows[8:] if r[2] == code))
    houses = series('BOGZ1FL155035013A')
    mobile = series('BOGZ1FL155012013A')
    fed_structures = series('BOGZ1FL155012665A')
    rows = []
    for year in (2006, 2007, 2008):
        stock = bea('FAAt501-A', 'k1r53101esoo', year)
        flow = bea('FAAt504-A', 'm1r53101esoo', year)
        household_stock = bea('FAAt501-A', 'k1r53105es00', year)
        household_flow = bea('FAAt504-A', 'm1r53105es00', year)
        value = houses[year] + mobile[year]
        land = 1 - stock / value
        rate = (flow / stock) * (1 - land)
        assert stock == fed_structures[year], (year, stock, fed_structures[year])
        assert 0 < land < 1 and 0 < rate < 1
        assert abs(rate - flow/value) < 1e-16
        rows.append(dict(year=year, owner_structures_stock_million=stock,
                         owner_depreciation_flow_million=flow,
                         owner_real_estate_ex_mobile_million=houses[year],
                         mobile_home_value_million=mobile[year],
                         housing_value_million=value, implied_land_share=land,
                         structures_rate=flow/stock, annual_depreciation=rate,
                         household_structures_rate=household_flow/household_stock,
                         household_rate_times_owner_structure_share=(household_flow/household_stock)*(1-land)))
    with (OUT / 'estimates.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    receipt = dict(status='calculated_and_adopted_for_next_calibration', reference_year=2007,
        units='Current dollars in millions; depreciation annual flow, stocks year-end',
        estimator='BEA owner-occupied depreciation / BEA owner-occupied structures stock times (1 - implied national land share); algebraically depreciation / national housing value',
        housing_value='Fed FL155035013 + FL155012013: owner-occupied real estate plus mobile homes, excluding separately held vacant land',
        sources=[dict(file='BEA_Section5.xlsx', url='https://apps.bea.gov/national/FixedAssets/Release/XLS/Section5All_xls.xlsx', vintage=book['FAAt501-A']['A6'].value),
                 *[dict(file=n+'.csv',url='https://fred.stlouisfed.org/graph/fredgraph.csv?id='+n) for n in ('BOGZ1FL155035013A','BOGZ1FL155012013A','BOGZ1FL155012665A')]],
        retrieved='2026-09-23', candidate=rows[1],
        checks=['Owner structures match independent Federal Reserve series exactly in all three years', 'Product formula equals direct flow/value within 1e-16', 'Positive rate and land share strictly between zero and one'],
        limitations=[
            'National owner-occupied stock proxy; applying one common rate to rental housing remains a model approximation.',
            'Residual land share combines replacement-cost structures with market-value housing; it is not a direct land appraisal.',
            'Unlike DUE, uses matched owner-occupied BEA row 11, and a 2007 national Fed residual land share rather than a later Bay Area tract estimate. Row 7 household-rate variant retained as sensitivity.',
            'BEA stocks include residential fixed assets; this calculation retains the Fed mapping to BEA owner-occupied row 11 rather than claiming parcel-by-parcel identity.',
            'Annual flow divided by year-end stock follows the selected convention. Adjacent years show sensitivity, not an estimated confidence interval.',
            'Official series are revised. Retained download hashes govern these values; cached web pages can contain earlier vintages.',
            'No model, target contract, objective, or paper edits; no new solves.'])
    for item in receipt['sources']:
        item['sha256']=hashlib.sha256((SRC/item['file']).read_bytes()).hexdigest()
    (OUT/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    tax_source = OUT.parent / 'property_tax/sources/us_estimates_sequence_0103.csv'
    with tax_source.open() as f:
        tax_row = next(csv.reader(f))
    assert tax_row[:6] == ['ACSSF', '2011e5', 'us', '000', '0103', '0000001']
    taxes, value = int(tax_row[120]), int(tax_row[55])
    assert sum(map(int, tax_row[121:123])) == taxes
    assert sum(map(int, tax_row[56:58])) == value
    results = [dict(parameter='housing_depreciation', annual_rate=rows[1]['annual_depreciation'],
                    reference_period='2007', status='adopted_for_next_calibration'),
               dict(parameter='property_tax', annual_rate=taxes/value,
                    reference_period='2007–2011 ACS five-year', status='adopted_for_next_calibration')]
    with (OUT.parent / 'national_inputs.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0])); writer.writeheader(); writer.writerows(results)
    print(json.dumps(results,indent=2))


if __name__ == '__main__':
    main()
