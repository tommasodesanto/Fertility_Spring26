"""Bounded source verification and uncertainty candidates; never adopts weights.

Reads only the two already indexed CPS partitions and the small NCHS cache.
No model imports, solves, raw-source edits, or external writes.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import statistics

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[4]
RECEIPT = ROOT / 'output/model/e5f_matched_pf_20260909a/parameter_target_audit/fertility/fertility_availability.json'
NCHS = ROOT / 'code/data/nchs_natality_timing'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def moments(rows):
    total = math.fsum(r['w'] for r in rows)
    moms = math.fsum(r['w'] for r in rows if r['n'] > 0)
    numerators = [math.fsum(r['w'] for r in rows if r['n'] == k) for k in range(3)]
    p3 = math.fsum(r['w'] for r in rows if r['n'] >= 3)
    return dict(n=len(rows), weight=total, mothers_weight=moms,
                childless=numerators[0]/total, exactly_one_among_mothers=numerators[1]/moms,
                parity1=numerators[1]/total, parity2=numerators[2]/total, parity3plus=p3/total,
                completed_capped5=math.fsum(r['w']*min(r['n'],5) for r in rows)/total,
                completed_uncapped=math.fsum(r['w']*r['n'] for r in rows)/total,
                top_parity_capped_mean=math.fsum(r['w']*min(r['n'],5) for r in rows if r['n']>=3)/p3)


def main():
    original = json.loads(RECEIPT.read_text())
    assert sha(original['schema_source']) == original['schema_sha256']
    raw = Path(original['cps_source'])
    assert raw.stat().st_size == original['raw_size_bytes']
    records, verification = {}, {}
    with raw.open('rb') as f:
        for year in [2004, 2006]:
            p = original['partitions'][str(year)]
            f.seek(p['byte_start'])
            data = f.read(p['byte_end_exclusive']-p['byte_start'])
            assert hashlib.sha256(data).hexdigest() == p['partition_sha256']
            assert len(data) % 261 == 0
            selected = []
            for offset in range(0, len(data), 261):
                row = data[offset:offset+261]
                assert int(row[:4]) == year and int(row[9:11]) == 6
                assert row.endswith(b'\n')
                if int(row[148:149]) != 2 or not 40 <= int(row[146:148]) <= 44:
                    continue
                n = int(row[238:241]); w = int(row[250:260])/10000
                assert n in range(21) or n == 999
                if n != 999 and w > 0:
                    selected.append(dict(age=int(row[146:148]), n=n, w=w))
            records[year] = selected
            verification[year] = dict(partition_sha256=p['partition_sha256'],
                                      records=len(selected), bytes_checked=len(data))
    cps = {str(y): moments(r) for y,r in records.items()}
    pooled = records[2004] + records[2006]
    cps['pooled'] = moments(pooled)
    target = original['cps_moments'][-1]
    for a,b in [('childless','childless_share'), ('exactly_one_among_mothers','exactly_one_given_mother'),
                ('completed_capped5','mean_children_ever_born_capped5')]:
        assert abs(cps['pooled'][a]-target[b]) < 2e-14
    cps_age = {str(age): moments([r for r in pooled if r['age']==age]) for age in range(40,45)}
    cps_age['40_41'] = moments([r for r in pooled if r['age']<=41])
    cps_age['42_44'] = moments([r for r in pooled if r['age']>=42])
    uncertainty = {}
    for key,denom in [('childless','weight'),('exactly_one_among_mothers','mothers_weight')]:
        annual = {}
        for y in ['2004','2006']:
            p = cps[y][key]; x = cps[y][denom]
            annual[y] = dict(value=p, denominator_population_units=x, b_parameter=2016,
                             approximate_gvf_se=math.sqrt(2016*p*(1-p)/x))
        alpha = cps['2004'][denom]/sum(cps[y][denom] for y in ['2004','2006'])
        a = alpha*annual['2004']['approximate_gvf_se']
        b = (1-alpha)*annual['2006']['approximate_gvf_se']
        uncertainty[key] = dict(annual=annual, base_share_2004=alpha,
            pooled_fixed_bases_zero_crossyear_covariance_se=math.hypot(a,b),
            pooled_fixed_bases_max_positive_correlation_se=a+b,
            qualification='Annual official generalized-variance approximation. Pooled candidates condition on observed base shares; cross-year covariance and random conditional-mother bases are not certified. Neither scale is adopted.')
    assert sha(NCHS/'first_birth_counts_year_age.csv') == original['nchs_sha256']
    counts = [{k:int(v) for k,v in r.items()} for r in csv.DictReader((NCHS/'first_birth_counts_year_age.csv').open())]
    rows = [r for r in counts if 2003<=r['year']<=2006]
    def midpoint(a): return 20 if a<22 else 44 if a>=42 else 20+4*((a-18)//4)
    def timing(rs):
        n = sum(r['n_first_births'] for r in rs)
        mu = sum(r['n_first_births']*midpoint(r['age']) for r in rs)/n
        p = sum(r['n_first_births'] for r in rs if r['age']>=30)/n
        vx = sum(r['n_first_births']*(midpoint(r['age'])-mu)**2 for r in rs)/n
        cov = sum(r['n_first_births']*(midpoint(r['age'])-mu)*((r['age']>=30)-p) for r in rs)/n
        return dict(n=n, midpoint_mean=mu, share30plus=p,
            conditional_multinomial_process_se_mean=math.sqrt(vx/n),
            conditional_multinomial_process_se_share=math.sqrt(p*(1-p)/n),
            conditional_multinomial_process_covariance=cov/n,
            first_births_below18=sum(r['n_first_births'] for r in rs if r['age']<18),
            first_births_above45=sum(r['n_first_births'] for r in rs if r['age']>45))
    nchs = dict(pooled=timing(rows), annual={str(y):timing([r for r in rows if r['year']==y]) for y in range(2003,2007)})
    assert nchs['pooled']['n'] == 6611269
    assert abs(nchs['pooled']['midpoint_mean']-25.976263860992496)<1e-13
    assert abs(nchs['pooled']['share30plus']-.2492780130410667)<1e-14
    nchs['leave_one_year_out'] = {str(y):timing([r for r in rows if r['year']!=y]) for y in range(2003,2007)}
    nchs['temporal_variation_not_sampling_se'] = {
        key: dict(annual_sample_sd=statistics.stdev(r[key] for r in nchs['annual'].values()),
                  annual_range=max(r[key] for r in nchs['annual'].values())-min(r[key] for r in nchs['annual'].values()),
                  maximum_leave_one_year_shift=max(abs(r[key]-nchs['pooled'][key]) for r in nchs['leave_one_year_out'].values()))
        for key in ['midpoint_mean','share30plus']}
    exclusions = [r for r in csv.DictReader((NCHS/'order_exclusion_shares_by_year.csv').open()) if 2003<=int(r['year'])<=2006]
    nchs['order_exclusions'] = exclusions
    u = sum(int(r['n_unknown_live_birth_order']) for r in exclusions)
    n = nchs['pooled']['n']; mu=nchs['pooled']['midpoint_mean']; p=nchs['pooled']['share30plus']
    nchs['extreme_unknown_order_bounds_not_se'] = dict(unknown_order_records=u,
        mean_lower=(n*mu+20*u)/(n+u), mean_upper=(n*mu+44*u)/(n+u),
        share_lower=n*p/(n+u), share_upper=(n*p+u)/(n+u),
        qualification='Conservative: every unknown-order record could be first order and at either midpoint extreme; no assumption that missingness is random.')
    cp=cps['pooled']; top=3.602359422009
    p3=(2.1-cp['parity1']-2*(1-cp['childless']-cp['parity1']))/(top-2)
    algebra=dict(normalization=2.1, inherited_top_bin_weight=top,
        required_p3plus_if_same_age_parity_distribution=p3,
        observed_p3plus=cp['parity3plus'], implied_p2=1-cp['childless']-cp['parity1']-p3,
        qualification='Algebraic conditional diagnostic only: actual 2.1 normalization and CPS age-window observers need not be the same distribution.')
    paths=[RECEIPT,NCHS/'build_first_birth_timing.R',NCHS/'first_birth_counts_year_age.csv',
           NCHS/'first_birth_counts_manifest.csv',NCHS/'order_exclusion_shares_by_year.csv',
           Path(original['schema_source']), Path(__file__)]
    payload=dict(status='source_checks_passed_weights_not_adopted', source_sha256={str(p):sha(p) for p in paths},
        cps_partition_checks=verification,cps=cps,cps_age_detail=cps_age,cps_uncertainty_candidates=uncertainty,
        nchs=nchs,normalization_algebra=algebra)
    (OUT/'verification.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(dict(status=payload['status'],cps_uncertainty=uncertainty,nchs_pooled=nchs['pooled'],normalization_algebra=algebra),indent=2))


if __name__ == '__main__':
    main()
