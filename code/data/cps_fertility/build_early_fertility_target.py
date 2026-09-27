"""CPS June 2004/2006 children ever born at exact age 25; Torch-only builder.

Reuses the certified fixed-width schema and June-partition hashes from the
September 10 fertility availability receipt. No target or weight is adopted.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import gzip
import json
import math
import numpy as np


def estimate(rows, seed, draws):
    x = np.asarray([r[2] for r in rows], dtype=float)
    w = np.asarray([r[3] for r in rows], dtype=float)
    years = np.asarray([r[0] for r in rows], dtype=int)
    sw = math.fsum(w)
    mean = math.fsum(x * w) / sw
    rng = np.random.default_rng(seed)
    bs = np.zeros(draws)
    bs3 = np.zeros(draws)
    strata = [np.flatnonzero(years == y) for y in (2004, 2006)]
    for b in range(draws):
        take = np.concatenate([rng.choice(s, len(s), replace=True) for s in strata if len(s)])
        bs[b] = np.dot(x[take], w[take]) / w[take].sum()
        bs3[b] = np.dot(np.minimum(x[take], 3), w[take]) / w[take].sum()
    return dict(n=len(rows), sum_supplement_weights=sw,
                mean_children_ever_born=mean,
                mean_children_ever_born_capped3=math.fsum(np.minimum(x, 3) * w) / sw,
                mean_children_ever_born_capped5=math.fsum(np.minimum(x, 5) * w) / sw,
                share_with_any_birth=math.fsum(w[x > 0]) / sw,
                childless_share=math.fsum(w[x == 0]) / sw,
                n_frever_above3=int((x > 3).sum()),
                weighted_share_frever_above3=math.fsum(w[x > 3]) / sw,
                n_frever_above5=int((x > 5).sum()),
                independent_person_year_stratified_bootstrap_se=float(bs.std(ddof=1)),
                capped3_bootstrap_se=float(bs3.std(ddof=1)),
                capped3_bootstrap_percentile_interval95=[float(v) for v in np.quantile(bs3, [.025, .975])],
                bootstrap_percentile_interval95=[float(v) for v in np.quantile(bs, [.025, .975])],
                annual={str(y): dict(n=int((years == y).sum()),
                    mean_children_ever_born=float(np.dot(x[years == y], w[years == y]) / w[years == y].sum()))
                    for y in (2004, 2006) if np.any(years == y)})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--receipt', type=Path, required=True)
    p.add_argument('--partitions', type=Path, required=True)
    p.add_argument('--compressed-source', type=Path)
    p.add_argument('--schema', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--seed', type=int, default=20260926)
    p.add_argument('--bootstrap-draws', type=int, default=2000)
    args = p.parse_args()
    assert args.bootstrap_draws >= 1000
    original = json.loads(args.receipt.read_text())
    assert hashlib.sha256(args.schema.read_bytes()).hexdigest() == original['schema_sha256']
    selected, completed, checks = [], [], {}
    raw_stream = gzip.open(args.compressed_source, 'rb') if args.compressed_source else None
    for year in (2004, 2006):
        info = original['partitions'][str(year)]
        source = args.partitions / f'june_{year}.dat'
        if raw_stream is None:
            data = source.read_bytes()
        else:
            raw_stream.seek(info['byte_start'])
            data = raw_stream.read(info['byte_end_exclusive'] - info['byte_start'])
            source.write_bytes(data)
        assert len(data) == info['byte_end_exclusive'] - info['byte_start']
        digest = hashlib.sha256(data).hexdigest()
        assert digest == info['partition_sha256']
        counts = {'all_records': 0, 'female24_26': 0, 'female24_26_invalid_frever': 0,
                  'female24_26_nonpositive_weight': 0}
        assert len(data) % 261 == 0
        for offset in range(0, len(data), 261):
            row = data[offset:offset + 261]
            assert int(row[:4]) == year and int(row[9:11]) == 6 and row.endswith(b'\n')
            counts['all_records'] += 1
            if int(row[148:149]) != 2:
                continue
            age, n, w = int(row[146:148]), int(row[238:241]), int(row[250:260]) / 10000
            if 24 <= age <= 26 or 40 <= age <= 44:
                assert n in range(21) or n == 999
            if 24 <= age <= 26:
                counts['female24_26'] += 1
                if n == 999:
                    counts['female24_26_invalid_frever'] += 1
                elif w <= 0:
                    counts['female24_26_nonpositive_weight'] += 1
                else:
                    selected.append((year, age, n, w))
            if 40 <= age <= 44 and n != 999 and w > 0:
                completed.append((n, w))
        checks[str(year)] = dict(path=str(source), bytes=len(data), sha256=digest, counts=counts)
    if raw_stream is not None:
        raw_stream.close()
    denom = math.fsum(w for n, w in completed)
    control = dict(n=len(completed), childless_share=math.fsum(w for n, w in completed if n == 0) / denom,
        exactly_one_given_mother=math.fsum(w for n, w in completed if n == 1) / math.fsum(w for n, w in completed if n > 0))
    reference = original['cps_moments'][-1]
    assert control['n'] == reference['n']
    for key in ('childless_share', 'exactly_one_given_mother'):
        assert abs(control[key] - reference[key]) < 2e-14
    exact = [r for r in selected if r[1] == 25]
    exact_weight = math.fsum(r[3] for r in exact)
    print(json.dumps(dict(stage='point_estimate_and_control_verified', n_age25=len(exact),
        mean_age25=math.fsum(r[2] * r[3] for r in exact) / exact_weight,
        capped3_mean_age25=math.fsum(min(r[2], 3) * r[3] for r in exact) / exact_weight,
        completed40_44_control=control)), flush=True)
    estimates = {str(age): estimate([r for r in selected if r[1] == age], args.seed + age, args.bootstrap_draws)
                 for age in (24, 25, 26)}
    estimates['24_26_diagnostic'] = estimate(selected, args.seed, args.bootstrap_draws)
    result = dict(status='verified_empirical_candidate_not_adopted',
        authoritative_candidate='25',
        estimator='sum(FRSUPPWT * FREVER) / sum(FRSUPPWT), pooled records; not average annual means',
        source='IPUMS CPS extract3; June 2004 and June 2006 partitions of cps_00003.dat',
        sample='Women SEX=2 at exact age25; FREVER0-20, positiveFRSUPPWT. Ages24/26 and24-26 are diagnostics.',
        timing='Completed integer age at June interview; cumulative live births as of interview, not exactly at25thbirthday.',
        weights='FRSUPPWT / 10000 per supplied loader; ratio invariant to common scaling',
        primary_topcoding='Model-matched candidate E[min(FREVER,3)] because model topstate is3+; raw and capped5means retained separately.',
        model_mapping='Exact observed age25 means ageinterval[25,26); lead specifies uniform birth-time observer interpolation weight0.875. Builder only estimates empiricalobject.',
        uncertainty='Person bootstrap stratified by survey year, retaining original weights and fixed annual samplecounts; NOT CPSdesign-consistent: no PSU/stratum/replicateweight correction or crossyear covariance. No calibrationweight adopted.',
        bootstrap=dict(seed=args.seed, draws=args.bootstrap_draws, numpy_version=np.__version__),
        partition_checks=checks, completed40_44_control=control, estimates=estimates,
        schema_sha256=original['schema_sha256'],
        builder_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'early_fertility_target.json').write_text(json.dumps(result, indent=2) + '\n')
    table = [dict(age=age, n=e['n'], mean_children_ever_born=e['mean_children_ever_born'],
                  mean_children_ever_born_capped3=e['mean_children_ever_born_capped3'],
                  capped3_bootstrap_se=e['capped3_bootstrap_se'],
                  mean_children_ever_born_capped5=e['mean_children_ever_born_capped5'],
                  bootstrap_se=e['independent_person_year_stratified_bootstrap_se'],
                  sum_supplement_weights=e['sum_supplement_weights']) for age, e in estimates.items()]
    with (args.output / 'early_fertility_target.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(table[0]))
        writer.writeheader(); writer.writerows(table)
    print(json.dumps(dict(estimates=table, control=control, status=result['status']), indent=2))


if __name__ == '__main__':
    main()
