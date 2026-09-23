"""Candidate aggregation of saved coefficients; no regression or target edit.

Holding cohort weights fixed across endpoints makes a within-cohort contrast
invariant to a cohort-specific additive normalization. It does not establish
causal identification, individual-panel balance, or a valid standard error.
"""
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[7]
OUT = Path(__file__).resolve().parent
SOURCE = ROOT / 'output/model/e5f_first_birth_measurement_review_20260905a/reference_cluster/full'
AUDIT = ROOT / 'code/data/psid_followup_mar2026/output/first_birth_correction_review/verification.json'


def read(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def main():
    names = ['coefficients_original.csv', 'input_support.csv', 'estimation_support_original.csv']
    expected = json.loads(AUDIT.read_text())['input_hashes']
    hashes = {}
    for name in names:
        path = SOURCE / name
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert digest == expected[str(path.relative_to(ROOT))], name
        hashes[name] = digest
    coeff = {(int(r['cohort']), r['event']): (float(r['coefficient']), float(r['variance']))
             for r in read(SOURCE / names[0])}
    def weights(filename):
        out = defaultdict(float)
        for row in read(SOURCE / filename):
            if int(row['never_treated']):
                continue
            g, k = int(float(row['first_child_year'])), int(float(row['K']))
            out[g, k] += float(row['weight_sum'])
        return out
    inp, fitted = weights(names[1]), weights(names[2])
    before = {g for (g, k), v in fitted.items() if k == -1 and v > 0}
    after = {g for (g, k), v in fitted.items() if k == 3 and v > 0}
    common = sorted(before & after)
    assert common
    rows = []
    endpoint_normalizations = []
    for g in common:
        b0, v0 = coeff[g, 'F1event']
        b1, v1 = coeff[g, 'L3event']
        assert all(math.isfinite(v) for v in (b0, v0, b1, v1))
        for event_name, value, variance in [('F1event', b0, v0), ('L3event', b1, v1)]:
            assert variance >= 0
            if variance == 0:
                # A supported endpoint can itself be the sole omitted baseline.
                # Check this explicitly rather than dropping that cohort.
                supported_events = {('F7event' if k <= -7 else 'L11event' if k >= 11
                                     else f'F{-k}event' if k < 0 else f'L{k}event')
                                    for (gg, k), w in fitted.items() if gg == g and w > 0}
                assert 'F2event' not in supported_events and value == 0
                zero_events = {e for e in supported_events if coeff[g, e] == (0.0, 0.0)}
                assert zero_events == {event_name}, (g, zero_events)
                endpoint_normalizations.append(dict(cohort=g, sole_supported_baseline=event_name))
        rows.append(dict(cohort=g, before=b0, after=b1, contrast=b1-b0,
                         fitted_weight_before=fitted[g, -1], fitted_weight_after=fitted[g, 3]))
    current = sum(w * coeff[g, 'L3event'][0] for (g, k), w in inp.items() if k == 3) / sum(w for (g,k),w in inp.items() if k == 3)
    current -= sum(w * coeff[g, 'F1event'][0] for (g,k),w in inp.items() if k == -1) / sum(w for (g,k),w in inp.items() if k == -1)
    assert abs(current - .7202462623815278) < 1e-10
    variants = []
    for endpoint, label in [('before', 'common_cohorts_fixed_prebirth_IW'), ('after', 'common_cohorts_fixed_postbirth_IW')]:
        denom = sum(r['fitted_weight_' + endpoint] for r in rows)
        estimate = sum(r['fitted_weight_' + endpoint] * r['contrast'] for r in rows) / denom
        # Deliberately large distinct arbitrary constants; compare every cohort.
        gaps = []
        shifted_total = 0.0
        for r in rows:
            shift = (r['cohort'] % 7 - 3) * .7
            delta = (r['after'] + shift) - (r['before'] + shift)
            gaps.append(abs(delta - r['contrast']))
            shifted_total += r['fitted_weight_' + endpoint] * delta / denom
        assert max(gaps) < 1e-12 and abs(shifted_total-estimate) < 1e-12
        variants.append(dict(name=label, estimate=estimate, standard_error=None,
                             shifted_estimate=shifted_total, maximum_cohort_invariance_gap=max(gaps)))
    receipt = dict(status='candidate saved-coefficient aggregation only; not adopted',
                   source_hashes=hashes, current_target_reproduced=current,
                   cohorts_before=len(before), cohorts_after=len(after), common_cohorts=len(common),
                   excluded_before=sorted(before-after), excluded_after=sorted(after-before),
                   before_weight_retained=sum(fitted[g,-1] for g in common)/sum(fitted[g,-1] for g in before),
                   after_weight_retained=sum(fitted[g,3] for g in common)/sum(fitted[g,3] for g in after),
                   variants=variants, no_new_regression=True,
                   endpoint_normalizations=endpoint_normalizations,
                   implementation_note='Initial positive-variance guard stopped at cohort 1971: its observed -1 is the sole omitted baseline. The corrected guard verifies this supported normalization explicitly instead of excluding the cohort.',
                   limitations=['Fixed cohort weights remove this additive-normalization dependence only.',
                                'Cohort composition and endpoint weighting differ from the frozen aggregate contrast.',
                                'People observed at both endpoints are not necessarily identical.',
                                'No full covariance matrix is saved, so no SE or calibration weight is supplied.',
                                'The author has not chosen this -1 to +3 contrast instead of a -2 reference window.',
                                'Parallel trends, controls, anticipation and the model counterfactual remain separate issues.'])
    (OUT/'common_cohort_contrast.json').write_text(json.dumps(receipt,indent=2)+'\n')
    with (OUT/'common_cohort_contrast.csv').open('w',newline='') as stream:
        w=csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n');w.writeheader();w.writerows(rows)
    print(json.dumps(receipt,indent=2))


if __name__ == '__main__':
    main()
