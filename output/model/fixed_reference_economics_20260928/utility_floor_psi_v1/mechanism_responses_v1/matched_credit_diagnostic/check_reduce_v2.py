"""Focused occupied-support checks plus original synthetic tests; zero solves."""
import json
import numpy as np
import check_reduce as original
from reduce_saved_v2 import reduce

w = original.w.copy()
p = original.p0.copy()
p[0,0,0,0,0,0] += .01
try:
    reduce(w,[p,original.p1],[original.f0,original.f1],original.k,[18,22],[True,False])
except ValueError as e:
    assert 'occupied matched' in str(e)
else:
    raise AssertionError('Occupied sum violation accepted')
w[0,0,0,0,0] = 0
r = reduce(w,[p,original.p1],[original.f0,original.f1],original.k,[18,22],[True,False])
d = r['probability_sum_diagnostics'][0]
assert d['global_bad_sum_count'] == 1 and d['occupied_bad_sum_count'] == 0
assert d['bad_sum_baseline_pre_mass'] == 0
assert d['global_maximum_nonzero_sum_error'] > .009
assert r['recovered_mass'] == 1 and r['excluded_mass'] == 2
assert r['probability_sum_tolerance'] == 1e-12
print(json.dumps({'status':'passed','model_calls':0,'checks':['original strict logit/weight/exclusion tests','positive-weight sum violation remains fatal','zero-weight violation excluded and globally reported','unchanged 1e-12 precision and recovery mask']},indent=2))
