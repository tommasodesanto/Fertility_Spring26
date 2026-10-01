"""Synthetic logit, masking and strict-gate checks; zero model calls."""
import json
import numpy as np
from reduce_saved import recover, reduce

k = .2
shape = (2, 1, 1, 2, 1)
wait = np.full(shape, 3.)
attempt = wait + np.array([-.4, .3])[:, None, None, None, None]
def fixture(a0, a1):
    actions = np.stack([a0, a1], axis=-1)
    mx = actions.max(axis=-1)
    x = np.exp((actions-mx[..., None])/k)
    probs = np.zeros(shape+(4,))
    probs[..., :2] = x/x.sum(axis=-1, keepdims=True)
    F = mx+k*np.log(x.sum(axis=-1))
    return probs, F
p0, f0 = fixture(wait, attempt)
p1, f1 = fixture(wait+.1, attempt+np.array([.2, -.2])[:, None, None, None, None])
valid, action, gap, sums = recover(p0, f0, k)
assert valid.all() and np.allclose(action[..., 0], wait) and np.allclose(action[..., 1], attempt)
assert np.allclose(gap, attempt-wait)
w = np.ones(shape)
r = reduce(w, [p0,p1], [f0,f1], k, [18,22], [True,False])
assert r['recovered_mass'] == 2 and r['excluded_mass'] == 2
assert np.isclose(r['rows'][0]['positive_gap_fraction'], .5)
assert np.isclose(r['rows'][0]['negative_gap_fraction'], .5)
assert np.isclose(r['rows'][0]['delta_wait_value'], .1)
assert np.isclose(r['rows'][0]['delta_try_value'], 0)
assert np.isclose(r['rows'][0]['delta_gap'], -.1)
corner = p1.copy(); corner[0,0,0,0,0,:2] = [0,1]
dead = f1.copy(); dead[1,0,0,0,0] = -1e10
r2 = reduce(w,[p0,corner],[f0,dead],k,[18,22],[True,False])
assert r2['recovered_mass'] == 0 and r2['weighted_excluded_fraction'] == 1
bad = p0.copy(); bad[0,0,0,0,0,0] += .01
try:
    recover(bad,f0,k)
except ValueError as e:
    assert 'sum' in str(e)
else:
    raise AssertionError('Mutated probability accepted')
bad = p0.copy(); bad[0,0,0,0,0,0] = np.nan
try:
    recover(bad,f0,k)
except ValueError:
    pass
else:
    raise AssertionError('Nonfinite probability accepted')
print(json.dumps(dict(status='passed',model_calls=0,checks=['interior action-level and gap inversion','matched weights and positive/negative fractions','zero-fecundity age exclusion','corner/dead-value weighted exclusion','mutated sum rejection','nonfinite probability rejection']),indent=2))
