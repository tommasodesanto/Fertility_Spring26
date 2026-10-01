"""Saved first-birth probability metadata only; no engine imports or solves."""
import hashlib
import json
from pathlib import Path
import resource
import time
import numpy as np

root = Path('/scratch/td2248/projects/utility_floor_mechanism_responses_v1_replacement2/results/responses')
started = time.monotonic()
with np.load(root/'q0_reference_inherited_states.npz', allow_pickle=False) as z:
    pre = z['g_pre']
assert hashlib.sha256(pre.tobytes()).hexdigest() == 'caf760cb04c6a941cb6890510a832d5c0b61b2744f4fa4c16b89e1c8f3317903'
w = pre[..., 0, 0].copy()
del pre
result = {'model_calls':0,'array_downloads':0,'cases':{},'baseline_pre_childless_mass':float(w.sum())}
for case in ('00_reference_p1.00','04_lifetime_repayment_only_p1.00'):
    with np.load(root/case/'solution_arrays.npz', allow_pickle=False) as z:
        p = z['fert_probs']
    assert p.shape[:-1] == w.shape
    s = p[..., :2].sum(axis=-1)
    total = p.sum(axis=-1)
    nonzero = s > 0
    error = np.abs(s-1)
    bad = nonzero & (error > 1e-12)
    occupied = w > 0
    rows = []
    for j in range(w.shape[3]):
        sl = (slice(None),slice(None),slice(None),j,slice(None))
        rows.append({'age':18+4*j,'childless_pre_mass':float(w[sl].sum()),
                     'bad_count':int(bad[sl].sum()),'occupied_bad_count':int((bad[sl]&occupied[sl]).sum()),
                     'bad_childless_pre_mass':float(w[sl][bad[sl]].sum()),
                     'zero_sum_childless_pre_mass':float(w[sl][~nonzero[sl]].sum()),
                     'endpoint_childless_pre_mass':float(w[sl][(p[sl][...,0]==0)|(p[sl][...,1]==0)|(p[sl][...,0]==1)|(p[sl][...,1]==1)].sum())})
    ix = np.argwhere(bad)
    examples=[]
    for mask,label in ((bad,'any'),(bad&occupied,'occupied')):
        if mask.any():
            index=np.unravel_index(np.argmax(np.where(mask,error,-1)),error.shape)
            examples.append({'kind':label,'index':[int(x) for x in index],'probabilities':p[index].tolist(),
                             'first_two_sum':float(s[index]),'absolute_error':float(error[index]),'pre_childless_mass':float(w[index])})
    result['cases'][case]={'dtype':str(p.dtype),'shape':list(p.shape),'finite':bool(np.isfinite(p).all()),
        'minimum_probability':float(p.min()),'maximum_probability':float(p.max()),
        'unused_actions_2_3_max':float(np.abs(p[...,2:]).max()),
        'first_two_vs_all_actions_max_difference':float(np.abs(total-s).max()),
        'nonzero_sum_count':int(nonzero.sum()),'zero_sum_count':int((~nonzero).sum()),
        'nonzero_sum_min':float(s[nonzero].min()),'nonzero_sum_max':float(s[nonzero].max()),
        'max_nonzero_sum_abs_error':float(error[nonzero].max()),
        'nonzero_error_quantiles':dict(zip(('0','50','90','99','99.9','100'),[float(x) for x in np.quantile(error[nonzero],[0,.5,.9,.99,.999,1])])),
        'errors_above_1e12_count':int(bad.sum()),'occupied_errors_above_1e12_count':int((bad&occupied).sum()),
        'occupied_bad_pre_childless_mass':float(w[bad].sum()),
        'maximum_occupied_nonzero_sum_error':float(error[occupied&nonzero].max(initial=0)),
        'bad_zero_or_one_action_count':int((bad&((p[...,0]==0)|(p[...,1]==0)|(p[...,0]==1)|(p[...,1]==1))).sum()),
        'by_age':rows,'representative_indices':examples}
result.update(status='completed_metadata_only',elapsed_seconds=time.monotonic()-started,
              maximum_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              script_sha256=EXECUTED_SOURCE_SHA256)
print(json.dumps(result,indent=2,allow_nan=False))
