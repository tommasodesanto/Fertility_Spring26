from pathlib import Path
import sys,pickle,json,csv,hashlib
import numpy as np
r=Path.cwd(); f=r/'tmp/e5f_overnight_20260905a/frozen_production/code/model';sys.path[:0]=[str(f/'tools'),str(f)]
with (r/'output/model/e5f_overnight_independent_verification_20260905a/numerical_smoke/dated_state.pkl').open('rb') as h:d=pickle.load(h)
p=d['parameters'];e=d['evaluation'];arr=[]
for field,value in vars(e.policy).items():
 if isinstance(value,np.ndarray):arr.append({'field':field,'shape':list(value.shape),'dtype':str(value.dtype),'bytes':int(value.nbytes)})
invalid=np.zeros((4,4),bool)
for n in range(4):invalid[n,n+1:]=True
badmass=float(e.g_pre[...,invalid].sum());assert badmass==0
checks={'saved_array_footprint':{'policy_arrays':arr,'policy_array_bytes':sum(x['bytes'] for x in arr),'continuation_probability_bytes':int(p._fert2_probs.nbytes),'valid_family_states':10,'allocated_family_states':16,'invalid_pre_fertility_mass':badmass,'grid_size':len(d['b_grid']),'location_count':p.I,'income_states':p.Nz,'uses_compiled_forward_distribution_flag':bool(getattr(p,'use_compiled_forward_distribution',False)),'uses_numba_scatter_flag':bool(getattr(p,'use_numba_scatter',False))}}
policyroot=r/'output/model/e5f_post2023_policy_mechanisms_eta063_kappa005_task010_production_20260904a';exposure={}
for ff in policyroot.glob('*/policy_path.csv'):
 rows=list(csv.DictReader(ff.open()));exposure[ff.parent.name]={'dates':len(rows),'fallback_dates':[x['calendar_year'] for x in rows if x['grid_resolution_fallback'].lower() in ['true','1']],'source_sha256':hashlib.sha256(ff.read_bytes()).hexdigest()}
checks['saved_policy_retry_exposure']=exposure
from intergen_eqscale_seq_optimized.calibration import extract_moments
from types import SimpleNamespace
par=np.array([0.,1.,0.,0.]);s=SimpleNamespace(mean_parity=1.,mean_completed_fertility=1.,parity_dist=par)
checks['verbose_fertility_units_probe']={'printed_legacy_TFR':2*s.mean_parity,'correct_literal_completed_fertility':extract_moments(s,SimpleNamespace(fertility_units='literal_topcode',tfr_top_bin_weight=3.602359422009))['tfr'],'example':'Every household has exactly one child'}
out={'python':sys.version,'numpy':np.__version__,'source':'frozen production','checks':checks}
(r/'output/model/e5f_code_architecture_audit_20260905a/saved_state_footprint.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
