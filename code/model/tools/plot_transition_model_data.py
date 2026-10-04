"""Reuse the standard assessment plots for the retained, dated 2023 transition."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name] = '1'
import csv, gzip, json, pickle, sys, importlib.util, subprocess
from pathlib import Path
from types import SimpleNamespace
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT/'code/model'), str(ROOT/'code/model/tools')]
from model_data_assessment import prepare_assessment, sha256
REF = ROOT/'output/model/transition_readiness_v1/current_baseline_20261003/retained_one_shock_v1'
OUT = REF/'model_data_assessment_2023'
OUT.mkdir(exist_ok=True)
DATA = OUT/'data'; DATA.mkdir(exist_ok=True)

def load(path):
    with gzip.open(path,'rb') as f: return pickle.load(f)

state_path = REF/'state_2023/actual_2023.pkl.gz'
state = load(state_path)
assert state['calendar_year'] == 2023
cache = OUT/'assessment_solution.pkl.gz'
if not cache.exists() or not hasattr(load(cache).solution,'V'):
    # Recover only this date's missing policies, using the original frozen engine
    # and saved continuation. No equilibrium or transition is re-estimated.
    sys.path.insert(0,str(REF/'source_snapshot/code/model/experiments/birth_count_choice'))
    from model.engine import solver as engine
    seed_path = ROOT/'output/model/transition_readiness_v1/current_baseline_20261003/local_continuation_v1/fit_job/run/candidate_0004/horizon_032/map_001/date_000/diagnostic_packet.pkl.gz'
    seed = load(seed_path)
    P = state['parameters']; b = state['b_grid']; price = np.array([state['forecast_prices'][0]])
    pre = P.birth_count_pre_distribution.copy(); post = P.birth_count_post_distribution.copy()
    first = P.birth_count_first_birth_tagged_distribution.sum(axis=(0,1,2,4,5,6))
    print('Recovering the 2023 policies at saved prices and continuation',flush=True)
    V,c,h,bp,choice,ten,loc,fert,fv,_ = engine.solve_bellman_full_markov_income(
        P.user_cost_rate*price + price-state['forecast_prices'][1],price,P,b,engine.precompute_shared(P,b),continuation_V=state['continuation_V'])
    error = float(np.max(np.abs(V-state['current_2023_V'])))
    np.testing.assert_allclose(V,state['current_2023_V'],rtol=2e-14,atol=1e-10,err_msg='2023 policy reproduction')
    import run_dynamic_population_transition as calendar
    maps = calendar.build_transition_maps(price,P,b,seed['shared'])
    g = engine.realize_current_cross_section(post,loc,choice,ten,maps.lmm_idx,maps.lmm_wt,maps.tmx_idx,maps.tmx_wt,
        use_compiled_scatter=bool(getattr(P,'use_numba_scatter',False)))
    s=SimpleNamespace(g=g,g_beginning_distribution=post,b_grid=b,p_eq=price,hR_pol=h,
        assessment_prebirth=pre,assessment_first_births=first, V=V,c_pol=c,bp_pol=bp,tenure_choice=choice,
        tenure_probs=ten,loc_probs=loc,fert_probs=fert,fert_value=fv,fert2_probs=getattr(P,'_fert2_probs',None),
        g_stay_distribution=getattr(P,'_g_stay_distribution',None))
    result=SimpleNamespace(P=P,solution=s)
    with gzip.open(cache,'wb') as f:pickle.dump(result,f)
    (OUT/'state_receipt.json').write_text(json.dumps(dict(checkpoint=str(state_path),checkpoint_sha256=sha256(state_path),
        calendar_year=2023,policy_value_max_abs_gap=error,transition_reestimated=False,policy_date_calls=1),indent=2)+'\n')
else:
    receipt=json.loads((OUT/'state_receipt.json').read_text())
    if receipt['checkpoint_sha256'] != sha256(state_path):raise RuntimeError('Saved 2023 checkpoint changed')
    result=load(cache)

psid=DATA/'psid_recent.csv'
if not psid.exists():
    subprocess.run(['/usr/local/bin/Rscript',str(ROOT/'code/data/psid_followup_mar2026/build_recent_assessment_extract.R'),str(DATA)],check=True)
acs=DATA/'housing_profile_by_age.csv'
if not acs.exists():
    spec=importlib.util.spec_from_file_location('housing_builder',ROOT/'code/empirical/housing/build_initial_housing_profile_diagnostic.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    data,fields,_=mod._open_sorted_source(); lo=mod._lower_bound(data,fields['year'],2023);hi=mod._lower_bound(data,fields['year'],2024)
    if hi==lo:raise RuntimeError('ACS source lacks 2023')
    acc={}
    for start in range(lo,hi,mod.CHUNK):
        block=data[start:min(start+mod.CHUNK,hi)]
        names=('year','sample','met2013','gq','pernum','relate','hhwt','age','ownershp','rooms','unitsstr','nchild','yngch','eldch')
        a={n:np.asarray(block[fields[n]]) for n in names}
        keep=(a['sample']==202301)&np.isin(a['gq'],(1,2))&(a['pernum']==1)&(a['relate']==1)&(a['hhwt']>0)&(a['age']>=18)&(a['age']<=85)&np.isin(a['ownershp'],(1,2))&(a['rooms']>0)
        frame=mod.pd.DataFrame({n:a[n][keep] for n in names})
        mod._aggregate_chunk(frame,set(),acc,None,None)
    mod._write_outputs(acc,DATA,dict(years=[2023],source=str(mod.SOURCE),source_size=mod.SOURCE.stat().st_size,source_mtime_ns=mod.SOURCE.stat().st_mtime_ns,
        builder=str(Path(mod.__file__)),builder_sha256=sha256(mod.__file__),target_system_changed=False))
prepare_assessment(result,OUT,output=OUT,data_year=2023,psid_path=psid,acs_path=acs)
print(OUT/'model_data_assessment.pdf')
