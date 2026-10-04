"""Bounded cached-only measurement audit. No solve, target edit or download."""
from pathlib import Path
import csv, hashlib, importlib.util, json, os
for key in ('NUMBA_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RAW = ROOT / 'code/data/cps_fertility/cache/jun24pub.csv'
REC = ROOT / 'output/model/credit_mechanism_20261004/evidence/preserved/fertility_by_income.json'
HELPER = HERE.parent / 'diagnostics/extract_common_states.py'
BASE = ROOT / 'output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/fable_analysis/credit_mechanism/credit_relaxation/phi_080'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
manifest=json.loads((ROOT/'code/data/cps_fertility/source_manifest.json').read_text())
assert sha(RAW)==manifest['sha256']
rows=[]
with RAW.open(newline='') as f:
    for r in csv.DictReader(f):
        if r['PESEX']!='2': continue
        age,n,w,inc,rel=(int(r['PRTAGE']),int(r['PTSF1']),float(r['PWSSWGT']),int(r['HEFAMINC']),int(r['PRFAMREL']))
        if (24<=age<=26 or 40<=age<=44) and 0<=n<=5 and w>0:
            rows.append((age,n,w,inc,rel,r['HRHHID']+':'+r['HRHHID2']))
A=np.array([r[:5] for r in rows],float)
cluster=np.unique([r[5] for r in rows],return_inverse=True)[1]

def third_alloc(m):
    """Fractionally split tied income categories, never sort within a tie."""
    m=np.asarray(m,float); total=m.sum(); lo=np.cumsum(m)-m; hi=np.cumsum(m)
    return np.array([np.maximum(0,np.minimum(hi,(k+1)*total/3)-np.maximum(lo,k*total/3))/np.maximum(m,1e-300) for k in range(3)])

def empirical(mask,kind):
    valid=mask&(A[:,3]>=1)&(A[:,3]<=16)
    w=A[:,2]*valid; cats=A[:,3].astype(int)
    if kind=='fixed_bins':
        allocation=np.array([cats<=10,(cats>=11)&(cats<=14),cats>=15],float)
    else:
        mass=np.array([w[cats==i].sum() for i in range(1,17)])
        fractions=third_alloc(mass)
        allocation=np.zeros((3,len(A)))
        for i in range(1,17): allocation[:,cats==i]=fractions[:,i-1,None]
    out=[]; scores=[]
    for k in range(3):
        wk=w*allocation[k]; denom=wk.sum(); x=np.minimum(A[:,1],3); mean=np.sum(wk*x)/denom
        score=np.bincount(cluster,weights=wk*(x-mean)/denom,minlength=cluster.max()+1)
        scores.append(score)
        out.append(dict(group=k+1,n_contributing=int((wk>0).sum()),weighted_share=float(denom/w.sum()),
                        mean_cap3=float(mean),mean_cap5=float(np.sum(wk*A[:,1])/denom),
                        childless=float(np.sum(wk*(A[:,1]==0))/denom),
                        two_plus=float(np.sum(wk*(A[:,1]>=2))/denom),
                        approximate_household_cluster_se=float(np.linalg.norm(score))))
    diff=out[2]['mean_cap3']-out[0]['mean_cap3']; se=float(np.linalg.norm(scores[2]-scores[0]))
    return dict(groups=out,high_minus_low=diff,approximate_contrast_se=se,
                approximate_interval95=[diff-1.96*se,diff+1.96*se])

emp={}
for label,mask in [('age24_26',(A[:,0]>=24)&(A[:,0]<=26)),('age40_44',A[:,0]>=40)]:
    for sample,sel in [('all_women',mask),('reference_or_spouse',mask&np.isin(A[:,4],[1,2]))]:
        for grouping in ('fixed_bins','fractional_weighted_thirds'):
            emp[f'{label}_{sample}_{grouping}']=empirical(sel,grouping)

spec=importlib.util.spec_from_file_location('cached_birth_map',HELPER); helper=importlib.util.module_from_spec(spec); spec.loader.exec_module(helper)
P,arrays,pre,rate,ages,fec,checks=helper.load(BASE)
post=arrays['g_beginning_distribution']
def zn(g,j): return g[:,:,:,j].sum(axis=(0,1,2,5))
def model_tab(q,kind):
    m=q.sum(1)
    if kind=='scratch_midpoint_thirds':
        group=np.digitize((np.cumsum(m)-m/2)/m.sum(),[1/3,2/3]); alloc=np.array([group==k for k in range(3)],float)
    else: alloc=third_alloc(m)
    out=[]
    for k in range(3):
        qk=(q*alloc[k,:,None]).sum(0); d=qk.sum()
        out.append(dict(group=k+1,weighted_share=float(d/m.sum()),mean_cap3=float(qk@np.arange(4)/d),childless=float(qk[0]/d)))
    return dict(groups=out,high_minus_low=out[-1]['mean_cap3']-out[0]['mean_cap3'])
model={}
for label,j in [('end22_25',1),('end42_45',6)]:
    for grouping in ('scratch_midpoint_thirds','fractional_weighted_thirds'):
        model[f'{label}_{grouping}']=model_tab(zn(post,j),grouping)
for label,age_values in [('uniform_interview_age24_26',range(24,27)),('uniform_interview_age40_44',range(40,45))]:
    qs=[]
    for age in age_values:
        midpoint=age+.5; j=int((midpoint-P['age_start'])//P['da']); f=(midpoint-ages[j])/P['da']
        q=(1-f)*zn(pre,j)+f*zn(post,j); qs.append(q/q.sum())
    for grouping in ('scratch_midpoint_thirds','fractional_weighted_thirds'):
        model[f'{label}_{grouping}']=model_tab(sum(qs)/len(qs),grouping)

original=json.loads(REC.read_text()); errors=[]
for label,orig_key in [('age24_26','age24_26_income'),('age40_44','age40_44_income')]:
    errors.extend(abs(x['mean_cap3']-y['mean_ceb_cap3']) for x,y in zip(emp[label+'_all_women_fixed_bins']['groups'],original['data_cps_june2024'][orig_key]))
for label,orig_key in [('end22_25','end_of_22_25_cell'),('end42_45','cell_42_45')]:
    errors.extend(abs(x['mean_cap3']-y['mean_ceb_cap3']) for x,y in zip(model[label+'_scratch_midpoint_thirds']['groups'],original['model_chain13'][orig_key]['by_tercile']))
payload=dict(status='cached_only_diagnostic_not_calibration',source=dict(cps=str(RAW),cps_sha256=sha(RAW),model=str(BASE),model_sha256=sha(BASE/'solution_arrays.npz'),helper_sha256=sha(HELPER),recovered_summary=str(REC)),
             checks=checks,maximum_pasted_mean_reproduction_error=float(max(errors)),empirical=emp,model=model,
             caveats=['Approximate household-cluster SE conditions on original weights; no CPS strata, PSU or replicate-weight design correction.',
                      'Weighted thirds fractionally split tied discrete income categories, using the same fraction for every woman in a tied category.',
                      'Uniform age diagnostic mixes exact prebirth/postbirth stock within each four-year cell; assumes uniform birth timing and equal exposure to each completed age.',
                      'Model income ranks current persistent Markov z; CPS ranks last-12-month combined family money income.'])
(HERE/'tables.json').write_text(json.dumps(payload,indent=2)+'\n')
print('pasted_max_error',payload['maximum_pasted_mean_reproduction_error'])
for domain,tables in [('CPS',emp),('model',model)]:
    for name,t in tables.items(): print(domain,name,[round(g['mean_cap3'],3) for g in t['groups']], 'H-L',round(t['high_minus_low'],3), 'SE',round(t.get('approximate_contrast_se',0),3))
