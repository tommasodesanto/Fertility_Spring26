"""Local diagnostic derivatives and trial proposals; never activate SMM targets.

Requires all 19 verified cases. Writes only proposals and local linear predictions;
a new nonlinear solve is required before any predicted improvement is a result.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.optimize import lsq_linear

HERE=Path(__file__).resolve().parent
NAMES=['beta_annual','kappa_fert','kappa_fert_continuation','chi','H0','theta0','theta1','first_birth_fixed_cost','h_P']

def value_at(data,path):
    for key in path.split('.'): data=data[key]
    return data

def transform(name,x): return math.log(-math.log(x)) if name=='beta_annual' else math.log(x)
def inverse(name,x): return math.exp(-math.exp(x)) if name=='beta_annual' else math.exp(x)
def write_csv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    raw=json.loads((HERE/'raw_case_summary.json').read_text())
    plan=json.loads((HERE/'cases.json').read_text())
    if raw['status']!='complete' or raw['case_count']!=19 or len(raw['cases'])!=19:
        raise RuntimeError('All 19 verified directions are required before a full Jacobian or joint proposal')
    cases={x['case_id']:x for x in raw['cases']}
    if set(cases)!={x['case_id'] for x in plan}:raise RuntimeError('Missing or duplicate planned case')
    for row in plan:
        c=cases[row['case_id']]
        if c['contract_sha256']!=row['contract_sha256'] or not c['numerical_gates_verified']:
            raise RuntimeError('Case provenance or numerical gate mismatch')
    fit=list(csv.DictReader((HERE.parent/'initial_fit_readout/target_fit.csv').open()))
    missing=[r['restriction_id'] for r in fit if r['role']!='normalization' and r['model_available'].lower()=='false']
    if missing!=['recent_parent_ownership']:raise RuntimeError('Unexpected missing empirical observations')
    scored=[r for r in fit if r['restriction_id'] not in ['initial_normalization','recent_parent_ownership']]
    if len(scored)!=11:raise RuntimeError('Explicit 11-observable diagnostic restriction set required')
    target=np.array([float(r['target']) for r in scored])
    params=list(csv.DictReader((HERE.parent/'initial_fit_readout/parameters.csv').open()))[:9]
    if [r['parameter'] for r in params]!=NAMES:raise RuntimeError('Unexpected coordinate order')
    base=cases['baseline'];u0=np.array([transform(n,base['parameters'][n]) for n in NAMES])
    def observe(c,projection):
        return np.array([value_at(c['early_measurement'],r['model_observation'].replace('uniform_birth_time',projection)) for r in scored],dtype=float)
    matrices={}; derivatives=[]
    for projection in ['uniform_birth_time','constant_post_cell']:
        m0=observe(base,projection); J=np.empty((11,9)); asymmetry=[]
        for k,n in enumerate(NAMES):
            minus,plus=cases[n+'_minus'],cases[n+'_plus']
            um,up=transform(n,minus['parameters'][n]),transform(n,plus['parameters'][n])
            mm,mp=observe(minus,projection),observe(plus,projection)
            J[:,k]=(mp-mm)/(up-um)
            jm,jp=(m0-mm)/(u0[k]-um),(mp-m0)/(up-u0[k])
            asymmetry.append(float(np.linalg.norm(jp-jm)/max(np.linalg.norm(J[:,k]),1e-12)))
            for i,row in enumerate(scored):
                derivatives.append(dict(projection=projection,restriction=row['restriction_id'],parameter=n,
                    coordinate='log annual discount rate' if n=='beta_annual' else 'log parameter',
                    centered_derivative=J[i,k],minus_derivative=jm[i],plus_derivative=jp[i],
                    base_moment=m0[i],target=target[i],base_gap=m0[i]-target[i]))
        if not np.isfinite(J).all():raise RuntimeError('Nonfinite derivative')
        matrices[projection]=dict(J=J,m0=m0,relative_column_asymmetry=asymmetry)
    references=list(csv.DictReader((HERE.parent/'initial_fit_readout/reference_precision_options.csv').open()))
    def reference_scale(restriction,option):
        matches=[r for r in references if r['restriction_id']==restriction and r['option']==option]
        if len(matches)!=1:raise RuntimeError('Reference scale must have one explicit source row')
        return float(matches[0]['scale'])
    fertility_reference=[reference_scale(r['restriction_id'],
        'pooled_GVF_correlation_one' if r['restriction_id'].startswith('cps_') else 'annual_temporal_dispersion')
        for r in scored[:4]]
    scales={
        'equal_relative_5pct':.05*np.abs(target),
        'reference_precision_diagnostic':np.array([
            *fertility_reference,
            *[float(r['empirical_standard_error'] or r['synthetic_scale']) for r in scored[4:]]])}
    if any(np.any(s<=0) or not np.isfinite(s).all() for s in scales.values()):raise RuntimeError('Invalid diagnostic scale')
    numerical_bounds=[]
    for n,r in zip(NAMES,params):
        endpoints=[transform(n,max(float(r['lower']),1e-8)),transform(n,float(r['upper']))]
        numerical_bounds.append((min(endpoints),max(endpoints)))
    trust=np.array([1.,.5,.5,.08,.5,.5,.5,.5,.5])
    lo=np.maximum(np.array([b[0] for b in numerical_bounds])-u0,-trust)
    hi=np.minimum(np.array([b[1] for b in numerical_bounds])-u0,trust)
    proposals=[]; spectra=[]
    for profile,scale in scales.items():
        for projections in [['uniform_birth_time'],['constant_post_cell'],['uniform_birth_time','constant_post_cell']]:
            A=np.vstack([matrices[p]['J']/scale[:,None] for p in projections])/math.sqrt(len(projections))
            r=np.concatenate([(matrices[p]['m0']-target)/scale for p in projections])/math.sqrt(len(projections))
            s=np.linalg.svd(A,compute_uv=False)
            spectra.append(dict(profile=profile,projections=projections,singular_values=s.tolist(),
                numerical_rank=int(np.linalg.matrix_rank(A)),condition_number=float(s[0]/s[-1]),
                caveat='Local numerical rank is not an identification or empirical-validity certificate'))
            for strength in [.1,1.,10.]:
                penalty=strength*float(np.median(s))**2
                solution=lsq_linear(np.vstack([A,np.eye(9)*math.sqrt(penalty)]),
                    np.r_[-r,np.zeros(9)],bounds=(lo,hi),tol=1e-11,max_iter=200)
                if not solution.success or not np.isfinite(solution.x).all():raise RuntimeError('Bounded linear diagnostic failed')
                for fraction in [.25,.5,1.]:
                    step=fraction*solution.x
                    theta={n:inverse(n,u0[k]+step[k]) for k,n in enumerate(NAMES)}
                    predictions={p:{row['restriction_id']:float(matrices[p]['m0'][i]+matrices[p]['J'][i]@step)
                        for i,row in enumerate(scored)} for p in matrices}
                    proposals.append(dict(profile=profile,projections=projections,ridge_strength=strength,
                        fraction=fraction,structural_candidate=theta,log_coordinate_step=step.tolist(),
                        local_base_score=float(r@r),linear_predicted_score=float(np.linalg.norm(r+A@step)**2),
                        predicted_moments=predictions,missing_recent_parent_model=None,
                        status='linear diagnostic proposal; requires nonlinear verification',calibrated_smm=False))
    result=dict(schema='e5f_initial_local_diagnostic_v1',calibrated_smm=False,actual_SMM_weights=None,
        observed_restrictions=11,proposed_scored_restrictions=12,structural_coordinates=9,
        missing_restrictions=missing,parameters=NAMES,spectra=spectra,
        asymmetry={p:v['relative_column_asymmetry'] for p,v in matrices.items()},
        diagnostic_scales={k:v.tolist() for k,v in scales.items()},
        scale_disclosure='Equal-relative profile uses synthetic 5% target scales. Reference profile uses unadopted CPS correlation-one approximations, NCHS annual temporal SD (not SE), other recorded SEs and synthetic bequest scale. Neither activates a target contract.',
        proposals=proposals,
        raw_panel_sha256=hashlib.sha256((HERE/'raw_case_summary.json').read_bytes()).hexdigest(),
        source_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/'local_analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    write_csv(HERE/'local_derivatives.csv',derivatives)
    print(json.dumps(dict(status='complete_diagnostic_analysis',proposals=len(proposals),spectra=spectra),indent=2))

if __name__=='__main__':main()
