"""Pure, zero-solve input construction for three experimental entry laws."""
from __future__ import annotations
import copy, hashlib, json, math
from pathlib import Path
from types import SimpleNamespace
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
GRID=ROOT/'output/model/publication_refactor_20260929/grid_resolution_v1'
PARAMETERS=('beta_annual','chi','first_birth_fixed_cost','kappa_fert','kappa_fert_continuation','theta0','delta_alpha_jump','child_benefit_curvature','tenure_choice_kappa')
ARMS={'empirical_credit':.25,'nonnegative_mean':0.}
PLAN=json.loads((HERE/'plan.json').read_text())
LANES=PLAN['lanes']
TABLE=ROOT/'output/model/fixed_reference_economics_20260928/entry_ratio_comparison_v1/collected/full/candidate_five_ratios/phase_b_ge/selected_root/parameters.csv'

def require(ok,message):
    if not ok:raise RuntimeError(message)

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def canonical(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def decode(v,a):
    if isinstance(v,list):return [decode(x,a) for x in v]
    if isinstance(v,dict):
        if '__ndarray__' in v:return a[v['__ndarray__']].copy()
        if '__npscalar__' in v:return np.dtype(v['__npscalar__']).type(v['value'])
        if '__tuple__' in v:return tuple(decode(x,a) for x in v['__tuple__'])
        if '__float__' in v:return float(v['__float__'])
        if '__dict__' in v:return {k:decode(x,a) for k,x in v['__dict__'].items()}
        raise RuntimeError('Unknown bundle tag')
    return v

def proposal(lane='nonnegative_mean_120x9'):
    require(lane in LANES,'Unknown lane')
    config=LANES[lane];source=PLAN['bundle_sources'][config['size']]
    folder=ROOT/source['folder'];meta=json.loads((folder/'bundle.json').read_text())
    require(sha(folder/'bundle.json')==source['bundle_sha256'],'Bundle metadata drift')
    require(meta['arrays_sha256']==source['arrays_sha256'] and sha(folder/'arrays.npz')==source['arrays_sha256'],'Grid arrays drift')
    # Fine lane loads the exact original bundle primitives; no interpolation of coarse inputs.
    with np.load(folder/'arrays.npz',allow_pickle=False) as a:
        P=SimpleNamespace(**{k:decode(v,a) for k,v in meta['parameters'].items() if not k.startswith('_')});grid=a['b_grid'].copy()
    require([P.Nb,P.Nz]==config['dimensions'] and grid.size==P.Nb,'Wrong lane grid')
    return P,grid

def linear_projection(grid,points,weights):
    """Exactly the native linear probability projection, including boundary clipping."""
    mass=np.zeros(len(grid))
    for value,weight in zip(points,weights):
        require(math.isfinite(float(value)) and weight>=0,'Invalid entry draw')
        if weight==0:continue
        b=float(np.clip(value,grid[0],grid[-1]));hi=int(np.searchsorted(grid,b,side='left'))
        if hi<=0:mass[0]+=float(weight)
        elif hi>=len(grid):mass[-1]+=float(weight)
        elif abs(grid[hi]-b)<=1e-14:mass[hi]+=float(weight)
        else:
            lo=hi-1;s=(b-grid[lo])/max(float(grid[hi]-grid[lo]),1e-12)
            mass[lo]+=float(weight)*(1-s);mass[hi]+=float(weight)*s
    require(mass.sum()>0,'Empty entry law')
    return mass/mass.sum()

def entry(P,grid,arm):
    require(arm in ARMS,'Unknown arm');P=copy.deepcopy(P)
    require(P.property_tax_lump_sum_transfer==0 and P.preference_spec=='eqscale','Entry cash formula contract changed')
    ratios=np.asarray(P.entry_wealth_ratio_nodes,float);weights=np.asarray(P.entry_wealth_ratio_weights,float);weights=weights/weights.sum()
    require(len(ratios)==5,'Exactly five empirical bin means required')
    original_ratios=ratios.copy();lam=float((ratios@weights)/(np.maximum(ratios,0)@weights))
    if arm=='zero_wealth':ratios=np.zeros_like(ratios)
    elif arm=='nonnegative_mean':ratios=lam*np.maximum(ratios,0)
    yperiod=np.asarray(P.income[0,0]*P.z_grid)
    annual=yperiod/float(P.period_years)/(1-float(P.tau_pay))
    points=ratios[:,None]*annual[None,:]
    C=np.column_stack([linear_projection(grid,points[:,z],weights) for z in range(P.Nz)])
    require(np.max(abs(C.sum(axis=0)-1))<2e-15,'Entry columns fail mass')
    zero=np.flatnonzero(grid==0.)
    require(len(zero)==1,'Exact zero wealth node required')
    if arm!='empirical_credit':require(float((C*(grid[:,None]<0)).sum())==0.,'Nonnegative arm has debt')
    if arm=='zero_wealth':require(np.array_equal(C[zero[0]],np.ones(P.Nz)),'Zero law spreads entry')
    P.fixed_reference_entry_conditional=C;P.native_fixed_reference_entry=True
    # Keep primitive nodes consistent with the operative conditional matrix.
    P.entry_wealth_ratio_nodes=ratios.copy()
    credit=ARMS[arm]
    slack=P.R_gross*grid[:,None]+yperiod[None,:]+credit
    occupied=C>0;minslack=float(slack[occupied].min())
    require(minslack>1e-6,'Necessary entrant current budget fails')
    draw=weights[:,None]*P.z_weights[None,:]
    clipped=(points<grid[0])|(points>grid[-1])
    rawmean=float((points*draw).sum());projected=float(grid@(C@P.z_weights))
    require(abs(rawmean-projected)<1e-12,'Entry mean changed by projection')
    originalmean=float((original_ratios@weights)*(annual@P.z_weights))
    if arm=='nonnegative_mean':require(abs(rawmean-originalmean)<1e-12,'Censor/scale did not preserve raw five-node mean')
    report=dict(arm=arm,credit=credit,dimensions=[int(P.Nb),int(P.Nz)],raw_ratio_nodes=original_ratios.tolist(),effective_ratio_nodes=ratios.tolist(),ratio_weights=weights.tolist(),lambda_positive=lam,mean_wealth=projected,raw_mean_wealth=rawmean,mean_annual_entry_income=float(annual@P.z_weights),negative_wealth_share=float((C@P.z_weights)[grid<0].sum()),zero_wealth_share=float((C@P.z_weights)[grid==0].sum()),clipped_draw_mass=float(draw[clipped].sum()),necessary_current_budget_minimum_slack=minslack,conditional_sha256=hashlib.sha256(C.tobytes()).hexdigest(),approximation='Existing five empirical quintile-bin means; censoring occurs after binning, not on raw survey observations.',lifecycle_solves=0)
    return P,report

def seed_and_bounds(lane='nonnegative_mean_120x9'):
    import csv
    rows=PLAN['reference_parameter_table']
    keyed={r['parameter']:r for r in rows}
    require(lane in LANES,'Unknown lane');seed=dict(LANES[lane]['seed'])
    bounds={k:[float(keyed[k]['lower']),float(keyed[k]['upper'])] for k in PARAMETERS}
    return seed,bounds,rows

def check_point(point,bounds):
    require(set(point)==set(PARAMETERS),'Wrong nine free coordinates')
    for key,value in point.items():require(math.isfinite(value) and bounds[key][0]<=value<=bounds[key][1],'Out of bounds '+key)

def bind(P,point,bounds):
    check_point(point,bounds);Q=copy.deepcopy(P)
    for key,value in point.items():
        if key!='beta_annual':setattr(Q,key,float(value))
    Q.beta=float(point['beta_annual'])**float(Q.period_years)
    Q.rho=1/Q.beta-1;Q.rho_hat=Q.rho
    Q.eps_fert=Q.kappa_fert
    # This direct binding intentionally never calls broad apply_overrides,
    # which would reconstruct fiscal/income objects outside this experiment.
    return Q
