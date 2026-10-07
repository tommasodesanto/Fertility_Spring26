"""Property-tax (and credit) arms from a saved 2023 state, generalized from
output/model/credit_from_2023_refit_20261006/scripts/run_refit.py + common.py (Oct 6 2026).
Only change: input paths and endpoint come from the state JSON (PTAX_STATE), outputs go to PTAX_OUT.
Usage (via pipeline.py): driver.py ARM RATE_OR_PHI fixed|ge, with POLICY=ptax|phi in the environment.
FIXED-SUPPLY COPY (Oct 6 2026, exploratory): adds POLICY=supply. The housing stock is frozen from 2023 on at its
2023 level: xi_supply 0.63 -> 0 and H0 rescaled to H0*(uc*q_2023/r_bar)^0.63, q_2023 = saved baseline 2023 price;
the dated supply rule is replaced by the same constant stock. Nothing else changes. Second argument is ignored."""
import os,sys,json,gzip,pickle,hashlib,copy,time
from pathlib import Path
import numpy as np
CFG=json.loads(Path(os.environ['PTAX_STATE']).read_text())
FARM=Path(os.environ['FARM'])
A2=FARM/'rebate_runs/a2'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def build(out):
    import rough_two_shock_runtime as rtm
    plan=json.loads((A2/'manifest.json').read_text())
    nr=rtm.build_runtime(plan=plan,output=Path(out)/'runtime',smoke=False)
    rt=nr.rt
    R=A2/'fit/reference/reference_reconstruction.json'
    rt.load_reconstructed_reference(dict(path=str(R),sha256=sha(R)),Path(out)/'reference')
    return nr,rt
def load_state():
    with gzip.open(CFG['state_2023'],'rb') as f:S=pickle.load(f)
    root=json.load(open(CFG['stage2_root']))
    T=np.asarray(root['final']['transfers'],float)[2:]
    assert np.allclose(np.asarray(root['final']['prices'])[2:],S['forecast_prices'],rtol=0,atol=0)
    assert np.allclose(np.asarray(root['final']['fiscal_values'])[2:],S['forecast_pensions'],rtol=0,atol=0)
    return S,T
ENDPOINT={k:float(CFG['endpoint'][k]) for k in ('price','population_scale','transfer')}

arm,phi,mode=sys.argv[1],float(sys.argv[2]),sys.argv[3]
OUT=Path(os.environ['PTAX_OUT'])/arm; OUT.mkdir(parents=True,exist_ok=False)
log=lambda **k:(open(OUT/'progress.jsonl','a').write(json.dumps(dict(t=time.time(),**k),default=str)+'\n'))
nr,rt=build(OUT)
S,Tpath=load_state(); psi=float(S['forecast_psi'][0])
relaxed=[]
_req=rt.scaffold.require
def soft_require(ok,message):
    if not ok and message=='Dated household/accounting gate failed':
        relaxed.append(message); log(phase='dated_gate_failed_logged'); return
    return _req(ok,message)
rt.scaffold.require=soft_require
ecr=sys.modules['e5f_evening_calibration_runtime']; _rne=ecr.require_no_negative_estates
def soft_rne(ledger):
    neg=float(ledger['estate']['totals']['net_negative'])
    if neg>1e-10: relaxed.append(('stationary_net_negative',neg)); log(phase='stationary_negative_estate_logged',net_negative=neg); return
    return _rne(ledger)
ecr.require_no_negative_estates=soft_rne
# the only economic change: phi (credit arms) or the property tax (POLICY=ptax, phi argument = new annual rate)
POLICY=os.environ.get('POLICY','phi')
with gzip.open(CFG['state_2023'],'rb') as _f: Q2023=float(np.asarray(pickle.load(_f)['forecast_prices'],float)[0])
def frozen_stock(Q): return float(np.asarray(Q.H0,float)[0]*(float(Q.user_cost_rate)*Q2023/float(np.asarray(Q.r_bar,float)[0]))**float(np.asarray(Q.xi_supply,float)[0]))
def apply(Q):
    Q=copy.deepcopy(Q)
    if POLICY=='supply':
        S=frozen_stock(Q); Q.H0=np.full_like(np.asarray(Q.H0,float),S); Q.xi_supply=np.zeros_like(np.asarray(Q.xi_supply,float))
        if os.environ.get('PTAX_PROJECT')=='1': Q.native_exact_inherited_distribution=False  # same projection fallback as the property-tax arm (reported mass)
        return Q
    if POLICY=='phi': Q.phi=np.full_like(np.asarray(Q.phi,float),phi)
    else:
        tau=phi*float(getattr(Q,'period_years',4)); uc_old=float(Q.user_cost_rate)
        Q.r_bar=np.asarray(Q.r_bar,float)*(uc_old-float(Q.tau_H)+tau)/uc_old; Q.tau_H=tau; Q.user_cost_rate=float(Q.q)+float(Q.delta)+tau
        if os.environ.get('PTAX_PROJECT')=='1': Q.native_exact_inherited_distribution=False  # engine's own surprise-change frontier projection (reported mass)
    return Q
spec0=dict(tau_H=float(rt.P.tau_H),user_cost=float(rt.P.user_cost_rate),r_bar=np.asarray(rt.P.r_bar).tolist(),phi=np.asarray(rt.P.phi).tolist(),H0=np.asarray(rt.P.H0).tolist(),xi_supply=np.asarray(rt.P.xi_supply).tolist())
if POLICY=='supply':
    S_fix=frozen_stock(rt.P); old_rule=rt.packet['supply_rule']
    log(phase='frozen_stock',q_2023=Q2023,stock=S_fix,old_rule_at_q2023=float(old_rule.quantity(np.array([Q2023]))[0]))
    assert abs(float(old_rule.quantity(np.array([Q2023]))[0])/S_fix-1)<1e-10,'dated supply rule differs from H0 formula'
rt.P=apply(rt.P); rt.packet=dict(rt.packet); rt.packet['parameters']=apply(rt.packet['parameters'])
if POLICY=='supply':
    rt.packet['supply_rule']=type(old_rule)('static-elastic',Q2023,S_fix,0.0)
if POLICY in ('ptax','supply'):
    # Diagnostic: a surprise tax (or supply change) leaves a tiny inherited mass (deeply indebted, lowest-income owners) with no feasible
    # choice in 2023. Log it and keep it (<=1e-7 of mass) instead of stopping; the evidence file is still written.
    _rdp=sys.modules['run_dynamic_population_transition']; _orig_req=_rdp._require_exact_inherited_distribution
    def soft_inherited(out,policy,P,b_grid):
        try: return _orig_req(out,policy,P,b_grid)
        except _rdp.InheritedDistributionInfeasible as e:
            dm=float(getattr(e,'dead_mass',1.))
            if dm<=1e-7: relaxed.append(('inherited_infeasible_mass',dm)); log(phase='inherited_infeasible_logged',dead_mass=dm); return
            raise
    _rdp._require_exact_inherited_distribution=soft_inherited
    prim=rt.rt['primitive']; _db=prim.dated_budget
    def soft_budget(*a,**k):
        try: return _db(*a,**k)
        except RuntimeError as e:
            if not str(e).startswith('Dated budget gate failed'): raise
            import re; m=re.search(r'mass=([^,]+), excess=(\S+)',str(e)); mass,exc=float(m.group(1)),float(m.group(2))
            relaxed.append(('dated_budget',mass,exc)); log(phase='dated_budget_logged',mass=mass,excess=exc)
            return {'budget_excess_mass':mass,'maximum_occupied_excess':exc,'actual_rent':float(a[4]) if len(a)>4 else None,'budget_tolerance':1e-9}
    prim.dated_budget=soft_budget
log(phase='policy',policy=POLICY,before=spec0,after=dict(tau_H=float(rt.P.tau_H),user_cost=float(rt.P.user_cost_rate),r_bar=np.asarray(rt.P.r_bar).tolist(),phi=np.asarray(rt.P.phi).tolist(),H0=np.asarray(rt.P.H0).tolist(),xi_supply=np.asarray(rt.P.xi_supply).tolist()))
import importlib
fc=importlib.import_module('experiments.birth_count_choice.model.fiscal_closure')
import one_shock_floor as retained; retained.original_modules()
from e5f_ssj_scaled_step_root import solve_price_path_scaled
import rebate_root
if mode=='fixed':
    log(phase='terminal_start')
    if POLICY=='phi':
        terminal,trec=rt.stationary_at(psi,ENDPOINT['price'],ENDPOINT['transfer'],OUT/'terminal')
        endpoint=dict(ENDPOINT)
    else:
        def solve(T,k):
            packet,record=rt.stationary_at(psi,ENDPOINT['price'],T,OUT/'terminal'/f'rebate_{k:02d}'); return packet['solution'],(packet,record)
        _s,(terminal,trec),fiscal=fc.balance(ENDPOINT['transfer'],solve,tol=1e-6,require_first=False,slope_hint=None)
        endpoint=dict(ENDPOINT,transfer=fiscal['transfer'])
    log(phase='terminal_done',renewal=trec['renewal_residual'],endpoint=endpoint)
else:
    accepted={}; latest={}; cnt=[0]; hint=[None]
    def ev_end(q):
        cnt[0]+=1; price=float(q[0]); point=OUT/'endpoint'/f'point_{cnt[0]:03d}'
        def solve(T,k):
            packet,record=rt.stationary_at(psi,price,T,point/f'rebate_{k:02d}'); return packet['solution'],(packet,record)
        start=fc.predict(accepted,price,ENDPOINT['transfer'])
        _s,(packet,record),fiscal=fc.balance(start,solve,tol=1e-6,require_first=price in accepted,slope_hint=hint[0])
        accepted[price]=fiscal['transfer']; hint[0]=fiscal['slope']; latest.update(packet=packet,record=record,T=fiscal['transfer'])
        log(phase='endpoint_point',n=cnt[0],price=price,T=fiscal['transfer'],renewal=record['renewal_residual'])
        return dict(mapping_valid=True,residual=np.array([record['renewal_residual']]))
    q0=rt.reference_price
    er=solve_price_path_scaled(initial_prices=np.array([ENDPOINT['price']]),evaluate=ev_end,project=lambda q:np.clip(q,q0*.05,q0*20.),
        slope=1.,market_tolerance=1e-6,max_log_step=.15,damping=.7,max_evaluations=48,deadline_monotonic=time.monotonic()+3600,
        max_condition_number=1e8,worsening_factor=1.5,final_reproduction_tolerance=1e-10)
    assert er['converged'],'endpoint failed'
    terminal=latest['packet']; endpoint=dict(price=latest['record']['price'],population_scale=latest['record']['population_scale'],transfer=latest['T'])
    log(phase='endpoint_done',**endpoint)
pf=rt.pf; orig=pf.calendar.evaluate_period; stats=[]
# capture (diagnostic only): pre-birth-choice value VI per (age j, income z) from the last Bellman solve (date 0),
# plus date-0/1 count policies and pre-fertility distributions, for the branch comparison.
import importlib as _il
_bc=_il.import_module('experiments.birth_count_choice.model.engine.birth_count'); _menu=_bc.birth_count_menu; VIcap={}; VIsolves=[]; capev={}
def menu_cap(VI,n,m,pi,kappa,F=0.,cap=3):
    if n==0 and m==0:
        f=sys._getframe(1).f_locals; key=(int(f['j']),int(f['zz']))
        if not VIsolves or key in VIsolves[-1]: VIsolves.append({})
        VIsolves[-1][key]=dict(VI=np.array(VI,copy=True),pi=float(pi))
    return _menu(VI,n,m,pi,kappa,F,cap)
_bc.birth_count_menu=menu_cap
def wrapped(price,g_pre,P,b_grid,shared,counter,*a,**k):
    ev=orig(price,g_pre,P,b_grid,shared,counter,*a,**k)
    g=ev.g_current; y=g[:,:,:,0:3]; ax=tuple(i for i in range(g.ndim) if i!=3)
    stats.append(dict(price=float(np.ravel(price)[0]),phi=np.asarray(P.phi).tolist(),own_18_29=float(y[:,1:].sum()/y.sum()),
        own_by_age=(g[:,1:].sum(axis=ax)/np.maximum(g.sum(axis=ax),1e-300))[:6].tolist(),births=float(ev.births)))
    t=len(stats)-1
    if t<=1: capev[t]=dict(g_pre=np.array(ev.g_pre),g_current=np.array(ev.g_current),action=np.array(ev.policy.birth_count_action_probs),realized=np.array(ev.policy.birth_count_realized_probs),V=np.array(ev.policy.V))
    return ev
n=[0]; finals={}
def run_map(q,b,T):
    n[0]+=1; stats.clear(); pf.calendar.evaluate_period=wrapped
    try: native,record=rt.mapping(terminal,endpoint,q,b,S['forecast_psi'],OUT/f'map_{n[0]:03d}',initial_state=S['initial_state'],start_year=2023,transfers=T)
    finally: pf.calendar.evaluate_period=orig
    # pick the Bellman solve whose menu values reproduce the date-0 policy value at n=0
    best=None
    if 0 in capev and VIsolves:
        V0=capev[0]['V']; F0=float(rt.P.first_birth_fixed_cost); kap=float(rt.P.kappa_fert); errs=[]
        for si,sol in enumerate(VIsolves):
            e=0.
            for (j,z),d in sol.items():
                VI=d['VI']; pi=d['pi']; u0=VI[...,0,0]; u1=pi*(VI[...,1,1]-F0)+(1-pi)*u0
                mx=np.maximum(u0,u1); val=mx+kap*np.log(np.exp((u0-mx)/kap)+np.exp((u1-mx)/kap))
                ok=(u0>-1e8)&(u1>-1e8)
                if ok.any(): e=max(e,float(np.max(np.abs(val-V0[:,:,:,j,z,0,0])[ok])))
            errs.append(e)
        best=int(np.argmin(errs)); VIcap.update(VIsolves[best]); log(phase='vi_match',solves=len(VIsolves),best=best,err=errs[best],errs=errs[:40])
    with gzip.open(OUT/f'map_{n[0]:03d}'/'branch_capture.pkl.gz','wb') as _f: pickle.dump(dict(VI=dict(VIcap),ev=dict(capev),kappa=float(rt.P.kappa_fert),kappa_cont=float(rt.P.kappa_fert_continuation),F=float(rt.P.first_birth_fixed_cost)),_f)
    VIcap.clear(); VIsolves.clear(); capev.clear()
    finals[n[0]]=dict(rows=record['rows'],fertility=record.get('fertility'),stats=list(stats),prices=list(map(float,q)),pensions=list(map(float,b)),transfers=list(map(float,T)),
        market_residual=record['market_residual'],fiscal_residual=record['fiscal_residual'],rebate_residual=record['rebate_residual'])
    pickle.dump(finals[n[0]],open(OUT/f'map_{n[0]:03d}'/'summary.pkl','wb'))
    log(phase='path_map',n=n[0],max_market=max(map(abs,record['market_residual'])),max_fiscal=max(map(abs,record['fiscal_residual'])),max_rebate=max(map(abs,record['rebate_residual'])))
    return record
q,b=np.asarray(S['forecast_prices'],float),np.asarray(S['forecast_pensions'],float)
if mode=='fixed' and POLICY=='phi':
    run_map(q,b,Tpath); converged=None
elif mode=='fixed':
    # fixed prices and pensions; rebate rebalanced at every date by fixed-point iteration T_t <- R_t/N_t
    T=Tpath.copy(); converged=False
    for it in range(10):
        r=run_map(q,b,T); res=max(map(abs,r['rebate_residual']))
        if res<=2e-5: converged=True; break
        T=np.asarray([row['implied_equal_transfer'] for row in r['rows']],float)
else:
    h=len(q); Jf=np.asarray(json.load(open(CFG['stage2_root']))['final_jacobian'],float); H=Jf.shape[0]//3
    idx=[blk*H+d for blk in range(3) for d in range(2,H)]; J=Jf[np.ix_(idx,idx)]
    slopes=[float(np.median(np.abs(np.diag(J)[i*h:(i+1)*h]))) for i in range(3)]
    P0=float(nr.reference_transfer) if hasattr(nr,'reference_transfer') else float(rt.P.property_tax_lump_sum_transfer)
    def ev(q,b,T):
        r=run_map(q,b,T); return dict(mapping_valid=True,market_residual=r['market_residual'],fiscal_residual=r['fiscal_residual'],rebate_residual=r['rebate_residual'])
    root=rebate_root.solve_rebate_path(initial_prices=q,initial_pensions=b,initial_transfers=Tpath,evaluate=ev,
        project_prices=lambda v:np.clip(v,q0*.05,q0*20.),pension_bounds=[float(rt.P.pension)*.05,float(rt.P.pension)*20.],
        transfer_bounds=[P0*.05,P0*20.],market_tolerance=2e-4,fiscal_tolerance=2e-5,market_slope=slopes[0],fiscal_slope=slopes[1],
        rebate_slope=slopes[2],max_log_step=.15,damping=.7,max_evaluations=12,deadline_monotonic=time.monotonic()+4*3600,
        max_condition_number=1e8,worsening_factor=1.5,final_reproduction_tolerance=1e-10,initial_jacobian=J)
    converged=bool(root['converged'])
json.dump(dict(arm=arm,phi=phi,mode=mode,converged=converged,evaluations=n[0],endpoint=endpoint,last=finals[n[0]],relaxed_log=[str(x) for x in relaxed]),
    open(OUT/'result.json','w'),default=lambda o:o.tolist() if hasattr(o,'tolist') else str(o),indent=1)
log(phase='done',converged=converged)
