"""Production stationary renewal closure with the exact native root and gates."""
from __future__ import annotations
import copy, json, time
from pathlib import Path
import numpy as np
from .inputs import validate_inputs
from . import native_phase_b as ge
from .engine import solver
from .reporting import build_context, historical_reporting_dependency

class Budget:
    """Native lifecycle/stage budget, without a calibration optimizer."""
    stage_deadline_seconds=300
    def __init__(self,out,deadline,max_lifecycle):
        self.out=Path(out);self.deadline_epoch=float(deadline);self.max_lifecycle=int(max_lifecycle)
        self.used_lifecycle=0;self.phase_a_lifecycle=0
        self.progress('initialized')
    @property
    def remaining_lifecycle(self): return self.max_lifecycle-self.used_lifecycle
    def claim_lifecycle(self,label):
        if self.used_lifecycle>=self.max_lifecycle: raise RuntimeError('Lifecycle cap reached')
        if time.time()+1>=self.deadline_epoch: raise RuntimeError('Global deadline reached')
        self.used_lifecycle+=1;self.progress('lifecycle_claimed',label=label)
    def progress(self,status,**extra):
        row=dict(status=status,time_epoch=time.time(),deadline_epoch=self.deadline_epoch,lifecycle_used=self.used_lifecycle,lifecycle_remaining=self.remaining_lifecycle,phase_a_lifecycle=0,**extra)
        self.out.mkdir(parents=True,exist_ok=True)
        (self.out/'latest.json').write_text(json.dumps(row,indent=2,default=str)+'\n')
        if status=='completed_case': (self.out/'latest_completed_case.json').write_text(json.dumps(row,indent=2,default=str)+'\n')

def solve_at_price(P,grid,price):
    """Explicit partial equilibrium at a supplied price, with canonical credit."""
    from . import credit
    price=float(price)
    if not np.isfinite(price) or price<=0: raise ValueError('positive finite price required')
    Q=copy.deepcopy(P);validate_inputs(Q,grid)
    credit.bind_engine_credit(Q,'corrected',float(Q.unsecured_credit_limit))
    sd=solver.precompute_shared(Q,np.asarray(grid))
    sol=solver.solve_markov_income_at_prices(np.array([price]),Q,np.asarray(grid),SD=sd,verbose=False,fast_stats=False)
    return {'solution':sol,'P':Q,'b_grid':np.asarray(grid).copy(),'price':price,'shared':sd}

def solve_stationary_ge(P,grid,*,out,price_start,budget_seconds=1800,max_lifecycle=32,closure='fixed_h0'):
    """Price clears birth renewal; both housing scales share this exact solve."""
    if closure not in ('fixed_h0','population_one'): raise ValueError('unknown closure')
    if not 2<=int(max_lifecycle)<=32: raise ValueError('native lifecycle limit must lie in [2,32]')
    if not np.isfinite(budget_seconds) or budget_seconds<=0: raise ValueError('positive time budget required')
    if not np.isfinite(price_start) or float(price_start)<=0: raise ValueError('positive price_start required')
    validate_inputs(P,grid)
    out=Path(out)
    if out.exists() and any(out.iterdir()): raise RuntimeError('refusing nonempty GE output directory')
    out.mkdir(parents=True,exist_ok=True)
    deadline=time.time()+float(budget_seconds)
    context=build_context(P,grid,out,price_start=price_start,deadline=deadline,max_lifecycle=max_lifecycle,closure=closure)
    budget=Budget(out,deadline,max_lifecycle)
    prior_observe=ge.observe_price;selected={}
    def capture(ctx,live,label,*,final=False):
        observed=prior_observe(ctx,live,label,final=final)
        if final and label=='selected_root': selected['live']=live
        return observed
    ge.observe_price=capture
    try: result=ge.run_phase_b(context,{'selected_d_bar':float(P.unsecured_credit_limit)},budget)
    finally: ge.observe_price=prior_observe
    if result.get('status')!='passed' or 'live' not in selected: raise RuntimeError('native GE acceptance failed: '+str(result.get('status')))
    live=selected['live'];cert=dict(result['selected']);price=float(result['selected_price'])
    factor=float((live['P'].user_cost_rate*price/live['P'].r_bar[0])**live['P'].xi_supply[0])
    implied=float(cert['normalized_housing_demand'])/factor
    cert.update(closure_mode=closure,implied_H0_at_population_one=implied,
        fixed_h0_population_scale=float(P.H0[0])/implied)
    report=out/'phase_b_ge/selected_root'
    # Native repeat verifies every numerical solution/shared array and both tables.
    # The native runner also requires byte-identical standard PNGs.
    other=out/'phase_b_ge/selected_repeat_final'
    import hashlib
    one={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (report/'standard_diagnostics').glob('*.png')}
    two={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (other/'standard_diagnostics').glob('*.png')}
    if len(one)!=17 or one!=two: raise RuntimeError('native selected repeat plot hashes differ')
    metadata=dict(context['reference_metadata'],reporting=historical_reporting_dependency(),renewal_root='adjusted_births/(2.1*entry)-1',target_fit_in_residual=False)
    return dict(solution=live['sol'],P=live['P'],b_grid=live['b_grid'].copy(),price=price,closure=cert,
        report_directory=str(report),reference_metadata=metadata,comparison_metadata={'exact_repeat':True,'standard_plot_hashes':one},
        lifecycle_solves=budget.used_lifecycle,shared=live['sd'])
