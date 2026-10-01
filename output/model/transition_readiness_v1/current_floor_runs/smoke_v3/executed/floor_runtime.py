"""Authenticated current-floor native bindings; no model solve at import/constructor.

Historical frozen objects supply observers only. The selected native state must be
rebuilt at its actual price and reproduced before a dated mapping is permitted.
"""
from __future__ import annotations
import contextlib, copy, csv, gzip, hashlib, importlib.util, json, pickle, sys, threading, time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
ROOT=Path(__file__).resolve().parents[4]
_LOCK=threading.RLock()

def require(ok,message):
    if not ok: raise RuntimeError(message)
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path): return json.loads(Path(path).read_text())
def write(path,value):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(value,indent=2,default=str,allow_nan=False)+'\n')
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module;spec.loader.exec_module(module);return module

def authenticate_handoff(pin):
    require(isinstance(pin,dict) and set(pin)=={'path','sha256'},'Pinned handoff required')
    path=Path(pin['path']);require(sha(path)==pin['sha256'],'Handoff hash differs')
    handoff=read(path);selected=handoff['selected_point']
    modern=selected.get('selected_repeat_status')=='exact_full_ge_repeat_passed'
    require((modern and selected.get('selected_postcheck_status')=='selected_numerically_verified') or selected.get('selected_repeat_status')=='selected_numerically_verified','Selected native repeat missing')
    provenance=handoff['checkpoint_and_sources']
    if 'source_pin_manifest_local' in provenance:
        sources=provenance['source_pin_manifest_local'];manifest=ROOT/sources['path'];digest=sources['sha256']
    else:
        manifest=ROOT/provenance['source_pins_manifest'];digest=provenance['source_pins_manifest_sha256']
    require(sha(manifest)==digest,'Source inventory hash differs')
    pins=read(manifest)
    for rel,digest in pins.items(): require(sha(ROOT/rel)==digest,'Executed floor source differs: '+rel)
    packet=ROOT/selected['local_packet']
    def resolved(value):
        path=Path(value);return path if path.is_absolute() else ROOT/path
    def report_path(key,names):
        if key in selected:return resolved(selected[key])
        found=[packet/name for name in names if (packet/name/'target_fit.csv').is_file()]
        require(len(found)==1,'Ambiguous or absent selected report: '+key)
        return found[0]
    if 'runtime_inputs' in handoff:
        report=resolved(handoff['runtime_inputs']['root_report_path']);repeated=resolved(handoff['runtime_inputs']['repeat_report_path'])
    else:
        report=report_path('root_local_path',('ROOT','root'));repeated=report_path('repeat_local_path',('REPEAT','repeat'))
    proof=handoff['selected_repeat_verification']
    plots=proof['standard_plot_sha256'] if modern else proof['plots']['sha256_by_filename']
    require(len(plots)==17,'Complete winner graph pins required')
    for directory in (report,repeated):
        for key,name,count in [('target_fit_csv','target_fit.csv',14),('parameters_csv','parameters.csv',31)]:
            require(sha(directory/name)==handoff['tables'][key]['sha256'],'Winner/repeat table differs: '+name)
            with (directory/name).open(newline='') as stream:require(len(list(csv.DictReader(stream)))==count,'Winner/repeat row count differs: '+name)
    for directory in ((repeated,) if modern else (report,repeated)):
        require({p.name for p in (directory/'standard_diagnostics').glob('*.png')}==set(plots),'Winner/repeat plot set differs')
        for name,digest in plots.items():require(sha(directory/'standard_diagnostics'/name)==digest,'Winner/repeat plot differs: '+name)
    root_closure=read(report/'closure.json');repeat_closure=read(repeated/'closure.json')
    if modern:
        closure_proof=proof['closure_receipts']
        for key in ('root_closure','repeat_closure','native_selected_repeat_receipt','selected_verification'):
            pin=closure_proof[key];require(sha(resolved(pin['path']))==pin['sha256'],'Native closure/verification pin differs: '+key)
        native_repeat=read(resolved(closure_proof['native_selected_repeat_receipt']['path']))
        require(native_repeat['status']=='exact_full_ge_repeat_passed' and native_repeat['standard_plot_hashes']==plots,'Native repeat receipt differs')
        root_flags={'exploration_unverified':True,'selected_price_repeat_performed':False,'standard_plots_deferred_to_final_selection':True}
        repeat_flags={'standard_plot_count':17,'standard_plot_supply_units':'physical supply divided by endogenous household population'}
        require(all(root_closure.get(k)==v for k,v in root_flags.items()) and all(repeat_closure.get(k)==v for k,v in repeat_flags.items()),'Selected search/report metadata differs')
        require({k:v for k,v in root_closure.items() if k not in root_flags}=={k:v for k,v in repeat_closure.items() if k not in repeat_flags},'Winner/repeat native economic closure differs')
        # Displayed handoff price may be rounded; the authenticated closure is exact.
        require(abs(float(repeat_closure['price'])-float(selected['price']))<=5e-10,'Winner closure price differs')
        report=repeated
    else:
        require(root_closure==repeat_closure,'Winner/repeat native closure differs')
        require(float(root_closure['price'])==float(selected['price']),'Winner closure price differs')
    return handoff,pins,report

class FloorRuntime:
    housing='static-elastic'
    exact_policy_cache_bytes=64*1024**3
    @classmethod
    def from_handoff(cls,pin,folder):
        self=cls();self.folder=Path(folder);self.folder.mkdir(parents=True,exist_ok=True)
        self.handoff,self.sourcepins,self.report=authenticate_handoff(pin);self.handoff_pin=dict(pin)
        old=ROOT/'output/model/fixed_reference_economics_20260928/utility_floor_round2_v1'
        sys.path[:0]=[str(old),str(ROOT/'code/model/tools'),str(ROOT/'code/model')]
        inputs=load('inputs',old/'inputs.py');self.runner=load('floor_native_runner',old/'runner.py')
        sys.path.insert(0,str(self.runner.BASE))
        base=load('floor_matched_workflow',self.runner.BASE/'run_comparison.py')
        self.ge=load('floor_native_observer',old/'phase_b_pilot.py')
        self.runner.install_reporter_on_authored(base.authored)
        ctx=base.authored.context_from_bundle(SimpleNamespace(bundle=ROOT/'output/model/publication_refactor_20260929/local_export_v1/inputs',reference_root=ROOT,out=self.folder))
        base.authored.authenticate_frozen(ctx)
        P,grid=inputs.proposal('floor_s0');P,entry=inputs.entry(P,grid,'nonnegative_mean')
        rows=self.runner.readtable(self.report/'parameters.csv');actual={r['parameter']:float(r['estimate']) for r in rows}
        bounds=self.handoff['contract']['bounds'];point={k:actual[k] for k in self.handoff['contract']['free_coordinates']}
        P=inputs.bind(P,point,bounds,'floor')
        from small_credit_lab import credit
        from small_credit_lab.engine import solver
        credit.bind_engine_credit(P,'corrected',0.)
        ctx.update(P=P,b_grid=grid,selected_d_bar=0.,reference_psi=float(P.psi_child),expected_parameters=actual,free_coordinates=list(point),out=self.folder,deadline_epoch=float('inf'))
        self.runner.install_observer_metadata(self.ge,'floor',bounds)
        for row in ctx['manifest']['full_parameter_table']:
            if row['parameter'] in bounds: row['lower'],row['upper']=map(str,bounds[row['parameter']])
        self.parameter_rows=rows
        self.ctx=ctx;self.P=P;self.grid=grid;self.model=solver;self.rt=ctx['prepared'].rt;self.pf=self.rt['primitive'].pf
        self.reference_price=float(self.handoff['selected_point']['selected_price']) if 'selected_price' in self.handoff['selected_point'] else float(self.handoff['selected_point']['price'])
        closure=read(self.report/'closure.json');self.reference_price=float(closure['price']);self.population_scale=float(closure['population_scale'])
        self.packet=None;self.initial_state=None;self.reference_verified=False;self.total_native_calls=0
        self.scaffold=load('floor_dated_scaffold',Path(__file__).parent/'pinned_tools/run_e5f_preference_transition.py')
        with self.native_bindings():
            got=ctx['fp'].actual_parameters(ctx['prepared'],P,grid)
            self.ge.validate_parameter_estimates(dict(expected_parameters=actual),rows,got)
            cohort=self.pf.calendar.entrant_cohort(np.array([1.]),P,grid)
            np.testing.assert_allclose(cohort.sum(axis=(1,2,4,5)),P.fixed_reference_entry_conditional*P.z_weights[None,:],rtol=0,atol=2e-16)
        write(self.folder/'constructor.json',dict(status='native_reconstruction_pending',policy_calls=0,entry=entry,identity=self.identity()))
        return self

    def identity(self):
        canonical=lambda value:hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        engine={k:v for k,v in self.sourcepins.items() if '/small_credit_lab/engine/' in k}
        return dict(reference_sha256=self.handoff_pin['sha256'],engine_sha256=canonical(engine),
                    entry_sha256=hashlib.sha256(self.P.fixed_reference_entry_conditional.tobytes()).hexdigest(),
                    grid_sha256=hashlib.sha256(self.grid.tobytes()).hexdigest(),
                    effective_parameters_sha256=canonical({r['parameter']:r['estimate'] for r in self.parameter_rows}),
                    source_pins=self.sourcepins)

    def native_status(self):
        return dict(reference_verified=self.reference_verified,exact_policy_cache_bytes=self.exact_policy_cache_bytes,reference_price=self.reference_price,
                    population_scale=self.population_scale,housing=self.housing,grid=[int(self.P.Nb),int(self.P.Nz)],
                    reconstruction_pending=not self.reference_verified)

    @contextlib.contextmanager
    def native_budget(self,deadline,remaining_calls):
        """Enforce a monotonic deadline and actual native-call allowance before entry.

        Nested guards all apply. An old cache swallowing an exception cannot
        start a fallback solve once the time/call allowance is exhausted.
        """
        require(isinstance(deadline,(int,float)) and np.isfinite(deadline),'Finite monotonic deadline required')
        require(type(remaining_calls) is int and remaining_calls>=0,'Nonnegative actual-call allowance required')
        if not hasattr(self,'_native_budget_stack'):self._native_budget_stack=[]
        budget=dict(deadline=float(deadline),remaining=remaining_calls,used=0)
        self._native_budget_stack.append(budget)
        try:yield budget
        finally:
            require(self._native_budget_stack[-1] is budget,'Native budget contexts exited out of order')
            self._native_budget_stack.pop()

    def _guard_native_call(self):
        budgets=getattr(self,'_native_budget_stack',())
        now=time.monotonic()
        for budget in budgets:
            if now>=budget['deadline']:raise TimeoutError('Native model-call deadline exhausted before solve')
            if budget['remaining']<=0:raise RuntimeError('Native model-call allowance exhausted before solve')
        for budget in budgets:
            budget['remaining']-=1;budget['used']+=1
        self.total_native_calls=getattr(self,'total_native_calls',0)+1

    @contextlib.contextmanager
    def native_api_bindings(self):
        """Expose only the genuine native utility omitted by the solver facade."""
        from small_credit_lab.engine import utils,parameters
        aliases=dict(interp_indices=utils.interp_indices,
            independent_child_maturation_active=parameters.independent_child_maturation_active,
            get_fecundity_by_age=parameters.get_fecundity_by_age,
            readiness_settled_state=parameters.readiness_settled_state,
            parent_age_maturation_active=parameters.parent_age_maturation_active)
        saved={name:(hasattr(self.model,name),getattr(self.model,name,None)) for name in aliases}
        require(all(not saved[name][0] or saved[name][1] is value for name,value in aliases.items()),'Conflicting native compatibility helper')
        try:
            for name,value in aliases.items():setattr(self.model,name,value)
            required=('__file__','interp_indices','birth_destination_child_state','precompute_shared','solve_bellman_full_markov_income',
                'solve_markov_income_at_prices','forward_distribution_markov_income',
                'advance_cohort_one_period_markov_income','property_tax_revenue_from_distribution',
                'build_forward_tenure_transition_maps','income_transition_values','entry_wealth_grid_weights',
                'DEAD_MASS_TOL','_censor_entry_dead_mass','_gate_dead_mass_at_age',
                'realize_current_cross_section','realize_stayer_cross_section','age_to_index',
                'apply_child_aging','bequest_utility_vec','compute_markov_statistics','get_completed_fertility',
                'income_at_state','owner_borrowing_floor','renter_borrowing_floor',
                'independent_child_maturation_active','get_fecundity_by_age','readiness_settled_state','parent_age_maturation_active')
            require(float(self.model.DEAD_MASS_TOL)==1e-12,'Native dead-mass tolerance changed')
            require(all(hasattr(self.model,name) for name in required),'Missing actual native callback API: '+','.join(name for name in required if not hasattr(self.model,name)))
            for name in required:
                if name not in ('DEAD_MASS_TOL','__file__'):
                    value=getattr(self.model,name)
                    require(callable(value) and getattr(value,'__module__','').startswith('small_credit_lab.engine.'),'Callback is not genuine native engine: '+name)
            yield {name:getattr(getattr(self.model,name),'__module__',self.model.__name__) for name in required}
        finally:
            for name,(existed,old) in saved.items():
                if existed:setattr(self.model,name,old)
                else:delattr(self.model,name)

    @classmethod
    def native_api_preflight(cls):
        source=ROOT/'output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed/source'
        sys.path.insert(0,str(source))
        from small_credit_lab.engine import solver,household,distribution,utils,parameters
        require(Path(solver.__file__).resolve().is_relative_to(source.resolve()) and Path(utils.__file__).resolve().is_relative_to(source.resolve()),'Native module source root differs')
        obj=cls();obj.model=solver
        require(solver.solve_bellman_full_markov_income is household.solve_bellman_full_markov_income,'Wrong native backward callee')
        require(solver.forward_distribution_markov_income is distribution.forward_distribution_markov_income,'Wrong native forward callee')
        with _LOCK,obj.native_api_bindings() as api:
            require(solver.interp_indices is utils.interp_indices,'Interpolation must be actual native utility')
            require(solver.independent_child_maturation_active is parameters.independent_child_maturation_active,'Maturation must be actual native helper')
            require(solver.get_fecundity_by_age is parameters.get_fecundity_by_age and solver.readiness_settled_state is parameters.readiness_settled_state and solver.parent_age_maturation_active is parameters.parent_age_maturation_active,'Indirect population helper identity differs')
            require(solver.birth_destination_child_state is household.birth_destination_child_state,'Birth destination helper differs')
            return dict(status='native_callback_api_passed',policy_calls=0,native_module=str(Path(solver.__file__).resolve()),native_utils=str(Path(utils.__file__).resolve()),native_parameters=str(Path(parameters.__file__).resolve()),source_hashes={str(Path(m.__file__).resolve()):sha(m.__file__) for m in (solver,utils,parameters)},callbacks=api)

    @contextlib.contextmanager
    def native_bindings(self):
        """Rebind BOTH actual callee chains, retaining unchanged calendar accounting."""
        from small_credit_lab.engine import household,distribution,kernels
        require(self.model.solve_bellman_full_markov_income is household.solve_bellman_full_markov_income,'Wrong backward callee')
        require(self.model.forward_distribution_markov_income is distribution.forward_distribution_markov_income,'Wrong stationary forward callee')
        require(self.model.advance_cohort_one_period_markov_income is distribution.advance_cohort_one_period_markov_income,'Wrong dated forward callee')
        for module in (household,distribution):
            for name,value in vars(module).items():
                if callable(value) and getattr(value,'__module__','').endswith('.kernels'):
                    require(value is getattr(kernels,name),'Wrong kernel callee: '+name)
        require(self.pf.transition.calendar is self.pf.calendar,'Population transition uses a different calendar module')
        targets=[(self.pf.calendar,'model',self.model),(self.rt['primitive'],'model',self.model),(self.rt['audit'],'model',self.model)]
        for name in ('apply_sequential_fertility','advance_sequential_calendar_distribution'):
            require(getattr(self.pf.transition,name).__globals__['calendar'] is self.pf.calendar,'Indirect population callback resolves a different calendar')
        configure=self.pf.transition.configure_sequential_model
        def forbidden_reconfigure(*args,**kwargs):
            raise RuntimeError('Legacy sequential reconfiguration forbidden inside native callbacks')
        with _LOCK,self.native_api_bindings():
            saved=[(obj,key,getattr(obj,key)) for obj,key,_ in targets];old_rt=self.rt['model']
            try:
                for obj,key,value in targets:setattr(obj,key,value)
                self.rt['model']=self.model;self.pf.transition.configure_sequential_model=forbidden_reconfigure
                require(self.pf.transition.calendar.model is self.model and self.pf.calendar.model is self.model,'Native forward callback chain differs')
                yield
            finally:
                changed=self.pf.transition.calendar.model is not self.model
                self.rt['model']=old_rt;self.pf.transition.configure_sequential_model=configure
                for obj,key,value in saved:setattr(obj,key,value)
                require(not changed,'Actual population callback replaced the genuine native engine')

    def stationary(self,psi,price,folder):
        """One counted native Bellman+KFE solve and unchanged stationary gates."""
        folder=Path(folder);folder.mkdir(parents=True,exist_ok=True);P=copy.deepcopy(self.P);P.psi_child=float(psi)
        ctx=dict(self.ctx,P=P,reference_psi=float(psi),out=folder)
        counter=self.pf.calendar.SolveCounter()
        with self.native_bindings():
            sd=self.model.precompute_shared(P,self.grid)
            native_solve=self.model.solve_markov_income_at_prices
            original_bellman=native_solve.__globals__['solve_bellman_full_markov_income'];calls=0
            def counted_stationary(*args,**kwargs):
                nonlocal calls
                self._guard_native_call()
                calls+=1;return original_bellman(*args,**kwargs)
            try:
                native_solve.__globals__['solve_bellman_full_markov_income']=counted_stationary
                sol=native_solve(np.array([price]),P,self.grid,SD=sd,verbose=False,fast_stats=False)
            finally:native_solve.__globals__['solve_bellman_full_markov_income']=original_bellman
            self.save_unverified_native_solve(P,self.grid,sd,sol,float(price),folder,calls)
            P._fert2_probs=sol.fert2_probs.copy();cal=self.pf.calendar
            policy=cal.policy_from_solution(sol,np.array([price]),P,self.grid,sd)
            pre,recon=cal.reconstruct_stationary_pre_fertility(sol,policy,P,self.grid,sd)
            require(recon['stationary_post_fertility_nesting_l1']<=5e-9 and recon['stationary_feasibility_projection_mass']==0.,'Native stationary reconstruction fails')
            supply=cal.HousingSupplyRule('static-elastic',float(price),float(P.H0[0]*(P.user_cost_rate*price/P.r_bar[0])**P.xi_supply[0]),float(P.xi_supply[0]))
            ev=cal.evaluate_period(np.array([price]),pre,P,self.grid,sd,counter,supply_rule=supply,supplied_policy=policy)
            packet=dict(parameters=P,b_grid=self.grid,shared=sd,solution=sol,evaluation=ev,stationary_g_pre=pre,supply_rule=supply,demographic_seed=self.ctx['reference'].get('demographic_seed'))
            gates=self.ctx['fp'].gates(packet,self.ctx['prepared'],folder,stationary=True)
            births=self.pf.transition.calendar_topcode_birth_accounting(ev.g_pre,ev.g_post_fertility,float(ev.births),P)['topcode_adjusted_birth_children']
        entry=float(sol.entry_rate);demand=float(np.sum(ev.demand_by_loc));supplyq=float(np.sum(ev.supply_by_loc));scale=supplyq/demand
        record=dict(price=float(price),psi_child=float(psi),pension=float(P.pension),population_scale=scale,renewal_residual=float(births/(2.1*entry)-1),gates=gates,accounting_valid=True,policy_calls=calls,total_native_calls=self.total_native_calls,lifecycle_calls=1,absolute_housing_demand=scale*demand,absolute_housing_supply=supplyq)
        write(folder/'stationary.json',record);return packet,record

    def save_unverified_native_solve(self,P,grid,shared,solution,price,folder,calls):
        """Save genuine solve evidence before observers; never a verified reference."""
        folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
        path=folder/'native_solve_unverified.pkl.gz';temporary=folder/'native_solve_unverified.pkl.gz.tmp'
        packet=dict(parameters=P,b_grid=grid,shared=shared,solution=solution,price=float(price),
            identity=self.identity(),total_native_calls=self.total_native_calls,policy_calls=calls,
            stage='native_solve_completed_reconstruction_pending',reference_verified=False)
        with gzip.open(temporary,'wb') as stream:pickle.dump(packet,stream)
        temporary.replace(path)
        receipt=dict(schema='current_floor_unverified_native_solve_v1',status='native_solve_completed_reconstruction_pending',
            checkpoint=dict(path=str(path),sha256=sha(path)),identity=self.identity(),price=float(price),
            psi_child=float(P.psi_child),policy_calls=calls,total_native_calls=self.total_native_calls,
            reference_verified=False,scientific_validation=False,production_ready=False)
        write(folder/'native_solve_unverified.json',receipt);return receipt

    def stationary_state(self,packet,population_scale=1.):
        scale=float(population_scale);require(np.isfinite(scale) and scale>0,'Invalid actual population')
        g=packet['stationary_g_pre']*scale
        return self.pf.stationary_initial_state(g,float(g[:,:,:,0].sum()),float(packet['evaluation'].births)*scale,packet['parameters'],1/2.1)

    def reconstruct_reference(self,folder):
        """Two fresh selected-price solves, full native 14/31/17 report comparison."""
        folder=Path(folder);results=[]
        for i in range(2):
            packet,record=self.stationary(self.P.psi_child,self.reference_price,folder/f'repeat_{i}')
            require(abs(record['renewal_residual'])<=1e-6,'Selected renewal differs')
            require(abs(record['population_scale']/self.population_scale-1)<=1e-10,'Selected population differs')
            live=dict(P=packet['parameters'],b_grid=self.grid,sd=packet['shared'],sol=packet['solution'],price=np.array([self.reference_price]))
            ctx=dict(self.ctx,out=folder/f'repeat_{i}',deadline_epoch=float('inf'))
            with self.native_bindings():self.ge.observe_price(ctx,live,'selected_root',final=True)
            report=folder/f'repeat_{i}/phase_b_ge/selected_root'
            self.runner.compare_repeated(self.report,report);results.append((packet,record))
        self.runner.compare_repeated(folder/'repeat_0/phase_b_ge/selected_root',folder/'repeat_1/phase_b_ge/selected_root')
        self.packet=results[-1][0];self.initial_state=self.stationary_state(self.packet,self.population_scale);self.reference_verified=True
        with gzip.open(folder/'selected_native_packet.pkl.gz','wb') as f:pickle.dump(self.packet,f)
        write(folder/'reference_reconstruction.json',dict(schema='current_floor_reference_reconstruction_v1',status='passed',policy_calls=sum(row[1]['policy_calls'] for row in results),lifecycle_calls=2,checkpoint=dict(path=str(folder/'selected_native_packet.pkl.gz'),sha256=sha(folder/'selected_native_packet.pkl.gz')),checkpoint_sha256=sha(folder/'selected_native_packet.pkl.gz'),identity=self.identity(),reports={str(report):{str(p.relative_to(report)):sha(p) for p in report.rglob('*') if p.is_file() and (p.name in ('parameters.csv','target_fit.csv','closure.json') or p.suffix=='.png')} for report in (folder/'repeat_0/phase_b_ge/selected_root',folder/'repeat_1/phase_b_ge/selected_root')}))
        return self.packet

    def load_reconstructed_reference(self,receipt_pin,folder):
        require(isinstance(receipt_pin,dict) and set(receipt_pin)=={'path','sha256'},'Pinned reconstruction receipt required')
        receipt_path=Path(receipt_pin['path'])
        require(sha(receipt_path)==receipt_pin['sha256'],'Reference reconstruction receipt hash differs')
        receipt=read(receipt_path)
        require(receipt.get('schema')=='current_floor_reference_reconstruction_v1' and receipt.get('status')=='passed','Reference reconstruction not verified')
        require(receipt['identity']==self.identity(),'Reconstructed reference numerical identity differs')
        require(receipt['lifecycle_calls']==2 and type(receipt['policy_calls']) is int and receipt['policy_calls']>=2,'Two native reference solves were not recorded')
        checkpoint=receipt['checkpoint'];require(set(checkpoint)=={'path','sha256'},'Pinned native checkpoint required')
        require(checkpoint['sha256']==receipt['checkpoint_sha256'] and sha(checkpoint['path'])==checkpoint['sha256'],'Reference checkpoint hash differs')
        reports=receipt['reports'];require(len(reports)==2,'Both reconstruction report inventories required')
        report_paths=[]
        for location,pins in reports.items():
            report=Path(location);report_paths.append(report)
            for relative,digest in pins.items():
                path=(report/relative).resolve()
                require(path.is_relative_to(report.resolve()) and sha(path)==digest,'Reconstructed report changed: '+relative)
            require(all(name in pins for name in ('target_fit.csv','parameters.csv','closure.json')),'Incomplete native reconstruction report')
            require(len([name for name in pins if name.endswith('.png')])==17,'Reconstructed17 graphs missing')
            self.runner.compare_repeated(self.report,report)
        self.runner.compare_repeated(*report_paths)
        with gzip.open(checkpoint['path'],'rb') as stream:packet=pickle.load(stream)
        require(np.array_equal(packet['b_grid'],self.grid),'Reconstructed reference grid differs')
        def encoded(value):
            if isinstance(value,np.ndarray):return ('array',value.dtype.str,value.shape,hashlib.sha256(value.tobytes()).hexdigest())
            if isinstance(value,np.generic):return encoded(value.item())
            if isinstance(value,float):return ('float',value.hex())
            if isinstance(value,dict):return {str(k):encoded(v) for k,v in value.items()}
            if isinstance(value,(list,tuple)):return (type(value).__name__,[encoded(v) for v in value])
            if isinstance(value,SimpleNamespace):return encoded(vars(value))
            return value
        def public(P):return {k:encoded(v) for k,v in vars(P).items() if not k.startswith('_') and k!='native_inherited_distribution_evidence_dir'}
        require(public(packet['parameters'])==public(self.P),'All public selected native primitives must match')
        require(float(packet['solution'].p_eq[0])==self.reference_price,'Reconstructed native price differs')
        require(np.array_equal(packet['stationary_g_pre'],packet['evaluation'].g_pre),'Reconstructed predecision state differs')
        require(abs(float(packet['stationary_g_pre'].sum())-1.)<=1e-9,'Saved checkpoint must contain unit stationary state')
        folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
        with self.native_bindings():
            actual=self.ctx['fp'].actual_parameters(self.ctx['prepared'],packet['parameters'],self.grid)
            self.ge.validate_parameter_estimates(dict(expected_parameters={r['parameter']:float(r['estimate']) for r in self.parameter_rows}),self.parameter_rows,actual)
            self.ctx['fp'].gates(packet,self.ctx['prepared'],folder,stationary=True)
            ev=packet['evaluation'];births=self.pf.transition.calendar_topcode_birth_accounting(ev.g_pre,ev.g_post_fertility,float(ev.births),packet['parameters'])['topcode_adjusted_birth_children']
        require(abs(float(births)/(2.1*float(packet['solution'].entry_rate))-1.)<=1e-6,'Restored native renewal gate fails')
        physical=float(packet['supply_rule'].quantity([self.reference_price])[0])
        expected=float(self.P.H0[0]*(self.P.user_cost_rate*self.reference_price/self.P.r_bar[0])**self.P.xi_supply[0])
        require(physical==expected and abs((physical/float(np.sum(ev.demand_by_loc)))/self.population_scale-1)<=1e-10,'Restored housing/population closure differs')
        self.packet=packet;self.initial_state=self.stationary_state(packet,self.population_scale);self.reference_verified=True
        write(folder/'reference_resume.json',dict(status='passed',identity=self.identity(),receipt=receipt_pin,checkpoint=checkpoint,policy_calls=0,lifecycle_calls=0,actual_initial_population=float(self.initial_state.g_pre.sum())))
        return packet

    def restore_reference(self,pin_checkpoint,folder):
        require(isinstance(pin_checkpoint,dict) and set(pin_checkpoint)=={'path','sha256','receipt'},'Checkpoint and reconstruction receipt pins required')
        require(sha(pin_checkpoint['path'])==pin_checkpoint['sha256'],'Pinned resume checkpoint differs')
        receipt=read(pin_checkpoint['receipt']['path'])
        require(receipt['checkpoint']=={k:pin_checkpoint[k] for k in ('path','sha256')},'Resume receipt names another checkpoint')
        return self.load_reconstructed_reference(pin_checkpoint['receipt'],folder)

    def mapping(self,terminal,endpoint,prices,pensions,psi_path,folder,initial_state=None,start_year=2007):
        require(self.reference_verified and self.packet is not None,'Fresh selected native reconstruction/repeat required')
        folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
        prices=np.asarray(prices,float);pensions=np.asarray(pensions,float);psi_path=np.asarray(psi_path,float)
        require(prices.shape==pensions.shape==psi_path.shape and prices.ndim==1 and len(prices)>0,'Dated path shape differs')
        require(all(np.isfinite(x).all() and (x>0).all() for x in (prices,pensions,psi_path)),'Invalid dated levels')
        state=copy.deepcopy(self.initial_state if initial_state is None else initial_state)
        pf=self.pf;clock=pf.transition.advance_adult_entry_clock;evaluate=pf.calendar.evaluate_period
        original_supply=self.scaffold.supply_rule;original_solve=self.model.solve_bellman_full_markov_income
        calls=0;queue_inputs=[];states={};current={}
        def counted(*args,**kwargs):
            nonlocal calls
            self._guard_native_call()
            calls+=1;return original_solve(*args,**kwargs)
        def observed_period(*args,**kwargs):
            ev=evaluate(*args,**kwargs);current.update(ev=ev,P=args[2] if len(args)>2 else kwargs['P']);return ev
        def observed_clock(queue,*args,**kwargs):
            queue_inputs.append(copy.deepcopy(queue));answer=clock(queue,*args,**kwargs)
            if len(queue_inputs)%2==0:
                t=len(queue_inputs)//2-1
                states[t]=dict(state=pf.PFInitialState(current['ev'].g_pre.copy(),queue_inputs[-2],queue_inputs[-1]),parameters=copy.deepcopy(current['P']))
            return answer
        with self.native_bindings():
            try:
                self.model.solve_bellman_full_markov_income=counted
                pf.calendar.evaluate_period=observed_period;pf.transition.advance_adult_entry_clock=observed_clock
                self.scaffold.supply_rule=lambda packet,_pf,housing:packet['supply_rule']
                native,record=self.scaffold.mapping(self.packet,self.ctx['prepared'],terminal,endpoint,prices,pensions,psi_path,'elastic_reference',folder,self.exact_policy_cache_bytes,capture=True,measure_fertility=True,initial_state=state,start_year=start_year)
            finally:
                self.model.solve_bellman_full_markov_income=original_solve
                pf.calendar.evaluate_period=evaluate;pf.transition.advance_adult_entry_clock=clock
                self.scaffold.supply_rule=original_supply
        require(len(queue_inputs)==2*len(prices) and len(states)==len(prices),'Both actual native queues must be captured each date')
        native.dated_states=states;native.floor_runtime_paths=dict(prices=prices,pensions=pensions,psi_path=psi_path,start_year=start_year)
        native.floor_runtime_record=record
        record.update(accounting_valid=all(record['gates'].values()),policy_calls=calls,total_native_calls=getattr(self,'total_native_calls',0),actual_initial_population=float(state.g_pre.sum()),housing=self.housing,exact_policy_cache_bytes=self.exact_policy_cache_bytes)
        write(folder/'mapping.json',record);return native,record

    def terminal_checks(self,packet,endpoint,native,psi_path,tolerance=1e-3,raw_queue_tolerance=1e-3):
        require(tolerance in (1e-6,1e-3) and raw_queue_tolerance in (1e-6,1e-3),'Original stationary or terminal tolerances required')
        state=self.stationary_state(packet,endpoint['population_scale']);P=packet['parameters'];entry=float(state.g_pre[:,:,:,0].sum())
        check=self.pf.terminal_convergence_diagnostics(evaluation=native,psi_path=psi_path,reference_state=state,reference_entry_flow=entry,reference_price=endpoint['price'],reference_psi=P.psi_child,base_parameters=P,tolerances={k:tolerance for k in self.pf.DEFAULT_TERMINAL_TOLERANCES})
        gap=float(np.max(np.abs(self.pf.birth_queue_values(native.terminal_state.scheduled_raw_entries)-self.pf.birth_queue_values(state.scheduled_raw_entries))))/max(entry,1e-15)
        check.update(raw_queue_maximum_relative_gap=gap,raw_queue_pass=gap<=raw_queue_tolerance)
        check['all_checks_pass']=check['all_checks_pass'] and check['raw_queue_pass'];return check

    def export_2023(self,candidate,folder):
        native=candidate['native_reply'];paths=native.floor_runtime_paths
        index=(2023-paths['start_year'])//4
        require(paths['start_year']+4*index==2023 and index in native.dated_states,'Actual 2023 clock capture missing')
        captured=native.dated_states[index];state=copy.deepcopy(captured['state'])
        packet=dict(schema='current_floor_actual_2023_v1',initial_state=state,parameters=captured['parameters'],b_grid=self.grid.copy(),current_2023_V=native.values[index].copy(),continuation_V=native.values[index+1].copy(),continuation_calendar_year=2027,terminal_V=native.values[-1].copy(),forecast_prices=paths['prices'][index:].copy(),forecast_pensions=paths['pensions'][index:].copy(),forecast_psi=paths['psi_path'][index:].copy(),calendar_year=2023,period=index,reference_identity=self.identity(),initial_population=float(state.g_pre.sum()),scheduled_entries=copy.deepcopy(state.scheduled_entries),scheduled_raw_entries=copy.deepcopy(state.scheduled_raw_entries))
        folder=Path(folder);folder.mkdir(parents=True,exist_ok=True);path=folder/'actual_2023.pkl.gz'
        with gzip.open(path,'wb') as f:pickle.dump(packet,f)
        receipt=dict(path=str(path),sha256=sha(path),calendar_year=2023,period=index,period_index=index,exact_native_state=True,reconstructed_or_rescaled=False,queue_lags=[16,20],forecast_and_continuation_saved=True,initial_population=packet['initial_population'],actual_inherited_state=True,no_rescaling=True,both_native_queues=True,reference_identity=self.identity())
        write(folder/'actual_2023.json',receipt);return receipt

    def render_standard(self,candidate,folder):
        native=candidate['native_reply'];folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
        with self.native_bindings():
            return self.scaffold.render_diagnostics(native.floor_runtime_record['diagnostic_packets'],folder,self.rt['audit'],self.ctx['manifest']['standard_diagnostic_names'])

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description='Zero-model-call current-floor authentication/import preflight')
    parser.add_argument('--handoff',type=Path);parser.add_argument('--handoff-sha256')
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--native-import',action='store_true');parser.add_argument('--native-api-only',action='store_true')
    args=parser.parse_args();pin=dict(path=str(args.handoff),sha256=args.handoff_sha256)
    if args.native_api_only:
        result=FloorRuntime.native_api_preflight();write(args.output/'native_api_preflight.json',result);print(json.dumps(result));sys.exit(0)
    require(args.handoff is not None and args.handoff_sha256 is not None,'Handoff path/SHA required')
    try:
        h,pins,report=authenticate_handoff(pin)
        result=dict(status='source_and_winner_evidence_passed',policy_calls=0,source_pins=len(pins),selected_report=str(report),native_setup_verified=False)
        if args.native_import:
            runtime=FloorRuntime.from_handoff(pin,args.output/'native_import')
            result.update(status='native_import_passed_reconstruction_pending',identity=runtime.identity(),native_setup_verified=True,native_status=runtime.native_status())
        write(args.output/'preflight.json',result);print(json.dumps(result))
    except Exception as exc:
        write(args.output/'preflight.json',dict(status='blocked',policy_calls=0,native_setup_verified=False,exception=type(exc).__name__,reason=str(exc)))
        raise
