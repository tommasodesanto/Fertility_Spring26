"""Authenticated current-floor native bindings; no model solve at import/constructor.

Historical frozen objects supply observers only. The selected native state must be
rebuilt at its actual price and reproduced before a dated mapping is permitted.
"""
from __future__ import annotations
import ast, contextlib, copy, csv, gzip, hashlib, importlib.util, inspect, json, pickle, sys, threading, time, types
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

def normalized_housing_contract(handoff):
    marker=handoff.get('normalized_housing_contract')
    if marker is None:return False
    expected=dict(schema='normalized_n0_fixed_h0_v1',N0=1.,H0_source='authenticated_parameter_table',counterfactual_H0_fixed=True)
    require(marker==expected and type(marker.get('N0')) in (int,float) and marker.get('counterfactual_H0_fixed') is True,
        'Explicit normalized fixed-H0 contract differs')
    return True

def validate_normalized_housing_report(handoff,closure,actual):
    if not normalized_housing_contract(handoff):
        require('H0_derived' not in closure,'Normalized H0 report requires explicit normalized contract')
        return
    require(closure.get('normalized_population')==1. and closure.get('population_scale')==1.,
        'Normalized reference population must be N0=1')
    H0=actual.get('H0')
    require(type(H0) in (int,float) and np.isfinite(H0) and .2<=H0<=80. and closure.get('H0_bounds')==[.2,80.] and H0==closure.get('H0_derived'),
        'Authenticated parameter H0 differs from normalized closure')

def bind_normalized_reference_housing(P,handoff,closure,actual):
    if not normalized_housing_contract(handoff):return P
    validate_normalized_housing_report(handoff,closure,actual)
    require(float(P.N_target)==1.,'Native normalized reference N_target must be 1')
    require(np.asarray(P.H0).shape==(1,),'Normalized housing requires the original one-market geometry')
    P.H0=np.array([float(actual['H0'])],dtype=float)
    return P

def compare_normalized_reference_reports(reference,current):
    """Numerical calibration equivalence with declared fixed-H0 reporting only."""
    def rows(path):
        with path.open(newline='') as stream:return list(csv.DictReader(stream))
    def equal(a,b,label):
        if isinstance(a,dict) and isinstance(b,dict):
            require(set(a)==set(b),'Normalized reference keys differ: '+label)
            for key in a:equal(a[key],b[key],label+'.'+key)
        elif isinstance(a,list) and isinstance(b,list):
            require(len(a)==len(b),'Normalized reference shape differs: '+label)
            for i,(x,y) in enumerate(zip(a,b)):equal(x,y,label+'.'+str(i))
        else:
            try:x,y=float(a),float(b)
            except (ValueError,TypeError):require(a==b,'Normalized reference metadata differs: '+label);return
            require(np.isfinite(x) and np.isfinite(y) and abs(x-y)<=1e-10,'Normalized reference numeric value differs: '+label)
    current=Path(current);reference=Path(reference)
    for name,count in (('target_fit.csv',14),('parameters.csv',31)):
        a,b=rows(reference/name),rows(current/name);require(len(a)==len(b)==count,'Normalized reference full row count differs')
        for x,y in zip(a,b):
            if name=='parameters.csv' and x.get('parameter')=='H0':
                require(y.get('parameter')=='H0' and x['status']=='derived calibrated housing supply coefficient at N0=1' and
                    y['status']=='fixed population scale; not identified by per-household targets',
                    'Normalized calibration H0 role differs')
                y.update(status='fixed calibrated housing supply coefficient at N0=1',lower=x['lower'],upper=x['upper'],near_bound=x['near_bound'])
                equal({k:v for k,v in x.items() if k!='status'},{k:v for k,v in y.items() if k!='status'},'H0')
            else:equal(x,y,name)
        if name=='parameters.csv':
            with (current/name).open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(b[0]));writer.writeheader();writer.writerows(b)
    a,b=read(reference/'closure.json'),read(current/'closure.json')
    require(a.get('standard_plot_supply_units')=='physical housing supply at normalized N0=1' and
        b.get('standard_plot_supply_units')=='physical supply divided by endogenous household population',
        'Normalized/fixed-H0 report supply units differ')
    declared={'normalized_population','H0_derived','H0_bounds','housing_supply_coefficient_role','standard_plot_supply_units'}
    require(a.get('housing_supply_coefficient_role')=='derived calibrated coefficient at N0=1','Normalized calibration housing role differs')
    equal({k:v for k,v in a.items() if k not in declared},{k:v for k,v in b.items() if k not in declared},'closure')
    aa={p.name:sha(p) for p in (reference/'standard_diagnostics').glob('*.png')}
    bb={p.name:sha(p) for p in (current/'standard_diagnostics').glob('*.png')}
    require(len(aa)==17 and aa==bb,'Normalized reference standard plots differ')
    receipt=dict(status='normalized_calibration_fixed_H0_reference_equivalence_passed',target_rows=14,parameter_rows=31,
        numeric_tolerance=1e-10,standard_plot_hashes=aa,H0_metadata_exception='Derived calibration coefficient becomes fixed calibrated transition coefficient',
        native_calls=0,source_reference=str(reference),fixed_H0_report=str(current))
    write(current/'normalized_reference_compatibility.json',receipt);return receipt

def authenticate_handoff(pin):
    require(isinstance(pin,dict) and set(pin)=={'path','sha256'},'Pinned handoff required')
    path=Path(pin['path']);require(sha(path)==pin['sha256'],'Handoff hash differs')
    handoff=read(path);selected=handoff['selected_point'];normalized=normalized_housing_contract(handoff)
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
    for directory,closure in ((report,root_closure),(repeated,repeat_closure)):
        with (directory/'parameters.csv').open(newline='') as stream:
            parameter_rows=list(csv.DictReader(stream));actual={row['parameter']:float(row['estimate']) for row in parameter_rows}
        if normalized:
            H0row=next(row for row in parameter_rows if row['parameter']=='H0')
            require(float(H0row['lower'])==.2 and float(H0row['upper'])==80.,'Normalized H0 external table bounds differ')
        validate_normalized_housing_report(handoff,closure,actual)
    if modern:
        closure_proof=proof['closure_receipts']
        for key in ('root_closure','repeat_closure','native_selected_repeat_receipt','selected_verification'):
            pin=closure_proof[key];require(sha(resolved(pin['path']))==pin['sha256'],'Native closure/verification pin differs: '+key)
        native_repeat=read(resolved(closure_proof['native_selected_repeat_receipt']['path']))
        require(native_repeat['status']=='exact_full_ge_repeat_passed' and native_repeat['standard_plot_hashes']==plots,'Native repeat receipt differs')
        root_flags={'exploration_unverified':True,'selected_price_repeat_performed':False,'standard_plots_deferred_to_final_selection':True}
        repeat_flags={'standard_plot_count':17,'standard_plot_supply_units':
            'physical housing supply at normalized N0=1' if normalized else 'physical supply divided by endogenous household population'}
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
        if normalized_housing_contract(self.handoff):
            P=bind_normalized_reference_housing(P,self.handoff,read(self.report/'closure.json'),actual)
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
        self.install_observer_adapters(self.folder/'observer_adapters')
        with self.native_bindings():
            got=ctx['fp'].actual_parameters(ctx['prepared'],P,grid)
            self.ge.validate_parameter_estimates(dict(expected_parameters=actual),rows,got)
            cohort=self.pf.calendar.entrant_cohort(np.array([1.]),P,grid)
            np.testing.assert_allclose(cohort.sum(axis=(1,2,4,5)),P.fixed_reference_entry_conditional*P.z_weights[None,:],rtol=0,atol=2e-16)
        write(self.folder/'constructor.json',dict(status='native_reconstruction_pending',policy_calls=0,entry=entry,identity=self.identity()))
        return self

    def install_observer_adapters(self,folder):
        """Inject genuine native dependencies; retain every observer computation."""
        from small_credit_lab.engine import parameters,utils
        from intergen_eqscale_seq_optimized import solver as legacy,utils as legacy_utils
        folder=Path(folder);folder.mkdir(parents=True,exist_ok=True);proof={}
        def normalized(fn):return ast.dump(ast.parse(inspect.getsource(fn)),include_attributes=False)
        for name in ('readiness_childless_states','readiness_gate_active','readiness_settled_state','get_fecundity_by_age','realize_current_cross_section','add_aggregate_wealth_bequest_flow_moments','annual_gross_income_at_state'):
            native=getattr(self.model,name,None) or getattr(parameters,name)
            require(normalized(native)==normalized(getattr(legacy,name)),'Observer helper semantics differ: '+name)
            proof[name]=dict(native_source=inspect.getsourcefile(native),native_sha256=sha(inspect.getsourcefile(native)),legacy_source=inspect.getsourcefile(getattr(legacy,name)),legacy_sha256=sha(inspect.getsourcefile(getattr(legacy,name))),normalized_ast_sha256=hashlib.sha256(normalized(native).encode()).hexdigest())
        require(normalized(utils.weighted_quantile)==normalized(legacy_utils.weighted_quantile),'Wealth quantile semantics differ')
        proof['weighted_quantile']=dict(native_source=utils.__file__,native_sha256=sha(utils.__file__),legacy_source=legacy_utils.__file__,legacy_sha256=sha(legacy_utils.__file__),normalized_ast_sha256=hashlib.sha256(normalized(utils.weighted_quantile).encode()).hexdigest())
        require(self.model.DEAD_MASS_TOL==legacy.DEAD_MASS_TOL and self.model.DEAD_VALUE_CUTOFF==legacy.DEAD_VALUE_CUTOFF,'Observer dead-tail constants differ')
        housing=self.rt['observe_initial_housing_wealth'].__globals__['__file__']
        recent_wrapper=self.rt['observe_recent_parent_flow']
        closure=inspect.getclosurevars(recent_wrapper)
        recent_modules=[v for v in closure.nonlocals.values() if isinstance(v,types.ModuleType) and hasattr(v,'observe_recent_parent_flow')]
        require(len(recent_modules)==1,'Original retained-tail observer closure ambiguous')
        recent=recent_modules[0].__file__
        modules={};receipts={}
        for label,path in (('housing',housing),('recent_parent',recent)):
            original=Path(path).read_text();tree=ast.parse(original)
            imports=[node for node in ast.walk(tree) if isinstance(node,ast.ImportFrom) and node.module in ('intergen_eqscale_seq_optimized','intergen_eqscale_seq_optimized.utils')]
            require(len(imports)==(2 if label=='housing' else 1),'Observer dependency imports changed')
            lines=original.splitlines(keepends=True);substitutions=[]
            for node in imports:
                require(node.lineno==node.end_lineno and len(node.names)==1,'Observer import is not the reviewed one-line dependency')
                symbol=node.names[0]
                if node.module=='intergen_eqscale_seq_optimized':
                    require(symbol.name=='solver' and symbol.asname=='model','Observer model import differs')
                    replacement='model = _CURRENT_FLOOR_NATIVE_MODEL'
                else:
                    require(symbol.name=='weighted_quantile' and symbol.asname is None,'Observer utility import differs')
                    replacement='weighted_quantile = _CURRENT_FLOOR_NATIVE_UTILS.weighted_quantile'
                line=lines[node.lineno-1];indent=line[:len(line)-len(line.lstrip())]
                lines[node.lineno-1]=indent+replacement+'\n';substitutions.append(dict(line=node.lineno,original=line.rstrip(),replacement=lines[node.lineno-1].rstrip()))
            transformed=''.join(lines);new_tree=ast.parse(transformed)
            # Prove normalized observer AST differs ONLY by the approved imports.
            replacements={n.lineno:n for n in imports}
            class RestoreImports(ast.NodeTransformer):
                def visit_Assign(self,node):return copy.deepcopy(replacements[node.lineno]) if node.lineno in replacements else self.generic_visit(node)
            restored=RestoreImports().visit(new_tree)
            require(ast.dump(restored,include_attributes=False)==ast.dump(tree,include_attributes=False),'Observer computation changed beyond dependency injection')
            target=folder/(label+'.py');target.write_text(transformed)
            module=types.ModuleType('current_floor_observer_'+label);module.__file__=str(target)
            module._CURRENT_FLOOR_NATIVE_MODEL=self.model;module._CURRENT_FLOOR_NATIVE_UTILS=utils
            exec(compile(transformed,str(target),'exec'),module.__dict__)
            modules[label]=module;receipts[label]=dict(original=dict(path=str(path),sha256=sha(path)),adapter=dict(path=str(target),sha256=sha(target)),substitutions=substitutions,computation_ast_unchanged=True)
        self.rt['observe_initial_housing_wealth']=modules['housing'].observe_initial_housing_wealth
        def observe_recent_with_retained_tail(*args,**kwargs):
            kwargs['diagnostic_allow_retained_dead_tail']=True
            return modules['recent_parent'].observe_recent_parent_flow(*args,**kwargs)
        self.rt['observe_recent_parent_flow']=observe_recent_with_retained_tail
        write(folder/'transformation.json',dict(schema='current_floor_observer_dependency_injection_v1',observers=receipts,helpers=proof,constants=dict(DEAD_MASS_TOL=self.model.DEAD_MASS_TOL,DEAD_VALUE_CUTOFF=self.model.DEAD_VALUE_CUTOFF),model_identity=self.model.__name__,no_model_name_masquerade=True))

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
            parent_age_maturation_active=parameters.parent_age_maturation_active,
            readiness_childless_states=parameters.readiness_childless_states,readiness_gate_active=parameters.readiness_gate_active)
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
                'independent_child_maturation_active','get_fecundity_by_age','readiness_settled_state','parent_age_maturation_active',
                'readiness_childless_states','readiness_gate_active','DEAD_VALUE_CUTOFF','add_aggregate_wealth_bequest_flow_moments','annual_gross_income_at_state')
            require(float(self.model.DEAD_MASS_TOL)==1e-12,'Native dead-mass tolerance changed')
            require(all(hasattr(self.model,name) for name in required),'Missing actual native callback API: '+','.join(name for name in required if not hasattr(self.model,name)))
            for name in required:
                if name not in ('DEAD_MASS_TOL','DEAD_VALUE_CUTOFF','__file__'):
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

    def debug_unverified_solve(self,receipt_pin,folder):
        """Reconstruct/observe one authentic saved solve, with zero new solves."""
        require(set(receipt_pin)=={'path','sha256'} and sha(receipt_pin['path'])==receipt_pin['sha256'],'Saved solve receipt hash differs')
        receipt=read(receipt_pin['path']);require(receipt['schema']=='current_floor_unverified_native_solve_v1' and receipt['status']=='native_solve_completed_reconstruction_pending','Not unverified native solve evidence')
        require(receipt['identity']==self.identity(),'Saved solve numerical identity differs')
        pin=receipt['checkpoint'];require(sha(pin['path'])==pin['sha256'],'Saved native solve bytes differ')
        with gzip.open(pin['path'],'rb') as stream:saved=pickle.load(stream)
        require(saved['identity']==self.identity() and np.array_equal(saved['b_grid'],self.grid),'Saved source/grid identity differs')
        P,sd,sol=saved['parameters'],saved['shared'],saved['solution'];price=float(saved['price'])
        require(price==self.reference_price and price==float(receipt['price']),'Saved selected price differs')
        public=lambda Q:{k:v for k,v in vars(Q).items() if not k.startswith('_') and k!='native_inherited_distribution_evidence_dir'}
        require(self.scaffold.serialized(public(P))==self.scaffold.serialized(public(self.P)),'Saved selected native primitives differ')
        folder=Path(folder);folder.mkdir(parents=True,exist_ok=True);before=self.total_native_calls
        outcome=dict(schema='current_floor_unverified_observer_debug_v1',scientific_validation=False,production_ready=False,reference_verified=False,original_native_calls=int(saved['total_native_calls']),new_native_calls=0,source_receipt=receipt_pin,source_checkpoint=pin)
        def no_new_solve(*args,**kwargs):raise RuntimeError('Debug forbids any new native solve')
        try:
            with self.native_bindings():
                original_bellman=self.model.solve_bellman_full_markov_income;original_stationary=self.model.solve_markov_income_at_prices
                try:
                    self.model.solve_bellman_full_markov_income=no_new_solve;self.model.solve_markov_income_at_prices=no_new_solve
                    P._fert2_probs=sol.fert2_probs.copy();cal=self.pf.calendar
                    policy=cal.policy_from_solution(sol,np.array([price]),P,self.grid,sd)
                    pre,recon=cal.reconstruct_stationary_pre_fertility(sol,policy,P,self.grid,sd)
                    require(recon['stationary_post_fertility_nesting_l1']<=5e-9 and recon['stationary_feasibility_projection_mass']==0.,'Saved native stationary reconstruction fails')
                    supply=cal.HousingSupplyRule('static-elastic',price,float(P.H0[0]*(P.user_cost_rate*price/P.r_bar[0])**P.xi_supply[0]),float(P.xi_supply[0]))
                    ev=cal.evaluate_period(np.array([price]),pre,P,self.grid,sd,cal.SolveCounter(),supply_rule=supply,supplied_policy=policy)
                    packet=dict(parameters=P,b_grid=self.grid,shared=sd,solution=sol,evaluation=ev,stationary_g_pre=pre,supply_rule=supply,demographic_seed=self.ctx['reference'].get('demographic_seed'))
                    gates=self.ctx['fp'].gates(packet,self.ctx['prepared'],folder,stationary=True)
                    live=dict(P=P,b_grid=self.grid,sd=sd,sol=sol,price=np.array([price]))
                    self.ge.observe_price(dict(self.ctx,out=folder,deadline_epoch=float('inf')),live,'selected_debug',final=True)
                    self.runner.compare_repeated(self.report,folder/'phase_b_ge/selected_debug')
                    outcome.update(status='debug_native_reconstruction_observers_report_passed',reconstruction=recon,full_target_rows=14,full_parameter_rows=31,standard_plots=17)
                finally:
                    self.model.solve_bellman_full_markov_income=original_bellman;self.model.solve_markov_income_at_prices=original_stationary
        except Exception as exc:
            outcome.update(status='FAILED',error_type=type(exc).__name__,error=str(exc));write(folder/'debug_receipt.json',outcome);raise
        require(self.total_native_calls==before,'Debug unexpectedly consumed native calls')
        write(folder/'debug_receipt.json',outcome);return outcome

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
            if normalized_housing_contract(self.handoff):compare_normalized_reference_reports(self.report,report)
            else:self.runner.compare_repeated(self.report,report)
            results.append((packet,record))
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
    parser.add_argument('--debug-unverified-receipt',type=Path);parser.add_argument('--debug-unverified-sha256');parser.add_argument('--output',type=Path,required=True);parser.add_argument('--native-import',action='store_true');parser.add_argument('--native-api-only',action='store_true')
    args=parser.parse_args();pin=dict(path=str(args.handoff),sha256=args.handoff_sha256)
    if args.native_api_only:
        result=FloorRuntime.native_api_preflight();write(args.output/'native_api_preflight.json',result);print(json.dumps(result));sys.exit(0)
    require(args.handoff is not None and args.handoff_sha256 is not None,'Handoff path/SHA required')
    try:
        h,pins,report=authenticate_handoff(pin)
        result=dict(status='source_and_winner_evidence_passed',policy_calls=0,source_pins=len(pins),selected_report=str(report),native_setup_verified=False)
        if args.native_import or args.debug_unverified_receipt:
            runtime=FloorRuntime.from_handoff(pin,args.output/'native_import')
            result.update(status='native_import_passed_reconstruction_pending',identity=runtime.identity(),native_setup_verified=True,native_status=runtime.native_status())
            if args.debug_unverified_receipt:
                require(args.debug_unverified_sha256 is not None,'Debug receipt SHA required')
                result['debug']=runtime.debug_unverified_solve(dict(path=str(args.debug_unverified_receipt),sha256=args.debug_unverified_sha256),args.output/'saved_solve_debug')
        write(args.output/'preflight.json',result);print(json.dumps(result))
    except Exception as exc:
        write(args.output/'preflight.json',dict(status='blocked',policy_calls=0,native_setup_verified=False,exception=type(exc).__name__,reason=str(exc)))
        raise
