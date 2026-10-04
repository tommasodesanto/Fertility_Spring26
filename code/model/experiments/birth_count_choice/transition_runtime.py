"""Current one-birth Estate-A bindings for the retained dated-transition controller.

Construction authenticates a safe saved case and builds observer context only.
The inherited runtime's numerical mappings, budgets, gates and queues are reused;
the old floor/small-credit constructor is never entered. No solve occurs at import.
"""
from __future__ import annotations

import contextlib
import copy
import csv
import gzip
import hashlib
import inspect
import json
import os
import platform
import pickle
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from ..transition_readiness import floor_runtime as retained
from .model import reporting, storage
from .model import native_phase_b
from .model.engine import adult_entry, distribution, household, kernels, parameters, solver, utils
from .model.observer_adapters import apply_count_fertility

ROOT = Path(__file__).resolve().parents[4]
SCHEMA = 'current_estate_a_transition_handoff_v1'
CASE = ROOT / 'output/model/experiments/birth_count_choice/estate_a_v1/single/cases/20261003T212605039706Z_a739edc3'
COUNT_FIELDS = ('birth_count_action_probs', 'birth_count_realized_probs')
BRIDGE_ABSOLUTE_TOLERANCE = 1e-10
ARRAY_COMPARISON_CRITERIA = dict(schema='current_estate_a_platform_array_criteria_v2',
    absolute_tolerance=1e-10, value_eps_multiplier=16., value_eps=float(np.finfo(np.float64).eps),
    living_control_relative_tolerance=float(np.sqrt(np.finfo(np.float64).eps)),
    occupied_mass_cutoff=1e-12, total_value_dead_cutoff=float(solver.DEAD_VALUE_CUTOFF),
    float32_probability_ulps=2., exact_float32_zero_support=True,
    dead_control_status='audited separately; no distance-equivalence claim',
    optimizer_identity_certified=False, native_numerical_gates_unchanged=True)
CONTROL_FIELDS = ('c_pol','bp_pol','c_pol_stay','bp_pol_stay')
require, sha, read, write = retained.require, retained.sha, retained.read, retained.write


def _path(value):
    value = Path(value)
    return value if value.is_absolute() else ROOT / value


def authenticate_handoff(pin):
    require(isinstance(pin, dict) and set(pin) == {'path', 'sha256'}, 'Pinned current-case handoff required')
    require(sha(pin['path']) == pin['sha256'], 'Current-case handoff hash differs')
    handoff = read(pin['path'])
    require(handoff.get('schema') == SCHEMA, 'Current Estate-A handoff schema differs')
    case = _path(handoff['saved_case']).resolve()
    pins = handoff['saved_files']
    required = {'metadata.json', 'native_result.npz', 'input_contract.json', 'target_fit.csv', 'parameters.csv',
                'native/phase_b_ge/selected_root/closure.json'}
    require(required.issubset(pins), 'Complete current-case input/report pins required')
    for relative, digest in pins.items():
        path = (case / relative).resolve()
        require(path.is_relative_to(case) and sha(path) == digest, 'Saved case file differs: ' + relative)
    sources = handoff['source_pins']
    require(sources and any('/birth_count_choice/model/engine/' in k for k in sources), 'Current engine source pins required')
    for relative, digest in sources.items():
        path = _path(relative).resolve()
        require(path.is_relative_to(ROOT) and sha(path) == digest, 'Current transition source differs: ' + relative)
    return handoff, sources, case


def validate_contract(P, closure, contract):
    from .model.estate_contract import experiment_flags
    flags = experiment_flags(1)
    require(contract.get('experiment_flags') == flags, 'Saved case is not the one-birth Estate-A contract')
    for key, value in flags.items():
        require(type(getattr(P, key, None)) is type(value) and getattr(P, key) == value,
                'Native Estate-A flag differs: ' + key)
    require(contract.get('closure') == 'fixed_h0' and closure.get('population_scale') is not None, 'Fixed baseline H0 required')
    require(np.isfinite(float(closure['population_scale'])) and float(closure['population_scale']) > 0,
            'Actual inherited population scale required')
    require(not bool(getattr(P, 'joint_nested_choice', False)), 'Joint-choice transition route is not authorized')
    require(getattr(P, 'child_state_mode', None) == 'independent_count' and
            not bool(getattr(P, 'readiness_gate_enabled', False)), 'Selected independent child-count contract differs')
    require(str(getattr(P, 'estate_receiver', 'none')) == 'none', 'Adult estate receiver is not authorized')
    require(getattr(P, 'unsecured_credit_limit', None) == 0., 'Selected renter credit capacity differs')
    require(not bool(getattr(P, 'experimental_natural_solvency', False)), 'Natural-solvency fallback is forbidden')
    require(float(P.psi) == .06 and np.asarray(P.H0).shape == (1,), 'Estate selling cost or one-market housing differs')
    # Native selected post-interest kernels map net sale/purchase into b before R.
    require(bool(getattr(P, 'native_purchase_income', False)) and bool(getattr(P, 'native_due_stayer_credit', False)),
            'Selected soft-credit purchase/stayer contract differs')


def compare_reports(reference, current):
    """Same effective primitives, complete native moments, and stable plot set."""
    reference, current = Path(reference), Path(current)
    for name, count in (('target_fit.csv', 14), ('parameters.csv', 31)):
        def rows(folder):
            with (folder / name).open(newline='') as stream:
                return list(csv.DictReader(stream))
        left, right = rows(reference), rows(current)
        require(len(left) == len(right) == count, 'Full report row count differs: ' + name)
        require(left == right, 'Current reference report fields differ: ' + name)
    left = {p.name: sha(p) for p in (reference / 'standard_diagnostics').glob('*.png')}
    right = {p.name: sha(p) for p in (current / 'standard_diagnostics').glob('*.png')}
    require(len(left) == 17 and left == right, 'Current reference standard diagnostic hashes differ')
    require(read(reference / 'closure.json') == read(current / 'closure.json'), 'Current reference closure differs')


def _numeric_gap(left, right, label):
    x, y = float(left), float(right)
    require(np.isfinite(x) == np.isfinite(y) and np.isnan(x) == np.isnan(y) and
            np.isposinf(x) == np.isposinf(y) and np.isneginf(x) == np.isneginf(y),
            'Nonfinite numeric masks differ: ' + label)
    gap = abs(x-y) if np.isfinite(x) else 0.
    require(gap <= BRIDGE_ABSOLUTE_TOLERANCE, 'Platform bridge numeric drift: ' + label)
    return gap


def verify_solution_arrays(saved, fresh):
    """Complete disclosed cross-platform comparison, not optimizer identity."""
    left = {key:value for key,value in vars(saved).items() if isinstance(value, np.ndarray)}
    right = {key:value for key,value in vars(fresh).items() if isinstance(value, np.ndarray)}
    require(left and set(left) == set(right), 'Complete saved/fresh ndarray inventory differs')
    controls_present = any(key in left for key in CONTROL_FIELDS)
    if controls_present:
        require(all(key in left for key in (*CONTROL_FIELDS,'V','birth_count_pre_distribution')),
                'Control comparison requires both complete policies and predecision distributions')
    dead = occupied = None
    if 'V' in left:
        require(left['V'].shape == right['V'].shape, 'Saved total-value shape differs')
        dead = left['V'] <= solver.DEAD_VALUE_CUTOFF
        require(np.array_equal(dead,right['V'] <= solver.DEAD_VALUE_CUTOFF), 'Total-value dead masks differ')
    if controls_present:
        require(left['birth_count_pre_distribution'].shape == right['birth_count_pre_distribution'].shape == dead.shape,
                'Occupied-control comparison state shape differs')
        occupied = (left['birth_count_pre_distribution'] > 1e-12) | (right['birth_count_pre_distribution'] > 1e-12)
    records = {}
    def errors(x,y):
        if not x.size:return dict(count=0,absolute_maximum_gap=0.,normalized_maximum_gap=0.)
        gap=np.abs(x.astype(np.float64)-y.astype(np.float64))
        scale=np.maximum(1.,np.maximum(np.abs(x.astype(np.float64)),np.abs(y.astype(np.float64))))
        return dict(count=int(x.size),absolute_maximum_gap=float(gap.max()),normalized_maximum_gap=float((gap/scale).max()))
    for key, x in left.items():
        y = right[key]
        require(x.shape == y.shape and x.dtype == y.dtype, 'Saved ndarray shape/dtype differs: ' + key)
        require(x.dtype.kind in 'biuf', 'Unsupported saved ndarray dtype: ' + key)
        finite = np.isfinite(x)
        require(np.array_equal(finite, np.isfinite(y)) and np.array_equal(np.isnan(x), np.isnan(y)) and
            np.array_equal(np.isposinf(x), np.isposinf(y)) and np.array_equal(np.isneginf(x), np.isneginf(y)),
            'Saved ndarray nonfinite masks differ: ' + key)
        if x.dtype.kind in 'biu':
            require(np.array_equal(x, y), 'Saved discrete ndarray differs: ' + key)
            metrics=errors(x,y); criterion='exact_discrete'; details={}
        else:
            a,b=x[finite].astype(np.float64),y[finite].astype(np.float64)
            metrics=errors(a,b); delta=np.abs(a-b);details={}
            if key == 'V':
                require(x.dtype == np.dtype(np.float64), 'Total value precision differs')
                bound=1e-10+16*np.finfo(np.float64).eps*np.maximum(np.abs(a),np.abs(b))
                criterion='value_absolute_plus_16_eps64_scale'
                details=dict(bound_formula='1e-10 + 16*eps64*max(abs(saved),abs(fresh))',
                    dead_count=int(dead.sum()),living_count=int((~dead).sum()),dead_masks_exact=True)
                require(np.all(delta <= bound),'Saved total-value numerical bound fails')
            elif key in CONTROL_FIELDS:
                require(x.shape == dead.shape and x.dtype == np.dtype(np.float64), 'Control precision/state shape differs: '+key)
                living=(~dead)&finite
                alive=errors(x[living],y[living]);occ=errors(x[occupied&finite],y[occupied&finite]);removed=errors(x[dead&finite],y[dead&finite])
                require(alive['normalized_maximum_gap'] <= np.sqrt(np.finfo(np.float64).eps),'Living control numerical bound fails: '+key)
                require(occ['absolute_maximum_gap'] <= 1e-10,'Occupied control absolute bound fails: '+key)
                criterion='living_control_scale_and_occupied_absolute'
                details=dict(living=alive,occupied_union=occ,dead_audit=removed,
                    bound_formula='living: sqrt(eps64)*max(1,abs(saved),abs(fresh)); occupied union: absolute<=1e-10',
                    dead_status=ARRAY_COMPARISON_CRITERIA['dead_control_status'],
                    dead_count=int(dead.sum()),living_count=int((~dead).sum()),occupied_union_count=int(occupied.sum()))
            elif x.dtype == np.dtype(np.float32):
                require(key.endswith('probs'), 'Float32 ULP criterion only supports probability arrays: '+key)
                require(np.array_equal(x==0,y==0),'Float32 probability zero support differs: '+key)
                magnitude=np.maximum(np.abs(x[finite]),np.abs(y[finite])).astype(np.float32)
                bound=2*np.spacing(magnitude).astype(np.float64)
                require(np.all(delta <= bound),'Float32 probability two-ULP bound fails: '+key)
                criterion='float32_probability_two_ulps'
                details=dict(bound_formula='2*spacing(float32(max(abs(saved),abs(fresh))))',zero_masks_exact=True,
                    zero_count=int((x==0).sum()),positive_support_count=int((x>0).sum()))
            else:
                require(metrics['absolute_maximum_gap'] <= 1e-10,'Saved ndarray numeric drift: '+key)
                criterion='float64_absolute_1e-10'
        records[key] = dict(shape=list(x.shape), dtype=str(x.dtype), **metrics,criterion=criterion,
            status='applicable_criteria_passed',nonfinite_count=int((~finite).sum()),details=details,
            saved_sha256=hashlib.sha256(x.tobytes()).hexdigest(), fresh_sha256=hashlib.sha256(y.tobytes()).hexdigest())
    return dict(array_count=len(records), criteria=ARRAY_COMPARISON_CRITERIA.copy(),all_applicable_criteria_passed=True,
        maximum_gap=max(record['absolute_maximum_gap'] for record in records.values()), arrays=records)


def compare_saved_platform_reports(reference, current):
    """Only derived report numerics and rendering may differ across platforms."""
    reference, current = Path(reference), Path(current)
    differences = []
    for name, count in (('target_fit.csv', 14), ('parameters.csv', 31)):
        def rows(folder):
            with (folder/name).open(newline='') as stream:
                return list(csv.DictReader(stream))
        left, right = rows(reference), rows(current)
        require(len(left) == len(right) == count, 'Full report row count differs: ' + name)
        for i, (x,y) in enumerate(zip(left,right)):
            require(set(x) == set(y), 'Platform bridge row structure differs: ' + name)
            for key in x:
                label = f'{name}.{i}.{key}'
                if name == 'target_fit.csv' and key in ('model','gap','loss_contribution') and x[key] != '':
                    gap = _numeric_gap(x[key],y[key],label)
                    if x[key] != y[key]: differences.append(dict(field=label, saved=x[key], fresh=y[key], absolute_gap=gap))
                else:
                    require(x[key] == y[key], 'Platform bridge exact report field differs: ' + label)
    def closure(x,y,label):
        require(type(x) is type(y), 'Platform bridge closure type differs: ' + label)
        if isinstance(x, dict):
            require(set(x) == set(y), 'Platform bridge closure keys differ: ' + label)
            for key in x: closure(x[key],y[key],label+'.'+key)
        elif isinstance(x, list):
            require(len(x) == len(y), 'Platform bridge closure shape differs: ' + label)
            for i,(a,b) in enumerate(zip(x,y)): closure(a,b,label+'.'+str(i))
        elif type(x) in (int,float):
            gap = _numeric_gap(x,y,label)
            if x != y: differences.append(dict(field=label,saved=x,fresh=y,absolute_gap=gap))
        else:
            require(x == y, 'Platform bridge closure metadata differs: ' + label)
    closure(read(reference/'closure.json'),read(current/'closure.json'),'closure')
    old = {p.name:sha(p) for p in (reference/'standard_diagnostics').glob('*.png')}
    new = {p.name:sha(p) for p in (current/'standard_diagnostics').glob('*.png')}
    require(len(old) == len(new) == 17 and set(old) == set(new), 'Platform bridge complete plot names differ')
    import matplotlib
    from matplotlib import ft2font
    return dict(derived_numeric_differences=differences,
        render_only_hash_differences={key:dict(saved=old[key],fresh=new[key]) for key in old if old[key] != new[key]},
        all_saved_plot_hashes=old, all_fresh_plot_hashes=new,
        fresh_render_versions=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__,
            platform=platform.platform(),freetype=ft2font.__freetype_version__),
        saved_render_versions='not recorded by the authenticated saved case; all saved hashes retained')


def attach_owned_count_policy(policy, P):
    """Capture full probabilities from this date's Bellman call before caching."""
    for key in COUNT_FIELDS:
        values = np.asarray(getattr(P, key, None))
        require(values.ndim == np.asarray(policy.V).ndim + 1 and values.shape[:-1] == policy.V.shape,
                'Date-owned full count policy unavailable: ' + key)
        setattr(policy, key, values.copy())
    return policy


class CurrentEstateARuntime(retained.FloorRuntime):
    housing = 'static-elastic'
    exact_policy_cache_bytes = 2 * 1024**3

    @classmethod
    def from_handoff(cls, pin, folder):
        self = cls()
        import numba
        requested_thread_text = os.environ.get('NUMBA_NUM_THREADS','1')
        requested_threads = int(requested_thread_text)
        require(1 <= requested_threads <= numba.config.NUMBA_NUM_THREADS,
                'Requested Numba threads exceed the initialized pool limit')
        self.folder = Path(folder); self.folder.mkdir(parents=True, exist_ok=True)
        self.handoff, self.sourcepins, case = authenticate_handoff(pin)
        self.handoff_pin = dict(pin)
        self.saved, self.case = storage.load_case(case)
        self.P = copy.deepcopy(self.saved.P); self.grid = self.saved.b_grid.copy()
        self.reference_price = float(self.saved.price)
        self.report = case / 'native/phase_b_ge/selected_root'
        closure = read(self.report / 'closure.json')
        validate_contract(self.P, closure, read(case / 'input_contract.json'))
        require(self.saved.closure['population_scale'] == closure['population_scale'], 'Saved/report population differs')
        self.population_scale = float(closure['population_scale'])
        with (self.report / 'parameters.csv').open(newline='') as stream:
            self.parameter_rows = list(csv.DictReader(stream))
        require(len(self.parameter_rows) == 31, 'Complete selected parameter table required')
        # Frozen reporting imports impose their historical one-thread env.
        # Restore the caller's explicit mask before any current kernel use;
        # never enlarge an already initialized Numba pool or silently fallback.
        try:
            self.ctx = reporting.build_context(self.P, self.grid, self.folder, price_start=self.reference_price,
                deadline=float('inf'), max_lifecycle=0, closure='fixed_h0')
        finally:
            os.environ['NUMBA_NUM_THREADS'] = requested_thread_text
            numba.set_num_threads(requested_threads)
        require(numba.get_num_threads() == requested_threads, 'Requested current Numba thread mask was not restored')
        self.rt = self.ctx['prepared'].rt; self.pf = self.rt['primitive'].pf
        # Reporting's recent-parent observer closes over this exact facade.
        # Reuse it so the retained caller/calendar identity guard remains valid.
        self.model = self.rt['model']
        self.ge = native_phase_b
        self.runner = SimpleNamespace(compare_repeated=self.compare_repeated)
        self.scaffold = retained.load('current_estate_a_dated_scaffold',
            Path(retained.__file__).parent / 'pinned_tools/run_e5f_preference_transition.py')
        self.packet = self.initial_state = None; self.reference_verified = False; self.total_native_calls = 0
        with self.native_bindings():
            actual = self.ctx['fp'].actual_parameters(self.ctx['prepared'], self.P, self.grid)
            self.ge.validate_parameter_estimates(dict(expected_parameters={r['parameter']:float(r['estimate'])
                for r in self.parameter_rows}), self.parameter_rows, actual)
            entrant = self.pf.calendar.entrant_cohort(np.array([1.]), self.P, self.grid)
            np.testing.assert_allclose(entrant.sum(axis=(1,2,4,5)),
                self.P.fixed_reference_entry_conditional*self.P.z_weights[None,:], rtol=0, atol=2e-16)
        write(self.folder / 'constructor.json', dict(status='current_estate_a_native_reconstruction_pending',
            policy_calls=0, population_scale=self.population_scale, identity=self.identity(),
            numba_threads=dict(requested=requested_threads, actual=numba.get_num_threads(),
                initialized_pool_limit=numba.config.NUMBA_NUM_THREADS, restored_after_frozen_reporting_import=True)))
        return self

    def stationary(self, psi, price, folder):
        packet, record = super().stationary(psi, price, folder)
        if float(psi) == float(self.P.psi_child) and float(price) == self.reference_price:
            arrays = verify_solution_arrays(self.saved.solution, packet['solution'])
            require(arrays['array_count'] == 78, 'Current selected case must verify all 78 solution ndarrays')
            native_receipt = read(Path(folder)/'native_solve_unverified.json')
            write(Path(folder)/'saved_array_verification.json', dict(schema='current_estate_a_saved_array_bridge_v1',
                status='passed',identity=self.identity(),saved_archive_sha256=sha(self.case/'native_result.npz'),
                current_report=str((Path(folder)/'phase_b_ge/selected_root').resolve()),
                native_checkpoint=native_receipt['checkpoint'], arrays=arrays))
        return packet, record

    def compare_repeated(self, reference, current):
        if Path(reference).resolve() != self.report.resolve():
            return compare_reports(reference,current)
        current = Path(current).resolve()
        proof_path = current.parent.parent/'saved_array_verification.json'
        require(proof_path.is_file(), 'Saved platform bridge requires verified solution arrays')
        proof = read(proof_path)
        require(proof.get('schema') == 'current_estate_a_saved_array_bridge_v1' and proof.get('status') == 'passed' and
            proof['identity'] == self.identity() and proof['saved_archive_sha256'] == sha(self.case/'native_result.npz') and
            proof['current_report'] == str(current), 'Saved platform array verification is stale or differently bound')
        checkpoint = proof['native_checkpoint']
        require(sha(checkpoint['path']) == checkpoint['sha256'], 'Verified native array checkpoint changed')
        arrays = proof['arrays']
        expected = {key for key,value in vars(self.saved.solution).items() if isinstance(value,np.ndarray)}
        require(arrays['array_count'] == 78 and set(arrays['arrays']) == expected and
            arrays['criteria'] == ARRAY_COMPARISON_CRITERIA and arrays['all_applicable_criteria_passed'] is True,
            'Complete saved array verification missing or criterion differs')
        # Recompute from the authenticated native checkpoint so a changed
        # criterion, status, gap, hash or support classification cannot pass.
        with gzip.open(checkpoint['path'],'rb') as stream:
            native = pickle.load(stream)
        require(native['identity'] == self.identity(), 'Native array checkpoint identity differs')
        recomputed = verify_solution_arrays(self.saved.solution,native['solution'])
        require(arrays == recomputed, 'Saved array verification evidence changed')
        comparison = compare_saved_platform_reports(reference,current)
        write(current/'saved_platform_bridge.json',dict(schema='current_estate_a_saved_platform_bridge_v1',
            status='passed',identity=self.identity(),array_verification=dict(path=str(proof_path),sha256=sha(proof_path)),
            absolute_tolerance=BRIDGE_ABSOLUTE_TOLERANCE, numerical_gates_unchanged=True,
            array_comparison_criteria=ARRAY_COMPARISON_CRITERIA, optimizer_identity_certified=False,
            same_host_repeat_remains_exact=True, **comparison))
        return comparison

    def identity(self):
        def digest(value):
            return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        engine = {key:value for key,value in self.sourcepins.items() if '/birth_count_choice/model/engine/' in key}
        return dict(reference_sha256=self.handoff_pin['sha256'], engine_sha256=digest(engine),
            entry_sha256=hashlib.sha256(self.P.fixed_reference_entry_conditional.tobytes()).hexdigest(),
            grid_sha256=hashlib.sha256(self.grid.tobytes()).hexdigest(),
            effective_parameters_sha256=digest({row['parameter']:row['estimate'] for row in self.parameter_rows}),
            source_pins=self.sourcepins)

    @contextlib.contextmanager
    def native_api_bindings(self):
        # This is the complete inherited compatibility API, flattened from the
        # same actual stage objects used by production_model_facade. No helper
        # is borrowed from the old floor or small-credit engine.
        required = ('__file__', 'interp_indices', 'birth_destination_child_state', 'precompute_shared',
            'solve_bellman_full_markov_income', 'solve_markov_income_at_prices',
            'forward_distribution_markov_income', 'advance_cohort_one_period_markov_income',
            'property_tax_revenue_from_distribution', 'build_forward_tenure_transition_maps',
            'income_transition_values', 'entry_wealth_grid_weights', 'DEAD_MASS_TOL',
            '_censor_entry_dead_mass', '_gate_dead_mass_at_age', 'realize_current_cross_section',
            'realize_stayer_cross_section', 'age_to_index', 'apply_child_aging', 'bequest_utility_vec',
            'compute_markov_statistics', 'get_completed_fertility', 'income_at_state',
            'owner_borrowing_floor', 'renter_borrowing_floor', 'independent_child_maturation_active',
            'get_fecundity_by_age', 'readiness_settled_state', 'parent_age_maturation_active',
            'readiness_childless_states', 'readiness_gate_active', 'DEAD_VALUE_CUTOFF',
            'add_aggregate_wealth_bequest_flow_moments', 'annual_gross_income_at_state')
        expected = {}
        for module in (parameters, utils, adult_entry, solver):
            expected.update({name:value for name,value in vars(module).items() if not name.startswith('__')})
        expected['__file__'] = solver.__file__
        for name in required:
            value = getattr(self.model, name, None)
            require(name in expected and value is expected[name], 'Native compatibility object differs: ' + name)
            if name in ('__file__', 'DEAD_MASS_TOL', 'DEAD_VALUE_CUTOFF'):
                continue
            require(callable(value) and value.__module__.startswith(solver.__package__ + '.'),
                    'Callback is not the selected current engine: ' + name)
        require(self.model.DEAD_MASS_TOL == 1e-12, 'Dead-mass gate differs')
        yield {name:getattr(getattr(self.model,name), '__module__', solver.__name__) for name in required}

    @contextlib.contextmanager
    def native_bindings(self):
        require(self.P.birth_count_choice_cap == 1 and self.P.child_state_mode == 'independent_count' and
            not self.P.readiness_gate_enabled, 'Dated one-birth independent-count contract differs')
        require(not bool(getattr(self.P, 'joint_nested_choice', False)), 'Legacy joint factorization is forbidden')
        require(self.model.solve_bellman_full_markov_income is household.solve_bellman_full_markov_income,
                'Wrong current backward callee')
        require(self.model.forward_distribution_markov_income is distribution.forward_distribution_markov_income,
                'Wrong current stationary forward callee')
        require(self.model.advance_cohort_one_period_markov_income is distribution.advance_cohort_one_period_markov_income,
                'Wrong current dated forward callee')
        for module in (household, distribution):
            for name,value in vars(module).items():
                if callable(value) and getattr(value, '__module__', '').endswith('.kernels'):
                    require(value is getattr(kernels, name), 'Wrong current kernel callee: ' + name)
        cal = self.pf.calendar; transition = self.pf.transition
        require(transition.calendar is cal, 'Dated calendar module differs')
        require(cal.apply_fertility is apply_count_fertility and transition.apply_sequential_fertility is apply_count_fertility,
                'Dated fertility callback must use the owned count kernel')
        factory = self.pf.policy_from_objects
        configure = transition.configure_sequential_model
        import e5f_overnight_estate_audit as estate
        from .model.estate_audit_adapter import adapt_audit
        old_estate = estate.audit
        estate_audit, estate_receipt = adapt_audit(getattr(estate, '_estate_a_original_audit', old_estate))
        def count_factory(objects, price, P, b_grid, shared):
            return attach_owned_count_policy(factory(objects, price, P, b_grid, shared), P)
        def forbidden(*args, **kwargs):
            raise RuntimeError('Legacy sequential reconfiguration forbidden inside current callbacks')
        targets = [(cal, 'model', self.model), (self.rt['primitive'], 'model', self.model),
                   (self.rt['audit'], 'model', self.model)]
        with retained._LOCK, self.native_api_bindings():
            saved = [(obj,key,getattr(obj,key)) for obj,key,_ in targets]; previous = self.rt['model']
            try:
                for obj,key,value in targets: setattr(obj,key,value)
                self.rt['model'] = self.model
                transition.configure_sequential_model = forbidden
                self.pf.policy_from_objects = count_factory
                estate.audit = estate_audit
                yield
            finally:
                changed = cal.model is not self.model
                self.pf.policy_from_objects = factory; transition.configure_sequential_model = configure
                estate.audit = old_estate
                self.rt['model'] = previous
                for obj,key,value in saved: setattr(obj,key,value)
                require(not changed, 'Dated callback replaced the selected native engine')

    def bootstrap_saved_reference(self, folder):
        """Reconstruct saved policy/distribution without a new Bellman or KFE solve.

        This is a bootstrap check, not permission for a dated mapping. The
        controller still requires the inherited two fresh selected-price repeats.
        """
        folder = Path(folder); folder.mkdir(parents=True, exist_ok=True)
        P = copy.deepcopy(self.P); sol = copy.deepcopy(self.saved.solution)
        with self.native_bindings():
            shared = self.model.precompute_shared(P, self.grid)
            P._fert2_probs = sol.fert2_probs.copy()
            policy = self.pf.calendar.policy_from_solution(sol, np.array([self.reference_price]), P, self.grid, shared)
            pre, reconstruction = self.pf.calendar.reconstruct_stationary_pre_fertility(sol, policy, P, self.grid, shared)
            require(reconstruction['stationary_post_fertility_nesting_l1'] <= 5e-9 and
                reconstruction['stationary_feasibility_projection_mass'] == 0., 'Saved stationary reconstruction fails')
            np.testing.assert_allclose(pre, sol.birth_count_pre_distribution, rtol=0, atol=5e-9)
            packet = dict(parameters=P, b_grid=self.grid, shared=shared, solution=sol,
                stationary_g_pre=pre, evaluation=SimpleNamespace(births=float(sol.birth_count_expected_children_by_age.sum())))
            state = self.stationary_state(packet, self.population_scale)
            require(abs(float(state.g_pre.sum())/self.population_scale - 1.) <= 1e-9, 'Bootstrap inherited population differs')
        receipt = dict(status='saved_policy_bootstrap_passed', policy_calls=0, reference_verified=False,
            reconstruction=reconstruction, actual_initial_population=float(state.g_pre.sum()),
            adjusted_queue=self.pf.birth_queue_values(state.scheduled_entries).tolist(),
            raw_queue=self.pf.birth_queue_values(state.scheduled_raw_entries).tolist(), identity=self.identity())
        write(folder / 'bootstrap.json', receipt)
        return receipt


# Compatibility seam: the existing controller resolves module.FloorRuntime.
FloorRuntime = CurrentEstateARuntime
