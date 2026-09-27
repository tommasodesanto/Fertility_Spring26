"""Opt-in, runtime-only credit/PF adapters; no imports or solves of the model.

Queue lists serialize the native three-plus-four waiting slots. Fiscal choices
are deliberately caller-owned: supply anticipated pensions and the inherited
fixed payroll tax to the native PF evaluator; no rebate is imposed here.
"""
from contextlib import nullcontext
import difflib
import hashlib
import inspect
import math
from pathlib import Path
import numpy as np


def _sha(source):
    return hashlib.sha256(source.encode()).hexdigest()


def _once(source, old, new):
    if source.count(old) != 1:
        raise ValueError('Source anchor must occur exactly once: ' + old[:80])
    return source.replace(old, new, 1)


def continuation_source(source):
    """Remove only the dated-continuation exclusion, keep exhaustive control."""
    return _once(source,
        'if continuation_V is not None or not bool(getattr(P, "exhaustive_saving_control", False)):',
        'if not bool(getattr(P, "exhaustive_saving_control", False)):')


def _save_patch(output, name, before, after):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    path = output / (name + '.generated.py')
    path.write_text(after)
    (output / (name + '.diff')).write_text(''.join(difflib.unified_diff(
        before.splitlines(True), after.splitlines(True), fromfile=name + '_before', tofile=name + '_after')))
    return path


def install_credit(model, helper, output, *, enabled=False, expected_helper_sha256=None):
    """Fresh-runtime only; preserve completed stationary benchmark unchanged."""
    if not enabled:
        return dict(enabled=False, baseline_noop=True)
    helper_path = Path(inspect.getfile(helper))
    actual = hashlib.sha256(helper_path.read_bytes()).hexdigest()
    if not expected_helper_sha256 or actual != expected_helper_sha256:
        raise ValueError('Explicit authenticated natural-credit helper hash required')
    installed = helper.install(model, Path(output) / 'stationary_credit', enabled=True)
    before = inspect.getsource(model.solve_bellman_full_markov_income)
    after = continuation_source(before)
    path = _save_patch(output, 'dated_credit_bellman', before, after)
    exec(compile(after, str(path), 'exec'), model.__dict__)
    return dict(enabled=True, helper_sha256=actual, stationary_installation=installed,
                before_sha256=_sha(before), after_sha256=_sha(after),
                death_valuation='Inherited current decision price; post-saving net liquidation solvency',
                limitations='Grid-node/native value-cutoff approximation; dynamic numerical validation outstanding')


def flatten_queue(queue):
    return list(queue.due_in_16) + list(queue.due_in_20)


def split_queue_step(slots, births, conversion, queue_type):
    if not math.isclose(float(conversion), 1 / 2.1, rel_tol=0, abs_tol=1e-15):
        raise ValueError('Retain birth-to-household conversion exactly once')
    values = tuple(float(x) for x in slots)
    if len(values) != 7:
        raise ValueError('Split queue requires exactly seven slots (three plus four)')
    queue = queue_type(values[:3], values[3:])
    due, future = queue.step(float(births))
    return due, flatten_queue(future)


def initialize_state(pf, queue_type, g_pre, adjusted_births, raw_births):
    """Stationary birth prehistory, without forcing births to actual entry mass."""
    distribution = np.asarray(g_pre, dtype=float)
    if distribution.ndim != 7 or not np.isfinite(distribution).all() or np.any(distribution < 0):
        raise ValueError('Exact finite nonnegative pre-fertility distribution required')
    return pf.PFInitialState(g_pre=distribution.copy(),
        scheduled_entries=flatten_queue(queue_type.constant_prehistory(adjusted_births)),
        scheduled_raw_entries=flatten_queue(queue_type.constant_prehistory(raw_births)))


def queue_source(source):
    source = _once(source, '    started = time.perf_counter()\n', '''    started = time.perf_counter()
    if historical_conditioning is not None:
        raise ValueError("Credit transition excludes historical age reweighting and migration")
    if getattr(base_parameters, "adult_entry_clock", None) != "split_birth_vintage":
        raise ValueError("Adopted split birth-vintage parameter required")
''')
    if source.count('transition.advance_birth_vintage_queue(') != 2:
        raise ValueError('Expected both raw and adjusted native queue calls')
    return source.replace('transition.advance_birth_vintage_queue(', '_credit_split_queue_step(')


def install_split_queue(pf, output, queue_type=None, *, enabled=True, expected_source_sha256=None):
    """Patch only PF queue calls; retain native household and fiscal mapping."""
    if not enabled:
        return dict(enabled=False, baseline_noop=True)
    if getattr(pf, '_credit_split_queue_installed', False):
        raise ValueError('Split PF queue already installed')
    before = inspect.getsource(pf.evaluate_path_at_prices)
    if expected_source_sha256 is not None and _sha(before) != expected_source_sha256:
        raise ValueError('Native PF function source hash differs')
    if queue_type is None:
        queue_type = pf.transition.SplitBirthEntryQueue
    queue_source_path = Path(inspect.getfile(queue_type))
    queue_file_sha256 = hashlib.sha256(queue_source_path.read_bytes()).hexdigest()
    # Check this is the native three/four-slot queue interface, not an ad hoc replacement.
    probe = queue_type.constant_prehistory(2.1)
    if len(probe.due_in_16) != 3 or len(probe.due_in_20) != 4:
        raise ValueError('Native queue interface mismatch')
    after = queue_source(before)
    path = _save_patch(output, 'split_queue_pf', before, after)
    pf.__dict__['_credit_split_queue_step'] = lambda q, b, c: split_queue_step(q, b, c, queue_type)
    exec(compile(after, str(path), 'exec'), pf.__dict__)
    pf._credit_split_queue_installed = True
    return dict(enabled=True, before_sha256=_sha(before), after_sha256=_sha(after),
                queue_source=str(queue_source_path), queue_source_sha256=queue_file_sha256,
                waiting_slots=[3, 4], entry_dates_after_birth=[4, 5], conversion=1 / 2.1,
                fiscal_closure='caller supplies inherited fixed tax and jointly balanced anticipated pensions')


def exact_cache_context(pf, cache_module, *, enabled=False, max_bytes=0):
    """Reuse existing bounded exact-call cache, never approximate continuation."""
    if not enabled:
        return nullcontext(None)
    if isinstance(max_bytes, bool) or not isinstance(max_bytes, int) or max_bytes <= 0:
        raise ValueError('Explicit positive cache byte budget required')
    if not callable(getattr(pf, 'solve_date_policy', None)):
        raise ValueError('Native dated policy function required')
    return cache_module.policy_cache(pf, max_bytes=max_bytes)


def install_continuation(model, output):
    """Enable dated continuation after authenticated stationary-credit install."""
    if not getattr(model, '_solvency_credit_benchmark_installed', False):
        raise ValueError('Authenticated stationary solvency adapter must be installed first')
    before = inspect.getsource(model.solve_bellman_full_markov_income)
    after = continuation_source(before)
    path = _save_patch(output, 'dated_credit_bellman', before, after)
    exec(compile(after, str(path), 'exec'), model.__dict__)
    return dict(before_sha256=_sha(before), after_sha256=_sha(after),
                death_valuation='Inherited current decision price',
                status='dynamic numerical regression required')


def initial_state(pf, packet, queue_type=None):
    """Preserve saved pre-fertility state and historically initialized entry prehistory."""
    if queue_type is None:
        queue_type = pf.transition.SplitBirthEntryQueue
    P = packet['parameters']
    if getattr(P, 'adult_entry_clock', None) != 'split_birth_vintage':
        raise ValueError('Packet must use adopted split birth-vintage entry')
    evaluation = packet['evaluation']
    accounting = pf.transition.calendar_topcode_birth_accounting(
        evaluation.g_pre, evaluation.g_post_fertility, float(evaluation.births), P)
    # Historical stationary prehistory is pinned to actual entrant mass.
    # Future queue updates still append actual endogenous births.
    entrant_mass = float(np.asarray(packet['stationary_g_pre'])[:, :, :, 0].sum())
    return initialize_state(pf, queue_type, packet['stationary_g_pre'],
        2.1 * entrant_mass, float(evaluation.births))
