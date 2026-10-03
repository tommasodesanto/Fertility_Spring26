"""Estate-A net death flow in the authenticated empirical wealth observer."""
from __future__ import annotations
import hashlib
import inspect
from pathlib import Path
from types import SimpleNamespace
from .estate_audit_adapter import _estate_a_allowed, IDENTITY

SOURCE_PIN = 'f711a61e8e177eb293d6ba58f8c6e7e6dc83f1eadb4af1cb1cf8522c055db435'
ANCHOR = '    model.add_aggregate_wealth_bequest_flow_moments(\n'
REPLACEMENT = '    _estate_a_aggregate_moments(model,\n'


def aggregate_moments(facade, legacy, stats, wealth_g, death_g, bp, P, bg, prices):
    # Keep the frozen living-stock/earnings calculation and all its outputs.
    legacy.add_aggregate_wealth_bequest_flow_moments(stats, wealth_g, death_g, bp, P, bg, prices)
    if _estate_a_allowed(P):
        net = SimpleNamespace()
        facade.add_aggregate_wealth_bequest_flow_moments(net, wealth_g, death_g, bp, P, bg, prices)
        stats.annual_bequest_flow = net.annual_bequest_flow
        stats.annual_bequest_flow_to_aggregate_wealth = net.annual_bequest_flow / max(stats.aggregate_wealth, 1e-12)


def adapt_observer(function, facade):
    """Clone its private wealth helper; retain the public observer identity."""
    helper = function.__globals__['_wealth_diagnostics']
    original = getattr(helper, '_estate_a_original_helper', helper)
    source_path = inspect.getsourcefile(original)
    file_pin = hashlib.sha256(Path(source_path).read_bytes()).hexdigest()
    if file_pin != SOURCE_PIN:
        raise RuntimeError('Housing/wealth observer source pin drift')
    before = inspect.getsource(original)
    if before.count(ANCHOR) != 1:
        raise RuntimeError('Housing/wealth observer balance-sheet call source drift')
    after = before.replace(ANCHOR, REPLACEMENT, 1)
    namespace = dict(original.__globals__)
    namespace['_estate_a_aggregate_moments'] = lambda *args: aggregate_moments(facade, *args)
    exec(compile(after, '<authenticated_estate_a_wealth_observer_clone>', 'exec'), namespace)
    clone = namespace[original.__name__]
    clone._estate_a_original_helper = original
    function.__globals__['_wealth_diagnostics'] = clone
    return dict(identity=IDENTITY, original_source_path=source_path,
        original_file_sha256=file_pin,
        original_function_sha256=hashlib.sha256(before.encode()).hexdigest(),
        clone_function_sha256=hashlib.sha256(after.encode()).hexdigest(),
        exact_replacement={'before':ANCHOR, 'after':REPLACEMENT},
        modified_statistics=['annual_bequest_flow','annual_bequest_flow_to_aggregate_wealth'],
        estate_formula='bp+(1-psi)*price*chosen_house; positive-part death-weighted and annualized once; no extra R',
        living_wealth_and_earnings_unchanged=True, age_geometry_and_state_weights_unchanged=True,
        targets_unchanged=True, frozen_source_mutation=False,
        net_flow_primitive_source=inspect.getsourcefile(facade.add_aggregate_wealth_bequest_flow_moments))


def install_estate_observer(context, P, facade, out):
    if not _estate_a_allowed(P):
        return
    function = context['prepared'].rt['observe_initial_housing_wealth']
    receipt = adapt_observer(function, facade)
    context['fp'].write(Path(out)/'estate_a_wealth_observer_receipt.json', receipt)
