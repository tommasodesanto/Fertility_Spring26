"""Authenticated in-memory estate-A exception to the donor-valuation guard."""
from __future__ import annotations
import hashlib
import inspect
from pathlib import Path

AUDIT_SOURCE_PINS = {
    '31bcdbeca73036da9306dcb2f3ef28a2c3ee43d878f60b75aa044e08aba5f962',
    'ce9e5c5e70c2e85fd86ebfcaa1a710e52a710277351b1c2cd651285ac0d7c8dc',
}
IDENTITY = 'estate_a_postsaving_net_selling_cost_v1'
ANCHOR = 'or bool(getattr(P, "bequest_net_of_selling_cost", False))'
REPLACEMENT = ANCHOR + ' and not _estate_a_allowed(P)'


def _estate_a_allowed(P):
    utility = getattr(P, 'bequest_net_of_selling_cost', False)
    flow = getattr(P, 'estate_flow_net_of_selling_cost', False)
    if not bool(utility) and not bool(flow):
        return False
    if (utility is not True or flow is not True
            or str(getattr(P, 'estate_receiver', 'none')) != 'none'
            or any(float(getattr(P, name, 0.)) != 0. for name in
                   ('estate_lump_sum_transfer', 'estate_probe_transfer', 'estate_tax_rate'))):
        raise ValueError('Estate-A audit requires both Boolean net-estate flags, receiver none, and zero estate transfers/tax')
    return True


def _label(ledger, receipt):
    ledger['donor_bequest_valuation_changes'] = True
    ledger['estate_a_identity'] = IDENTITY
    ledger['donor_estate_valuation'] = 'post-saving bp+(1-psi)*price*chosen_house; no extra interest'
    ledger['estate_a_audit_receipt'] = receipt
    old = 'Availability valuation is provisional; donor utility remains unchanged.'
    if ledger['caveats'].count(old) != 1:
        raise RuntimeError('Estate audit caveat source drift')
    ledger['caveats'] = [
        'Availability valuation is provisional; experimental estate-A donor utility uses the same net housing valuation.'
        if item == old else item for item in ledger['caveats']]
    return ledger


def adapt_audit(original):
    """Clone only one guard clause; every ledger calculation/gate is inherited."""
    source_path = inspect.getsourcefile(original)
    file_pin = hashlib.sha256(Path(source_path).read_bytes()).hexdigest()
    if file_pin not in AUDIT_SOURCE_PINS:
        raise RuntimeError('Estate audit source is not an authenticated supported pin')
    before = inspect.getsource(original)
    if before.count(ANCHOR) != 1:
        raise RuntimeError('Estate audit donor guard source drift')
    after = before.replace(ANCHOR, REPLACEMENT, 1)
    namespace = dict(original.__globals__, _estate_a_allowed=_estate_a_allowed)
    exec(compile(after, '<authenticated_estate_a_audit_clone>', 'exec'), namespace)
    clone = namespace[original.__name__]
    receipt = dict(original_source_path=source_path, original_file_sha256=file_pin,
        original_function_sha256=hashlib.sha256(before.encode()).hexdigest(),
        clone_function_sha256=hashlib.sha256(after.encode()).hexdigest(),
        exact_replacement={'before': ANCHOR, 'after': REPLACEMENT},
        identity=IDENTITY, all_accounting_arithmetic_and_numerical_gates_unchanged=True,
        parameter_mutation=False, frozen_source_mutation=False)
    shortfall = original.__globals__['EstateFundingShortfall']
    def audit(evaluation, P, b_grid, **kwargs):
        active = _estate_a_allowed(P)
        if not active:
            return original(evaluation, P, b_grid, **kwargs)
        try:
            ledger = clone(evaluation, P, b_grid, **kwargs)
        except shortfall as exc:
            _label(exc.audit, receipt)
            raise
        return _label(ledger, receipt)
    return audit, receipt


def install_estate_audit(context, P, out):
    """Called only after reporting's frozen authentication has succeeded."""
    if not _estate_a_allowed(P):
        return
    prepared = context['prepared']
    contract = prepared.estate
    module = getattr(contract, 'module', contract)
    original = getattr(module, '_estate_a_original_audit', module.audit)
    adapted, receipt = adapt_audit(original)
    module._estate_a_original_audit = original
    module.audit = adapted
    context['fp'].write(Path(out)/'estate_a_audit_adapter_receipt.json', receipt)
