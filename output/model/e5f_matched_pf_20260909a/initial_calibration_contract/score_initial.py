"""Pure twelve-moment diagonal minimum-distance score; no I/O or model solves.

``score_initial`` consumes unchanged observation dictionaries and a supplied,
externally pinned contract. No weights, samples, target values, or certification
are inferred here. ``fingerprint`` hashes canonical JSON, not original file bytes.
The caller must supply the independently expected contract hash, not recompute it
from an untrusted replacement. A receipt binds each observation and parameter
packet to one checkpoint and the contract's complete source-fingerprint map.
This checks receipt consistency; it does not independently re-solve or certify
an equilibrium, empirical provenance, or local/global identification.

Contract schema ``e5f_initial_minimum_distance_v1`` requires:
* contract_id, objective_name; source_fingerprints (nonempty name:SHA256 map);
* cps_projection, explicitly uniform_birth_time or constant_post_cell;
* target_rows: all thirteen named rows. Preserve empirical_builder, sample,
  definition, empirical_provenance_contract_id, empirical_record_id and
  empirical_source_path. The twelve scored rows require positive actual_weight,
  model_observation and role='proposed_scored_restriction'; normalization has
  its existing separate role, target 2.1 and actual_weight=None;
* parameter_restrictions: the nine named coordinates, lower/upper/transform;
  near_bound_fraction in [0,.5], normalization_tolerance in (0,5e-4];
* recent_parent_approximation: maintained=True, exact_acs=False, approximation_id,
  declaration, observer_id, moment, snapshot, age_projection and
  diagnostic_allow_residence_proxy=True. This explicitly connects the preserved
  ACS target to a distinct maintained model approximation.

Evaluation receipt schema ``e5f_initial_score_receipt_v1`` requires status='verified',
source_fingerprints, checkpoint_sha256, numerical_gates_verified=True,
recent_parent_certified=True, recent_parent_approximation_id, and input_sha256
with hashes of early_measurement, recent_parent_observation, normalization and
parameters. ``parameters`` is the full list of supplied parameter-table rows;
numeric CSV cells are accepted there alone. Every structural row must carry the
supplied bounds and transform; extra fixed/derived rows are retained in full.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
from collections.abc import Mapping

SCHEMA = 'e5f_initial_minimum_distance_v1'
RECENT_MOMENT = 'ownership_current_birth_from_empty_dependent_home_minus_current_empty_home_30_55'
RECENT_OBSERVER = 'e5f_recent_parent_flow_diagnostic_v1'
RECENT_TARGET = 0.16289550916123285
PARAMETERS = ('beta_annual', 'kappa_fert', 'kappa_fert_continuation', 'chi', 'H0',
              'theta0', 'theta1', 'first_birth_fixed_cost', 'h_P')
NORMALIZATION_ID = 'initial_normalization'
MOMENT_PATHS = {
    'cps_childlessness': 'fertility.{projection}.moments.childless_rate_40_44',
    'cps_exactly_one': 'fertility.{projection}.moments.exactly_one_among_mothers_40_44',
    'nchs_mean_age': 'fertility.{projection}.moments.period_mean_age_first_birth',
    'nchs_share30': 'fertility.{projection}.moments.period_share_first_births_age30plus',
    'wealth_earnings': 'housing_wealth.moments.aggregate_wealth_to_annual_gross_labor_earnings',
    'bequest_wealth': 'housing_wealth.moments.annual_bequest_flow_to_aggregate_wealth',
    'old_dispersion': 'housing_wealth.moments.old_total_wealth_to_annual_income_p90_p50_7684',
    'mean_rooms': 'housing_wealth.moments.aggregate_mean_occupied_rooms_capped9_18_85',
    'ownership_30_55': 'housing_wealth.moments.own_rate_30_55',
    'first_birth_rooms': 'housing_wealth.moments.housing_increment_0to1',
    'family_rooms': 'housing_wealth.moments.prime30_55_model_dependent_3plus_minus_1to2_rooms_capped9',
    'recent_parent_ownership': RECENT_MOMENT,
}
PROVENANCE_FIELDS = ('empirical_builder', 'sample', 'definition',
                     'empirical_provenance_contract_id', 'empirical_record_id',
                     'empirical_source_path')


class ScoreContractError(ValueError):
    """A required contract, observation, parameter or receipt is inadmissible."""


def _require(condition, message):
    if not condition:
        raise ScoreContractError(message)


def fingerprint(value):
    """SHA256 of strict canonical JSON; reject all nonfinite/unserializable data."""
    try:
        payload = json.dumps(value, sort_keys=True, separators=(',', ':'),
                             ensure_ascii=True, allow_nan=False).encode('utf-8')
    except (TypeError, ValueError) as exc:
        raise ScoreContractError('Cannot fingerprint nonfinite or non-JSON data') from exc
    return hashlib.sha256(payload).hexdigest()


def _hash(value, label):
    _require(isinstance(value, str) and len(value) == 64
             and all(c in '0123456789abcdef' for c in value), label + ' must be SHA256')
    return value


def _text(value, label):
    _require(isinstance(value, str) and bool(value.strip()), label + ' is required')
    return value


def _number(value, label, csv_cell=False):
    if csv_cell and isinstance(value, str):
        try:
            value = float(value)
        except ValueError as exc:
            raise ScoreContractError(label + ' is not numeric') from exc
    _require(type(value) in (int, float) and math.isfinite(value), label + ' must be finite')
    return float(value)


def _rows(rows, key, label):
    _require(isinstance(rows, list) and all(isinstance(x, Mapping) for x in rows),
             label + ' must be a list of rows')
    indexed = {}
    for row in rows:
        name = _text(row.get(key), label + '.' + key)
        _require(name not in indexed, 'Duplicate ' + label + ': ' + name)
        indexed[name] = row
    return indexed


def _lookup(data, path):
    for key in path.split('.'):
        _require(isinstance(data, Mapping) and key in data, 'Missing model observation: ' + path)
        data = data[key]
    return _number(data, 'Model observation ' + path)


def _validate_recent(packet, approximation, checkpoint):
    _require(isinstance(approximation, Mapping), 'Recent-parent approximation declaration required')
    _require(approximation.get('maintained') is True and approximation.get('exact_acs') is False,
             'Recent-parent approximation must be explicitly maintained and not exact ACS')
    for name in ('approximation_id', 'declaration'):
        _text(approximation.get(name), 'recent_parent_approximation.' + name)
    expected = {'observer_id': RECENT_OBSERVER, 'moment': RECENT_MOMENT,
                'snapshot': 'synchronized_post_fertility_snapshot',
                'age_projection': 'uniform_within_age_cell',
                'diagnostic_allow_residence_proxy': True}
    for key, value in expected.items():
        _require(approximation.get(key) == value, 'Unsupported recent-parent approximation: ' + key)
    _require(packet.get('observer_id') == RECENT_OBSERVER and packet.get('moment') == RECENT_MOMENT,
             'Wrong recent-parent model observation')
    _require(packet.get('available') is True, 'Recent-parent observation unavailable')
    _require(packet.get('status') == 'diagnostic_only_measurement_approximations_unresolved',
             'Expected unchanged recent-parent diagnostic packet')
    _require(packet.get('production_eligible') is False and packet.get('target_contract_activated') is False,
             'Do not relabel the original recent-parent observer as an exact activated target')
    metadata = packet.get('metadata', {})
    for key in ('snapshot', 'age_projection', 'diagnostic_allow_residence_proxy'):
        _require(metadata.get(key) == expected[key], 'Recent-parent packet convention mismatch: ' + key)
    _require(metadata.get('age_interval') == [30., 56.], 'Recent-parent age interval mismatch')
    _require(metadata.get('policy_input_provenance', {}).get('checkpoint_sha256') == checkpoint,
             'Recent-parent observation belongs to a different checkpoint')
    value = _number(packet.get('model_value'), 'recent-parent model_value')
    rates = []
    for name in ('selected_birth', 'current_empty'):
        group = packet.get('groups', {}).get(name, {})
        den = _number(group.get('denominator'), name + ' denominator')
        num = _number(group.get('owner_numerator'), name + ' numerator')
        rate = _number(group.get('ownership_rate'), name + ' ownership_rate')
        _require(den > 0 and 0 <= num <= den and 0 <= rate <= 1,
                 'Invalid recent-parent group: ' + name)
        _require(math.isclose(num / den, rate, rel_tol=0, abs_tol=1e-12),
                 'Recent-parent group ratio mismatch: ' + name)
        rates.append(rate)
    _require(math.isclose(value, rates[0] - rates[1], rel_tol=0, abs_tol=1e-12),
             'Recent-parent ownership contrast mismatch')
    return value


def _parameter_table(parameters, restrictions, near_fraction):
    supplied = _rows(parameters, 'parameter', 'parameter')
    bounds = _rows(restrictions, 'parameter', 'parameter restriction')
    _require(set(bounds) == set(PARAMETERS), 'Exactly nine named structural restrictions are required')
    _require(set(PARAMETERS).issubset(supplied), 'Missing structural parameter')
    result = []
    for original in parameters:
        row = copy.deepcopy(original)
        name = row['parameter']
        value = _number(row.get('estimate'), name + ' estimate', csv_cell=True)
        row['estimate'] = value
        if name in bounds:
            restriction = bounds[name]
            low = _number(restriction.get('lower'), name + ' contract lower')
            high = _number(restriction.get('upper'), name + ' contract upper')
            _require(low < high and low <= value <= high, 'Parameter outside its bounds: ' + name)
            _require(_number(row.get('lower'), name + ' lower', True) == low
                     and _number(row.get('upper'), name + ' upper', True) == high,
                     'Parameter bound mismatch: ' + name)
            transform = _text(restriction.get('transform'), name + ' contract transform')
            _require(row.get('transform') == transform, 'Parameter transform mismatch: ' + name)
            row.update(lower=low, upper=high, near_bound=min(value-low, high-value) <= near_fraction*(high-low),
                       structural_coordinate=True)
        else:
            _text(row.get('status'), name + ' fixed/derived restriction status')
            low, high = row.get('lower'), row.get('upper')
            if low not in (None, '') or high not in (None, ''):
                low = _number(low, name + ' lower', True); high = _number(high, name + ' upper', True)
                _require(low <= value <= high, 'Fixed/derived parameter outside supplied bounds: ' + name)
            row['structural_coordinate'] = False
        result.append(row)
    return result


def score_initial(contract, *, expected_contract_sha256, early_measurement,
                  recent_parent_observation, normalization, parameters, evaluation_receipt):
    """Return twelve weighted squared gaps, all thirteen rows and full parameters.

    ``loss = fsum(actual_weight * (model - target)**2)`` over the twelve named
    scored restrictions only. The 2.1 normalization must pass its own tolerance;
    it has null weight/contribution and never affects the objective. The two CPS
    projections require separate explicit contracts; averaging is unsupported.
    No supplied input is modified. Invalid or incomplete inputs raise
    ``ScoreContractError`` rather than returning a penalty or dropping rows.
    """
    _require(isinstance(contract, Mapping), 'Score contract must be a mapping')
    for name, packet in (('early_measurement', early_measurement),
                         ('recent_parent_observation', recent_parent_observation),
                         ('normalization', normalization), ('evaluation_receipt', evaluation_receipt)):
        _require(isinstance(packet, Mapping), name + ' must be a mapping')
    _hash(expected_contract_sha256, 'expected_contract_sha256')
    actual_hash = fingerprint(contract)
    _require(actual_hash == expected_contract_sha256, 'Contract fingerprint mismatch')
    _require(contract.get('schema') == SCHEMA, 'Wrong score contract schema')
    for name in ('contract_id', 'objective_name'):
        _text(contract.get(name), name)
    projection = contract.get('cps_projection')
    _require(projection in ('uniform_birth_time', 'constant_post_cell'),
             'Choose exactly one explicit CPS projection; no implicit averaging')
    targets = _rows(contract.get('target_rows'), 'restriction_id', 'target')
    _require(len(targets) - (NORMALIZATION_ID in targets) >= len(PARAMETERS),
             'Fewer scored rows than free parameters')
    _require(set(targets) == set(MOMENT_PATHS) | {NORMALIZATION_ID},
             'All twelve named moments plus separate normalization are required')
    for name, row in targets.items():
        for field in PROVENANCE_FIELDS:
            _text(row.get(field), name + '.' + field)
        _number(row.get('target'), name + ' target')
    normrow = targets[NORMALIZATION_ID]
    _require(normrow['target'] == 2.1 and normrow.get('actual_weight') is None
             and normrow.get('role') == 'normalization_separate_from_scored_objective'
             and normrow.get('model_observation') == 'final.normalization.completed_fertility',
             'The 2.1 normalization must remain separate and unweighted')
    _require(targets['recent_parent_ownership']['target'] == RECENT_TARGET,
             'Preserve the authoritative ACS recent-parent target')
    sources = contract.get('source_fingerprints')
    _require(isinstance(sources, Mapping) and sources, 'Source fingerprint map required')
    for name, pin in sources.items():
        _text(name, 'source fingerprint name'); _hash(pin, name)
    receipt = evaluation_receipt
    _require(receipt.get('schema') == 'e5f_initial_score_receipt_v1'
             and receipt.get('status') == 'verified'
             and receipt.get('numerical_gates_verified') is True
             and receipt.get('recent_parent_certified') is True,
             'Verified numerical and recent-parent certification receipt required')
    _require(receipt.get('source_fingerprints') == sources, 'Source fingerprint mismatch')
    checkpoint = _hash(receipt.get('checkpoint_sha256'), 'checkpoint_sha256')
    inputs = {'early_measurement': early_measurement, 'recent_parent_observation': recent_parent_observation,
              'normalization': normalization, 'parameters': parameters}
    _require(receipt.get('input_sha256') == {k: fingerprint(v) for k, v in inputs.items()},
             'Evaluation input fingerprint mismatch')
    approximation = contract.get('recent_parent_approximation', {})
    _require(isinstance(approximation, Mapping), 'Recent-parent approximation declaration required')
    _require(receipt.get('recent_parent_approximation_id') == approximation.get('approximation_id'),
             'Recent-parent certification names a different approximation')
    recent = _validate_recent(recent_parent_observation, approximation, checkpoint)
    _require(early_measurement.get('weights_assigned') is False
             and early_measurement.get('empirical_target_contract_activated') is False
             and early_measurement.get('calibrated_smm') is False,
             'Expected unchanged early measurement packet')
    fertility = early_measurement.get('fertility', {})
    _require(projection in fertility and fertility[projection].get('metadata', {}).get('age_projection') == projection,
             'CPS observer projection metadata mismatch')
    tolerance = _number(contract.get('normalization_tolerance'), 'normalization_tolerance')
    _require(0 < tolerance <= 5e-4, 'Invalid normalization tolerance')
    normtarget = _number(normalization.get('target'), 'normalization target')
    normvalue = _number(normalization.get('completed_fertility'), 'normalized completed fertility')
    gap = normvalue - 2.1
    _require(normtarget == 2.1 and abs(gap) <= tolerance, 'Initial normalization failed')
    _require(math.isclose(_number(normalization.get('absolute_gap'), 'normalization absolute_gap'),
                         abs(gap), rel_tol=0, abs_tol=1e-12), 'Normalization gap receipt mismatch')
    near = _number(contract.get('near_bound_fraction'), 'near_bound_fraction')
    _require(0 <= near <= .5, 'Invalid near-bound fraction')
    parameter_table = _parameter_table(parameters, contract.get('parameter_restrictions'), near)
    normalized_parameter = next((x for x in parameter_table if x['parameter'] == 'psi_child'), None)
    if normalized_parameter is not None:
        psi = _number(normalization.get('psi_child'), 'normalization psi_child')
        _require(math.isclose(normalized_parameter['estimate'], psi, rel_tol=0, abs_tol=1e-12),
                 'Normalized preference differs from parameter receipt')
    fits = []
    for original in contract['target_rows']:
        row = copy.deepcopy(original); name = row['restriction_id']
        # Do not carry the baseline CSV's stale model availability or file path
        # into the current evaluation. Empirical provenance remains untouched.
        row.update(model_available=True, model_checkpoint_sha256=checkpoint, model_source_path=None, calibrated_smm=False)
        if name == NORMALIZATION_ID:
            row.update(model=normvalue, gap=gap, actual_weight=None, loss_contribution=None, scored=False,
                       evaluation_input='normalization', model_source_json_location='completed_fertility')
        else:
            expected_path = MOMENT_PATHS[name].format(projection=projection)
            _require(row.get('model_observation') == expected_path, 'Wrong model observation mapping: ' + name)
            _require(row.get('role') == 'proposed_scored_restriction', 'Wrong scored role: ' + name)
            weight = _number(row.get('actual_weight'), name + ' actual_weight')
            _require(weight > 0, 'Working weight must be strictly positive: ' + name)
            value = recent if name == 'recent_parent_ownership' else _lookup(early_measurement, expected_path)
            difference = value - row['target']
            try:
                contribution = weight * difference ** 2
            except OverflowError as exc:
                raise ScoreContractError('Nonfinite score contribution: ' + name) from exc
            _number(contribution, name + ' loss contribution')
            row.update(model=value, gap=difference, actual_weight=weight,
                       loss_contribution=contribution, scored=True,
                       evaluation_input='recent_parent_observation' if name == 'recent_parent_ownership' else 'early_measurement',
                       model_source_json_location='model_value' if name == 'recent_parent_ownership' else expected_path)
            if name == 'recent_parent_ownership':
                row['maintained_model_approximation'] = copy.deepcopy(approximation)
        fits.append(row)
    try:
        loss = math.fsum(row['loss_contribution'] for row in fits if row['scored'])
    except OverflowError as exc:
        raise ScoreContractError('Nonfinite total score') from exc
    _number(loss, 'total score')
    return {'schema': 'e5f_initial_minimum_distance_result_v1', 'objective_name': contract['objective_name'],
            'contract_id': contract['contract_id'], 'contract_sha256': actual_hash,
            'evaluation_receipt_sha256': fingerprint(receipt), 'checkpoint_sha256': checkpoint,
            'source_fingerprints': copy.deepcopy(sources), 'cps_projection': projection,
            'loss': loss, 'scored_moment_count': 12, 'restriction_count': 13,
            'free_parameter_count': 9, 'target_fit': fits, 'parameters': parameter_table,
            'normalization': copy.deepcopy(normalization), 'calibrated_smm': False,
            'benchmark_certified': False, 'identification_established': False,
            'identification_note': 'Moment count alone is not an identification proof.',
            'weight_interpretation': 'Supplied frozen working weights; no optimal-weight claim.'}
