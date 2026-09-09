"""Unchanged twelve-row measurement on a supplied historical PF path.

This observer does not solve a model or approve the empirical parent/date
mapping. It retains the 2019--2023 matched branch and old-cohort prehistory.
"""
from __future__ import annotations

import copy
import math
from types import SimpleNamespace
from typing import Any

import numpy as np

import run_e5f_transition_calibration as measurement


PARAMETER_NAMES = (
    'beta_annual', 'kappa_fert', 'kappa_fert_continuation', 'chi', 'H0',
    'theta0', 'theta1', 'hbar_child_rooms', 'first_birth_fixed_cost',
    'hbar_first_child_jump', 'psi_child_change_2023',
)
FIXED_FIELDS = ('J', 'age_start', 'da', 'period_years', 'n_parity',
    'n_child_states', 'tfr_top_bin_weight', 'tenure_choice_kappa',
    'kappa_fert', 'kappa_fert_continuation', 'beta', 'chi', 'H0', 'theta0',
    'theta1', 'hbar_child_rooms', 'first_birth_fixed_cost', 'hbar_first_child_jump')
CHOICE_FLAGS = ('joint_nested_choice', 'fertility_nest_choice',
                'two_shock_choice', 'exhaustive_saving_control')


def _parameter_snapshot(parameters):
    # Measurement needs coordinates and age geometry, not mutable Bellman
    # caches on P; avoid retaining or copying full joint probability tensors.
    values = {name: copy.deepcopy(getattr(parameters, name))
              for name in (*FIXED_FIELDS, 'psi_child')}
    values.update({name: bool(getattr(parameters, name, False)) for name in CHOICE_FLAGS})
    return SimpleNamespace(**values)


class HistoricalMomentObserver:
    """Observe global periods 0--4 (2007--2023); later tail dates are ignored.

    old_first_birth_accounting must come from the supplied normalized old
    stationary evaluation, before reweighting its age margins to 2007.
    old_normalization_tolerance is explicit; this class cannot normalize a seed.
    """

    def __init__(self, *, target_system, expected_target_fingerprint: str,
                 old_parameters, old_first_birth_accounting: dict[str, Any],
                 old_normalization: dict[str, Any], old_normalization_tolerance: float):
        canonical = measurement.e5_target_system_for_profile('baseline')
        if (target_system.count != 12 or target_system.name != canonical.name
                or target_system.fingerprint != expected_target_fingerprint
                or target_system.fingerprint != canonical.fingerprint):
            raise ValueError('Expected the pinned unchanged baseline twelve-target system')
        target_system.require_identified(11)
        if (not math.isfinite(old_normalization_tolerance) or not 0 < old_normalization_tolerance <= 5e-4
                or float(old_normalization.get('target', math.nan)) != 2.1
                or old_normalization.get('status') not in ('derived_intercept', 'normalized_at_initial_guess')
                or not math.isfinite(float(old_normalization.get('completed_fertility', math.nan)))
                or abs(float(old_normalization['completed_fertility']) - 2.1) > old_normalization_tolerance
                or not math.isclose(float(old_normalization.get('psi_child', math.nan)),
                                    float(old_parameters.psi_child), rel_tol=0, abs_tol=1e-12)):
            raise ValueError('Old stationary 2.1 normalization must be supplied and verified')
        if float(old_parameters.period_years) != 4. or int(old_parameters.J) != 17:
            raise ValueError('Historical target adapter requires the existing four-year, 17-age model')
        self.target_system = target_system
        self.fingerprint = expected_target_fingerprint
        self.old_parameters = _parameter_snapshot(old_parameters)
        self.normalization = copy.deepcopy(old_normalization)
        self.normalization_tolerance = float(old_normalization_tolerance)
        self.old_accounting = {}
        for key in ('at_risk', 'flow', 'hazard'):
            values = np.asarray(old_first_birth_accounting[key], dtype=float)
            if values.shape != (17,) or not np.isfinite(values).all() or np.any(values < 0):
                raise ValueError(f'Invalid old first-birth accounting: {key}')
            self.old_accounting[key] = values.copy()
        implied = np.divide(self.old_accounting['flow'], self.old_accounting['at_risk'],
            out=np.zeros(17), where=self.old_accounting['at_risk'] > 1e-15)
        if (np.any(self.old_accounting['hazard'] > 1 + 1e-12)
                or not np.allclose(implied, self.old_accounting['hazard'], rtol=0, atol=2e-12)):
            raise ValueError('Old first-birth flow/risk/hazard accounting is inconsistent')
        self.period_records = {}
        self.pending_branch = None
        self.moments = None
        self.timing = None
        self.housing_measurement = None
        self.parameters_2023 = None
        self.last_period = -1

    def __call__(self, period, evaluation, parameters, grid, shared):
        if type(period) is not int or period != self.last_period + 1:
            raise ValueError('Historical moment dates must arrive once in chronological order from period zero')
        if self.target_system.fingerprint != self.fingerprint:
            raise RuntimeError('Target fingerprint changed during observation')
        if period > 4:
            if self.moments is None:
                raise RuntimeError('Person tail arrived before the 2023 measurement')
            self.last_period = period
            return
        for name in FIXED_FIELDS:
            if not np.array_equal(np.asarray(getattr(parameters, name)),
                                  np.asarray(getattr(self.old_parameters, name))):
                raise ValueError(f'Historical measurement parameter changed across dates: {name}')
        for name in CHOICE_FLAGS:
            if bool(getattr(parameters, name, False)) != bool(getattr(self.old_parameters, name, False)):
                raise ValueError(f'Historical choice arm changed: {name}')
        self.period_records[period] = {'first_birth_accounting_by_age':
            measurement.first_birth_accounting_by_age(evaluation, parameters)}
        if period == 3:
            self.pending_branch = measurement.begin_dated_first_birth_housing_branch(
                evaluation, parameters, grid, shared, origin_period=3)
        if period == 4:
            if self.pending_branch is None:
                raise RuntimeError('Missing the 2019 first-birth matched branch')
            housing = measurement.finish_dated_first_birth_housing_branch(
                self.pending_branch, evaluation, parameters, grid, shared, destination_period=4)
            if housing.get('census_age_bridge_applied') is not False:
                raise RuntimeError('The fixed birth/control branch must not receive an aggregate age bridge')
            moments = measurement.transition_cross_section_moments(
                evaluation, parameters, grid, shared, self.target_system.moment_names,
                housing_increment_override=float(housing['housing_response']))
            timing = measurement.cohort_timing_moments(
                self.old_accounting, self.period_records, self.old_parameters,
                terminal_period=4, terminal_childless_probability=float(moments['childless_rate']))
            moments['mean_age_first_birth'] = float(timing['mean_age_first_birth'])
            moments['share_first_births_age30plus'] = float(timing['share_first_births_age30plus'])
            if abs(float(timing['synthetic_childless_probability']) - float(moments['childless_rate'])) > 2e-10:
                raise RuntimeError('Synthetic first-birth cohort misses the 2023 childless stock')
            if set(moments) != set(self.target_system.moment_names) or not all(math.isfinite(float(v)) for v in moments.values()):
                raise RuntimeError('The complete twelve-moment vector is missing or nonfinite')
            self.moments, self.timing, self.housing_measurement = moments, timing, housing
            self.parameters_2023 = _parameter_snapshot(parameters)
            self.pending_branch = None
        self.last_period = period

    def result(self, candidate: str) -> dict[str, Any]:
        """Return all twelve rows and the existing weighted-square loss."""
        if self.moments is None or sorted(self.period_records) != list(range(5)):
            raise RuntimeError('Historical moment measurement is incomplete')
        if self.target_system.fingerprint != self.fingerprint:
            raise RuntimeError('Target fingerprint changed before reporting')
        fit_rows, loss = measurement.target_fit_rows(candidate, self.moments, self.target_system)
        return dict(status='complete_conditional_pf_moments_diagnostic_only',
            target_set=self.target_system.name, target_fingerprint=self.fingerprint,
            target_count=12, measurement_year=2023, moments=copy.deepcopy(self.moments),
            target_fit_rows=fit_rows, loss=loss, timing=copy.deepcopy(self.timing),
            first_birth_housing_did=copy.deepcopy(self.housing_measurement),
            old_normalization=copy.deepcopy(self.normalization),
            old_normalization_tolerance=self.normalization_tolerance,
            equilibrium_verified=False, production_promoted=False,
            outstanding=['ownership parent/control model-data alignment',
                         'empirical date-window alignment',
                         'market and terminal equilibrium certification'],
            measurement_scope='unchanged targets and groups on supplied PF choices; no target promotion')

    def report(self, candidate: str, *, theta: dict[str, float], parameter_domain,
               supply_rule) -> dict[str, Any]:
        """Add the existing parameter report under an explicitly configured domain.

        calibration.parameter_rows reads its module's configured domain. Fail
        rather than silently falling back to that module's older import default.
        """
        result = self.result(candidate)
        domain = tuple(tuple(row) for row in parameter_domain)
        if (tuple(row[0] for row in domain) != PARAMETER_NAMES
                or domain != tuple(measurement.TRANSITION_SEARCH_DOMAIN)):
            raise ValueError('Explicit eleven-parameter domain must match the configured parameter_rows helper')
        for name, lower, upper, transform in domain:
            if not (math.isfinite(lower) and math.isfinite(upper) and lower < upper):
                raise ValueError(f'Invalid parameter bounds: {name}')
            if name == 'psi_child_change_2023':
                continue
            key = 'beta' if name == 'beta_annual' else name
            if not np.allclose(theta[key], getattr(self.parameters_2023, key), rtol=0, atol=1e-12):
                raise ValueError(f'Reported parameter differs from observed PF path: {key}')
        if float(supply_rule.elasticity) != .63 or float(self.parameters_2023.tenure_choice_kappa) != .005:
            raise ValueError('Unexpected external supply elasticity or tenure scale')
        rows = measurement.parameter_rows(theta, {'new_psi_child': self.parameters_2023.psi_child},
            old_psi_child=self.old_parameters.psi_child, joint_panel=False)
        for row in rows:
            if row['parameter'] == 'psi_child_2007':
                row['value'] = float(self.old_parameters.psi_child)
            elif row['parameter'] == 'psi_child_2023':
                row['value'] = float(self.parameters_2023.psi_child)
        for name, value in (('tenure_choice_kappa', .005), ('housing_supply_elasticity', .63)):
            rows.append(dict(parameter=name, value=value, lower_bound=math.nan, upper_bound=math.nan,
                transform='external_profile', is_free_parameter=False,
                status='externally_fixed_profile_not_estimated', near_bound=False))
        result['parameter_rows'] = rows
        return result
