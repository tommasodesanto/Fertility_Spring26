"""Default-off mechanism switches for the stationary-state sandbox.

These switches are NOT implemented by editing the production package
(`code/model/intergen_eqscale_seq_optimized/`). They are pure-Python runtime
overrides that live entirely in this file: at solve time, `sandbox_context()`
temporarily monkeypatches `solver.precompute_shared` with a wrapper defined
here, and restores the original function object when the `with` block exits.
No package file is ever written to.

Two of the switches described in the sandbox task have a clean runtime
hook and are implemented below:

  (a) child_benefit_form in {"linear", "log", "power"} with
      child_benefit_curvature -- replaces psi*m by psi*log(1+m) or
      psi*((1+m)**(1-eps)-1)/(1-eps) in the child-benefit term that production
      builds at solver.py:2254 (`psi_v[nn, cs] = P.psi_child * nk`).

  (b) scale_weighting in {"deflate", "multiply"} -- changes the household
      equivalence-scale exponent applied at solver.py:2278 from sigma-1 to
      sigma, for the "power"/"sqrt" eqscale_form branches (the branches that
      actually apply a sigma-dependent exponent in the live code). The
      "linear"/gamma_e branch (solver.py:2284, `escale = 1 + gamma_e*nk`)
      applies no sigma exponent at all in production, at any sigma; there is
      therefore no well-defined "multiply" analogue for it without inventing
      a new functional form the author has not specified. When
      eqscale_form == "linear", scale_weighting == "multiply" is a no-op
      (deflate and multiply coincide) and this is logged in summary.md.

A third former switch, (c) child_earnings_penalty, is NOT implemented here
because it needs no sandbox hook at all: the production package implements
it directly (parameters.py:393-404, `apply_overrides` coerces a scalar or
4-vector into P.child_earnings_penalty with one entry per m = 0, 1, 2, 3+).
Like the other native package switches (mortgage_origination_only,
mortgage_amortization, rental_wedge_intercept/slope/knee, estate_receiver,
bequest_net_of_selling_cost, child_maturation_mode, mu_young, a_rise,
a_full), it passes through the sandbox's override path untouched into the
package's own apply_overrides (see _SANDBOX_ONLY_OVERRIDE_KEYS below, which
lists the only keys stripped out before that call).

Bitwise-off guarantee: when a spec leaves child_benefit_form == "linear" and
scale_weighting == "deflate" (the defaults), `sandboxed_precompute_shared`
returns the ORIGINAL package function's output completely unmodified -- it
never rebuilds or mutates the arrays. This makes the "off" path trivially
bitwise-identical to unpatched production, which is what
tests/test_mechanisms.py checks.
"""
from __future__ import annotations

import contextlib
from types import SimpleNamespace
from typing import Any, Iterator

import numpy as np

from intergen_eqscale_seq_optimized import solver as _solver

# Captured once, at import time, before sandbox_context() ever patches
# _solver.precompute_shared -- so sandboxed_precompute_shared always calls
# the real package function, never itself, no matter how many times
# sandbox_context() is entered/exited or re-entered.
_ORIGINAL_PRECOMPUTE_SHARED = _solver.precompute_shared

# Formerly: child_earnings_penalty was listed here as NOT_IMPLEMENTED on the
# theory that no sandbox-only hook existed (income_at_state takes no
# child-state argument). That theory was wrong: the production package
# implements the penalty natively (parameters.py:393-404), so no sandbox hook
# is needed -- the override passes straight through to the package. The
# mapping is kept (empty) so that check_switches_supported() remains a valid
# call for existing callers.
NOT_IMPLEMENTED: dict[str, str] = {}


def _child_benefit_value(psi_child: float, nk: np.ndarray, form: str, curvature: float) -> np.ndarray:
    """psi*f(m) for f in {identity, log(1+m), CRRA-style power of (1+m)}."""
    if form == "linear":
        return psi_child * nk
    if form == "log":
        return psi_child * np.log1p(nk)
    if form == "power":
        eps = float(curvature)
        if eps == 1.0:
            return psi_child * np.log1p(nk)
        return psi_child * (np.power(1.0 + nk, 1.0 - eps) - 1.0) / (1.0 - eps)
    raise ValueError(f"Unknown child_benefit_form: {form!r} (expected linear, log, power)")


def _nk_kp_grids(P: SimpleNamespace) -> tuple[np.ndarray, np.ndarray]:
    """Reproduce the (nn, cs) -> (nk, kp) mapping at solver.py:2235-2247 exactly.

    nk is "number of children currently at home" (or "children ever born" in
    shared-clock mode) and kp marks which (nn, cs) cells are child-having
    cells (as opposed to childless baseline cells). This is read-only
    bookkeeping duplicated from the package's own predicates
    (independent_child_maturation_active, readiness_gate_active), imported
    rather than reimplemented, so it tracks the package if those predicates'
    *inputs* change; only the small loop shape is copied.
    """
    K = P.n_child_stages
    csm1 = K + 1
    nk = np.zeros((P.n_parity, P.n_child_states), dtype=float)
    kp = np.zeros((P.n_parity, P.n_child_states), dtype=bool)
    independent = _solver.independent_child_maturation_active(P)
    readiness = _solver.readiness_gate_active(P)
    for nn in range(P.n_parity):
        for cs in range(P.n_child_states):
            if independent:
                this_nk = cs if cs <= nn else 0
                this_kp = this_nk > 0
            else:
                this_nk = nn
                this_kp = (cs >= 1) and (cs < csm1)
            if readiness and nn == 0 and cs == 1:
                this_kp = False
            nk[nn, cs] = this_nk
            kp[nn, cs] = this_kp
    return nk, kp


def sandboxed_precompute_shared(P: SimpleNamespace, b_grid: np.ndarray) -> SimpleNamespace:
    """Drop-in replacement for solver.precompute_shared with two extra switches.

    When both switches are at their default ("linear", "deflate") this calls
    straight through to the original function and returns its result
    untouched -- see the module docstring for why that guarantees bitwise
    identity at the off value.
    """
    form = str(getattr(P, "child_benefit_form", "linear")).lower()
    weighting = str(getattr(P, "scale_weighting", "deflate")).lower()
    if form == "linear" and weighting == "deflate":
        return _ORIGINAL_PRECOMPUTE_SHARED(P, b_grid)

    if form not in {"linear", "log", "power"}:
        raise ValueError(f"Unknown child_benefit_form: {form!r}")
    if weighting not in {"deflate", "multiply"}:
        raise ValueError(f"Unknown scale_weighting: {weighting!r}")

    SD = _ORIGINAL_PRECOMPUTE_SHARED(P, b_grid)
    nk, kp = _nk_kp_grids(P)
    curvature = float(getattr(P, "child_benefit_curvature", 0.5))

    psi_v = SD.psi_v.copy()
    if form != "linear":
        new_values = _child_benefit_value(float(P.psi_child), nk, form, curvature)
        psi_v = np.where(kp, new_values, 0.0)

    escale = SD.escale_flat.reshape(P.n_parity, P.n_child_states).copy()
    eqscale_form = str(getattr(P, "eqscale_form", "linear")).lower()
    if weighting == "multiply" and eqscale_form in {"power", "sqrt"}:
        sigma = float(P.sigma)
        if eqscale_form == "power":
            base = np.power((2.0 + 0.7 * nk) / 2.0, 0.7)
        else:  # sqrt
            base = np.power((2.0 + nk) / 2.0, 0.5)
        multiply_escale = np.power(base, sigma)
        escale = np.where(kp, multiply_escale, escale)
    # else: multiply has no effect for eqscale_form == "linear" (no sigma
    # exponent exists there in production to shift from sigma-1 to sigma).

    if form == "linear" and np.array_equal(psi_v, SD.psi_v) and np.array_equal(escale, SD.escale_flat.reshape(P.n_parity, P.n_child_states)):
        return SD  # nothing actually changed (e.g. multiply requested but eqscale_form == linear)

    nc = SD.nc
    c_bar = SD.c_bar
    h_bar = SD.h_bar
    g_bar = SD.g_bar
    alpha_bar = SD.alpha_flat.reshape(P.n_parity, P.n_child_states)
    triples = np.column_stack([
        c_bar.reshape(-1, order="F"),
        h_bar.reshape(-1, order="F"),
        psi_v.reshape(-1, order="F"),
    ])
    unique_triples, type_map = np.unique(triples, axis=0, return_inverse=True)

    return SimpleNamespace(
        c_bar=c_bar,
        h_bar=h_bar,
        psi_v=psi_v,
        g_bar=g_bar,
        cb_flat=c_bar.reshape(1, nc, order="F"),
        hb_flat=h_bar.reshape(1, nc, order="F"),
        psi_flat=psi_v.reshape(1, nc, order="F"),
        gb_flat=g_bar.reshape(1, nc, order="F"),
        alpha_flat=alpha_bar.reshape(1, nc, order="F"),
        escale_flat=escale.reshape(1, nc, order="F"),
        nc=nc,
        b=SD.b,
        bp=SD.bp,
        phi_state=SD.phi_state,
        phi_choice=SD.phi_choice,
        n_types=unique_triples.shape[0],
        type_map=type_map,
        type_cb=unique_triples[:, 0],
        type_hb=unique_triples[:, 1],
        type_psi=unique_triples[:, 2],
        birth_dp=SD.birth_dp,
        birth_entry_grant=SD.birth_entry_grant,
    )


def check_switches_supported(overrides: dict[str, Any]) -> None:
    """Raise NotImplementedError early for any switch this sandbox can't honor."""
    for name, reason in NOT_IMPLEMENTED.items():
        value = overrides.get(name)
        if value not in (None, 0, 0.0, False, {}):
            raise NotImplementedError(f"sandbox switch '{name}' is not implemented: {reason}")


# Sandbox-invented override keys with no corresponding attribute on the
# production parameter object. parameters.apply_overrides (imported into
# solver.py's namespace as `apply_overrides`, parameters.py:313) rejects any
# override key not already in `vars(P)` or its SUPPORTED_DYNAMIC_OVERRIDE_KEYS
# allowlist, so these three keys must be stripped out of the override dict
# before it reaches the real apply_overrides and set as plain attributes on
# the P object it returns instead -- otherwise every spec that uses a
# mechanism switch would fail with "Unknown parameter override(s)".
_SANDBOX_ONLY_OVERRIDE_KEYS = ("child_benefit_form", "child_benefit_curvature", "scale_weighting")
_ORIGINAL_APPLY_OVERRIDES = _solver.apply_overrides


def _sandboxed_apply_overrides(P: SimpleNamespace, overrides: Any) -> SimpleNamespace:
    if isinstance(overrides, dict):
        stripped = {k: v for k, v in overrides.items() if k not in _SANDBOX_ONLY_OVERRIDE_KEYS}
        extra = {k: v for k, v in overrides.items() if k in _SANDBOX_ONLY_OVERRIDE_KEYS}
    else:
        stripped, extra = overrides, {}
    P = _ORIGINAL_APPLY_OVERRIDES(P, stripped)
    for key, value in extra.items():
        setattr(P, key, value)
    return P


@contextlib.contextmanager
def sandbox_context() -> Iterator[None]:
    """Monkeypatch solver.precompute_shared and solver.apply_overrides for one solve.

    Both patches are process-local and always reverted in `finally`, including
    on exception. Neither ever touches any file under
    code/model/intergen_eqscale_seq_optimized/.
    """
    original = _solver.precompute_shared
    original_apply = _solver.apply_overrides
    _solver.precompute_shared = sandboxed_precompute_shared
    _solver.apply_overrides = _sandboxed_apply_overrides
    try:
        yield
    finally:
        _solver.precompute_shared = original
        _solver.apply_overrides = original_apply
