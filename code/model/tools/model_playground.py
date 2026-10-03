"""Interactive fixed-price access through the canonical production engine.

The historical authenticated loader remains available under an explicit legacy
name. Importing this module loads no model state and performs no lifecycle solve.
"""
from __future__ import annotations

import copy
import hashlib
import importlib
import json
import os
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

for _thread_var in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_thread_var] = "1"


ROOT = Path(__file__).resolve().parents[3]
EXPLORER_CASES = ROOT / "output/model/fixed_reference_economics_20260928/soft_timing_review_v1/explorer_cases.json"
SELECTION = ROOT / "output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_selected.json"
REFERENCE_CLOSURE = ROOT / "output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_postcheck/selected_postcheck/phase_b_ge/selected_repeat_final/closure.json"
V2_DIR = ROOT / "output/model/fixed_reference_economics_20260928/normalized_calibration_v2"
ORIGINAL_ENGINE = ROOT / "output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed/source"

PARAMETER_ORDER = (
    "beta_annual", "chi", "first_birth_fixed_cost", "kappa_fert",
    "kappa_fert_continuation", "theta0", "h_P", "child_benefit_curvature",
    "tenure_choice_kappa", "psi_child",
)
PARAMETER_DESCRIPTIONS = {
    "beta_annual": "Annual discount factor; native period beta is beta_annual ** period_years.",
    "chi": "Owner housing-service premium.",
    "first_birth_fixed_cost": "Fixed utility cost of the first birth.",
    "kappa_fert": "First-birth choice shock/logit scale.",
    "kappa_fert_continuation": "Later-birth attempt choice shock/logit scale.",
    "theta0": "Bequest utility scale.",
    "h_P": "Physical room floor added at the first child.",
    "child_benefit_curvature": "Curvature of the child benefit by children at home.",
    "tenure_choice_kappa": "Tenure-choice logit scale.",
    "psi_child": "Child benefit scale.",
}
_PARAMETER_NATIVE_FIELDS = {
    "beta_annual": {"beta", "rho", "rho_hat"},
    "chi": {"chi"},
    "first_birth_fixed_cost": {"first_birth_fixed_cost"},
    "kappa_fert": {"kappa_fert", "eps_fert"},
    "kappa_fert_continuation": {"kappa_fert_continuation"},
    "theta0": {"theta0"},
    "h_P": {"hbar_first_child_jump"},
    "child_benefit_curvature": {"child_benefit_curvature"},
    "tenure_choice_kappa": {"tenure_choice_kappa"},
    "psi_child": {"psi_child"},
}
_DISPLAYED_NATIVE_FIELDS = set().union(*_PARAMETER_NATIVE_FIELDS.values())
_STRUCTURAL_FIELDS = {
    "Nb", "Nz", "wealth_grid_nodes", "income_states", "b_min", "b_max",
    "b_grid_power", "earnings_transaction_grid", "period_years", "J_R",
    "z_grid", "z_weights", "Pi_z", "entry_wealth_ratio_nodes",
    "entry_wealth_ratio_weights", "fixed_reference_entry_conditional",
    "fixed_reference_entry_grid",
}
_RECOMPUTED_FIELDS = {"q", "user_cost_rate", "rho", "rho_hat", "eps_fert",
                      "pension", "pension_by_loc", "income"}


def _same_value(left, right) -> bool:
    import numpy as np
    try:
        return bool(np.array_equal(np.asarray(left), np.asarray(right), equal_nan=True))
    except TypeError:
        return bool(np.array_equal(np.asarray(left), np.asarray(right)))


def _policy_tools():
    try:
        from . import model_policy_tools
    except ImportError:
        import model_policy_tools
    return model_policy_tools


def bind_parameters(P: SimpleNamespace, parameters: dict[str, float]) -> SimpleNamespace:
    """Copy P and apply the ten native calibration coordinates without clamping."""
    import numpy as np
    if set(parameters) != set(PARAMETER_ORDER):
        missing = sorted(set(PARAMETER_ORDER) - set(parameters))
        extra = sorted(set(parameters) - set(PARAMETER_ORDER))
        raise ValueError(f"Expected the ten selected parameters; missing={missing}, extra={extra}")
    values = {name: float(parameters[name]) for name in PARAMETER_ORDER}
    if not np.isfinite(list(values.values())).all():
        raise ValueError("Parameter values must be finite")
    if not 0.0 < values["beta_annual"] < 1.0:
        raise ValueError("beta_annual must lie strictly between zero and one")
    if values["chi"] <= 0.0 or values["h_P"] < 0.0 or values["tenure_choice_kappa"] < 0.0:
        raise ValueError("chi must be positive; h_P and tenure_choice_kappa must be nonnegative")
    if not 0.0 <= values["child_benefit_curvature"] < 1.0:
        raise ValueError("child_benefit_curvature must be in [0, 1)")
    Q = copy.deepcopy(P)
    for name, value in values.items():
        if name not in {"beta_annual", "h_P"}:
            setattr(Q, name, value)
    # Native calibration mapping: annual discounting is compounded to the
    # model period and the associated discount-rate fields move with it.
    Q.beta = values["beta_annual"] ** float(Q.period_years)
    Q.rho = 1.0 / Q.beta - 1.0
    Q.rho_hat = Q.rho
    Q.eps_fert = values["kappa_fert"]
    # h_P is the physical first-child room floor; later-child floors stay zero.
    Q.child_room_floor = True
    Q.hbar_first_child_jump = values["h_P"]
    Q.hbar_child_rooms = 0.0
    return Q


def _selected_parameters() -> dict[str, float]:
    selected = json.loads(SELECTION.read_text())["selected"]
    point = selected["parameters"]
    if tuple(point.keys()) != PARAMETER_ORDER and set(point) != set(PARAMETER_ORDER):
        raise RuntimeError("Selected soft record does not contain the expected ten parameters")
    return {name: float(point[name]) for name in PARAMETER_ORDER}


def load_saved_solution(case: str = "soft") -> SimpleNamespace:
    """Load a saved explorer case after checking its recorded SHA-256; no solve."""
    import numpy as np
    config = json.loads(EXPLORER_CASES.read_text())
    spec = next((item for item in config["cases"] if item["id"] == case), None)
    if spec is None:
        raise KeyError(f"Unknown saved case {case!r}; available: {[x['id'] for x in config['cases']]}")
    path = Path(spec["arrays"])
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != spec["sha256"]:
        raise RuntimeError(f"Saved array fingerprint differs for {case}: {path}")
    with np.load(path, allow_pickle=False) as data:
        values = {key: data[key].copy() for key in data.files}
    values.update(price=float(spec["price"]), timing=spec["timing"], case_id=case)
    return SimpleNamespace(**values)


def _install_read_only_overlay() -> None:
    """Install the existing hash-checked, read-only local source mapping once."""
    marker = "_fertility_model_playground_overlay"
    overlay = ROOT / "output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/local_runtime/bootstrap.py"
    installed = getattr(sys, marker, None)
    if installed:
        if installed != str(overlay):
            raise RuntimeError("A different frozen-source overlay is already installed")
        return
    source = overlay.read_text()
    prefix = source.split("if '--preflight-context' in sys.argv:", 1)[0]
    if "DIGESTS=" not in prefix or "Frozen overlay write forbidden" not in prefix:
        raise RuntimeError("The authenticated read-only overlay has an unexpected shape")
    scope = {"__file__": str(overlay), "__name__": "authenticated_local_overlay"}
    exec(compile(prefix, str(overlay), "exec"), scope)
    setattr(sys, marker, str(overlay))


def load_reference_model() -> tuple[SimpleNamespace, np.ndarray, SimpleNamespace]:
    """Build the authenticated soft reference inputs and return ``P, b_grid, solver``.

    This only constructs inputs and shared model objects. It does not solve.
    The returned price is in ``P.reference_price``; the closure's derived H0 is
    retained in ``P.H0``.
    """
    import numpy as np
    _install_read_only_overlay()
    sys.path.insert(0, str(V2_DIR))
    import run_psi as v2
    import numba
    numba.set_num_threads(1)

    v2.native.verify_sources()
    selected = json.loads(SELECTION.read_text())
    source = ROOT / selected["source"]
    if hashlib.sha256(source.read_bytes()).hexdigest() != selected["source_sha256"]:
        raise RuntimeError("Selected soft parameter source hash differs")
    source_record = json.loads(source.read_text())
    if source_record[selected["source_key"]]["best"] != selected["selected"]:
        raise RuntimeError("Selected soft parameter record differs from its source")
    point = selected["selected"]["parameters"]
    if set(point) != set(PARAMETER_ORDER):
        raise RuntimeError("Selected soft record does not contain the expected ten parameters")
    if selected["selected"]["weight_fingerprint"] != v2.weight_fingerprint({}):
        raise RuntimeError("Soft reference weight fingerprint differs")
    if v2.native.target_identity(selected["selected"]["target_fit"]) != v2.CONFIG["base_target_contract"]:
        raise RuntimeError("Soft reference target contract differs")

    import copy
    original_lane = copy.deepcopy(v2.inputs.LANES["floor_s0"])
    _, bounds, _ = v2.inputs.seed_and_bounds("floor_s0")
    bounds = {name: tuple(value) for name, value in bounds.items()}
    if bounds["h_P"] != (0.1, 2.3):
        raise RuntimeError("Expected base h_P bound changed")
    bounds.update(h_P=(0.1, 2.6), psi_child=tuple(v2.CONFIG["psi_bounds"]))
    try:
        v2.inputs.LANES["floor_s0"].update(seed=dict(point), bounds=bounds, free_coordinates=list(point))
        P, b_grid = v2.inputs.proposal("floor_s0")
        P, _ = v2.inputs.entry(P, b_grid, "nonnegative_mean")
        P = bind_parameters(P, {name: float(point[name]) for name in PARAMETER_ORDER})
        with tempfile.TemporaryDirectory(prefix="model-playground-check-") as temp:
            P = v2.native.utility_checks(P, b_grid, "floor_s0", Path(temp))
    finally:
        v2.inputs.LANES["floor_s0"] = original_lane

    # Match the fixed-price experiment: corrected zero-credit contract.
    sys.path.insert(0, str(ORIGINAL_ENGINE))
    from small_credit_lab import credit
    from small_credit_lab.engine import solver
    expected_package = (ORIGINAL_ENGINE / "small_credit_lab").resolve()
    if Path(credit.__file__).resolve() != expected_package / "credit.py":
        raise RuntimeError(f"Unexpected credit module origin: {credit.__file__}")
    if Path(solver.__file__).resolve() != expected_package / "engine/solver.py":
        raise RuntimeError(f"Unexpected solver module origin: {solver.__file__}")
    credit.bind_engine_credit(P, "corrected", 0.0)
    closure = json.loads(REFERENCE_CLOSURE.read_text())
    if closure["normalized_population"] != 1.0 or P.N_target != 1.0:
        raise RuntimeError("Soft reference population normalization changed")
    P.H0 = np.array([closure["H0_derived"]])
    P.reference_price = float(closure["price"])
    P.native_inherited_distribution_evidence_dir = tempfile.mkdtemp(prefix="model-playground-inherited-state-")
    return P, b_grid, solver


# Explicit name for old oracle and diagnosis callers; ordinary playground use
# goes through CanonicalModelPlayground and the production input/engine API.
load_legacy_authenticated_reference_model = load_reference_model


def solve_at_price(P: SimpleNamespace, b_grid: np.ndarray, solver: SimpleNamespace, price: float) -> SimpleNamespace:
    """Run one explicit one-core stationary solve at ``price``; no market root."""
    import numpy as np
    if not np.isfinite(price) or price <= 0:
        raise ValueError("price must be finite and positive")
    shared = solver.precompute_shared(P, b_grid)
    return solver.solve_markov_income_at_prices(
        np.array([float(price)]), P, b_grid, SD=shared, fast_stats=False
    )


class ModelResult:
    """A native solution with convenient plots and aggregate summaries."""

    def __init__(self, solution, P, parameters, price, *, label="fixed-price experiment"):
        self.solution = solution
        self.P = copy.deepcopy(P)
        self.parameters = dict(parameters)
        self.price = float(price)
        self.label = label

    def __getattr__(self, name):
        return getattr(self.solution, name)

    def aggregates(self):
        import numpy as np
        aggregate_solution = _policy_tools().aggregate_solution
        return aggregate_solution(
            self.solution, houses=np.asarray(self.P.H_own, dtype=float),
            age_start=int(self.P.age_start), period_years=float(self.P.period_years),
        )

    def plot_policy(self, **kwargs):
        import numpy as np
        plot_policy = _policy_tools().plot_policy
        options = dict(houses=np.asarray(self.P.H_own, dtype=float),
                       age_start=int(self.P.age_start),
                       period_years=float(self.P.period_years))
        options.update(kwargs)
        return plot_policy(self.solution, **options)

    def plot_aggregates(self, *, wealth_range="central"):
        import numpy as np
        plot_aggregates = _policy_tools().plot_aggregates
        return plot_aggregates(
            self.solution, houses=np.asarray(self.P.H_own, dtype=float),
            age_start=int(self.P.age_start), period_years=float(self.P.period_years),
            wealth_range=wealth_range,
        )


class AuthenticatedLegacyModelPlayground:
    """Historical authenticated soft reference for explicit oracle callers."""

    def __init__(self):
        P, b_grid, solver = load_reference_model()
        self.P = P
        self.b_grid = b_grid.copy()
        self.solver = solver
        self._reference_params = _selected_parameters()
        self.params = dict(self._reference_params)
        self.reference_price = float(P.reference_price)
        self._authenticated_base_P = copy.deepcopy(P)
        self.last_result = None

    def solve(self, *, price=None, overrides=None):
        """Solve once at a fixed price; this does not solve a market-price root."""
        if overrides:
            unknown = set(overrides) - set(PARAMETER_ORDER)
            if unknown:
                raise KeyError(f"Unknown model parameters: {sorted(unknown)}")
            self.params.update({key: float(value) for key, value in overrides.items()})
        # Keep advanced edits to non-estimated native P fields. The visible
        # ten-parameter dictionary deliberately takes precedence on overlaps.
        P = bind_parameters(self.P, self.params)
        solve_price = self.reference_price if price is None else float(price)
        solution = solve_at_price(P, self.b_grid, self.solver, solve_price)
        result = ModelResult(solution, P, self.params, solve_price)
        self.last_result = result
        return result

    def saved_result(self, case="soft"):
        """Wrap the hash-checked saved solution as the comparison baseline."""
        solution = load_saved_solution(case)
        price = float(solution.price)
        return ModelResult(solution, self._authenticated_base_P, self._reference_params, price,
                           label=f"saved {case} solution")

    def reset_parameters(self):
        """Restore the ten estimated coordinates to the authenticated selection."""
        self.params.clear()
        self.params.update(self._reference_params)

    def reset(self):
        """Restore the full initialized P object and all ten selected values."""
        restored = copy.deepcopy(self._authenticated_base_P)
        vars(self.P).clear()
        vars(self.P).update(vars(restored))
        self.reset_parameters()

    def compare(self, baseline, changed):
        """Return comparable aggregate levels and changed-minus-baseline gaps."""
        left, right = baseline.aggregates(), changed.aggregates()
        fields = tuple(left["overall"])
        overall = {
            key: {"baseline": left["overall"][key], "changed": right["overall"][key],
                  "difference": right["overall"][key] - left["overall"][key]}
            for key in fields
        }
        by_age = []
        for old, new in zip(left["by_age"], right["by_age"]):
            row = {"age": old["age"]}
            for key in fields:
                a, b = old[key], new[key]
                row[key] = {"baseline": a, "changed": b,
                            "difference": None if a is None or b is None else b - a}
            by_age.append(row)
        return {"units": left["units"], "overall": overall, "by_age": by_age,
                "baseline_price": baseline.price, "changed_price": changed.price,
                "changed_parameters": dict(changed.parameters)}

    def show_parameters(self):
        """Print the editable selected primitives and their native interpretation."""
        print("parameter                    value            meaning")
        for name in PARAMETER_ORDER:
            print(f"{name:28s} {self.params[name]:<16.10g} {PARAMETER_DESCRIPTIONS[name]}")


class CanonicalModelPlayground:
    """Editable partial-equilibrium access through the production inputs/engine."""

    def __init__(self, parameter_file=None):
        import numpy as np
        from production import equilibrium
        from production.inputs import load_inputs
        from production.parameter_files import load_parameter_file
        from run_model import PARAMETER_FILE
        self.parameter_file = PARAMETER_FILE if parameter_file is None else parameter_file
        self._preset = load_parameter_file(self.parameter_file)
        self.parameters = copy.deepcopy(self._preset["parameters"])
        self.params = self.parameters
        self.external_inputs = copy.deepcopy(self._preset["external_inputs"])
        self.native_overrides = copy.deepcopy(self._preset["native_overrides"])
        self.P, self.b_grid = load_inputs(self.parameters, self.external_inputs,
                                         self.native_overrides)
        self._initial_P = copy.deepcopy(self.P)
        self._initial_grid = np.asarray(self.b_grid).copy()
        self.solver = equilibrium
        self.reference_price = float(self._preset["price_guess"])
        self.last_result = None

    def _direct_native_edits(self):
        if not _same_value(self.b_grid, self._initial_grid):
            raise ValueError("Wealth-grid edits are unsupported; use the authenticated 120-node grid.")
        initial = vars(self._initial_P)
        current = vars(self.P)
        added = set(current) - set(initial)
        if added:
            raise ValueError(f"Unknown native P fields cannot be edited: {sorted(added)}")
        changed = {name for name, value in current.items()
                   if not _same_value(value, initial[name])}
        structural = changed & _STRUCTURAL_FIELDS
        if structural:
            raise ValueError(f"Structural grid or entry-law edits are unsupported: {sorted(structural)}")
        displayed = changed & _DISPLAYED_NATIVE_FIELDS
        if displayed:
            raise ValueError("Edit displayed parameters through model.params instead of model.P: "
                             + ", ".join(sorted(displayed)))
        derived = changed & _RECOMPUTED_FIELDS
        if derived:
            raise ValueError("These P fields are derived by the production input loader; edit their "
                             "source primitives instead: " + ", ".join(sorted(derived)))
        return {name: copy.deepcopy(current[name]) for name in changed}

    def solve(self, *, price=None, overrides=None, external_inputs=None,
              native_overrides=None):
        """Run one fixed-price solve through the canonical production engine."""
        import numba
        numba.set_num_threads(1)
        from production.inputs import load_inputs
        unknown = set(overrides or {}) - set(PARAMETER_ORDER)
        if unknown:
            raise KeyError(f"Unknown model parameters: {sorted(unknown)}")
        self.parameters.update({key: float(value) for key, value in (overrides or {}).items()})
        parameters = copy.deepcopy(self.parameters)
        external = copy.deepcopy(self.external_inputs)
        external.update(external_inputs or {})
        native = copy.deepcopy(self.native_overrides)
        native.update(native_overrides or {})
        for name, value in self._direct_native_edits().items():
            # A direct edit replaces unchanged preset defaults. Independent
            # dictionary or per-call edits remain explicit conflicting controls.
            for controls, preset, supplied in (
                    (external, self._preset["external_inputs"], external_inputs or {}),
                    (native, self._preset["native_overrides"], native_overrides or {})):
                if name in controls and not _same_value(controls[name], value):
                    unchanged_default = (name in preset
                                         and _same_value(controls[name], preset[name])
                                         and name not in supplied)
                    if not unchanged_default:
                        raise ValueError(f"Conflicting input and direct P edit for {name}")
                    del controls[name]
            native[name] = value
        P, grid = load_inputs(parameters, external, native)
        solve_price = self.reference_price if price is None else float(price)
        outcome = self.solver.solve_at_price(P, grid, solve_price)
        result = ModelResult(outcome["solution"], outcome["P"], parameters,
                             solve_price, label="production fixed-price experiment")
        self.last_result = result
        return result

    def reset_parameters(self):
        self.parameters.clear()
        self.parameters.update(copy.deepcopy(self._preset["parameters"]))

    def reset(self):
        import numpy as np
        from production.inputs import load_inputs
        self.reset_parameters()
        self.external_inputs.clear()
        self.external_inputs.update(copy.deepcopy(self._preset["external_inputs"]))
        self.native_overrides.clear()
        self.native_overrides.update(copy.deepcopy(self._preset["native_overrides"]))
        self.P, self.b_grid = load_inputs(self.parameters, self.external_inputs,
                                        self.native_overrides)
        self._initial_P = copy.deepcopy(self.P)
        self._initial_grid = np.asarray(self.b_grid).copy()
        self.reference_price = float(self._preset["price_guess"])

    def saved_result(self, case="soft"):
        solution = load_saved_solution(case)
        return ModelResult(solution, self.P, self.parameters, solution.price,
                           label=f"historical saved {case} solution")

    def compare(self, baseline, changed):
        left, right = baseline.aggregates(), changed.aggregates()
        fields = tuple(left["overall"])
        overall = {key: {"baseline": left["overall"][key],
                         "changed": right["overall"][key],
                         "difference": right["overall"][key] - left["overall"][key]}
                   for key in fields}
        by_age = []
        for old, new in zip(left["by_age"], right["by_age"]):
            row = {"age": old["age"]}
            for key in fields:
                a, b = old[key], new[key]
                row[key] = {"baseline": a, "changed": b,
                            "difference": None if a is None or b is None else b-a}
            by_age.append(row)
        return {"units": left["units"], "overall": overall, "by_age": by_age,
                "baseline_price": baseline.price, "changed_price": changed.price,
                "changed_parameters": dict(changed.parameters)}

    def show_parameters(self):
        print("Production input                     value")
        for name, value in self.parameters.items():
            print(f"{name:36s} {value}")


LegacyModelPlayground = AuthenticatedLegacyModelPlayground
ModelPlayground = CanonicalModelPlayground


def main(argv=None) -> None:
    global sol, model, P
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--params", help="Parameter file for inputs and saved cache")
    args = parser.parse_args(argv)
    model_dir = str(ROOT / "code/model")
    if model_dir not in sys.path:
        sys.path.insert(0, model_dir)
    from run_model import PARAMETER_FILE
    from production.parameter_files import output_root_for, describe_saved_case
    from production.storage import load_latest
    selected = PARAMETER_FILE if args.params is None else args.params
    try:
        cached, case_dir = load_latest(output_root_for(selected))
    except FileNotFoundError:
        sol = None
        print(f"No saved solution yet for {selected}; ready for an explicit fixed-price solve.")
    else:
        sol = cached.solution
        print(describe_saved_case(case_dir, selected))
        print(f"Saved price: {cached.price:.15g}")
        print(f"V shape={sol.V.shape}; wealth grid nodes={len(cached.b_grid)}; arrays loaded with no solve")
    model = CanonicalModelPlayground(parameter_file=selected)
    P = model.P
    model.show_parameters()
    print("Ready: model.solve(), model.params['beta_annual']=0.98, result.aggregates(), result.plot_policy(), result.plot_aggregates()")


if __name__ == "__main__":
    main()
