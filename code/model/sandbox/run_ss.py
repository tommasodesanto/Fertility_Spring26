#!/usr/bin/env python3
"""Solve the 2007 general-equilibrium stationary state for one spec.

    PYTHONPATH=. python sandbox/run_ss.py --spec baseline
    PYTHONPATH=. python sandbox/run_ss.py --spec kappa_h_zero --fast

Reuses the production package's own solver, closure roots and moment
extractor (code/model/intergen_eqscale_seq_optimized/, and the closure/root
helpers in code/model/tools/audit_closed_reproductive_closure.py and
code/model/tools/run_e5f_transition_calibration.py) rather than
reimplementing any of them. See sandbox/README.md for the full contract.

Writes exactly four files to --out (default output/model/sandbox/<spec>/):
summary.md, moments.csv, parameters.csv, graphs.pdf.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

SANDBOX_ROOT = Path(__file__).resolve().parent
MODEL_ROOT = SANDBOX_ROOT.parent
REPO_ROOT = MODEL_ROOT.parents[1]
TOOLS_ROOT = MODEL_ROOT / "tools"
RETAINED_DIR = REPO_ROOT / "output/model/e5f_final_night_20260913/corrected_initial"
DEFAULT_OUT_ROOT = REPO_ROOT / "output/model/sandbox"
PROFILE = "e5f-floor"
FERTILITY_TARGET = 2.1
FERTILITY_TOLERANCE = 1.0e-6
FULL_NB = 120
FAST_NB = 40

sys.path[:0] = [str(MODEL_ROOT), str(TOOLS_ROOT)]

import mechanisms  # noqa: E402  (sandbox/mechanisms.py)
import target_table  # noqa: E402  (sandbox/target_table.py)
from spec_io import load_spec  # noqa: E402


PERIOD_YEARS = 4.0  # code/model/intergen_eqscale_seq_optimized/run_e1_chain.py: PERIOD_YEARS


def load_retained_theta() -> tuple[dict[str, float], float, dict[str, Any]]:
    candidate = json.loads((RETAINED_DIR / "candidate_result.json").read_text())
    theta = dict(candidate["proposal"]["parameters"])
    psi_child = float(candidate["score"]["normalization"]["psi_child"])
    return theta, psi_child, candidate


def display_theta_to_overrides(display_theta: dict[str, float]) -> dict[str, float]:
    """Undo the reporting-only renames in parameters.csv / candidate_result.json.

    Both files report search-domain coordinates under DISPLAY names, not the
    literal override keys `run_model_cp_dt` accepts:

      - "beta_annual" reports `theta["beta"] ** (1/PERIOD_YEARS)`
        (run_e1_chain.py:346-361, parameter_rows() in
        run_e5f_transition_calibration.py:1511-1512). Inverting: override
        "beta" = beta_annual ** PERIOD_YEARS.
      - "h_P" is a reparametrization of the first-child room jump that
        replaced the (jump, slope) pair once hbar_child_rooms (the slope) was
        fixed at zero: h_P = hbar_first_child_jump + hbar_child_rooms
        (code/model/tools/build_e5f_utility_review_packet.py:87-95,
        "old_jump = h_P - hbar_child_rooms"). Inverting with
        hbar_child_rooms == 0 (this candidate's zero restriction): override
        "hbar_first_child_jump" = h_P.

    All other keys (kappa_fert, kappa_fert_continuation, chi, H0, theta0,
    theta1, first_birth_fixed_cost) are literal override names already.
    """
    overrides = dict(display_theta)
    if "beta_annual" in overrides:
        overrides["beta"] = overrides.pop("beta_annual") ** PERIOD_YEARS
    if "h_P" in overrides:
        overrides["hbar_first_child_jump"] = overrides.pop("h_P") - float(overrides.get("hbar_child_rooms", 0.0))
    return overrides


def build_overrides(closure: Any, chain: Any, nb: int, theta: dict[str, float]) -> dict[str, float]:
    """Exactly mirror run_e5f_transition_calibration.py's REPAIRED_MODEL_PROFILE
    construction (lines ~1768-1773), which is what actually produced the
    retained candidate (its 9 free names -- including first_birth_fixed_cost
    and h_P/hbar_first_child_jump -- match E5F_INCOME_ENTRY_DOMAIN, not the
    plainer E5F_DOMAIN that profile="e5f-floor" alone gives you):

        base = closure.make_overrides(chain, theta, nb=nb, profile="e5f-floor")
        base.update(profile_overrides)   # = e5f_income_entry_overrides(), the
                                          # 15-state e6b permanent-income grid
        base.update(theta)               # theta wins over the profile's own
                                          # first_birth_fixed_cost=0.0 default

    The first build of this sandbox omitted the `e5f_income_entry_overrides()`
    merge entirely (no income-entry profile, no tenure_choice_kappa=0.005),
    which is why --spec baseline landed at a materially different equilibrium
    (psi=0.115 vs retained 0.149, and moments off by up to 2.4).
    """
    from intergen_eqscale_seq_optimized.e5f_income_entry_profile import e5f_income_entry_overrides

    raw_theta = display_theta_to_overrides(theta)
    overrides = closure.make_overrides(chain, raw_theta, nb=nb, profile=PROFILE)
    overrides.update(e5f_income_entry_overrides())
    overrides.update(raw_theta)
    overrides["hbar_child_rooms"] = 0.0  # retained calibration: "zero restriction"
    overrides["tenure_choice_kappa"] = 0.005  # retained: externally fixed (E5_FIXED default is 0.0)
    return overrides


EXPECTED_ATTRIBUTES = {
    "fecundity_omega1": 0.02, "fecundity_omega2": 0.134,
    "delta": 1.0 - (1.0 - 0.011) ** PERIOD_YEARS, "q": (1.0 + 0.02) ** PERIOD_YEARS - 1.0,
    "tenure_choice_kappa": 0.005, "alpha_cons": 0.733, "hbar_child_rooms": 0.0,
    "child_state_mode": "independent_count", "child_room_floor": True,
}


def print_attribute_diff(P: Any) -> None:
    """Diagnostic: compare the solved parameter object against the coordinator's
    checklist (survival, fecundity, independent_count, delta, q,
    tenure_choice_kappa, alpha_cons, hbar_child_rooms, and the 15-state income
    grid from e5f_income_entry_profile/e6b_profile)."""
    print("ATTRIBUTE_DIFF_CHECK")
    for name, expected in EXPECTED_ATTRIBUTES.items():
        actual = getattr(P, name, "<missing>")
        ok = "OK" if actual == expected or (isinstance(expected, float) and isinstance(actual, (int, float))
                                             and math.isclose(actual, expected, rel_tol=0, abs_tol=1e-12)) else "MISMATCH"
        print(f"  {name}: expected={expected!r} actual={actual!r} [{ok}]")
    z_grid = np.asarray(getattr(P, "z_grid", []))
    print(f"  z_grid (income states): count={z_grid.size} (expected 15 under e5f_income_entry_overrides)")
    survival = np.asarray(getattr(P, "survival_probs", []))
    print(f"  survival_probs: shape={survival.shape}")


def apply_spec(overrides: dict[str, float], spec: dict[str, Any]) -> tuple[dict[str, float], dict[str, Any]]:
    """Merge a spec's parameter overrides on top of the retained calibration.

    Returns (merged_overrides, sandbox_switch_dict) where the second element
    holds only the three mechanism-switch keys (for reporting).
    """
    merged = dict(overrides)
    switch_keys = ("child_benefit_form", "child_benefit_curvature", "scale_weighting", "child_earnings_penalty")
    switches = {}
    for key, value in spec.get("overrides", {}).items():
        if key in switch_keys:
            switches[key] = value
        merged[key] = value
    mechanisms.check_switches_supported(switches)
    return merged, switches


def warm_start_price(overrides: dict[str, float], baseline_out: Path | None) -> np.ndarray | None:
    if baseline_out is None:
        return None
    parameters_csv = baseline_out / "parameters.csv"
    if not parameters_csv.exists():
        return None
    with parameters_csv.open() as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        if row["parameter"] == "_solved_price":
            price = np.asarray(json.loads(row["estimate"]), dtype=float)
            overrides["p_init_override"] = price
            return price
    return None


def solve_stationary_state(
    chain: Any,
    calib: Any,
    overrides: dict[str, float],
    *,
    initial_psi: float,
    fix_psi: bool,
    psi_mode: str = "root",
    trace_path: Path | None = None,
) -> tuple[Any, Any, np.ndarray, dict[str, Any], int]:
    with mechanisms.sandbox_context():
        if psi_mode == "joint":
            import joint_psi

            sol, P, price, seconds, diagnostics = joint_psi.solve_old_steady_state_joint(
                chain,
                overrides,
                initial_psi=initial_psi,
                completed_fertility_target=FERTILITY_TARGET,
                completed_fertility_tolerance=FERTILITY_TOLERANCE,
                normalize=not fix_psi,
                trace_path=trace_path,
            )
        else:
            sol, P, price, seconds, diagnostics = calib.solve_old_steady_state(
                chain,
                overrides,
                initial_psi=initial_psi,
                completed_fertility_target=FERTILITY_TARGET,
                completed_fertility_tolerance=FERTILITY_TOLERANCE,
                normalize=not fix_psi,
            )
    return sol, P, price, diagnostics, int(diagnostics["stationary_solves"])


def parameter_rows(theta: dict[str, float], P: Any, price: np.ndarray, spec: dict[str, Any]) -> list[dict[str, Any]]:
    bounds = {
        "beta_annual": (0.94, 0.99),
        "kappa_fert": (0.02, 50.0),
        "kappa_fert_continuation": (0.02, 50.0),
        "chi": (0.1, 5.0),
        "H0": (0.2, 80.0),
        "theta0": (0.0, 8.0),
        "theta1": (0.02, 16.0),
        "first_birth_fixed_cost": (0.0, 8.0),
        "h_P": (0.1, 2.3),
    }
    rows = []
    for name, value in theta.items():
        lower, upper = bounds.get(name, (None, None))
        near_bound = False
        if lower is not None:
            span = max(upper - lower, 1e-12)
            near_bound = min(abs(value - lower), abs(upper - value)) <= 0.02 * span
        rows.append({
            "parameter": name,
            "estimate": value,
            "lower": lower,
            "upper": upper,
            "near_bound": near_bound,
            "status": "structural coordinate (retained, possibly spec-overridden)",
        })
    rows.append({"parameter": "psi_child", "estimate": float(P.psi_child), "lower": None, "upper": None,
                 "near_bound": False, "status": "normalized to completed fertility 2.1 unless spec fixes psi"})
    # payroll_tax and housing_supply_elasticity are display names from the dated
    # transition reporting (run_e5f_transition_calibration.py); the plain
    # stationary solve exposes the first as P.tau_pay directly, and has no
    # "supply_rule.elasticity" object at all (that belongs to the dated
    # transition's supply rule, calendar.normalize_date0_housing_supply,
    # lines ~2006-2030 -- not part of a bare run_model_cp_dt stationary call).
    rows.append({"parameter": "payroll_tax", "estimate": getattr(P, "tau_pay", None), "lower": None,
                 "upper": None, "near_bound": False, "status": "externally fixed (P.tau_pay)"})
    rows.append({"parameter": "housing_supply_elasticity", "estimate": None, "lower": None, "upper": None,
                 "near_bound": False,
                 "status": "not applicable: this is a dated-transition supply_rule.elasticity value "
                           "(0.63 in the retained report); the bare stationary solve has no such object"})
    for name in ("tenure_choice_kappa", "alpha_cons", "sigma",
                 "child_benefit_form", "child_benefit_curvature", "scale_weighting"):
        rows.append({"parameter": name, "estimate": getattr(P, name, None), "lower": None, "upper": None,
                     "near_bound": False, "status": "externally fixed or sandbox switch"})
    rows.append({"parameter": "_solved_price", "estimate": json.dumps(list(np.asarray(price).ravel())),
                 "lower": None, "upper": None, "near_bound": False, "status": "warm-start cache (not a model parameter)"})
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty csv: {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_graphs(sol: Any, P: Any, out_dir: Path) -> None:
    """Reuse the package's own 17-graph diagnostic packet, plus policy-function
    plots over wealth by age/children/tenure, as one multi-page PDF."""
    import shutil

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from intergen_eqscale_seq_optimized.diagnostics import write_diagnostics

    png_dir = out_dir / "_standard_diagnostics_pngs"
    if png_dir.exists():
        shutil.rmtree(png_dir)
    png_dir.mkdir(parents=True)
    try:
        write_diagnostics(sol, P, png_dir)  # writes 17 PNGs plus summary.json into png_dir
        png_paths = sorted(png_dir.glob("*.png"))

        with PdfPages(out_dir / "graphs.pdf") as pdf:
            for png_path in png_paths:
                image = plt.imread(png_path)
                fig, ax = plt.subplots(figsize=(11, 8.5))
                ax.imshow(image)
                ax.axis("off")
                ax.set_title(png_path.stem.replace("_", " "))
                pdf.savefig(fig)
                plt.close(fig)
            for fig in policy_function_figures(sol, P):
                pdf.savefig(fig)
                plt.close(fig)
    finally:
        shutil.rmtree(png_dir, ignore_errors=True)


def policy_function_figures(sol: Any, P: Any) -> list[Any]:
    """Policy functions over wealth, by age, children-at-home and tenure.

    These are supplemental to (not a replacement for) the standard 17-graph
    packet, per the project's diagnostic-graph-set stability rule. Solved
    policy arrays are shaped (Nb, tenure_types, I, J, [Nz,] n_parity,
    n_child_states) depending on the income process
    (solver.py:2397/2885) -- rather than hard-code one layout, this locates
    the wealth axis (matches b_grid's length) and the age axis (length P.J)
    by shape and is defensive about everything else: any failure here is
    caught and simply yields fewer supplemental pages, since the mandatory
    17-graph packet (written separately, above) is unaffected either way.
    """
    import matplotlib.pyplot as plt

    figures: list[Any] = []
    try:
        b_grid = np.asarray(getattr(sol, "b_grid", None))
        hR_pol = np.asarray(getattr(sol, "hR_pol", None))
        if b_grid.ndim != 1 or hR_pol.ndim < 2:
            return figures
        wealth_axis = next(ax for ax, size in enumerate(hR_pol.shape) if size == b_grid.size)
        age_axis = next(ax for ax, size in enumerate(hR_pol.shape) if size == int(P.J) and ax != wealth_axis)
        parity_axis = next(
            (ax for ax, size in enumerate(hR_pol.shape) if size == int(P.n_parity) and ax not in (wealth_axis, age_axis)),
            None,
        )
        ages_to_plot = sorted({0, int(P.J) // 2, int(P.J) - 1})
        fig, axes = plt.subplots(1, len(ages_to_plot), figsize=(5 * len(ages_to_plot), 4), squeeze=False)
        for column, age_index in enumerate(ages_to_plot):
            ax = axes[0][column]
            n_parity_values = hR_pol.shape[parity_axis] if parity_axis is not None else 1
            for nn in range(n_parity_values):
                indexer: list[Any] = [0] * hR_pol.ndim
                indexer[wealth_axis] = slice(None)
                indexer[age_axis] = age_index
                if parity_axis is not None:
                    indexer[parity_axis] = nn
                series = hR_pol[tuple(indexer)]
                if series.ndim != 1 or series.shape[0] != b_grid.size:
                    continue
                ax.plot(b_grid, series, label=f"children={nn}")
            ax.set_title(f"Renter housing policy, age index {age_index}")
            ax.set_xlabel("liquid wealth b")
            ax.set_ylabel("rented rooms hR")
            ax.legend(fontsize=7)
        fig.suptitle("Policy functions over wealth by age and children at home (renter rooms; first tenure/income slice)")
        figures.append(fig)
    except (StopIteration, AttributeError, TypeError, ValueError, IndexError):
        return figures
    return figures


def write_summary(
    out_dir: Path,
    spec_name: str,
    spec: dict[str, Any],
    switches: dict[str, Any],
    diagnostics: dict[str, Any],
    target_rows: list[dict[str, Any]],
    loss: float,
    parameter_table: list[dict[str, Any]],
    timings: dict[str, Any],
    warm_start_note: str,
) -> None:
    lines = [f"# Stationary-state sandbox run: {spec_name}", ""]
    lines.append(f"Spec file: `sandbox/specs/{spec_name}.yaml`")
    lines.append(f"Overrides applied: {json.dumps(spec.get('overrides', {}))}")
    lines.append(f"Mechanism switches: {json.dumps(switches) if switches else '(none; all at default-off)'}")
    lines.append("")
    lines.append("## Timing")
    lines.append(f"- Wall time per Bellman/root evaluation: {timings['per_evaluation_seconds']:.2f} s (mean over "
                 f"{timings['evaluations']} evaluation(s))")
    lines.append(f"- Total wall time: {timings['total_seconds']:.2f} s")
    lines.append(f"- Number of root evaluations (psi-normalization GE solves): {timings['evaluations']}")
    lines.append(f"- Warm start: {warm_start_note}")
    lines.append(f"- psi-normalization status: {diagnostics['status']}")
    lines.append(f"- Grid: {'fast (coarse Nb)' if timings['fast'] else 'full production Nb'} "
                 f"(Nb={timings['nb']}, J=17)")
    lines.append("")
    lines.append(f"## Loss: {loss:.6f}")
    lines.append("")
    lines.append("## 13-row target table")
    lines.append("| Moment | Target | Model | Gap | Weight | Loss contribution | Source |")
    lines.append("|---|---:|---:|---:|---:|---:|---|")
    for row in target_rows:
        lines.append(
            f"| {row['label']} | {row['target']:.6g} | {row['model']:.6g} | {row['gap']:.4g} | "
            f"{row['weight'] if row['weight'] is not None else '-'} | "
            f"{row['loss_contribution'] if row['loss_contribution'] is not None else '-'} | {row['source']} |"
        )
    lines.append("")
    lines.append("**Caveat on the 4 fertility-timing rows** (childless women 40-44, exactly-one-child "
                  "mothers 40-44, mean first-birth age, first births at 30+): the retained "
                  "target_fit.csv computes these from a specialized dated period/cohort-timing "
                  "pipeline (`transition_cross_section_moments` / `cohort_timing_moments` in "
                  "run_e5f_transition_calibration.py), not from the plain stationary "
                  "`extract_moments()` output this sandbox reuses. The sandbox reports the closest "
                  "available stationary analogue and flags it 'approximate' below; do not treat "
                  "these four rows as exact reproductions.")
    lines.append("")
    lines.append("## Parameter table")
    lines.append("| Parameter | Estimate | Lower | Upper | Near bound | Status |")
    lines.append("|---|---:|---:|---:|---|---|")
    for row in parameter_table:
        if row["parameter"] == "_solved_price":
            continue
        lines.append(f"| {row['parameter']} | {row['estimate']} | {row['lower']} | {row['upper']} | "
                     f"{row['near_bound']} | {row['status']} |")
    lines.append("")
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", required=True)
    parser.add_argument("--out", default=None)
    parser.add_argument("--fast", action="store_true", help="Coarser wealth grid (Nb=%d instead of %d)." % (FAST_NB, FULL_NB))
    parser.add_argument("--warm-start-from", default=None,
                        help="Output directory of a prior run to warm-start the price root from.")
    parser.add_argument("--package-root", default=None,
                        help="Directory containing an alternate intergen_eqscale_seq_optimized/ "
                             "(e.g. a fetched pre-main snapshot) to import instead of this repo's "
                             "code/model/. Overrides a spec's package_root key. Default: none, "
                             "i.e. this sandbox imports the package on main as it always has.")
    args = parser.parse_args()

    spec = load_spec(SANDBOX_ROOT / "specs" / f"{args.spec}.yaml")
    out_dir = Path(args.out) if args.out else DEFAULT_OUT_ROOT / args.spec
    out_dir.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(TOOLS_ROOT))

    package_root = args.package_root or spec.get("package_root")
    if package_root:
        # Inserted ahead of MODEL_ROOT/TOOLS_ROOT (both already in sys.path)
        # so BOTH `import intergen_eqscale_seq_optimized` and the tools
        # modules below (audit_closed_reproductive_closure,
        # run_e5f_transition_calibration) resolve to this snapshot instead
        # of main's package. Needed to reproduce the September 13
        # corrected_initial state, which was solved on branch
        # codex/balanced-social-security at commit 70abd4a8 plus a corrected
        # solver.py -- main differs from that snapshot in 8 files. See
        # sandbox/README.md's package-version notice.
        package_root_path = Path(package_root)
        if not package_root_path.is_absolute():
            package_root_path = (REPO_ROOT / package_root_path).resolve()
        sys.path.insert(0, str(package_root_path / "tools"))
        sys.path.insert(0, str(package_root_path))
        # `mechanisms` (imported at this file's top, before args/spec are
        # known) already triggered `from intergen_eqscale_seq_optimized
        # import solver`, caching main's copy in sys.modules -- inserting
        # package_root into sys.path *after* that does nothing on its own.
        # Purge the cached package and re-execute mechanisms.py so its
        # `_solver`/`_ORIGINAL_PRECOMPUTE_SHARED` bindings pick up the
        # snapshot's solver module instead.
        import importlib

        for name in list(sys.modules):
            if name == "intergen_eqscale_seq_optimized" or name.startswith("intergen_eqscale_seq_optimized."):
                del sys.modules[name]
        importlib.reload(mechanisms)

    import audit_closed_reproductive_closure as closure
    import run_e5f_transition_calibration as calib
    from intergen_eqscale_seq_optimized import solver as _resolved_solver

    def _sha256(path: Path) -> str:
        import hashlib
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()

    resolved_solver_path = Path(_resolved_solver.__file__).resolve()
    provenance = {
        "package_root_requested": package_root,
        "resolved_solver_path": str(resolved_solver_path),
        "resolved_solver_sha256": _sha256(resolved_solver_path),
        "main_solver_sha256": _sha256(MODEL_ROOT / "intergen_eqscale_seq_optimized/solver.py"),
        "resolved_closure_path": str(Path(closure.__file__).resolve()),
        "resolved_calib_path": str(Path(calib.__file__).resolve()),
    }
    print("PACKAGE_PROVENANCE " + json.dumps(provenance))
    if package_root and resolved_solver_path == (MODEL_ROOT / "intergen_eqscale_seq_optimized/solver.py").resolve():
        raise RuntimeError(
            "--package-root was given but intergen_eqscale_seq_optimized.solver still "
            f"resolved to main's copy ({resolved_solver_path}); sys.path priority failed."
        )

    chain = closure.load_chain(profile=PROFILE)
    theta, retained_psi, candidate = load_retained_theta()

    nb = FAST_NB if args.fast else FULL_NB
    overrides = build_overrides(closure, chain, nb, theta)
    overrides, switches = apply_spec(overrides, spec)

    psi_mode = str(spec.get("psi_mode", "root")).lower()
    if psi_mode not in ("root", "fixed", "joint"):
        raise ValueError(f"psi_mode must be 'root', 'fixed', or 'joint', got {psi_mode!r}")
    fix_psi = psi_mode == "fixed" or bool(spec.get("fix_psi", False))
    initial_psi = float(spec.get("psi_child", retained_psi)) if fix_psi else retained_psi

    warm_start_note = "none available (no baseline run found)"
    baseline_dir = Path(args.warm_start_from) if args.warm_start_from else (DEFAULT_OUT_ROOT / "baseline")
    if args.spec != "baseline":
        price = warm_start_price(overrides, baseline_dir if baseline_dir.exists() else None)
        if price is not None:
            warm_start_note = f"price root warm-started from {baseline_dir}/parameters.csv (saved price {list(price)})"
    if warm_start_note.startswith("none"):
        warm_start_note += (
            ". Bellman/value-function warm start is NOT implemented: "
            "run_model_cp_dt (solver.py:1064) has no V-init override hook "
            "analogous to p_init_override, so the value function is always "
            "cold-started inside the package; only the price root can be "
            "warm-started from here."
        )

    t0 = time.perf_counter()
    trace_path = (out_dir / "residual_trace.csv") if psi_mode == "joint" else None
    sol, P, price, diagnostics, evaluations = solve_stationary_state(
        chain, calib, overrides, initial_psi=initial_psi, fix_psi=fix_psi, psi_mode=psi_mode,
        trace_path=trace_path,
    )
    total_seconds = time.perf_counter() - t0
    print_attribute_diff(P)

    moments = chain.extract_moments(sol, P)
    target_rows, loss = target_table.build_target_rows(moments, diagnostics)

    parameter_table = parameter_rows(theta, P, price, spec)
    write_csv(out_dir / "parameters.csv", parameter_table)

    moment_rows = [{"moment": row["label"], "model": row["model"], "source": row["source"]} for row in target_rows]
    write_csv(out_dir / "moments.csv", moment_rows)

    write_graphs(sol, P, out_dir)

    timings = dict(
        per_evaluation_seconds=diagnostics["stationary_solve_seconds"] / max(evaluations, 1),
        evaluations=evaluations,
        total_seconds=total_seconds,
        fast=args.fast,
        nb=nb,
    )
    write_summary(out_dir, args.spec, spec, switches, diagnostics, target_rows, loss, parameter_table, timings, warm_start_note)

    print(f"Wrote {out_dir}/{{summary.md,moments.csv,parameters.csv,graphs.pdf}}")
    print(f"loss={loss:.6f} evaluations={evaluations} total_seconds={total_seconds:.1f}")


if __name__ == "__main__":
    main()
