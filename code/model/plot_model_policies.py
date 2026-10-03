"""Create editable policy plots from the latest saved model run.

Run this script after ``run_model.py``. It only loads a saved run; it never
initializes the reference model or runs a solve. Edit the settings below, then
run the file again to write figures under that run's ``policy_plots`` folder.
"""
from __future__ import annotations

import sys
import textwrap
from pathlib import Path


# ---------------------------------------------------------------------------
# Author settings. Ages are model ages in years; income states and owner-size
# indices are human-counted from 1. Child counts n and m use model definitions.
# ---------------------------------------------------------------------------
RUN_DIRECTORY = None  # None selects the latest saved author run.
AGE_YEARS = [26, 30, 34]
INCOME_STATES = [3, 5, 7]  # Human 1-based productivity-state indices.
CHILDREN_EVER_BORN = 0  # n: children ever born.
CHILDREN_AT_HOME = 0  # m: children currently at home.
FAMILY_STATES = [(0, 0), (1, 1), (2, 1), (2, 2)]  # (n, m); impossible pairs are skipped.

POLICY_HOUSING_BRANCH = "renter"  # "renter" or "owner".
POLICY_OWNER_SIZE_INDEX = 1  # Human 1-based index in P.H_own when branch is owner.
OWNER_POLICY = "buying"  # "buying" or "staying" for a selected owner branch.
INHERITED_TENURE = "renter"  # Used for ownership and birth-attempt probability objects.
INHERITED_OWNER_SIZE_INDEX = 1  # Human 1-based index in P.H_own when inherited tenure is owner.

AXIS = "assets"  # "assets", "earnings" (same saved b values), or "node_index".
WEALTH_LIMITS = "central"  # Or an explicit pair such as (-10, 40); None uses the full grid.
NODE_INDEX_LIMITS = None  # Optional node-index limits when AXIS == "node_index".
SHOW_FIGURES = False

ROOT = Path(__file__).resolve().parents[2]
TOOLS_DIR = Path(__file__).resolve().parent / "tools"
def _tenure_index(kind: str, owner_size_index: int, houses) -> int:
    if kind == "renter":
        return 0
    if kind != "owner":
        raise ValueError("Tenure setting must be 'renter' or 'owner'.")
    if not 1 <= int(owner_size_index) <= len(houses):
        raise ValueError(f"Owner-size index must be between 1 and {len(houses)}.")
    return int(owner_size_index)


def _age_index(age: int, age_start: int, period_years: int, n_ages: int) -> int:
    step = (float(age) - age_start) / period_years
    if not step.is_integer():
        raise ValueError(f"Age {age} is not a model decision age.")
    result = int(step)
    if not 0 <= result < n_ages:
        raise ValueError(f"Age {age} is outside the saved model age range.")
    return result


def _family_state(n: int, m: int, shape) -> None:
    if n < 0 or n >= shape[5] or m < 0 or m >= shape[6]:
        raise ValueError(f"Family state (n={n}, m={m}) is outside the saved state space.")
    if m > n:
        raise ValueError(f"Impossible family state: m={m} children at home exceeds n={n} ever born.")


def _axis_values(b):
    import numpy as np

    if AXIS == "assets":
        return b, "Assets b (model units)"
    if AXIS == "earnings":
        return b, "Financial assets / mean annual gross earnings"
    if AXIS == "node_index":
        return np.arange(1, b.size + 1, dtype=float), "Saved asset-grid node index (1-based)"
    raise ValueError("AXIS must be 'assets', 'earnings', or 'node_index'.")


def _limits(x, b, solution):
    limits = NODE_INDEX_LIMITS if AXIS == "node_index" else WEALTH_LIMITS
    if limits is None:
        return float(x[0]), float(x[-1])
    if limits == "central" and AXIS != "node_index":
        import numpy as np

        beginning = np.asarray(solution.g_beginning_distribution, dtype=float)
        pooled = beginning.sum(axis=tuple(range(1, beginning.ndim)))
        total = float(pooled.sum())
        if total <= 0.0 or not np.isfinite(total):
            raise ValueError("Cannot derive the central wealth range: saved population mass is empty.")
        cdf = np.cumsum(pooled) / total
        lower = max(0, int(np.searchsorted(cdf, 0.0001)) - 1)
        upper = min(b.size - 1, int(np.searchsorted(cdf, 0.995)) + 1)
        return float(b[lower]), float(b[upper])
    if isinstance(limits, str):
        raise ValueError("WEALTH_LIMITS must be 'central', None, or an increasing pair.")
    if len(limits) != 2 or not limits[0] < limits[1]:
        raise ValueError("Axis limits must be an increasing pair.")
    return float(limits[0]), float(limits[1])


def _state_value(values, valid):
    import numpy as np

    y = np.asarray(values, dtype=float)
    mask = np.asarray(valid, dtype=bool) & np.isfinite(y)
    return y, mask


def _draw_policy(ax, x, values, valid, label):
    """Draw finite saved values at original grid nodes without feasibility claims."""
    import numpy as np

    y, mask = _state_value(values, valid)
    ax.plot(x, np.where(mask, y, np.nan), marker=".", markersize=5, label=label)
    return y, mask


def _finish_axis(ax, x, plotted_values, xlabel, ylabel, title, limits):
    """Set exact grid x endpoints and y limits from visible finite values."""
    import numpy as np

    lo, hi = limits
    ax.set_xlim(lo, hi)
    visible = (x >= lo) & (x <= hi)
    y_parts = []
    for values, valid in plotted_values:
        y = np.asarray(values, dtype=float)
        use = np.asarray(valid, dtype=bool) & visible & np.isfinite(y)
        if np.any(use):
            y_parts.append(y[use])
    if not y_parts:
        raise ValueError(f"No finite plotted values are visible for: {title}")
    joined = np.concatenate(y_parts)
    ymin, ymax = float(np.min(joined)), float(np.max(joined))
    pad = 0.05 * (ymax - ymin) if ymax > ymin else max(abs(ymin) * 0.05, 0.05)
    ax.set_ylim(ymin - pad, ymax + pad)
    if " · " in title:
        title = title.replace(" · ", " ·\n", 1)
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)
    if title.startswith("Raw "):
        note = (
            "Raw arrays can contain initialized entries at infeasible states; plotted nodes do not certify "
            "feasibility or occupancy. Owner-buying branches use saved b-grid coordinates and can reflect "
            "transaction-shifted wealth."
        )
        ax.figure.text(0.5, 0.012, textwrap.fill(note, width=92),
                       ha="center", va="bottom", fontsize=7)
    elif "probability" in title.lower():
        note = "State-specific probability only; plotted asset nodes do not certify feasibility or population occupancy."
        ax.figure.text(0.5, 0.012, textwrap.fill(note, width=92),
                       ha="center", va="bottom", fontsize=7)


def _plot_run(result, run_dir: Path) -> list[Path]:
    import numpy as np
    import matplotlib.pyplot as plt

    sol, P = result.solution, result.P
    b = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    houses = np.asarray(P.H_own, dtype=float).reshape(-1)
    V = np.asarray(sol.V)
    shape = V.shape
    if len(shape) != 7 or shape[0] != b.size:
        raise ValueError(f"Expected saved seven-dimensional policy arrays; V has shape {shape}.")
    period_years = int(P.period_years)
    age_start = int(P.age_start)
    if OWNER_POLICY not in {"buying", "staying"}:
        raise ValueError("OWNER_POLICY must be 'buying' or 'staying'.")
    if not AGE_YEARS or not INCOME_STATES:
        raise ValueError("AGE_YEARS and INCOME_STATES must each include at least one value.")
    for age in AGE_YEARS:
        _age_index(age, age_start, period_years, shape[3])
    for state in INCOME_STATES:
        if not 1 <= int(state) <= shape[4]:
            raise ValueError(f"Income state must be between 1 and {shape[4]}.")
    _family_state(CHILDREN_EVER_BORN, CHILDREN_AT_HOME, shape)

    policy_tenure = _tenure_index(POLICY_HOUSING_BRANCH, POLICY_OWNER_SIZE_INDEX, houses)
    inherited_tenure = _tenure_index(INHERITED_TENURE, INHERITED_OWNER_SIZE_INDEX, houses)
    if policy_tenure == 0 and OWNER_POLICY == "staying":
        raise ValueError("OWNER_POLICY='staying' requires POLICY_HOUSING_BRANCH='owner'.")
    x, xlabel = _axis_values(b)
    axis_limits = _limits(x, b, sol)
    output = run_dir / "policy_plots"
    output.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    figures_to_show = []

    def save(fig, stem):
        path = output / f"{stem}.png"
        fig.tight_layout(rect=(0, 0.12, 1, 1))
        fig.savefig(path, dpi=160, bbox_inches="tight")
        written.append(path)
        if SHOW_FIGURES:
            plt.show(block=False)
            figures_to_show.append(fig)
        else:
            plt.close(fig)

    def age_income_family(age, income_state, n=CHILDREN_EVER_BORN, m=CHILDREN_AT_HOME):
        _family_state(n, m, shape)
        j = _age_index(age, age_start, period_years, shape[3])
        z = int(income_state) - 1
        ix = (slice(None), policy_tenure, 0, j, z, n, m)
        return j, z, ix

    # Figure 1: raw consumption policies at the selected housing branch, by age.
    fig, ax = plt.subplots(figsize=(7, 4.5))
    plotted = []
    for age in AGE_YEARS:
        _, _, ix = age_income_family(age, INCOME_STATES[0])
        field = sol.c_pol_stay if policy_tenure > 0 and OWNER_POLICY == "staying" else sol.c_pol
        y, mask = _draw_policy(ax, x, field[ix], np.ones(b.size, dtype=bool), f"Age {age}")
        plotted.append((y, mask))
    branch_name = ("renter" if policy_tenure == 0 else
                   f"owner, {houses[policy_tenure - 1]:g} rooms, {OWNER_POLICY} branch")
    _finish_axis(ax, x, plotted, xlabel,
        f"Consumption per {period_years}-year model period\n(mean annual gross-earnings units)",
        f"Raw conditional consumption · {branch_name} · income state {INCOME_STATES[0]}, n={CHILDREN_EVER_BORN}, m={CHILDREN_AT_HOME}", axis_limits)
    save(fig, "01_consumption_by_age")

    # Figure 2: next-period financial assets from the same raw conditional branch.
    fig, ax = plt.subplots(figsize=(7, 4.5))
    plotted = []
    for age in AGE_YEARS:
        _, _, ix = age_income_family(age, INCOME_STATES[0])
        field = sol.bp_pol_stay if policy_tenure > 0 and OWNER_POLICY == "staying" else sol.bp_pol
        y, mask = _draw_policy(ax, x, field[ix], np.ones(b.size, dtype=bool), f"Age {age}")
        plotted.append((y, mask))
    _finish_axis(ax, x, plotted, xlabel, "Next-period financial assets\n(mean annual gross-earnings units)",
        f"Raw conditional next assets · {branch_name} · income state {INCOME_STATES[0]}, n={CHILDREN_EVER_BORN}, m={CHILDREN_AT_HOME}", axis_limits)
    save(fig, "02_next_assets_by_age")

    # Figure 3: renter housing services, with no owner-size policy substituted.
    fig, ax = plt.subplots(figsize=(7, 4.5))
    plotted = []
    for age in AGE_YEARS:
        j = _age_index(age, age_start, period_years, shape[3])
        z = INCOME_STATES[0] - 1
        ix = (slice(None), 0, 0, j, z, CHILDREN_EVER_BORN, CHILDREN_AT_HOME)
        y, mask = _draw_policy(ax, x, sol.hR_pol[ix], np.ones(b.size, dtype=bool), f"Age {age}")
        plotted.append((y, mask))
    _finish_axis(ax, x, plotted, xlabel, "Renter housing services (rooms)",
        f"Raw renter housing policy · income state {INCOME_STATES[0]}, n={CHILDREN_EVER_BORN}, m={CHILDREN_AT_HOME}", axis_limits)
    save(fig, "03_renter_housing_by_age")

    # Figure 4: probability of choosing ownership now, conditional on inherited tenure.
    fig, ax = plt.subplots(figsize=(7, 4.5))
    plotted = []
    for age in AGE_YEARS:
        j = _age_index(age, age_start, period_years, shape[3])
        z = INCOME_STATES[0] - 1
        ix = (slice(None), inherited_tenure, 0, j, z, CHILDREN_EVER_BORN, CHILDREN_AT_HOME)
        probabilities = np.asarray(sol.tenure_probs[ix], dtype=float)
        values = probabilities[:, 1:].sum(axis=1)
        normalized = (np.isfinite(probabilities).all(axis=1)
                      & (probabilities.sum(axis=1) > 0.99)
                      & (probabilities.sum(axis=1) < 1.01))
        y, mask = _draw_policy(ax, x, values, normalized, f"Age {age}")
        plotted.append((y, mask))
    inherited_name = ("renter" if inherited_tenure == 0 else
                      f"owner, {houses[inherited_tenure - 1]:g} rooms")
    _finish_axis(ax, x, plotted, xlabel, "Probability",
        f"Ownership probability · conditional on inherited {inherited_name}, income state {INCOME_STATES[0]}, n={CHILDREN_EVER_BORN}, m={CHILDREN_AT_HOME}", axis_limits)
    save(fig, "04_ownership_probability_by_age")

    # Figure 5: first-birth attempt probability, explicitly among childless households.
    fig, ax = plt.subplots(figsize=(7, 4.5))
    plotted = []
    for age in AGE_YEARS:
        j = _age_index(age, age_start, period_years, shape[3])
        z = INCOME_STATES[0] - 1
        values = np.asarray(sol.fert_probs[
            (slice(None), inherited_tenure, 0, j, z, 1)], dtype=float)
        y, mask = _draw_policy(ax, x, values, np.isfinite(values), f"Age {age}")
        plotted.append((y, mask))
    _finish_axis(ax, x, plotted, xlabel, "Probability",
        f"First-birth attempt probability · childless, m=0, inherited {inherited_name}, income state {INCOME_STATES[0]}", axis_limits)
    save(fig, "05_birth_attempt_probability_by_age")

    # Figure 6: consumption policies across the requested human-counted income states.
    fig, ax = plt.subplots(figsize=(7, 4.5))
    plotted = []
    age = AGE_YEARS[0]
    for income_state in INCOME_STATES:
        _, _, ix = age_income_family(age, income_state)
        field = sol.c_pol_stay if policy_tenure > 0 and OWNER_POLICY == "staying" else sol.c_pol
        y, mask = _draw_policy(ax, x, field[ix], np.ones(b.size, dtype=bool), f"Income state {income_state}")
        plotted.append((y, mask))
    _finish_axis(ax, x, plotted, xlabel,
        f"Consumption per {period_years}-year model period\n(mean annual gross-earnings units)",
        f"Raw conditional consumption · age {age}, {branch_name}, n={CHILDREN_EVER_BORN}, m={CHILDREN_AT_HOME}", axis_limits)
    save(fig, "06_consumption_by_income_state")

    # Figure 7: next-assets policies across the requested income states.
    fig, ax = plt.subplots(figsize=(7, 4.5))
    plotted = []
    for income_state in INCOME_STATES:
        _, _, ix = age_income_family(age, income_state)
        field = sol.bp_pol_stay if policy_tenure > 0 and OWNER_POLICY == "staying" else sol.bp_pol
        y, mask = _draw_policy(ax, x, field[ix], np.ones(b.size, dtype=bool), f"Income state {income_state}")
        plotted.append((y, mask))
    _finish_axis(ax, x, plotted, xlabel, "Next-period financial assets\n(mean annual gross-earnings units)",
        f"Raw conditional next assets · age {age}, {branch_name}, n={CHILDREN_EVER_BORN}, m={CHILDREN_AT_HOME}", axis_limits)
    save(fig, "07_next_assets_by_income_state")

    # Figure 8: consumption across selected family states; impossible pairs are skipped.
    fig, ax = plt.subplots(figsize=(7, 4.5))
    plotted = []
    for n, m in FAMILY_STATES:
        if m > n:
            print(f"Skipping impossible child state (n={n}, m={m}).")
            continue
        try:
            _family_state(n, m, shape)
        except ValueError as exc:
            print(f"Skipping family state (n={n}, m={m}): {exc}")
            continue
        _, _, ix = age_income_family(age, INCOME_STATES[0], n, m)
        field = sol.c_pol_stay if policy_tenure > 0 and OWNER_POLICY == "staying" else sol.c_pol
        y, mask = _draw_policy(ax, x, field[ix], np.ones(b.size, dtype=bool), f"n={n}, m={m}")
        plotted.append((y, mask))
    if not plotted:
        plt.close(fig)
        raise ValueError("FAMILY_STATES contains no valid family states to plot.")
    _finish_axis(ax, x, plotted, xlabel,
        f"Consumption per {period_years}-year model period\n(mean annual gross-earnings units)",
        f"Raw conditional consumption · age {age}, income state {INCOME_STATES[0]}, {branch_name}", axis_limits)
    save(fig, "08_consumption_by_family_state")

    if SHOW_FIGURES:
        plt.show()
        for fig in figures_to_show:
            plt.close(fig)
    return written


def main() -> list[Path]:
    """Load the selected saved run and write eight editable policy figures."""
    if str(TOOLS_DIR) not in sys.path:
        sys.path.insert(0, str(TOOLS_DIR))
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    try:
        if RUN_DIRECTORY is None:
            from production.storage import load_latest
            result, run_dir = load_latest(ROOT / "output/model/local_solution")
        else:
            from production.storage import load_case
            result, run_dir = load_case(RUN_DIRECTORY)
    except ImportError as exc:
        raise RuntimeError(
            "Saved-run loading is unavailable. Save a completed run with "
            "code/model/run_model.py, then run this plotter again."
        ) from exc
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "No completed local GE case was found. Run code/model/run_model.py first, "
            "then run this script; it never solves implicitly."
        ) from exc
    paths = _plot_run(result, Path(run_dir))
    print(f"Saved {len(paths)} policy figures to {Path(run_dir) / 'policy_plots'}")
    return paths


if __name__ == "__main__":
    main()
