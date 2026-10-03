"""Read-only aggregates and raw policy plots for saved model solutions.

The aggregate contract uses ``g`` as the realized, post-tenure current
cross-section. ``g_beginning_distribution`` is the post-fertility,
pre-tenure distribution and is used only for inherited assets and its pooled
wealth distribution. Consumption is a flow per model period; with the default
four-year period it is four-year consumption. Financial assets are stocks in
mean annual gross-earnings units. The change in financial assets includes the
housing transaction mapping and is neither national-account saving nor
next-period financial assets.
"""
from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np


_POLICY_FIELDS = (
    "c_pol", "bp_pol", "c_pol_stay", "bp_pol_stay", "hR_pol",
)
_MASS_TOL = 1.0e-10
_VALUE_INFEASIBLE = -1.0e9


def _get(sol: Any, name: str, default: Any = None) -> Any:
    return sol.get(name, default) if isinstance(sol, Mapping) else getattr(sol, name, default)


def _required(sol: Any, name: str) -> np.ndarray:
    value = _get(sol, name)
    if value is None:
        raise ValueError(f"saved solution is missing required array {name!r}")
    return np.asarray(value)


def _houses(sol: Any, houses: Any, tenure_count: int) -> np.ndarray:
    metadata = _get(sol, "houses", _get(sol, "house_sizes", None))
    if houses is None:
        houses = metadata if metadata is not None else (2, 4, 6, 8, 10)
    result = np.asarray(houses, dtype=float).reshape(-1)
    if result.size + 1 != tenure_count:
        raise ValueError(
            f"houses describes {result.size} owner tenures but solution has {tenure_count} tenures"
        )
    if metadata is not None and not np.array_equal(result, np.asarray(metadata, dtype=float).reshape(-1)):
        raise ValueError("houses argument conflicts with saved solution metadata")
    if not np.all(np.isfinite(result)) or np.any(result <= 0):
        raise ValueError("house sizes must be finite and positive")
    return result


def _check_mass(name: str, mass: np.ndarray) -> None:
    if not np.all(np.isfinite(mass)):
        raise ValueError(f"{name} contains nonfinite distribution mass")
    if np.any(mass < -_MASS_TOL):
        raise ValueError(f"{name} contains negative distribution mass")


def _profiles(
    sol: Any,
    houses: Any,
    age_start: int,
    period_years: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    b = _required(sol, "b_grid").astype(float, copy=False).reshape(-1)
    g = _required(sol, "g").astype(float, copy=False)
    stay = _required(sol, "g_stay_distribution").astype(float, copy=False)
    beginning = _required(sol, "g_beginning_distribution").astype(float, copy=False)
    if g.ndim != 7:
        raise ValueError(f"g must be seven-dimensional, got shape {g.shape}")
    if g.shape[0] != b.size or any(x.shape != g.shape for x in (stay, beginning)):
        raise ValueError("b_grid, g, g_stay_distribution, and g_beginning_distribution shapes conflict")
    if not np.all(np.isfinite(b)) or np.any(np.diff(b) <= 0):
        raise ValueError("b_grid must be finite and strictly increasing")
    for name, arr in (("g", g), ("g_stay_distribution", stay),
                      ("g_beginning_distribution", beginning)):
        _check_mass(name, arr)
    total = float(g.sum())
    if total <= 0:
        raise ValueError("realized distribution g has no population mass")
    if not np.isclose(float(beginning.sum()), total, rtol=0.0, atol=_MASS_TOL):
        raise ValueError("g_beginning_distribution and realized g have different population mass")
    if not np.allclose(beginning.sum(axis=(0, 1, 2, 4, 5, 6)),
                       g.sum(axis=(0, 1, 2, 4, 5, 6)), rtol=0.0, atol=_MASS_TOL):
        raise ValueError("g_beginning_distribution and realized g have different per-age mass")
    if np.any(stay > g + _MASS_TOL):
        raise ValueError("g_stay_distribution is not a subset of realized g")
    if np.any(stay[:, 0, ...] > _MASS_TOL):
        raise ValueError("renter stayer mass must be zero")

    size = g.shape
    hh = _houses(sol, houses, size[1])
    policies = {name: _required(sol, name).astype(float, copy=False) for name in _POLICY_FIELDS}
    if any(x.shape != size for x in policies.values()):
        raise ValueError("a policy array shape differs from the saved distribution shape")
    values = _get(sol, "V")
    V = None if values is None else np.asarray(values)
    if V is not None and V.shape != size:
        raise ValueError("V shape differs from the saved distribution shape")

    buyer_mass = g - stay
    if np.any(buyer_mass < -_MASS_TOL):
        raise ValueError("stayer mass exceeds realized mass beyond tolerance")
    # Round only tolerated subtraction noise to zero.
    buyer_mass = np.maximum(buyer_mass, 0.0)
    if V is not None:
        occupied_buyer = buyer_mass > 0
        occupied_stay = stay > 0
        if np.any(~np.isfinite(V[occupied_buyer])) or np.any(V[occupied_buyer] <= _VALUE_INFEASIBLE):
            raise ValueError("buyer-policy mass occupies a nonfinite or initialized-infeasible value entry")
        if np.any(~np.isfinite(V[occupied_stay])) or np.any(V[occupied_stay] <= _VALUE_INFEASIBLE):
            raise ValueError("stayer-policy mass occupies a nonfinite or initialized-infeasible value entry")

    c_flow = buyer_mass * 0.0
    bp_flow = buyer_mass * 0.0
    # Never multiply a zero-mass state by an initialized NaN/Inf policy value.
    buyer = buyer_mass > 0
    stayer = stay > 0
    c_flow[buyer] = policies["c_pol"][buyer] * buyer_mass[buyer]
    bp_flow[buyer] = policies["bp_pol"][buyer] * buyer_mass[buyer]
    c_flow[stayer] += policies["c_pol_stay"][stayer] * stay[stayer]
    bp_flow[stayer] += policies["bp_pol_stay"][stayer] * stay[stayer]
    if (not np.all(np.isfinite(policies["c_pol"][buyer]))
            or not np.all(np.isfinite(policies["bp_pol"][buyer]))
            or not np.all(np.isfinite(policies["c_pol_stay"][stayer]))
            or not np.all(np.isfinite(policies["bp_pol_stay"][stayer]))):
        raise ValueError("occupied consumption or next-assets policy entries are nonfinite")

    # g is indexed by the realized tenure. Owners occupy fixed house sizes;
    # renter services come from the conditional saved renter housing policy.
    renter = g[:, 0, ...] > 0
    renter_h = policies["hR_pol"][:, 0, ...]
    if not np.all(np.isfinite(renter_h[renter])):
        raise ValueError("occupied renter housing policy entries are nonfinite")

    bshape = (b.size,) + (1,) * (g.ndim - 1)
    b_grid = b.reshape(bshape)
    inherited_assets_mass = beginning * b_grid
    # b is inherited wealth in the post-fertility, pre-tenure distribution.
    # Subtracting it from b' retains any current housing-transaction wealth
    # shift as part of the financial-asset change.
    change_mass = bp_flow - inherited_assets_mass
    room_mass = np.zeros_like(g)
    room_mass[:, 0, ...][renter] = renter_h[renter] * g[:, 0, ...][renter]
    for ten, room_count in enumerate(hh, start=1):
        room_mass[:, ten, ...] = g[:, ten, ...] * room_count
    own_mass = g.copy()
    own_mass[:, 0, ...] = 0.0

    by_age: list[dict[str, Any]] = []
    for age_i in range(size[3]):
        age_mass = float(np.sum(g[:, :, :, age_i, :, :, :]))
        denominator_assets = float(np.sum(beginning[:, :, :, age_i, :, :, :]))
        row: dict[str, Any] = {
            "age": int(age_start + period_years * age_i),
            "population_mass": age_mass,
            "mean_consumption": None if age_mass <= 0 else float(c_flow[:, :, :, age_i, :, :, :].sum() / age_mass),
            "mean_inherited_assets": None if denominator_assets <= 0 else float(inherited_assets_mass[:, :, :, age_i, :, :, :].sum() / denominator_assets),
            "mean_next_assets": None if age_mass <= 0 else float(bp_flow[:, :, :, age_i, :, :, :].sum() / age_mass),
            "mean_asset_change": None if age_mass <= 0 else float(change_mass[:, :, :, age_i, :, :, :].sum() / age_mass),
            "mean_rooms": None if age_mass <= 0 else float(room_mass[:, :, :, age_i, :, :, :].sum() / age_mass),
            "ownership_rate": None if age_mass <= 0 else float(own_mass[:, :, :, age_i, :, :, :].sum() / age_mass),
        }
        by_age.append(row)

    inherited_pooled = beginning.sum(axis=tuple(range(1, beginning.ndim)))
    pooled_mass = float(inherited_pooled.sum())
    # All returned values are finite JSON primitives/lists; empty age bins use null.
    overall = {
        "mean_consumption": float(c_flow.sum() / total),
        "mean_inherited_assets": float(inherited_assets_mass.sum() / pooled_mass),
        "mean_next_assets": float(bp_flow.sum() / total),
        "mean_asset_change": float(change_mass.sum() / total),
        "mean_rooms": float(room_mass.sum() / total),
        "ownership_rate": float(own_mass.sum() / total),
    }
    result = {
        "units": {
            "consumption": f"flow per {period_years}-year model period, in mean annual gross-earnings units",
            "financial_assets": "stock in mean annual gross-earnings units",
            "rooms": "rooms",
            "ownership_rate": "fraction of realized household mass",
            "population_mass": "saved stationary household mass",
            "asset_change_definition": "mean next-period b' minus inherited b from the post-fertility, pre-tenure distribution; includes housing transactions and is not national-account saving",
        },
        "population_mass": total,
        "overall": overall,
        "by_age": by_age,
        "inherited_asset_distribution": {
            "asset_nodes": b.tolist(),
            "mass": inherited_pooled.tolist(),
            "probability": (inherited_pooled / pooled_mass).tolist(),
        },
    }
    return result


def aggregate_solution(
    sol: Any,
    houses: Any = None,
    age_start: int = 18,
    period_years: int = 4,
) -> dict[str, Any]:
    """Return realized population aggregates and inherited wealth distribution.

    Required native arrays are stored either as ``sol`` attributes or mapping
    keys. ``g`` is realized after location/tenure transactions, while
    ``g_beginning_distribution`` is post-fertility and pre-tenure. Buyer
    policies apply to ``g - g_stay_distribution``; owner-stayer policies apply
    to the stayer subset. Empty age cells are encoded as JSON ``null``.
    """
    if period_years <= 0:
        raise ValueError("period_years must be positive")
    return _profiles(sol, houses, int(age_start), int(period_years))


def _selector_index(sol: Any, age: int, income: int, branch: int,
                    children: int, at_home: int, age_start: int,
                    period_years: int) -> tuple[np.ndarray, int]:
    b = _required(sol, "b_grid").astype(float, copy=False).reshape(-1)
    shape = _required(sol, "g").shape
    age_index_float = (float(age) - age_start) / period_years
    if not age_index_float.is_integer():
        raise ValueError("age must match a model decision date")
    age_i = int(age_index_float)
    checks = ((branch, shape[1], "branch/tenure"), (income, shape[4], "income"),
              (children, shape[5], "children"), (at_home, shape[6], "at_home"),
              (age_i, shape[3], "age"))
    for value, length, label in checks:
        if not 0 <= value < length:
            raise ValueError(f"{label} selector {value} is outside [0, {length})")
    return b, age_i


def plot_policy(
    sol: Any,
    variable: str = "consumption",
    age: int = 30,
    income: int = 4,
    branch: int = 0,
    children: int = 0,
    at_home: int = 0,
    owner_policy: str = "buying",
    axis: str = "assets",
    xlim: tuple[float, float] | None = None,
    houses: Any = None,
    age_start: int = 18,
    period_years: int = 4,
):
    """Plot one raw conditional policy slice and return its Matplotlib Axes.

    ``branch`` is the saved tenure index (0 renter, 1+ owner house sizes).
    Policies are displayed at original grid nodes without interpolation,
    tenure averaging, or support masking.
    """
    import matplotlib.pyplot as plt

    b, age_i = _selector_index(
        sol, age, income, branch, children, at_home,
        age_start, period_years,
    )
    if owner_policy not in {"buying", "staying"}:
        raise ValueError("owner_policy must be 'buying' or 'staying'")
    if axis not in {"assets", "node_index", "earnings"}:
        raise ValueError("axis must be 'assets', 'node_index', or 'earnings'")
    ix = (slice(None), branch, 0, age_i, income, children, at_home)
    if variable == "consumption":
        field = "c_pol_stay" if branch > 0 and owner_policy == "staying" else "c_pol"
        y = _required(sol, field)[ix]
        ylabel = f"Consumption per {period_years}-year model period / mean annual gross earnings"
    elif variable == "next_assets":
        field = "bp_pol_stay" if branch > 0 and owner_policy == "staying" else "bp_pol"
        y, ylabel = _required(sol, field)[ix], "Next-period financial assets"
    elif variable == "housing":
        hh = _houses(sol, houses, _required(sol, "g").shape[1])
        y = (_required(sol, "hR_pol")[ix] if branch == 0
             else np.full(b.size, hh[branch - 1]))
        ylabel = "Housing services (rooms)"
    elif variable == "assets_distribution":
        y = _required(sol, "g_beginning_distribution")[ix]
        ylabel = "Mass at inherited-asset node"
    else:
        raise ValueError("variable must be consumption, next_assets, housing, or assets_distribution")
    if axis == "node_index":
        x = np.arange(1, b.size + 1)
        xlabel = "Saved asset-grid node index"
    else:
        x = b.copy()
        xlabel = ("Financial assets / mean annual gross earnings" if axis == "earnings"
                  else "Assets b (model units)")
    fig, ax = plt.subplots()
    ax.plot(x, y, marker=".")
    x_note = ("b is in mean annual gross-earnings units" if axis == "assets" else "")
    ax.set(xlabel=xlabel, ylabel=ylabel,
           title=f"{variable.replace('_', ' ').title()} · age {age}, income state {income + 1}, tenure {branch}"
                 + (f" · {x_note}" if x_note else ""))
    if xlim is not None:
        if xlim[0] > xlim[1]:
            raise ValueError("xlim must be an increasing pair")
        ax.set_xlim(*xlim)
        visible = np.isfinite(x) & np.isfinite(y) & (x >= xlim[0]) & (x <= xlim[1])
        if not np.any(visible):
            raise ValueError("xlim contains no finite saved policy nodes")
        y_min, y_max = float(np.min(y[visible])), float(np.max(y[visible]))
        pad = 0.05 * (y_max - y_min) if y_max > y_min else max(abs(y_min) * 0.05, 0.05)
        ax.set_ylim(y_min - pad, y_max + pad)
    return ax


def plot_aggregates(sol: Any, houses: Any = None, age_start: int = 18,
                    period_years: int = 4, wealth_range: str = "central"):
    """Draw per-age aggregates and the pooled inherited-asset distribution."""
    import matplotlib.pyplot as plt

    data = aggregate_solution(sol, houses, age_start, period_years)
    if wealth_range not in {"central", "all"}:
        raise ValueError("wealth_range must be 'central' or 'all'")
    rows = data["by_age"]
    ages = np.asarray([row["age"] for row in rows])
    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    for ax, field, label in (
        (axes[0, 0], "mean_consumption", "Consumption per model period"),
        (axes[0, 1], "mean_next_assets", "Next-period financial assets"),
        (axes[0, 2], "mean_rooms", "Mean housing services (rooms)"),
        (axes[1, 0], "ownership_rate", "Ownership rate"),
    ):
        vals = [row[field] for row in rows]
        ax.plot(ages, [np.nan if v is None else v for v in vals], marker="o")
        ax.set(xlabel="Age", ylabel=label)
    distribution = data["inherited_asset_distribution"]
    wealth_x = np.asarray(distribution["asset_nodes"], dtype=float)
    wealth_mass = np.asarray(distribution["mass"], dtype=float)
    if wealth_range == "central":
        cdf = np.cumsum(wealth_mass / wealth_mass.sum())
        lower = max(0, int(np.searchsorted(cdf, 0.0001)) - 1)
        upper = min(wealth_x.size - 1, int(np.searchsorted(cdf, 0.995)) + 1)
        visible = (wealth_x >= wealth_x[lower]) & (wealth_x <= wealth_x[upper])
        axes[1, 1].set_xlim(wealth_x[lower], wealth_x[upper])
    else:
        visible = np.ones(wealth_x.size, dtype=bool)
    axes[1, 1].plot(wealth_x[visible], wealth_mass[visible], marker=".")
    mass_min, mass_max = float(np.min(wealth_mass[visible])), float(np.max(wealth_mass[visible]))
    mass_pad = 0.05 * (mass_max - mass_min) if mass_max > mass_min else max(abs(mass_min) * 0.05, 0.05)
    axes[1, 1].set_ylim(mass_min - mass_pad, mass_max + mass_pad)
    axes[1, 1].set(xlabel="Assets b (model units); mean annual gross-earnings units",
                   ylabel="Population mass at node")
    axes[1, 2].axis("off")
    axes[1, 2].text(0, 0.8, f"Population mass: {data['population_mass']:.8g}", va="top")
    axes[1, 2].text(0, 0.6, f"Mean inherited assets: {data['overall']['mean_inherited_assets']:.6g}", va="top")
    axes[1, 2].text(0, 0.4, f"Mean asset change: {data['overall']['mean_asset_change']:.6g}", va="top")
    fig.tight_layout()
    return fig, axes
