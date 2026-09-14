"""Plot a frozen native successive-surprise candidate; never solve or edit slides."""
from pathlib import Path
import argparse
import csv
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    root = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "output/model/e5f_original_queue_20260913a/long_successive_refit/readout/frozen_big_run.json")
    args = parser.parse_args()
    packet = json.loads(args.input.read_text())
    out = args.input.parent
    get = lambda key: packet[key]["content"]
    rows = get("round_01/mappings/mapping_07/rows.json")
    fert = get("round_01/mappings/mapping_07/fertility.json")
    best = get("round_01/best_so_far.json")
    receipt = get("terminal/root_receipt.json")
    verified = get("terminal/verified_endpoint.json")
    initial = get("initial_reference")["initial"]
    reference = initial["quantities"]
    terminal = receipt["endpoint_reference"]
    years = np.array([r["calendar_year"] for r in rows])
    np.testing.assert_array_equal(years, 2007 + 4 * np.arange(104))
    np.testing.assert_array_equal(years, [f["calendar_year"] for f in fert])
    coords = np.r_[[r["asset_price"] for r in rows],
                   [r["pension_period_units"] for r in rows],
                   [r["equal_transfer_period_units"] for r in rows]]
    np.testing.assert_array_equal(coords, best["prices"])
    assert best["evaluation"] == 7 and best["mapping_valid"]
    assert receipt["stationary_endpoint_verified"] and verified["verified"]
    assert all(verified["audit"]["checks"].values())
    assert all(r["psi_child"] == verified["psi"] == best["psi"] for r in rows)
    np.testing.assert_allclose(fert[0]["age_mass"], initial["fertility"]["age_mass"], atol=2e-12, rtol=0)
    np.testing.assert_allclose(rows[0]["adult_population"], reference["adult_population"], atol=2e-12, rtol=0)
    tfr = np.array([f["period_tfr_topcode_adjusted"] for f in fert])
    for row, f in zip(rows, fert):
        np.testing.assert_allclose(sum(f["birth_flow_topcode_adjusted"]), row["birth_children_topcode_adjusted"], atol=2e-12, rtol=0)
        np.testing.assert_allclose(sum(f["age_specific_birth_rate_topcode_adjusted"]), f["period_tfr_topcode_adjusted"], atol=2e-12, rtol=0)
    initial_tfr = initial["fertility"]["period_tfr_topcode_adjusted"]
    # There is no mortality over fertile ages in this closure: each stationary
    # fertile-age cell has entrant mass E, so sum_j(B_j/E) = total births / E.
    terminal_tfr = terminal["birth_children_topcode_adjusted"] / terminal["entry_flow"]
    np.testing.assert_allclose(terminal_tfr, 2.1 * terminal["renewal_ratio"], atol=2e-7, rtol=0)
    hh = np.array([r["adult_population"] for r in rows]) / reference["adult_population"] * 100
    housing = np.array([r["housing_demand"] for r in rows]) / reference["housing_demand"] * 100
    price = np.array([r["asset_price"] for r in rows]) / reference["asset_price"] * 100
    terminal_values = [terminal_tfr, terminal["population_households"] / reference["adult_population"] * 100,
                       terminal["housing_demand"] / reference["housing_demand"] * 100,
                       terminal["asset_price"] / reference["asset_price"] * 100]
    series = [(tfr, initial_tfr, "Period fertility", "Births per woman"),
              (hh, 100, "Household heads", "Initial steady state = 100"),
              (housing, 100, "Total housing", "Initial steady state = 100"),
              (price, 100, "House price", "Initial steady state = 100")]
    plt.rcParams.update({"font.size": 12, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 7.8))
    x = np.r_[2003, 2007, years]
    for ax, (values, pre, title, ylabel), endpoint in zip(axes.flat, series, terminal_values):
        y = np.r_[pre, pre, values]
        artist, = ax.plot(x, y, color="#d97815", lw=2.3, label="Current forecast")
        np.testing.assert_array_equal(artist.get_xdata(), x)
        np.testing.assert_array_equal(artist.get_ydata(), y)
        ax.axhline(pre, color="#777777", lw=.8, ls=":")
        ax.axvline(2007, color=".65", lw=.8, ls=":")
        ax.axhline(endpoint, color="#c43c39", lw=1.1, ls="--", label="Verified new steady state")
        ax.plot(2423, endpoint, "D", color="#c43c39", ms=5)
        ax.set(title=title, ylabel=ylabel, xlabel="Start of four-year period", xlim=(2003, 2428),
               xticks=[2007, 2103, 2203, 2303, 2423])
        ax.grid(axis="y", alpha=.16)
    axes[0, 0].legend(frameon=False, fontsize=10, loc="lower right")
    fig.subplots_adjust(left=.075, right=.98, top=.94, bottom=.11, wspace=.25, hspace=.43)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"big_run_transition.{ext}", dpi=160, facecolor="white")
    plt.close(fig)
    table = [dict(year=int(y), period_fertility=float(tfr[i]), households_index=float(hh[i]),
                  housing_index=float(housing[i]), house_price_index=float(price[i])) for i, y in enumerate(years)]
    with (out / "big_run_transition.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(table[0])); writer.writeheader(); writer.writerows(table)
    result = dict(model_solves=0, status="First candidate of first surprise; not a fitted four-shock history",
                  stage=0, psi=best["psi"], mapping=7, dates=len(rows), artist_arrays_verified=True,
                  native_best_coordinates_match_exactly=True, terminal_steady_state_verified=True,
                  first_target=get("../contract.json")["target"], first_model=float(tfr[0]),
                  max_housing_gap_percent=100*max(abs(r["relative_market_residual"]) for r in rows),
                  max_paygo_gap_percent=100*max(abs(r["scaled_pension_budget_residual"]) for r in rows),
                  max_rebate_gap_percent=100*max(abs(r["scaled_government_budget_residual"]) for r in rows),
                  terminal_values=terminal_values, terminal_distance=get("round_01/mappings/mapping_07/terminal_distance.json"),
                  terminal_tfr_definition="Stationary adjusted births divided by entrants; equal fertile-age masses and no fertile-age mortality",
                  selected_rows=[table[i] for i in (0, 4, 14, 103)],
                  frozen_input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest())
    (out / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
