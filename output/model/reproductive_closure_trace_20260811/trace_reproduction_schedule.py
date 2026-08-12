#!/usr/bin/env python3
"""Fixed-price trace of the reproduction and demand schedules (smoke grade).

At each fixed asset price P, solve the household problem and normalized
stationary distribution once (no equilibrium loop), holding the funded
baseline policy fixed (1% property tax, lump-sum transfer at its baseline
value). Record the demographic and housing aggregates that the reproductive
closure needs:

    entry_rate            E  (replacement requirement per household)
    entrants_mature_total B  (mature entrant households per household)
    total_births_kfe      births per household per period
    housing_demand        D  (physical rooms per household)

Smoke settings (Nb=40, loose tolerance, J=17 preserved). Diagnostic only.
"""

from __future__ import annotations

import csv
import sys
import time
import traceback
from pathlib import Path

import numpy as np

ROOT = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
MODEL_DIR = ROOT / "code/model"
TOOLS_DIR = MODEL_DIR / "tools"
for p in (str(TOOLS_DIR), str(MODEL_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

SOURCE = (
    ROOT
    / "output/model/intergen_e5f_child_room_floor_psinneg_extended_20260806/report/results.json"
)
BASELINE_PRICE = 0.776082707252744
BASELINE_TRANSFER = 0.2091064453125
ANNUAL_TAX = 0.01
OUT_CSV = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("trace_points.csv")
PRICES = (
    [float(x) for x in sys.argv[2:]]
    if len(sys.argv) > 2
    else [
        0.25,
        0.33,
        0.40,
        0.47,
        0.54,
        0.62,
        0.70,
        BASELINE_PRICE,
        0.86,
        0.95,
        1.08,
        1.25,
    ]
)

import run_e5_repaired_policy_with_entry as drv  # noqa: E402

theta, base_overrides, target_system, arm = drv.load_repaired_profile(SOURCE, smoke=True)
import run_intergen_funded_policy_with_entry as closure  # noqa: E402
from run_intergen_funded_property_tax_test import case_overrides  # noqa: E402

FIELDS = [
    "price",
    "rent",
    "entry_rate_E",
    "mature_entrants_B",
    "births_per_household",
    "housing_demand_D",
    "housing_supply_HS",
    "B_over_E",
    "lifetime_births_nu",
    "clearing_scale_S",
    "status",
    "seconds",
]

write_header = not OUT_CSV.exists()
print(f"arm={arm}; {len(PRICES)} price points -> {OUT_CSV}", flush=True)
for price in PRICES:
    started = time.time()
    row = {name: "" for name in FIELDS}
    row["price"] = price
    try:
        overrides = case_overrides(
            base_overrides, ANNUAL_TAX, False, BASELINE_TRANSFER, True
        )
        overrides.update(solve_mode="pe", p_fixed=np.array([float(price)], dtype=float))
        solution, parameters, _ = closure.run_model_cp_dt(overrides, verbose=False)
        user_cost = float(parameters.user_cost_rate)
        entry = float(solution.entry_rate)
        mature = float(solution.entrants_mature_total)
        births = float(solution.total_births_kfe)
        demand = float(np.asarray(solution.housing_demand, dtype=float).reshape(-1)[0])
        supply = float(parameters.H0[0]) * (
            user_cost * price / float(parameters.r_bar[0])
        ) ** float(parameters.xi_supply[0])
        row.update(
            rent=user_cost * price,
            entry_rate_E=entry,
            mature_entrants_B=mature,
            births_per_household=births,
            housing_demand_D=demand,
            housing_supply_HS=supply,
            B_over_E=mature / entry,
            lifetime_births_nu=births / entry,
            clearing_scale_S=supply / demand,
            status="ok",
        )
    except Exception as error:  # noqa: BLE001 - diagnostic sweep records failures
        row["status"] = f"failed: {type(error).__name__}: {error}"
        traceback.print_exc()
    row["seconds"] = round(time.time() - started, 1)
    with OUT_CSV.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        if write_header:
            writer.writeheader()
            write_header = False
        writer.writerow(row)
    print(
        f"P={price:.4f}: {row['status']} ({row['seconds']}s) "
        f"B/E={row['B_over_E'] if row['B_over_E'] != '' else 'n/a'}",
        flush=True,
    )
print("done", flush=True)
