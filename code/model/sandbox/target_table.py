"""Build the 13-row target-fit table from a solved stationary state.

Reuses the retained calibration's own target/weight values (from
output/model/e5f_final_night_20260913/corrected_initial/target_fit.csv) and
the package's own `chain.extract_moments(sol, P)` dict for model values. It
does not reimplement any moment definition.

Important limitation, stated here and repeated in every summary.md: the
retained target_fit.csv's 12 scored rows were produced by a specialized dated
period/cohort-timing measurement pipeline
(`transition_cross_section_moments` / `cohort_timing_moments` in
code/model/tools/run_e5f_transition_calibration.py), which differs from the
plain stationary `extract_moments()` this sandbox reuses for two of the four
fertility-timing rows and for the room-count rows (which use uncapped
analogues, since extract_moments has no rooms-capped-at-9 measure). Four rows
(1 normalization + 3 with a reasonable extract_moments analogue) are close
approximations; one row (exactly-one-child-among-mothers-40-44) has no
analogue in extract_moments at all. See MOMENT_KEY_MAP below for the row-by
-row mapping and MOMENT_KEY_MAP[*].source for what each row actually is.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[3]
RETAINED_TARGET_FIT = REPO_ROOT / "output/model/e5f_final_night_20260913/corrected_initial/target_fit.csv"

EXACT = "exact_stationary_key"
APPROX = "approximate_stationary_analogue"
UNAVAILABLE = "unavailable_in_extract_moments"

# label (as it appears in the retained target_fit.csv) -> (extract_moments key or None, source)
MOMENT_KEY_MAP: dict[str, tuple[str | None, str]] = {
    "Initial model completed fertility": ("tfr", EXACT),
    "Childless women, ages 40–44": ("childless_rate", APPROX),
    "Exactly one child among mothers, ages 40–44": (None, UNAVAILABLE),
    "Period mean first-birth age": ("mean_age_first_birth", APPROX),
    "First births at age 30+": ("share_first_births_age30plus", APPROX),
    "Wealth / annual gross labor earnings": ("aggregate_wealth_to_annual_gross_labor_earnings", EXACT),
    "Annual bequests / aggregate wealth": ("annual_bequest_flow_to_aggregate_wealth", EXACT),
    "Old wealth/income p90 / median, ages 76–84": ("old_total_wealth_to_annual_income_p90_p50_7684", EXACT),
    "Mean occupied rooms, capped at 9": ("aggregate_mean_occupied_rooms_18_85", APPROX),
    "Ownership, heads 30–55": ("own_rate_3055", EXACT),
    "First-birth room response, −1 to +3": ("housing_increment_0to1", APPROX),
    "Rooms: 3+ versus 1–2 resident children (model dependent proxy)": (
        "prime30_55_parent_3plus_minus_1to2_mean_rooms", APPROX,
    ),
    "Recent-parent ownership gap": ("own_gap_newparent_nonparent_3055", EXACT),
}


def load_retained_targets() -> list[dict[str, str]]:
    with RETAINED_TARGET_FIT.open() as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 13:
        raise ValueError(f"Expected 13 retained target rows, found {len(rows)}")
    return rows


def build_target_rows(moments: dict[str, Any], psi_diagnostics: dict[str, Any]) -> tuple[list[dict[str, Any]], float]:
    retained = load_retained_targets()
    rows: list[dict[str, Any]] = []
    loss = 0.0
    for retained_row in retained:
        label = retained_row["label"]
        key, source = MOMENT_KEY_MAP.get(label, (None, UNAVAILABLE))
        target = float(retained_row["target"])
        scored = retained_row["scored"] == "True"
        weight = float(retained_row["actual_weight"]) if scored and retained_row["actual_weight"] else None
        if label == "Initial model completed fertility":
            model = float(psi_diagnostics["completed_fertility"])
        elif key is not None and key in moments and math.isfinite(float(moments[key])):
            model = float(moments[key])
        else:
            model = float("nan")
        gap = model - target if math.isfinite(model) else float("nan")
        contribution = weight * gap * gap if (weight is not None and math.isfinite(gap)) else None
        if contribution is not None:
            loss += contribution
        rows.append(dict(
            label=label, target=target, model=model, gap=gap, weight=weight,
            loss_contribution=contribution, scored=scored, source=source,
            extract_moments_key=key,
        ))
    return rows, loss
