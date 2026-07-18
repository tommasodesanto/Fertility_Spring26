"""Small diagnostic calibration driver for the one-market model scaffold."""

from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np

from .solver import income_at_state, run_model_cp_dt


PERIOD_YEARS = 4.0
AGE_START = 18.0
FERTILITY_START_AGE = 18.0
FERTILITY_END_AGE = 42.0
RETIREMENT_AGE = 66.0
CHILD_MATURITY_AGE = 18.0


CORE_TARGETS = {
    "own_rate": 0.575,
    "young_owner_rate": 0.25,
    "old_owner_rate": 0.764,
    "mean_completed_fertility": 0.95,
    "childless_rate": 0.20,
}


PSID_ENTRY_WEALTH_RATIO_NODES_2535 = np.array(
    [-2.51940697, -0.07907025, 0.10228762, 0.35287169, 3.03955200],
    dtype=float,
)
PSID_ENTRY_WEALTH_RATIO_WEIGHTS_2535 = np.array(
    [0.2000027, 0.2000966, 0.2000024, 0.1998877, 0.2000107],
    dtype=float,
)
PSID_ENTRY_WEALTH_RATIO_SOURCE_2535 = (
    "PSID PSIDSHELF_MOBILITY, young childless renters ages 25-35, "
    "weighted quintile-bin means of NETWORTH2R/INCFAMR; mean 0.17922556, "
    "weighted median 0.09996729."
)


def external_entry_wealth_overrides() -> dict[str, Any]:
    """Externally calibrated entrant wealth distribution in annual-income units."""
    return {
        "entry_wealth_mode": "income_ratio_distribution",
        "entry_wealth_ratio_nodes": PSID_ENTRY_WEALTH_RATIO_NODES_2535.copy(),
        "entry_wealth_ratio_weights": PSID_ENTRY_WEALTH_RATIO_WEIGHTS_2535.copy(),
        "entry_wealth_ratio_source": PSID_ENTRY_WEALTH_RATIO_SOURCE_2535,
        "entry_wealth_spread_nodes": 1,
    }


PSID_ENTRY_WEALTH_RATIO_NODES_1824 = np.array(
    [
        -2.22252859121604,
        -0.0526403178090525,
        0.104244689233362,
        0.351572998737395,
        3.1034542945567,
    ],
    dtype=float,
)
PSID_ENTRY_WEALTH_RATIO_WEIGHTS_1824 = np.array(
    [
        0.199425762696792,
        0.200377472937802,
        0.199863308508567,
        0.200087860612093,
        0.200245595244747,
    ],
    dtype=float,
)
PSID_ENTRY_WEALTH_RATIO_SOURCE_1824 = (
    "PSID reference-person childless renters ages 18-24, primary rows of "
    "code/data/psid_followup_mar2026/output/intergen_income_entry_targets_20260716/"
    "block2_entry_wealth_18_24.csv: weighted quintile-bin means of wealth over "
    "annual family income; mean 0.258855836903436 (person-bootstrap SE "
    "0.125965630995927, 499 replications), weighted median 0.0985177806342346."
)


def external_entry_wealth_overrides_1824() -> dict[str, Any]:
    """Model-entry-age (18-24) externally calibrated entrant wealth distribution."""
    return {
        "entry_wealth_mode": "income_ratio_distribution",
        "entry_wealth_ratio_nodes": PSID_ENTRY_WEALTH_RATIO_NODES_1824.copy(),
        "entry_wealth_ratio_weights": PSID_ENTRY_WEALTH_RATIO_WEIGHTS_1824.copy(),
        "entry_wealth_ratio_source": PSID_ENTRY_WEALTH_RATIO_SOURCE_1824,
        "entry_wealth_spread_nodes": 1,
    }


OLD_NONLOCATION_TARGETS = {
    "tfr": 1.70,
    "childless_rate": 0.15,
    "mean_age_first_birth": 26.0,
    "own_rate": 0.57547241,
    "own_family_gap": 0.16766167,
    "housing_increment_0to1": 0.66443467,
    "housing_increment_1to2": 0.56581378,
    "young_liquid_wealth_to_income": 0.60,
    "old_age_own_rate": 0.76426097,
    "old_age_parent_childless_gap": 0.070,
}


OLD_NONLOCATION_NO_TIMING_TARGETS = {
    k: v for k, v in OLD_NONLOCATION_TARGETS.items() if k != "mean_age_first_birth"
}


CANDIDATE_NO_TIMING_V0_TARGETS = {
    **OLD_NONLOCATION_NO_TIMING_TARGETS,
    "liquid_wealth_to_income": 1.20,
    "housing_user_cost_share": 0.24,
    "prime_childless_renter_median_rooms": 4.0,
    "prime_childless_owner_median_rooms": 6.0,
}


CANDIDATE_REPLACEMENT_V1_TARGETS = {
    "tfr": 1.70,
    "childless_rate": 0.15,
    "own_rate": 0.57547241,
    "own_family_gap": 0.16766167,
    "housing_increment_0to1": 0.66443467,
    "housing_increment_1to2": 0.48803143,
    "young_liquid_wealth_to_income": 0.17922556,
    "old_nonhousing_wealth_to_income_6575": 6.41854289,
    "old_parent_childless_nonhousing_wealth_to_income_gap_6575": 1.00744952,
    "prime30_55_childless_renter_mean_rooms": 3.80528810,
    "prime30_55_childless_owner_mean_rooms": 6.22404983,
    "prime30_55_childless_renter_share_rooms_ge6": 0.13765221,
    "prime30_55_childless_owner_share_rooms_ge6": 0.59613112,
}


CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_TARGETS = {
    **{k: v for k, v in CANDIDATE_REPLACEMENT_V1_TARGETS.items() if k != "old_nonhousing_wealth_to_income_6575"},
    "old_nonhousing_wealth_to_income_median_6575": 2.23046078,
}


CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_V1_TARGETS = {
    **CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_TARGETS,
    "old_minus_young_owner_rate_6575_2534": 0.42309487,
}


CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_OWNGAP_V1_TARGETS = {
    **{
        k: v
        for k, v in CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_TARGETS.items()
        if k != "old_parent_childless_nonhousing_wealth_to_income_gap_6575"
    },
    "old_age_parent_childless_gap": 0.08254406,
    "old_minus_young_owner_rate_6575_2534": 0.42309487,
}


CANDIDATE_REPLACEMENT_OLD_RETENTION_V1_TARGETS = {
    **CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_TARGETS,
    "old_age_own_rate": 0.76426097,
    "old_age_parent_childless_gap": 0.08254406,
}


CANDIDATE_REPLACEMENT_YOUNG_OLD_OWN_V1_TARGETS = {
    **CANDIDATE_REPLACEMENT_OLD_RETENTION_V1_TARGETS,
    "own_rate_2534": 0.34116609,
}


CANDIDATE_REPLACEMENT_YOUNG_OLD_ROOMGAP_V1_TARGETS = {
    **CANDIDATE_REPLACEMENT_YOUNG_OLD_OWN_V1_TARGETS,
    "prime30_55_childless_owner_minus_renter_mean_rooms": 2.41876173,
}


CANDIDATE_REPLACEMENT_ROOMGAP_14MOMENT_V1_TARGETS = {
    **{
        k: v
        for k, v in CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_TARGETS.items()
        if k
        not in {
            "prime30_55_childless_owner_mean_rooms",
            "prime30_55_childless_renter_share_rooms_ge6",
            "young_liquid_wealth_to_income",
        }
    },
    "young_childless_renter_liquid_wealth_to_annual_gross_income_2535": 0.17922556,
    "prime30_55_childless_owner_minus_renter_mean_rooms": 2.41876173,
    "old_age_own_rate": 0.76426097,
    "own_rate_2534": 0.34116609,
}


CANDIDATE_REPLACEMENT_ROOMGAP_14MOMENT_TFR192_V1_TARGETS = {
    **CANDIDATE_REPLACEMENT_ROOMGAP_14MOMENT_V1_TARGETS,
    "tfr": 1.92,
}


CANDIDATE_REPLACEMENT_POST_AUDIT_V1_TARGETS = {
    **{
        k: v
        for k, v in CANDIDATE_REPLACEMENT_ROOMGAP_14MOMENT_V1_TARGETS.items()
        if k != "housing_increment_1to2"
    },
    "tfr": 1.918,
    "childless_rate": 0.188,
    "prime30_55_parent_3plus_minus_1to2_mean_rooms": 0.36769955881,
}


CANDIDATE_REPLACEMENT_TOTAL_MEDIAN_V1_TARGETS = {
    **{
        k: v
        for k, v in CANDIDATE_REPLACEMENT_V1_TARGETS.items()
        if k
        not in {
            "old_nonhousing_wealth_to_income_6575",
            "old_parent_childless_nonhousing_wealth_to_income_gap_6575",
        }
    },
    "old_total_wealth_to_income_median_6575": 5.26375360,
    "old_parent_childless_total_wealth_to_income_gap_6575": 1.28413467,
}


CANDIDATE_REPLACEMENT_TOTAL_DUE_LIFECYCLE_V1_TARGETS = {
    **CANDIDATE_REPLACEMENT_TOTAL_MEDIAN_V1_TARGETS,
    "old_minus_young_owner_rate_6575_2534": 0.42309487,
}


# Minimal internally calibrated bequest block. These three PSID targets replace
# the legacy nonhousing level, parent-childless nonhousing gap, and old-age
# ownership moment while leaving the other 12 active moments unchanged.
CANDIDATE_REPLACEMENT_BEQUEST_INTERNAL_V1_TARGETS = {
    **{
        k: v
        for k, v in CANDIDATE_REPLACEMENT_POST_AUDIT_V1_TARGETS.items()
        if k
        not in {
            "old_nonhousing_wealth_to_income_median_6575",
            "old_parent_childless_nonhousing_wealth_to_income_gap_6575",
            "old_age_own_rate",
        }
    },
    "old_total_estate_wealth_to_annual_income_median_7684": 6.50131577436537,
    "old_total_estate_wealth_to_annual_income_p90_p50_7684": 3.44811075444552,
    "old_2plus_minus_1_total_estate_wealth_to_annual_income_median_gap_6575": 0.101108088567873,
}


NONHOUSING_MEDIAN_TARGET = 1.90821154211154
NONHOUSING_MEDIAN_WEIGHT = 83.74916751466371


# M4 standard-bequest recalibration set. Relative to the internal bequest set,
# this keeps the late-life estate median level, drops the p90/p50 tail ratio
# and the 2+-minus-1 family-size estate gap, and restores the legacy
# nonhousing composition median so late-life wealth composition stays
# disciplined while theta0 is the only free bequest parameter.
CANDIDATE_REPLACEMENT_BEQUEST_MEDIAN_COMPOSITION_V1_TARGETS = {
    **{
        k: v
        for k, v in CANDIDATE_REPLACEMENT_BEQUEST_INTERNAL_V1_TARGETS.items()
        if k
        not in {
            "old_total_estate_wealth_to_annual_income_p90_p50_7684",
            "old_2plus_minus_1_total_estate_wealth_to_annual_income_median_gap_6575",
        }
    },
    "old_nonhousing_wealth_to_income_median_6575": NONHOUSING_MEDIAN_TARGET,
}


# PSID share of 65-75 reference persons holding at least one year of annual
# family income in nonhousing net worth. Headline row of block3 in
# code/data/psid_followup_mar2026/output/intergen_income_entry_targets_20260716/
# (block3_composition_share_6575.csv); the weight is the inverse person-cluster
# bootstrap variance 1/SE^2 with SE 0.0102949617302331 (499 replications).
OLD_NONHOUSING_GE1X_SHARE_TARGET = 0.608333139649131
OLD_NONHOUSING_GE1X_SHARE_SE = 0.0102949617302331
OLD_NONHOUSING_GE1X_SHARE_WEIGHT = 9435.18732291246


# M5 income-disciplined recalibration set. Relative to the M4
# median-composition set, this drops the old-age nonhousing median level,
# disciplines the late-life balance sheet with the PSID nonhousing >= 1x
# annual-income share, and restores the old-age ownership rate at the legacy
# ACS-sourced value and weight (ACS-consistent with the other ownership
# targets; the PSID reference-person alternative 0.834 is documented in
# block4_oldage_ownership_6575.csv of the same 20260716 output folder).
CANDIDATE_REPLACEMENT_INCOME_DISCIPLINED_V1_TARGETS = {
    **{
        k: v
        for k, v in CANDIDATE_REPLACEMENT_BEQUEST_MEDIAN_COMPOSITION_V1_TARGETS.items()
        if k != "old_nonhousing_wealth_to_income_median_6575"
    },
    "old_nonhousing_ge_1x_income_share_6575": OLD_NONHOUSING_GE1X_SHARE_TARGET,
    "old_age_own_rate": 0.76426097,
}


# Non-behavioral target-object ledger. Keep this near the target values so
# calibration changes document the model statistic, empirical object, and known
# measurement caveats in the same edit.
TARGET_MOMENT_OBJECTS: dict[str, dict[str, str]] = {
    "tfr": {
        "model": "2 * mean_completed_fertility among post-fertility ages.",
        "data": "Completed-fertility-equivalent target, not a period TFR flow.",
        "status": "target-choice",
        "issue": "The label is legacy shorthand; reports should call this completed-fertility-equivalent.",
    },
    "childless_rate": {
        "model": "parity_dist[0], measured after the fertility window.",
        "data": "Completed childlessness / ever-childless object.",
        "status": "needs-source-pin",
        "issue": "Do not substitute ACS no-child-in-household for completed childlessness.",
    },
    "own_rate": {
        "model": "own_rate_3055: owner mass divided by total mass, ages 30-55.",
        "data": "ACS household heads, DUE-housing sample, ages 30-55.",
        "status": "clean",
        "issue": "Earlier aggregate-owner mapping was fixed; keep own_rate mapped to own_rate_3055.",
    },
    "own_family_gap": {
        "model": "own_gap_newparent_nonparent_3055: new-parent ownership minus nonparent ownership, ages 30-55.",
        "data": "ACS household heads, DUE-housing sample, new-parent minus no-child household ownership, ages 30-55.",
        "status": "mostly-clean",
        "issue": "Data no-child household and model nonparent/never-parent states are close but not identical.",
    },
    "housing_increment_0to1": {
        "model": "housing_increment_0to1_eventstudy_t3: birth cohort minus no-birth control, default horizon 0 in 4-year periods.",
        "data": "PSID first-birth rooms event-study, about 3 years post-birth.",
        "status": "clean-after-fix",
        "issue": "June 25 fix removed the old raw 12-year lifecycle-drift window.",
    },
    "housing_increment_1to2": {
        "model": "housing_increment_1to2_proxy_t3: two-plus-child birth-state housing minus one-child birth-state housing.",
        "data": "PSID quick second-birth room change, post-3.",
        "status": "proxy",
        "issue": "This is an additional-child housing-demand target, not a sequential second-birth hazard.",
    },
    "prime30_55_parent_3plus_minus_1to2_mean_rooms": {
        "model": "mean realized housing services for high-child current-parent states minus one-child current-parent states, ages 30-55.",
        "data": "ACS/MMS household heads, ages 30-55, parents with nchild > 2 minus parents with 1 <= nchild <= 2, requiring a child-under-18 signal.",
        "status": "candidate-replacement",
        "issue": "Replaces the old PSID second-birth pre/post target in the post-audit target set; maps the coarse model high-child state to the empirical 3+ child bin.",
    },
    "young_liquid_wealth_to_income": {
        "model": "young childless renter liquid wealth divided by 4-year period after-tax income.",
        "data": "PSID young childless renter NETWORTH2R / INCFAMR, annual family-income denominator.",
        "status": "measurement-error",
        "issue": "Replace in serious target sets with an explicit annual-income statistic such as young_childless_renter_liquid_wealth_to_annual_gross_income_2535.",
    },
    "young_childless_renter_liquid_wealth_to_annual_gross_income_2535": {
        "model": "young childless renter liquid wealth divided by annual gross-normalized income, ages 25-35.",
        "data": "PSID young childless renter NETWORTH2R / INCFAMR, ages 25-35.",
        "status": "candidate-replacement",
        "issue": "Mean target can use 0.17922556; median should be re-extracted as weighted before hard targeting.",
    },
    "young_childless_renter_liquid_wealth_to_annual_gross_income_2535_median": {
        "model": "weighted median of young childless renter liquid wealth over annual gross-normalized income, ages 25-35.",
        "data": "PSID median of young childless renter NETWORTH2R / INCFAMR, ages 25-35.",
        "status": "candidate-replacement",
        "issue": "Weighted PSID median is 0.099967; existing CSV median 0.061550 is unweighted.",
    },
    "old_parent_childless_nonhousing_wealth_to_income_gap_6575": {
        "model": "parent minus childless ratio-of-sums: nonhousing wealth divided by model period income, ages 65-75.",
        "data": "PSID parent minus childless mean nonhousing net worth / annual income, ages 65-75.",
        "status": "needs-fix",
        "issue": "Period-vs-annual denominator and ratio-of-sums vs mean-of-ratios are not identical.",
    },
    "old_nonhousing_wealth_to_income_median_6575": {
        "model": "weighted median of nonhousing wealth divided by annual gross income, ages 65-75.",
        "data": "PSID reference persons, weighted median nonhousing net worth / annual family income, ages 65-75.",
        "status": "internally-calibrated-balance-sheet-target",
        "issue": "Target 1.90821154211154 and inverse-variance weight 83.74916751466371 use the same person-cluster bootstrap as the estate targets (499 replications).",
    },
    "old_total_estate_wealth_to_annual_income_median_7684": {
        "model": "weighted median of (b + pH) / annual gross income, ages 76-84.",
        "data": "PSID reference persons, weighted median NETWORTHR / INCFAMR, ages 76-84.",
        "status": "internally-calibrated-bequest-target",
        "issue": "Total estate wealth includes gross housing value because liquid b is net of secured debt.",
    },
    "old_total_estate_wealth_to_annual_income_p90_p50_7684": {
        "model": "weighted p90 divided by weighted median of (b + pH) / annual gross income, ages 76-84.",
        "data": "PSID reference persons, p90/median of NETWORTHR / INCFAMR, ages 76-84.",
        "status": "internally-calibrated-bequest-target",
        "issue": "Upper-tail target identifies the luxury shift jointly with the old-age median.",
    },
    "old_2plus_minus_1_total_estate_wealth_to_annual_income_median_gap_6575": {
        "model": "2+-child minus 1-child weighted median of (b + pH) / annual gross income, ages 65-75.",
        "data": "PSID reference persons, completed children 2+ minus 1, median NETWORTHR / INCFAMR, ages 65-75.",
        "status": "internally-calibrated-bequest-target-weak",
        "issue": "The 2+ bin matches the live 0/1/2+ parity state; person-bootstrap SE is 0.563.",
    },
    "old_nonhousing_ge_1x_income_share_6575": {
        "model": "asset-distribution mass share of ages 65-75 with liquid b / annual gross income >= 1, all tenures, positive-income states.",
        "data": "PSID reference persons ages 65-75, share with nonhousing net worth of at least one year of annual family income.",
        "status": "internally-calibrated-balance-sheet-target",
        "issue": "Target 0.608333139649131 and inverse-variance weight 9435.18732291246 use the block3 person-cluster bootstrap (SE 0.0102949617302331, 499 replications).",
    },
    "prime30_55_childless_renter_mean_rooms": {
        "model": "mass-weighted mean realized renter housing services for childless-in-household ages 30-55.",
        "data": "ACS household heads, childless renters, mean ROOMS, ages 30-55.",
        "status": "clean",
        "issue": "Means are not affected by the Markov nonlinear-threshold bug.",
    },
    "prime30_55_childless_owner_mean_rooms": {
        "model": "mass-weighted mean owner housing rung for childless-in-household ages 30-55.",
        "data": "ACS household heads, childless owners, mean ROOMS, ages 30-55.",
        "status": "diagnostic-or-demoted",
        "issue": "Correlated with the owner-renter mean gap; not in the active 14-moment roomgap target.",
    },
    "prime30_55_childless_renter_share_rooms_ge6": {
        "model": "renter mass share with h_R >= 6 for childless-in-household ages 30-55.",
        "data": "ACS household heads, childless renters, share ROOMS >= 6, ages 30-55.",
        "status": "clean-after-fix",
        "issue": "Must be computed before income-state policy collapse; fixed by markov_renter_room_moments.",
    },
    "prime30_55_childless_owner_share_rooms_ge6": {
        "model": "owner mass share with owner rung H >= 6 for childless-in-household ages 30-55.",
        "data": "ACS household heads, childless owners, share ROOMS >= 6, ages 30-55.",
        "status": "clean",
        "issue": "If missed, treat as economics/product-menu problem, not measurement.",
    },
    "prime30_55_childless_owner_minus_renter_mean_rooms": {
        "model": "childless owner mean housing services minus childless renter mean housing services, ages 30-55.",
        "data": "ACS household heads, childless owner-renter mean ROOMS gap, ages 30-55.",
        "status": "clean",
        "issue": "Preferred smooth room-segmentation target in the active 14-moment roomgap system.",
    },
    "old_age_own_rate": {
        "model": "old_age_own_rate_6575: owner mass divided by total mass, ages 65-75.",
        "data": "ACS household heads, DUE-housing sample, ownership ages 65-75.",
        "status": "clean",
        "issue": "Current mismatch is structural/economic, not an obvious measurement bug.",
    },
    "old_age_parent_childless_gap": {
        "model": "old_age_parent_childless_gap_6575: parent old ownership minus childless old ownership.",
        "data": "PSID/ACS old parent-minus-childless ownership gap, ages 65-75 depending on target variant.",
        "status": "demoted-or-variant",
        "issue": "Not in active 14-moment roomgap target; source/sample must be pinned before reuse.",
    },
    "own_rate_2534": {
        "model": "owner mass divided by total mass, ages 25-34.",
        "data": "ACS household heads, DUE-housing sample, ownership ages 25-34.",
        "status": "clean",
        "issue": "Clean measurement object; current miss is economic/value-ranking or grid/product mapping.",
    },
}


CORE_FEASIBILITY_V1_TARGETS = {
    k: CANDIDATE_NO_TIMING_V0_TARGETS[k]
    for k in [
        "tfr",
        "childless_rate",
        "own_rate",
        "own_family_gap",
        "housing_increment_0to1",
        "housing_increment_1to2",
        "prime_childless_renter_median_rooms",
        "prime_childless_owner_median_rooms",
    ]
}


COST_TEST_V1_TARGETS = {
    **CORE_FEASIBILITY_V1_TARGETS,
    "housing_user_cost_share": CANDIDATE_NO_TIMING_V0_TARGETS["housing_user_cost_share"],
}


OLDAGE_TEST_V1_TARGETS = {
    **CORE_FEASIBILITY_V1_TARGETS,
    "old_age_own_rate": CANDIDATE_NO_TIMING_V0_TARGETS["old_age_own_rate"],
    "old_age_parent_childless_gap": CANDIDATE_NO_TIMING_V0_TARGETS["old_age_parent_childless_gap"],
}


ROOMCOST_TEST_V1_TARGETS = {
    k: CANDIDATE_NO_TIMING_V0_TARGETS[k]
    for k in [
        "tfr",
        "childless_rate",
        "own_rate",
        "own_family_gap",
        "old_age_own_rate",
        "old_age_parent_childless_gap",
        "housing_user_cost_share",
        "housing_increment_0to1",
        "housing_increment_1to2",
    ]
}


CORE_WEIGHTS = {
    "own_rate": 8.0,
    "young_owner_rate": 8.0,
    "old_owner_rate": 6.0,
    "mean_completed_fertility": 10.0,
    "childless_rate": 3.0,
}


OLD_NONLOCATION_WEIGHTS = {
    "tfr": 12.0,
    "childless_rate": 12.0,
    "mean_age_first_birth": 12.0,
    "own_rate": 12.0,
    "own_family_gap": 10.0,
    "housing_increment_0to1": 8.0,
    "housing_increment_1to2": 2.0,
    "young_liquid_wealth_to_income": 12.0,
    "old_age_own_rate": 10.0,
    "old_age_parent_childless_gap": 10.0,
}


OLD_NONLOCATION_NO_TIMING_WEIGHTS = {
    k: v for k, v in OLD_NONLOCATION_WEIGHTS.items() if k != "mean_age_first_birth"
}


CANDIDATE_NO_TIMING_V0_WEIGHTS = {
    **OLD_NONLOCATION_NO_TIMING_WEIGHTS,
    "liquid_wealth_to_income": 12.0,
    "housing_user_cost_share": 250.0,
    "prime_childless_renter_median_rooms": 10.0,
    "prime_childless_owner_median_rooms": 10.0,
}


CANDIDATE_REPLACEMENT_V1_WEIGHTS = {
    "tfr": 20.0,
    "childless_rate": 20.0,
    "own_rate": 100.0,
    "own_family_gap": 45.0,
    "housing_increment_0to1": 14.0,
    "housing_increment_1to2": 8.0,
    "young_liquid_wealth_to_income": 12.0,
    "old_nonhousing_wealth_to_income_6575": 0.8,
    "old_parent_childless_nonhousing_wealth_to_income_gap_6575": 2.0,
    "prime30_55_childless_renter_mean_rooms": 6.0,
    "prime30_55_childless_owner_mean_rooms": 6.0,
    "prime30_55_childless_renter_share_rooms_ge6": 25.0,
    "prime30_55_childless_owner_share_rooms_ge6": 25.0,
}


CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_WEIGHTS = {
    **{k: v for k, v in CANDIDATE_REPLACEMENT_V1_WEIGHTS.items() if k != "old_nonhousing_wealth_to_income_6575"},
    "old_nonhousing_wealth_to_income_median_6575": 0.8,
}


CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_V1_WEIGHTS = {
    **CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_WEIGHTS,
    "old_minus_young_owner_rate_6575_2534": 80.0,
}


CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_SOFT_V1_WEIGHTS = {
    **CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_WEIGHTS,
    "old_minus_young_owner_rate_6575_2534": 30.0,
}


CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_OWNGAP_V1_WEIGHTS = {
    **{
        k: v
        for k, v in CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_WEIGHTS.items()
        if k != "old_parent_childless_nonhousing_wealth_to_income_gap_6575"
    },
    "old_age_parent_childless_gap": 40.0,
    "old_minus_young_owner_rate_6575_2534": 80.0,
}


CANDIDATE_REPLACEMENT_OLD_RETENTION_V1_WEIGHTS = {
    **CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_WEIGHTS,
    "old_age_own_rate": 160.0,
    "old_age_parent_childless_gap": 40.0,
}


CANDIDATE_REPLACEMENT_YOUNG_OLD_OWN_V1_WEIGHTS = {
    **CANDIDATE_REPLACEMENT_OLD_RETENTION_V1_WEIGHTS,
    "own_rate_2534": 80.0,
}


CANDIDATE_REPLACEMENT_YOUNG_OLD_ROOMGAP_V1_WEIGHTS = {
    **CANDIDATE_REPLACEMENT_YOUNG_OLD_OWN_V1_WEIGHTS,
    "prime30_55_childless_owner_minus_renter_mean_rooms": 12.0,
}


CANDIDATE_REPLACEMENT_ROOMGAP_14MOMENT_V1_WEIGHTS = {
    **{
        k: v
        for k, v in CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_WEIGHTS.items()
        if k
        not in {
            "prime30_55_childless_owner_mean_rooms",
            "prime30_55_childless_renter_share_rooms_ge6",
            "young_liquid_wealth_to_income",
        }
    },
    "young_childless_renter_liquid_wealth_to_annual_gross_income_2535": 12.0,
    "prime30_55_childless_owner_minus_renter_mean_rooms": 12.0,
    "old_age_own_rate": 160.0,
    "own_rate_2534": 80.0,
}


CANDIDATE_REPLACEMENT_POST_AUDIT_V1_WEIGHTS = {
    **{
        k: v
        for k, v in CANDIDATE_REPLACEMENT_ROOMGAP_14MOMENT_V1_WEIGHTS.items()
        if k != "housing_increment_1to2"
    },
    "prime30_55_parent_3plus_minus_1to2_mean_rooms": CANDIDATE_REPLACEMENT_ROOMGAP_14MOMENT_V1_WEIGHTS[
        "housing_increment_1to2"
    ],
}

CANDIDATE_REPLACEMENT_POST_AUDIT_WEALTHSTRESS_V1_WEIGHTS = {
    **CANDIDATE_REPLACEMENT_POST_AUDIT_V1_WEIGHTS,
    "young_childless_renter_liquid_wealth_to_annual_gross_income_2535": 120.0,
}


CANDIDATE_REPLACEMENT_TOTAL_MEDIAN_V1_WEIGHTS = {
    **{
        k: v
        for k, v in CANDIDATE_REPLACEMENT_V1_WEIGHTS.items()
        if k
        not in {
            "old_nonhousing_wealth_to_income_6575",
            "old_parent_childless_nonhousing_wealth_to_income_gap_6575",
        }
    },
    "old_total_wealth_to_income_median_6575": 0.8,
    "old_parent_childless_total_wealth_to_income_gap_6575": 2.0,
}


CANDIDATE_REPLACEMENT_TOTAL_DUE_LIFECYCLE_V1_WEIGHTS = {
    **CANDIDATE_REPLACEMENT_TOTAL_MEDIAN_V1_WEIGHTS,
    "old_minus_young_owner_rate_6575_2534": 80.0,
}


CANDIDATE_REPLACEMENT_BEQUEST_INTERNAL_V1_WEIGHTS = {
    **{
        k: v
        for k, v in CANDIDATE_REPLACEMENT_POST_AUDIT_V1_WEIGHTS.items()
        if k
        not in {
            "old_nonhousing_wealth_to_income_median_6575",
            "old_parent_childless_nonhousing_wealth_to_income_gap_6575",
            "old_age_own_rate",
        }
    },
    # Diagonal inverse person-bootstrap variances from
    # bequest_calibration_targets.csv. The full covariance remains reported.
    "old_total_estate_wealth_to_annual_income_median_7684": 18.585767349158665,
    "old_total_estate_wealth_to_annual_income_p90_p50_7684": 56.97941874555749,
    "old_2plus_minus_1_total_estate_wealth_to_annual_income_median_gap_6575": 3.154769728751094,
}


CANDIDATE_REPLACEMENT_BEQUEST_MEDIAN_COMPOSITION_V1_WEIGHTS = {
    **{
        k: v
        for k, v in CANDIDATE_REPLACEMENT_BEQUEST_INTERNAL_V1_WEIGHTS.items()
        if k
        not in {
            "old_total_estate_wealth_to_annual_income_p90_p50_7684",
            "old_2plus_minus_1_total_estate_wealth_to_annual_income_median_gap_6575",
        }
    },
    "old_nonhousing_wealth_to_income_median_6575": NONHOUSING_MEDIAN_WEIGHT,
}


CANDIDATE_REPLACEMENT_INCOME_DISCIPLINED_V1_WEIGHTS = {
    **{
        k: v
        for k, v in CANDIDATE_REPLACEMENT_BEQUEST_MEDIAN_COMPOSITION_V1_WEIGHTS.items()
        if k != "old_nonhousing_wealth_to_income_median_6575"
    },
    "old_nonhousing_ge_1x_income_share_6575": OLD_NONHOUSING_GE1X_SHARE_WEIGHT,
    "old_age_own_rate": 160.0,
}


CORE_FEASIBILITY_V1_WEIGHTS = {
    "tfr": 20.0,
    "childless_rate": 20.0,
    "own_rate": 90.0,
    "own_family_gap": 45.0,
    "housing_increment_0to1": 14.0,
    "housing_increment_1to2": 8.0,
    "prime_childless_renter_median_rooms": 8.0,
    "prime_childless_owner_median_rooms": 8.0,
}


COST_TEST_V1_WEIGHTS = {
    **CORE_FEASIBILITY_V1_WEIGHTS,
    "housing_user_cost_share": 320.0,
}


OLDAGE_TEST_V1_WEIGHTS = {
    **CORE_FEASIBILITY_V1_WEIGHTS,
    "old_age_own_rate": 160.0,
    "old_age_parent_childless_gap": 40.0,
}


ROOMCOST_TEST_V1_WEIGHTS = {
    "tfr": 20.0,
    "childless_rate": 20.0,
    "own_rate": 100.0,
    "own_family_gap": 45.0,
    "old_age_own_rate": 160.0,
    "old_age_parent_childless_gap": 40.0,
    "housing_user_cost_share": 320.0,
    "housing_increment_0to1": 14.0,
    "housing_increment_1to2": 8.0,
}


CANDIDATE_NO_TIMING_OWNHEAVY_V1_WEIGHTS = {
    **CANDIDATE_NO_TIMING_V0_WEIGHTS,
    "own_rate": 120.0,
    "own_family_gap": 50.0,
    "old_age_own_rate": 80.0,
    "old_age_parent_childless_gap": 20.0,
    "housing_user_cost_share": 120.0,
    "liquid_wealth_to_income": 8.0,
    "young_liquid_wealth_to_income": 8.0,
    "housing_increment_0to1": 6.0,
    "prime_childless_renter_median_rooms": 5.0,
    "prime_childless_owner_median_rooms": 5.0,
}


CANDIDATE_NO_TIMING_REFINEMENT_V1_WEIGHTS = {
    **CANDIDATE_NO_TIMING_V0_WEIGHTS,
    "tfr": 20.0,
    "childless_rate": 20.0,
    "own_rate": 110.0,
    "own_family_gap": 45.0,
    "old_age_own_rate": 160.0,
    "old_age_parent_childless_gap": 40.0,
    "housing_user_cost_share": 320.0,
    "liquid_wealth_to_income": 20.0,
    "young_liquid_wealth_to_income": 14.0,
    "housing_increment_0to1": 14.0,
    "housing_increment_1to2": 8.0,
    "prime_childless_renter_median_rooms": 6.0,
    "prime_childless_owner_median_rooms": 6.0,
}


TARGET_SETS = {
    "core": (CORE_TARGETS, CORE_WEIGHTS),
    "old_nonlocation": (OLD_NONLOCATION_TARGETS, OLD_NONLOCATION_WEIGHTS),
    "old_nonlocation_no_timing": (OLD_NONLOCATION_NO_TIMING_TARGETS, OLD_NONLOCATION_NO_TIMING_WEIGHTS),
    "candidate_no_timing_v0": (CANDIDATE_NO_TIMING_V0_TARGETS, CANDIDATE_NO_TIMING_V0_WEIGHTS),
    "candidate_replacement_v1": (CANDIDATE_REPLACEMENT_V1_TARGETS, CANDIDATE_REPLACEMENT_V1_WEIGHTS),
    "candidate_replacement_nh_median_v1": (
        CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_TARGETS,
        CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_WEIGHTS,
    ),
    "candidate_replacement_due_lifecycle_v1": (
        CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_V1_TARGETS,
        CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_V1_WEIGHTS,
    ),
    "candidate_replacement_due_lifecycle_soft_v1": (
        CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_V1_TARGETS,
        CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_SOFT_V1_WEIGHTS,
    ),
    "candidate_replacement_due_lifecycle_owngap_v1": (
        CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_OWNGAP_V1_TARGETS,
        CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_OWNGAP_V1_WEIGHTS,
    ),
    "candidate_replacement_old_retention_v1": (
        CANDIDATE_REPLACEMENT_OLD_RETENTION_V1_TARGETS,
        CANDIDATE_REPLACEMENT_OLD_RETENTION_V1_WEIGHTS,
    ),
    "candidate_replacement_young_old_own_v1": (
        CANDIDATE_REPLACEMENT_YOUNG_OLD_OWN_V1_TARGETS,
        CANDIDATE_REPLACEMENT_YOUNG_OLD_OWN_V1_WEIGHTS,
    ),
    "candidate_replacement_young_old_roomgap_v1": (
        CANDIDATE_REPLACEMENT_YOUNG_OLD_ROOMGAP_V1_TARGETS,
        CANDIDATE_REPLACEMENT_YOUNG_OLD_ROOMGAP_V1_WEIGHTS,
    ),
    "candidate_replacement_roomgap_14moment_v1": (
        CANDIDATE_REPLACEMENT_ROOMGAP_14MOMENT_V1_TARGETS,
        CANDIDATE_REPLACEMENT_ROOMGAP_14MOMENT_V1_WEIGHTS,
    ),
    "candidate_replacement_roomgap_14moment_tfr192_v1": (
        CANDIDATE_REPLACEMENT_ROOMGAP_14MOMENT_TFR192_V1_TARGETS,
        CANDIDATE_REPLACEMENT_ROOMGAP_14MOMENT_V1_WEIGHTS,
    ),
    "candidate_replacement_post_audit_v1": (
        CANDIDATE_REPLACEMENT_POST_AUDIT_V1_TARGETS,
        CANDIDATE_REPLACEMENT_POST_AUDIT_V1_WEIGHTS,
    ),
    "candidate_replacement_post_audit_wealthstress_v1": (
        CANDIDATE_REPLACEMENT_POST_AUDIT_V1_TARGETS,
        CANDIDATE_REPLACEMENT_POST_AUDIT_WEALTHSTRESS_V1_WEIGHTS,
    ),
    "candidate_replacement_total_median_v1": (
        CANDIDATE_REPLACEMENT_TOTAL_MEDIAN_V1_TARGETS,
        CANDIDATE_REPLACEMENT_TOTAL_MEDIAN_V1_WEIGHTS,
    ),
    "candidate_replacement_total_due_lifecycle_v1": (
        CANDIDATE_REPLACEMENT_TOTAL_DUE_LIFECYCLE_V1_TARGETS,
        CANDIDATE_REPLACEMENT_TOTAL_DUE_LIFECYCLE_V1_WEIGHTS,
    ),
    "candidate_replacement_bequest_internal_v1": (
        CANDIDATE_REPLACEMENT_BEQUEST_INTERNAL_V1_TARGETS,
        CANDIDATE_REPLACEMENT_BEQUEST_INTERNAL_V1_WEIGHTS,
    ),
    "candidate_replacement_bequest_median_composition_v1": (
        CANDIDATE_REPLACEMENT_BEQUEST_MEDIAN_COMPOSITION_V1_TARGETS,
        CANDIDATE_REPLACEMENT_BEQUEST_MEDIAN_COMPOSITION_V1_WEIGHTS,
    ),
    "candidate_replacement_income_disciplined_v1": (
        CANDIDATE_REPLACEMENT_INCOME_DISCIPLINED_V1_TARGETS,
        CANDIDATE_REPLACEMENT_INCOME_DISCIPLINED_V1_WEIGHTS,
    ),
    "candidate_no_timing_core_feasibility_v1": (
        CORE_FEASIBILITY_V1_TARGETS,
        CORE_FEASIBILITY_V1_WEIGHTS,
    ),
    "candidate_no_timing_cost_test_v1": (
        COST_TEST_V1_TARGETS,
        COST_TEST_V1_WEIGHTS,
    ),
    "candidate_no_timing_oldage_test_v1": (
        OLDAGE_TEST_V1_TARGETS,
        OLDAGE_TEST_V1_WEIGHTS,
    ),
    "candidate_no_timing_roomcost_test_v1": (
        ROOMCOST_TEST_V1_TARGETS,
        ROOMCOST_TEST_V1_WEIGHTS,
    ),
    "candidate_no_timing_ownheavy_v1": (
        CANDIDATE_NO_TIMING_V0_TARGETS,
        CANDIDATE_NO_TIMING_OWNHEAVY_V1_WEIGHTS,
    ),
    "candidate_no_timing_refinement_v1": (
        CANDIDATE_NO_TIMING_V0_TARGETS,
        CANDIDATE_NO_TIMING_REFINEMENT_V1_WEIGHTS,
    ),
}


def run_small_calibration(
    outdir: Path,
    *,
    n_cases: int = 24,
    seed: int = 1234,
    J: int = 12,
    Nb: int = 40,
    n_house: int = 4,
    max_iter_eq: int = 35,
    progress: bool = True,
    target_set: str = "old_nonlocation",
) -> dict[str, Any]:
    """Run a small random-search diagnostic calibration.

    This is intentionally not a formal SMM routine. Its purpose is to find a
    locally more interpretable parameter point for inspecting the scaffold.
    """

    rng = np.random.default_rng(seed)
    targets, weights = get_target_set(target_set)
    outdir.mkdir(parents=True, exist_ok=True)
    cases_path = outdir / "cases.jsonl"
    best_path = outdir / "best.json"
    meta = {
        "target_set": str(target_set),
        "targets": targets,
        "weights": weights,
        "n_cases": int(n_cases),
        "seed": int(seed),
        "J": int(J),
        "period_years": PERIOD_YEARS,
        "age_start": AGE_START,
        "fertility_start_age": FERTILITY_START_AGE,
        "fertility_end_age": FERTILITY_END_AGE,
        "retirement_age": RETIREMENT_AGE,
        "child_maturity_age": CHILD_MATURITY_AGE,
        "stage_durations": [CHILD_MATURITY_AGE / PERIOD_YEARS],
        "expected_child_years": [CHILD_MATURITY_AGE],
        "Nb": int(Nb),
        "n_house": int(n_house),
        "max_iter_eq": int(max_iter_eq),
        "progress": bool(progress),
        "status": "diagnostic_only_not_formal_smm",
    }
    (outdir / "metadata.json").write_text(json.dumps(meta, indent=2, sort_keys=True))

    base = base_overrides(J=J, Nb=Nb, n_house=n_house, max_iter_eq=max_iter_eq)
    best: dict[str, Any] | None = None
    start = time.perf_counter()
    cases_path.write_text("")
    for idx in range(n_cases):
        theta = draw_candidate(rng, idx)
        overrides = {**base, **theta}
        t0 = time.perf_counter()
        try:
            sol, P, p_eq = run_model_cp_dt(overrides, verbose=False)
            moments = extract_moments(sol, P)
            loss = diagnostic_loss(moments, targets=targets, weights=weights)
            status = "ok"
            err = float(getattr(sol, "best_max_abs_rel_excess", np.nan))
            timings = getattr(sol, "timings", {})
        except Exception as exc:  # noqa: BLE001 - calibration should checkpoint failures.
            moments = {}
            loss = math.inf
            status = f"failed: {type(exc).__name__}: {exc}"
            p_eq = np.array([np.nan])
            err = math.inf
            timings = {}
        elapsed = time.perf_counter() - t0
        record = {
            "case": int(idx),
            "status": status,
            "loss": float(loss),
            "theta": jsonable(theta),
            "moments": jsonable(moments),
            "p_eq": jsonable(p_eq),
            "market_residual": float(err),
            "elapsed_sec": float(elapsed),
            "timings": jsonable(timings),
        }
        with cases_path.open("a") as fh:
            fh.write(json.dumps(record, sort_keys=True) + "\n")
        if best is None or record["loss"] < best["loss"]:
            best = record
            best_path.write_text(json.dumps(best, indent=2, sort_keys=True))
        if progress:
            best_loss = float(best["loss"]) if best is not None else math.inf
            print(
                f"case {idx + 1}/{n_cases}: loss={record['loss']:.4g}, "
                f"resid={record['market_residual']:.2e}, best={best_loss:.4g}, "
                f"elapsed={elapsed:.1f}s",
                flush=True,
            )

    summary = {
        "best": best,
        "elapsed_sec": float(time.perf_counter() - start),
        "metadata": meta,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    return summary


def run_informed_smoke(
    outdir: Path,
    *,
    J: int = 16,
    Nb: int = 40,
    n_house: int = 6,
    max_iter_eq: int = 25,
    labels: list[str] | None = None,
    case_limit: int | None = None,
    fixed_price: float | None = None,
    progress: bool = True,
) -> dict[str, Any]:
    """Run a deterministic production-parameter smoke panel.

    This is a parameter-led smoke check, not a formal calibration. It varies
    economically meaningful internal blocks around the current scaffold and
    keeps finance, supply, and menu objects fixed at their benchmark inputs.
    """

    targets, weights = get_target_set("old_nonlocation")
    rank_targets = {k: v for k, v in targets.items() if k != "mean_age_first_birth"}
    rank_weights = {k: v for k, v in weights.items() if k != "mean_age_first_birth"}
    outdir.mkdir(parents=True, exist_ok=True)
    cases_path = outdir / "cases.jsonl"
    best_path = outdir / "best.json"
    base = base_overrides(J=J, Nb=Nb, n_house=n_house, max_iter_eq=max_iter_eq)
    if fixed_price is not None:
        base.update(
            {
                "solve_mode": "pe",
                "p_fixed": np.array([float(fixed_price)]),
                "w_fixed": np.array([1.0]),
                "entry_shares_fixed": np.array([1.0]),
            }
        )
    candidates = informed_smoke_candidates()
    if labels:
        wanted = {str(label) for label in labels}
        candidates = [candidate for candidate in candidates if str(candidate["label"]) in wanted]
        missing = sorted(wanted - {str(candidate["label"]) for candidate in candidates})
        if missing:
            raise ValueError(f"Unknown informed-smoke labels: {missing}")
    if case_limit is not None:
        candidates = candidates[: max(0, int(case_limit))]
    meta = {
        "status": "deterministic_informed_smoke_not_formal_calibration",
        "case_count": len(candidates),
        "labels": [str(candidate["label"]) for candidate in candidates],
        "fixed_price": None if fixed_price is None else float(fixed_price),
        "J": int(J),
        "Nb": int(Nb),
        "n_house": int(n_house),
        "max_iter_eq": int(max_iter_eq),
        "fixed_external_or_first_stage": [
            "entry_wealth_mode",
            "entry_wealth_ratio_nodes",
            "entry_wealth_ratio_weights",
            "q",
            "delta",
            "tau_H",
            "phi",
            "pti_limit",
            "psi",
            "income_age_profile",
            "z_process",
            "H_own",
            "hR_max",
            "H0",
            "r_bar",
            "eta_supply",
            "theta1",
        ],
        "varied_blocks": [
            "beta",
            "alpha_cons",
            "c_bar_0",
            "c_bar_n",
            "h_bar_0",
            "h_bar_jump",
            "h_bar_n",
            "psi_child",
            "kappa_fert",
            "chi",
            "theta0",
            "theta_n",
        ],
        "missing_internal_parameter": "lambda_a_first_birth_timing_shifter",
        "ranking_loss": "old_nonlocation_excluding_mean_age_first_birth",
        "targets": targets,
        "weights": weights,
        "ranking_targets": rank_targets,
        "ranking_weights": rank_weights,
    }
    (outdir / "metadata.json").write_text(json.dumps(meta, indent=2, sort_keys=True))
    cases_path.write_text("")
    best: dict[str, Any] | None = None
    start = time.perf_counter()
    for idx, candidate in enumerate(candidates):
        label = str(candidate["label"])
        theta = dict(candidate["theta"])
        overrides = {**base, **theta}
        t0 = time.perf_counter()
        try:
            sol, P, p_eq = run_model_cp_dt(overrides, verbose=False)
            moments = extract_moments(sol, P)
            full_loss = diagnostic_loss(moments, targets=targets, weights=weights)
            rank_loss = diagnostic_loss(moments, targets=rank_targets, weights=rank_weights)
            status = "ok"
            err = float(getattr(sol, "best_max_abs_rel_excess", np.nan))
            timings = getattr(sol, "timings", {})
        except Exception as exc:  # noqa: BLE001 - smoke panel should checkpoint failures.
            moments = {}
            full_loss = math.inf
            rank_loss = math.inf
            status = f"failed: {type(exc).__name__}: {exc}"
            p_eq = np.array([np.nan])
            err = math.inf
            timings = {}
        elapsed = time.perf_counter() - t0
        record = {
            "case": int(idx),
            "label": label,
            "status": status,
            "rank_loss": float(rank_loss),
            "full_old_nonlocation_loss": float(full_loss),
            "theta": jsonable(theta),
            "moments": jsonable(moments),
            "p_eq": jsonable(p_eq),
            "market_residual": float(err),
            "elapsed_sec": float(elapsed),
            "timings": jsonable(timings),
        }
        with cases_path.open("a") as fh:
            fh.write(json.dumps(record, sort_keys=True) + "\n")
        if best is None or record["rank_loss"] < best["rank_loss"]:
            best = record
            best_path.write_text(json.dumps(best, indent=2, sort_keys=True))
        if progress:
            print(
                f"case {idx + 1}/{len(candidates)} {label}: "
                f"rank={record['rank_loss']:.4g}, full={record['full_old_nonlocation_loss']:.4g}, "
                f"resid={record['market_residual']:.2e}, elapsed={elapsed:.1f}s",
                flush=True,
            )

    summary = {
        "best": best,
        "elapsed_sec": float(time.perf_counter() - start),
        "metadata": meta,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    return summary


def informed_smoke_candidates() -> list[dict[str, Any]]:
    beta_hi = 0.98**PERIOD_YEARS
    beta_lo = 0.94**PERIOD_YEARS
    return [
        {"label": "baseline", "theta": {}},
        {"label": "patient_beta", "theta": {"beta": beta_hi}},
        {"label": "impatient_beta", "theta": {"beta": beta_lo}},
        {
            "label": "more_housing_share",
            "theta": {"alpha_cons": 0.60, "h_bar_0": 3.5, "c_bar_0": 0.08 * PERIOD_YEARS},
        },
        {
            "label": "less_housing_share",
            "theta": {"alpha_cons": 0.78, "h_bar_0": 4.5, "c_bar_0": 0.12 * PERIOD_YEARS},
        },
        {
            "label": "higher_family_space_need",
            "theta": {"h_bar_jump": 1.10, "h_bar_n": 0.90, "h_bar_0": 3.6},
        },
        {
            "label": "lower_child_cost_higher_taste",
            "theta": {"c_bar_n": 0.32, "psi_child": 0.10, "kappa_fert": 5.5},
        },
        {
            "label": "higher_child_cost_lower_taste",
            "theta": {"c_bar_n": 0.72, "psi_child": 0.045, "kappa_fert": 4.0},
        },
        {"label": "higher_owner_premium", "theta": {"chi": 1.30}},
        {
            "label": "stronger_bequest_parent_tilt",
            "theta": {"theta0": 0.80, "theta_n": 0.45},
        },
        {
            "label": "combined_internal_push",
            "theta": {
                "beta": beta_hi,
                "alpha_cons": 0.62,
                "c_bar_0": 0.08 * PERIOD_YEARS,
                "c_bar_n": 0.40,
                "h_bar_0": 3.5,
                "h_bar_jump": 1.00,
                "h_bar_n": 0.80,
                "psi_child": 0.09,
                "kappa_fert": 5.5,
                "chi": 1.25,
                "theta0": 0.75,
                "theta_n": 0.40,
            },
        },
    ]


def base_overrides(*, J: int, Nb: int, n_house: int, max_iter_eq: int) -> dict[str, Any]:
    lifecycle = lifecycle_overrides(J)
    return {
        "I": 1,
        "w_hat": np.array([1.0]),
        "E_loc": np.array([0.0]),
        "N_0": np.array([1.0]),
        "entry_shares": np.array([1.0]),
        "entry_by_loc": np.array([1.0 / J]),
        "r_bar": np.array([0.16]),
        "H0": np.array([4.0]),
        "eta_supply": np.array([1.0]),
        **lifecycle,
        "J": int(J),
        "Nb": int(Nb),
        "n_house": int(n_house),
        "H_own": np.linspace(2.0, 10.0, int(n_house)),
        "hR_max": 6.0,
        "max_iter_eq": int(max_iter_eq),
        "tol_eq": 1e-4,
        "use_pti_constraint": False,
        "scalar_market_refine": True,
        "scalar_market_refine_method": "brent",
        "scalar_market_refine_iter": 16,
        "scalar_market_refine_max_expand": 8,
        **external_entry_wealth_overrides(),
    }


def lifecycle_overrides(J: int) -> dict[str, Any]:
    j = int(J)
    period_years = PERIOD_YEARS
    fertile_start = age_to_period_number(FERTILITY_START_AGE, period_years)
    fertile_end = age_to_period_number(FERTILITY_END_AGE, period_years)
    retirement_idx = int(round((RETIREMENT_AGE - AGE_START) / period_years))
    return {
        "period_years": period_years,
        "da": period_years,
        "age_start": AGE_START,
        "A_m": CHILD_MATURITY_AGE,
        "stage_durations": np.array([CHILD_MATURITY_AGE / period_years]),
        "J_R": max(2, min(j - 1, retirement_idx)),
        "A_f_start": max(1, min(j, fertile_start)),
        "A_f_end": max(1, min(j, fertile_end)),
    }


def age_to_period_number(age: float, period_years: float) -> int:
    return int(round((float(age) - AGE_START) / period_years)) + 1


def draw_candidate(rng: np.random.Generator, idx: int) -> dict[str, Any]:
    if idx == 0:
        return {}
    if idx == 1:
        return {"hR_max": 4.0}
    if idx == 2:
        return {"hR_max": 5.0}
    if idx == 3:
        return {"phi": np.array([0.95, 0.95, 0.95])}
    if idx == 4:
        return {"hR_max": 4.0, "H0": np.array([3.0])}
    phi = rng.uniform(0.82, 0.97)
    return {
        "alpha_cons": rng.uniform(0.62, 0.82),
        "phi": np.array([phi, phi, phi]),
        "c_bar_n": rng.uniform(0.24, 0.72),
        "chi": rng.uniform(0.75, 1.50),
        "kappa_fert": rng.uniform(3.0, 7.0),
        "psi": rng.uniform(0.00, 0.08),
        "psi_child": rng.uniform(0.03, 0.12),
        "theta0": rng.uniform(0.12, 0.85),
        "theta_n": rng.uniform(0.00, 0.50),
        "theta1": rng.uniform(-0.05, 0.08),
        "h_bar_jump": rng.uniform(0.35, 1.10),
        "h_bar_n": rng.uniform(0.25, 0.90),
        "hR_max": rng.uniform(3.6, 6.5),
        "owner_h_bar_scale": rng.uniform(0.75, 1.15),
        "owner_size_cost": rng.uniform(0.0, 0.05),
        "owner_size_cost_ref": rng.uniform(2.0, 8.0),
        "H0": np.array([rng.uniform(3.0, 5.8)]),
    }


def extract_moments(sol: Any, P: Any | None = None) -> dict[str, Any]:
    household_parity = float(getattr(sol, "mean_completed_fertility", np.nan))
    parity = np.asarray(getattr(sol, "parity_dist", np.array([np.nan])), dtype=float).reshape(-1)
    renter_rooms = float(getattr(sol, "prime_childless_renter_median_rooms", np.nan))
    owner_rooms = float(getattr(sol, "prime_childless_owner_median_rooms", np.nan))
    own_rate_2534 = float(getattr(sol, "own_rate_2534", np.nan))
    old_age_own_rate = float(getattr(sol, "old_age_own_rate_6575", np.nan))
    total_mass = float(getattr(sol, "total_mass", 1.0))
    return {
        "aggregate_mean_occupied_rooms_18_85": float(
            getattr(sol, "aggregate_housing_demand", np.nan)
        ) / max(total_mass, 1e-12),
        "tfr": 2.0 * household_parity,
        "own_rate": float(getattr(sol, "own_rate_3055", np.nan)),
        "aggregate_own_rate": float(getattr(sol, "own_rate", np.nan)),
        "young_owner_rate": float(getattr(sol, "young_owner_rate", np.nan)),
        "old_owner_rate": float(getattr(sol, "old_owner_rate", np.nan)),
        "old_age_own_rate": old_age_own_rate,
        "mean_completed_fertility": household_parity,
        "childless_rate": float(getattr(sol, "childless_rate", np.nan)),
        "own_rate_3055": float(getattr(sol, "own_rate_3055", np.nan)),
        "own_rate_2534": own_rate_2534,
        "old_minus_young_owner_rate_6575_2534": old_age_own_rate - own_rate_2534,
        "own_rate_3544": float(getattr(sol, "own_rate_3544", np.nan)),
        "own_rate_nonparents_3055": float(getattr(sol, "own_rate_nonparents_3055", np.nan)),
        "own_rate_newparents_3055": float(getattr(sol, "own_rate_newparents_3055", np.nan)),
        "own_gap_newparent_nonparent_3055": float(getattr(sol, "own_gap_newparent_nonparent_3055", np.nan)),
        "own_family_gap": float(
            getattr(sol, "own_gap_newparent_nonparent_3055", getattr(sol, "own_family_gap", np.nan))
        ),
        "old_age_parent_childless_gap_6575": float(getattr(sol, "old_age_parent_childless_gap_6575", np.nan)),
        "old_age_parent_childless_gap": float(getattr(sol, "old_age_parent_childless_gap_6575", np.nan)),
        "old_nonhousing_wealth_to_income_6575": float(getattr(sol, "old_nonhousing_wealth_to_income_6575", np.nan)),
        "old_total_wealth_to_income_6575": float(getattr(sol, "old_total_wealth_to_income_6575", np.nan)),
        "old_nonhousing_wealth_to_income_median_6575": float(
            getattr(sol, "old_nonhousing_wealth_to_income_median_6575", np.nan)
        ),
        "old_total_wealth_to_income_median_6575": float(
            getattr(sol, "old_total_wealth_to_income_median_6575", np.nan)
        ),
        "old_nonhousing_ge_1x_income_share_6575": float(
            getattr(sol, "old_nonhousing_ge_1x_income_share_6575", np.nan)
        ),
        "old_parent_childless_nonhousing_wealth_to_income_gap_6575": float(
            getattr(sol, "old_parent_childless_nonhousing_wealth_to_income_gap_6575", np.nan)
        ),
        "old_parent_childless_total_wealth_to_income_gap_6575": float(
            getattr(sol, "old_parent_childless_total_wealth_to_income_gap_6575", np.nan)
        ),
        "old_parent_childless_nonhousing_wealth_to_income_median_gap_6575": float(
            getattr(sol, "old_parent_childless_nonhousing_wealth_to_income_median_gap_6575", np.nan)
        ),
        "old_parent_childless_total_wealth_to_income_median_gap_6575": float(
            getattr(sol, "old_parent_childless_total_wealth_to_income_median_gap_6575", np.nan)
        ),
        "old_total_estate_wealth_to_annual_income_median_7684": float(
            getattr(sol, "old_total_estate_wealth_to_annual_income_median_7684", np.nan)
        ),
        "old_total_estate_wealth_to_annual_income_p90_p50_7684": float(
            getattr(sol, "old_total_estate_wealth_to_annual_income_p90_p50_7684", np.nan)
        ),
        "old_2plus_minus_1_total_estate_wealth_to_annual_income_median_gap_6575": float(
            getattr(
                sol,
                "old_2plus_minus_1_total_estate_wealth_to_annual_income_median_gap_6575",
                np.nan,
            )
        ),
        "mean_age_first_birth": float(getattr(sol, "mean_age_first_birth", np.nan)),
        "attempt_hazard_by_age": np.asarray(getattr(sol, "attempt_hazard_by_age", np.array([np.nan])), dtype=float),
        "first_birth_hazard_by_age": np.asarray(getattr(sol, "first_birth_hazard_by_age", np.array([np.nan])), dtype=float),
        "first_birth_age_distribution": np.asarray(
            getattr(sol, "first_birth_age_distribution", np.array([np.nan])), dtype=float
        ),
        "share_first_births_age30plus": float(getattr(sol, "share_first_births_age30plus", np.nan)),
        "childless_chosen_45": float(getattr(sol, "childless_chosen_45", np.nan)),
        "childless_clock_45": float(getattr(sol, "childless_clock_45", np.nan)),
        "housing_increment_0to1": float(getattr(sol, "housing_increment_0to1_eventstudy_t3", np.nan)),
        "housing_increment_1to2": float(
            getattr(sol, "housing_increment_1to2_proxy_t3", getattr(sol, "housing_increment_1to2", np.nan))
        ),
        "prime30_55_parent_3plus_minus_1to2_mean_rooms": float(
            getattr(sol, "prime30_55_parent_3plus_minus_1to2_mean_rooms", np.nan)
        ),
        "young_liquid_wealth_to_income": float(getattr(sol, "young_liquid_wealth_to_income", np.nan)),
        "liquid_wealth_to_income": float(getattr(sol, "liquid_wealth_to_income", np.nan)),
        "wealth_to_income": float(getattr(sol, "wealth_to_income", np.nan)),
        "young_all_liquid_wealth_to_annual_gross_income_2530": float(
            getattr(sol, "young_all_liquid_wealth_to_annual_gross_income_2530", np.nan)
        ),
        "young_all_liquid_wealth_to_annual_gross_income_2530_median": float(
            getattr(sol, "young_all_liquid_wealth_to_annual_gross_income_2530_median", np.nan)
        ),
        "young_childless_liquid_wealth_to_annual_gross_income_2535": float(
            getattr(sol, "young_childless_liquid_wealth_to_annual_gross_income_2535", np.nan)
        ),
        "young_childless_liquid_wealth_to_annual_gross_income_2535_median": float(
            getattr(sol, "young_childless_liquid_wealth_to_annual_gross_income_2535_median", np.nan)
        ),
        "young_childless_renter_liquid_wealth_to_annual_gross_income_2535": float(
            getattr(sol, "young_childless_renter_liquid_wealth_to_annual_gross_income_2535", np.nan)
        ),
        "young_childless_renter_liquid_wealth_to_annual_gross_income_2535_median": float(
            getattr(sol, "young_childless_renter_liquid_wealth_to_annual_gross_income_2535_median", np.nan)
        ),
        "parity_share_0": float(parity[0]) if parity.size > 0 else np.nan,
        "parity_share_1": float(parity[1]) if parity.size > 1 else np.nan,
        "parity_share_2plus": float(np.sum(parity[2:])) if parity.size > 2 else np.nan,
        "prime_childless_renter_median_rooms": renter_rooms,
        "prime_childless_owner_median_rooms": owner_rooms,
        "prime_childless_owner_minus_renter_rooms": owner_rooms - renter_rooms,
        "prime30_55_childless_renter_mean_rooms": float(
            getattr(sol, "prime30_55_childless_renter_mean_rooms", np.nan)
        ),
        "prime30_55_childless_owner_mean_rooms": float(
            getattr(sol, "prime30_55_childless_owner_mean_rooms", np.nan)
        ),
        "prime30_55_childless_owner_minus_renter_mean_rooms": float(
            getattr(sol, "prime30_55_childless_owner_minus_renter_mean_rooms", np.nan)
        ),
        "prime30_55_childless_renter_share_rooms_ge6": float(
            getattr(sol, "prime30_55_childless_renter_share_rooms_ge6", np.nan)
        ),
        "prime30_55_childless_owner_share_rooms_ge6": float(
            getattr(sol, "prime30_55_childless_owner_share_rooms_ge6", np.nan)
        ),
        "prime30_55_parent_owner_minus_renter_mean_rooms": float(
            getattr(sol, "prime30_55_parent_owner_minus_renter_mean_rooms", np.nan)
        ),
        "housing_user_cost_share": housing_user_cost_share(sol, P),
        "renter25_45_all_cap_share": float(getattr(sol, "renter25_45_all_cap_share", np.nan)),
        "owner_neg_liquid_share_2534": float(getattr(sol, "owner_neg_liquid_share_2534", np.nan)),
        "market_residual": float(getattr(sol, "best_max_abs_rel_excess", np.nan)),
    }


def housing_user_cost_share(sol: Any, P: Any | None) -> float:
    if P is None or not hasattr(sol, "g") or not hasattr(sol, "hR_pol"):
        return np.nan
    g = np.asarray(sol.g, dtype=float)
    hR = np.asarray(sol.hR_pol, dtype=float)
    user_cost = np.asarray(getattr(sol, "owner_user_cost", np.nan), dtype=float).reshape(-1)
    if g.ndim != 7 or hR.shape != g.shape or user_cost.size != int(P.I):
        return np.nan
    z_grid = np.asarray(getattr(P, "z_grid", [1.0]), dtype=float).reshape(-1)
    total_income = 0.0
    total_housing = 0.0
    for i in range(P.I):
        for j in range(P.J):
            for zz, z_value in enumerate(z_grid):
                yj = income_at_state(P, i, j, float(z_value))
                mass_ijz = float(np.sum(g[:, :, i, j, zz, :, :]))
                total_income += yj * mass_ijz
                total_housing += float(np.sum(g[:, 0, i, j, zz, :, :] * hR[:, 0, i, j, zz, :, :])) * user_cost[i]
                for ten in range(1, 1 + int(P.n_house)):
                    mass = float(np.sum(g[:, ten, i, j, zz, :, :]))
                    total_housing += mass * user_cost[i] * float(P.H_own[ten - 1])
    return float(total_housing / max(total_income, 1e-12))


def get_target_set(target_set: str) -> tuple[dict[str, float], dict[str, float]]:
    key = str(target_set).lower()
    if key not in TARGET_SETS:
        raise ValueError(f"Unknown target_set {target_set!r}; choose from {sorted(TARGET_SETS)}")
    targets, weights = TARGET_SETS[key]
    return dict(targets), dict(weights)


def get_target_moment_objects(target_set: str | None = None) -> dict[str, dict[str, str]]:
    """Return documented model/data objects for target moments.

    This helper is deliberately non-behavioral: it does not affect the loss or
    target values. It keeps the measurement ledger available to calibration
    reports and future target revisions.
    """
    if target_set is None:
        return {name: dict(meta) for name, meta in TARGET_MOMENT_OBJECTS.items()}
    targets, _ = get_target_set(target_set)
    missing = {
        "model": "undocumented",
        "data": "undocumented",
        "status": "needs-documentation",
        "issue": "Add this moment to TARGET_MOMENT_OBJECTS before treating the target set as production.",
    }
    return {name: dict(TARGET_MOMENT_OBJECTS.get(name, missing)) for name in targets}


def diagnostic_loss(
    moments: dict[str, float],
    *,
    targets: dict[str, float] | None = None,
    weights: dict[str, float] | None = None,
) -> float:
    if targets is None or weights is None:
        targets, weights = get_target_set("core")
    loss = 0.0
    for name, target in targets.items():
        value = float(moments.get(name, np.nan))
        if not np.isfinite(value):
            return math.inf
        loss += float(weights.get(name, 1.0)) * (value - target) ** 2
    residual = float(moments.get("market_residual", np.nan))
    if not np.isfinite(residual) or residual > 5e-3:
        loss += 100.0
    return float(loss)


def jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value
