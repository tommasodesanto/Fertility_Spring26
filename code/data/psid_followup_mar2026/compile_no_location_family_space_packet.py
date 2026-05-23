#!/usr/bin/env python3
"""Compile immediate no-location PSID family-space facts from existing outputs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output" / "no_location_family_space_packet_20260523"


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def weighted_mean(x: pd.Series, w: pd.Series) -> float:
    ok = x.notna() & w.notna() & (w > 0)
    if not ok.any():
        return float("nan")
    return float(np.average(x[ok], weights=w[ok]))


def by_group_summary(df: pd.DataFrame, group: str) -> pd.DataFrame:
    rows = []
    for value, g in df.groupby(group, dropna=False):
        w = g["iw_pre"] if "iw_pre" in g else pd.Series(np.ones(len(g)), index=g.index)
        rows.append(
            {
                group: value,
                "n": int(len(g)),
                "weight": float(w.sum(skipna=True)),
                "own_post3": weighted_mean(g["own_post3"], w),
                "moved_to_own_post3": weighted_mean(g["moved_to_own_post3"], w),
                "move_post3": weighted_mean(g["move_post3"], w),
                "rooms_pre": weighted_mean(g["rooms_pre"], w),
                "rooms_change_post3": weighted_mean(g["rooms_change_post3"], w),
                "incfamr_pre": weighted_mean(g["incfamr_pre"], w),
                "networth2r_pre": weighted_mean(g["networth2r_pre"], w),
            }
        )
    return pd.DataFrame(rows)


def fmt(x: float, digits: int = 3) -> str:
    if pd.isna(x):
        return "NA"
    return f"{x:.{digits}f}"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    rooms_first = read_csv(
        ROOT / "output" / "sa_rooms_first_birth_one_variant_v1" / "rooms_f_c_y_all_summary.csv"
    )
    rooms_second = read_csv(
        ROOT
        / "output"
        / "sa_rooms_second_birth_lasttreated_controls_v1"
        / "rooms_s_c_y_all_summary.csv"
    )
    moving_windows = read_csv(
        ROOT / "output" / "moving_first_birth_window_summary_v1" / "moving_first_birth_window_summary_v1.csv"
    )
    wealth_tercile_existing = read_csv(ROOT / "output" / "wealth_v1" / "wealth_tercile_outcomes_v1.csv")
    parenthood_quintile = read_csv(
        ROOT / "output" / "fertility_wealth_v1" / "parenthood_by_quintile_v12.csv"
    )

    event_path = ROOT / "output" / "wealth_v1" / "wealth_event_sample_pre_renters_v1.dta"
    event = pd.read_stata(event_path)
    event = event.copy()
    numeric_cols = [
        "incfamr_pre",
        "networth2r_pre",
        "iw_pre",
        "own_post3",
        "moved_to_own_post3",
        "move_post3",
        "rooms_pre",
        "rooms_change_post3",
        "nw2_tercile",
    ]
    for col in numeric_cols:
        event[col] = pd.to_numeric(event[col], errors="coerce")
    event = event[
        event["rooms_pre"].between(1, 12)
        & event["rooms_change_post3"].between(-10, 10)
        & event["incfamr_pre"].notna()
    ].copy()
    event["income_tercile"] = pd.qcut(event["incfamr_pre"], 3, labels=[1, 2, 3], duplicates="drop")
    event["rooms_slack_bin"] = np.select(
        [event["rooms_pre"] <= 4, event["rooms_pre"].between(5, 6), event["rooms_pre"] >= 7],
        ["low_slack_1_4", "mid_5_6", "high_7plus"],
        default=None,
    )
    event["wealth_tercile"] = event["nw2_tercile"]

    income_event = by_group_summary(event.dropna(subset=["income_tercile"]), "income_tercile")
    slack_event = by_group_summary(event.dropna(subset=["rooms_slack_bin"]), "rooms_slack_bin")
    wealth_event = by_group_summary(event.dropna(subset=["wealth_tercile"]), "wealth_tercile")
    income_slack_event = (
        event.dropna(subset=["income_tercile", "rooms_slack_bin"])
        .groupby(["income_tercile", "rooms_slack_bin"], dropna=False)
        .apply(
            lambda g: pd.Series(
                {
                    "n": int(len(g)),
                    "weight": float(g["iw_pre"].sum(skipna=True)),
                    "own_post3": weighted_mean(g["own_post3"], g["iw_pre"]),
                    "moved_to_own_post3": weighted_mean(g["moved_to_own_post3"], g["iw_pre"]),
                    "move_post3": weighted_mean(g["move_post3"], g["iw_pre"]),
                    "rooms_pre": weighted_mean(g["rooms_pre"], g["iw_pre"]),
                    "rooms_change_post3": weighted_mean(g["rooms_change_post3"], g["iw_pre"]),
                }
            )
        )
        .reset_index()
    )

    rooms_event_summary = pd.DataFrame(
        [
            {
                "source": "psid_first_birth_rooms",
                "pre_event_mean": float(rooms_first.loc[0, "pre_event_mean"]),
                "coef_p3": float(rooms_first.loc[0, "coef_p3"]),
                "se_p3": float(rooms_first.loc[0, "se_p3"]),
                "coef_p5": float(rooms_first.loc[0, "coef_p5"]),
                "se_p5": float(rooms_first.loc[0, "se_p5"]),
                "sample_obs": int(rooms_first.loc[0, "sample_obs"]),
                "sample_ids": int(rooms_first.loc[0, "sample_ids"]),
            },
            {
                "source": "psid_second_birth_rooms",
                "pre_event_mean": float(rooms_second.loc[0, "pre_event_mean"]),
                "coef_p3": float(rooms_second.loc[0, "coef_p3"]),
                "se_p3": float(rooms_second.loc[0, "se_p3"]),
                "coef_p5": float(rooms_second.loc[0, "coef_p5"]),
                "se_p5": float(rooms_second.loc[0, "se_p5"]),
                "sample_obs": int(rooms_second.loc[0, "sample_obs"]),
                "sample_ids": int(rooms_second.loc[0, "sample_ids"]),
            },
        ]
    )

    moving_key = moving_windows[moving_windows["window"].isin(["pre2", "post0to3", "post0to5"])].copy()

    rooms_event_summary.to_csv(OUT / "psid_birth_rooms_event_summary.csv", index=False)
    moving_key.to_csv(OUT / "psid_first_birth_moving_windows.csv", index=False)
    income_event.to_csv(OUT / "psid_prebirth_renter_outcomes_by_income_tercile.csv", index=False)
    slack_event.to_csv(OUT / "psid_prebirth_renter_outcomes_by_rooms_slack.csv", index=False)
    wealth_event.to_csv(OUT / "psid_prebirth_renter_outcomes_by_wealth_tercile.csv", index=False)
    income_slack_event.to_csv(OUT / "psid_prebirth_renter_outcomes_by_income_x_slack.csv", index=False)
    wealth_tercile_existing.to_csv(OUT / "psid_existing_wealth_tercile_outcomes.csv", index=False)
    parenthood_quintile.to_csv(OUT / "psid_parenthood_by_wealth_quintile.csv", index=False)

    first = rooms_event_summary.loc[rooms_event_summary["source"] == "psid_first_birth_rooms"].iloc[0]
    second = rooms_event_summary.loc[rooms_event_summary["source"] == "psid_second_birth_rooms"].iloc[0]
    post3_size = moving_key[(moving_key["outcome"] == "moved_for_size") & (moving_key["window"] == "post0to3")]
    pre2_size = moving_key[(moving_key["outcome"] == "moved_for_size") & (moving_key["window"] == "pre2")]
    low_inc = income_event.sort_values("income_tercile").iloc[0]
    high_inc = income_event.sort_values("income_tercile").iloc[-1]
    low_slack = slack_event[slack_event["rooms_slack_bin"] == "low_slack_1_4"].iloc[0]
    high_slack = slack_event[slack_event["rooms_slack_bin"] == "high_7plus"].iloc[0]

    lines = [
        "# PSID No-Location Family-Space Packet",
        "",
        "This packet compiles existing PSID outputs and a fresh split of the saved pre-birth-renter event sample. It does not use location-specific PSID.",
        "",
        "## Core Dynamic Housing Facts",
        "",
        f"- First birth raises rooms by `{fmt(first['coef_p3'])}` rooms by +3 years (SE `{fmt(first['se_p3'])}`); pre-event mean `{fmt(first['pre_event_mean'])}`.",
        f"- Second birth raises rooms by `{fmt(second['coef_p3'])}` rooms by +3 years (SE `{fmt(second['se_p3'])}`); pre-event mean `{fmt(second['pre_event_mean'])}`.",
        f"- Moved-for-size rate rises from pre-birth `t=-2` `{fmt(float(pre2_size['mean'].iloc[0]))}` to post-birth years 0-3 `{fmt(float(post3_size['mean'].iloc[0]))}` in the existing window summary.",
        "",
        "## Pre-Birth Renter Heterogeneity",
        "",
        f"- Top pre-birth income tercile: own by +3 `{fmt(high_inc['own_post3'])}`, moved-to-own by +3 `{fmt(high_inc['moved_to_own_post3'])}`, rooms change `{fmt(high_inc['rooms_change_post3'])}`.",
        f"- Bottom pre-birth income tercile: own by +3 `{fmt(low_inc['own_post3'])}`, moved-to-own by +3 `{fmt(low_inc['moved_to_own_post3'])}`, rooms change `{fmt(low_inc['rooms_change_post3'])}`.",
        f"- Low-slack renters (1-4 rooms pre-birth): own by +3 `{fmt(low_slack['own_post3'])}`, move by +3 `{fmt(low_slack['move_post3'])}`, rooms change `{fmt(low_slack['rooms_change_post3'])}`.",
        f"- High-slack renters (7+ rooms pre-birth): own by +3 `{fmt(high_slack['own_post3'])}`, move by +3 `{fmt(high_slack['move_post3'])}`, rooms change `{fmt(high_slack['rooms_change_post3'])}`.",
        "",
        "## Wealth-Fertility Composition Fact",
        "",
        "Young-adult wealth and parenthood are non-monotone in the existing PSID output: parenthood by age 35 across wealth quintiles is "
        + ", ".join(fmt(x) for x in parenthood_quintile["parent_by_35"].tolist())
        + ".",
        "",
        "## Outputs",
        "",
        "- `psid_birth_rooms_event_summary.csv`",
        "- `psid_first_birth_moving_windows.csv`",
        "- `psid_prebirth_renter_outcomes_by_income_tercile.csv`",
        "- `psid_prebirth_renter_outcomes_by_rooms_slack.csv`",
        "- `psid_prebirth_renter_outcomes_by_wealth_tercile.csv`",
        "- `psid_prebirth_renter_outcomes_by_income_x_slack.csv`",
        "- `psid_parenthood_by_wealth_quintile.csv`",
        "",
        "## Interpretation Guardrails",
        "",
        "These are no-location PSID facts. They support the family-space transition and heterogeneity by income, wealth, and pre-birth space slack. They do not identify center-to-periphery moves or metro family-size-premium effects.",
    ]
    (OUT / "PSID_NO_LOCATION_FAMILY_SPACE_PACKET.md").write_text("\n".join(lines) + "\n")
    print(f"Wrote PSID no-location family-space packet to {OUT}")


if __name__ == "__main__":
    main()
