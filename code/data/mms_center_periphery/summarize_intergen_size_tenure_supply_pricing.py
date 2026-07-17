#!/usr/bin/env python3
"""Summarize existing ACS/MMS size-by-tenure supply packets for intergen audits.

This script deliberately does not read the raw IPUMS extract. It aggregates
the precomputed ACS/MMS supply CSVs created by:

- analyze_family_size_supply_menu.R
- analyze_couillard_bedroom_supply_panel.R
- audit_intergen_quality_adjusted_room_targets.R
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = ROOT / "code" / "data" / "mms_center_periphery"
OUT_DIR = DATA_DIR / "output_intergen_size_tenure_supply_pricing"


def weighted_from_cells(
    df: pd.DataFrame,
    group_cols: list[str],
    weight_col: str,
    value_cols: list[str],
) -> pd.DataFrame:
    rows = []
    for keys, group in df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(group_cols, keys))
        row["weight"] = group[weight_col].sum()
        row["n_cells"] = len(group)
        for col in value_cols:
            valid = group[col].notna() & group[weight_col].notna() & (group[weight_col] > 0)
            if valid.any():
                denom = group.loc[valid, weight_col].sum()
                row[col] = (group.loc[valid, col] * group.loc[valid, weight_col]).sum() / denom
            else:
                row[col] = pd.NA
        rows.append(row)
    out = pd.DataFrame(rows)
    first = group_cols[:-1]
    if first:
        out["share_within_" + "_".join(first)] = out["weight"] / out.groupby(first)["weight"].transform("sum")
    return out


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def summarize_room_bins() -> pd.DataFrame:
    cells = read_csv(DATA_DIR / "output_family_size_supply" / "acs_family_size_supply_cells.csv")
    out = weighted_from_cells(
        cells,
        ["tenure", "size_bin"],
        "hh_weight",
        [
            "mean_rooms",
            "mean_rent",
            "median_rent",
            "mean_rent_per_room",
            "median_rent_per_room",
        ],
    )
    out["share_within_tenure"] = out["weight"] / out.groupby("tenure")["weight"].transform("sum")
    return out.sort_values(["tenure", "size_bin"])


def summarize_bedroom_bins() -> pd.DataFrame:
    cells = read_csv(DATA_DIR / "output_couillard_bedroom_supply" / "acs_bedroom_supply_cells.csv")
    out = weighted_from_cells(
        cells,
        ["tenure", "bedroom_bin"],
        "hh_weight",
        ["mean_bedrooms", "mean_rooms", "median_rent", "median_rent_per_bedroom"],
    )
    out["share_within_tenure"] = out["weight"] / out.groupby("tenure")["weight"].transform("sum")
    return out.sort_values(["tenure", "bedroom_bin"])


def summarize_structure() -> pd.DataFrame:
    path = (
        DATA_DIR
        / "output_intergen_quality_adjusted_room_targets"
        / "intergen_quality_adjusted_room_structure_summary.csv"
    )
    df = read_csv(path)
    return df.sort_values(["tenure", "weight_share_within_tenure", "unitsstr_label"], ascending=[True, False, True])


def summarize_structure_groups(structure: pd.DataFrame) -> pd.DataFrame:
    grouped = structure.copy()
    grouped["structure_group"] = "other"
    grouped.loc[grouped["unitsstr_label"].eq("1-family house, detached"), "structure_group"] = "single_family_detached"
    grouped.loc[grouped["unitsstr_label"].eq("1-family house, attached"), "structure_group"] = "single_family_attached"
    grouped.loc[grouped["unitsstr_label"].str.contains("family building", na=False), "structure_group"] = "multifamily_2plus"
    grouped.loc[grouped["unitsstr_label"].eq("Mobile home or trailer"), "structure_group"] = "mobile_home"

    rows = []
    for (tenure, structure_group), group in grouped.groupby(["tenure", "structure_group"]):
        weight = group["weight"].sum()
        row = {
            "tenure": tenure,
            "structure_group": structure_group,
            "weight": weight,
            "n": group["n"].sum(),
        }
        for col in ["mean_rooms", "mean_bedrooms", "mean_rent_to_income"]:
            valid = group[col].notna() & group["weight"].notna() & (group["weight"] > 0)
            if valid.any():
                row[col] = (group.loc[valid, col] * group.loc[valid, "weight"]).sum() / group.loc[valid, "weight"].sum()
            else:
                row[col] = pd.NA
        rows.append(row)
    out = pd.DataFrame(rows)
    out["share_within_tenure"] = out["weight"] / out.groupby("tenure")["weight"].transform("sum")
    return out.sort_values(["tenure", "share_within_tenure"], ascending=[True, False])


def pick(df: pd.DataFrame, **kwargs) -> pd.Series:
    mask = pd.Series(True, index=df.index)
    for key, value in kwargs.items():
        mask &= df[key].eq(value)
    if not mask.any():
        raise KeyError(kwargs)
    return df.loc[mask].iloc[0]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    room = summarize_room_bins()
    bedroom = summarize_bedroom_bins()
    structure = summarize_structure()
    structure_group = summarize_structure_groups(structure)

    room.to_csv(OUT_DIR / "room_bin_tenure_summary.csv", index=False)
    bedroom.to_csv(OUT_DIR / "bedroom_bin_tenure_summary.csv", index=False)
    structure.to_csv(OUT_DIR / "childless_30_55_structure_tenure_summary.csv", index=False)
    structure_group.to_csv(OUT_DIR / "childless_30_55_structure_group_tenure_summary.csv", index=False)

    owner_large = pick(room, tenure="owner", size_bin="L_7plus")
    renter_large = pick(room, tenure="renter", size_bin="L_7plus")
    owner_small = pick(room, tenure="owner", size_bin="S_1_4")
    renter_small = pick(room, tenure="renter", size_bin="S_1_4")
    renter_mid = pick(room, tenure="renter", size_bin="M_5_6")
    owner_b3 = pick(bedroom, tenure="owner", bedroom_bin="B_3plus")
    renter_b3 = pick(bedroom, tenure="renter", bedroom_bin="B_3plus")

    det_owner = pick(structure, tenure="owner", unitsstr_label="1-family house, detached")
    det_renter = pick(structure, tenure="renter", unitsstr_label="1-family house, detached")
    multi_owner = pick(structure_group, tenure="owner", structure_group="multifamily_2plus")
    multi_renter = pick(structure_group, tenure="renter", structure_group="multifamily_2plus")

    facts = {
        "source_packets": [
            "output_family_size_supply/acs_family_size_supply_cells.csv",
            "output_couillard_bedroom_supply/acs_bedroom_supply_cells.csv",
            "output_intergen_quality_adjusted_room_targets/intergen_quality_adjusted_room_structure_summary.csv",
        ],
        "room_bins_all_heads": {
            "owner_large_7plus_share": float(owner_large["share_within_tenure"]),
            "renter_large_7plus_share": float(renter_large["share_within_tenure"]),
            "owner_small_1to4_share": float(owner_small["share_within_tenure"]),
            "renter_small_1to4_share": float(renter_small["share_within_tenure"]),
            "renter_mid_5to6_share": float(renter_mid["share_within_tenure"]),
            "renter_mean_rent_1to4": float(renter_small["mean_rent"]),
            "renter_mean_rent_5to6": float(renter_mid["mean_rent"]),
            "renter_mean_rent_7plus": float(renter_large["mean_rent"]),
            "renter_mean_rent_per_room_1to4": float(renter_small["mean_rent_per_room"]),
            "renter_mean_rent_per_room_7plus": float(renter_large["mean_rent_per_room"]),
        },
        "bedroom_bins_all_heads": {
            "owner_3plus_bedroom_share": float(owner_b3["share_within_tenure"]),
            "renter_3plus_bedroom_share": float(renter_b3["share_within_tenure"]),
        },
        "structure_childless_heads_30_55": {
            "owner_detached_share": float(det_owner["weight_share_within_tenure"]),
            "renter_detached_share": float(det_renter["weight_share_within_tenure"]),
            "owner_detached_mean_rooms": float(det_owner["mean_rooms"]),
            "renter_detached_mean_rooms": float(det_renter["mean_rooms"]),
            "owner_multifamily_share": float(multi_owner["share_within_tenure"]),
            "renter_multifamily_share": float(multi_renter["share_within_tenure"]),
            "owner_multifamily_mean_rooms": float(multi_owner["mean_rooms"]),
            "renter_multifamily_mean_rooms": float(multi_renter["mean_rooms"]),
        },
    }
    (OUT_DIR / "key_facts.json").write_text(json.dumps(facts, indent=2) + "\n")

    readme = f"""# Intergen Size-Tenure Supply And Pricing Summary

Generated from existing ACS/MMS packet CSVs. This folder does not read the raw
IPUMS extract and does not change model/calibration code.

## Key Facts

- All heads: owner 7+ room share {facts['room_bins_all_heads']['owner_large_7plus_share']:.3f}; renter 7+ room share {facts['room_bins_all_heads']['renter_large_7plus_share']:.3f}.
- All heads: renter 1-4 room share {facts['room_bins_all_heads']['renter_small_1to4_share']:.3f}; owner 1-4 room share {facts['room_bins_all_heads']['owner_small_1to4_share']:.3f}.
- All-head renters: mean monthly rent is {facts['room_bins_all_heads']['renter_mean_rent_1to4']:.0f} in 1-4 room units, {facts['room_bins_all_heads']['renter_mean_rent_5to6']:.0f} in 5-6 room units, and {facts['room_bins_all_heads']['renter_mean_rent_7plus']:.0f} in 7+ room units.
- All-head renters: mean rent per room is {facts['room_bins_all_heads']['renter_mean_rent_per_room_1to4']:.0f} in 1-4 room units and {facts['room_bins_all_heads']['renter_mean_rent_per_room_7plus']:.0f} in 7+ room units.
- All heads: owner 3+ bedroom share {facts['bedroom_bins_all_heads']['owner_3plus_bedroom_share']:.3f}; renter 3+ bedroom share {facts['bedroom_bins_all_heads']['renter_3plus_bedroom_share']:.3f}.
- Childless heads 30-55: owner detached single-family share {facts['structure_childless_heads_30_55']['owner_detached_share']:.3f}; renter detached single-family share {facts['structure_childless_heads_30_55']['renter_detached_share']:.3f}.

## Outputs

- `room_bin_tenure_summary.csv`
- `bedroom_bin_tenure_summary.csv`
- `childless_30_55_structure_tenure_summary.csv`
- `childless_30_55_structure_group_tenure_summary.csv`
- `key_facts.json`

## Guardrail

These are descriptive stock/pricing facts. They support a model-menu/supply
diagnosis, but they do not by themselves identify a structural supply elasticity
or a final SMM target.
"""
    (OUT_DIR / "README.md").write_text(readme)


if __name__ == "__main__":
    main()
