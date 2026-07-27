#!/usr/bin/env python3
"""Build CPS childlessness and completed-fertility gradients for women 40--44."""

from __future__ import annotations

import csv
import math
from itertools import combinations
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
SOURCE_PATH = REPO_ROOT / "code" / "data" / "cps_fertility" / "cache" / "jun24pub.csv"
OUTPUT_PATH = SCRIPT_DIR / "output" / "childlessness_gradient.csv"

EDUCATION_GROUPS = (
    ("less_than_high_school", range(31, 38)),
    ("high_school", range(38, 40)),
    ("some_college_associate", range(40, 43)),
    ("bachelor", range(43, 44)),
    ("advanced", range(44, 47)),
)


def weighted_quantile_code_bins(frame: pd.DataFrame) -> list[list[int]]:
    """Partition ordered HEFAMINC codes into four contiguous, near-equal weighted bins."""
    code_weights = frame.groupby("HEFAMINC", sort=True)["PWSSWGT"].sum()
    codes = [int(code) for code in code_weights.index]
    weights = [float(code_weights.loc[code]) for code in codes]
    if len(codes) < 4:
        raise RuntimeError(f"Fewer than four observed HEFAMINC codes: {codes}")
    total = sum(weights)
    cumulative = [0.0]
    for weight in weights:
        cumulative.append(cumulative[-1] + weight)
    # There are only 16 observed codes in this CPS file. Exhaustively select
    # the three cut points that minimize squared deviations from equal weighted
    # mass, preserving the ordinal and categorical nature of HEFAMINC.
    best_cuts: tuple[int, int, int] | None = None
    best_loss = math.inf
    for cuts in combinations(range(1, len(codes)), 3):
        boundaries = (0, *cuts, len(codes))
        masses = [
            cumulative[boundaries[index + 1]] - cumulative[boundaries[index]]
            for index in range(4)
        ]
        loss = sum((mass - total / 4) ** 2 for mass in masses)
        if loss < best_loss:
            best_loss = loss
            best_cuts = cuts
    if best_cuts is None:
        raise RuntimeError("Could not identify HEFAMINC cut points.")
    boundaries = (0, *best_cuts, len(codes))
    bins = [
        codes[boundaries[index] : boundaries[index + 1]]
        for index in range(4)
    ]
    if len(bins) != 4 or any(not item for item in bins):
        raise RuntimeError(f"Could not construct four nonempty HEFAMINC bins: {bins}")
    return bins


def summary_row(grouping: str, group: str, codes: list[int], subset: pd.DataFrame) -> dict[str, object]:
    weight = subset["PWSSWGT"].sum()
    if len(subset) == 0 or not math.isfinite(weight) or weight <= 0:
        raise RuntimeError(f"No positive-weight valid observations for {grouping}: {group}")
    childless_share = (subset["PWSSWGT"] * (subset["PTSF1"] == 0)).sum() / weight
    completed_fertility = (subset["PWSSWGT"] * subset["PTSF1"]).sum() / weight
    return {
        "grouping": grouping,
        "group": group,
        "codes": ",".join(str(code) for code in codes),
        "unweighted_n": len(subset),
        "weighted_population": weight,
        "weighted_share_childless": childless_share,
        "approx_binomial_se_childless": math.sqrt(childless_share * (1 - childless_share) / len(subset)),
        "weighted_mean_children_ever_born": completed_fertility,
        "age_range": "40-44",
        "weight_variable": "PWSSWGT",
        "fertility_universe": "PTSF1 valid codes 0-5; approximate binomial SE ignores survey design",
    }


def main() -> None:
    if not SOURCE_PATH.exists():
        raise FileNotFoundError(f"Missing CPS cache: {SOURCE_PATH}")
    required = {"PESEX", "PRTAGE", "PTSF1", "PWSSWGT", "PEEDUCA", "HEFAMINC"}
    header = set(pd.read_csv(SOURCE_PATH, nrows=0, encoding="utf-8-sig").columns)
    missing = sorted(required - header)
    if missing:
        raise RuntimeError(f"Missing required CPS columns: {', '.join(missing)}")

    frame = pd.read_csv(SOURCE_PATH, usecols=sorted(required), encoding="utf-8-sig")
    for column in required:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame[
        (frame["PESEX"] == 2)
        & frame["PRTAGE"].between(40, 44)
        & frame["PTSF1"].between(0, 5)
        & (frame["PWSSWGT"] > 0)
        & frame["PEEDUCA"].notna()
        & frame["HEFAMINC"].notna()
    ].copy()
    if frame.empty:
        raise RuntimeError("The CPS women 40-44 valid-fertility sample is empty.")
    frame["PEEDUCA"] = frame["PEEDUCA"].astype(int)
    frame["HEFAMINC"] = frame["HEFAMINC"].astype(int)
    frame["PTSF1"] = frame["PTSF1"].astype(int)

    observed_education = sorted(frame["PEEDUCA"].unique().tolist())
    expected_education = list(range(31, 47))
    if observed_education != expected_education:
        raise RuntimeError(
            "Observed PEEDUCA codes among valid women 40-44 do not match the requested "
            f"31--46 grouping. Observed code distribution: {frame['PEEDUCA'].value_counts().sort_index().to_dict()}"
        )

    rows: list[dict[str, object]] = []
    for label, code_range in EDUCATION_GROUPS:
        codes = list(code_range)
        rows.append(summary_row("education", label, codes, frame[frame["PEEDUCA"].isin(codes)]))

    income_bins = weighted_quantile_code_bins(frame)
    for index, codes in enumerate(income_bins, start=1):
        label = f"income_bin_{index}"
        rows.append(summary_row("family_income", label, codes, frame[frame["HEFAMINC"].isin(codes)]))

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    columns = list(rows[0])
    with OUTPUT_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    result = pd.DataFrame(rows)
    print("CPS childlessness and completed-fertility gradient")
    print(result[["grouping", "group", "codes", "unweighted_n", "weighted_share_childless", "approx_binomial_se_childless", "weighted_mean_children_ever_born"]].to_string(index=False))
    print(f"Wrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
