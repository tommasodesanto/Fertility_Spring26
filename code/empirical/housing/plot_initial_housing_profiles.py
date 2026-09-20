#!/usr/bin/env python3
"""Plot empirical age profiles from the completed housing sufficient statistics."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _profile_driver_hash() -> str | None:
    driver = Path(__file__).with_name("build_initial_housing_profile_diagnostic.py")
    return _sha(driver) if driver.exists() else None


def _build_model_overlay(source: Path, model_dir: Path, output_dir: Path) -> dict:
    emp = pd.read_csv(source)
    emp = emp[(emp["geography_scope"] == "active42") & (emp["sample"] == "all_structures") & (emp["age_kind"] == "four_year")].copy()
    emp = emp.groupby(["age_lower", "age_upper"], as_index=False)[["hhwt", "rooms_capped9_sum", "owner_hhwt"]].sum()
    emp["mean_rooms_capped9"] = emp["rooms_capped9_sum"] / emp["hhwt"]
    emp["ownership_rate"] = emp["owner_hhwt"] / emp["hhwt"]
    emp = emp[["age_lower", "age_upper", "mean_rooms_capped9", "ownership_rate"]]
    emp["source"] = "ACS 42 metros"
    frames = [emp]
    model_hashes = {}
    for path in sorted(model_dir.glob("*.json")):
        if path.name in {"launch_manifest.json", "submission.json"} or path.name.endswith("_summary.json"):
            continue
        d = json.loads(path.read_text())
        if d.get("status") != "completed":
            raise ValueError(f"model artifact not completed: {path}")
        p = d["profiles"]
        ages = [18 + 4 * i for i in p["age_cells"]]
        frames.append(pd.DataFrame({"age_lower": ages, "age_upper": [a + 3 for a in ages],
                                    "mean_rooms_capped9": p["mean_rooms_capped9_by_age"],
                                    "ownership_rate": p["ownership_by_age"], "source": d["family"]}))
        model_hashes[d["family"]] = _sha(path)
    series = pd.concat(frames, ignore_index=True)
    series.to_csv(output_dir / "model_overlay_series.csv", index=False)
    colors = {"ACS 42 metros": "#111111", "original": "#c8553d", "stationary_new_income": "#577590", "refit_new_income": "#43aa8b"}
    labels = {"ACS 42 metros": "ACS 42 metros", "original": "Original", "stationary_new_income": "New income", "refit_new_income": "Refit new income"}
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    for src in labels:
        sub = series[series.source == src]; ax.plot(sub.age_lower, sub.mean_rooms_capped9, marker="o", color=colors[src], label=labels[src])
    ax.set(xlabel="Physical age cell: 18 + 4j", ylabel="Mean rooms, capped at 9", title="Physical housing size by age"); ax.grid(alpha=.25); ax.legend(frameon=False)
    fig.savefig(output_dir / "model_overlay_mean_rooms_capped9.png", dpi=160); plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    for src in labels:
        sub = series[series.source == src]; ax.plot(sub.age_lower, 100 * sub.ownership_rate, marker="o", color=colors[src], label=labels[src])
    ax.set(xlabel="Physical age cell: 18 + 4j", ylabel="Ownership rate (%)", title="Ownership by age"); ax.grid(alpha=.25); ax.legend(frameon=False)
    fig.savefig(output_dir / "model_overlay_ownership.png", dpi=160); plt.close(fig)
    return {"model_dir": str(model_dir), "model_artifact_sha256": model_hashes,
            "model_measurement": "physical rooms: renter hR_pol and owner H_own, each capped at 9 before g_current weighting",
            "ownership_note": "model has no DUE classification; overlay is all-household/all-structure only",
            "comparison_note": "checkpoint family prices and preferences may differ; overlay is descriptive and not a causal effect"}


def build(input_dir: Path, output_dir: Path, model_dir: Path | None = None) -> None:
    output_dir.mkdir(parents=True, exist_ok=False)
    source = input_dir / "housing_profile_by_age.csv"
    provenance = json.loads((input_dir / "provenance.json").read_text())
    df = pd.read_csv(source)
    df = df[(df["age_kind"] == "four_year") & (df["sample"].isin(["all_structures", "DUE"]))].copy()
    group_cols = ["geography_scope", "sample", "age_lower", "age_upper"]
    sums = df.groupby(group_cols, as_index=False)[["hhwt", "rooms_capped9_sum", "owner_hhwt"]].sum()
    sums["mean_rooms_capped9"] = sums["rooms_capped9_sum"] / sums["hhwt"]
    sums["ownership_rate"] = sums["owner_hhwt"] / sums["hhwt"]
    sums.to_csv(output_dir / "initial_housing_profile_series.csv", index=False)

    colors = {"active42": "#c8553d", "national": "#577590"}
    scope_labels = {"active42": "42 metros", "national": "National"}
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    for scope in ("active42", "national"):
        sub = sums[(sums["geography_scope"] == scope) & (sums["sample"] == "all_structures")]
        ax.plot(sub["age_lower"], sub["mean_rooms_capped9"], marker="o", color=colors[scope], label=scope_labels[scope])
    ax.set(xlabel="Householder age, four-year cell lower bound", ylabel="Mean rooms, capped at 9", title="ACS housing size by age (all structures)")
    ax.legend(frameon=False); ax.grid(alpha=.25)
    fig.savefig(output_dir / "mean_rooms_capped9_by_age.png", dpi=160); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True, sharey=True)
    for jj, sample in enumerate(("all_structures", "DUE")):
        ax = axes[jj]
        for scope in ("active42", "national"):
            sub = sums[(sums["geography_scope"] == scope) & (sums["sample"] == sample)]
            ax.plot(sub["age_lower"], 100 * sub["ownership_rate"], marker="o", color=colors[scope], label=scope_labels[scope])
        sample_label = "All structures" if sample == "all_structures" else "DUE (UNITSSTR 3–10)"
        ax.set_title(f"Ownership by age: {sample_label}"); ax.set_xlabel("Age-cell lower bound"); ax.grid(alpha=.25)
    axes[0].set_ylabel("Ownership rate (%)"); axes[1].legend(frameon=False)
    fig.savefig(output_dir / "ownership_by_age_active42_vs_national.png", dpi=160); plt.close(fig)

    manifest = {"source_csv": str(source), "source_csv_sha256": _sha(source),
                "plot_helper_sha256": _sha(Path(__file__)),
                "input_provenance_source_sha256": provenance.get("source_sha256_from_existing_canonical_receipt"),
                "input_profile_driver_hash": _profile_driver_hash(),
                "input_reader_hash": provenance.get("reader_sha256"),
                "aggregation": "sum HHWT, capped rooms and owner HHWT across tenure/current-child cells; ratios after aggregation",
                "age_definition": "four-year bins beginning at 18; no equal-age weighting",
                "model_data_used": model_dir is not None}
    if model_dir is not None:
        manifest["model_overlay"] = _build_model_overlay(source, model_dir, output_dir)
    (output_dir / "plot_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, default=None)
    args = parser.parse_args()
    build(args.input, args.output, args.model_dir)


if __name__ == "__main__":
    main()
