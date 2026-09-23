#!/usr/bin/env python3
"""Reproduce available saved PSID room-estimand receipts; run no regression."""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[7]
DATA = ROOT / "code/data/psid_followup_mar2026/output"
MAIN = DATA / "sa_rooms_first_birth_household_aligned_v1"
REVIEW = DATA / "first_birth_correction_review"
OUT = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def rows(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


receipt = rows(MAIN / "target_receipt.csv")[0]
events = {int(r["relative_time"]): r for r in rows(MAIN / "event_study_estimates.csv")}
beta3 = float(receipt["component_l3"])
beta_m1 = float(receipt["component_f1"])
cov3m1 = float(receipt["covariance_l3_f1"])
se3 = float(events[3]["se"])
se4 = float(events[4]["se"])
se_m1 = float(events[-1]["se"])
contrast_var_from_components = se3**2 + se_m1**2 - 2 * cov3m1
beta2 = float(events[2]["b"])
beta4 = float(events[4]["b"])
mean234 = (beta2 + beta3 + beta4) / 3

result = [
    dict(estimand="beta(+3)", estimate=beta3, standard_error=se3,
         status="available_marginal_SE", formula="saved coefficient and marginal event-study SE",
         note="SE uses saved event-study table precision; covariance not needed."),
    dict(estimand="beta(+3)-beta(-1)", estimate=beta3-beta_m1,
         standard_error=float(receipt["standard_error"]),
         status="available_full_covariance", formula="sqrt(V33 + V(-1,-1) - 2 V(3,-1))",
         note=(f"Receipt covariance={cov3m1:.16g}; component reconstruction with rounded event CSV "
               f"SEs gives SE={math.sqrt(contrast_var_from_components):.12f}.")),
    dict(estimand="beta(+4)", estimate=beta4, standard_error=se4,
         status="available_marginal_SE", formula="saved coefficient and marginal event-study SE",
         note="SE uses saved event-study table precision; covariance not needed."),
    dict(estimand="mean(beta(+2), beta(+3), beta(+4))", estimate=mean234,
         standard_error=None, status="unavailable_missing_pairwise_covariances",
         formula="(beta2 + beta3 + beta4)/3",
         note="Requires all pairwise covariance terms V22,V33,V44,V23,V24,V34; only V(3,-1) is stored in the target receipt."),
]

with (OUT / "room_estimand_results.csv").open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=result[0].keys())
    w.writeheader()
    w.writerows(result)

source_files = [
    MAIN / "target_receipt.csv", MAIN / "event_study_estimates.csv",
    MAIN / "metadata.json", MAIN / "sa_rooms_first_birth_household_aligned_v1.log",
    REVIEW / "README.md", REVIEW / "event_curve_audit.csv",
    REVIEW / "all_wave_validation_receipt.json",
    ROOT / "code/model/tools/e5f_initial_housing_observer.py",
    ROOT / "code/model/tools/run_e5f_transition_calibration.py",
    ROOT / "output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/final_readout/B/selected/native_summary.json",
    OUT / "reconstruct_second_birth_shares.R",
]
hashes = {str(p.relative_to(ROOT)): sha256(p) for p in source_files if p.is_file()}
(OUT / "source_hashes.json").write_text(json.dumps(hashes, indent=2) + "\n", encoding="utf-8")

expected = float(receipt["estimate"])
assert math.isclose(beta3-beta_m1, expected, rel_tol=0, abs_tol=5e-13)
assert abs(math.sqrt(contrast_var_from_components)-float(receipt["standard_error"])) < 2e-7
print(json.dumps({"target": expected, "target_se": float(receipt["standard_error"]),
                  "contrast_se_from_rounded_components": math.sqrt(contrast_var_from_components),
                  "mean_beta_2_4": mean234, "hash_count": len(hashes)}, indent=2))
