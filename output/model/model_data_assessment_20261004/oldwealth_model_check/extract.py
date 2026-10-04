"""Read-only old-wealth audit of the adopted cached stationary solution."""
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "code/model"))
from production.storage import load_case
from production.engine.shared import annual_gross_income_at_state
from production.engine.utils import weighted_quantile

OUT = Path(__file__).resolve().parent
CASE = ROOT / "output/model/local_solution/cases/20261003T175652812716Z_b1c72f13"
PSID = ROOT / "code/data/psid_followup_mar2026/output/model_assessment/psid_2005_2007.csv"
PSID_META = PSID.with_name("metadata.json")
ASSESSMENT_META = CASE / "aggregate_plots/model_data_assessment/metadata.json"
def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()

result, case = load_case(CASE)
s, P = result.solution, result.P
g = np.asarray(s.g_beginning_distribution, float)
bg = np.asarray(result.b_grid, float)
ph = np.atleast_1d(np.asarray(s.owner_asset_price, float))
z = np.atleast_1d(np.asarray(P.z_grid, float))
assert g.ndim == 7 and g.shape[0] == bg.size
assert g.shape[4] == z.size

rows = []
all_vals = []
all_wts = []
observer_vals = []
observer_wts = []
observer_wealth_vals = []
for j in range(int(P.J)):
    age = float(P.age_start + j * P.da)
    mass = float(g[:, :, :, j].sum())
    fm = hm = own = debt_own = neg_financial = neg_total = 0.0
    vals, wts = [], []
    for i in range(int(P.I)):
        for ten in range(g.shape[1]):
            house = ph[i] * float(P.H_own[ten - 1]) if ten else 0.0
            wealth = bg + house
            for zz in range(g.shape[4]):
                a = g[:, ten, i, j, zz].sum(axis=(1, 2))
                fm += float(np.dot(a, bg))
                hm += float(a.sum()) * house
                own += float(a.sum()) * int(ten > 0)
                debt_own += float(a[bg < 0].sum()) * int(ten > 0)
                neg_financial += float(a[bg < 0].sum())
                neg_total += float(a[wealth < 0].sum())
                overlap = max(0.0, min(age + float(P.da), 85.0) - max(age, 76.0)) / float(P.da)
                if overlap > 0:
                    income = annual_gross_income_at_state(P, i, j, float(z[zz]))
                    observer_vals.append(wealth / income)
                    observer_wts.append(overlap * a)
                    observer_wealth_vals.append(wealth)
                    if 76 <= age <= 84:
                        vals.append(wealth / income)
                        wts.append(a)
    rows.append(dict(age=int(age), mass=mass, financial_mean=fm/mass,
                     home_value_mean=hm/mass, total_wealth_mean=(fm+hm)/mass,
                     owner_share=own/mass, owner_debt_share=debt_own/own if own else None,
                     negative_financial_share=neg_financial/mass,
                     negative_total_share=neg_total/mass))
    if vals:
        all_vals += vals
        all_wts += wts

v = np.concatenate(all_vals)
w = np.concatenate(all_wts)
q50 = float(weighted_quantile(v, w, .5))
q90 = float(weighted_quantile(v, w, .9))
ov = np.concatenate(observer_vals)
ow = np.concatenate(observer_wts)
op50 = float(weighted_quantile(ov, ow, .5))
op90 = float(weighted_quantile(ov, ow, .9))
wealth_values = np.concatenate(observer_wealth_vals)
model_wealth_p50 = float(weighted_quantile(wealth_values, ow, .5))
model_wealth_p90 = float(weighted_quantile(wealth_values, ow, .9))
with PSID.open(newline="") as stream:
    psid_rows = list(csv.DictReader(stream))
def number(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return float("nan")
data_ages = np.array([number(row["age"]) for row in psid_rows])
data_wealth = np.array([number(row["total_net_wealth"]) for row in psid_rows])
data_weights = np.array([number(row["weight"]) for row in psid_rows])
keep = ((data_ages >= 76) & (data_ages <= 84) & np.isfinite(data_wealth)
        & np.isfinite(data_weights) & (data_weights > 0))
data_p50, data_p90 = (float(x) for x in weighted_quantile(data_wealth[keep], data_weights[keep], [.5, .9]))
normalizer = float(json.loads(ASSESSMENT_META.read_text())["normalizers"]["psid_annual_gross_labor_age18_65"])
summary = {
    "case": str(case.relative_to(ROOT)),
    "state_shape": list(g.shape),
    "old_model_age_cells": [r["age"] for r in rows if 76 <= r["age"] <= 84],
    "old_ratio_p50": q50,
    "old_ratio_p90": q90,
    "old_ratio_p90_p50": q90/q50,
    "saved_old_ratio_p90_p50": float(s.old_total_wealth_to_annual_income_p90_p50_7684),
    "observer_age_cell_weights": {"74": 0.5, "78": 1.0, "82": 0.75},
    "observer_old_ratio_p50": op50,
    "observer_old_ratio_p90": op90,
    "observer_old_ratio_p90_p50": op90/op50,
    "wealth_only_diagnostic_not_target": {
        "definition": "PSID finite NETWORTHR ages 76–84; model beginning b+pH with survey age-cell overlap; no family-income or child-history filter",
        "psid_n": int(keep.sum()),
        "psid_weighted_p50_2022usd": data_p50,
        "psid_weighted_p90_2022usd": data_p90,
        "psid_p90_p50": data_p90/data_p50,
        "psid_working_gross_earnings_normalizer_2022usd": normalizer,
        "psid_p50_normalized": data_p50/normalizer,
        "psid_p90_normalized": data_p90/normalizer,
        "model_p50_native_units": model_wealth_p50,
        "model_p90_native_units": model_wealth_p90,
        "model_p90_p50": model_wealth_p90/model_wealth_p50,
    },
    "source_sha256": {
        "case_native_result_npz": sha256(CASE / "native_result.npz"),
        "psid_selected_cache": sha256(PSID),
        "psid_cache_metadata": sha256(PSID_META),
        "assessment_metadata": sha256(ASSESSMENT_META),
    },
    "source_paths": {
        "case": str(CASE.relative_to(ROOT)),
        "psid_cache": str(PSID.relative_to(ROOT)),
        "psid_cache_metadata": str(PSID_META.relative_to(ROOT)),
        "assessment_metadata": str(ASSESSMENT_META.relative_to(ROOT)),
    },
    "realized_owner_share_last_cell": float(np.asarray(s.own_by_age)[-1]),
    "pension": float(P.pension),
    "bequest_spec": P.bequest_spec,
    "normalize_bequest_utility": bool(P.normalize_bequest_utility),
}
with (OUT / "age_profile.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
(OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
for row in rows:
    if row["age"] >= 66:
        print(row)
