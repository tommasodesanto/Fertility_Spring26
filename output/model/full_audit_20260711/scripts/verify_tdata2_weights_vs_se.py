"""Verifier check for finding T-DATA-2: used weights vs inverse-variance from moment_se.csv."""
import csv

from intergen_housing_fertility import calibration as cal

targets, weights = cal.get_target_set("candidate_replacement_post_audit_v1")
weights = dict(weights)
weights["aggregate_mean_occupied_rooms_18_85"] = 6.0  # overnight additional weight

se = {}
path = "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/moment_standard_errors/output/moment_se.csv"
with open(path) as f:
    for row in csv.DictReader(f):
        se[row["key"]] = row["boot_se"]

print(f"{'moment':55s} {'w_used':>8s} {'boot_se':>10s} {'1/se^2':>12s} {'used/inv':>10s}")
ratios = {}
for k, w in sorted(weights.items()):
    s = se.get(k, "ABSENT")
    if s in ("NA", "ABSENT"):
        print(f"{k:55s} {w:8.1f} {s:>10s} {'--':>12s} {'--':>10s}")
        continue
    inv = 1.0 / float(s) ** 2
    ratios[k] = w / inv
    print(f"{k:55s} {w:8.1f} {float(s):10.5f} {inv:12.1f} {w / inv:10.5f}")

mx = max(ratios, key=ratios.get)
mn = min(ratios, key=ratios.get)
print(f"\nmax ratio {mx} = {ratios[mx]:.4f}; min ratio {mn} = {ratios[mn]:.6f}")
print(f"spread max/min = {ratios[mx] / ratios[mn]:.1f}x")
print(f"n moments with SE: {len(ratios)}; total weighted moments: {len(weights)}")
