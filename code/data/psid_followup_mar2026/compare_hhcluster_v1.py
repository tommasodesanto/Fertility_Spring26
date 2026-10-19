"""Compare household- and person-clustered SEs for the A2h first-birth event study.

Reads output/sa_first_birth_hhcluster_v1/coefficients_household.csv and the
person-clustered runs with identical point estimates; writes comparison.csv.
"""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent / "output"
REPO = Path(__file__).resolve().parents[3]
PERSON = {
    "rooms_full": HERE / "sa_rooms_first_birth_v2/A2h/coefficients.csv",
    "own": HERE / "sa_first_birth_outcomes_v3/A2h_own/coefficients.csv",
    "moved": HERE / "sa_first_birth_outcomes_v3/A2h_moved/coefficients.csv",
    "rooms_cs": REPO / "output/model/first_birth_rooms_event_20261010/measurement_bridge/original/rooms_coefficients.csv",
}
WINDOWS = ["Wleft", "Wm2", "Wm1", "Wp1", "Wp2", "Wp3", "Wp4", "Wp5", "Wp6", "Wright"]

hh = {}
for r in csv.DictReader(open(HERE / "sa_first_birth_hhcluster_v1/coefficients_household.csv")):
    hh[(r["run"], r["window"])] = (float(r["estimate"]), float(r["se_household"]))

rows, maxdiff = [], 0.0
for run, path in PERSON.items():
    pers = {r["coefficient"]: (float(r["estimate"]), float(r["standard_error"])) for r in csv.DictReader(open(path))}
    for w in WINDOWS:
        b_h, s_h = hh[(run, w)]
        b_p, s_p = pers[w]
        maxdiff = max(maxdiff, abs(b_h - b_p))
        rows.append({"run": run, "window": w, "estimate": b_h, "estimate_person_run": b_p,
                     "se_person": s_p, "se_household": s_h, "ratio": s_h / s_p})
with open(HERE / "sa_first_birth_hhcluster_v1/comparison.csv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0]))
    wr.writeheader()
    wr.writerows(rows)
print(f"max |estimate difference| = {maxdiff:.3e}")
for r in rows:
    print(f"{r['run']:11s} {r['window']:7s} {r['estimate']:8.4f} {r['se_person']:.4f} {r['se_household']:.4f} {r['ratio']:.3f}")
