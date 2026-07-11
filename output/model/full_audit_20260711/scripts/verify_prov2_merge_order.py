import numpy as np
from intergen_housing_fertility.local_panel import income_process_overrides, income_process_fingerprint
from intergen_housing_fertility.production_profile import production_profile_overrides

RHO = 0.9601845894041878
SD = 0.06453733259357768

income = income_process_overrides(5, "rouwenhorst", SD, RHO)
profile = production_profile_overrides()

collide = sorted(set(income) & set(profile))
print("colliding keys:", collide)
for k in ["z_grid", "z_weights", "income_shock_persistence"]:
    print(k, "rouwenhorst:", np.round(np.asarray(income[k], dtype=float), 6))
    print(k, "profile    :", np.round(np.asarray(profile[k], dtype=float), 6))
print("Pi_z rouwenhorst:\n", np.round(income["Pi_z"], 4))
print("Pi_z profile    :\n", np.round(profile["Pi_z"], 4))
print("fingerprint rouwenhorst:", income_process_fingerprint(income))
prof_income = {k: profile[k] for k in ["z_grid", "z_weights", "Pi_z", "income_shock_persistence"]}
print("fingerprint profile    :", income_process_fingerprint(prof_income))

# Simulate the two merge orders exactly as in run_local_panel_case
extra = dict(profile)
extra["normalize_bequest_utility"] = True
dirty = {**extra, **income}     # dirty: income last among the two
head = {**income, **extra}      # HEAD: extra_overrides last
print("dirty-order z_grid  :", np.round(np.asarray(dirty["z_grid"], dtype=float), 6))
print("head-order  z_grid  :", np.round(np.asarray(head["z_grid"], dtype=float), 6))
print("orders equal on z_grid:", np.allclose(np.asarray(dirty["z_grid"], float), np.asarray(head["z_grid"], float)))
