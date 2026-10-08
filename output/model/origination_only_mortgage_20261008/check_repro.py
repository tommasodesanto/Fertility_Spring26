"""R80 (switch off) vs the saved 14.402 baseline packet and vs rental-menu cell A6; R95 vs cell A6_LTV95."""
import numpy as np
from pathlib import Path
M = Path('../rental_menu_precaution_20261007/cells')
REF = Path('../jmp_draft_deck_20261005/refit_best_14p40_mac_r3_chain0/cases_v2/baseline/solution_arrays.npz')
KEYS = ['V', 'g', 'birth_count_action_probs', 'birth_count_pre_distribution', 'hR_pol', 'bp_pol', 'c_pol', 'tenure_probs']
def cmp(a, r, label):
    a, r = np.load(a), np.load(r); worst = 0.0; same = True
    for k in KEYS:
        if k in a.files and k in r.files:
            x, y = a[k].astype(float), r[k].astype(float)
            same &= bool(np.array_equal(x, y, equal_nan=True)); worst = max(worst, float(np.nanmax(np.abs(x - y))) if x.size else 0.0)
    print(f'{label}: bitwise identical on {len([k for k in KEYS if k in a.files and k in r.files])} arrays = {same}, max abs diff {worst:.3g}')
cmp('cells/R80/solution_arrays.npz', REF, 'R80 vs saved 14.402 packet')
cmp('cells/R80/solution_arrays.npz', M / 'A6/solution_arrays.npz', 'R80 vs rental-menu A6')
cmp('cells/R95/solution_arrays.npz', M / 'A6_LTV95/solution_arrays.npz', 'R95 vs rental-menu A6_LTV95')
