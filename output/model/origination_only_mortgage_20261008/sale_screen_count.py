"""Sale-to-rent screen R b + (1-psi) p H >= 0 (kernels.py:341/419, active because unsecured_credit_limit = 0.0 is set).
Count childless households aged 22-33 (start of period, origin tenure owned) that the screen bars from selling into renting
although the sale-to-rent budget is feasible with current income: R b + S + y - r h_bar(branch) - c_min > 0 at b' = 0.
Waiting branch: childless renter (h_bar = 0). Birth branch: one child at home (h_bar = h_P = 2.596 rooms at rent r per room).
Weights: start-of-period distribution (birth_count_pre_distribution). Also weighted by the first-birth attempt probability p
(birth branch) and 1 - p (waiting branch)."""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE)); import run_cells as rc
P, grid, price, T, _ = rc.build('R80')
from model.engine.shared import income_at_state, income_transition_values
z = np.asarray(income_transition_values(P)[0], float)
R, psi, H, cmin = float(P.R_gross), float(P.psi), np.asarray(P.H_own, float), float(P.c_min)
r = float(P.user_cost_rate) * price; hP = float(P.hbar_first_child_jump); b = np.asarray(grid, float).reshape(-1)
js = [1, 2, 3]
Y = np.array([[income_at_state(P, 0, j, float(zz)) for zz in z] for j in js])     # (3, Nz)
rows = []
for c in rc.CELLS:
    S = np.load(HERE / 'cells' / c / 'solution_arrays.npz')
    pre = S['birth_count_pre_distribution'][:, 1:, 0, 1:4, :, 0, 0].astype(float)      # (b, owned rung, j, z)
    p = S['birth_count_action_probs'][:, 1:, 0, 1:4, :, 0, 0, 1].astype(float)
    Sv = (1 - psi) * price * H[None, :, None, None]; B = b[:, None, None, None]; y = Y[None, None, :, :]
    fail = R * B + Sv < 0
    feas_wait = R * B + Sv + y - cmin > 0
    feas_birth = R * B + Sv + y - r * hP - cmin > 0
    ltv = -B / (price * H[None, :, None, None])
    tot_owners = pre.sum(); tot_childless = S['birth_count_pre_distribution'][:, :, 0, 1:4, :, 0, 0].sum()
    m = lambda sel, w=pre: float((w * sel).sum())
    rows.append(dict(cell=c, owners_share_of_childless=tot_owners / tot_childless, screen_fails=m(fail) / tot_owners,
                     barred_feasible_wait=m(fail & feas_wait) / tot_owners, barred_feasible_birth=m(fail & feas_birth) / tot_owners,
                     barred_feasible_wait_x_1mp=m(fail & feas_wait, pre * (1 - p)) / tot_owners, barred_feasible_birth_x_p=m(fail & feas_birth, pre * p) / tot_owners,
                     above_8684=m(ltv > (1 - psi) / R) / tot_owners, max_ltv=float(np.broadcast_to(ltv, pre.shape)[pre > 1e-12].max()),
                     barred_mass_of_population=m(fail & feas_wait)))
hdr = ['cell', 'owners / childless 22-33', 'screen fails', 'barred, feasible (wait)', 'barred, feasible (birth)', 'x (1-p)', 'x p', 'LTV > 86.84%', 'max entering LTV']
md = '| ' + ' | '.join(hdr) + ' |\n|' + '---|' * len(hdr) + '\n' + ''.join(
    f"| {d['cell']} | {100*d['owners_share_of_childless']:.1f}% | {100*d['screen_fails']:.2f}% | {100*d['barred_feasible_wait']:.2f}% | {100*d['barred_feasible_birth']:.2f}% | "
    f"{100*d['barred_feasible_wait_x_1mp']:.2f}% | {100*d['barred_feasible_birth_x_p']:.2f}% | {100*d['above_8684']:.2f}% | {100*d['max_ltv']:.1f}% |\n" for d in rows)
note = (f'Shares are of childless owners aged 22-33 (start of period). r = {r:.4f} per room-period, h_P = {hP:.3f}, c_min = {cmin}, '
        f'cutoff (1-psi)/R = {(1-psi)/R:.4f}. "x p" / "x (1-p)" weight each state by its first-birth attempt / waiting probability.\n')
(HERE / 'sale_screen_count.md').write_text(note + '\n' + md); print(note); print(md)
json.dump(rows, open(HERE / 'sale_screen_count.json', 'w'), indent=1)
