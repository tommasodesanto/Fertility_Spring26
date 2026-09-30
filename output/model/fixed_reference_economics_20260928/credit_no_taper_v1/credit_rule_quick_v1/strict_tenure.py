"""Native tenure kernels with only the strict owner-sale gate added."""
import numpy as np
from numba import njit
from intergen_eqscale_seq_optimized.kernels import _interp_with_clip

@njit(cache=True)
def tenure_choice_kernel(
    Vd,                  # (Nb, nt, I, npar, ncs)
    b_grid,              # (Nb,)
    heq,                 # (I, nt)
    hcost,               # (I, nt)
    dp_arr,              # (I, nt, npar, ncs)
    bmo,                 # (I, nt, npar, ncs)
    birth_dp,            # (npar, ncs, nt, nt) bool
    birth_entry_grant,   # (I, nt, npar, ncs)
    Vd_stay,             # (Nb, nt, I, npar, ncs) values read at to == tn
    strict_interpolated_support=False,
    transaction_support=False,
):
    # Discrete tenure-choice argmax over `tn` given conditional values Vd
    # for each (origin tenure `to`, location, parity, child-state, b).
    # Three branches handle: stay/move-as-renter (tn=0), buy-on-entry
    # (to=0 -> tn>=1, with optional birth grant or entry grant), and
    # sell-then-rebuy (to>=1 -> tn != to). Infeasible states (below
    # down-payment threshold dp or below borrowing limit bm) get -1e10.
    # A stayer (to == tn) reads Vd_stay, which carries the stayer mortgage
    # floors when that switch is on and matches Vd bit for bit otherwise.
    Nb, nt, I, npar, ncs = Vd.shape
    VH = np.empty((Nb, nt, I, npar, ncs))
    tcj = np.empty((Nb, nt, I, npar, ncs), dtype=np.int16)
    NEG_INF = -1e10
    for id_ in range(I):
        for to in range(nt):
            sp = heq[id_, to] if to > 0 else 0.0
            for nn in range(npar):
                for cs in range(ncs):
                    dpn0_unused = 0.0  # placeholder; per-tn dp/bm pulled inside tn loop
                    for b in range(Nb):
                        bg_b = b_grid[b]
                        best_v = NEG_INF
                        best_tn = 0
                        # tn = 0 (renter)
                        if to == 0:
                            v0 = Vd[b, 0, id_, nn, cs]
                        else:
                            ba = bg_b + sp
                            v0 = _interp_with_clip(b_grid, Vd[:, 0, id_, nn, cs], ba, strict_interpolated_support, transaction_support)
                            if bg_b + heq[id_, to] < 0.0:
                                v0 = NEG_INF
                        if v0 > best_v:
                            best_v = v0
                            best_tn = 0
                        # tn >= 1 (owner tenures)
                        for tn in range(1, nt):
                            hc = hcost[id_, tn]
                            dpn = dp_arr[id_, tn, nn, cs]
                            bmn = bmo[id_, tn, nn, cs]
                            if to == tn:
                                v_tn = Vd_stay[b, tn, id_, nn, cs]
                            elif to == 0:
                                bab = bg_b - hc
                                if birth_dp[nn, cs, to, tn]:
                                    bag = bab if bab > bmn else bmn
                                    v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], bag, strict_interpolated_support, transaction_support)
                                elif birth_entry_grant[id_, tn, nn, cs] > 0:
                                    gfix = birth_entry_grant[id_, tn, nn, cs]
                                    babg = bab + gfix
                                    v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], babg, strict_interpolated_support, transaction_support)
                                    if (bg_b + gfix) < dpn or babg < bmn:
                                        v_tn = NEG_INF
                                else:
                                    v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], bab, strict_interpolated_support, transaction_support)
                                    if bg_b < dpn or bab < bmn:
                                        v_tn = NEG_INF
                            else:
                                bar = bg_b + sp - hc
                                v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], bar, strict_interpolated_support, transaction_support)
                                dpc = dpn - sp
                                if bg_b < dpc or bar < bmn:
                                    v_tn = NEG_INF
                            if v_tn > best_v:
                                best_v = v_tn
                                best_tn = tn
                        VH[b, to, id_, nn, cs] = best_v
                        tcj[b, to, id_, nn, cs] = best_tn
    return VH, tcj

@njit(cache=True)
def tenure_logit_kernel(
    Vd,                  # (Nb, nt, I, npar, ncs)
    b_grid,              # (Nb,)
    heq,                 # (I, nt)
    hcost,               # (I, nt)
    dp_arr,              # (I, nt, npar, ncs)
    bmo,                 # (I, nt, npar, ncs)
    birth_dp,            # (npar, ncs, nt, nt) bool
    birth_entry_grant,   # (I, nt, npar, ncs)
    kappa,               # taste-shock scale
    Vd_stay,             # (Nb, nt, I, npar, ncs) values read at to == tn
    transaction_support=False,
):
    Nb, nt, I, npar, ncs = Vd.shape
    VH = np.empty((Nb, nt, I, npar, ncs))
    tcj = np.empty((Nb, nt, I, npar, ncs), dtype=np.int16)
    probs = np.zeros((Nb, nt, I, npar, ncs, nt), dtype=np.float32)
    NEG_INF = -1e10
    tiny_kappa = kappa if kappa > 1e-12 else 1e-12
    vals = np.empty(nt)
    for id_ in range(I):
        for to in range(nt):
            sp = heq[id_, to] if to > 0 else 0.0
            for nn in range(npar):
                for cs in range(ncs):
                    for b in range(Nb):
                        bg_b = b_grid[b]
                        best_v = NEG_INF
                        best_tn = 0
                        if to == 0:
                            v0 = Vd[b, 0, id_, nn, cs]
                        else:
                            ba = bg_b + sp
                            v0 = _interp_with_clip(b_grid, Vd[:, 0, id_, nn, cs], ba, False, transaction_support)
                            if bg_b + heq[id_, to] < 0.0:
                                v0 = NEG_INF
                        vals[0] = v0
                        if v0 > best_v:
                            best_v = v0
                            best_tn = 0
                        for tn in range(1, nt):
                            hc = hcost[id_, tn]
                            dpn = dp_arr[id_, tn, nn, cs]
                            bmn = bmo[id_, tn, nn, cs]
                            if to == tn:
                                v_tn = Vd_stay[b, tn, id_, nn, cs]
                            elif to == 0:
                                bab = bg_b - hc
                                if birth_dp[nn, cs, to, tn]:
                                    bag = bab if bab > bmn else bmn
                                    v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], bag, False, transaction_support)
                                elif birth_entry_grant[id_, tn, nn, cs] > 0:
                                    gfix = birth_entry_grant[id_, tn, nn, cs]
                                    babg = bab + gfix
                                    v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], babg, False, transaction_support)
                                    if (bg_b + gfix) < dpn or babg < bmn:
                                        v_tn = NEG_INF
                                else:
                                    v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], bab, False, transaction_support)
                                    if bg_b < dpn or bab < bmn:
                                        v_tn = NEG_INF
                            else:
                                bar = bg_b + sp - hc
                                v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], bar, False, transaction_support)
                                dpc = dpn - sp
                                if bg_b < dpc or bar < bmn:
                                    v_tn = NEG_INF
                            vals[tn] = v_tn
                            if v_tn > best_v:
                                best_v = v_tn
                                best_tn = tn
                        tcj[b, to, id_, nn, cs] = best_tn
                        if best_v <= -1e9:
                            VH[b, to, id_, nn, cs] = best_v
                            for tn in range(nt):
                                probs[b, to, id_, nn, cs, tn] = 0.0
                        else:
                            denom = 0.0
                            for tn in range(nt):
                                ev = np.exp((vals[tn] - best_v) / tiny_kappa)
                                probs[b, to, id_, nn, cs, tn] = ev
                                denom += ev
                            if denom <= 0.0:
                                VH[b, to, id_, nn, cs] = best_v
                                probs[b, to, id_, nn, cs, best_tn] = 1.0
                            else:
                                for tn in range(nt):
                                    probs[b, to, id_, nn, cs, tn] = probs[b, to, id_, nn, cs, tn] / denom
                                VH[b, to, id_, nn, cs] = best_v + tiny_kappa * np.log(denom)
    return VH, tcj, probs
