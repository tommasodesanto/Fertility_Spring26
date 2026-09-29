"""Exact-match transformations for the isolated two-birth solver experiment.

This module edits strings only. The lead applies it to the authenticated source
copy on Torch; it must never be applied to the shared working solver in place.
"""

SOLVER = 'code/model/intergen_eqscale_seq_optimized/solver.py'
SOURCE_PATHS = (SOLVER,)


HELPERS = '''def two_births_active(P):
    """Validate the narrow, default-off two-opportunity architecture."""
    if not bool(getattr(P, "two_births_per_period", False)):
        return False
    if (not bool(getattr(P, "sequential_births", False))
            or bool(getattr(P, "joint_nested_choice", False))
            or not independent_child_maturation_active(P)
            or parent_age_maturation_active(P) or readiness_gate_active(P)
            or int(P.n_parity) != 4 or int(P.n_child_states) != 4):
        raise ValueError("Two births require sequential independent counts, four child states, and no joint/readiness/parent-age mode")
    return True


def two_birth_last_opportunity(wait_value, success_value, pi, kappa):
    """The sole remaining attempt after one successful birth this period.

    Its inclusive value includes an additional independent Gumbel opportunity;
    this is an explicit experimental preference/choice-tree change.
    """
    if not np.isfinite(kappa) or float(kappa) <= 0.0:
        raise ValueError("Continuation birth scale must be positive")
    if not np.isfinite(pi) or not 0.0 <= float(pi) <= 1.0:
        raise ValueError("Fecundity must lie in [0, 1]")
    menu = np.stack((wait_value, pi * success_value + (1.0 - pi) * wait_value), axis=-1)
    scaled = menu / kappa
    inclusive, _ = logsumexp(scaled, axis=menu.ndim - 1)
    shifted_weights = np.exp(scaled - np.max(scaled, axis=-1, keepdims=True))
    probabilities = shifted_weights / shifted_weights.sum(axis=-1, keepdims=True)
    probabilities[np.max(menu, axis=-1) <= DEAD_VALUE_CUTOFF, :] = 0.0
    return kappa * inclusive, probabilities


def two_birth_fertility_step(pre, first_try, continuation, extra, pi):
    """Route at most two births from original source pools.

    State arrays end in (children ever born, children at home). Continuation
    and extra arrays end in (action, original child-count slot, at-home count).
    The extra slots are original n=0/1, unlike continuation slots n-1=0/1.
    Returns post mass, births/attempts/risk by boundary, and first-birth-tagged
    post mass. The three flow arrays retain all leading state axes.
    """
    if pre.shape[-2:] != (4, 4):
        raise ValueError("Two-birth flow needs four ever-born and at-home states")
    leading = pre.shape[:-2]
    if (first_try.shape != leading or continuation.shape != leading + (2, 2, 4)
            or extra is None or extra.shape != leading + (2, 2, 4)):
        raise ValueError("Two-birth policy arrays do not match the source mass")
    if not np.isfinite(pi) or not 0.0 <= float(pi) <= 1.0:
        raise ValueError("Fecundity must lie in [0, 1]")
    post = pre.copy()
    births = np.zeros(leading + (3,))
    attempts = np.zeros_like(births)
    risk = np.zeros_like(births)
    first_post = np.zeros_like(pre)
    for nn in range(3):
        for cs in range(nn + 1):
            source = pre[..., nn, cs]
            probability = first_try if nn == 0 else continuation[..., 1, nn - 1, cs]
            attempted = source * probability
            success = pi * attempted
            extra_attempted = (success * extra[..., 1, nn, cs]
                               if nn <= 1 else np.zeros_like(success))
            extra_success = pi * extra_attempted
            post[..., nn, cs] -= success
            post[..., nn + 1, cs + 1] += success - extra_success
            births[..., nn] += success
            attempts[..., nn] += attempted
            risk[..., nn] += source
            if nn <= 1:
                post[..., nn + 2, cs + 2] += extra_success
                births[..., nn + 1] += extra_success
                attempts[..., nn + 1] += extra_attempted
                risk[..., nn + 1] += success
            if nn == 0:
                first_post[..., 1, 1] += success - extra_success
                first_post[..., 2, 2] += extra_success
    return post, births, attempts, risk, first_post


def snapshot_extra_birth_probabilities(P):
    """Own a price-specific cache; reject absent caches when the flag is on."""
    if not two_births_active(P):
        return None
    probabilities = getattr(P, "_fert_extra_probs", None)
    if probabilities is None:
        raise RuntimeError("Two-birth solution lacks its extra-opportunity probabilities")
    return probabilities.copy()


'''


def _replace(text, old, new, label):
    count = text.count(old)
    if count != 1:
        raise ValueError(f'{label}: expected one exact source match, found {count}')
    return text.replace(old, new, 1)


def patch_text(text):
    if 'def two_births_active(' in text:
        raise ValueError('Refusing to patch an already modified solver')
    text = _replace(text, 'def solve_bellman_full_markov_income(\n',
                    HELPERS + 'def solve_bellman_full_markov_income(\n', 'helpers')
    text = _replace(text,
        '    t0 = time.perf_counter()\n    natural_credit = validate_native_solvency_mode(P)\n',
        '    t0 = time.perf_counter()\n    two_births = two_births_active(P)\n    natural_credit = validate_native_solvency_mode(P)\n',
        'Bellman architecture gate')
    text = _replace(text,
        '    fert_value = np.zeros((Nb, nt, I, J, Nz))\n\n    ctx = _build_housing_stage_ctx',
        '    fert_extra_probs = (np.zeros((Nb, nt, I, J, Nz, 2, 2, ncs))\n'
        '                        if two_births else None)\n'
        '    fert_value = np.zeros((Nb, nt, I, J, Nz))\n\n    ctx = _build_housing_stage_ctx',
        'extra cache allocation')
    text = _replace(text,
        '                if bool(getattr(P, "sequential_births", False)):\n                    Vfa = np.empty((Nb, nt, I, 2))\n',
        '                if bool(getattr(P, "sequential_births", False)):\n'
        '                    kf_extra_raw = getattr(P, "kappa_fert_continuation", None)\n'
        '                    kf_extra = float(P.kappa_fert) if kf_extra_raw is None else float(kf_extra_raw)\n'
        '                    Vfa = np.empty((Nb, nt, I, 2))\n',
        'first branch continuation scale')
    text = _replace(text,
        '                        first_dest = VI[:, :, :, 1, 1]\n                    Vfa[:, :, :, 1] = (\n',
        '                        first_dest = VI[:, :, :, 1, 1]\n'
        '                    if two_births:\n'
        '                        first_dest, extra_first = two_birth_last_opportunity(\n'
        '                            first_dest, VI[:, :, :, 2, 2], pi_j, kf_extra)\n'
        '                        fert_extra_probs[:, :, :, j, zz, :, 0, 0] = extra_first\n'
        '                    Vfa[:, :, :, 1] = (\n',
        'first birth extra inclusive value')
    text = _replace(text,
        '                                cont_dest = VI[:, :, :, nn + 1, destination_cs]\n                            V2[:, :, :, 1] = (\n',
        '                                cont_dest = VI[:, :, :, nn + 1, destination_cs]\n'
        '                            if two_births and nn <= 1:\n'
        '                                cont_dest, extra_later = two_birth_last_opportunity(\n'
        '                                    cont_dest, VI[:, :, :, nn + 2, destination_cs + 1], pi_j, kf_cont)\n'
        '                                fert_extra_probs[:, :, :, j, zz, :, nn, cs] = extra_later\n'
        '                            V2[:, :, :, 1] = (\n',
        'later birth extra inclusive value')
    text = _replace(text,
        '    P._fert2_probs = fert2_probs\n    P._joint_choice = joint\n',
        '    P._fert2_probs = fert2_probs\n    P._fert_extra_probs = fert_extra_probs\n    P._joint_choice = joint\n',
        'Bellman cache export')

    # Keep the old payload ordering and append one owned cache. Flag-off arrays
    # and arithmetic are untouched; its appended value is explicitly None.
    text = _replace(text,
        '        sol._model_payload = (V, c_pol, hR_pol, bp_pol, tc, tp, lp_j, fp, fv, r, p, P._fert2_probs.copy())\n',
        '        sol._model_payload = (V, c_pol, hR_pol, bp_pol, tc, tp, lp_j, fp, fv, r, p, P._fert2_probs.copy(), snapshot_extra_birth_probabilities(P))\n'
        '        sol.fert_extra_probs = sol._model_payload[-1]\n',
        'fast price cache ownership')
    text = _replace(text,
        '    V, c_pol, hR_pol, bp_pol, tc, tp, lp_j, fp, fv, r, p, fert2_probs = payload\n',
        '    V, c_pol, hR_pol, bp_pol, tc, tp, lp_j, fp, fv, r, p, fert2_probs, fert_extra_probs = payload\n'
        '    P._fert_extra_probs = None if fert_extra_probs is None else fert_extra_probs.copy()\n',
        'accepted cache restoration')
    text = _replace(text,
        '        fert2_probs=getattr(P, "_fert2_probs", None),\n        joint_choice=getattr(P, "_joint_choice", None),\n',
        '        fert2_probs=getattr(P, "_fert2_probs", None),\n'
        '        fert_extra_probs=snapshot_extra_birth_probabilities(P),\n'
        '        joint_choice=getattr(P, "_joint_choice", None),\n',
        'full solution ownership')

    # Restrict edits to this function so repeated legacy patterns cannot match.
    start = text.index('def forward_distribution_markov_income(\n')
    end = text.index('\ndef collapse_markov_policy(', start)
    forward = text[start:end]
    forward = _replace(forward,
        '    fec = get_fecundity_by_age(P)\n',
        '    two_births = two_births_active(P)\n'
        '    extra_policy = getattr(P, "_fert_extra_probs", None)\n'
        '    if two_births and extra_policy is None:\n'
        '        raise RuntimeError("Two-birth KFE lacks its matching extra probabilities")\n'
        '    fec = get_fecundity_by_age(P)\n',
        'KFE architecture and cache gate')
    begin = forward.index('                    # Snapshot every upward at-risk pool BEFORE any birth flow\n')
    finish = forward.index('                    birth_mass = realized1\n', begin)
    original = forward[begin:finish]
    flagged = '''                    first_birth_post = None
                    if two_births:
                        source = g[:, :, :, j, zz, :, :].copy()
                        post, born, attempted, risk, first_birth_post = two_birth_fertility_step(
                            source, pa[:, :, :, 1], P._fert2_probs[:, :, :, j, zz],
                            extra_policy[:, :, :, j, zz], pi_j)
                        g[:, :, :, j, zz, :, :] = post
                        realized1 = first_birth_post.sum(axis=(-2, -1))
                        born_total = born.sum(axis=(0, 1, 2))
                        attempted_total = attempted.sum(axis=(0, 1, 2))
                        risk_total = risk.sum(axis=(0, 1, 2))
                        first_births_by_age[j] += float(born_total[0])
                        second_births_by_age[j] += float(born_total[1])
                        third_births_by_age[j] += float(born_total[2])
                        second_attempts_by_age[j] += float(attempted_total[1])
                        third_attempts_by_age[j] += float(attempted_total[2])
                        second_at_risk_by_age[j] += float(risk_total[1])
                        third_at_risk_by_age[j] += float(risk_total[2])
                        total_births += float(born_total.sum())
                        births_by_loc += born.sum(axis=(0, 1, 3))
                    else:
'''
    flagged += ''.join('    ' + line if line.strip() else line for line in original.splitlines(keepends=True))
    forward = forward[:begin] + flagged + forward[finish:]
    forward = _replace(forward,
        '                        birth_cohort[:, :, :, zz, 1, 1] = realized1\n',
        '                        if two_births:\n'
        '                            birth_cohort[:, :, :, zz, :, :] = first_birth_post\n'
        '                        else:\n'
        '                            birth_cohort[:, :, :, zz, 1, 1] = realized1\n',
        'first-birth event cohort destinations')
    text = text[:start] + forward + text[end:]
    return text


def patch_sources(sources):
    """Return transformed source strings; preserve all other mapping entries."""
    result = dict(sources)
    result[SOLVER] = patch_text(result[SOLVER])
    return result
