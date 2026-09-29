"""Exact-match, isolated deterministic sequential fertility overlay."""

SOLVER = 'code/model/intergen_eqscale_seq_optimized/solver.py'
SOURCE_PATHS = (SOLVER,)


def replace_once(source, old, new, label):
    count = source.count(old)
    if count != 1:
        raise RuntimeError(f'{label}: expected one source match, found {count}')
    return source.replace(old, new, 1)


HELPER = '''def exact_zero_birth_choice(menu):
    """Maximum value and deterministic wait-first action; dead menus choose neither."""
    value = np.max(menu, axis=3)
    probability = np.zeros_like(menu)
    live = value > DEAD_VALUE_CUTOFF
    probability[..., 0] = live & (menu[..., 0] >= menu[..., 1])
    probability[..., 1] = live & (menu[..., 1] > menu[..., 0])
    return value, probability


'''


def patch_text(source):
    if 'def exact_zero_birth_choice(' in source:
        raise RuntimeError('Already patched')
    source = replace_once(source, 'def solve_bellman_full_markov_income(\n',
                          HELPER + 'def solve_bellman_full_markov_income(\n', 'helper')
    source = replace_once(source,
        '                    lf = Vfa / P.kappa_fert\n'
        '                    ls, pr = logsumexp(lf, axis=3)\n'
        '                    pr[np.max(Vfa, axis=3) <= DEAD_VALUE_CUTOFF, :] = 0.0\n'
        '                    fert_probs[:, :, :, j, zz, :2] = pr\n'
        '                    fert_value[:, :, :, j, zz] = P.kappa_fert * ls\n',
        '                    if P.kappa_fert == 0.0:\n'
        '                        choice_value, pr = exact_zero_birth_choice(Vfa)\n'
        '                    else:\n'
        '                        lf = Vfa / P.kappa_fert\n'
        '                        ls, pr = logsumexp(lf, axis=3)\n'
        '                        pr[np.max(Vfa, axis=3) <= DEAD_VALUE_CUTOFF, :] = 0.0\n'
        '                        choice_value = P.kappa_fert * ls\n'
        '                    fert_probs[:, :, :, j, zz, :2] = pr\n'
        '                    fert_value[:, :, :, j, zz] = choice_value\n', 'first choice')
    source = replace_once(source,
        '                            l2, p2 = logsumexp(V2 / kf_cont, axis=3)\n'
        '                            p2[np.max(V2, axis=3) <= DEAD_VALUE_CUTOFF, :] = 0.0\n',
        '                            if kf_cont == 0.0:\n'
        '                                continuation_value, p2 = exact_zero_birth_choice(V2)\n'
        '                            else:\n'
        '                                l2, p2 = logsumexp(V2 / kf_cont, axis=3)\n'
        '                                p2[np.max(V2, axis=3) <= DEAD_VALUE_CUTOFF, :] = 0.0\n'
        '                                continuation_value = kf_cont * l2\n', 'later choice')
    source = replace_once(source,
        '                            V[:, :, :, j, zz, nn, cs] = kf_cont * l2\n',
        '                            V[:, :, :, j, zz, nn, cs] = continuation_value\n', 'later value')
    return source


def patch_sources(sources):
    result = dict(sources)
    result[SOLVER] = patch_text(result[SOLVER])
    return result
