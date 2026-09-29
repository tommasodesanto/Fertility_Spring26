"""Build isolated two-birth calendar/observer sources; never edit active code.

Caller authenticates input source pins and installs returned top-level classes /
functions only in an isolated Torch process. Every textual anchor is unique;
changed or already-patched input fails rather than silently skipping work.
The flag-off numerical path is retained. Flag-on uses policy-owned probabilities.
"""
from __future__ import annotations


PREFIX = "code/model/tools/"
SOURCE_FILES = tuple(PREFIX + name for name in (
    "run_dynamic_population_transition.py",
    "run_e5f_open_population_transition.py",
    "run_e5f_transition_calibration.py",
    "e5f_initial_fertility_observer.py",
    "e5f_recent_parent_flow_observer.py",
))
SOURCE_PATHS = SOURCE_FILES


def _once(source: str, before: str, after: str, label: str) -> str:
    count = source.count(before)
    if count != 1:
        raise ValueError(f"{label}: expected one source anchor, found {count}")
    return source.replace(before, after, 1)


EXTRA_ACCESSOR = '''def policy_extra_birth_probs(
    policy: PolicyBundle, P: SimpleNamespace,
) -> np.ndarray | None:
    """Owned conditional extra opportunity; never read mutable parameter caches."""
    if not bool(getattr(P, "two_births_per_period", False)):
        return None
    probabilities = getattr(policy, "fert_extra_probs", None)
    if probabilities is None:
        raise RuntimeError("Two-birth policy lacks owned extra-birth probabilities")
    expected = np.shape(policy.fert_probs)[:-1] + (2, 2, int(P.n_child_states))
    if (np.shape(probabilities) != expected or not np.isfinite(probabilities).all()
            or np.any(probabilities < 0.0) or np.any(probabilities > 1.0)):
        raise RuntimeError("Two-birth policy has invalid extra-birth probabilities")
    return probabilities


'''


TWO_BIRTH_CALENDAR = '''    if bool(getattr(P, "two_births_per_period", False)):
        if (str(getattr(P, "child_state_mode", "")) != "independent_count"
                or model.readiness_gate_active(P)
                or model.parent_age_maturation_active(P)
                or any(bool(getattr(P, key, False)) for key in
                       ("joint_nested_choice", "fertility_nest_choice", "two_shock_choice"))):
            raise RuntimeError("Two-birth calendar requires the isolated independent-count architecture")
        if fert_extra_probs is None:
            raise RuntimeError("Two-birth calendar requires explicit policy-owned extra probabilities")
        expected = g_pre.shape[:5] + (2, 2, int(P.n_child_states))
        if (np.shape(fert_extra_probs) != expected
                or not np.isfinite(fert_extra_probs).all()
                or np.any(fert_extra_probs < 0.0) or np.any(fert_extra_probs > 1.0)):
            raise RuntimeError("Invalid two-birth calendar probability array")
        for j in range(int(P.J)):
            if not (int(P.A_f_start) <= j + 1 <= int(P.A_f_end)):
                continue
            cell_post, flow, _, _, _ = model.two_birth_fertility_step(
                g_pre[:, :, :, j], fert_probs[:, :, :, j, :, 1],
                continuation[:, :, :, j], fert_extra_probs[:, :, :, j],
                float(fecundity[j]),
            )
            out[:, :, :, j] = cell_post
            births += float(np.sum(flow))
            births_by_loc += np.sum(flow, axis=(0, 1, 3, 4))
        if not np.isfinite(out).all() or float(np.min(out)) < -1e-13:
            raise RuntimeError("Two-birth fertility produced invalid population mass")
        if abs(float(np.sum(out) - np.sum(g_pre))) > 2e-10:
            raise RuntimeError("Two-birth fertility failed mass conservation")
        return out, births, births_by_loc
'''


def patch_sources(sources: dict[str, str]) -> dict[str, str]:
    """Return five patched source texts, keyed by full repository-relative path."""
    missing = set(SOURCE_FILES) - set(sources)
    if missing:
        raise ValueError(f"Missing reporting sources: {sorted(missing)}")
    result = dict(sources)
    name = SOURCE_FILES[0]
    text = sources[name]
    text = _once(text, "    c_pol_stay: np.ndarray | None = None\n",
                 "    c_pol_stay: np.ndarray | None = None\n    fert_extra_probs: np.ndarray | None = None\n", name + " field")
    text = _once(text, "            self.fert2_probs = self.fert2_probs.copy()\n",
                 "            self.fert2_probs = self.fert2_probs.copy()\n"
                 "        if self.fert_extra_probs is not None:\n"
                 "            self.fert_extra_probs = self.fert_extra_probs.copy()\n", name + " copy")
    text = _once(text, "@dataclass\nclass PeriodEvaluation:", EXTRA_ACCESSOR + "@dataclass\nclass PeriodEvaluation:", name + " accessor")
    for obj, field in (("solution", "fert_extra_probs"), ("P", "_fert_extra_probs")):
        anchor = f'        fert2_probs=getattr({obj}, "' + ("fert2_probs" if obj == "solution" else "_fert2_probs") + '", None),\n'
        text = _once(text, anchor, anchor + f'        fert_extra_probs=getattr({obj}, "{field}", None),\n', name + " policy " + obj)
    text = _once(text,
        "    fert2_probs: np.ndarray | None = None,\n) -> tuple[np.ndarray, float, np.ndarray]:\n    out = g_pre.copy()",
        "    fert2_probs: np.ndarray | None = None,\n    *, fert_extra_probs: np.ndarray | None = None,\n) -> tuple[np.ndarray, float, np.ndarray]:\n"
        "    if bool(getattr(P, 'two_births_per_period', False)):\n"
        "        raise RuntimeError('Configure the sequential fertility operator for two-birth policies')\n"
        "    out = g_pre.copy()", name + " generic interface")
    text = _once(text, "            gated, policy.fert_probs, P, continuation\n",
        "            gated, policy.fert_probs, P, continuation,\n"
        "            fert_extra_probs=policy_extra_birth_probs(policy, P),\n", name + " evaluate")
    text = _once(text, "            g_pre, policy.fert_probs, P, policy_continuation_birth_probs(policy, P)\n",
        "            g_pre, policy.fert_probs, P, policy_continuation_birth_probs(policy, P),\n"
        "            fert_extra_probs=policy_extra_birth_probs(policy, P),\n", name + " reconstruct")
    result[name] = text

    name = SOURCE_FILES[1]
    text = sources[name]
    text = _once(text,
        '    fert2_probs: np.ndarray | None = None,\n) -> tuple[np.ndarray, float, np.ndarray]:\n    """Apply at most one sequential birth per household during the period."""',
        '    fert2_probs: np.ndarray | None = None,\n    *, fert_extra_probs: np.ndarray | None = None,\n) -> tuple[np.ndarray, float, np.ndarray]:\n'
        '    """Apply original births, or the explicit experimental two-opportunity kernel."""', name + " interface")
    text = _once(text, "    births_by_loc = np.zeros(int(P.I))\n    for j in range(int(P.J)):\n",
        "    births_by_loc = np.zeros(int(P.I))\n" + TWO_BIRTH_CALENDAR + "    for j in range(int(P.J)):\n", name + " kernel")
    result[name] = text

    name = SOURCE_FILES[2]
    text = sources[name]
    text = _once(text,
        "    flows = np.zeros((int(P.J), n_progressions), dtype=float)\n",
        "    flows = np.zeros((int(P.J), n_progressions), dtype=float)\n"
        "    two_births = bool(getattr(P, 'two_births_per_period', False))\n"
        "    extra = calendar.policy_extra_birth_probs(evaluation.policy, P)\n"
        "    expanded_risk = np.zeros_like(flows)\n", name + " flow setup")
    text = _once(text,
        "        pi_j = float(fecundity[j])\n        for zz in range(evaluation.g_pre.shape[4]):\n",
        "        pi_j = float(fecundity[j])\n"
        "        if two_births:\n"
        "            post_cell, cell_flows, _, cell_risk, _ = model.two_birth_fertility_step(\n"
        "                evaluation.g_pre[:, :, :, j], evaluation.policy.fert_probs[:, :, :, j, :, 1],\n"
        "                continuation[:, :, :, j], extra[:, :, :, j], pi_j)\n"
        "            if np.max(np.abs(post_cell - evaluation.g_post_fertility[:, :, :, j])) > 2e-10:\n"
        "                raise RuntimeError('Two-birth observer does not replay dated population')\n"
        "            flows[j] = np.sum(cell_flows, axis=(0, 1, 2, 3))\n"
        "            expanded_risk[j] = np.sum(cell_risk, axis=(0, 1, 2, 3))\n"
        "            continue\n"
        "        for zz in range(evaluation.g_pre.shape[4]):\n", name + " flow computation")
    text = _once(text, '        "birth_flow_first": flows[:, 0],\n',
        '        "birth_order_risk": expanded_risk if two_births else None,\n'
        '        "birth_flow_first": flows[:, 0],\n', name + " flow output")
    text = _once(text,
        "    policy = evaluation.policy\n    fecundity = model.get_fecundity_by_age(P)\n",
        "    policy = evaluation.policy\n"
        "    extra = calendar.policy_extra_birth_probs(policy, P)\n"
        "    fecundity = model.get_fecundity_by_age(P)\n", name + " same-policy owned cache")
    text = _once(text,
        "    joint = calendar.joint_nested_enabled(P)\n    fecundity = model.get_fecundity_by_age(P)\n",
        "    joint = calendar.joint_nested_enabled(P)\n"
        "    extra = calendar.policy_extra_birth_probs(policy, P)\n"
        "    fecundity = model.get_fecundity_by_age(P)\n", name + " dated owned cache")
    text = _once(text,
        "            birth_cohort[:, :, :, zz, 1, 1] = realized\n",
        "            birth_cohort[:, :, :, zz, 1, 1] = realized\n"
        "            if bool(getattr(P, 'two_births_per_period', False)):\n"
        "                second = realized * extra[:, :, :, j, zz, 1, 0, settled] * float(fecundity[j])\n"
        "                birth_cohort[:, :, :, zz, 1, 1] -= second\n"
        "                birth_cohort[:, :, :, zz, 2, 2] = second\n", name + " same-policy event")
    text = _once(text,
        "                treated[:, :, :, j, zz, 1, 1] = realized\n",
        "                treated[:, :, :, j, zz, 1, 1] = realized\n"
        "                if bool(getattr(P, 'two_births_per_period', False)):\n"
        "                    second = realized * extra[:, :, :, j, zz, 1, 0, settled] * float(fecundity[j])\n"
        "                    treated[:, :, :, j, zz, 1, 1] -= second\n"
        "                    treated[:, :, :, j, zz, 2, 2] = second\n", name + " dated event")
    text = _once(text,
        "            calendar.policy_continuation_birth_probs(evaluation.policy, P),\n        )\n        control_post = control_pre",
        "            calendar.policy_continuation_birth_probs(evaluation.policy, P),\n"
        "            fert_extra_probs=calendar.policy_extra_birth_probs(evaluation.policy, P),\n"
        "        )\n        control_post = control_pre", name + " destination event")
    text = _once(text,
        "    current temporary-equilibrium policy fixed and deliberately does not add a\n    second birth inside the four-year window.\n",
        "    current temporary-equilibrium policy fixed. In the isolated two-birth\n"
        "    experiment the origin first-birth families carry their actual n=1/n=2\n"
        "    destination mixture; the confirmed-childless control is unchanged.\n", name + " event doc")
    result[name] = text

    name = SOURCE_FILES[3]
    text = sources[name]
    text = _once(text,
        "    if (np.any(at_risk > pre_parity[:, 0] + FLOW_ATOL)\n",
        "    boundary_risk = pre_parity[:, :3]\n"
        "    if bool(getattr(P, 'two_births_per_period', False)):\n"
        "        boundary_risk = np.asarray(period['birth_order_risk'], dtype=float)\n"
        "        if (boundary_risk.shape != flows.shape or not np.isfinite(boundary_risk).all()\n"
        "                or np.any(boundary_risk < 0.0)):\n"
        "            raise RuntimeError('Invalid expanded within-period birth risk sets')\n"
        "    if (np.any(at_risk > pre_parity[:, 0] + FLOW_ATOL)\n", name + " risk input")
    text = _once(text, "            or np.any(flows > pre_parity[:, :3] + FLOW_ATOL)):\n",
        "            or np.any(flows > boundary_risk + FLOW_ATOL)):\n", name + " risk assertion")
    text = _once(text, '                "age_projection": "uniform_birth_time",\n',
        '                "age_projection": "uniform_birth_time",\n'
        '                "two_birth_projection_convention": (\n'
        '                    "common-event pre/post interpolation; both within-period births share the inherited event time proxy; ordered dates unavailable"\n'
        '                    if bool(getattr(P, "two_births_per_period", False)) else None),\n', name + " projection disclosure")
    text = _once(text,
        '            "birth_time_assumption": ("one parity transition per cell, uniformly distributed within the four-year age interval"\n',
        '            "birth_time_assumption": (("common-event pre/post interpolation for up to two births; both share the inherited event-time proxy; ordered within-cell dates unavailable"\n'
        '                                       if bool(getattr(P, "two_births_per_period", False)) else\n'
        '                                       "one parity transition per cell, uniformly distributed within the four-year age interval")\n', name + " full projection disclosure")
    result[name] = text

    name = SOURCE_FILES[4]
    text = sources[name]
    text = _once(text,
        '    _array(model.get_fecundity_by_age(P), "fecundity", (ages,), probability=True)\n',
        '    two_births = bool(getattr(P, "two_births_per_period", False))\n'
        '    extra = transition.calendar.policy_extra_birth_probs(policy, P)\n'
        '    _array(model.get_fecundity_by_age(P), "fecundity", (ages,), probability=True)\n', name + " owned cache")
    text = _once(text,
        "        return transition.apply_sequential_fertility(mass, first, P, continuation)\n",
        "        return transition.apply_sequential_fertility(\n"
        "            mass, first, P, continuation, fert_extra_probs=extra)\n", name + " kernel call")
    text = _once(text,
        "    # The official kernel is linear at fixed probabilities. Its original risk\n"
        "    # pools prevent a new first birth receiving a second birth in this call.\n"
        "    empty_pre = pre * empty\n"
        "    empty_post, selected_births, _ = fertility(empty_pre)\n"
        "    birth_post = empty_post * ((parity > 0) & (child_state == 1))\n"
        "    control_post = post * empty\n",
        "    # Tag original empty-dependent homes, then count successful households\n"
        "    # once even when one household receives two children in the period.\n"
        "    empty_pre = pre * empty\n"
        "    empty_post, selected_births, _ = fertility(empty_pre)\n"
        "    birth_post = empty_post * ((parity > 0) & (child_state > 0 if two_births else child_state == 1))\n"
        "    control_post = post * empty\n"
        "    selected_households = (float(empty_pre.sum()) - float((empty_post * empty).sum())\n"
        "                           if two_births else selected_births)\n", name + " household selection")
    text = _once(text,
        "        selected_birth_flow_error=abs(float(birth_post.sum()) - selected_births),\n",
        "        selected_birth_flow_error=abs(float(birth_post.sum()) - selected_households),\n", name + " household accounting")
    text = _once(text,
        "    weights = uniform_age_cell_overlap(P, 30., 56.)\n    groups = dict(\n",
        "    first_current = continuation_current = None\n"
        "    if two_births:\n"
        "        first_post, _, _ = fertility(pre * never)\n"
        "        continuation_post, _, _ = fertility(pre * former)\n"
        "        first_current = transport(first_post * ((parity > 0) & (child_state > 0)))\n"
        "        continuation_current = transport(continuation_post * ((parity > 0) & (child_state > 0)))\n"
        "        accounting['origin_tag_partition_error'] = _error(\n"
        "            first_current + continuation_current, birth_current)\n"
        "        _check(accounting['origin_tag_partition_error'], 'origin-tagged birth households')\n"
        "        accounting['empty_home_birth_households'] = float(selected_households)\n"
        "    weights = uniform_age_cell_overlap(P, 30., 56.)\n    groups = dict(\n", name + " origin tags")
    text = _once(text,
        "        first_birth=_group(birth_current, weights, parity == 1),\n"
        "        continuation_birth=_group(birth_current, weights, parity >= 2),\n",
        "        first_birth=(_group(first_current, weights) if two_births\n"
        "                     else _group(birth_current, weights, parity == 1)),\n"
        "        continuation_birth=(_group(continuation_current, weights) if two_births\n"
        "                            else _group(birth_current, weights, parity >= 2)),\n", name + " subgroup definitions")
    text = _once(text,
        '            birth_probability_source="supplied policy.fert_probs and owned policy.fert2_probs",\n',
        '            birth_probability_source=("supplied policy.fert_probs, owned policy.fert2_probs and owned policy.fert_extra_probs"\n'
        '                                      if two_births else "supplied policy.fert_probs and owned policy.fert2_probs"),\n', name + " cache metadata")
    text = _once(text,
        '                "At most one explicit birth per four-year period; top-code representative does not scale household mass",\n',
        '                ("At most two explicit births per period; selected birth households counted once and groups tagged by original parent status"\n'
        '                 if two_births else "At most one explicit birth per four-year period; top-code representative does not scale household mass"),\n', name + " timing warning")
    text = _once(text,
        "be None if that subgroup is empty. Attribution uses destination parity after\n"
        "the real one-birth kernel, which current housing transport preserves.\n",
        "be None if that subgroup is empty. Under the two-birth diagnostic, attribution\n"
        "uses original parent status because first-birth families can reach n=2.\n", name + " doc")
    result[name] = text
    return result


# Earlier integration messages used this spelling; both return the same mapping.
patch_texts = patch_sources
