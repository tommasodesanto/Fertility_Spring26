# Normalized CES-limit Cobb--Douglas shares (isolated experiment)

This adapter is an experimental, runtime-only change.  It does not edit
`code/model/production/`, does not change targets, finance, timing, entry,
wealth, fertility architecture, gross-estate semantics, or the canonical
renter/owner housing-service mapping through `chi`.

For every family state, including childless states, the normalized material
composite denominator applies; it evaluates the composite as

\[
Q(c,s;a)=\frac{c^a s^{1-a}}{a^a(1-a)^{1-a}},\qquad
M(m)=e(m)^{\sigma-1}[a(m)^{a(m)}(1-a(m))^{1-a(m)}]^{\sigma-1}.
\]

The raw utility-unit switch is intentional: the normalized composite makes
the equal-price optimum equal total material expenditure, but the CRRA value
multiplier is `M(m)` above.  There is no `alpha0` numerator and no reference
rent (`rstar`). The author-authorized rule is `alpha(m)=.733` for `m=0` and
`alpha(m)=clip(.733-delta_alpha_jump-delta_alpha*m,.05,.95)` for parents.
Both share parameters are free on `[0,.25]`. Child benefits retain their
supplied curvature and the physical housing floor is fixed at `h_P=0`.

`adapter.make_evaluator` has the canonical evaluator call shape.  It applies
the hook only while a solver call runs, restoring the two active imported
lookups: `production.engine.shared.apply_child_preferences` and
`production.equilibrium.build_context`.  `preflight_contexts(out)` is an
actual-context check only in the isolated staged dependency snapshot
(`CES_NORMALIZED_SHARES_STAGED_CONTEXT=1`); otherwise it writes one explicit
deferred-dependency receipt and performs no context build. The existing
`delta_alpha_jump` and `delta_alpha` rows record the two free share parameters;
`h_P` remains a fixed-zero row documenting removal of the floor. The 31-row
manifest stays unique; `alpha_cons` remains explicit. The legacy reference-rent input is labeled inactive.
The experiment identity and multiplier are
recorded in the context; they do not merely alter the returned
optimizer receipt.

Each reporting-context directory and each completed selected-root/repeat
report directory keeps native `target_fit.csv` unchanged, then independently
writes `target_fit_experimental.csv` and `utility_contract.json`. The new
contract, `ces_normalized_jump_slope_family_rooms_v1`, promotes only
`family_rooms` to scored, with target `0.38509964969278165` and weight
`280.52808370152104` (11 scored rows of 14 for 11 free coordinates). Its
ACS target is `E[min(ROOMS,9)|NCHILD>=3,YNGCH<18]-E[min(ROOMS,9)|NCHILD=1-2,YNGCH<18]`,
using 2005/06 national heads aged 30--55, positive HHWT, no DUE restriction,
and a weighted mean difference without fixed effects. National uncertainty is
unavailable; the weight is inherited from a 42-metro bootstrap, not national
uncertainty. The native model-dependent-child observer remains a proxy. The
full target and weight fingerprints are pinned; 11 moments for 11 free
coordinates do not establish informative rank.

## Verification status

The local adapter suite reports 8 passed and 1 skipped; all four mock-chain
checks pass. The v1 and v2 zero-solve preflights failed while packaging
historical dependencies. Attempt 3's v3 package is being prepared, not launched,
so no actual mounted reporting preflight, numerical smoke, calibration, or
adoption is established. Once launched, the native selected-point postcheck
must verify all 11 coordinates, the complete 14-row experimental target CSV,
the 31-row parameter table, 17 diagnostic plots, and an exact repeat.
