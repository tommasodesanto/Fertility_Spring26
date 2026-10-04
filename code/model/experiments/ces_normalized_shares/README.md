# Normalized CES-limit Cobb--Douglas shares (isolated experiment)

This adapter is an experimental, runtime-only change.  It does not edit
`code/model/production/` or change finance, timing, entry, wealth, fertility
architecture, gross-estate semantics, or the canonical renter/owner
housing-service mapping through `chi`. Its separate experimental score adds
`family_rooms` as a scored moment; it leaves the native target table intact.

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

V4 and V5 smoke jobs passed the full native numerical gates and collectors, then Slurm reported failure in the EXIT-receipt bookkeeping. Production uses the separately verified launcher repair.

The local adapter suite reports 11 passed and 1 skipped; all four mock-chain
checks pass. V1 and v2 zero-solve preflights failed while packaging historical
dependencies. V3 established the three-context preflight. V4 smoke job 19132298
completed two GEs and passed its native postcheck and collector, including an
exact repeat, but Slurm reported `FAILED 1:0` afterward because the EXIT receipt
function lacked an `os` import. V5 smoke job 19132940 passed the full native
gate and collector; its 17 diagnostic plot hashes match v4. The final reviewed
launcher fixes the EXIT receipt function's import, quoting, and valid-JSON
newline only; model, search, target, and budget source stayed unchanged. Four
local and Torch EXIT fixtures passed, preserving valid JSON and the original exit code.
The four-chain production array 19133352 (`0-3%4`) was submitted October 3 at
22:31:16 EDT; all four tasks were verified `RUNNING` at 22:32:52 EDT. No
calibrated result or adoption is established. Evidence is linked in the
[experiment packet](../../../../output/model/experiments/ces_normalized_shares/overnight_v1/README.md).

V5 retains `first_birth_fixed_cost=1.9` as a free coordinate in `[0,8]`; its
original seed is unbracketed. A native GE may start only with more than 2,700
seconds left: a 900-second minimum search-GE window plus the 1,800-second
final-verification reserve. An authenticated numerical candidate that cannot
be bracketed receives a `1e12` penalty; at a bounded-budget stop, the best
candidate is postchecked. Other errors remain terminal without retries or
fallback.
