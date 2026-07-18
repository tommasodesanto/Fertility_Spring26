# HANDOFF: the income process — situation briefing (no plan, no recommendation)

Audience: a fresh Fable session. Your first job is ONLY to familiarize:
read this file, then the artifacts it points to, and confirm your
understanding back to Tommaso. Do not propose a plan, do not rank options,
do not start work. Tommaso will decide the direction together with you.
General session context and operational gotchas:
`docs/prompts/HANDOFF_fable_next_session_20260718.md`.

## The situation

The model's household income process is an annual AR(1) in logs
(persistence 0.960, innovation s.d. 0.0645), sampled at the 4-year model
frequency, discretized with a 5-state Rouwenhorst chain. Two facts about it:

1. The persistence is inside the range PSID studies report (0.88–0.96).
   The innovation s.d. is roughly one third of the range those studies
   report (0.12–0.25). Historically the value was reverse-engineered from
   an older ad hoc grid; it survives because larger values break the model
   (next point). The paper currently states this honestly: the variance is
   "provisional," with a footnote explaining the feasibility conflict, and
   the calibration section was circulated to advisors in that form on
   2026-07-17.

2. Why larger risk breaks the model: households must afford a subsistence
   bundle every period (Stone–Geary consumption floor c_bar_0 plus a
   minimum rental dwelling; together ~32% of mean annual income at the
   current calibrated c_bar_0 ≈ 31%). The model has no safety net, no
   default, no family transfers, no labor-supply margin. With innovation
   s.d. at 0.20, the lowest persistent income state pays ~20% of mean
   income at young ages — below the bundle — and a household stuck there
   has an empty budget set: no feasible choice exists. The solver's
   feasibility gate (installed 2026-07-11 after a real bug) then rejects
   the whole parameter vector. This was verified by direct probes,
   including with all entry debt removed. A related fact: the current
   calibration sits exactly ON the feasibility frontier — shaving the
   bottom income state by ~2% (e.g., through grid renormalization) already
   breaks it.

There is a second, related weakness at the other end of the distribution:
the model has no late-life wealth tail (estate p90/p50 = 1.75 vs 3.45 in
the PSID; untargeted diagnostic). More income risk was one hoped-for fix
for that too.

## What has been tried, and what the evidence showed

These are computational/mathematical results, not opinions. Artifacts:
probes and censuses summarized in
`~/Desktop/income_risk_advice_20260717/2_PROBLEM_STATEMENT.md` (items 1–8)
and `10_topstate_probe_results.txt`; the full verified literature survey and
an outside (ChatGPT-Pro) assessment in
`docs/model/intergen_income_risk_feasibility_decision_memo_20260717.md`.

- Giving entrants more initial wealth does not fix it: the bottom state is
  a flow shortfall (income < required spending every period), so any finite
  buffer only delays the failure; even zero-debt entry breaks at age 22.
- Shifting the mean of the innovation does nothing: a normal's support is
  unbounded, and after mean-one normalization the shift is a no-op.
- Truncating the income support from below and then raising sigma: the
  bound level is either the measured safety net (~0.15 of mean income,
  which is below the bundle, so infeasibility remains) or a level chosen to
  equal the bundle (circular); and with the bottom pinned, the added
  variance is all upside, so measured downside risk does not actually
  increase.
- A rare Castaneda-style high-income state (<1% of households) moves the
  estate p90/50 barely (+0.06 of the +1.70 gap): thin tails sit above the
  90th percentile and cannot drag it.
- A prevalent-and-persistent high-income configuration (~12% of households
  at ~2.2x mean income, near-permanent) moved the estate p90/50 from 1.75
  to 2.74 at fixed parameters, with fertility, young liquid wealth, and
  aggregate ownership also moving toward their targets, and the estate
  median and composition share moving away (a refit would have to recover
  them). Tommaso is skeptical of adding heterogeneity as a "solution" —
  treat this as a measured fact about that configuration, not a
  recommendation.
- A shifted-lognormal process, log(e−s) AR(1) (a friend's suggestion), is
  the mathematically correct way to bound the support at s (estimable,
  smooth); it inherits the same unresolved question as truncation: what
  justifies the level of s. The outside assessment adds a units caution
  (compare post-transfer resource guarantees, not raw ratios) and a
  mean-preserving intercept formula.
- A means-tested transfer floor (Hubbard–Skinner–Zeldes 1995; De Nardi–
  French–Jones 2010, eq. 10 — both verified verbatim from the papers) is
  the literature's standard device and works mechanically. At the measured
  US guarantee for childless households (~10–20% of mean income) it does
  NOT clear our ~32% bundle, so it forces c_bar_0 down and reopens the
  fertility calibration; at a guarantee equal to our bundle it becomes a
  large program (and HSZ's ~$7,000 / ~31%-of-median figure is a bundle for
  a mother with two children — verified — not a precedent for a universal
  childless floor at 30%).
- Equivalence scales (replacing the Stone–Geary intercepts) would remove
  the feasibility problem at its root (no intercept, no empty budget sets)
  but re-identify the entire fertility block — a preference redesign, not
  an income-process fix.
- Own PSID estimation exists and is reusable: on transfer-inclusive family
  income, persistent+transitory gives rho 0.975 / sigma_eta 0.177 /
  sigma_transitory 0.413; AR(1)-only gives 0.944 / 0.285 (pipeline:
  `code/data/psid_followup_mar2026/estimate_intergen_income_entry_targets.R`,
  block 1; person-bootstrap SEs). These are close to Boar–Gorea–Midrigan's
  published disposable-income process.

## The underlying tension, stated once

Every route to realistic income risk runs through the same fork: the
subsistence bundle (~32% of mean income, driven by the calibrated c_bar_0)
is larger than any measured safety-net guarantee for childless households.
So either (i) something in the model supplies resources at the bottom
(floor/transfers — with a level that must be justified), or (ii) the bottom
of the income support is bounded (with a level that must be justified), or
(iii) c_bar_0 comes down and the fertility block re-identifies, or (iv) the
preference structure changes so no hard bundle exists, or (v) the current
low-variance process stays, defended as-is (the current shipping position).
Heterogeneity at the top is a separate question that touches the wealth
tail, not the bottom feasibility.

## Key artifacts (read in this order)

1. `docs/model/intergen_income_risk_feasibility_decision_memo_20260717.md`
   — verified survey of floors/levels/precedents + outside assessment.
2. `~/Desktop/income_risk_advice_20260717/2_PROBLEM_STATEMENT.md` and
   `10_topstate_probe_results.txt` — the falsifications and both probe
   rounds, with numbers.
3. `latex/quantification_rewrite_review.tex`, income paragraph + footnote —
   what the paper currently says (circulated to advisors).
4. `CALIBRATION_STATUS.md` top sections — the M5 calibration this all sits
   on (loss 9.044; established block 3.166; old-age ownership 0.954 vs
   0.764 is a separate disclosed miss, not part of this task).
5. `code/data/psid_followup_mar2026/estimate_intergen_income_entry_targets.R`
   — the estimation pipeline (block 1 = income process).

## Your first message back to Tommaso

Summarize: the problem in your own words, the fork above, and the evidence
status of each attempted route. Ask nothing except whether your
understanding is right. Then Tommaso decides the direction with you.
