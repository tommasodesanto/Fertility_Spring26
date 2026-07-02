> Tracked mirror of the Phase-7 final report. Full evidence set (STATUS, COMMANDS, NOTES, CODEX_TASKS, HANDOFF_RESTART, SIZE_MECHANICS_MAP, SUPPORT_AUDIT_READOUT, MEASUREMENT_MAP_AUDIT, PRICE_SCHEDULE_DEFORMATION_READOUT, PROPERTY_TAX_SIZE_RECONCILIATION, SHADOW_REBATE_ACCOUNTING) lives in output/model/fable_size_mapping_audit_20260701/ (gitignored) and in the isolated copy ~/Desktop/Projects/Fertility/Fertility_Spring26_fable_size_mapping_audit_20260701 (own git repo, final commit e94b14a).

# FINAL_FABLE_SIZE_MAPPING_DIAGNOSIS — can the current one-market model deliver the fertility direction through size access?

2026-07-01. Fixed theta = current best (loss 6.0014895774, p 0.6759498581,
`candidate_replacement_post_audit_v1`), all work in the isolated copy
`Fertility_Spring26_fable_size_mapping_audit_20260701` (source @4ca392f,
dirty tree; baseline reproduced exactly). Evidence: Phases 0-6 deliverables
in this folder; ~80 GE solves total; drivers + run logs committed in the copy.

## Executive diagnosis (10 bullets)

1. **The answer to the central question is NO for the size-access channel and
   a qualified YES for the direction.** Δτ_H>0 ⇒ Δp<0, ΔyoungOwn>0, ΔTFR>0
   holds robustly in the canonical model (three of the four desired signs);
   ΔPr(H≥6 | fertility-marginal young)>0 fails under every defensible
   support, measurement map, and cost-unit mapping we could construct without
   changing model equations.
2. Capitalization is hardwired: the supply curve is written in the rental
   rate, so a property-tax rise lowers the asset price ~one-for-one at
   unchanged rent (−12.2%/−21.9% across all 18 tested worlds). The tax
   relaxes down payments without touching flow costs — that is the entire
   transmission.
3. Tax-induced young owners always buy the cheapest rung on any menu, and the
   cheapest rung is always below the parity-1 family floor (4.52 rooms at
   this theta): **tax-induced young ownership is fertility-dead starter
   ownership by floor geometry**, under every support tried (canonical, fine,
   ratio-2.5, ratio-3.3, renter-cap-8).
4. The model's fertility margin is not where the mechanism assumes: young
   childless renters carry a derivative mass of only ~0.004; the margin
   concentrates at older fertile ages where Pr(R≥6|birth) is already 0.84
   (0.97 in ACS-quantile units — satiated). The marginal size action is a
   6→8 upgrade (empirical top-third), where tax pass-through is positive but
   tiny (+0.002-0.004 quantile-mapped).
5. The July-1 deformation "revival" does not survive scrutiny: (i) the hook
   only priced the Bellman — the forward distribution moved wealth at base
   prices (bias: TFR understated by 0.16-0.19 at λ=0.5/1); (ii) after the
   consistency patch, ALL of the fertility gain comes from the `anchor2`
   normalization, i.e. an unfunded absolute subsidy to rungs ≥4. Neutral
   flattenings (average- or family-rung-cost-preserving) LOWER fertility and
   collapse young ownership (0.215→0.055).
6. The H mapping IS distorted — the quantile bridge says model rung 4 ≈ ACS
   4-5 rooms, rung 6 spans ACS 5-7+, rungs 8/10 are both ACS 7+ — but
   correcting units is not the blocked pipe: the R≥6 no-pass-through result
   is only mildly a measurement artifact (threshold placement), and no
   re-mapping makes the young size-access channel material.
7. The single strongest structural fertility lever found anywhere is the
   RENTER cap: hR_max 6→8 raises baseline TFR by +0.071 at fixed theta —
   an order of magnitude above any policy effect — because the parity-2
   floor (5.888) is cap-pinned. But it collapses the family-ownership block
   (own_family_gap 0.296→0.035 vs target 0.168; own_rate 0.525→0.411):
   **family-size rental supply and the ownership fit trade off one-for-one**
   in a one-market model with rent = user_cost·p.
8. Parent LTV is confirmed placebo across every world: 95-100% of its
   ownership response sits in zero-birth-probability states (Phases 2-4),
   and the Phase-6 accounting explains why — fertility-marginal young
   renters already afford the starter down payment; their gap is the family
   rung, which LTV relief on any rung does not target as well as cash.
9. The missing budget object is quantitatively real and correctly targeted
   only if young-targeted: rebating the model's own extra tax revenue to
   young (or young-renter) households bridges the H=6 down-payment gap for
   29-36% of the fertility-marginal band (88-91% of crossings at the family
   rung, not the starter rung). Uniform rebates are too dilute (zero
   crossings at tax2).
10. Two accounting bugs/artifacts discovered in passing: (a) the
    Bellman-only implementation of BOTH reduced-form hooks (price multiplier,
    old-owner sale wedge) — the prior sale-wedge failure is unreliable
    evidence; (b) stochastic child aging matures 22% of children within one
    period, creating age-26 "ACS-childless" empty-nest owners that flatter
    the childless-owner room moments.

## The ten required questions

**1. What worked without model changes.** Exact baseline reproduction; the
directionally correct triple (Δp<0, ΔyoungOwn>0, ΔTFR>0) for the property
tax at canonical config; threshold re-measurement (R≥7/top-third) revealing
small positive marginal size pass-through; the quantile units bridge; the
shadow-rebate accounting showing a well-targeted margin exists.

**2. What did not work.** Every attempt to make the tax (or parent LTV) move
derivative-weighted H≥6 access for young households: five supports, two
measurement maps plus threshold family, thirteen cost-unit deformations. The
neutral deformations actively hurt. The one 4/4 cell (anchor2 λ=0.25) is an
unfunded subsidy world with a 4× worse fit and a destroyed birth-event
housing response.

**3. Is the current H mapping likely wrong?** Partly. The triple identity
(services = rooms = cost units) is exactly what the code imposes (verified,
no bridge anywhere), and the quantile map shows real distortion at the top of
the support. But the audit shows fixing the mapping alone cannot deliver the
mechanism; the mapping error mostly misleads INTERPRETATION (which rung is
"family-size", where the top-half cut sits), not behavior.

**4. Is linear pH defensible once H is correctly mapped?** Yes — with the
support re-expressed in properly bridged units. DUE prices size linearly
within location; the ACS rent schedule's concavity (ε=0.3187) is about
unit-quality composition across bins, not evidence for concave asset pricing.
The audit's neutral deformations show concave pricing is not even useful
here. Keep pH linear; fix the units the support is stated in, and consider a
support whose top matches the empirical top (rungs 8/10 are currently both
"7+" — one of them is redundant in room units and absorbs cost mass).

**5. Is the multiplier hook enough for the next serious calibration?** As of
the consistency patch it is internally coherent (Bellman + KFE + wealth
stats), so it can carry unit re-normalizations of the SUPPORT during
calibration experiments. It is NOT enough for policy results: its wedge
revenue accrues to no one and its supply side does not exist. If the paper
needs size-dependent prices, they must be externally disciplined AND funded.

**6. Can property-tax capitalization generate fertility here after size
mapping is corrected?** Direction yes, mechanism no. TFR rises ~+0.02-0.03
via cheaper ownership options at birth states across all fertile ages and a
0.3-0.5pp fall in completed childlessness — not via young households moving
into family-size homes. Correcting the size mapping does not change this;
the young size-access leg fails everywhere defensible.

**7. Parent LTV.** Keep it strictly as placebo/contrast. Its footprint is
ownership in zero-birth states in every configuration tested; the Phase-6
targeting logic explains structurally why it cannot become the main policy.

**8. Government budget: necessary but not sufficient?** Exactly that. It is
necessary for interpretation (three sinks already leak resources: ψ, τ_H,
estate tax; the capitalization experiment is a pure levy), and the revenue
is large enough to matter (young-renter rebate ≈ 18% of eligible income at
tax2; ≈ the full H=6 down payment at tax3). It is not sufficient: without
targeting, crossings are ≈0; and crossing the dp threshold is an option-value
statement, not yet a birth.

**9. What would force a true model change?** (i) If a funded young-targeted
transfer (the minimal addition) still fails to raise TFR through the family
rung — that would indict the one-shot-completed-fertility valuation of the
birth-state option and only then justify architecture work; (ii) any paper
claim requiring rents and prices to move separately (they cannot: r ≡
user_cost·p) — that forces either segmented rental/owner markets or an
independent rental supply curve; (iii) any old-owner release/lock-in story —
requires basis/realization taxation, which does not exist and cannot be
faked Bellman-only (the wedge experiment is now known to be unreliable).
Sequential fertility is NOT indicated by anything in this audit: direct
large-house cost cuts move completed fertility strongly, so the one-shot
architecture is not the block (rule E).

**10. Exact recommended next run.** In order:
   (a) *Fixed-theta diagnostic (no model change, half a day):* re-run the
   tax + LTV suite at the canonical config with the fertility-margin
   decomposition split by AGE and PARITY-ALTERNATIVE to pin exactly which
   fertile ages/parities produce the +0.021 TFR — this sharpens the paper
   story about what capitalization does move.
   (b) *Minimal model change (one small object):* implement the funded
   transfer — budget identity G = τ_H·p·E[H_own·m] rebated per eligible
   state (young-renter default), wealth increment on the KFE AND Bellman
   sides symmetrically. Run tax2/tax3 with rebate on/off at fixed theta.
   The Phase-6 accounting predicts this is the first configuration with a
   fighting chance of Δτ_H ⇒ ΔTFR through the family rung. Budget: 6 GE
   solves plus code; smoke-test the eligibility mask against Phase-6 masses.
   (c) Only if (b) succeeds qualitatively: a disciplined recalibration in
   quantile-bridged units — replace the R≥6 target statement by the
   empirical top-third object (Phase 3), reconsider the redundant rung 10,
   and re-examine hR_max as an estimated (not imposed) object given its
   load-bearing role — with the identification ledger updated accordingly.
   No blind recalibration before (b).

## Decision rules

- **A (support/mapping fixes it): REJECTED.** No support or mapping
  correction makes the tax move derivative-weighted H≥6 access.
- **B (one-market scalar price too blunt for the size margin): ACCEPTED**,
  with precision: the scalar price moves the starter margin powerfully and
  the family-size margin not at all; the fertility-relevant size constraint
  in this model is the rental cap and the family-rung down payment, neither
  of which a uniform price level can target.
- **C (only extreme concave schedules work): PARTIALLY** — only the
  *subsidy-like* normalization works, only at mild λ, and it must be
  empirically disciplined as a funded policy to be usable; as a schedule
  claim it fails (extreme λ flips the sign).
- **D (nothing works): REJECTED** — the direction is deliverable today; the
  mechanism through young size access is not.
- **E (sequential fertility): NOT INDICATED.** The one-shot architecture
  responds strongly to direct size-cost variation; the block is targeting,
  not timing structure.

## Honest limitations

All results are fixed-theta; a recalibration under any variant would shift
levels (not, we judge, the targeting geometry, which is driven by floors vs
menu minima and by where the derivative mass sits). The derivative-weighted
objects treat the multinomial birth choice as a grouped binary logit (exact
only for common shifts). The shadow-rebate crossings are option-value
accounting, not behavior. The J17/Nb60/max_iter_eq=10 numerics leave ~1e-4
residuals; nothing here hinges on effects below ~5e-3.
