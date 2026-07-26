# E5b calibration autopsy: what is actually failing

2026-07-26. Diagnosis at the certified E5b winner (loss 385.1). Extraction:
`output/model/eqscale_seq_e5b_autopsy_20260726/` (verification gate passed at
machine zero on all twelve moments; three GE solves). Read with the probe
verdicts of 2026-07-25.

## Summary

The four visible failures are not four problems. They collapse into two
upstream defects:

1. THE MODEL HAS NO PERMANENT HETEROGENEITY. This one defect produces the
   childlessness miss, the wealth-dispersion miss, and most of the aggregate
   wealth miss -- 84 percent of the remaining loss.
2. THE ENTRY DECISION HAS ALMOST NO PRICE EXPOSURE AND NO READINESS
   DYNAMICS. This produces the timing-shape miss and the zero policy
   elasticity.

## Defect 1: no permanent heterogeneity

Evidence.
- Model cross-sectional gross labor income p90/p50 = 1.61. The Floden-Linde
  process is a RESIDUAL wage process -- estimated after controlling for
  observables -- and the model uses it as the ONLY source of income
  differences. There is no fixed-effect layer underneath it. (Stationary
  log-sd of the continuous process is 0.415, implying p90/p50 about 1.7;
  the five-state chain delivers 1.61, so discretization is NOT the issue --
  the process itself is this compressed. Empirical earnings p90/p50
  including permanent heterogeneity is far above this; measure the exact
  PSID number, but the direction is not in doubt.)
- Wealth p90/p50 at 76-84 is 2.07 against 3.45 -- and the model's own
  internal ratio (5.61/2.67 = 2.10) shows wealth dispersion is already
  amplifying income dispersion. No bequest-curvature parameter can stretch a
  1.61 income distribution into a 3.45 wealth distribution; theta1 was
  convicted of a crime committed upstream.
- Aggregate wealth 5.45 vs 6.87: the missing top tail is missing aggregate
  wealth, so the level miss and the dispersion miss are partly the same
  miss.
- Childlessness is FLAT across income states (0.0807 to 0.0797, first to
  fifth state) and nearly flat in wealth (0.067 to 0.091). With ex-ante
  identical households, childlessness can only be noise plus biology: chosen
  0.030 + clock 0.051 = 0.080, ceiling ~0.12 on the whole kappa_E frontier,
  against 0.188 in the data with strong permanent gradients.

Implication. A permanent type layer (income fixed effects, optionally with
correlated psi) is not a fertility patch -- it is the missing bottom layer
of the income process, and it attacks childlessness, p90/p50, and aggregate
wealth simultaneously. This also resolves the author's Floden-Linde doubt
precisely: FL is a defensible estimate of the WITHIN component; the model
omits the BETWEEN component entirely. The planned own-PSID estimation should
measure both components (ledger item, income-process review).

## Defect 2: the entry margin

Evidence.
- Static decomposition (median young renter budget): the utility increments
  of a first child are psi +0.585, scale -0.083, tilt +0.014. A 10 percent
  rent increase moves the total by -0.0033 -- about half a percent of the
  child surplus, against taste-noise scales of 1.25 and 3.55. Housing enters
  the first-birth decision two to three orders of magnitude too weakly for
  prices to matter.
- GE confirmation: H0 +-7 percent moves equilibrium prices +-2.6 percent and
  ownership by +-1.6 points; tfr moves +-0.0006 and mean first-birth age by
  0.005 years. Tenure responds to prices; fertility does not.
- Timing shape: the model's attempt hazard RISES monotonically with age
  (0.23 at 18 to 0.57 at 42) and the realized first-birth distribution
  declines monotonically from 18. The data rise to a peak at 19, dip, and
  form a second plateau at 25-30 -- a readiness hump the model cannot
  produce with a single logit scale. Bin by bin: model has MORE first
  births than the data at 18-25, FEWER at 26-33, and a late tail three
  times the data (9.7 vs 3.0 percent of first births at 38+).
- The 38+ tail is priced by fecundity: our schedule still gives pi = 0.50
  at 42, while the Sommer/Trussell-Wilson object is ~97 percent permanently
  infertile by 45. The unfitted late-fecundity constants (open item #37)
  are not innocuous after all -- they are the direct source of the late
  tail. My 2026-07-25 conclusion that the fecundity refit "moves nothing"
  addressed the 30/35/40 levels and missed the tail; retracted for ages
  38-45.

Implication. The timing fix has two named parts: a readiness/arrival
process for the missing 25-30 hump, and a corrected late-fecundity tail for
the 38+ excess. The elasticity fix is separate and quantified: any housing
gate must raise the price-sensitive share of the entry decision by roughly
two orders of magnitude; the tilt channel cannot be scaled up to do this
(the 0.66-room increment caps it), so the gate must run through a discrete
margin (tenure premium with children, or requirement), to be designed only
after defect 1 and the timing parts are in.

## Loss forensics

Childlessness contributes 200 of 385; p90/p50 contributes 108; share30+ 61;
aggregate wealth 13; everything else under 2 each. Under both counterfactual
reweightings (share30+ SE doubled; timing rows at 5 percent of target) the
ranking is unchanged -- childlessness and p90/p50 dominate regardless. The
declared-SE worry does NOT drive the verdict; the model does. The weights
discussion can wait until defect 1 is addressed.

## Proposed order of work (nothing approved yet)

1. Measure the between/within income split in the PSID (extends the
   registered income-process review) and the empirical childlessness
   gradient by income/education (CPS; M-side or E-side builder).
2. Permanent income types (2-3 fixed-effect levels re-using the existing
   state machinery), refit. Expect movement on childlessness level and
   gradient, p90/p50, aggregate wealth.
3. Fecundity tail fit (task #37, now with a specific job: kill the 38+
   excess) plus readiness arrival for the 25-30 hump; refit.
4. Only then: the housing-gate design, judged against the two-orders-of-
   magnitude requirement, with Stone-Geary re-entry as the stated fallback.
