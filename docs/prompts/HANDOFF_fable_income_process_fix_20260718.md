# HANDOFF: fix the income process properly (single-task session)

Audience: a fresh Fable session working with Tommaso. This is the one task.
Read `CALIBRATION_STATUS.md` top sections first; general session context is
in `docs/prompts/HANDOFF_fable_next_session_20260718.md` (operational
gotchas at the bottom of that file apply here too).

## The problem, in four sentences

The model's income process (annual AR(1), rho 0.960, sigma 0.0645, 5-state
Rouwenhorst) is provisional: the persistence is literature-consistent but
the innovation variance is ~1/3 of PSID estimates, retained because larger
values are infeasible — under the Stone–Geary minima (c_bar_0 + rent ~32%
of mean income) with no safety net, households at the bottom persistent
state cannot afford the required bundle (verified by direct probes; the
paper currently defends this honestly and flags it provisional). Separately,
the income process carries no permanent heterogeneity, which is why the
model has no late-life wealth tail (estate p90/50 = 1.75 vs 3.45,
untargeted). Fixing the process means fixing both ends: real heterogeneity
at the top, and an honest treatment of risk at the bottom. Everything
below is evidence-backed; do not re-derive it.

## What is established (artifacts; do not re-litigate)

- Falsified: richer entry endowments (flow problem — zero-debt entry still
  breaks at 22); innovation-mean shifts (support unbounded; no-op after
  normalization); bottom truncation alone (bound level is either measured
  ~0.15 — below the bundle, fixes nothing — or circular at 0.32); rare
  Castaneda top states (<1% mass moves estate 90/50 by +0.06 at most).
- Validated: prevalent-persistent high earners — ~12% mass at ~2.2x mean
  income, near-permanent — move estate 90/50 from 1.75 to 2.74 AT FIXED
  PARAMETERS, with TFR, young liquid wealth, and aggregate ownership all
  moving toward targets (estate median and composition share move away;
  recoverable in refit). Probes + all falsifications:
  `~/Desktop/income_risk_advice_20260717/` files 2 and 10;
  full survey and outside advice:
  `docs/model/intergen_income_risk_feasibility_decision_memo_20260717.md`.
- Own PSID estimates (all-transfer-inclusive family income; pipeline
  `code/data/psid_followup_mar2026/estimate_intergen_income_entry_targets.R`,
  block 1): persistent+transitory rho 0.9749 / sigma_eta 0.1769 /
  sigma_transitory 0.4128; AR(1)-only 0.9436 / 0.2846. These bracket the
  literature and independently match Boar–Gorea–Midrigan (0.964, 0.150).
- Feasibility mechanics: the current calibration sits exactly ON the
  feasibility frontier (a 2% shave of the bottom state breaks it); entry
  censoring exists (`entry_wealth_censor_to_frontier`, age-18 only) and the
  interior gate is untouchable (July-11 rule).

## The plan (staged; Tommaso decides at each gate)

**Stage A — first-stage estimation (data, ~1 day, no model code).**
Estimate a two-type income structure from PSID, by education (college /
non-college), in the transfer-inclusive family-income concept the model
uses: (i) type shares; (ii) type-specific age–income profiles (normalized
so the aggregate mean profile is unchanged); (iii) within-type residual
AR(1) (rho_g, sigma_g) by the same autocovariance method already coded in
block 1. Template: HSZ (1995) estimate income processes by education
(published Table 2: persistence ~0.955/0.946/0.955, innovation variances
~0.033/0.025/0.016 — cite the published table, not the WP). Rationale for
types over big shocks: Guvenen's heterogeneous-profiles point — permanent
profile differences masquerade as ultra-persistent shocks. Deliverable: the
estimated objects with person-bootstrap SEs + a one-page summary for
Tommaso's sign-off. Extend the existing R pipeline; person-cluster
bootstrap, seed 20260715, same conventions.

**Stage B — specification choice (Tommaso's call, with these defaults).**
Default: 2 types x 5 within-type states (10-state chain; ~2x solve cost),
types permanent for life. Budget fallback: 6-state chain with a
near-permanent top block matched to the college (level, mass, within-type
persistence) — cheaper but exposes the p90-inside-the-atom artifact, so if
chosen, report estate p85/p90/p95 and masses either side of 10%. Within-type
sigmas: use the Stage-A estimates IF the feasibility probe passes at the
bottom type; if not, bound the non-college sigma at the computed feasibility
frontier and DOCUMENT the bound (this is the honest interim treatment of
the bottom — the full bottom fix, floor-plus-c_bar_0-re-estimation, is a
separate decision that should wait until the type refit shows where c_bar_0
lands). Renormalization must hold the bottom state's LEVEL fixed (the
mean-one renormalization trap broke two probe cells; fix is in the probe-3
construction).

**Stage C — one recalibration (overnight).**
Same 15 moments / 14 parameters; income process FROZEN from Stage A (zero
new estimated parameters — this is the identification discipline; NEVER let
the estate ratio identify the income objects, that is circular). Nested
seed, 8 chains x 3:55, strict collector with all gates machine-encoded
(clone the M5 scripts; the M5 collector is the template). Predeclared gates
BEFORE launch, including: established-12 within tolerance of 3.166; both
late-life wealth targets recovered; young liquid wealth within band; and
the counterfactual-stability rule — recompute the two policy experiments
(grant; grant+tax) and require the headline effects to keep sign and order
of magnitude (M5 reference: +2.1%/+0.90%; +2.7%/−18.4%). Full free/fixed
parameter table in front of Tommaso before launch.

**Stage D — verdict and paper language.**
If gates pass: promote as the paper's income process ("household income
combines two permanent education types with within-type persistent risk,
estimated from the PSID"), keep the feasibility footnote for the bottom
(now scoped to within-type risk), and report the estate tail diagnostic —
expected to be materially closer. If the estate tail improves but gates on
established moments fail: report as robustness, keep M5 as baseline. If the
policy effects flip sign: STOP — that is a finding about the paper's
headline results, and the transfer-floor/c_bar_0 redesign moves from parked
to necessary.

**Also in this session (cheap, independent): the feasibility appendix.**
One script over a (rho, sigma) grid recording dead-node mass and
first-failure age; turns the paper's feasibility claim into a figure/table.
~1 hour; the probe machinery is described in the decision memo (scratchpad
scripts are gone — rebuild from the recipes there; the M5 contract template
is in `probe`-style: args -> arm_contract -> common_overrides -> solve).

## What Tommaso decides, explicitly

1. Stage A sign-off: are the estimated type objects (shares, premia,
   within-type risk) the right empirical anchors?
2. Stage B: 10-state (right object, 2x cost) vs 6-state (cheap, atom
   caveat); and accept the within-type bottom bound if the frontier binds.
3. Stage C: the predeclared gates and the launch.
4. Stage D: promote vs robustness vs stop.

## Hard rules carried over

Nb=120 standard (no 240 talk). Deterministic replays and probes run
locally; the calibration search runs on Torch (needs Tommaso's live ssh
session — check before promising). Present plans and WAIT before launching.
Document at every step: run contract before launch, CALIBRATION_STATUS at
launch and collect, commit protectively. The word "parity" is banned. The
full parameter table goes in front of Tommaso before any launch.
