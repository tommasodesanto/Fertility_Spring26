# HANDOFF: E-series strand — full state and marching orders (2026-07-23)

For the next Fable session. Launch Claude from
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26` (the
previous session was launched from the Job Market folder by accident; all
work still landed here, but memory keying and relative links misbehaved).

## 0. Division of labor — read this first

Two agents work this repository in parallel:

- **The other agent owns the M strand**: the one-shot Stone–Geary production
  model (`code/model/intergen_housing_fertility_optimized/`), its one-shot
  14-moment ledger, and its repairs. On July 23 it committed a wealth-timing
  repair (`d0dce5e`, documented in `a85ee17`, figure in `2638678`) and
  overnight calibration waves (`5b9e3ae`, `fa55d50`, `ed4cb12`).
- **This session owns the E strand**: the equivalence-scale / sequential-
  fertility architecture in `code/model/intergen_eqscale_seq_optimized/`
  (an isolated fork; production untouched).

Coordination rules: do not edit the M packages or their calibration
contracts; before editing shared files (`CALIBRATION_STATUS.md`, data
scripts), run `git log --oneline -10` and `git status -sb` to see what the
other agent changed; write dated sections rather than rewriting theirs.

## 1. Mandatory startup (per AGENTS.md)

Read in order: `memory/AGENT_MEMORY.md`; latest `memory/daily/*.md`;
`CALIBRATION_STATUS.md` (the July 23 "consolidated state of play" section
is the recap of everything below); then this strand's plan of record:
`docs/model/eqscale_calibration_reconciliation_20260722.md` (audit of the
old ledger, the reconciled target system, the July-23 decision record in
its Section 3.6, and the phased implementation plan).

## 2. What the E-series is, and what was built

Architecture (all verified in code): preferences are Cobb–Douglas CRRA
(sigma=2) with NO Stone–Geary intercepts — children act through (i) an
expenditure-share tilt `alpha(n) = clip(alpha0 - (delta_jump +
delta_alpha*n), .05, .95)`, (ii) a needs scale multiplying flow utility,
(iii) additive per-child utility `psi*n` — all only while a child is at
home. Fertility is sequential: childless choose {wait, try}; parity-n
parents with a child at home choose {stop, try for n+1}; attempts succeed
with external age-declining fecundity `pi_j = clip(1 - 0.02*exp(0.134*
(age-18)),0,1)`, zero from 45 (within 1-3pp of Leridon anchors; formal LS
fit is a pending half-day). Income risk is honest (no feasibility floors
needed). One housing market, owner rungs {2,4,6,8,10}, renter cap 6 rooms,
down-payment phi=0.80, supply elasticity 1.75.

**Literal parity (L4), built July 22 by this strand** — gated, default-off,
bitwise-nested (golden V/g hash check), 102 tests incl. 9 new in
`tests/test_literal_parity.py`:
- `n_parity=4` literal birth states 0/1/2/3+ (override plumbing existed;
  maturation matrix already routes nn>=2 correctly);
- generalized per-parity attempt margin (side-channel slot nn-1 of
  `fert2_probs`; shape unchanged up to n_parity=4; payload contract intact);
- KFE third-birth split with all at-risk pools snapshotted BEFORE any birth
  flow lands (no same-period chaining, tested);
- `entrant_conversion_factor` (0.5 = two matured children form one entrant
  household; the old x2 tfr convention carried this implicitly, so
  demographics are unchanged);
- `fertility_units="literal_topcode"` with `tfr_top_bin_weight` (top-coded
  3+ bin weighted at its empirical mean, 3.4 provisional pending a CPS
  measurement); `child_bin_high_cutoff=3` literal room bins;
- `get_completed_fertility` reads parity for 3+ while preserving the legacy
  mapping verbatim on nn<=2 cells (bitwise nesting requirement — a max()
  rewrite broke V-hash equality via zero-mass bequest-table cells; the
  lesson: nesting means bitwise, including unreachable cells).

## 3. Results so far (with their exact status)

- **E2** (bin-parity, collected July 20): strict loss 5.486 on the 15-moment
  system, all interior. **E3** (literal parity, July 22 overnight, 8/8
  chains): strict 2.806. **E3b** (continuation): 2.294. Torch runs:
  E3 array `14583403`, collector fix in `654ab11` (the collector's contract
  check now accepts arm `E3_L4`), continuation `14620597/14620598`.
  CAVEAT: all these losses contain wealth rows measured with the invalid
  hybrid timing (Section 4) — rankings are suggestive, nothing is certified
  until re-measured. Also the E3/E3b winners hit completed fertility by
  gutting childlessness (0.03-0.04 vs 0.188) with psi ~ +2.7-2.9, kappa_f ~
  13-16, gamma_e ~ 0: fertility carried by logit noise — not acceptable
  economics regardless of loss.
- **Fertility frontier scan** (147 cells, psi x kappa_f x gamma_e at the E3
  winner, L4 on; `run_fertility_frontier_scan.py`, output
  `output/model/eqscale_fertility_frontier_20260722/`): NO cell reaches
  completed fertility 1.918 AND childlessness 0.188 jointly; best joint
  cells sit at CF~1.26, p0~0.17 with 3+ share collapsing. Diagnosis: one
  psi plus iid per-period taste noise ties the entry margin to the
  progression margin; a LINEAR needs scale cannot separate them. This
  finding is computed from mass objects — unaffected by the timing bug.
- **Family-gap scan** (25 cells over delta_jump x delta_alpha at the E2
  winner; `run_family_gap_scan.py`, output
  `output/model/eqscale_family_gap_scan_20260722/`): the parent/nonparent
  ownership gap rises monotonically in both tilts from 0.043 to 0.160
  (target 0.168) — the mechanism exists without a committed-housing device;
  costs appear at the alpha-clip. Also mass-based, clean.

## 4. The wealth-timing bug (audited and confirmed July 23)

The wealth moments were computed on a hybrid distribution
(`asset_current` = `assign_current_cross_section_to_beginning_assets`):
mass reassigned to the NEWLY chosen tenure while retaining INHERITED
beginning-of-period liquid wealth b. A buyer keeps pre-purchase b next to
the new house (double count); a seller keeps the pre-sale mortgage next to
a renter label. Worst at old ages (mass selling). At the one-shot winner:
estate p90/p50 3.55 hybrid vs 1.91/2.15 coherent; bequest flow 0.009 vs
0.0126/0.0120; wealth/earnings 6.20 vs 6.01/5.93. Affected in the E-series:
rows `young_childless_renter_liquid_wealth...` (8),
`old_total_estate_wealth..._median_7684` (12),
`old_nonhousing_ge_1x_income_share_6575` (13) — and hence beta, theta0,
theta1 in every calibration that targeted them.

**Adopted conventions (author-confirmed):** cross-sectional wealth stocks
at BEGINNING of period (b_t + p*H_t with incoming labels; in a stationary
population this equals the end-of-period cross-section — the survey
analogue); bequest flow and any decedent-based estate object AT DEATH
(post-saving b' + p*H_held, death-probability-weighted, terminal age
included); the living-PSID p90/p50 (3.448) stays cross-sectional; housing
in estates stays gross of sale costs (stated convention). Principle: one
rule — match the instant at which the data counterpart is measured.

**FIRST TASK OF THE NEXT SESSION:** audit the other agent's repair commit
`d0dce5e` (`git show d0dce5e --stat`) and determine whether it covered ONLY
`intergen_housing_fertility_optimized` or also the E-package. The E-package
(`intergen_eqscale_seq_optimized/solver.py`) has the identical
`asset_g=asset_current` pattern at ~lines 4153/4179/4741. If not covered,
port the repair to the E-package under the same gates (golden bitwise
check where applicable, tests, tiny-config verification), then recompute
the E2/E3 target tables at their winners so the wealth rows are honest.

## 5. Decisions made July 23 (all author-confirmed; provenance in
`docs/model/eqscale_calibration_reconciliation_20260722.md` Section 3.6)

1. **Equivalence scale: imposed, not estimated.** e(n) = ((2+0.7n)/2)^0.7
   relative to a childless couple — the scale of Scholz–Seshadri–
   Khitatrakun (2006, JPE 114(4), p.619: "Equivalence scale.—This is
   obtained from Citro and Michael (1995) and takes the form (A_j +
   0.7K_j)^0.7"), also used by Borella–De Nardi–Yang (2023 ReStud; NBER
   w26097 pp.10 and 63). Both SET it; neither calibrates it. PDFs:
   `docs/reference/scholz_seshadri_khitatrakun_2006_jpe.pdf`,
   `docs/reference/borella_denardi_yang_nber_w26097.pdf`. Consequence:
   gamma_e leaves the free list (12 -> 11 free); the two parity-binned CEX
   increments (0 vs 1-2: available; 1-2 vs 3+: from the July-22 pass)
   become overidentification checks of the imposed curvature. Square-root
   scale = declared robustness alternative. NOTE the adaptation: BDY feed
   the scale an age profile of average children; we feed actual parity n —
   say so in the text. Concavity rationale: first child most expensive =>
   entry margin stays costly (childlessness survives) while continuation
   is cheap (2-3-kid families) — the decoupling the frontier scan showed
   the linear scale cannot deliver.
2. **Income process: imported pair.** Floden–Linde (2001, RED 4(2),
   406-437; PDF read directly): US rho=0.9136 (se 0.0090), innovation
   variance 0.0426 => sigma_eps=0.206 annual — passed through the HSV
   (2017, QJE 132(4)) log-linear tax function, which preserves the AR(1)
   exactly and scales sigma by (1 - tau_p). HSV WP gives tau_p in
   0.155-0.185; PIN THE PUBLISHED QJE BENCHMARK VALUE when wiring (do not
   trust memory; the commonly cited number is 0.181 but verify the
   published table). Result: rho=0.9136, sigma~0.17, 5-state Rouwenhorst
   (annual-to-4-year aggregation in `income_process_overrides` verified
   correct: rho^4, sigma*sqrt(sum rho^{2k})). Rationale: the model's flat
   payroll tax shifts levels only — a flat tax cannot compress log risk;
   the HSV wedge stands in for the progressivity the model lacks. This
   replaces the incoherent current pairing (M5 persistence 0.9602 with SS
   sigma 0.20).
3. **Family ownership gap (0.168) enters the hard loss** as the single
   ex-ante overidentifying row.
4. **Tenure moment:** deterministic state-conditional E[p(1-p)] analogue
   for the next diagnostic round (it LOWER-bounds the data object — kappa_T
   biased up — disclosed); the seeded simulated cross-fitted PSID
   replication is a hard gate before any paper calibration. Data target:
   model-feasible Brier 0.1176 (re-measured July 22 without married/year).

## 6. The target ledger after these decisions

**11 free parameters:** beta, theta0, theta1, alpha0, delta_jump,
delta_alpha, psi, kappa_f, chi, kappa_T, H0.
**12 hard rows:** wealth/after-tax-earnings 6.90 (borrowed DNY — remeasure
before paper); bequest flow/wealth 0.0088 (borrowed Gale–Scholz); estate
p90/p50 living 76-84 = 3.448 (PSID, cross-sectional timing); LES
expenditure slope -> alpha0 = 0.733 (CEX); first-child rooms 0.664
(horizon-0 diff-in-diff — the active Markov analogue; the 12-year variant
is legacy-path-only); 3+ vs 1-2 rooms gap 0.368 (literal bins); completed
fertility 1.918 (rename the moment from "tfr"; = parity mean with top bin
at its measured CPS weight); childlessness 0.188; ownership 30-55 0.5755;
tenure Brier 0.1176; aggregate rooms 5.780; family ownership gap 0.168.
**Validation (out of loss):** parity-binned CEX consumption increments
(now checks of the imposed scale), LES intercept + price coefficient
(homotheticity test the spec is allowed to fail), young liquid wealth
0.179 (expect a miss under honest risk — disclosed), owner-renter rooms
2.419, childless renter rooms 3.805, owner >=6 share 0.596, own_rate_2534
0.341, old-age ownership 0.764, old nonhousing share 0.608, estate median
6.501, timing battery (mean age at first birth, share 30+, hazard
profiles, PP 1->2 ~0.77 literal, PP 2->3), tenure switch rate,
fertility-income gradient (to build: polices noise-driven fertility).
**Externals:** theta_n=0; sigma=2; SSK scale; FL x HSV income; fecundity
(0.02, 0.134, 45) pending LS fit; SSA survival from 66; phi=0.80; sale
cost 6%; depreciation 1.1%/yr; property tax 1%/yr; r=2%/yr; menus
{2,4,6,8,10}/hR_max=6; eta=1.75; PSID 18-24 entrant wealth; 18-yr child
duration; L4 conventions (n_parity=4, entrant factor 0.5, literal bins,
top-bin weight measured).

## 7. Critical path (in order; no step may be skipped)

0. Verify/port the wealth-timing repair to the E-package (Section 4).
1. Implement the SSK scale as a gated `eqscale_form` option (default
   linear, bitwise-nested; "power" activates ((2+0.7n)/2)^0.7 — note at
   sigma=2 the utility multiplier equals the scale itself; document the
   sigma-dependence once).
2. Wire the FL x HSV income external (pin published HSV tau first).
3. RERUN the 147-cell fertility frontier under the imposed scale — the
   decisive test: does concavity unlock CF 1.918 + childlessness 0.188
   jointly? Pre-registered: if yes -> proceed; if no -> the problem is the
   preference structure, and the (rejected for now) options are per-parity
   psi or permanent taste heterogeneity — the author explicitly dislikes
   parameter-adding; bring evidence, not knobs.
4. Data passes: CPS top-bin weight; A6 timing-moment targets
   (cohort-consistent, NCHS/NSFG/PSID); fecundity LS fit; the CEX 1-2 vs
   3+ increment on literal bins if not already in the July-22 outputs.
5. Freeze the ledger (SE-based weights where measured; declared synthetic
   SEs for borrowed rows), Jacobian smoke at the seed (11 columns, rank
   must be 11 at 1e-3), reachability precheck at the seed (any row >3x off
   triggers review), THEN one 8-chain Torch recalibration under Long-Run
   Search Safety (exact-loop smoke first, strict repeats, collectors).

## 8. Standing user directives (repeated failures cost trust — obey)

- **Explain like to an advisor**: define every object in one plain sentence
  before naming it; no session shorthand; one worked example; don't
  inflate pedigree (a report is a report); state established vs
  conjectured. Now codified in AGENTS.md/CLAUDE.md "Communication
  Standard". Absolute paths in chat, always.
- **Do not add parameters to fix fit problems.** The author rejected a
  psi-split proposal on these grounds; the concave scale won because it
  costs zero parameters. Off-the-shelf, famous-paper imports beat
  home-brew estimation wherever possible.
- **Completed fertility is PARAMOUNT.** Never demote it. (And it equals
  cohort TFR only in steady state; the data target 1.918 is the CPS
  children-ever-born cohort stock, never the period TFR ~1.6.)
- **Verify sources against the actual paper PDF** before citing numbers
  (this session: wrong NBER WP number caught only by opening the file;
  CALIBRATION_STATUS carried a stale c_bar_n bound claim; the July-20 spec
  note's "gamma_e reproduces OECD" claim was units-inflated).
- Loose search losses are never results; only strict twice-repeated
  collector winners count. Never compare losses across target systems.
- Report full target-fit tables (every moment, target, model, gap, weight)
  and every parameter with bounds when reporting any calibration.

## 9. Where things are (absolute paths)

- Repo: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`
- E-package: `code/model/intergen_eqscale_seq_optimized/` (L4 tests in
  `tests/test_literal_parity.py`; chain runner `run_e1_chain.py` with env
  gate `E3_L4=1`; scan drivers `run_family_gap_scan.py`,
  `run_fertility_frontier_scan.py`; collector `collect_e1.py`).
- Cluster: `code/cluster/submit_intergen_e3_l4.sh`,
  `submit_intergen_e3b_continuation.sh`, `submit_intergen_e3_gap_scan.sh`,
  `submit_intergen_e3_fert_frontier.sh`; wrapper `code/cluster/torch.sh`
  (account `torch_pr_570_general`, remote
  `$SCRATCH/projects/Fertility_Spring26`, sync via rsync of exact files —
  the remote is NOT a git clone). Slurm chains: afterok on an array
  requires ALL tasks to succeed; a failed collector auto-cancels
  dependents (bit us once: fixed in `654ab11`).
- E-series outputs: `output/model/eqscale_seq_l4_recalibration_20260722/`
  (production + report + continuation + report_continuation),
  `output/model/eqscale_fertility_frontier_20260722/`,
  `output/model/eqscale_family_gap_scan_20260722/`.
- Reference PDFs: `docs/reference/scholz_seshadri_khitatrakun_2006_jpe.pdf`
  (scale on p.619), `docs/reference/borella_denardi_yang_nber_w26097.pdf`
  (pp.10, 63). Floden–Linde and HSV are cited with links in the
  reconciliation doc Section 3.6.
- Plan of record: `docs/model/eqscale_calibration_reconciliation_20260722.md`.
- This handoff: `docs/prompts/HANDOFF_fable_eseries_20260723.md`.

## 10. Known tensions to expect (do not be surprised; do not hide)

- After the timing repair, the estate-dispersion tension returns (coherent
  p90/p50 ~1.9-2.1 vs 3.448 at the audited one-shot winner). The E-series'
  honest income risk is the structural answer (direction robust, levels to
  be re-measured); if it falls short, the nonhousing-saving margin
  discussion from `docs/model/intergen_bequest_balance_sheet_fable_review_
  20260716.md` reopens.
- Young liquid wealth will overshoot under honest risk (disclosed
  validation miss; FL x HSV's sigma~0.17 < 0.20 trims it somewhat).
- The renamed completed-fertility row plus childlessness must bind
  TOGETHER; watch for optimizer configurations with psi>0 large and
  kappa_f>10 (noise-driven fertility) — gate acceptance on the fertility
  block's economics, not just the loss.
