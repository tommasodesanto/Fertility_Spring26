# M4 standard-bequest calibration: independent audit (2026-07-16, evening)

Verdict: **C — provisional.** Numerically clean; scientifically conditional.
Auditor: Fable, per `docs/prompts/HANDOFF_claude_evaluate_m4_calibration_20260716.md`.
Five audit strands (contract/code, artifacts/chronology, tenure-kappa history,
authorization transcript, income-process provenance), all findings verified
against local artifacts; per-evaluation Torch records were unreachable (SSH
timeout) — three items rest on collector-gate logic plus local mirrors.

## Headline result (verified bit-level)

Winner: loss `4.401508692086581`, residual `1.50e-5`, `theta0=0.247497728306924`,
`theta1=0.4217108770366293`, established-12 loss `3.6497` (M1 comparable:
`3.841`), estate median `6.5046` vs `6.5013`, nonhousing median `2.0029` vs
`1.9082`. All five machine-checked gates pass; two strict repeats bit-identical;
one parameter vector across final_report / identification / diagnostics.
Artifacts: `output/model/intergen_standard_bequest_recalibration_20260716/`.

## Why it is provisional, not final

1. **theta1 is set-identified, not point-identified.** Weakest Jacobian
   direction loads 0.996 on theta1 alone; its column norm is 8.6–14x below
   theta0's; its derivative on the nonhousing median is exactly 0 and on the
   estate median −0.0067 — what little curvature it has comes from prime-age
   rooms moments. The reported value is bit-inherited from a profile cell at
   the *edge* of a 5-cell ladder centered on the failed six-chain winner
   (0.8229); **no evaluation below theta1=0.4217 is verifiable anywhere**, and
   the final "13-dim polish" winner was literally its first initial-simplex
   vertex (a beta-only perturbation, RNG-reconstructed). Profile: 0.4217 →
   4.4687, 0.5891 → 4.5715, then steep (0.823 → 12.65, 1.61 → 80.98).
2. **Six-chain search failure was budget, not geometry.** All chains
   time-starved (~215 of 1,000 evals in 75 min); two chains *starting* at
   theta1=0.25 finished eligible at losses 229.9 / 28.2. The nested theta0=0
   reference (307.68) fails through one row: estate median collapses to 2.463
   (contribution 303.1) — i.e., **theta0 is cleanly identified**; theta1 is not.
3. **tenure_choice_kappa=0 was never surfaced for M4 by name.** Real empirical
   trail (estimated ~0 repeatedly June 28–30; production [0,0.12] search froze
   it at 0; July-13 gradient points into the deterministic face; full-rank
   restricted system), plus the historical eps_ten design rule. But: never in
   the approving pre-launch exchange; no free/fixed table shown pre-launch;
   `CALIBRATION_STATUS` said "only theta_n=0 is externally restricted" (false);
   and a deeper basin at kappa=0.018 (6.401 vs 6.99, July-14 chain_c) was
   excluded on design grounds. Sweep evidence says positive kappa *lowers*
   aggregate and old-age ownership — the direction of M4's two biggest misses
   (own_rate 0.715 vs 0.575; old-age own 0.964, path acceptance False).
   Honest label: "M4 conditional on deterministic tenure choice."
4. **The income process is a numerical construction** (user finding, fully
   confirmed): `(rho, sigma) = (0.9601845894041878, 0.06453733259357768)` is the
   exact inversion of the provisional stay-or-redraw grid (rho4 = 0.85 by
   construction; stationary log var 0.05337 bit-identical), built 2026-07-10 as
   a discretization control (`tmp/rouwenhorst_matched_nb120_short_SOURCE_SNAPSHOT.txt`)
   and promoted to production the same night without external anchor, git
   history, or status entry. It carries ~25% of the Sommer–Sullivan (0.90,
   0.20) stationary variance and ~10% of the innovation variance. The July-10
   SS-parameter runs were abandoned after the young liquid-wealth moment
   exploded (1.90 vs 0.179, contribution 35.5) — the calibration selected the
   income process. The s_R "income dispersion rejected" diagnostic and the
   reachability "tail unreachable" verdict are both conditional on this
   process (working-life risk was never varied). Repo issue INC-1 already
   flags this "Open; urgent".
5. **Entrant wealth is partly circular** (user finding, confirmed): the age
   25–35 childless-renter quintile distribution is injected at age 18
   (`calibration.py:33-45`) and its mean `0.17922556` is bit-identical to the
   hard 25–35 target. Not an accounting identity — the moment ranged over
   [−0.18, +1.90] across parameter points with the same injection, so behavior
   dominates — but the initial condition decays slowest exactly under the
   current low-risk income process, and the July-10 "cannot absorb realistic
   income risk" verdict was contaminated by entrants starting too rich. The
   two repairs (real income risk up, honest 18–24 entry wealth down) are
   offsetting on this moment.
6. **The reinstated nonhousing median is grid-discrete**: model value
   2.0029448702957864 = a b-grid node ratio, identical at three different
   parameter vectors; node spacing ≈ 3 bootstrap SEs; Jacobian row exactly
   zero. It must be replaced by a smooth composition moment (standing rule:
   never hard-target grid-discrete medians).

## Untargeted diagnostics at the M4 winner

p90/p50 estate 1.763 (PSID 5.291 raw / 3.448 ratio); 76–84 nonhousing
p50/p75/p90 = 0/0/0 (PSID 2.37/10.34/26.27); top-decile housing share 0.996
(PSID 0.186); **family-size estate gap 4.22 vs 0.10** even with child-blind
bequests (the housing channel alone overshoots the family wealth gradient
~40x); ownership path acceptance False (+36pp climb ages 34–42 vs ACS; 95%+
ownership from the mid-50s); retirement income degenerate.

## Process failures recorded

Six-chain budget starvation; profile grid silently centered on a failed
winner; two mid-run launches in an unattended window against the
present-and-wait rule; acceptance criterion 6 (identification) never encoded
in `acceptance_criteria.csv`; profile README wrong about its task-4 cell;
`CALIBRATION_STATUS.md` never updated post-launch (canonical file denied the
run existed); all M4 code untracked until 2026-07-16 late (committed with this
memo); no Nb=240 verification of the winner; `code/model/README.md`
regeneration commands point at the superseded 12.654 winner.

Also corrected: the "29.7%" loss mis-framing (self-corrected in-session);
Nakajima–Telyukova provenance in `bequest_specification_memo_20260714.tex`
(gamma=0.43 / zeta≈$19,600 appear in NO version of that paper; published
values gamma=20.534, zeta=$7,619/yr in 2000 dollars, Table I Panel B — the
0.25-of-income anchor survives via their implied $30,000 retiree income);
Kvaerner (2023) title in the July-15 note.

## Disposition

M4 is the warm start and baseline for M5, not the paper object. The M5
contract (`docs/model/m5_recalibration_contract_20260716.md`) replaces the
income process and entry condition with directly estimated PSID objects,
swaps the grid-discrete median for a smooth composition share, and settles
theta1 (two-sided profile) and tenure kappa (sweep, then decision with an
ownership-path moment) with predeclared gates.
