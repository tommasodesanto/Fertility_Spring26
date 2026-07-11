# FABLE CROSS-CHECK PROMPT — combined-spec audit adjudication (2026-07-11)

You are an independent adjudicator. A full audit of the intergenerational housing–fertility
model (repo `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`, HEAD `fe092ca99030e54fb256206725e66d9694dc539d`,
DIRTY working tree = the audited object) reached the conclusions below. Your job is to try to
FALSIFY the consequential ones, not to summarize them. Read the cited lines yourself; rerun the
small commands; disagree where the evidence disagrees. Do not modify any file outside
`output/model/full_audit_20260711/crosscheck/` (create it). The model venv is
`code/model/.venv/bin/python`; run from `code/model` with
`env NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=$PWD`.

## Candidates under audit

| Object | Value |
|---|---|
| current_bound_best | loss 14.780020699972585, residual 2.906086971430419e-05, arm global_current, scratch `.../overnight_20260710_wave2_lateral/global_current/task_1/polish/cases.jsonl` line 362 |
| housing_relaxed_best | loss 13.90346531862862, residual 3.140213892938831e-05, arm global_housing, `.../global_housing/task_17/polish/cases.jsonl` line 207 (DIAGNOSTIC bounds) |
| Local extracts | `output/model/combined_recalibration/overnight_20260710_report/{current_bound_best,housing_relaxed_best,summary}.json` |
| Paper's frozen candidate | loss 16.13191908615311, `output/model/intergen_circulation_candidate_20260711/` |

## File hashes (SHA256) — confirm you are reading the audited bits

| File | SHA256 |
| `code/model/intergen_housing_fertility/calibration.py` | `3756da5da88e3123c74572b204835305f2d06cadbee6451cb0ef5f03763afde6` |
| `code/model/intergen_housing_fertility/kernels.py` | `c94cc6bede09c356bd58aff7acf02fc1f6249c1c88288f0f404118637687246e` |
| `code/model/intergen_housing_fertility/local_panel.py` | `2df99277f78f02804561aa395dd848a66c847be3f341cb9a07c63fa6c199676c` |
| `code/model/intergen_housing_fertility/parameters.py` | `d986b6b67d44f62023668ba18e25d0862fffccdc119286f850e257890ef1f829` |
| `code/model/intergen_housing_fertility/production_profile.py` | `922e0b947ceefb3072647ad83d8fae3c2f40a24abbb5916481da702c6126cc3f` |
| `code/model/intergen_housing_fertility/solver.py` | `3ef00d1d1c58e813f1636425ba42654009e28e4f2f5262a3c1ecb24e8a73339e` |
| `code/model/tools/run_intergen_combined_recalibration.py` | `afcd5755350f374324d49ac300e588b94c71a22f75355e2d3f6116ce79583d3c` |
| `tmp/overnight_combined_20260710/reduce_results.py` | `969dec9a550573062a8f3867547600714d599f5fd7fab42ba0e1d1fa28d79ed6` |
| `tmp/overnight_combined_20260710/run_overnight_combined.py` | `aa15e99026a6cdc1c54a85b7a8db30b6b06407245f3af8d42b56070c91651376` |
| `tmp/overnight_combined_20260710/run_wave2_lateral.py` | `62c635ed70dd5ae9a247022705690545345afe55077997084dc7341d5e11b9c4` |
| `output/model/combined_recalibration/overnight_20260710_report/current_bound_best.json` | `77651d72e41344acccc8f068695e0062b7e2663f126d291e1f13845a1adf16f5` |
| `output/model/combined_recalibration/overnight_20260710_report/housing_relaxed_best.json` | `e4f6d4b0b0835f307b854604510d322cb0575d6e90b227cd5ac3a1e6716241c8` |
| `output/model/combined_recalibration/overnight_20260710_report/summary.json` | `d6cb34be7d5dbf958a5ce84381384e7d58d3afa9d75ea75f0e51ad5f49e7c074` |

## Minimal file list to inspect

1. `code/model/intergen_housing_fertility/solver.py` — lines 2305-2312, 3602-3679, 3710-3716, 4614-4616, 4826-4877, 5604-5620
2. `code/model/intergen_housing_fertility/local_panel.py` — lines 510-700 (run_local_polish), 834-894 (run_local_panel_case), 1031-1098 (income), 1151-1167 (selection)
3. `code/model/intergen_housing_fertility/parameters.py` — lines 77-160, 207-340, 408-471
4. `code/model/intergen_housing_fertility/production_profile.py` (entire, 231 lines)
5. `code/model/intergen_housing_fertility/kernels.py` — lines 895-950
6. `tmp/overnight_combined_20260710/run_overnight_combined.py`, `run_wave2_lateral.py`, `reduce_results.py`, `run_array.sh`, `run_wave2.sh`
7. `output/model/combined_recalibration/overnight_20260710_report/summary.json`
8. `output/model/full_audit_20260711/` — LEAD_WORKSHEET.md, repro/*.json, workflow_findings.json, FULL_AUDIT_20260711.md
9. `latex/intergenerational_housing_fertility_paper_draft.tex` — lines 470-530, 1000-1040
10. `code/data/moment_standard_errors/` (SE harness) and `output/model/intergen_circulation_candidate_20260711/README.md`

## Claims to adjudicate (each with its falsification test)

### C1 (FATAL). Uniform fertility on infeasible nodes contaminates ~1/5 of births.
Claim: at nodes where all fertility options carry the −1e10 sentinel, `lf = Vfa/P.kappa_fert; ls, pr = logsumexp(lf, axis=3)` (solver.py:2309-2311) yields pr=(1/3,1/3,1/3); the KFE injects 20.4% of each entering cohort onto such nodes (solver.py:3710-3716, entry_wealth_mode income_ratio_distribution, c_bar_0=1.28, entry income 0.650); j=0 invalid-node mass 0.011986 splits into EXACT thirds by parity; total cohort births ≈0.0583 → ≈21% mechanical (gross upper bound); TFR net of the artifact ≈1.57 vs reported 1.983.
Falsify by: running `output/model/full_audit_20260711/scripts/audit_diag_packet.py` (or the probe scripts referenced in diagnostics/DIAGNOSTICS_REPORT.md) at the current-bound theta and checking (a) fert_probs rows at V≤−1e9 nodes, (b) mass-by-parity equal thirds at j=0, (c) the artifact-birth share arithmetic. Judgment required: is FATAL right for the fertility block given the upper-bound caveat, and does it also invalidate childless_rate-based claims in the paper's 16.13 candidate?

### C2 (MAJOR). Equilibrium-acceptance slack ≈0.27 and Nb=240 drift +1.8.
Claim: at the same theta, max_iter_eq=30 (and tol 2.5e-5) give loss 15.051 vs reported 14.780; Nb=240 gives 16.587 (16.340 tight); own_rate_2534 falls 0.2506→0.1956 at 240.
Falsify by: re-running `scripts/reproduce_candidate.py` with `--max-iter-eq 30`, then `--nb 240` (each <10 min locally) and comparing to `repro/current_eq30.json` / `repro/current_nb240.json`. Judgment: does 0.27 slack + 1.8 drift void the 14.78-vs-13.90 ranking and the "search at 120 without a 240 verification" reporting?

### C3 (MAJOR). Both winners descend from ONE wave-1 seed, probably from a diagnostic arm, clipped onto the production box.
Claim: beta=0.8411423842686492, c_bar_n=0.4559787014242588, theta0=0.1318301350511569 identical to 16 digits across arms; coordinate diffs decompose into pattern steps (ΔH0=2·0.08·9=1.44); h̄₀=1.0 sits exactly at the production lower bound while the housing winner is at 0.921<1 → the shared seed was representable only in the relaxed box and `global_unit_from_theta`+np.clip projected it.
Falsify by: reading summary.json theta values; checking `lp.global_unit_from_theta`/`theta_from_global_unit` clip behavior (local_panel.py:930-963); OR (decisive, needs Torch SSH) reading wave-1 `overnight_launch.json` source_seed_path and wave-2 `lateral_summary.json`. Judgment: does this void the "admissible current-bound best" framing, or is clip-then-362-eval-polish a legitimate box-constrained estimate?

### C4 (MAJOR). H0's [1,10] bound mechanically caps the rooms target.
Claim: supply=H0·(p·user_cost_rate/0.16)^1.75=5.551 at p_eq=0.68976 with H0=9.9997; target 5.780 needs H0≈10.4.
Falsify by: arithmetic from parameters.py:146 (r_bar=0.16), the recomputed user_cost_rate=0.16571, and solver.py:1471; or one solve with H0 forced to 10.4 (diagnostic). Judgment: MAJOR (defeats the stated identification of the new target) vs MINOR (2.1% of loss)?

### C5 (MAJOR). SE-based weights exist, unused; ranking flips under them.
Claim: `code/data/moment_standard_errors/` (July 5-6) reproduces 12/14 targets to <5e-7 and implies inverse-variance weights differing from the hand weights by up to ~500× (relative); under SE weights housing_relaxed_best ranks WORSE than current_bound_best (the old-age wealth-gap and own_family_gap re-weighting dominates).
Falsify by: recomputing both candidates' losses under 1/SE² weights from that folder's outputs against the stored moment vectors (`scripts/recompute_loss.py` pattern). Judgment: does this justify withholding all cross-candidate ranking claims?

### C6 (MAJOR). Not reproducible from committed source; stale-defaults trap family.
Claim: HEAD lacks the run_local_polish kwargs (crash) AND has the pre-flip merge order (income before profile → old 5-point z-grid would silently win); parameters.py defaults are the OLD economy (q=0.169859, δ=0.077632, old income process, normalize_bequest_utility=False); production_profile_overrides() carries none of the combined spec; eight tools + CLI replay records under the wrong economy.
Falsify by: `git stash` is FORBIDDEN — instead `git show HEAD:code/model/intergen_housing_fertility/local_panel.py | grep -n "income" | head` and compare merge order; check parameters.py:80-82; run any tool (e.g. build_intergen_mechanics_packet.py --help) and inspect which overrides it can pass. Judgment: MAJOR-process vs blocking-for-circulation?

### C7 (MAJOR). Paper quantitative content = frozen 16.13 candidate; policy magnitudes fixed-theta with 2.4× closure sensitivity.
Claim: all tables/figures match `output/model/intergen_circulation_candidate_20260711` digit-for-digit (figures byte-identical); +3.4%/+4.9% births numbers come from the phase9b battery at that theta; closure variants in the same table give price effects differing up to 2.4×; four estimates sit on undisclosed bounds.
Falsify by: spot-checking three table numbers against the frozen packet's target_fit.csv, `cmp` on the two repaired figures, and reading the battery README + its closure columns. Judgment: is "withhold policy magnitudes" right, or does a disclosed fixed-theta diagnostic framing suffice for a workshop draft?

### C8 (MAJOR, provenance of a target). Rooms target is a 42-metro MMS object from an uncommitted builder edit, self-flagged "candidate".
Falsify by: `git diff HEAD -- code/data/mms_center_periphery/build_intergen_one_market_housing_targets.R` and reading the output CSV's status field + sample filters. Judgment: acceptable target with relabeled description, or re-audit before next estimation?

### C9 (MAJOR). Income process dispersion is a re-parameterization of the old ad hoc grid.
Claim: ρ_a=0.85^(1/4) exactly; σ chosen to hit stationary log-sd 0.23102 = the old grid's value exactly; half the Sommer–Sullivan dispersion (0.459) that the repo's own comparison tool cites as the anchor.
Falsify by: arithmetic (0.9601845894^4; the weighted log-variance of [0.6,0.8,1,1.2,1.4] at [.1,.2,.4,.2,.1]); read `code/model/tools/compare_intergen_income_processes.py` docstring. Judgment: is this defensible as disclosed continuity, or must the paper's income-calibration claim change / the process be re-estimated?

### C10 (disputed — lead overrode one auditor). Owner services in the paper.
One auditor claimed the draft displays χ·H−h̄ (finance F1); the lead read tex:489 `\chi_O[h-\bar h(n,s)]` and withdrew the finding; a verifier A/B showed the alternative form would shift ownership moments 9-12pp (why the display matters).
Falsify by: reading latex/intergenerational_housing_fertility_paper_draft.tex lines 485-495 in the DIRTY tree. Judgment: confirm the withdrawal or reinstate.

### C11 (MAJOR). Grid-discrete hard-targeted median.
Claim: old_nonhousing_wealth_to_income_median_6575 is an uninterpolated b-grid value; one node move shifts the objective 0.54-1.24, exceeding the 0.88 candidate gap; violates the standing no-hard-median rule.
Falsify by: locating the median computation in solver.py (weighted_median on grid values) and checking step size at the candidate. Judgment: replace-before-next-estimation vs disclose?

## Open questions where independent judgment is genuinely required

1. Severity calibration of C1: FATAL-for-fertility-block vs FATAL-for-the-whole-calibration (the fertility parameters co-move with the housing block through h̄_n, c̄_n).
2. Is the honest headline number 16.3-16.6 (Nb=240 tight/production) — and should ALL paper numbers be re-solved at 240 despite the 120-search protocol?
3. Does the single-basin finding (C3) demand a genuinely global re-search after the fixes, or is seed-descent acceptable given the 15-moment system's noise floor?
4. Which of the five persistent misses (old own, old liquid, young wealth, room gap, renter rooms) are worth model surgery (exit margin, transfers) vs re-weighting under SEs (C5) BEFORE the next search is bought?
5. The blocked cluster items (§12 of FULL_AUDIT_20260711.md): rank which of the six matter enough to gate circulation.

## Ground rules

Cite file:line for every verdict. Quote what you read. Where you disagree, say what evidence
would settle it. The lead's own worksheets (LEAD_WORKSHEET.md, ECONOMIC_ASSESSMENT.md) are
evidence, not authority.
