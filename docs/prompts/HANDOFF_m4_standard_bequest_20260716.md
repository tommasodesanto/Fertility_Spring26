# HANDOFF: M4 standard bequest recalibration (2026-07-16)

Paste this whole file into a fresh chat. It is self-contained. Goal: finish
wiring the M4 arm, verify it, run it on Torch, fill the advisor note, send.

> **SUPERSEDED CONTRACT NOTICE (user clarification, later July 16).** The user
> clarified before launch that M4 must estimate both `theta0` and `theta1`
> internally; `theta_n=0` alone is external. The live contract is therefore 14
> moments / 13 free parameters, with `theta1=0.25` used only as one dispersed
> optimizer start. The corrected implementation, bounds, starts, run budget,
> and identification gate are recorded in
> `output/model/intergen_standard_bequest_recalibration_20260716/README.md` and
> the top section of `CALIBRATION_STATUS.md`. Any 12-parameter or externally
> fixed-`theta1` instruction below is historical and must not be followed.

Repo root: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`

---

## 0. FIRST ACTION — establish the unknown

A subagent was mid-way through wiring the M4 arm when the session hit its
usage limit. Its last words were "Now verify: py_compile, bash -n, and the
test suite" — meaning it had probably **written edits it never verified**.
Nobody has looked at those edits. Run this first:

```bash
cd /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
git status -s code/model/intergen_housing_fertility/ code/model/tools/ code/cluster/
rg -n "M4|median_composition_v1|NONHOUSING_MEDIAN" \
  code/model/intergen_housing_fertility/calibration.py \
  code/model/tools/run_intergen_bequest_exit_chain.py
ls code/cluster/ | grep -i standard_bequest
```

Then decide: keep and verify the partial edits, or `git checkout` those two
files and redo the wiring cleanly (spec in §3). **The two model files are
objective-critical — the lead must read the diff line by line against §3
before anything runs.** Note both files also carry pre-existing uncommitted
changes from earlier sessions; do not blanket-revert.

---

## 1. What was decided (with Tommaso, this session)

Return the bequest block to the standard package used by the structural
savings / housing-retirement literatures. No new mechanism, no revolution.

**Specification (M4):**

$$B(W)=\theta_0\,\frac{(\theta_1+W)^{1-\sigma}-\theta_1^{1-\sigma}}{1-\sigma},
\qquad W=\max\{b+pH,0\},\quad \sigma=2$$

- De Nardi (2004) child-blind luxury warm glow. In code: existing
  `linear_child_scale` spec with `theta_n = 0`, `normalize_bequest_utility =
  True`. **No new code path.**
- `theta0` — the ONLY estimated bequest parameter.
- `theta1 = 0.25` — external (Nakajima–Telyukova anchor, see §2).
- `theta_n = 0` — external (Hurd 1989 null; Kopczuk–Lupton 2007; DFJM 2025;
  our own PSID gap 0.101 with SE 0.563).

**Targets: 14 moments, 12 free parameters** (11 clean-frontier + `theta0`).
The 12 established moments unchanged; the two late-life wealth targets are
LEVELS (DFJ/NT style):

| Moment | Target | Weight (1/SE²) |
|---|---:|---:|
| `old_total_estate_wealth_to_annual_income_median_7684` | 6.50131577436537 | 18.585767349158665 |
| `old_nonhousing_wealth_to_income_median_6575` | **1.90821154211154** | **83.74916751466371** |

Demoted to reported-but-untargeted diagnostics: estate p90/p50, the
2+-minus-1-child gap, the housing/nonhousing decomposition, ownership path.

**Why:** p90/p50 was 84% of the M3 loss (139.93 of 166.65) and is
structurally unreachable — it demands a financial/business wealth tail the
model has no state for. The family gap is statistically zero (t≈0.18) and a
marital-composition cancellation (−0.886 married / +0.884 nonmarried) in a
model without marriage. Chasing them wrecked the established block: summed
loss on the 12 shared moments 3.84 (M1) → 16.92 (M3), and the young renter
liquid-wealth moment flipped sign (+0.19 → −0.18).

Full diagnosis, verification table, and the identification argument:
`docs/model/intergen_bequest_balance_sheet_fable_review_20260716.md`
(committed `c65dcb4`, amended by the literature-anchor commit, pushed).

---

## 2. Verified inputs — DO NOT re-derive, these cost real time

### R bootstrap (COMPLETE; lead verified the diff and CSV firsthand)

`old_nonhousing_wealth_to_income_median_6575 = 1.90821154211154`,
bootstrap SE `0.109272216032637`, weight `1/SE² = 83.74916751466371`.
499 reps, person-cluster, seed 20260715, same draws as the other three.

- Modified (UNCOMMITTED): `code/data/psid_followup_mar2026/audit_intergen_bequest_family_size_targets.R`
  — new moment computed inside `bequest_target_values()` (deterministic, no
  new RNG calls), same filters as the 65–75 estate targets.
- Regenerated: `.../output/intergen_bequest_family_size_audit/bequest_calibration_targets.csv`
  (now 4 rows) and the 4×4 covariance CSV.
- Regression gate PASSED: the three legacy rows are bit-identical
  (6.50131577436537/0.2319582116443; 3.44811075444552/0.132477154936326;
  0.101108088567873/0.563010076238571).

**KEY FINDING:** the legacy stored `2.23046078` (`calibration.py:106`) came
from the OLD all-person sample. The corrected reference-person sample gives
**1.908** (−14%). The M1 winner's model value was **2.0029 — within one
bootstrap SE of the corrected target.** The composition target should not
fight the established block.

### Literature (COMPLETE; verified from primary sources, not abstracts)

- **Nakajima–Telyukova (2017), JF 72(2), published, Table I Panel B:**
  `v(a) = γ(a+ζ)^{1-σ}/(1-σ)`, **γ = 20.534, ζ = $7,619/year in 2000
  dollars**, σ = 2.006.
  **The project's recalled γ=0.43, ζ≈$19,600 are WRONG — they appear in NO
  version of that paper.** (`bequest_specification_memo_20260714.tex` ~line
  110 must be corrected.) Provenance of the version history: 2013 draft
  γ=3.92; WP 14-27 γ=6.576, ζ=15,354 *biennial*; published γ=20.534,
  ζ=7,619 *annual*.
  **The 0.25 anchor SURVIVES**: the paper's own transfer figures ("$252 =
  0.84% of median annual after-tax income") imply retired-homeowner median
  after-tax income of $30,000, so ζ/income = 7,619/30,000 ≈ **0.25**.
- **De Nardi (2004), ReStud, eq. (8) p.750, Table 4 p.752:**
  `φ(b) = φ₁(1+b/φ₂)^{1-σ}`, φ₁ = −9.5, φ₂ = 11.6 (recall correct).
  **BUT φ₂ is in model units and the paper states NO dollar mapping** — it
  CANNOT anchor a conversion. Any "De Nardi-scale sensitivity" framing is
  unsupported; drop it.
- **De Nardi–French–Jones (2010), JPE, eq. (3) p.44, Table 3 p.57:**
  `φ(e) = θ(e+k)^{1-ν}/(1-ν)`, θ = 2,360, **k = $273,000 (1998 dollars)**,
  ν = 3.84; both bequest params statistically insignificant (p.62).
  k/top-quintile annuity income ≈ 13.7. **This is the right high-scale
  sensitivity anchor** (replaces the unsupported De Nardi one).
- **Kværner (2023):** correct title is "How Large Are Bequest Motives?
  Estimates Based on Health Shocks", RFS 36(8), 3382–3422. The July-15
  note's bibitem title is wrong.
- PDFs were downloaded to a session scratchpad that is likely gone; re-fetch
  from De Nardi's NBER page / Wiley if you need to re-check.

---

## 3. M4 wiring spec (what the dead agent was doing)

**Files:** `code/model/intergen_housing_fertility/calibration.py`,
`code/model/tools/run_intergen_bequest_exit_chain.py`, the M3 collector
(find via `code/cluster/submit_intergen_internal_bequest_recalibration.sh`),
new cluster scripts, tests. Do NOT touch `solver.py`, `parameters.py`,
`kernels.py`, `production_profile.py`.

1. **calibration.py** — new target set
   `candidate_replacement_bequest_median_composition_v1` = current
   `candidate_replacement_bequest_internal_v1` (defined `:206-220`) MINUS
   `old_total_estate_wealth_to_annual_income_p90_p50_7684` MINUS
   `old_2plus_minus_1_total_estate_wealth_to_annual_income_median_gap_6575`
   PLUS `old_nonhousing_wealth_to_income_median_6575` (moment already exists,
   `solver.py:5374`). = 13 base moments; driver appends the rooms moment → 14.
   Register in `TARGET_SETS`. Constants: `NONHOUSING_MEDIAN_TARGET =
   1.90821154211154`, `NONHOUSING_MEDIAN_WEIGHT = 83.74916751466371`
   (these are FINAL — no placeholders). Estate median/weight unchanged.
   Provenance entry mirroring `:299-314`.
2. **run_intergen_bequest_exit_chain.py** — arm `M4`, contract IDENTICAL to
   M2: `active = BASE_DOMAIN + THETA0_DOMAIN`; fixed `theta_n=0.0`,
   `tenure_choice_kappa=0.0`, `theta1=float(args.theta1)` (default 0.25).
   Mechanism overrides as M0/M1/M3 (`linear_child_scale`,
   `normalize_bequest_utility=True`, `owner_ltv_taper=False`),
   `use_age_survival=True` — add M4 to every `{M1,M2,M3}` gate including
   `survival_schedule`. `target_system()` assertion expects **14** for M4
   (15 otherwise). Nested seeds, reusing M2's injection mechanism:
   (a) strict M1 winner from
   `output/model/intergen_mortality_recalibration_20260715/report/results.json`
   at `theta0=0` (exact nested submodel); (b) M3 winner's 11 shared
   coordinates from
   `output/model/intergen_internal_bequest_recalibration_20260715/report/results.json`
   with `theta0=0.31024680236397595`. **Collector gate: FAIL if the free-θ₀
   winner's tight loss exceeds the θ₀=0 nested seed's tight loss** (this is
   the July-15 A3 failure that must not recur).
3. **Cluster scripts** — clone the M3 submit/smoke/collector into
   `submit_intergen_standard_bequest_recalibration*.sh`: `--arm M4`, 6
   chains, 75 min or 1000 evals per chain, search `(max_iter_eq,tol_eq) =
   (10,1e-4)`, two tight winner solves at `(40,2.5e-5)`, outdir
   `output/model/intergen_standard_bequest_recalibration_20260716/`, account
   `torch_pr_570_general`.
4. **Tests** — M4 set has 13 base entries, contains the nonhousing + estate
   medians, excludes p90_p50 and the family gap; arm M4 has 12 active params.
   Run: `cd code/model && NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1
   MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=$PWD .venv/bin/python
   -m pytest intergen_housing_fertility/tests -q` (43 tests passed as of
   2026-07-14).

---

## 4. Run sequence

1. Lead verifies the wiring diff line-by-line against §3.
2. Tests pass → local exact-loop smoke (mirror how the M3 smoke was invoked
   locally; keep it under ~10 min).
3. **Sync to Torch — the remote is NOT a git clone**; copy changed files
   individually to `/scratch/td2248/projects/Fertility_Spring26/...`
   (`ssh torch` works; queue was near-empty). `code/cluster/torch.sh` has no
   family for this battery — the M3 pattern submits the script directly over
   ssh.
4. Torch smoke → 6 production chains → dependent strict collector.
5. **Pass criteria (ex ante, from the review memo §6):**
   1. Strict, exactly repeated tight winner (bit-identical two solves).
   2. Estate median within 1 bootstrap SE (|gap| ≤ 0.232).
   3. Nonhousing median within 2 SEs of 1.908 (≤ 0.219) and strictly positive.
   4. Summed loss on the 12 established moments ≤ 1.15 × M1's 3.84 = **4.42**.
   5. Free-θ₀ winner weakly below the injected θ₀=0 seed.
   - If (2) and (3) cannot hold jointly: STOP. That is real information —
     escalate the mechanism decision (leading candidate: out-of-pocket
     medical/LTC expense risk after retirement, DFJ-style, calibrated
     externally from HRS/MEPS so it consumes no moments).
6. Fill the `[PENDING]` blocks in the advisor note, recompile, send.

---

## 5. Advisor note state

`docs/model/bequest_standard_recalibration_note_20260716.tex` — UNCOMMITTED,
compiles clean (`pdflatex` ×2, BUILD_OK, PDF exists). Structure: executive
conclusion, diagnosis (PSID/model decomposition table), the standard
specification with the corrected NT provenance, the child-blind defense with
citations, the M4 design, then `\pending{...}` markers for: the 14-moment
target-fit table, the 12-parameter table, untargeted diagnostics, and the
nested-zero verdict. **It was edited after the last successful build (NT
provenance paragraph + nonhousing target 1.908) and has NOT been recompiled.**

**The referee answer to "your paper has children — how can bequests be
child-blind?"** (drafted, in the note, §"Why a child-blind motive"): children
direct resources through HOUSING while alive — they raise required space,
parents buy bigger houses, transaction costs make them sticky into old age,
and housing IS the estate. So parents leave bigger bequests with no
child-directed taste for posthumous giving. Hurd (1989): the elderly with
children decumulate no differently. Kopczuk–Lupton (2007): children affect
motive *incidence*, not scale. DFJM (2025): children not even a conditioning
variable. Kværner (2023): models the motive child-blind, puts children in
inter vivos transfers (McGarry 1999). Our PSID gap: 0.101, SE 0.563. Tommaso
wants direct page/table references here rather than paraphrase.
**Bonus if M4 behaves:** with θ_n = 0 the model should still generate a small
positive family-size estate gap through housing alone — matching the data as
an untargeted prediction.

---

## 6. Corrections owed elsewhere (not yet done)

- `docs/model/bequest_specification_memo_20260714.tex` ~line 110: γ=0.43 /
  ζ≈$19,600 → γ=20.534 / ζ=$7,619 (2000$, published JF Table I Panel B). The
  ζ/income ≈ 0.25 claim is right but for the wrong reason.
- `docs/model/bequest_exit_note_20260715.tex`: Kværner bibitem title wrong;
  "age-82 forces terminal liquidation" is FALSE for production/M3 (verified:
  terminal continuation is `Vnr = Vbq`, tenure-indexed, no forced sale —
  `solver.py:2189-2190`). It described only the rejected owner-LTV taper arms.
- `CALIBRATION_STATUS.md` July-15 section: same forced-liquidation error. The
  file has a July-16 pointer section added at the top (uncommitted) plus ~464
  lines of pre-existing uncommitted changes from earlier sessions — do not
  blanket-revert.
- `parameters.py:67`: `normalize_bequest_utility` defaults to `False` while
  every M-arm overrides it to `True`. Latent trap; flip the default or assert.

---

## 7. Environment gotchas

- **Codex CLI is NOT installed** — `npm install -g @openai/codex`, then
  `/codex:setup`. Both coding tasks bounced off this and had to be rerouted
  to built-in subagents. Fix this first if you want the standard routing.
- Torch: authenticated, `ssh torch` OK, account `torch_pr_570_general`,
  remote dir `/scratch/td2248/projects/Fertility_Spring26/` (not a git repo).
- Local venv: `code/model/.venv` — imports OK (numba 0.58.1, numpy 1.24.3).
- `memory/` is a symlink into the detached nightly store — files there are
  NOT committable. `memory/daily/2026-07-16.md` was written for context.

---

## 8. Honest status

Not delivered: the M4 run and therefore the advisor-ready note. The session
burned its budget on the six-agent audit and lost the wiring agent to the
usage limit mid-edit.

Delivered and durable: the diagnosis (why the last week failed, with
line-level evidence), the corrected literature anchors (one of which
overturns a wrong number that has been sitting in the July-14 memo), the
corrected nonhousing target with its bootstrap SE, and a note that needs only
run numbers. Remaining critical path is roughly: verify/redo wiring (~30 min)
→ tests + smoke (~20 min) → Torch smoke + 6 chains (~90 min wall) →
collector + fill note (~30 min).
