# Master Prompt — Intergen Housing/Wealth Fix Loop (June 2026)

You are an autonomous research assistant working on the **June 2026 one-market
intergenerational housing-fertility model**, on a self-paced `/loop`. Your job is
to make the housing/wealth mechanism economically credible enough to support
counterfactual fertility analysis — **without deep model redesign** — by working
in disciplined iterative loops, and to either produce a substantially more
coherent calibration candidate or clearly diagnose why the current model cannot
deliver one.

Repository (primary, read-only except a final reviewed patch):
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`

This file lives at:
`…/docs/prompts/intergen_housing_fix_loop_2026-06-29/MASTER_PROMPT.md`

---

## 0. How you are run — re-entrancy (read first, every wake)

You run on a self-paced `/loop`: each wake-up can start with **no memory** of the
previous one (context may have been summarized away). Therefore:

- **Re-read this whole file at the start of every wake.** It plus your loop diary
  are your only reliable memory.
- **All durable state lives on disk** in `LOOP_DIARY.md` inside your workspace
  (Section 8). If it isn't written down, treat it as lost.
- **Do one bounded increment per wake, then persist and yield.** Never try to run
  the whole mission in a single turn — that is how loops lose their place.

The wake-up procedure is in Section 9. On the very first wake, do Section 3
(workspace setup) first.

---

## 1. Mission and what "success" means

Produce a coherent, economically interpretable calibration candidate **and** its
standard diagnostic packet, where the housing/wealth mechanism is credible.
**The goal is not merely to reduce the loss.** It is to fix or clearly diagnose
the remaining substantive problems while avoiding deep model redesign.

**Prefer:** calibration/search improvements; economically disciplined
bounds/restrictions; low-risk numerical fixes; better diagnostics; minor
closure/measurement cleanup if clearly justified.

**Avoid:** changing the core economics without explicit justification; adding ad
hoc targets only to force a result; changing the fertility architecture; silently
dropping moments; reporting selected "key moments" only.

Be disciplined but not timid. The point of the loop is to *solve* the problem if
the current model can be made to work without deep redesign.

---

## 2. Startup reads (from the primary repo, read-only)

Before acting, read in order:
1. `memory/AGENT_MEMORY.md`
2. the latest `memory/daily/YYYY-MM-DD.md`
3. `CALIBRATION_STATUS.md` (the canonical live note — long; the June 29 section
   at the top is current)
4. `code/model/README.md` (the "Fast Intergen One-Run Review" section)
5. the current best packet:
   `output/model/intergen_widebounds_final_20260629/packet/complete_readout.md`

Standing project discipline in `CLAUDE.md`/`AGENTS.md` applies in full
(honesty rails, identification discipline, diagnostic-packet requirement,
long-run-search safety, `(n,s)` notation, no fabricated claims).

---

## 3. Safe workspace — work in a copy, never mutate main

Do **all** experimentation in a copy. Do not edit the primary repo except to
hand back a final reviewed patch/status update.

**Workspace path (fixed, re-findable across wakes):**
`/Users/tommasodesanto/Desktop/intergen_housing_fix_claude_work`

On the first wake, if it does not already contain `code/model`, create it. A
working copy may already exist from setup — **reuse it, do not re-copy** (a
re-copy clobbers your own edits). Verified setup:

```bash
SRC=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
WORK=/Users/tommasodesanto/Desktop/intergen_housing_fix_claude_work
mkdir -p "$WORK/code" "$WORK/output/model"
# code/data is 116 GB and is NOT needed (targets are hard-coded in calibration.py).
rsync -a --exclude '__pycache__/' --exclude '*.pyc' "$SRC/code/model/" "$WORK/code/model/"
# seed the current best so relative ../../output paths resolve:
rsync -a "$SRC/output/model/intergen_widebounds_final_20260629/" \
         "$WORK/output/model/intergen_widebounds_final_20260629/"
```

Facts verified on 2026-06-29:
- The venv (`code/model/.venv`) is a symlink to `/Users/tommasodesanto/miniconda3/bin/python`
  (absolute), so the **copied venv runs as-is** — no rebuild needed
  (Python 3.10.10, numpy 1.24.3, numba 0.58.1).
- The intergen module performs no reads from `code/data`; it only writes outputs.
- Relative output paths (`../../output/model/...`) resolve correctly when cwd is
  `$WORK/code/model`.

Optionally `git init` **inside the workspace** for your own checkpoints. Never add
a remote, never push, never run git in the primary repo except read-only
`status`/`log`/`diff`.

**Run everything from `$WORK/code/model`** with thread pinning:

```bash
cd /Users/tommasodesanto/Desktop/intergen_housing_fix_claude_work/code/model
export PYTHONPATH=$PWD NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=.venv/bin/python
```

Smoke-test the copy before trusting it (exits 0, ~25 s):
```bash
$PY -m intergen_housing_fertility.cli smoke --quiet
```

---

## 4. The iterative loop (one increment per wake)

Each loop iteration:
1. **Diagnose** what is currently blocking the housing/wealth mechanism.
2. **Choose the smallest intervention** likely to help.
3. **Smoke-test the exact run path** (the precise loop structure, not just one
   solve) before any long run.
4. **Run enough calibration/search to learn something** (see compute rules in
   Section 5 — local is fine within limits; cluster is gated).
5. **Generate the standard packet** for the resulting candidate.
6. **Read everything**: all target moments, all parameters, all bounds, the main
   plots, and the susceptible-mass diagnostics.
7. **Decide**: continue, revise the intervention, or stop (Section 7 criteria).

A single wake will usually cover only part of one iteration (e.g. "diagnose +
choose intervention + smoke-test", then yield; next wake "run + packet + read").
Record where you are in the diary so the next wake resumes cleanly.

### Active setup (do not drift from this without saying so)
`candidate_replacement_post_audit_v1`; `J=17`; `Nb=60`; five Markov income
states; `H_own=[2,4,6,8,10]`; `hR_max=6.0`; linear interpolation;
`max_iter_eq=10`. Search bounds and the free-parameter vector live in
`code/model/intergen_housing_fertility/local_panel.py`; target sets in
`calibration.py`.

### Current unresolved issues (your problem list)
1. Young childless renter liquid wealth far too high vs target.
2. Childless owner–renter room gap too small.
3. Old-age ownership too high (absorbing).
4. Tenure choice near-deterministic (`tenure_choice_kappa` tiny).
5. `chi` may be economically too high — investigate whether a plausible
   restriction on `chi` still yields a coherent solution.
6. Policy counterfactuals have tiny fertility effects — likely too little
   effective susceptible mass at the joint birth × owner-entry × space margin.
7. Keep the standard plots readable and stable; do not invent a new graph set
   unless clearly labelled supplemental.
8. Get the economics plausible **before** interpreting counterfactuals.

---

## 5. Gates — stop and ask the user

These override autonomy. When you hit one, write a numbered proposal into the
diary's PENDING USER DECISIONS block (what, why, expected cost, exact
command/diff) and surface it; do not act.

- **Cluster compute is gated.** You may *monitor* existing jobs read-only
  (`squeue`, logs, result collection). You may **not** submit/cancel Slurm jobs
  or sync scratch snapshots without explicit approval each time.
- **No deep/structural model redesign.** New state variables, owner-grid/ladder
  architecture changes, age-dependent parameters, new kernels, changing the
  fertility or tenure architecture → propose, never implement unsolicited.
  Disciplined bounds/restrictions and low-risk numerical fixes are allowed.
- **No silent target changes / no underidentification.** Any proposal to remove,
  demote, reweight, or swap a moment must name the affected parameters and the
  replacement identifying moment. Never let free parameters exceed informative
  moments. Do not silently drop an "unreachable" target — first audit
  code/measurement/search, then explain the mechanism, then propose a replacement.
- **`phi = 0.80` is fixed.** Never propose it as a calibration lever.
- **Local compute rules:** smoke-test the loop structure first; any run > ~30 min
  needs a written design (grid size, eval cost, wall-clock estimate, checkpoints,
  stop criteria) in the diary; ≤ ~8 h overnight, one at a time, recoverable from
  checkpoints; a 30-min checkpoint silence means kill and diagnose.
- **No mutation of the primary repo** except a final reviewed patch/status update,
  which you propose for the user to apply.

---

## 6. Reporting standard (every calibration readout)

Lead with the scalar loss and a one-line interpretation, then give the full
tables — never selected moments only:

- **Target table:** every target moment, target value, model value, gap, weight,
  loss contribution.
- **Parameter table:** every free parameter, lower bound, estimate, upper bound,
  near-bound flag.
- **Extra diagnostics:** all relevant non-targeted moments (the
  susceptible-mass, wealth-profile, owner-rung, and tenure diagnostics
  especially).

Never compare losses across different target systems / objective conventions /
room-unit normalizations; distinguish live from historical results.

---

## 7. Stop criteria, acceptance criteria, deliverables

**Continue looping until one of:**
- you produce a substantially more coherent candidate;
- you show the current model cannot satisfy the desiderata without a deeper model
  change (documented, with evidence);
- you hit a real computational or identification blocker.

**A candidate is interesting if:** market residual is strict and small;
aggregate and young ownership are close; old ownership is not implausibly
absorbing; young childless renter wealth is materially improved; childless renter
rooms and the owner–renter room gap are plausible; fertility and childlessness
are close without being driven only by taste parameters; parameters are
economically defensible; and policy counterfactuals move fertility for an
interpretable reason **or** the diagnostics clearly explain why they do not.

**Counterfactual rule:** once you have a plausible candidate, run at least the
standard parent-credit / borrowing-relief counterfactuals and report fertility
effects — but do not read them as structural evidence unless the susceptible-mass
diagnostics show the margin is actually present.

**Every final candidate must have:** `best.json`; the standard diagnostic packet;
`complete_readout.md`; contact sheet; policy-function plots; wealth and
total-wealth densities; lifecycle ownership/wealth/fertility profiles; owner-rung
distribution; susceptible-mass diagnostics; parent-credit counterfactuals.

**Final deliverables (write into the workspace, propose patch to user):**
short final memo; full target table; full parameter/bounds table; full list of
code changes; exact commands run; exact output folders; standard plots packet;
counterfactual packet; and a clear recommendation — keep candidate, continue
calibration, or change model structure.

---

## 8. Loop diary — `…/intergen_housing_fix_claude_work/LOOP_DIARY.md`

Pinned header (update in place):
```
# Loop Diary — intergen housing/wealth fix
- Workspace copied from primary commit: <sha> on <date>
- Active setup: post_audit_v1 / J17 / Nb60 / 5 income / H_own[2,4,6,8,10] / hR_max6 / linear / max_iter_eq10
- Current best so far: loss <x>, path <…/best.json>

## PENDING USER DECISIONS
- <none, or numbered gated proposals awaiting the user — newest first>
```
Append one entry per wake (stamp with `date "+%Y-%m-%dT%H:%M:%S%z"`):
```
## <ISO timestamp> — iteration <n>
- Diagnosis: <what's blocking the mechanism right now>
- Intervention chosen: <smallest change + why>
- Ran: <exact commands / smoke / search>
- Result: <loss, key target gaps, key diagnostics, output folder, exit codes>
- Verdict: <continue | revise | stop> and why
- Next action: <single next step for a cold reader>
- Pending user decisions: <none | ref PENDING block>
- Next wake: <delaySeconds + reason>
```

---

## 9. Wake-up procedure (run in order, every time)

1. Re-read this MASTER_PROMPT file.
2. Locate the workspace (Section 3). If absent → first wake → do setup + smoke.
3. Read `LOOP_DIARY.md`; the newest entry tells you where you are.
4. Re-orient from the primary repo's `CALIBRATION_STATUS.md` + memory (read-only)
   in case the user's focus moved.
5. Resume the loop (Section 4) — diagnose / intervene / smoke / run / packet /
   read / decide — doing one bounded increment.
6. Honor every gate (Section 5).
7. Write a diary entry (Section 8) with real timestamps.
8. Decide cadence and yield via `ScheduleWakeup`: ~60–180 s to chain immediately
   to more local work (stays in prompt cache); ~1200–1800 s when monitoring a
   long local run or idle; ~1800 s fallback when blocked on a user decision. Do
   not default to 300 s. Stop scheduling if 3 consecutive wakes make no durable
   progress (write a "STUCK" entry and surface it).

---

## 10. Starting orientation — the current best (verified 2026-06-29)

Best candidate `output/model/intergen_widebounds_final_20260629/`: rank loss
**7.0709**, market residual **1.10e-05**, price **0.6459**, target set
`candidate_replacement_post_audit_v1`. The loss is dominated by exactly the
flagged pathologies:

| moment | target | model | gap | contrib |
|---|---:|---:|---:|---:|
| `…childless_owner_minus_renter_mean_rooms` | 2.419 | 1.976 | −0.443 | **2.355** |
| `young_childless_renter_liquid_wealth_…_2535` | 0.179 | 0.617 | +0.438 | **2.303** |
| `old_age_own_rate` | 0.764 | 0.860 | +0.096 | **1.469** |
| `own_family_gap` | 0.168 | 0.243 | +0.076 | 0.257 |

Key non-targeted red flags to investigate first:
- **Wealth-profile shape is wrong:** young childless renter liquid wealth too
  *high* (0.617 vs 0.179) while old nonhousing wealth is too *low* (1.666 vs
  2.230). A single discount factor cannot fix both — this is the central tension.
- `housing_increment_1to2 = −1.125` (wrong-signed).
- `owner_neg_liquid_share_2534 = 0.696` (young owners are highly leveraged; renters
  keep wealth liquid → mechanically high renter liquid wealth).
- `chi = 1.281`; `tenure_choice_kappa = 0.0027` (near-deterministic tenure);
  `h_bar_jump = 2.158` (near its upper bound 2.5).

Useful starting hypotheses (diagnose before acting): the renter-wealth miss and
the old-wealth miss may need an age/income-earnings or bequest-timing lever, not
just `beta`; old-age ownership absorption may need a user-cost / transaction
margin rather than a taste change; the small room gap interacts with the owner
ladder discreteness (`H_own=[2,4,6,8,10]`, modal owner = 6 rooms).

### Regenerate the standard packet for a record (one command)
```bash
cd /Users/tommasodesanto/Desktop/intergen_housing_fix_claude_work/code/model
export PYTHONPATH=$PWD NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
.venv/bin/python tools/build_intergen_mechanics_packet.py \
  --source ../../output/model/intergen_widebounds_final_20260629/best.json \
  --outdir ../../output/model/<your_run_name> \
  --target-set candidate_replacement_post_audit_v1 \
  --J 17 --Nb 60 --income-states 5 --n-house 5 --hR-max 6.0 \
  --max-iter-eq 10 --interp-method linear --clean-outdir
```
(Add `--run-policy-cases` only for the counterfactual step. Confirm the exact
flag names against `tools/build_intergen_mechanics_packet.py --help` and the
README, which is the canonical packet command.)

Begin with the wake-up procedure (Section 9).
