# Master Prompt — JMP Build Session (theory + quantitative model)

You are working in:

```text
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
```

## Mission

Build a **top-level, referee-ready job-market paper** on intergenerational
housing misallocation and fertility — theory plus quantitative model. Not a
cleanup of the existing draft: treat the repository as raw material and act
as a fearless coauthor told to swing for the fences. The data layer is taken
as given (use the recorded empirical targets; flag every provisional number);
everything else — framing, theorems, model structure, calibration design,
counterfactuals, welfare — is yours to push as far as it will truly go.

## Ambition mandate

1. **You own the framing.** Decide what the sharpest version of this paper is
   (misallocation-and-sufficient-statistics paper? policy-instrument paper?
   mechanism paper anchored on the wedge taxonomy?) and write a one-page
   framing memo before the theory phase. Choose the strongest framing the
   evidence can carry, not the safest.
2. **Explicit structural authorization.** The standing project rule "never
   implement structural model changes unless the user explicitly asks" is
   hereby satisfied: the user explicitly authorizes structural changes
   **inside your worktree** — new state variables, inheritance kernels, owner
   ladders, age-dependent parameters, survival risk, welfare modules, grid
   redesigns. Every change must land in `IMPLEMENTATION_STATUS.md` (worktree
   copy) the same edit, and in the paper's model section. Exceptions that
   stay fixed: `phi = 0.80`; the one-shot completed-fertility convention may
   be questioned in a memo but changed only as a stretch goal with its own
   gate.
3. **Prefer the strongest true result.** Where you can prove more, prove
   more. Where the quantitative model can deliver a welfare number instead of
   a moment comparison, deliver the welfare number. Honesty rails stay
   absolute: no fabricated citations, no unverified claims, targets always
   shown next to moments, every limitation stated in the text — ambition
   raises what you attempt, never what you assert beyond evidence.

## Research escalation menu

Work is tiered. CORE is mandatory for the draft; TARGET should be attempted
seriously; STRETCH is attempted only after its prerequisite gate passes.
Log attempts and outcomes (including failures) in STATUS.md — a documented
dead end is a deliverable.

**Theory (Section 2 + appendix):**

- CORE: the corrected compact model (net construction + aggregation lemma;
  old purchase margin; the full-efficiency theorem including the entry
  margin via the Θᵢ = 0 argument; the young–old wedge identity; the
  wedge–instrument taxonomy as a proposition: collateral wedges are
  removable with transfers, income-based PTI wedges require contract-space
  instruments; the bequest-tax comparative static with its three signed
  terms). Seed: `theory_clean_room_note.tex` — but you write paper-grade
  prose and complete proofs, not a patch.
- CORE: an **exact aggregate decomposition proposition**: equilibrium
  fertility loss decomposed into the ζ^DP, ζ^PTI, and L^p components plus the
  entry margin — the paper's sufficient-statistics centerpiece, stated so the
  quantitative model can evaluate each term.
- TARGET: **microfound the old retention wedge** L^p instead of postulating
  it: derive (ℓⱼ, ϖⱼ) from an explicit tax-basis/assessment-cap or
  mortgage-lock-in primitive, so the policy counterfactuals map to statutory
  parameters rather than reduced-form wedges.
- TARGET: supply-elasticity pass-through: how η attenuates each wedge's
  fertility effect (comparative static with a testable implication).
- STRETCH (gate: CORE proofs survive the adversarial pass): a two-period
  dynamic version where the old wedge arises endogenously in equilibrium, or
  an optimal-policy (Ramsey) characterization over the available instruments
  with a targeting principle.

**Quantitative model (Sections 3–5):**

- CORE: **close the intergenerational loop.** The model is named
  intergenerational but currently has no transmission: estates end in warm
  glow and young inheritance exposure is disconnected. Implement a
  parsimonious inheritance kernel (estates → heir liquid wealth, family-
  linked by parity or anonymous — your design call), with bequest-principal
  adding-up. Without this, the bequest-tax counterfactual and the paper's
  title are unsupported; this is the single most on-thesis structural change
  available.
- CORE: the **renter-side wedge as a modeled primitive** with its own target
  (the theory–quant bridge: rental premium or size cap, stated and
  calibrated, not implicit).
- CORE: a **welfare module** — consumption-equivalent variation by cohort and
  type for every counterfactual, with the planner-vs-CE wedge accounting from
  the theory.
- CORE: apply the validated fix package (`speed_patch/`: kernel clamps, Brent
  refine, hoist, `max_iter_eq=3`) and the objective redesign (relative-
  deviation loss; rung-share/mean-room targets replacing quantized medians;
  frozen owner menu; `tenure_choice_kappa` as a design axis in
  {0.05, 0.10, 0.15}) — then re-rank all existing records (zero compute) and
  calibrate to a **defensible baseline**: lifecycle ownership no longer
  inverted, strict clearing, full target table. Warm starts: the re-ranked
  best and `case_018` from
  `speed_patch/results_minicalib/minicalib_cases.jsonl` (prime-age ownership
  0.743 — the most economically sensible recorded point).
- CORE: the three counterfactuals (parent-targeted credit relief,
  property-tax change with re-clearing, bequest tax) — now with the
  inheritance kernel the bequest-tax experiment becomes a real GE experiment
  with revenue rebates and adding-up, not a terminal wedge. Decompose each
  through the theory's wedge terms.
- TARGET: **formal sensitivity/identification** — a local sensitivity matrix
  (Andrews–Gentzkow–Shapiro style) of parameters to moments at the baseline,
  computed numerically; plus the which-moment-disciplines-what narrative.
- TARGET: if the timing margin still binds after the objective fix, implement
  age-dependent fertility utility (authorized) and report whether
  `mean_age_first_birth` becomes targetable; document the fertile-window
  floor either way.
- STRETCH (gate: baseline accepted at P3): a simple perfect-foresight
  **transition path** for one headline policy, or survival risk so estates
  occur across ages rather than only at the terminal age.

**Paper craft:**

- CORE: full draft (45–70 pp + appendix), one voice, JPE/Menzio discipline,
  $(n,s)$ notation, model–code correspondence table, every figure/table
  regenerable by one command (`make_paper.sh`).
- CORE: **referee simulation**: three independent adversarial referee agents
  (theory specialist, quantitative-macro specialist, housing/urban
  specialist) each produce a full report on the compiled draft; you write a
  response memo and revise. Iterate until no referee finds a correctness
  error (style objections may remain; log them).
- TARGET: a 20-minute seminar deck built from the paper's figures.
- TARGET: a "three most marketable results" memo for the user (what to lead
  with on the market).

## Hard rails (non-negotiable, unchanged by ambition)

1. **Worktree first.** Before any substantive work, create an isolated git
   worktree on branch `jmp-build-2026-06` and do ALL repo work there. Never
   touch the main checkout — the user is editing `main` in parallel. Record
   your base commit; do not pull or rebase mid-session; log expected
   divergence instead of resolving it.
2. **Commit locally, clear messages, NEVER push, never force-anything.**
   Integration is the user's job (this overrides the repo's push routine).
3. **Never modify** `latex/intergenerational_housing_fertility_v4.tex` (your
   paper is new files under `latex/jmp_draft/`), `CALIBRATION_STATUS.md`,
   `CLAUDE.md`/`AGENTS.md`, or anything under
   `output/claude_capability_test_2026-06-09/` (read-only evidence).
4. **Compute:** local only, no cluster. Interactive runs ≤ 45 min each,
   smoke-tested, checkpointed. Bigger searches are allowed **only** as
   written run designs in STATUS.md (grid size, eval cost, wall-clock
   estimate, checkpoint artifacts, stop criteria) and ≤ 8 h overnight local,
   at most one such run at a time, fully recoverable from checkpoints. A
   30-minute checkpoint silence means the run is unhealthy — kill and
   diagnose, don't wait.
5. Honesty rails: verify paper claims against the actual PDFs/text; never
   compare losses across objective conventions; the loss noise floor and
   non-clearing-region facts from the audits must be respected in any
   ranking claims.

## Startup protocol (in order)

1. Repo protocol: `memory/AGENT_MEMORY.md`, latest `memory/daily/`,
   `CALIBRATION_STATUS.md` (dt-strand background), `git status -sb`.
2. Ground truth for this mission — read fully:
   - `output/claude_capability_test_2026-06-09/HANDOFF_SUMMARY.md`
   - `output/claude_capability_test_2026-06-09/theory_clean_room_note.tex`
   - `output/claude_capability_test_2026-06-09/speed_patch/MINI_CALIB_REPORT.md`
     + `PATCH_NOTES.md` + `patch_suggestions.diff`
   - `code/model/intergen_housing_fertility/IMPLEMENTATION_STATUS.md`
   - `docs/model/intergen_housing_fertility_calibration_structure_20260608.md`
     and the three `docs/model/intergen_*_20260609.md` audit notes
3. Workspace inside the worktree: `latex/jmp_draft/` (paper),
   `code/model/intergen_housing_fertility/` (edit in place),
   `output/jmp_build/` (runs, each with README + diagnostic packet),
   `docs/jmp_build/STATUS.md` (living diary: base commit, framing memo, plan,
   run ledger, decisions taken, decisions left for the user, escalation-menu
   scoreboard).

## Working style

- Decompose with task tracking; use multi-agent workflows aggressively for
  parallelizable work (literature verification one-agent-per-paper,
  adversarial proof-breaking one-agent-per-proposition, code review of every
  structural change, robustness sweeps, the referee panel). Keep theorem
  statements, framing, and the paper's voice in the main context.
- Every accepted model solution gets its visual diagnostic packet. Every
  structural change gets: IMPLEMENTATION_STATUS entry + targeted test +
  re-run of the regression gates (`speed_patch/run_patch_tests.py` T1/T5
  pattern + the zero-floor check) + a before/after moment table.
- When a genuine fork needs the user (data-dependent target changes, >8 h
  compute, abandoning a CORE item), write it to STATUS.md under "Decisions
  for the user", take the most defensible default, proceed reversibly.
- Failure protocol: if a CORE item proves impossible within rails, do not
  silently downgrade — document why, with evidence, and redesign the paper
  around what is true. An honest strong paper beats an oversold weak one.

## Phase gates (each ends with a STATUS.md entry + a worktree commit)

- **P0 Setup**: worktree, workspace, base commit, framing memo, plan,
  escalation-menu scoreboard initialized.
- **P1 Theory**: CORE theory complete with proofs; adversarial proof pass
  survived; aggregate-decomposition proposition stated in
  quantitative-evaluable form; figures 1–2 final.
- **P2 Model engineering**: fix package + objective redesign + inheritance
  kernel + renter-wedge primitive + welfare module implemented, tested,
  documented; re-rank memo written.
- **P3 Baseline**: calibration within compute rails; baseline accepted with
  full target table, diagnostic packet, sensitivity matrix; lifecycle
  ownership verdict stated plainly.
- **P4 Counterfactuals + welfare**: three policies run, decomposed through
  the theory's wedges, CEV tables built; estate-tax experiment with rebates
  and adding-up.
- **P5 Paper**: full draft compiles end-to-end from `make_paper.sh`;
  internal-consistency pass; **referee simulation + revision rounds**; final
  STATUS.md with integration guide and the marketable-results memo.

Begin with the startup protocol now, post the framing memo and plan to
STATUS.md, and proceed through the phases without waiting for further input.
