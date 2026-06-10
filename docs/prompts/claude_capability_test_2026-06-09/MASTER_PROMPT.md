# Master Prompt For Claude

You are working locally in:

```text
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
```

Use your strongest available reasoning mode. Take your time. I do not want a
fast generic review. I want a serious research-engineering pass over a
quantitative macro/urban/fertility project, with the standards of a top-five or
top-field economics journal.

Your task is to evaluate the project along four margins:

1. the current theory artifact;
2. the full active codebase;
3. the calibration system and current status;
4. the project as a whole.

Do not rewrite active source files or paper files unless I later approve it.
You may create scratch outputs under:

```text
output/claude_capability_test_2026-06-09/
```

## Startup Protocol

First read `CLAUDE.md`, `AGENTS.md`, and
`docs/prompts/claude_capability_test_2026-06-09/CONTEXT_MANIFEST.md`.

Then follow the repository startup protocol exactly:

1. read `memory/AGENT_MEMORY.md`;
2. read the latest `memory/daily/YYYY-MM-DD.md`;
3. read `CALIBRATION_STATUS.md`;
4. read `SESSION_DIARY.md` only if historical context is needed.

Run `git status -sb` before doing substantive work. The worktree may be dirty.
Do not revert or overwrite user changes.

## Working Style

Be technical, skeptical, and precise. Do not produce a warmed-over summary.
Every important claim should be tied to one of:

- a file and line reference;
- an equation, definition, proposition, or proof step in the TeX;
- a command you ran;
- a calibration output whose comparability you have verified.

If you are uncertain, say exactly what would have to be checked.

Use clean economics style: concise primitives, clear timing, explicit
equilibrium definitions, precise propositions, and proofs that distinguish
private optimality, resource feasibility, market clearing, instruments, and
externalities. Avoid odd naming, tangents, and decorative prose. For exposition,
aim for the discipline and clarity of JPE-style macro theory and papers by
Guido Menzio, even though the field here is different.

## Theory Pass

Read the current theory artifact:

- `latex/intergenerational_housing_fertility_v4.tex`
- `latex/intergenerational_housing_fertility_v4.pdf` if useful

Use the archived v3 only for comparison:

- `latex/archive/intergenerational_housing_fertility_v3_20260609/`

The current goal is not to restart from scratch. The project already has a
working draft. Assess whether the current compact analytical model is clean and
correct, and how it should be sharpened.

Required theory checks:

1. State the primitives, timing, agents, choice sets, constraints, and market
   clearing conditions as they currently appear.
2. Derive the young household FOCs and the old household FOCs from the stated
   problems. Do not assume the formulas are right.
3. Derive the constrained planner problem and its FOCs. Be explicit about what
   the planner can and cannot do.
4. Compare the competitive equilibrium and planner FOCs.
5. Evaluate the efficiency proposition(s). Are they true as written? Are the
   assumptions sufficient? Are any objects being smuggled in as primitives?
6. Be especially careful about the housing wedge:
   - If the down-payment/friction is an individual borrowing constraint, can a
     constrained planner eliminate it, or only if it has transfers/instruments?
   - Is the old-age property-tax capitalization term a social wedge, a private
     wedge, or an incidence term?
   - Which wedge affects the intensive fertility choice, and which affects
     entry through the outside option?
7. Check that warm-glow bequests and the outside option are still present and
   integrated cleanly, not merely appended.
8. Assess the two candidate figures in `latex/figures/fig6_*` and
   `latex/figures/fig7_*`. Are they conceptually correct? What should be
   changed before they enter the paper?

If the current theory is flawed, give the smallest clean correction. If it is
basically right, give a sharper theorem/proposition statement and proof outline.

## Codebase Pass

Audit the active Python model, beginning with the files named in
`CONTEXT_MANIFEST.md`. Do not rely on archived MATLAB or old outputs unless
needed for reference.

Required code checks:

1. Map the active code path: parameter construction, solver, objective,
   calibration driver, result collection, diagnostics.
2. Link the most important model equations or economic objects to code
   locations.
3. Identify any clear bugs, stale assumptions, silent normalizations, target
   mismatches, or missing diagnostics. Use file/line references.
4. Run the smallest safe smoke checks you can identify from local documentation.
   Do not launch long searches. Anything expected to exceed 30 minutes requires
   a written run design and explicit approval.
5. Propose focused tests or assertions that would catch the highest-risk
   failures.

## Calibration Pass

Treat `CALIBRATION_STATUS.md` as the starting source of truth. Then verify what
you can against the code and outputs.

Required calibration checks:

1. Extract the current target system from active code and compare it to
   `CALIBRATION_STATUS.md`.
2. Identify the active calibrated parameters and objective weights.
3. If cheaply available, report current model moments next to target values. If
   not cheaply available, explain exactly what command or data is missing.
4. Do not compare losses across incompatible target systems, room scales,
   closures, or geographies.
5. Diagnose identification: which moments discipline which parameters, which
   parameters are likely substitutable, and which margins are weakly pinned
   down.
6. Propose the next calibration experiment with:
   - run size;
   - smoke test;
   - checkpoint artifacts;
   - stop criteria;
   - expected wall-clock risk.

## Project-Level Pass

Assess whether the project is organized so that the next serious work session
can proceed without losing state. Look for:

- source-of-truth conflicts;
- stale documentation;
- untracked files that appear important;
- generated output mixed with active source;
- missing README/status notes;
- fragile workflows.

Do not turn this into cosmetic cleanup. Focus on issues that could cause wrong
economics, wrong calibration, or lost reproducibility.

## Deliverable

Write your main report to:

```text
output/claude_capability_test_2026-06-09/claude_full_pass_report.md
```

Follow the structure in:

```text
docs/prompts/claude_capability_test_2026-06-09/OUTPUT_SCHEMA.md
```

Also read:

```text
docs/prompts/claude_capability_test_2026-06-09/SCORING_RUBRIC.md
```

At the end of your response to me, include:

- the path to the report;
- commands/tests run;
- the three most important findings;
- the one next action you recommend.

Do not claim a check was run unless it was run. Do not hide uncertainty.
