# Launch guide — JMP build session

How to run the autonomous JMP session defined in `MASTER_PROMPT.md`, in
parallel with your own work on `main`, without the two interfering.

## 1. Launch

Open a **new session** in the repo root and send, as the first message:

```text
ultracode read docs/prompts/jmp_autonomy_2026-06-10/MASTER_PROMPT.md — that is
your prompt. Execute it. Work only in your own git worktree; my checkout of
main is in active use in parallel.
```

Why each piece matters:

- **`ultracode`** — the keyword opts the session into multi-agent Workflow
  orchestration for its whole duration, so it will fan out subagents for
  literature verification, adversarial proof-breaking, code review of every
  structural change, robustness sweeps, and the three-referee simulation by
  default. Without it, the session must ask before each multi-agent run.
- Append a **token budget directive** like `+2m` if you want the fan-outs
  scaled generously — recommended for this mission, which is multi-day by
  design (referee rounds, overnight calibration runs).
- The worktree sentence is redundant with the prompt's hard rails — keep it
  anyway; it is the one safety property you care most about.

Scope expectations: this is the *super-ambitious* variant. The prompt
explicitly authorizes structural model changes inside the worktree
(inheritance kernel to actually close the intergenerational loop, renter-wedge
primitive, welfare/CEV module, age-dependent fertility utility, owner-ladder
redesign), new theory (aggregate decomposition proposition, microfounded
retention wedge, optimal-policy stretch goal), formal sensitivity analysis,
and a referee-simulation phase. It also allows **one overnight local run at a
time (≤8 h, written design + checkpoints)** on top of the ≤45-min interactive
runs — so don't be surprised to find a calibration grinding at 3am; its run
ledger in STATUS.md is the place to check. `phi=0.80` and the no-push /
no-main-edits rails are unchanged. Progress is tracked on an escalation-menu
scoreboard (CORE / TARGET / STRETCH) in STATUS.md, and a documented dead end
counts as a deliverable — expect honest failure reports, not silent
downgrades.

## 2. Let it run long

After the first turn (it will set up the worktree, read the ground-truth
files, and post its plan to `docs/jmp_build/STATUS.md` in the worktree),
type:

```text
/loop
```

with no interval. That puts the session in self-paced mode: it schedules its
own wake-ups, continues through the phase gates between your check-ins, and
keeps going when calibration runs finish. You can close the window; it
resumes itself. To stop it, just interrupt or send any message.

Notes:

- There is no `/goal` skill in this environment; the master prompt + `/loop`
  is the mechanism that plays that role (a standing objective with phase
  gates, plus autonomous continuation).
- You do NOT need to invoke brainstorming/planning skills yourself — the
  session triggers them on its own where they apply.
- `/workflows` shows live multi-agent progress if you are watching.

## 3. Monitor without interfering

- Its diary: `<worktree>/docs/jmp_build/STATUS.md` (base commit, plan, run
  ledger, decisions it took, decisions it left for you).
- Its outputs: `<worktree>/latex/jmp_draft/`, `<worktree>/output/jmp_build/`.
- The worktree path will be printed in its first status entry (a sibling
  directory created by `git worktree add`, branch `jmp-build-2026-06`).
- It commits to its branch only and never pushes. Your `main` checkout is
  untouched by construction; you can keep applying the production fixes from
  `HANDOFF_SUMMARY.md` concurrently.

## 4. Safety properties you can rely on

- All repo edits in the worktree branch; no pushes; no force operations.
- `latex/intergenerational_housing_fertility_v4.tex`, `CALIBRATION_STATUS.md`,
  CLAUDE/AGENTS, and the capability-test evidence folder are declared
  read-only in its rails.
- Compute is bounded (≤45 min per run, ~4 h/day, checkpointed, no cluster).
  Anything bigger lands as a written run design in STATUS.md for you.
- `phi=0.80` fixed; `(n,s)` notation; no fabricated citations; targets always
  shown next to moments.

## 5. When it finishes (or when you want to integrate)

Its final STATUS.md will contain an integration guide. Review on the branch:

```bash
git -C <worktree> log --oneline main..HEAD
git -C <worktree> diff main --stat
```

Then merge/cherry-pick from your side at your own pace. If you applied
production fixes to `main` in the meantime, the merge is yours to drive — the
session was told not to rebase onto your moving main.

## 6. If something goes wrong

- Interrupt any time; the worktree keeps all state; `/loop` can be restarted.
- If a calibration run hangs: its own rules treat a 30-minute silent
  checkpoint gap as unhealthy — but you can also just ask "status?" and it
  will report from its run ledger.
- Worst case, the whole experiment is one `git worktree remove` away from
  gone, with `main` never having been touched.
