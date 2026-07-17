# Launch guide — active-focus `/loop` session

A standing, self-paced session that keeps making safe, useful progress on the
project's active focus between your check-ins. It works in an isolated **full
copy** of the code, asks before any cluster compute, and never pushes.

## What `/loop` actually does

- `/loop 30m <prompt>` — **interval mode**: the harness re-sends `<prompt>` every
  30 minutes. Good for dumb polling.
- `/loop` with **no interval** — **self-paced mode**: the session schedules its
  own wake-ups, continues between them, and resumes itself even if you close the
  window. This is what we use. The session decides its own cadence.

The trick that makes it reliable: each wake-up can start with no memory of the
last one, so the prompt forces the session to re-read its instructions and its
diary from disk every time. That's why the prompt lives in a file.

## Launch (two steps)

1. Open a **new session** in the repo root and send, as the first message:

   ```text
   Read docs/prompts/active_focus_loop_2026-06-29/MASTER_PROMPT.md — that is your
   prompt. Execute it now: set up your isolated workspace copy and do the first
   increment. My checkout of main is in active use in parallel; do not touch it.
   ```

   It will create the workspace copy, smoke-test it, post its first diary entry,
   and schedule its first wake-up.

2. Then send:

   ```text
   /loop
   ```

   with no interval. That puts it in self-paced mode. You can close the window.

Optional: prefix the first message with `ultracode` only if you want it to use
multi-agent fan-outs by default — for a careful, mostly-local active-focus loop
you usually don't need it, and it costs more tokens unattended.

## Where to watch it

- Its memory/ledger: `../Fertility_Spring26_loop_workspace/LOOP_DIARY.md` (a
  sibling folder of your repo). The pinned top block shows the current focus and
  any **PENDING USER DECISIONS** waiting on you.
- Its work and runs: everything under `../Fertility_Spring26_loop_workspace/`.
- Nothing it does writes into your primary checkout. Reading your live
  `CALIBRATION_STATUS.md` / memory is the only contact, and it's read-only.

## How it stays safe

- **Isolated copy, not a worktree** — a full `rsync` copy with `.git` excluded,
  so there is no remote and no way to push. Your `main` checkout is never edited.
- **Cluster is gated** — it can monitor `squeue` and collect finished results,
  but it will stop and ask before submitting/cancelling any job or spending
  scratch compute. Approve in chat when you want it to launch.
- **Structural model changes, target-system changes, and `phi`** are all gated —
  it proposes, you approve. Standing project rules still hold.
- **Bounded turns + a no-progress guard** — one increment per wake, and it stops
  spinning if three wake-ups make no real progress.

## Driving it

- **Approve a gated request:** just answer the PENDING USER DECISIONS item in
  chat (e.g. "yes, launch that 4-hour search"). It picks up the approval on its
  next wake (or sooner if you message it).
- **Redirect it:** send a normal message ("switch focus to the old-age ownership
  pathology"). Any message interrupts the loop; it re-orients on the next turn.
- **Stop it:** send any message, or just don't restart it. All state is on disk
  in the workspace, so you can resume later by sending `/loop` again.
- **Ask for status:** "status?" — it will summarize from the diary.

## Integrating its work

Because it works in a copy, nothing reaches `main` until you pull it over.
Review with a diff against your checkout when you're ready:

```bash
diff -ru --exclude='.git' --exclude='output' \
  ~/Desktop/Projects/Fertility/Fertility_Spring26 \
  ~/Desktop/Projects/Fertility/Fertility_Spring26_loop_workspace | less
```

Then copy over only what you want. Worst case, the whole experiment is one
`rm -rf ../Fertility_Spring26_loop_workspace` away from gone, with `main` never
having been touched.
