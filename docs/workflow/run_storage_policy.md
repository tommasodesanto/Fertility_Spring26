# Run storage policy (Mac)

Read before launching any local run, calibration chain, battery or package
build. Written 2026-10-08 after the Mac fell to 18 GB free during overnight
chains. Rule of thumb: save only what a later step reads.

## Why the disk fills

Measured 2026-10-08 (repo total 358 GB on a 926 GB disk):

| Source | Size | Cause |
|---|---|---|
| `.git/objects` loose, unreachable | 82.7 GB | outputs hashed by `git add -A` / `git stash -u` and then dropped; `gc.auto=0`, so nothing ever prunes them. `tmp/` is not git-ignored. |
| Per-evaluation `solution_arrays.npz` | ~90 GB | every evaluation of a search keeps a 25–80 MB full solution (CES 46 GB, cap55 12 GB, CES composite 11 GB) |
| `diagnostic_packet.pkl.gz` / `native_solve_unverified.pkl.gz` | ~27 GB in `transition_readiness_v1` alone | 1,130 + 373 copies across smoke, test and local runs |
| Package builds already staged on Torch | ~5 GB | Mac copies kept after upload |
| `executed_P.json` | ~20 GB in `credit_mechanism_20261004` | dense arrays written as JSON |
| Code-tree copies per evaluation | 323k `.py` files in one folder | each case stages its own copy of the code |

## Rules

1. **Per evaluation**: write only the parameter vector, the moment/fit table,
   the loss and a status line (JSON/CSV, under 1 MB). No solution arrays, no
   diagnostic packets, no code copies.
2. **Solution arrays**: keep one `best_so_far` (overwritten in place) per
   chain, plus the final verified postcheck. Anything else is regenerable
   from the parameter vector and the code commit.
3. **Provenance**: record the git commit and a hash of any uncommitted diff in
   the run manifest instead of copying the code tree into each case.
4. **Dense arrays**: `.npz` (compressed), never JSON.
5. **Smoke and test runs**: write to the scratchpad or a folder named
   `smoke_*`, and delete it once the real launch is confirmed.
6. **Torch hand-offs**: after a package is staged on Torch and its checksum
   verified there, delete the Mac copy and note the Torch path in the README.
7. **Budget before launch**: state the expected write (evaluations x size per
   evaluation). A local launch needs free space of at least
   guard (15 GB) + expected write + 10 GB headroom. Above 5 GB of expected
   writes, say so to Tommaso before launching.
8. **Disk guard**: every local launcher stops cleanly below 15 GB free, as the
   `local_mac_calibrate_r1_disk15.py` chains already do.
9. **On completion**: the thread that launched a run deletes its per-case
   folders after the summary, best point and postcheck are saved, and writes
   the final size in the experiment README.
10. **Git**: never `git add -A`, `git add .` or `git stash -u` in this repo;
    add files by path. `tmp/` and `output/` hold no tracked files.
    Never run `git prune`, `git gc` or `git stash drop` to free space: deleted
    objects are unrecoverable. A prune needs a verified external copy of
    `.git` and Tommaso's word first (Oct 9 2026: an emergency prune removed
    82 GB before the agreed Torch backup finished, including three unmerged
    stashes). In a disk emergency, pause new writers instead.

## Outside the repo

`~/fertility_runs/` (18 GB) follows the same rules. Personal iCloud data
(WhatsApp backup tars, about 24 GB rewritten nightly) is outside the
project's control; turning on "Optimize Mac Storage" for iCloud Drive keeps
it in the cloud.

The whole Desktop, including this repo, is under iCloud Desktop sync
(`com.apple.icloud.desktop`). iCloud churn can take tens of GB overnight and
leaves `<name> 2` conflict copies inside `.git/objects`.
