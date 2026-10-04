# One-birth Estate-A global exploration, stage 1

This stage changes the **search design only**. Relative to the verified recovery
chain 1, earnings, initial wealth and income distributions, timing, transfers and
floors, preferences, estate accounting, birth cap, solver gates, and the 14 target
values and weights do not change. The pinned target fingerprint is
`c7a3d185668122e508a6c322bc5ef0715ebb0ecb23948c8d9b184ee25d1cde70`;
the weight fingerprint is
`f762ebb5684ab30487b3b8b64fc10977fda396b520035d91c0c5c803255f88e4`.
The verified same-contract incumbent loss is 21.275413361071312; the canceled
nearby-start continuation's 19.14485687996558 is provisional.

Stage 1 draws a scrambled Sobol sequence with `d=10`, `m=6` (64 points), and
seed `20261004`. It covers the full fixed search bounds below. Let `u` be a
coordinate in `[0,1]`: linear coordinates use `lo+(hi-lo)*u`; positive
logarithmic coordinates use `exp(log(lo)+(log(hi)-log(lo))*u)`.

| Coordinate | Bounds | Transform |
|---|---:|---|
| `beta_annual` | 0.93–0.99 | linear |
| `chi` | 0.1–5 | log |
| `child_benefit_curvature` | 0–0.8 | linear |
| `first_birth_fixed_cost` | 0–8 | linear |
| `h_P` | 0.1–2.6 | linear |
| `kappa_fert` | 0.02–50 | log |
| `kappa_fert_continuation` | 0.02–50 | log |
| `psi_child` | 0.01–0.5 | linear |
| `tenure_choice_kappa` | 0.001–0.1 | log |
| `theta0` | 0–8 | linear |

The plan checks 64 distinct points and one point in each of 64 marginal Sobol
strata for each coordinate. This is dispersed exploration, not a guarantee of
finding a global optimum in ten dimensions. The verified incumbent is an
independent control, evaluated once in the native preflight and not counted
among the 64 exploratory points.

The initial production array has 16 one-core tasks with four sequential points
per task, at most 16 concurrent. Each task has a 90-minute Slurm wall limit;
each point has a 20-minute native evaluation allowance and its own saved result.
At the observed 3–5 minutes for an ordinary completed GE, this is roughly
15–30 minutes wall and 4–6 core-hours; hard points may consume the full limit.
No task retries, automatic extensions, or stage-2 submission are configured.
The production array uses a five-minute delayed start so its exactly-once
submission receipt is durable before workers can run.
The exact two-case smoke uses a mock residual and checks loop/checkpoint/output
execution only, under 5 minutes and zero model solves. A separate one-case
incumbent preflight uses the unchanged native GE and reporting gates with a
25-minute limit. Both must pass before production submission.

For each full native result, the runner records the full result and a case
checkpoint, latest-completed summary, best-so-far summary, and heartbeat.
`inadmissible_numerical` and case `budget_exhausted` retain a reason and no
finite loss; an unknown evaluator status, exception, source drift, target drift,
or reporting-gate failure terminates that task. After stage 1, select up to
eight **distinct transformed-space basins** among feasible points plus the
verified incumbent and the canceled continuation's 19.14485687996558
provisional point (subject to fresh native validation) for a separately
reviewed bounded local refinement, with at most 100 calls and ten hours per
one-core job. Boundary coordinates of the incumbent remain eligible even
though Sobol points lie in the interior. The refinement stage is not
automatically launched. A loss below 13 counts only after the unchanged fresh
native selected-point and exact-repeat gates pass. Otherwise stop at the finite
budget and report the full 14-row fit and 10 free-parameter table for any
verified winner.

After the array is terminal, copy `control/plan.json` and each task's
`run/{cases,completed}.json` and launcher receipt into the same relative layout
as the Torch result directory.
Run `python3 collect_stage.py --root <copied-stage-root> --out
<collection-directory>`. The collector checks plan and target pins, summarizes
every task, writes the best provisional fit and parameter tables, and proposes
distinct transformed-space starts for review. It does not launch refinement
or declare a verified winner.
