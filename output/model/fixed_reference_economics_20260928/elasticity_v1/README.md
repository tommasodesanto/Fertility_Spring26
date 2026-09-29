# Local housing-price elasticities: union-grid successor

**September 29 recovery authorization:** after the lead proposed completing
only the two +1% cases with a new two-solve, 15-minute-per-case, 35-minute-total
budget, the author instructed the lead to proceed. The cheaper workers completed
the saving repair and renderer; the lead reviewed the source diff and verified
zero-solve Torch smoke **18814979** passed. Recovery production job **18815133**
is submitted (pending priority at 15:10 EDT), with dependent zero-solve rendering
job **18815186**. Hourly monitoring is active. See [recovery_v1](recovery_v1/README.md)
for immutable sources, input pins, and the separate recovery contract. No old
deadline is reset and no passed case is rerun. The stopped-run status below is
historical; it does not describe the newly authorized recovery.

The new recovery retains the same 602-node grid, frozen preferences including
psi, earnings, inherited distribution, entry, fiscal inputs, supply primitives,
repayment rules and all scientific gates. The two experimental cases prescribe
house price and mapped rent at 1.01 times reference; the second also retains the
already-specified removal of artificial borrowing limits. Checkpoint writing
will be atomic and precede plotting, with explicit saving progress. Completed
0.99 and 1.00 cases remain read-only inputs. The intended result is a verified
three-price comparison; ±2% robustness remains outside this recovery budget.

**September 29, 12:17 EDT: numerical work stopped; full elasticity comparison
incomplete.** Job 18801318 completed all four q0 controls and stopped on its
forecast. Continuation 18803216 completed `reference_990` and `credit_990`, then
hit the unchanged 600-second case limit in `reference_1010`. Its last progress
label is `standard_17_plot_rendering`, which also covers subsequent provenance
checks and checkpoint/receipt writing. A checkpoint exists, but no completion
receipt certifies that case. The bounded Torch audit found all 17 plots and
both gate files, but its 185,708,331-byte gzip checkpoint is truncated (EOFError
before the end-of-stream marker). Its write time immediately precedes the
timeout: checkpoint finalization, rather than model nonconvergence or a
documented scientific gate failure, prevented completion. Do not treat it as
a passed result or load the partial checkpoint for a continuation.

Six cases have passing receipts; a seventh lifecycle case was attempted.
The +1% credit case and all ±2% cases were not attempted. No centered elasticity
or full price-response figure is available. The original absolute deadline
and all solve/case limits remain unchanged; no retry or new solve was launched.
The hourly monitor is paused after the terminal failure. A bounded zero-solve
audit has authenticated the passed cases and diagnosed the incomplete output.
The [partial readout](partial_readout_v1/README.md) contains the one-sided
0.99-to-1 log elasticities, full tables and four unique 17-plot packets. Impact
birth elasticities are −0.436 under baseline credit rules and −0.468 with
solvency-only credit; cohort completed-fertility elasticities are −0.532 and
−0.551. The lead checked the extraction formula, common inherited-state
definition and receipt identities. These are one-sided prescribed-price
responses; symmetry, step-size robustness and GE impact remain unverified.
The earlier verified borrowing GE and supply endpoints
are unaffected. Launch descriptions below are historical records, not current
running status or instructions to resubmit.

**Conditional continuation submitted: Torch job 18803216**, dependent on clean
completion of 18801318. It uses `source_v4/`, logs to `continue_v4.log`, and
writes `results_v2/solve_v4/` under the same remote root. It can run only after
the four successful q0 cases and a forecast-only stop. No numerical case is
repeated and the original absolute deadline remains **12:43:42 EDT on
September 29** (epoch `1790700222.6872504`). V3 was never launched.
Do not resubmit either job. The hourly heartbeat now monitors this workflow;
the earlier borrowing GE remains complete.

The lead reviewed the v3 controller and the v4 metadata/identity diff, and
matched both v4 source hashes to Torch before submission. Torch syntax and
synthetic partial-formatter tests passed (84 comparison and 28 elasticity
rows). Controller SHA256 is
`a2f66adaa6331ec4217f7e5fe07e64978b6c85fcf3b046642f2d87f5d595d0a5`;
launcher SHA256 is
`6d042d9d7e5dc97f1336c234a13132b8d27ed5217c186cda76abfdcf56e0965e`.
At submission, the larger-grid reference control and its fresh repeat passed
in 527.622 and 591.879 seconds. The repeat matches 113 arrays, all 14 fit rows,
31 parameter estimates and 17 standard plot hashes exactly. Completed
fertility is 2.1008623 on the 602-node grid versus 2.1007842 on the earlier
262-node grid. Only the numerical grid changes; the frozen calibration is
not replaced. The first credit control was running.

**Running: Torch job 18801318**, using the frozen `source_v2/` and [recorded plan](plan_v2.json). The [launch receipt](launch_v2.json) records the original 90-minute deadline; no retry or extension is authorized. Remote results: `/scratch/td2248/projects/fixed_reference_elasticity_v2_20260929/results_v2/solve_v2/`; log: `/scratch/td2248/projects/fixed_reference_elasticity_v2_20260929/solve_v2.log`. Plan SHA256: `6354bafd069fbfc4eda21271fc91747a22e3aeb8c90a57293f7dfde01badf28f`.

**Do not launch the ten-case `source/` v1 controller or submit a duplicate v2 job.** The v1 zero-solve preflight found severe price-dependent floor tightening on the fixed 262-node grid. That source and receipt are preserved as a rejected numerical design.

Reference: **2007 stationary reference — block0506, September 28 verified export**.
Frozen manifest SHA-256: `147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4`.
The only proposed economic changes are prescribed house prices with rents mapped by the unchanged user-cost rate, and the solvency-only credit regime. Preferences including child benefit `psi_child=0.1355551166583114`, earnings, entry endowments, survival, fiscal/estate rules, and the housing-supply function are retained. No fertility normalization occurs. The prescribed-price cases do not certify housing-market clearing, general equilibrium, or demographic renewal.

V2 constructs one common grid from the union of the original 160 nodes and natural-credit boundary nodes at all five factors 0.98, 0.99, 1, 1.01 and 1.02. Every price and regime uses that identical grid, preserving entrant and inherited point masses exactly. Four q0 cases first solve reference and credit, each with a fresh exact repeat. The old 262-node q0 checkpoints serve only as a numerical-refinement comparison. Eight shocks follow as paired ±1% cases, then paired ±2% cases. Each case is a fresh process with one lifecycle solve, 600-second cap; twelve solves and 5,400 seconds are the total caps. One CPU/16 GiB, no retries. After the q0 controls, a conservative forecast from observed case wall times stops before another case if the remaining budget is insufficient.

The [v2 preflight receipt](preflight_v2/preflight.json) from Torch job `18801218` passed in 36 seconds, with zero solves. Its [common grid](preflight_v2/common_grid.json) has **602 nodes**; each individual price grid has 262. Across all five prices, inherited and entrant unsupported mass is zero, the independent solvency recurrence passes, and the maximum union-grid floor tightening equals the candidate-specific maximum of 0.0016005. The elementwise extra tightening is exactly zero. [Frozen source hashes](source_v2.sha256) identify the staged v2 bundle. Preflight receipt SHA-256 is `17b4030fb5e1923f3117a0394624c2777ea6c6364e0ae9af5405dcc3a6b784b8`; common-grid SHA-256 is `cc0c21e5540c155dcf09ac1f0a6a2a6ee12d48c1426de33620b97241bda291e1`; driver SHA-256 is `5d461d4999c428c77438ee43dd63db0e69ee0dd10356b0e09d65c62880adc5ca`. A linear extrapolation from an approximately 140-second 262-node solve gives roughly 322 seconds per 602-node solve, or 3,864 seconds for twelve solves before audits and plots. The actual q0 wall times and conservative controller forecast determine whether the 5,400-second total can finish.

The zero-solve preflight checks the independent minimum-over-tenures solvency recurrence at all prices, original atom embedding, entrant/inherited economic and discrete support, and the elementwise numerical-floor difference between the fixed grid and a candidate-specific grid. Candidate-specific grids are diagnostics only. The preflight refuses any price with unsupported inherited or entrant mass. It does not establish policy convergence or the occupied saving-floor incidence; the latter comes from each solved credit audit. Torch job `18800762` completed in 39 seconds with zero solves. [Its receipt](preflight_v1/preflight.json) passes the recurrence, exact q0 grid, and all five entrant/inherited support checks. The source hash list is [source_v1.sha256](source_v1.sha256).

The fixed grid causes potentially material numerical floor tightening away from q0. Maximum tightening relative to the economic floor is 0.449435 at 0.98q0, 0.398416 at 0.99q0, 0.001601 at q0, 0.075448 at 1.01q0, and 0.138157 at 1.02q0. Candidate-specific grids would have maximum 0.001601 at each price. Worst elementwise extra tightening from holding the grid fixed is 0.448835, 0.398316, 0, 0.074248, and 0.136657 respectively; the largest differences occur in owner tenure index 5. These values are numerical support restrictions, not occupied policy responses. They may affect derivative interpretation despite zero unsupported inherited or entrant mass. The ten-case lifecycle run is still unlaunched.

The rejected v1 Torch source remains frozen at `/scratch/td2248/projects/fixed_reference_elasticity_20260929/source_v1/`. The large reference and old credit checkpoints remain on Torch. The successor source and results belong under `/scratch/td2248/projects/fixed_reference_elasticity_v2_20260929/`.

The following command was executed once for job 18801318; it is recorded for provenance, not resubmission:

```sh
ssh torch 'sbatch --output=/scratch/td2248/projects/fixed_reference_elasticity_v2_20260929/solve_v2.log --error=/scratch/td2248/projects/fixed_reference_elasticity_v2_20260929/solve_v2.log /scratch/td2248/projects/fixed_reference_elasticity_v2_20260929/source_v2/launch.sh'
```

The launcher checks the preflight, old checkpoint hashes and frozen sources, creates a hash-pinned plan, then runs at most twelve sequential cases. `latest_completed.json`, `best_so_far.json` and `progress.json` are updated per case. A budget stop writes `stopped.json` and retains partial cases without calling the elasticity table complete. `comparison.csv` separates impact (identical inherited distribution) from recomputed normalized-cohort outcomes. `elasticities.csv` reports centered log elasticities at ±1% and ±2%, plus lower and upper one-sided slopes, only for strictly positive outcomes. Each case retains complete target fits, parameter table, gates, checkpoint on Torch and the standard 17 plots. The best-so-far criterion is only smallest absolute prescribed-price housing excess, never GE certification.

The lead submitted the v2 main job as Torch `18801318`; log `/scratch/td2248/projects/fixed_reference_elasticity_v2_20260929/solve_v2.log`. At the bounded preparation read, the first reference q0 case passed in 527.622 seconds and the fresh reference repeat was running. The lead owns monitoring and collection.

`source_v3/` is an unlaunched first continuation draft, preserved for review. **Use `source_v4/` for any continuation.** It prepares a controller-only continuation if and only if job `18801318` completes cleanly after the v2 conservative wall-time forecast stops it immediately after all four q0 controls. It hard-pins the original plan SHA-256 `6354bafd069fbfc4eda21271fc91747a22e3aeb8c90a57293f7dfde01badf28f`, job ID, original absolute deadline `1790700222.6872504`, and frozen v2/v4 source bytes. It does not change the v2 solver, plan, scientific gates, 600-second case cap, or twelve-solve total. It verifies the four receipts, checkpoints, 14/31 tables and 17 plots; symlinks those completed case directories through a read-only bind; and calls only the unchanged v2 child driver for remaining cases. A regime-specific forecast reserves 1.25 times each regime's maximum observed case wall time, capped at 600 seconds. It starts a paired ±1% group only if all four cases fit the remaining original budget; paired ±2% cases are optional under the same test. If only ±1% completes, it writes explicitly partial ±1% centered and one-sided elasticity tables with v2-compatible columns, source/receipt pins and CSV hashes, and no five-price completion claim. Full `completed.json` carries the reference label and v2-compatible completion metadata. [V4 source hashes](source_v4.sha256) identify this unlaunched preparation. It refuses scientific failure, a nonterminal v2 job, changed pins, or any other stop reason. Torch syntax and synthetic partial-formatter checks passed: 84 comparison rows, 28 elasticity rows, with verified columns and hashes.

Only if the lead reviews the four actual q0 costs and the v2 run ends in that specific forecast stop, the prepared command is:

```sh
ssh torch 'sbatch --dependency=afterok:18801318 --kill-on-invalid-dep=yes --output=/scratch/td2248/projects/fixed_reference_elasticity_v2_20260929/continue_v4.log --error=/scratch/td2248/projects/fixed_reference_elasticity_v2_20260929/continue_v4.log /scratch/td2248/projects/fixed_reference_elasticity_v2_20260929/source_v4/launch_continue.sh'
```
