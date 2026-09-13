# Codex worker task

## Goal

Build an isolated, minimal prototype applying the official sequence-space
Jacobian toolkit to the current original-household-queue transition. The author
explicitly requested an agent to try this alongside the existing transition
jobs. Deliver implementation and a small validation, not another broad methods
essay. Stop after the worker_fast profile's30-minute wall-clock budget and
return the best available result; no silent restart or scope expansion.

## Scope

Exclusive source ownership: `code/model/tools/e5f_sequence_space_prototype.py`
and, only if needed for substantive checks,
`code/model/tools/test_e5f_sequence_space_prototype.py`.
Exclusive note: `docs/model/e5f_sequence_space_prototype.md`.
Scratch/results: `output/model/e5f_sequence_space_prototype_20260913/`.
Do not edit existing scientific modules. The lead will verify the adapter
against the actual residual equations before trusting or promoting results.

## Context

Full project startup required, using bounded excerpts and canonical current
status rather than historical broad scans. Read `docs/model/transition_terminal_method_review.md`
and the leading entries of `CALIBRATION_STATUS.md`.

The model has a verified original2007 stationary household distribution and a
verified new stationary endpoint after a permanent preference decline from
0.1489153145785918 to0.09221854783921073. No immigration. Original four-vintage
household-entry queue, adjusted births/2.1, same population law throughout.
Housing supply curve retains absolute scale and elasticity0.63. Payroll tax
0.179, endogenous balanced pension,1% annual property tax equally rebated.
Do not reset population mass to one along a transition. Unit-mass stationary
distribution is only a way to derive averages and the equilibrium scale.

Current root unknowns are all dated log house prices, pensions and rebates.
Residuals are housing relative imbalance and200 times each fiscal relative
imbalance. The retained finite-root tolerance is2e-4 on this scaled vector,
fresh numerical replay2e-10; keep all native household/feasibility/mass gates.
The100-period stationary-endpoint root takes approximately28 minutes per full
mapping; two mappings have reduced its score44.55 to14.29. Ten-period native
mapping about3 minutes; an exact ten-period no-shock loop just passed.

Verified frozen base on Torch:
`/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/afternoon_original_queue_20260913a/spec.json`.
Use `run_e5f_original_queue_experiments.load_context(spec)` from that batch's
`source/` folder. It verifies all original scientific hashes and cache proof.
Existing interface examples locally:
`code/model/tools/e5f_original_queue_experiment.py`,
`code/model/tools/run_e5f_original_queue_experiments.py` (fixed_terminal_path),
`code/cluster/run_e5f_fixed_terminal_horizon.py`.
Frozen verified endpoint:
`/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/afternoon_original_queue_20260913a_terminal_restart_v1/endpoint/terminal.pkl.gz`
with sibling `root_receipt.json`. Local tensor-free receipt:
`output/model/e5f_original_queue_20260913a/terminal_restart_v1/verified_endpoint_receipt.json`.

Inspect official package source, especially nonlinear path updates, different
initial/terminal steady states, and block/Jacobian interfaces:
https://github.com/shade-econ/sequence-jacobian
https://raw.githubusercontent.com/shade-econ/sequence-jacobian/master/src/sequence_jacobian/blocks/block.py
https://straub.scholars.harvard.edu/sites/g/files/omnuum7751/files/straub/files/sequence_space_jacobian.pdf

Choose the smallest truthful integration. A custom block using the exact native
household/population operator is acceptable. Preserve entry queues and changing
population levels explicitly. Do not pretend ordinary finite-difference
Jacobians automatically obtain the package's fast-news performance. Distinguish
the practical native bridge from a full fast Jacobian implementation, and name
the missing derivative objects if the latter is not feasible within this task.

## Do not touch

No modifications to original solver, utility, targets, parameters, weights,
entry rules, fiscal closure, acceptance gates, frozen sources, running jobs,
presentation or manuscript files. No commits or pushes. No global package or
environment changes. If installing the official package is useful, use an
isolated temporary target/environment and record version or commit. No broad
literature/codebase audit or automatic subordinate agents.

## Required output

1. Minimal adapter/prototype with equations and data-shape contract documented.
2. A zero-shock check and derivative/check comparison where feasible; distinguish
   native checks from toy/interface-only checks.
3. A recoverable ten-period benchmark design: solve count, parallelism, runtime
   estimate, preserved convergence gates and the exact next command.
4. Compressed final report: files, what works, tests/evidence, what is missing,
   whether a larger benchmark is worth running. No claim of speedup without a
   measured comparison, or of equilibrium merely because a linear solve works.

## Verification

Pure bookkeeping/interface tests locally are fine. Do not run the native model
on this Intel Mac or on the cluster login node. You may submit at most ONE
small native directional-derivative/smoke job, limited to10 minutes,1CPU,24GiB,
at most six full evaluations and horizon2, using a separately named scratch
batch, all hashes pinned. Torch account torch_pr_570_general; numerical Python
/share/apps/anaconda3/2025.06/bin/python; all numerical-library threads1.
No full ten-period Jacobian grid, calibration, or long root solve before lead
review. If time is insufficient, prepare the exact native smoke rather than
claiming it ran. Never relax a gate to make the prototype pass.

## Stop and report if

The original population law cannot be represented by the chosen package block;
the adapter changes the native residual or normalization; the correct frozen
files/package cannot be accessed; a package integration needs a broad rewrite;
or the30-minute worker deadline is reached. Report the concrete bottleneck and
the smaller implementation that remains possible.
