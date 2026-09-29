# Codex worker task

## Goal

Implement an isolated, reviewable removal of arbitrary renter debt taper.
Author explicitly selected removal. Do not change shared active model code.
Keep zero NEW unsecured borrowing as comparison control; positive credit is
an open author choice, not adopted. Retain lifetime repayment/nonnegative estates.

## Scope

Exclusive write ownership ONLY `output/model/fixed_reference_economics_20260928/credit_no_taper_v1/`.
Prepare one patch/application driver, targeted tests and compact README. Use
frozen source via existing Torch mount rather than copying repo/checkpoint.
Optional patch files/full changed source may remain in this folder; never edit
`code/model/` or another worker folder. Read native `parameters.py`, `solver.py`,
`kernels.py` in intergen_eqscale_seq_optimized and runtime named by frozen contract.

## Context

Full startup required perAGENTS. Reference **2007 stationary reference — block0506,
September 28 verified export**. Manifest at output/model/fertility_identification_20260928/fixed_reference_manifest.json.
Original frozen Torch root `/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project`.
solverSHA b637a655a9344b63f4461ee0fa4796c04bd98188477c4e6ace2c48ae0fc8aec1;
parametersSHA66f86697c2c58ca3864305bf13dd2be71a008905b2beb573f1a4ebafabef5464;
kernelsSHA639c9a21797dbc9f2a0e9a891f283c115353c2edfcb89c959a7fe9f32b86ca27.
No baseline switch; preserve frozen source/checkpoint and comparison permanently.

Current renter floor `min(s[j+1]*min(b,0),-D[j+1])`; lambda_d0 means no new
credit but existing debt may roll over. Arbitrary weights taper42–62. Allnative
renter kernel paths use these scalar weights/floor. Simply weights1 until
terminal0 is NOT sufficient: native baseline renter branch lacks a separate
positive-mortality estate mask. Natural-credit mode has a separate strict
solvency solver but enabling it would remove OTHER credit restrictions, forbidden.

## Mathematical specification owned by lead

Renter saving b' before death realization: floor=min(b,0) when death impossible,
floor=0 when current decision has positive death probability or terminal death.
Death condition exactly native `j==J-1 or(use_age_survival and survival_probs[j]<1)`.
Keep lambda_d0 and reject positive credit in this isolated mode until specified.
No age42–62 schedule. Owner buyer/incumbent/estate credit and purchase cash
conditions remain byte-for-byte or mathematically identical. Retain all other
parameters incl psi_child, entry/earnings/fiscal/supply/targets/grids/timing.

Suggested minimal path, if exact: add explicit optional renter rule flag to
isolated `build_debt_caps`, default absent preserves frozen outputs. In newmode
build weights[j+1]=1 on no-death decisions and0 on possible-death decisions,
caps0; existing floor/kernel then implements exact rule on all branch paths.
This is a binary estate bound, NOT a postponed arbitrary taper. Verify flag
survives native parameter-rebuild calls and mortality arrays available at
builder time; do not guess. Another minimal implementation is allowed if exact
and all paths align. Do not bypass solvency or enable natural_credit globally.

## Do not touch

No shared files or checkpoint mutations/downloads; no model imports/numerical
tests on Mac; no source hashing/bulk copies locally; no production solves,
Slurm jobs, recalibration, parameter changes, Google writes or Gitcommit/push.
Do not restart failed accounting18832560; unrelated, budget closed. No subagents.
Twenty-minute cap, stop if assumptions, unsupportedpath or identity conflict.

## Required output

Minimal reproducible patch driver validating expected originalsource hashes
before applying only to a new versioned directory; explicit runtime activation
recipe; regression/targeted fixture tests for zero debt, negative debt, all ages,
first mortality, terminal, owner boundaries and absent-flag baseline identity.
Tests must exercise actual changed function/path, not mirrorimplementation.
Include a Torch-only verification launcher prepared for leadreview, NOT submitted:
1CPU/16GiB/5min, zero lifecycle/GE solves, no retries, compact passed/failure
receipt, exact source/test hash manifest. Do not stage new immutableversion until
leadreview. Sourcefile reads/AST parsing can be lightweight local; no imports.
Short README status with every economic change and all uncomputed steps.

## Verification

Lead will inspect model-critical diff line by line and run only exact targeted
verification on Torch after review. Preserve all17standard diagnostic plots.
No inference about behavior from passing fixtures. Future fixedparameter GE
recompute has a separate explicit solve/time budget, not authorized by smoke.

## Stop and report if

Any premise cannot be met without economic changes, cross-path mismatch, or
time limit. Return specific diff and uncertainties; no silent retry or extension.
