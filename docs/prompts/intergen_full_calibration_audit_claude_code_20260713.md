# Claude Code prompt: full adversarial audit of the repaired calibration

Use this entire document as the task. Do not ask me to restate the context.

Repository:
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`

## Role and standard

Act as the lead independent reviewer of a quantitative-macroeconomics
calibration intended for a top economics journal. Perform a complete,
adversarial audit of the live one-market intergenerational housing--fertility
calibration after the July 2026 entrant-feasibility repair.

This is a falsification exercise, not an editing exercise. Do not trust status
notes, reported losses, optimizer labels, or prior audit conclusions until you
have reconstructed them from code and raw records. Be willing to conclude that
a reported improvement is a search artifact, a numerical artifact, outside the
admissible model, weakly identified, or economically meaningless.

Use parallel subagents only for genuinely independent scopes such as target
provenance, search forensics, and code-vs-math review. The lead reviewer must
personally verify the objective, target table, parameter mapping, selected
candidate, and identification conclusion before reporting them.

## Non-negotiable boundaries

- Read `AGENTS.md` and follow the repository startup protocol.
- Read, in order: `memory/AGENT_MEMORY.md`, the latest
  `memory/daily/YYYY-MM-DD.md`, `CALIBRATION_STATUS.md`, and
  `code/model/README.md`.
- This audit is read-only with respect to model code, targets, calibration
  status, paper, policies, and canonical result files. Do not implement fixes,
  change bounds, reweight moments, alter targets, commit, push, or open a PR.
- You may write audit outputs only under
  `output/model/full_calibration_audit_20260713_claude/`.
- Put temporary files under that audit folder or `tmp/claude_calibration_audit/`.
- Production standard is `J=17`, `Nb=120`. Do not demand or substitute
  `Nb=240`; the author explicitly chose `Nb=120` as production.
- Do not compare current losses with pre-repair, pre-combined, pre-room-scale,
  or different-target-system losses.
- Never interpret children-at-home as current fertility, and never reinterpret
  the one-shot fertility architecture as a sequential birth-hazard model.
- Do not finalize while the focused profile jobs below are still running. Audit
  static code and completed rounds first, then poll the jobs without cancelling
  or modifying them. If Torch is inaccessible, state that limitation and leave
  all focused-run claims `NOT TESTED` rather than guessing.

## Canonical code and specifications to inspect in full

Active package:
`code/model/intergen_housing_fertility/`

At minimum read and cross-reference:

- `parameters.py`
- `production_profile.py`
- `calibration.py`
- `kernels.py`
- `solver.py`
- `local_panel.py`
- `diagnostics.py`
- all tests under `intergen_housing_fertility/tests/`
- `code/model/tools/run_intergen_feasibility_fix_diagnostic.py`
- `code/model/tools/run_intergen_combined_recalibration.py`
- `code/model/tools/run_intergen_combined_wide_discovery.py`
- `code/model/tools/run_intergen_combined_wide_polish.py`
- `code/model/tools/run_intergen_repaired_surrogate.py`
- `code/model/tools/report_intergen_calibration_refinement.py`
- `code/model/tools/collect_intergen_focused_profile.py`
- the corresponding submit scripts under `code/cluster/`

Settled feasibility specification and earlier audit evidence:

- `output/model/full_audit_20260711/FEASIBILITY_FIX_SPEC_FOR_CODEX.md`
- `output/model/full_audit_20260711/FULL_AUDIT_20260711.md`
- `output/model/full_audit_20260711/diagnostics/jacobian_current_bound.json`

Paper/math comparison, only where needed to test implementation consistency:

- `latex/intergenerational_housing_fertility_paper_draft.tex`
- `docs/model/intergen_one_market_identification_ledger_20260618.md`
- `docs/model/intergen_sensitivity_jacobian_audit_20260618.md`

## Result trees

### Completed repaired production round

- `output/model/intergen_refinement_policy_20260712/calibration/report/`
- Reported incumbent claim: loss `8.925945490698217`, strict residual
  `6.0097429391390856e-05`.

### Completed overnight three-algorithm round

- `output/model/intergen_overnight_calibration_20260713/report_classic/`
- `output/model/intergen_overnight_calibration_20260713/report/`
- `output/model/intergen_overnight_calibration_20260713/README.md`

Claims to falsify:

1. Best production-box point is loss `8.895656229998638`, residual
   `4.3899301847149713e-05`.
2. Best relaxed point independently reproduces near loss `8.415008941173383`,
   residual `4.4298659299111956e-05`.
3. The relaxed point is outside the production box at approximately
   `h_bar_0=0.5300` and `theta_n=2.4546`, with `theta0` almost zero.
4. Projecting only `h_bar_0` and `theta_n` to the production bounds makes that
   vector entrant-infeasible.
5. ExtraTrees proposed 64 points; 49 were entrant-infeasible, 15 solved, 14
   were strict, and none beat the incumbent.
6. Across both wide waves there were 14,355 strict solves. Higher `chi` did not
   help in the sampled records: the best strict point with `chi>1.15` was much
   worse than the best near `chi=1.1`.

Raw completed results are on Torch under:

`/scratch/td2248/projects/Fertility_Spring26_20260711_feasibility_recal/output/model/overnight_calibration_20260713/`

### Focused relaxed polish and h0-by-chi profile (may still be running)

Local run ledger:
`output/model/intergen_focused_relaxed_profile_20260713/README.md`

Remote result root:

`/scratch/td2248/projects/Fertility_Spring26_20260711_feasibility_recal/output/model/focused_relaxed_profile_20260713/`

Jobs:

- Free transformed Nelder--Mead polishes: `13480103`--`13480109`.
- 25-cell profile array: `13480110`.
- Dependent collector: `13480164`.

Profile design:

- deterministic tenure and no bequest (`tenure_choice_kappa=0`,
  `theta0=theta_n=0`);
- fixed `h_bar_0 in {0.25,0.40,0.55,0.75,1.00}`;
- fixed `chi in {0.90,1.00,1.10,1.25,1.50}`;
- reoptimize all remaining active coordinates within every cell.

The very early provisional claims (`8.3764`, then `8.2285`) are hypotheses,
not accepted results. Verify from final raw records and fresh re-solves.

Torch access, if needed:

```bash
ssh torch
squeue -u "$USER"
```

Use account `torch_pr_570_general`; never use project 571.

## Audit questions

### A. Exact calibration contract

1. Reconstruct the actual live objective directly from code. List all 15
   moments, targets, weights, units, and exact loss formula.
2. Trace every target to its empirical source or construction script. Flag
   stale targets, sample mismatches, means-versus-medians, household concepts,
   room-unit inconsistencies, and targets whose code does not match their
   labels.
3. Verify that the production driver, relaxed DE driver, relaxed polish driver,
   surrogate validator, reporting code, and all cluster launchers use the same:
   `J`, `Nb`, income process, `q`, `delta`, `eta_supply`, bequest normalization,
   feasibility repair, debt taper, housing menu, target set, weights, and market
   tolerance.
4. Verify the complete parameter mapping. In particular check
   `beta_annual -> beta_annual**4`, `H0`, fixed versus searched coordinates,
   transformed relaxed coordinates, and every nested arm.
5. State the correct parameter and informative-moment count for each branch:
   production, deterministic tenure, no bequest, and each fixed profile cell.

### B. Entrant feasibility and economic admissibility

1. Re-audit the implemented repair line by line against
   `FEASIBILITY_FIX_SPEC_FOR_CODEX.md`, including the owner unsecured-position
   accounting, seller post-liquidation case, current-iterate-price convention,
   post-sale grid clips, and retirement taper.
2. Confirm that every calibration path catches `InfeasibleThetaError` and never
   ranks a dead-node, placeholder-policy, non-strict, or non-clearing solve.
3. Check whether relaxed `h_bar_0<1` is mathematically and economically valid in
   the implemented Stone--Geary/housing menu, or whether it exploits a grid,
   interpolation, clipping, or utility-domain artifact.
4. Inspect poorest/richest, entrant, renter/owner, parent/nonparent, and oldest
   boundary states at the final relaxed candidate. Verify budget identities,
   borrowing floors, down-payment `(1-phi)` convention, positive-mass
   feasibility, and post-62 debt restrictions.

### C. Search integrity and completeness

1. Inventory every task, wave, seed, arm, fixed coordinate, transform, solve
   count, strict count, infeasible count, exit state, wall time, and best loss.
   Work from raw `cases.jsonl`, metadata, Slurm accounting, and stderr--not
   summary prose.
2. Detect missing/truncated tasks, overwritten result directories, duplicated
   points, stale seeds, configuration drift, non-strict winners, and records
   accidentally counted more than once.
3. Audit the transformed DE and Nelder--Mead implementations line by line:
   simplex construction, clipping, contraction/shrink logic, boundary
   degeneracy, loss ordering, time/evaluation stopping, and checkpoint safety.
4. Decide whether the broad relaxed search genuinely covered the domain or was
   too sparse in 14 dimensions. Quantify coverage rather than saying merely
   "many evaluations."
5. For the focused run, independently reconstruct the 5-by-5 `h_bar_0 x chi`
   surface from raw strict records. Confirm that every cell fixed the intended
   values and genuinely reoptimized the remaining coordinates.
6. Determine whether the final improvement comes from `h_bar_0`, `chi`, the
   bequest ridge, another coordinate, or a numerical interaction. Use nested
   arms and controlled profiles, not simple correlations alone.

### D. Fresh numerical reproduction

Perform fresh, cache-free re-solves in new audit-owned directories for:

1. the July 12 incumbent;
2. the best production-box point;
3. the completed relaxed `8.415` point;
4. the final best focused-polish point;
5. the best fixed profile cell;
6. at least one infeasible projection/counterexample used in the narrative.

For every solve record the exact theta, code hash, overrides, elapsed time,
price, residual, strict flag, feasibility census, all moments, reported loss,
and independently recomputed loss. Repeat the final selected solve at least
twice to test determinism.

Then re-evaluate the principal finalists under a tighter equilibrium evaluator
(`max_iter_eq=40`, `tol_eq=2.5e-5`, same `Nb=120`) without changing the model or
objective. Quantify ranking drift and winner's curse. Do not silently replace
the production loss with the tight loss; report both.

### E. Identification

1. Recompute a weighted moment Jacobian at the final production-box point and
   at the final preferred relaxed/focused point using fresh finite differences
   and strict equilibrium solves.
2. Report singular values, rank at relative tolerances `1e-2`, `1e-3`, and
   `1e-4`, condition number, near-null vectors, and parameter-specific moment
   loadings.
3. Test the prior findings that `psi_child` is near-null and that
   `theta_n` trades against the child-space block. At candidates with
   `theta0≈0`, explicitly assess whether `theta_n` has any economically
   meaningful identification.
4. Separate external/nested restrictions from estimated parameters. Do not call
   the system identified merely because 15 exceeds 14.
5. Explain which variation identifies `h_bar_0`, `h_bar_jump`, `h_bar_n`,
   `chi`, `H0`, fertility tastes, tenure smoothing, and the bequest block.

### F. Fit, attainable set, and economic interpretation

1. Produce the full target-fit table for every finalist: target, model, gap,
   weight, and loss contribution. Recompute totals independently.
2. Produce every free/estimated parameter, production restriction, relaxed
   restriction, normalized bound position, and whether it is near or outside a
   production bound.
3. Explain why loss remains far from zero. Quantify how much is due to old-age
   ownership, aggregate ownership, large-home share, room gaps, fertility, and
   mean rooms.
4. Determine whether the remaining loss reflects incomplete search, numerical
   noise, weak identification, internally inconsistent targets, or a structural
   inability of the model to hit the moments jointly.
5. In particular, diagnose why lowering `h_bar_0` can improve fertility/room
   moments while old-age ownership worsens. Decide whether this is a coherent
   economic trade-off or suspicious implementation behavior.
6. Give a precise verdict on `chi`: previously near a bound, now apparently
   preferred near 1.1. Use the controlled profile to distinguish a changed
   optimum from poor search.
7. Loss zero is not presumed attainable. State what a defensible benchmark is,
   given 15 moments, the weights, sampling uncertainty, grid noise, and only
   approximately 11--12 informative directions.

### G. Tests and diagnostics

1. Run the complete intergenerational unittest suite and report every test,
   skip, failure, and runtime.
2. Run the inert-diff, entrant-feasibility, borrowing-floor, debt-taper,
   down-payment, target/weight, and population-closure regression checks.
3. Regenerate the standard diagnostic packet for the final selected candidate.
   Inspect it visually; do not infer correctness from file existence.
4. Check monotonicity in wealth, choice-probability bounds and aggregation,
   market clearing, boundary mass, renter-cap mass, owner-rung use, lifecycle
   ownership, fertility/child-at-home profiles, and value/policy discontinuities.

## Required verdicts

For each claim use exactly one of:

- `CONFIRMED`
- `CONFIRMED WITH QUALIFICATION`
- `REFUTED`
- `NOT TESTED`

At minimum give explicit verdicts on:

1. objective and parameter contract consistency;
2. feasibility repair correctness in every calibration path;
3. July 12 incumbent reproduction;
4. production-box improvement;
5. relaxed improvement;
6. focused-polish improvement;
7. `h_bar_0` interpretation and admissibility;
8. `chi` not being the useful relaxed margin;
9. bequest-ridge interpretation;
10. search completeness;
11. numerical robustness/tight-evaluator ranking;
12. identification;
13. whether the preferred point is publishable as a production calibration;
14. whether the remaining high loss is primarily a search problem or a model/
    target-system problem.

## Deliverables

Write everything under:
`output/model/full_calibration_audit_20260713_claude/`

Required files:

1. `FULL_CALIBRATION_AUDIT_20260713.md`
   - ten-line executive summary first;
   - verdict table;
   - findings register ordered `FATAL`, `MAJOR`, `MODERATE`, `MINOR`;
   - evidence and reproduction commands for every material claim;
   - clear separation of production, relaxed-frontier, no-bequest, and profile
     estimands.
2. `CLAIM_VERIFICATION.csv`
   - claim, verdict, evidence file, command, numerical result, caveat.
3. `search_inventory.csv`
   - every task/wave/arm/config/count/state/best.
4. `target_fit/`
   - complete target tables for all finalists.
5. `parameters/`
   - complete parameter/restriction tables for all finalists.
6. `reproduction/`
   - fresh solve records, hashes, configs, and loss reconstructions.
7. `identification/`
   - Jacobians, singular values/vectors, rank tables, and interpretation.
8. `diagnostics/`
   - standard visual packet plus a short visual inspection log.
9. `profile/`
   - independently reconstructed `h_bar_0 x chi` table/figure and nested-arm
     comparison.
10. `CODEX_CROSSCHECK_PROMPT.md`
    - a standalone falsification prompt naming the smallest file set and exact
      commands needed for another agent to challenge the audit.

Also provide a concise final chat response containing the ten-line executive
summary and links to the main report, claim table, search inventory, and
cross-check prompt.

## Audit discipline

- Cite file paths and line numbers for code claims.
- Record commands exactly and preserve stdout/stderr under the audit folder.
- Never convert an optimizer failure into evidence of economic infeasibility
  without reading the actual exception/census.
- Never call a surrogate prediction a model result.
- Never infer search completeness from wall time or evaluation count alone.
- Never hide a target miss behind the scalar loss.
- If a requested computation cannot be completed, mark it `NOT TESTED`, explain
  why, and do not weaken the corresponding conclusion silently.
- Do not implement any fixes. End with a ranked recommendation list separating:
  (i) numerical/search follow-up, (ii) identification/target work,
  (iii) model changes, and (iv) paper disclosure.
