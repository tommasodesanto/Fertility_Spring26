**Terminal outcome:** all three corrected jobs stopped because renewal was unbracketed in the inherited price interval. No new equilibrium or fit was obtained. See [full results and tables](RESULTS.md). The original twelve searches are unchanged.

# Fixed-parameter utility diagnostics

Experimental, author-authorized diagnostic; no adoption or recalibration. Existing twelve searches are untouched. All three new arms use the authenticated pilot-selected `nonnegative_mean` nine-coordinate vector, corrected scalar credit and owner sale solvency, mortality, 2% annual interest, D=0, 120×9 grid, five-bin nonnegative mean-preserving entry, fixed H0 and child-benefit scale. Price clears actual birth renewal and population clears physical housing supply. Target/weight contract is unchanged (14 rows, 10 scored).

- `floor`: historical parenthood-only physical housing requirement h_P=1.8900476600128304 rooms, constant childless consumption share, compensation off, nonlinear child benefit retained.
- `no_A`: compensation off, inherited first-child expenditure-share change retained, no floor, nonlinear child benefit retained.
- `constant_alpha`: compensation off, constant childless consumption share, no floor, nonlinear child benefit retained.

Children currently at home determine the floor, equivalence scale and benefit. The floor has no subsequent-child slope. Renters use residual physical rooms h−h_P; owners use χ(h−h_P), retaining the native owner service premium after the physical-room subtraction. This is the historical convention; it is not χh−h_P. The reference rent is retained as an inactive input with compensation off. Owner products with h≤h_P are infeasible under the native strict floor gate.

`authenticated_control.json` pins the common selected control, its full 14-target/31-parameter tables, closure and 17 actual standard PNG hashes. Both independent full-GE repeats have identical tables/closure and native recorded PNG hashes. Repeat PNG bytes were not separately downloaded. Control PNGs remain in the original pilot selected report; compact control tables are copied into `reference_control/`.

Runner changes only process-local input bindings and the observer's reported h_P, which the inherited observer otherwise hardcodes to zero. Model engine sources, numerical gates and objective are unchanged. The native reporter check uses its own authentication directory. Each arm requires a native full GE and an independent full-GE repeat, in addition to each GE's selected-price repeat. At most two full GEs per arm; one CPU and 24 GiB each. Deadline is 2026-09-30 23:41 UTC (epoch1790811660), including authentication and preflight. No retries or optimizer. Credible one-round estimation requires a larger explicit budget.

Preparation: all three exact shared-utility array and numerical-formula checks pass locally with zero lifecycle solves. Local native authentication is blocked by pre-existing drift in active `code/model/tools/e5f_exact_policy_cache.py` relative to the frozen runtime contract; Torch's pinned frozen project must pass the exact native reporter and mocked GE preflight before any solve. Runtime results are not yet available.

Invocation: `runner.py --mode preflight|run --arm floor|no_A|constant_alpha --out NEW_DIRECTORY --deadline-seconds REMAINING_SECONDS`. `deployment/launch_torch.sh` dispatches the three arms with the common deadline. The lead must review and authorize submission.

## Cluster launch receipt

Array18901961 launched three fixed-parameter diagnostics (floor, no_A, constant_alpha), each one CPU/24 GiB and two full GE evaluations including a repeat, with common hard stop2026-09-30T23:41:00Z. Existing twelve calibration jobs were preserved. All three exact CLI preflights passed with zero lifecycle solves. Native launch failed before solving in floor/no_A with `Fresh interpreter required: model already imported`; no automatic retry was performed. See `launch.json` and `deployment/cluster_receipts/` for the recorded evidence. These are diagnostics, not refits.

## Authentication sequencing repair

Initial array18901961 passed all three native reporter/mocked-loop preflights, then stopped native execution before any lifecycle solve because the same interpreter authenticated twice. The revised runner retains full reporter authentication in the separate preflight process. Run mode installs only the process-local reporting adapter, then authenticates exactly once through the inherited native evaluator. That evaluator still checks every effective parameter before its first lifecycle call. Economics, inputs, gates and the 23:41 UTC deadline are unchanged; original failed receipts are preserved.

### Authorized orchestration repair

The first array failed before any native lifecycle solve because its zero-solve verification imported the model before the native loader required a fresh interpreter. Its archive and receipts are retained in `deployment/attempt1/` and the original remote path. Reviewed attempt2, array18902151 at `/scratch/td2248/projects/parenthood_floor_quick_v2`, separates those two authentication paths. It retains the same three economic variants and hard23:41Z deadline; source and solver checks are unchanged. This retry is an orchestration correction, not an economic change.
