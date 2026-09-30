# Proposed 120 wealth × 9 income numerical comparison

**Numerical follow-up:** Torch job 18879780 completed the 160×15 full-GE control in 477.756 seconds, but rejected the first 120×9 fixed-price evaluation because the inherited feasibility routine relocated occupied entrant mass to higher wealth. The external guard rejected this change; there is no accepted 120×9 equilibrium, speed comparison or calibration. The failed proposal is preserved without retries or credit changes. See [runner evidence](runner/README.md).

The preparation record below describes input construction, which did not relocate entry mass. This must be distinguished from the subsequent forward-distribution routine. The proposed discretization is experimental and has not replaced the reference.

## Verified facts

- Reference: **2007 stationary reference — block0506, September 28 verified
  export**. Bundle SHA256 `427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7`;
  checkpoint `b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d`.
- Effective reference dimensions are 160 wealth nodes and 15 income states,
  authenticated from arrays and `make_grid`, rather than constructor defaults.
  Income is a **single** 15-state Rouwenhorst chain. Historical 3×5 income
  grouping arrays are inactive (`permanent_income_levels_enabled=False`).
- Period persistence is 0.7345934905942886 and stationary log-income standard
  deviation 0.7130810449881093. The period is four years. The implied annual
  persistence and innovation standard deviation are 0.9257884726460368 and
  0.2695745373962748; these are recovered from authenticated inputs, not newly
  estimated values. Reconstructing the reference chain matches income nodes
  within 1.78e-15 and its transition within floating precision.
- Entry uses the fixed 160×15 conditional matrix, with 169 occupied cells and
  50 occupied wealth nodes. Historical PSID ratio atoms do not override it.
  The flag `entry_wealth_censor_to_frontier=True` is retained. Although fixed-entry input construction bypasses its historical ratio-based branch, the forward-distribution routine also uses this flag. That later routine attempted wealth relocation in the 120×9 trial, causing rejection. The earlier description of the flag as dormant throughout the solve was incorrect.
- The displayed contract has 14 target rows: ten positively weighted rows,
  three zero-weight validation rows and completed-fertility normalization 2.1.
  The 31 displayed parameter rows include ten searched coordinates, separately
  normalized `psi_child`, and externally fixed/inherited restrictions; they are
  not 31 free parameters. Target/weight fingerprint:
  `20e531855075b807886da8e936930dca9090b989b2ce381de39042cf206efdbd`.
  [All reference target fits](reference_target_fit.csv) and
  [all reference estimates and restrictions](reference_parameters.csv) are
  copied byte-for-byte, explicitly as reference evidence, not new-grid results.

## Proposed immutable construction

The [one preparation driver](prepare.py) authenticates the frozen bundle, then
writes a separately identified [proposed bundle](proposed_120x9/bundle.json).
It has a new schema and cannot pass the reference loader as an exact reference.
Its SHA is recorded independently in `preflight.json`. It is not a runnable
production input or an adopted entry specification.

The wealth grid is a deterministic 120-node subset of the original 160 nodes.
All 50 occupied entry wealth atoms, zero and both original endpoints (-12 and
3000) are compulsory; remaining nodes fill the largest index gaps. Thus entry
wealth itself is not interpolated, truncated, forgiven or censored. This is a
transparent numerical grid selection, not an optimal grid claim.

The nine-state income chain uses the same period persistence and stationary
log variance, with mean income multiplier normalized to one. Its discrete
tail support and higher moments differ, as with any smaller Rouwenhorst
approximation. Cumulative-probability-bin overlap transports the old
conditional wealth distributions into the nine new income bins. Every old
income bin's mass is assigned by overlap; no mass is deleted and no node is
clipped. This preserves the entire original wealth marginal and new stationary
income weights. It approximates the original rank association; it does not
preserve the exact original discrete joint wealth-income law.

| Joint entry statistic | Reference 160×15 | Proposed 120×9 |
|---|---:|---:|
| Total mass | 1 | 1 |
| Mean wealth | 0.1865196792468183 | 0.1865196792468183 |
| Wealth second moment | 2.616827638490699 | 2.616827638490699 |
| Negative-wealth mass | 0.2625085160126798 | 0.2625085160126798 |
| Mean income multiplier | 1 | 1 |
| Income second moment | 1.628765744226251 | 1.605446747506460 |
| Wealth-income covariance | 0.1193721443144067 | 0.1122969087784427 |

Covariance falls by about 5.93%. This is visible approximation error requiring
lead review before a resolution experiment, not a claim of exact entry-law
preservation. The largest original wealth-atom marginal mass difference is
6.94e-18. Simple nonnegative linear interpolation in income would fail for
0.18310546875% of the old income mass outside the smaller income support; this
alternative was diagnosed and not applied. Two inactive historical 3×5 index
arrays are omitted in the proposal, with the permanent-group feature still
disabled. All other primitive fields remain unchanged.

Credit must be selected explicitly for a matched pair: either the original
reference contract (scalar unset) or corrected diagnostic scalar D=0.14 in both
arms. D=0.14 is experimental, not empirically estimated or adopted. This packet
does not choose a new magnitude or resolve the upstream borrowing decision.

## Exact integration points and remaining work

For the initial grid comparison, reuse the indexed matched-credit route in
`../small_credit_replication_v1/arms/indexed/`: `single_price.solve_fixed_price`
and `phase_b_ge.run_phase_b`. The latter solves prices for actual birth renewal,
then population for absolute housing supply. Default normalized-population
`refactor_lab` GE must not replace it. Hold reference `psi_child` and every
economic primitive fixed across dimensions; do not normalize it to offset
grid error. Original-credit and corrected-credit pairs are separate comparisons.

Minimal isolated adapter changes needed before an executable loop smoke:

1. Replace `driver.context_from_bundle` only in an isolated adapter with a
   proposal loader pinned to `preflight.json`, checking every unchanged field,
   only the declared changed fields and effective array dimensions. The
   reference control continues through the exact frozen loader.
2. First run `authenticate_frozen_observers` with the untouched reference
   context, including its exact reference wealth-grid equality gate. Then
   explicitly replace live context P/grid with the authenticated proposal,
   recording both identities and changes. `phase_b_ge.observe_price` already
   supplies live P/grid to native calendar, account, observer and diagnostic
   calls; verify these shape-general paths with fixtures. A wholesale observer
   rewrite or weakened reference authentication is not needed. Do not
   impersonate the frozen checkpoint.
3. `phase_b_ge.observe_price` currently requires a float scalar D and audits
   corrected-credit renter floors. For an original-contract pair, use a
   separately explicit original-rule audit/observer branch; do not turn unset
   credit into D=0. Corrected D=.14 can reuse the existing audit.
4. Preserve all 14 fit rows, all 31 estimates/restrictions, all 17 named
   diagnostics, raw entry-feasibility, estate, PAYGO and calendar checks, actual
   renewal price residual and population-scaled housing closure. Saved native
   array comparisons across grids require moment/policy interpolation receipts,
   not array equality. Exact repeats within each grid still require equality.
5. Run the actual controller's synthetic exact-loop smoke and then one fresh
   paired GE plus repeat only after review. This packet has zero-solve input
   tests; it **has not passed an exact-loop smoke**.

For later full recalibration, the existing source-pinned controller is
`output/model/fertility_identification_20260928/two_stream_overnight_v1/search.py`.
Its worker invokes the native `e5f_evening_calibration_runtime.setup` factory
through the driver lineage. `EveningObjective.bind` copies prepared parameters;
`StationaryObjective.normalize` passes `selected['b_grid']` and
`self.rt['model']` to `solve_balanced_initial_equilibrium`, separately adjusting
`psi_child` to fertility 2.1. It does **not** currently consume this typed input
bundle or the refactor engine. A separately pinned factory adapter must bind
the proposed parameters/grid and provide the cheaper engine's compatible
entry point; copied source/config and fresh full-evaluator scientific replay
are needed. No direct import substitution is certified here. This normalization
and market-closure contract differs from the fixed-psi renewal-price/population
comparison, and must be reconciled before production calibration.

The ten coordinates and bounds remain the existing contract: H0 [.2,80],
annual beta [.94,.99], chi [.1,5], first-birth fixed cost [0,8], first- and
later-birth choice scales [.02,50], theta0 [0,8], housing-share jump [0,.25],
child-benefit curvature [0,.8] and tenure-choice scale [.001,.1]. Ten scored
moments plus the normalization equality provide the count restriction; local
numerical rank and statistical identification at new-grid estimates are
unverified. Do not drop or reweight targets to obtain a run. Full recalibration
also remains dependent on the upstream conceptual borrowing decision.

## Verification, budgets and smallest next computation

The preparation script passed on one CPU in 0.31 seconds with zero solves.
Checks cover source input authentication, 160×15 reconstruction, effective
120×9 round-trip dimensions, strictly ordered wealth nodes, all occupied entry
atoms retained, conditional mass sums, valid entry indices, nonnegative
stochastic transition rows and invariant distribution. Reference bundle hash
remains unchanged. Regenerate with a NumPy-enabled Python:

Final probability-contraction check passed in 0.245 seconds with Python
`-W error`, with no warnings. The tiny conditional-mass contraction uses
`np.einsum(..., optimize=False)` and agrees exactly with independent ordered
broadcast summation; all entries are finite and nonnegative. CDF input weights
are finite, positive and sum to one. An earlier Apple BLAS matrix-product path
emitted floating-point warnings despite finite results; this path was replaced,
not its warnings suppressed. Final wealth-atom mass error remains 6.94e-18.
The proposed bundle is now frozen; its final pins are in `preflight.json`.
Regeneration creates a new NPZ container hash and requires a new proposal pin;
do not regenerate underneath an authenticated runner.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
python3 output/model/publication_refactor_20260929/grid_resolution_v1/prepare.py
```

The smallest next computation is the isolated adapter's mocked **exact**
controller loop, zero lifecycle solves, with a 300-second cap. After that,
proposed numerical budget is one 160×15 GE, one 120×9 GE and one exact 120×9
repeat: at most 18 lifecycle evaluations, 300 seconds per evaluation and 2400
seconds total, one CPU. No credit-cap search. Checkpoint, heartbeat and latest
completed/best summary must appear after every price trial. Stop on source or
target drift, entry/estate/fiscal/closure failure, missing outputs, or any cap;
no automatic extension or retry. These are proposed budgets, not launch approval.

Matched indexed 160×15 workflow previously took 417 seconds. State count falls
to 45% and a wealth-quadratic work proxy to 33.75%; actual timing is unknown.
Approximately 7–10 minutes per full workflow is a conservative preparation
estimate, not a speedup claim. The 40-minute envelope is finite and may yield
an incomplete comparison if any stage reaches its cap.
