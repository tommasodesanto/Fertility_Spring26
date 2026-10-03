# Independent evidence check of Fable memo (October 3, 2026)

Scope: read only, no model solve. Reviewed `fable_analysis/ECONOMIC_MEMO.md`,
`analysis/native_grid_analysis.py`, saved JSON, the exact-repeat arrays for
original chain 15 and alternative chain 13, the executed timing engine, and the
historical hard/quarter engine. All paths below are relative to the repository
root. This receipt does not adopt a calibration or revise the memo.

## Grounded results

- The saved-array identities are sound for the three tested distributions.
  `analysis/native_grid_analysis.py:79-93` compares child-count shares from
  `g`, `g_beginning_distribution`, and
  `g_cross_sectional_wealth_distribution` with the observer's pre/post
  accounting. Maximum post differences are below 2e-15; maximum pre
  differences are about 0.0165. The NPZ lists no pre-fertility distribution.
  This conclusion applies to these selected saved objects, not every possible
  engine output. Source `code/model/experiments/purchase_timing_sandbox/source/refactor_lab/engine/distribution.py:704-735`
  subtracts realized first births before saving `g`; `:1293-1295` copies that
  post-fertility `g` into both named distributions.
- For the current selected points, reversing first-birth selection with
  \(g_{0,\mathrm{post}}/(1-\pi_j a_j)\) reproduces the observer's pre-fertility
  childless mass to about 1e-15 at ages 18, 22, and 26. The code and source
  are `analysis/native_grid_analysis.py:149-160` and the engine lines above.
  Conditioning the policy slices on childless inherited renters and weighting
  by this reconstructed mass is appropriate. The feasibility groups also
  differ in income, so their average attempt rates are descriptive.
- The age-25 observer is reproduced exactly: original 0.5330788222,
  alternative 0.5338046821 (`analysis/native_grid_analysis.py:95-110`; saved
  `analysis/out/native_grid_analysis.json`). In the original point, any birth
  is 0.441959 versus CPS 0.457254; children per mother are 1.20617 versus
  1.77041; two or more are 0.091120 versus data bounds 0.176137–0.323110.
  The 0.875 post-cell weight corresponds to age 25.5 within the model's
  [22,26) cell. CPS completed age 25 means interviews in [25,26); this is an
  approximation based on uniform within-cell timing.
- A separate zero-solve composition check across *all* inherited tenures
  confirms the model's upper-income concentration: original-arm first-birth
  flow shares in income states 6–9 are 88.8% at ages 18–21 and 65.8% at
  ages 22–25; states 1–4 contribute 0.02% and 0.89%. The corresponding
  alternative-arm shares are 90.3% and 66.5%. This supports the model-side
  descriptive statement more directly than the inherited-renter-only slices.
  The contrasting empirical claim about US income groups was not rechecked.
- All 14 snapshot hashes in `input_snapshot/manifest.json` match both their
  snapshot and their currently named source file. The two 51 MB NPZ arrays
  are selected by chain-specific paths in `analysis/native_grid_analysis.py:25-28`
  but are not in that manifest. Current SHA-256: original
  `2170a4bd2e1ced398fa9000cb96d0d2129b61c2ad8324b9da748874b5977fc80`,
  alternative `029c85e0b2d7bc4f7d4680713807f4939bf183b231d26931b040fd67d8302168`.

## Corrections needed before relying on the memo

1. **Closing-screen formula and shares.** `analysis/native_grid_analysis.py:179`
   uses \(b+y/R\ge(1-\phi)Q\). The executed soft engine sets
   `dp_choice=(dp_arr-income_for_purchase)/Rg` with
   `dp_arr=(1-phi)*Q` and then compares `bg_b<dpn`
   (`code/model/experiments/purchase_timing_sandbox/source/refactor_lab/engine/household.py:188-190,894-901`;
   `engine/kernels.py:315-329`). Thus its screen is
   \(Rb+y\ge(1-\phi)Q\), equivalently
   \(b+y/R\ge(1-\phi)Q/R\). Zero-solve recomputation from the same arrays:
   original pass shares ages 18/22/26 are 0.971/0.951/0.992 instead of
   0.971/0.951/0.933; alternative are 0.971/0.952/0.933 instead of
   0.965/0.949/0.933. The ending-floor inequalities in script lines 175–178
   are consistent with the stated budgets. Revise memo lines 99–125 and figure
   references after correcting the script.
2. **Historical closing tests.** Memo lines 204–207 say hard/quarter rules only
   tighten the ending floor and that no closing-binding constraint was tested.
   The historical plan explicitly says hard closing resources
   \(A\ge(1-\phi)Q\) and quarter
   \(A+0.25S\ge(1-\phi)Q\)
   (`purchase_rules_overnight_v1/local_runtime/local_plan.json:148-150,3075-3077`).
   The historical engine rejects renter buyers below its `dp_choice`
   (`purchase_rules_overnight_v1/engines/quarter/refactor_lab/engine/kernels.py:314-328`),
   and adds a saving floor for the quarter arm (`:912-918`). With native
   purchase income, `household.py:905-912` uses current income in closing
   resources. A strict *beginning liquid wealth alone* test may still be
   untested, but the memo must acknowledge the historical hard/quarter closing
   tests and their older contracts.
3. **Conditional, not universal, fertility ceiling.** The 0.676 original
   ceiling in `analysis/native_grid_analysis.py:123-131` fixes the observed
   18–21 first-birth flow and the 22–25 first-birth flow, while setting second
   births among early mothers to certainty. It is not a bound on all feasible
   parameters or policies. The calculated 0.349 early-first-birth share needed
   to reach 0.810 holds the 22–25 first-birth hazard fixed; it implies an
   age-25 any-birth share about 0.504, above CPS 0.457. If the empirical
   any-birth share were held fixed, the required two-child share would be
   0.352, implying an early-first-birth share at least 0.403 with uniform
   interpolation and certain second births. These are accounting scenarios,
   not a proof the target is unreachable. The mean-first-birth-age fit (memo
   line 165) does not forbid earlier births because later births could offset
   the mean. Figure F1's title at script line 268 must say the displayed
   ceiling conditions on current first-birth timing/flows.
4. **Causal and identification language.** Memo lines 29–36, 124–136,
   177–196, and 215–220 use “not a credit phenomenon at all,” “binds nobody,”
   “income, not credit,” “not identified by the financial block,” and “null is
   structural.” The selected slices show weak first-birth movement under the
   compared financial/timing contracts and a strong income gradient; they do
   not establish zero credit effects in all specifications, a causal income
   mechanism, or rank identification. The memo's own original-arm table has
   14–24% of young childless renters failing its ending-floor feasibility
   screen, and 4–9% of young realized owners at the debt floor. Call the
   evidence local/descriptive and state the contracts. An owner policy at the
   floor is not by itself proof that the floor causes the chosen amount.
5. **Other literal corrections.** Memo line 59 calls temporary financing paths
   “tax-shock”; verify and label the shock correctly. Line 194 calls a
   6.927-to-4.458 wealth target change a “halving”; it is a 35.6% reduction.
   Line 250 extends an earlier grid-refinement diagnostic to the new winners
   without a matched test. Lines 131–132 say owner floor contact is under 2%
   after age 38; the JSON shows 30.5% original and 70.9% alternative at the
   terminal age-82 cell (ages 38–78 are under 2%). “Near-bound” fertility
   scales alone do not prove probabilities are nearly deterministic; the
   observed policy gradient supports a narrower descriptive statement.

## Suggested restrained conclusion

At the two selected soft timing points, the age-25 children-ever-born miss is
mainly the number of children among mothers, while the share with any birth is
close to the data. Held-coordinate timing and older financing comparisons show
small first-birth responses under their particular contracts despite sizable
ownership responses. The evidence does not yet establish that borrowing never
matters for fertility, that the empirical age-25 target is infeasible, or that
the financial parameter block lacks identification. The hard and quarter
experiments test stricter closing rules at older calibrated points; a matched
policy comparison under the current accepted soft contract remains outstanding.
