# Three entrant-wealth calibration pilots

All three pilots completed with two exact selected full-GE repeats each.
See [final comparison and full tables](RESULTS.md). No production adoption.

The author-approved parallel one-hour experiments are described below. Each array arm uses one CPU,
24 GiB, one computational thread, the 120-asset/9-income grid, the indexed
corrected-credit engine and the full birth-renewal-price/population-housing
closure from the verified `credit053_v2` comparison.

| Arm | Initial wealth | Renter next-period floor |
|---|---|---|
| empirical_credit | Five empirical wealth/income bin means independently times current annual gross income | -0.25 |
| zero_wealth | All entrants at zero | 0 |
| nonnegative_mean | Negative bin means set to zero; positive means multiplied by 0.363238516 to preserve original mean wealth | 0 |

The common annual saving and borrowing rate is 2%; homeowner financing is
unchanged. The .25 allowance adopts the KMV2018 ratio of unsecured capacity to
average annual labor earnings, without their borrowing premium. Censoring the
five-bin approximation is distinct from censoring raw survey observations.
Inputs, tails and projection are reported per arm. No feasibility relocation
is allowed. Necessary input-budget checks are not lifecycle certificates.

`plan.json` embeds all 31 reference parameters/bounds and the complete 14-row
moment/weight/role contract (10 scored rows). H0 and psi_child are fixed: under
this closure H0 scales population rather than per-household moments, while
price clears birth renewal. The other nine historical coordinates, including
tenure-choice dispersion, genuinely bind the native parameter object. Beta is
converted from annual to period units, its discount-rate aliases are updated,
and the fertility-dispersion alias is retained. Native observers verify every
actual parameter; annual-beta round trips alone allow 2e-12 floating tolerance.

`runner.py --mode preflight --arm ARM --out NEW --deadline-seconds 300` runs the
same mocked price and search loops with zero lifecycle claims. `test_pilot.py`
checks entry laws, all nine parameter bindings, all31 actual-parameter drift
rejections, the D=0 price-loop path, baseline failure, failed probes and deficient
Jacobian rank. The preparation tests and all three native baseline and selected-repeat validations pass.

`runner.py --mode run --arm ARM --out NEW --deadline-seconds 3600` is restricted
to single-core Torch Slurm. Source pins cover the reused81-file stage and all
new runtime files. The launcher authenticates the old and new stage inventories.
It never calls the default normalized-population equilibrium CLI, recalibrates
the child-benefit level, revises targets or restores legacy credit tapering.

The exact upper-bound loop is two full baseline GEs, nine one-sided derivative
probes, two damped Gauss–Newton proposals, and two selected full-GE repeats:
15 full GEs, each capped at20 lifecycle calls and300 seconds per lifecycle call.
Every GE itself includes its native selected-price repeat. Each native repeated
GE checks arrays, complete14/31 tables, closure and actual17 standard PNG hashes.
Search begins only after the separately repeated baseline passes all gates.
Accounting, source, parameter or feasibility failures halt the arm; a genuine
unbracketed exploratory root is recorded as an inadmissible proposal, never a
zero derivative. Gauss–Newton is skipped if the baseline Jacobian lacks rank9.

The hard one-hour launcher clock includes preflight, authentication, compilation,
search and reporting. Search preserves a1000-second reserve; final verification
also respects the inherited700-second native-start guard. The15-point loop is
an upper bound, not a promise: incomplete derivatives or missing final repeats
remain explicitly provisional. No automatic retry or extension. Prior observed
cold full workflows took187 seconds at120×9, so15 workflows suggest47 minutes;
changed candidates can cost more. The deadline takes precedence.

Outputs include input_contract.json, initial_entry_distribution.npz, latest.json,
latest_completed.json, best_so_far.json, cases.json, identification.json and
completed.json. Each successful point contains full target/parameter tables and
17 standard diagnostics under `POINT/phase_b_ge/selected_root/`. Best-so-far is
provisional until two final full-GE repeats. The derivative report is at the
baseline, not a fresh identification certificate at the selected point.
No reference promotion, transition or fine-grid convergence follows this pilot.
# Launch — September 30, 17:01 New York

Submitted Torch array **18895422**, tasks0/1/2 for empirical_credit,
zero_wealth and nonnegative_mean. Each task has one CPU,24GiB and a hard
one-hour clock. Lead independently passed8 tests, all3 exact300-second
preflight CLIs,85 runtime/input pins and the8-file overlay archive hashes.
This paragraph records the launch state; terminal evidence is in `RESULTS.md`
and `final_verification.json`. See `launch.json` and `lead_verification.json`.
