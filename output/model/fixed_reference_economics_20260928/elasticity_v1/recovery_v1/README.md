# Bounded +1% recovery, September 29

**Production submitted: Torch job 18815133**, pending priority at 15:10 EDT on
September 29. The lead reviewed the unchanged scientific block, final saving
logic and source identities and independently verified final smoke 18814979
completed with exit 0:0. No production solve has started as of this check.
Do not resubmit. Its clock begins on entering `launch.sh`; `solve_v1/launch.json`
will record the fixed 2100-second deadline. Hourly monitoring is active.
Zero-solve rendering job **18815186** has dependency `afterok:18815133` and
automatic cancellation if that dependency fails. The renderer additionally
requires the passed three-price receipt and authenticates both input CSVs.
See `../../slide_inputs_v1/recovery_render_v1/README.md` for the figure/table
paths. Preparation statements below record the worker's handback before the
lead's production submission.

This is a distinct, author-approved recovery contract for the frozen **2007 stationary reference — block0506, September 28 verified export** (manifest SHA-256 `147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4`). It does not reopen the old 5400-second deadline or turn the failed `reference_1010` artifact into a passed case. The prior v2/v4 history has seven lifecycle attempts, six certified passes and one failed +1% reference attempt. The truncated failed checkpoint is excluded.

The new run may make **two** fresh lifecycle attempts, `reference_1010` then `credit_1010`, on the identical 602-node common union grid at prescribed factor 1.01. The reference run retains original borrowing limits; the credit run removes artificial renter, purchaser and incumbent-owner limits subject to natural lifetime repayment and net-estate solvency. The only price change is +1% in house price and its mapped rent; preferences, earnings, entry/inherited distributions, transfers, fiscal objects and supply remain frozen. There is no fertility normalization and no general-equilibrium certification. Every case retains all applicable household/cohort/impact/credit/estate gates, a full 14-row target fit, 31-parameter table, 17 standard plots and the complete checkpoint packet.

The first two +1% cases are **new cases**, with 900 seconds each including reporting, a single attempt each, and 2100 seconds total from entry into `launch.sh`. One CPU and 16 GiB are requested. The six earlier passed cases (four q0 controls and two −1% cases) are read-only comparison inputs, pinned by receipt, checkpoint, fit table, parameter table and all 17 plot SHA-256 hashes. Their 602-node q0 fresh repeats had exact array/table/plot checks. The original v2 plan SHA-256 is `6354bafd069fbfc4eda21271fc91747a22e3aeb8c90a57293f7dfde01badf28f`.

The recovery driver copies the v2 scientific solve/gate segment byte-for-byte (61 lines from `rt = prepared.rt` through `impact_summary = credit.aggregates(...)`). It adds a pre-checkpoint identity check, writes the **full** checkpoint through a `.tmp` gzip file, closes and `fsync`s it, then atomically renames it before standard plots. The final path cannot contain an interrupted pickle stream. A durable `checkpoint_saved.json` records its hash and time even if later plotting fails. The source/input identity check runs again after plots and the case is passed only when all 17 plots and its receipt exist. A saved checkpoint alone is not a passed case. No duplicate checkpoint write occurs.

The resulting comparison covers factors 0.99, 1.00 and 1.01 in both regimes. `comparison.csv` has `regime,scope,price_factor,outcome,value,unit,prescribed_price,mapped_rent`; scope is `impact` for the fixed inherited distribution and `cohort` for recomputed normalized-entry composition. `elasticities.csv` has `regime,scope,outcome,step,central_log_elasticity,lower_one_sided_log_elasticity,upper_one_sided_log_elasticity`, with `step=0.01`. The central measure is `[ln Y(1.01q0) − ln Y(.99q0)]/[ln(1.01) − ln(.99)]`. Completed fertility and mean first-birth age are cohort measures only. No ±2% outcome or five-price completion is claimed. Successful completion has `status="passed"`, `complete_three_price=true`, `complete_five_price=false`, the factors and step sizes, both CSV hashes, two new solves, six reused cases and one historical failed attempt.

The staged Torch root is `/scratch/td2248/projects/fixed_reference_elasticity_recovery_20260929`. Final production plan SHA-256: `b12c59a89b1259494de94c071b2e1176b5e22afaf8503ea36d7e47b3b30c0962`; driver SHA-256: `fafec584517acdbe0872296711cfe2b5be11a242d2b2de0179e3720e055ab012`; launcher SHA-256: `4ee6eb93d17e6d33f299a1bb7446caf7847b5ab3a425998ce82135de9886da52`. The first zero-solve smoke job, `18814847`, failed **only** in its final test-status print after its assertions passed: the synthetic temporary plan had been deleted. Its source, plan and logs are preserved remotely as `source_failed_smoke_v1`, `plan_failed_smoke_v1.json`, `smoke_v1.log` and `smoke_v1.err`. Corrected zero-solve smoke job `18814888` **completed successfully in seven seconds**. Its passing source, plan and logs are preserved as `source_passed_smoke_v2`, `plan_passed_smoke_v2.json`, `smoke_v2.log` and `smoke_v2.err`. Final zero-solve smoke job `18814979` **completed successfully in nine seconds** with the final source/plan pins. It authenticated the prior six cases and tested exact two-case controller flow, no-retry failed child, budget and case-timeout branches, full atomic checkpoint and interrupted atomic write. It called no lifecycle solver. The final staged source and plan are read-only; production has not been submitted.

After lead review, the **only** production submission command is:

```bash
ssh -o BatchMode=yes torch 'sbatch /scratch/td2248/projects/fixed_reference_elasticity_recovery_20260929/source/launch.sh'
```

This worker has not submitted that production job. The launcher writes `recovery.log` and `recovery.err`; each case writes its own log, progress, full checkpoint, tables, plots, gates and passed receipt. `latest_completed.json`, `best_so_far.json` and a final `completed.json` are in `solve_v1`. A failure leaves an atomic `failure.json` and no automatic retry.
