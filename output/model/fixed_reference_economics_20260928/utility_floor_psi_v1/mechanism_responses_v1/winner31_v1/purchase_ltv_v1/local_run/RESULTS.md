# Mortgage financing experiments: fixed prices and stationary GE

**Earlier fixed-price status:** completed preliminary comparison at the saved reference price. Those cells passed household, purchase-accounting, estate, fiscal, and policy-array gates. They are prescribed-price diagnostics, not general-equilibrium roots; later GE results are separate below. Production adoption is false.

The accepted baseline is the original verified 31.284 point with (q_0=0.719168368828958). Both 90% cases reuse the same saved baseline PRE distribution (SHA-256 `459ab9229e7ce58376bda8c933b583ec059104ede31c57627ce7368c41e5a51b`). Mortgage cases are purchase LTV .90 / stayer collateral share .80, then .90 / .90. Original targets, weights, timing, earnings, entry, grids, and household preferences remain fixed. The .90/.80 and .90/.90 results each have 14 target rows, 31 parameter rows, and the unchanged set of 17 standard diagnostic plots.

| Measure | Baseline .80/.80 | Purchase .90 / stayer .80 | Both .90/.90 |
|---|---:|---:|---:|
| Births per household, common baseline PRE | 0.115253846 | 0.115237704 (−0.000016142) | 0.115266648 (+0.000012802) |
| First-birth flow, common PRE | 0.049912693 | 0.049894082 (−0.000018611) | 0.049905588 (−0.000007105) |
| Second-birth flow, common PRE | 0.041457721 | 0.041459271 (+0.000001550) | 0.041470373 (+0.000012653) |
| Third-bin entry flow, common PRE | 0.023883432 | 0.023884351 (+0.000000919) | 0.023890686 (+0.000007254) |
| Ownership, common PRE | 0.659373600 | 0.664250109 (+0.004876509) | 0.675949055 (+0.016575455) |
| Rooms per household, common PRE | 5.977125401 | 5.983666917 (+0.006541516) | 6.005350742 (+0.028225341) |
| Completed fertility, saved cohort summary | 2.099999968 | 2.097688986 (−0.002310982) | 2.097382457 (−0.002617510) |
| Childless share, CPS ages 40–44 observer | 0.202532984 | 0.203299541 (+0.000766557) | 0.203454140 (+0.000921156) |
| Mean first-birth age | 25.956378034 | 25.959859413 (+0.003481379) | 25.959080785 (+0.002702752) |
| Share of first births at ages 30+ | 0.229310813 | 0.229771336 (+0.000460523) | 0.229791268 (+0.000480455) |
| Mean children ever born, capped at 3, age 25 | 0.531217550 | 0.530581197 (−0.000636354) | 0.530641300 (−0.000576250) |

Childlessness, first-birth age and the age-30+ share are the saved `uniform_birth_time` observer measures; completed fertility is the saved cohort-summary measure. The target-fit CSVs contain all 14 target/model/gap/weight/contribution rows, with exact target, weight, role, and row-order identity across the three cells. The 31-row parameter CSVs are retained for all cells; the only differing parameter row is `financed_share` (`.80` in baseline and purchase-only, `.90` in both); purchase origination LTV is recorded separately in each closure receipt.

## Full packets

- Baseline: [14-row target fit](retry5/results/baseline_80_80/target_fit.csv), [31-row parameters](retry5/results/baseline_80_80/parameters.csv), [17 diagnostic plots and summary](retry5/results/baseline_80_80/standard_diagnostics/summary.json), [closure](retry5/results/baseline_80_80/closure.json), [numeric acceptance receipt](retry5/results/baseline_80_80/numeric_acceptance.json).
- Purchase LTV .90 / stayer .80: [14-row target fit](retry7/results/purchase_90_stayer_80/target_fit.csv), [31-row parameters](retry7/results/purchase_90_stayer_80/parameters.csv), [17 diagnostic plots and summary](retry7/results/purchase_90_stayer_80/standard_diagnostics/summary.json), [closure](retry7/results/purchase_90_stayer_80/closure.json), [gate ledger](retry7/results/purchase_90_stayer_80/gates.json).
- Both .90/.90: [14-row target fit](retry7/results/both_90_90/target_fit.csv), [31-row parameters](retry7/results/both_90_90/parameters.csv), [17 diagnostic plots and summary](retry7/results/both_90_90/standard_diagnostics/summary.json), [closure](retry7/results/both_90_90/closure.json), [gate ledger](retry7/results/both_90_90/gates.json).
- Search/run record: [completed two-cell receipt](retry7/results/completed.json), [latest case record](retry7/results/latest_completed.json), [local source receipt](retry7/results/source_receipt.json), [runner log](retry7/runner.log).

## Source identity and limits

The immutable winner source-binding SHA-256 is `47b809b4c1d2a18373c2828b2d90a0b23fed498ea64a925eea8b935ba5dc2c93` ([complete binding](../../source_binding.json)); target-contract SHA-256 is `db60605ef444b747b2ed7f9482b6b8c90ef5e8b1d75c49a3ac1ff6f4275c9ba1`. The local execution receipt pins the household engine `2a34f5f28c0c63ca3d24d9e759cd1b5b92aadaddb5fe04795e65abbc0b448082`, household override `ecd712fe036324acfc3020bd4d9d77d4f250115d214ec29b5889f1241d587ae8`, original purchase audit `1289173378b978a51e10f85fffaba34b5840bbda6f3d02c9e1c959e7f318f161`, and path-adjusted local audit wrapper `7f30eaa55839f8d941a604331a7e3e9bd96031d57b54f26adcc3527e93cac1b3`; that wrapper changes only the inherited audit import path and retains its strict inherited-contract SHA check (`71eebc1d42e28ccd304f921535f0d154f3a994451b5c53bcb8097233ce18669a`). The executed local driver SHA-256 is `82ca484a050e27f5d8f86fd23ab46a8954ea86c95fad1bc7acb02dc9ea0e6549`.

The earlier 100% purchase-LTV / .80-stayer attempt failed the production negative-estate gate and is preserved at [its failure receipt](retry6/results/failure.json); it is **not** a valid 100% estimate and is not evidence of global infeasibility. The subsequently identified buyer death-estate-floor implementation gap is documented separately; the .90 results do not resolve it or imply general-equilibrium effects.

Lifecycle accounting across retries: five case-level lifecycle starts are recorded (two baseline attempts, one 100%/.80 attempt, and these two .90 attempts); three complete gated result packets are retained (accepted baseline and both .90 cases). The missing-audit baseline attempt and failed 100%/.80 attempt did not produce accepted case packets. Four earlier failures occurred before a lifecycle start. No 100% retry or baseline repeat was run in this two-cell pass.

## Local stationary GE mortgage-LTV follow-up

This bounded one-core exercise solves price and population closure at the unchanged 31.284 winner point. It is a comparative-static experiment, not a recalibration, global identification result, or adoption proposal. The four cases vary buyer origination LTV and owner-stayer collateral share. A reviewed buyer death-estate floor correction is active. Earnings, entry, preferences, target contract, and the 120×9 grid remain fixed. Every selected point retains the full 14-moment fit, 31-parameter table, and standard 17 diagnostic plots.

The existing closure chooses population (N=S/D) and solves the adjusted-birth/entry renewal condition. The comparator is the saved .80/.80 candidate at q0=0.71916837. It is a verified GE after population scaling: the canonical closure reports N=0.9222615666500361, zero absolute housing residual, and renewal residual −1.53366778166e−8. Its 0.08429 relative residual is the normalized cohort’s unscaled housing-market residual and does not invalidate the population-scaled GE. Changes below are against its stationary q0 cohort and are not fixed-PRE impact estimates.

| Buyer / stayer LTV | GE q | Rent | N | Renewal residual | PAYGO residual | Ownership (Δ q0) | Births/HH (Δ q0) | Rooms/HH (Δ q0) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| .90 / .80 | .71770027 | .12938657 | .91912454 | 1.41e−11 | 3.34e−14 | .6678165 (+.0084429) | .115251085 (−.000002761) | 5.9898095 (+.0126841) |
| .90 / .90 | .71751282 | .12935278 | .91138002 | 8.55e−10 | 2.84e−14 | .6896946 (+.0303210) | .115249252 (−.000004594) | 6.0397144 (+.0625890) |
| 1.00 / .80 | .71617704 | .12911196 | .91751919 | 1.11e−11 | 3.42e−14 | .6804383 (+.0210647) | .115242202 (−.000011644) | 5.9922636 (+.0151382) |
| 1.00 / 1.00 | .71226701 | .12840707 | .90601584 | −4.63e−08 | 2.38e−14 | .7326132 (+.0732396) | .115225513 (−.000028333) | 6.0474517 (+.0703263) |

All selected roots passed the unchanged household, purchase, estate, and fiscal gates; housing clears under N, renewal is within 1e−6, and the native scaled 16/20 entry-queue/distribution check passed. Selected and same-point-repeat tables match for all four cases: exact target/weight/role/order, model/gap/loss within 2e−12, and all 31 estimates within 2e−12. Completed fertility is near 2.1 because replacement is imposed by the renewal closure; it is not evidence of an independently recovered fertility response.

The frozen q0 PRE remains unchanged (SHA-256 `459ab9229e7ce58376bda8c933b583ec059104ede31c57627ce7368c41e5a51b`). Its optional impact at a new price is unavailable: the first .95-price evaluation found 6.45718872382e−12 positive mass at infeasible inherited states under the exact zero-projection gate. No mass was projected or replaced. Earlier fixed-price/common-PRE results remain separate from these stationary GE cohorts.

| Case | GE closure | Full 14-row fit | Full 31-row parameters | Standard 17 plots | Same-point repeat fit / parameters |
|---|---|---|---|---|---|
| .90 / .80 | [closure](ge_retry4/purchase_90_stayer_80/root_6_0p997959/closure.json) | [fit](ge_retry4/purchase_90_stayer_80/root_6_0p997959/target_fit.csv) | [parameters](ge_retry4/purchase_90_stayer_80/root_6_0p997959/parameters.csv) | [summary](ge_retry4/purchase_90_stayer_80/root_6_0p997959/standard_diagnostics/summary.json) | [fit](ge_retry4/purchase_90_stayer_80/selected_repeat_7_0p997959/target_fit.csv), [parameters](ge_retry4/purchase_90_stayer_80/selected_repeat_7_0p997959/parameters.csv) |
| .90 / .90 | [closure](ge_retry5/both_90_90/root_5_0p997698/closure.json) | [fit](ge_retry5/both_90_90/root_5_0p997698/target_fit.csv) | [parameters](ge_retry5/both_90_90/root_5_0p997698/parameters.csv) | [summary](ge_retry5/both_90_90/root_5_0p997698/standard_diagnostics/summary.json) | [fit](ge_retry5/both_90_90/selected_repeat_6_0p997698/target_fit.csv), [parameters](ge_retry5/both_90_90/selected_repeat_6_0p997698/parameters.csv) |
| 1.00 / .80 | [closure](ge_retry6/purchase_100_stayer_80/root_5_0p995841/closure.json) | [fit](ge_retry6/purchase_100_stayer_80/root_5_0p995841/target_fit.csv) | [parameters](ge_retry6/purchase_100_stayer_80/root_5_0p995841/parameters.csv) | [summary](ge_retry6/purchase_100_stayer_80/root_5_0p995841/standard_diagnostics/summary.json) | [fit](ge_retry6/purchase_100_stayer_80/selected_repeat_6_0p995841/target_fit.csv), [parameters](ge_retry6/purchase_100_stayer_80/selected_repeat_6_0p995841/parameters.csv) |
| 1.00 / 1.00 | [closure](ge_retry7/both_100_100/root_4_0p990404/closure.json) | [fit](ge_retry7/both_100_100/root_4_0p990404/target_fit.csv) | [parameters](ge_retry7/both_100_100/root_4_0p990404/parameters.csv) | [summary](ge_retry7/both_100_100/root_4_0p990404/standard_diagnostics/summary.json) | [fit](ge_retry8/both_100_100/selected_repeat_1_0p990404/target_fit.csv), [parameters](ge_retry8/both_100_100/selected_repeat_1_0p990404/parameters.csv) |

Same-q0 fixed-price results are separate. At q0, .90/.90 has births/HH .115266648, ownership .675949055, rooms/HH 6.005350742; 1.00/.80 has .115221425, .668716002, 5.983263355; 1.00/1.00 has .115181114, .692838064, 6.000115817. Full packets: [both .90](retry7/results/both_90_90/closure.json), [buyer 1.00 / stayer .80](retry8/results/purchase_100_stayer_80/closure.json), [both 1.00](retry9/results/both_100_100/closure.json).

Adapter and provenance: [GE driver](run_purchase_ltv_ge.py); native mortgage gates from [fixed_price_responses.py](../../fixed_price_responses.py); closure/root/scaled-step math from [credit_ge_v1/run_ge.py](../../../../../credit_ge_v1/run_ge.py). Per-case certification receipts: [purchase .90/stayer .80](ge_retry4/purchase_90_stayer_80/ge_completion.json), [both .90](ge_retry5/both_90_90/ge_completion.json), [purchase 1.00/stayer .80](ge_retry6/purchase_100_stayer_80/ge_completion.json), [both 1.00](ge_retry7/both_100_100/ge_completion.json). The original .95 inherited-PRE failure is [here](ge_retry3/failure.json). One earlier multi-case process ended at fresh-interpreter re-authentication after its first selected point; that selected case is preserved at [ge_retry4 case receipt](ge_retry4/purchase_90_stayer_80/ge_completion.json). The both-100 repeat tables passed; a post-check writer counter error is retained at [ge_retry8 failure](ge_retry8/failure.json), with the passed repeat cell in [its receipt](ge_retry8/both_100_100/completed.json).
