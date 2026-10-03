# Calibration and model status

**Reconciled October 3, 2026.** This is the consolidated current-state note. The
complete previous status, previous memory, reconstruction evidence and preservation
manifest are in [the context archive](calibration_archive/context_refresh_20261002/README.md).
Read historical chronology there, in dated daily notes or in the named experiment
packets when needed.

The current working reference remains the verified **soft purchase-constraint
candidate, loss 23.078309**. It is experimental, not a certified paper baseline.
The matched original versus alternative timing search is now terminal. Verified
as of **October 3, 2026, 06:00 New York**, arrays **19086987 and 19087556**
have 48 terminal chains: 46 passed fresh native selected-point checks and two
(chain 23 in each arm) found no admissible candidate. The lowest verified
original-timing loss is **18.44530519407432** (chain 15); the experimental
alternative-timing loss is **13.771131463467462** (chain 13). These searches
do not certify optimizer convergence or adopt either result as a new reference.
An isolated post-interest-timing recalibration with a narrower experimental
PSID wealth numerator has also finished: ten of ten chains passed native
verification; its lowest new-contract loss is **48.170707377609034** (chain 2).
The target contracts differ, so these native losses cannot be ranked directly.

## Reference identities and navigation

| Object | Identity and role | Authoritative evidence |
|---|---|---|
| Frozen paper reference | `paper-baseline-2026-09-14`; retained checkout `tmp/paper_baseline_sep14/`. Preserve source and results. | `PAPER_BASELINE.md` in that checkout; [baseline checker](code/model/tools/check_paper_baseline.py) |
| September 28 fixed-economics reference | Older equilibrium/normalization and utility objects. A refactor oracle, not interchangeable with the current soft calibration. | [refactor report](output/model/publication_refactor_20260929/REPORT.md), [refactor runtime](code/model/refactor_lab/README.md) |
| Current soft selected point | Original timing; chain 16 / case 0046; verified loss 23.078309. | [selection](output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_selected.json), [verification](output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_verification.json) |
| Matched timing comparison | Soft constraint, original versus experimental post-interest transaction timing; 24 matched starts per arm, 48 terminal chains, 46 verified. Lowest verified losses 18.445305 and 13.771131. Neither fit is adopted. | [48-chain receipt](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/collection.json), [complete fit and parameter readout](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/RESULTS.md), [driver plan](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/driver_plan.json) |
| Experimental wealth-numerator comparison | Post-interest timing with one new PSID aggregate wealth/earnings target; ten verified chains, lowest new-contract loss 48.170707. No adoption. | [ten-chain verification](output/model/fixed_reference_economics_20260928/alternative_wealth_local_20261003_v1/collection/verification.json), [three-arm comparison and full tables](output/model/fixed_reference_economics_20260928/alternative_wealth_local_20261003_v1/comparison/COMPARISON.md) |
| Historical purchase-rule comparison | Fresh hard and quarter results: 88.588403 and 48.319938. Older policy exercises use earlier points with losses 97.011220 and 51.556036. | [fresh calibration readout](output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/fresh_calibration_v1/README.md) |
| Historical transition initializer | Normalized-v1 chain 20 / case 0028_nm, loss 30.371888. Its one-shock transition attempt failed acceptance. | [transition deployment/readout](output/model/transition_readiness_v1/normalized_restart_v1/resume_preparation/deployment/v2/README.md) |

Active model code is under [`code/model/`](code/model/README.md). The timing
experiment and calibration wrapper are under
`code/model/experiments/purchase_timing_sandbox/`
and `code/cluster/soft_timing_calibration/`. The interactive Python interface is
documented in [MODEL_PLAYGROUND.md](code/model/tools/MODEL_PLAYGROUND.md).
Check executed saved inputs and the deployment manifest before assuming a generic
model command reproduces the current normalized calibration.

## Working economic and accounting contract

These are the objects executed at the selected soft point. Retaining them in this
comparison does not resolve every empirical or publication issue.

**Units and household choices.** One period is four years. Households choose
consumption, housing, tenure, saving and fertility over the lifecycle. Children
ever born, `n`, and children currently at home, `m`, are different states.
The highest birth-count state represents three or more children. Housing is
measured in rooms; owner products have 2, 4, 6, 8 or 10 rooms, and rental
housing has the retained six-room cap. The active solution has 120 asset nodes and
nine income states. Legacy labels containing “B15” do not establish 15 executed
income states.

**Utility.** Retained risk aversion is \(\sigma=2\), the consumption share is
\(\alpha=0.733\), and there is no adopted child-dependent consumption-share
change. The equivalence scale is
\[
e(m)=\left(\frac{2+0.7m}{2}\right)^{0.7}.
\]
The parenthood housing floor is a physical-room object applied before the housing
taste parameter \(\chi\). The child benefit is \(\psi m^{1-\gamma}\), with zero
benefit at \(m=0\); the first-birth cost and bequest motive are separate objects.
Historical lower-\(A(m)\), renter-borrowing and other utility probes are
diagnostics, not silently adopted components of this reference.

**Earnings and entry.** Earnings use a deterministic age profile and a single
persistent four-year process: \(\rho=0.7345934906\), innovation standard deviation
\(0.4838308245\). There is no separately adopted permanent-type plus transitory
decomposition. See the [external earnings estimate](output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/single_process_external_estimate.json).

The selected input uses the provisional `nonnegative_mean` entry-wealth mapping:
negative five-bin mean wealth ratios are floored at zero after binning, then
positive means are rescaled to preserve the retained mean. Effective ratios are
\([0,0,0.037866,0.127705,1.127294]\). Mean entry wealth is \(0.186520\);
mean annual entrant income is \(0.720554\); the zero-wealth share is \(0.666948\).
This is not a recoding of every raw survey observation. Tommaso provisionally
preferred this option; production paper adoption and donor/recipient accounting
remain separate decisions. Earlier negative-entry-cell infeasibility is not a
current blocker at this point. The executed entry contract is in
[the selected input contract](output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_postcheck/input_contract.json).

**Finance and timing.** Annual real interest is 0.02, giving
\(R=1.02^4=1.08243216\). The financed share is \(\phi=0.8\), so a down payment
uses \(1-\phi\). Renter debt capacity is zero. Selling cost is 0.06.
Annual depreciation \(0.0141614372\) is compounded to a four-year rate
\(0.0554537908\); annual property tax \(0.0105983608\) is multiplied by four to
\(0.0423934431\). These conventions are explicit, not interchangeable
annualizations.

Let \(b\) be beginning net financial wealth, \(S\) net sale proceeds, \(Q\) the
purchase cost, \(y\) current income, \(c\) consumption and \(K\) other costs.
The original convention is
\[
b'=R(b+S-Q)+y-c-K.
\]
The alternative is
\[
b'=Rb+S-Q+y-c-K.
\]
The change moves current net transaction financing outside the interest factor.
At held choices the difference is \((R-1)(Q-S)\). Sale, forward-budget and
solvency maps must use the chosen convention consistently. Existing debt is
already in \(b\); do not subtract it again. The applicable owner ending-debt
floor remains \(b'\geq-\phi Q\). A Pro recommendation motivated the comparison;
it did not authorize replacing the reference timing. Hard versus soft purchase
constraints also remain a distinct specification choice.

**Demography and closure.** The selected soft equilibrium fixes household scale
\(N=1\), solves price using stationary birth renewal, and derives the housing
supply scale \(H_0=6.757074\) from market clearing. \(\psi\) is jointly free
in the ten-parameter calibration; it is not independently normalized back to
2.1 at each proposal. The unscored normalization row checks adjusted births
divided by entry.

The top-count correction adds \((w-3)\) times the flow into the three-or-more
state to recorded births, where \(w\) is its representative count. Potential
entry is adjusted births divided by 2.1, once. A half-age-16 / half-age-20
entry queue is implemented, approximating mean entry age 18. Child departure
uses the retained \(2/9\) rule, without a newborn exemption; parental death
removes dependency without an additional entry flow. See
[adult-entry accounting](code/model/intergen_eqscale_seq_optimized/adult_entry.py).
This is a working household demography approximation, not explicit tracking of
individual offspring genealogies or certified resident-person accounting.
The October 2 check found no extra factor of two in this path; different
three-plus weighting and age windows still require care when comparing fertility
statistics.

At the selected point, entry is \(0.0617334562\), adjusted births
\(0.1296402579\), the birth-renewal residual is approximately
\(-7.08\times10^{-10}\), the housing residual is zero, and the PAYGO residual
is \(2.71\times10^{-14}\). Payroll tax is \(0.0802807096\) and the period
pension is \(0.9177840475\). The accepted CPS ASEC pension/gross-earnings ratio
\(0.2294460119\) underlies the retained PAYGO normalization; a transition
holds the adopted payroll tax and balances period pensions, rather than adding
unrecorded government spending. See
[verified closure](output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_postcheck/selected_postcheck/phase_b_ge/selected_repeat_final/closure.json).

Housing supply elasticity 0.63 is externally fixed provisionally.
The current equilibrium is a closed single housing market, with no adopted
outside-entry/geographic migration valve. For a policy endpoint under the
retained supply function, the reported population scale uses
\(N=H_{0,\mathrm{base}}/H_{0,\mathrm{derived}}\). Holding \(H_0\) fixes the supply
schedule, not the quantity at every price. This population-scale effect must
not be described as a total fertility effect.

## Verified current fit and parameters

The local selected-point repeat passed, reproducing loss **23.07830929416065**.
The largest saved-moment discrepancy was \(5.33\times10^{-15}\). The packet has
14 target rows, 31 parameter records and the established 17 plots. A first
cross-platform exact-double comparison was rejected for roundoff; the accepted
local repeat did not relax economic acceptance gates. This is selected-point
verification, not proof of optimizer convergence.

The gap below is model minus target. Blank weights identify the separate
normalization; zero weights identify validation rows. Rounding is for display;
the linked CSVs retain full precision.

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Stationary birth-renewal normalization | 2.100 | 2.100 | -1.486e-09 | — | — |
| Childless women, ages 40–44 | 0.198 | 0.202 | 0.004 | 35532.304 | 0.553 |
| One child among mothers, ages 40–44 | 0.214 | 0.213 | -0.001 | 26952.821 | 0.028 |
| Mean age at first birth | 25.976 | 26.032 | 0.056 | 139.828 | 0.436 |
| First births at ages 30+ (validation) | 0.249 | 0.234 | -0.015 | 0 | 0 |
| Net wealth / labor earnings | 6.927 | 6.551 | -0.376 | 7.595 | 1.073 |
| Child-directed bequest flow / wealth | 0.007 | 0.007 | -4.449e-04 | 5165289.256 | 1.023 |
| Older-household wealth dispersion (validation) | 3.516 | 2.968 | -0.548 | 0 | 0 |
| Mean occupied rooms | 5.729 | 5.957 | 0.228 | 128.021 | 6.639 |
| Ownership, heads ages 30–55 | 0.676 | 0.665 | -0.011 | 2339.362 | 0.282 |
| First-birth room response | 1.465 | 1.324 | -0.141 | 137.565 | 2.726 |
| Rooms: 3+ versus 1–2 children (validation) | 0.385 | 0.304 | -0.081 | 0 | 0 |
| Recent-parent ownership gap | 0.128 | 0.118 | -0.009 | 27055.823 | 2.371 |
| Children ever born, capped at 3, age 25 | 0.810 | 0.528 | -0.282 | 100.000 | 7.947 |

Source: [complete target-fit CSV](output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_target_fit.csv).
Early fertility and mean rooms account for approximately 63% of loss. Exactly
one child is conditional on being a mother, not an unconditional share of women.
The ten positive-weight moments match the ten free parameters in count;
informative rank and weak identification at this point have not been certified.

| Free parameter | Estimate | Lower | Upper | Near bound |
|---|---:|---:|---:|---|
| `beta_annual` | 0.967 | 0.940 | 0.990 | no |
| `chi` | 1.097 | 0.100 | 5.000 | no |
| `first_birth_fixed_cost` | 0.352 | 0 | 8.000 | no |
| `kappa_fert` | 0.124 | 0.020 | 50.000 | yes |
| `kappa_fert_continuation` | 0.363 | 0.020 | 50.000 | yes |
| `theta0` | 0.106 | 0 | 8.000 | no |
| `child_benefit_curvature` | 0.101 | 0 | 0.800 | no |
| `tenure_choice_kappa` | 0.014 | 0.001 | 0.100 | no |
| `psi_child` | 0.178 | 0.010 | 0.500 | no |
| `h_P` | 2.504 | 0.100 | 2.600 | no |

Source: [all 31 parameter records](output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_parameters.csv).
“Near” follows the saved screen: within 1% of the search interval, not actual
endpoint contact. The derived \(H_0\) has admissibility range \([0.2,80]\);
\(\theta_1=0.0081930841\) remains fixed. Fixed objects and the executed input
contract are linked above.

## Empirical target provenance and measurement limits

The row-by-row [provenance review](calibration_archive/context_refresh_20261002/target_provenance_review.json)
preserves builders, source records, estimators, samples, controls/clustering,
uncertainty where available, current weights, observer definitions and warnings.
It is a reconstruction of the executed contract, not a new target registry or
a certification of statistical design. Governing reviewed sources are
[the September 27 measurement review](docs/model/e5f_target_measurement_review_20260927.md)
and [accepted-input reconciliation](docs/model/accepted_input_reconciliation_20260926.md).

- **CPS fertility stocks:** pooled June 2004/2006 women ages 40–44, valid
  `FREVER` 0–20 and positive supplement weights. Childlessness is unconditional;
  exactly one child is conditional on mothers. Official annual generalized
  variance approximations exist, but pooled covariance is not design-certified.
  The model uses uniform within-period birth timing and a reproductive household
  member proxy, not a demonstrated female-exposure reconstruction.
- **Early fertility:** the same supplements at exact age 25, with children ever
  born capped at three. The recorded person bootstrap standard error is 0.028035,
  stratified by year; it is not CPS design-consistent. Weight 100 is a working
  calibration choice, not inverse survey variance. The model's uniform
  age interpolation is an approximation.
- **NCHS timing:** first births in 2003–2006, ages 12–49, mapped to retained
  four-year age-cell representative values. Annual dispersion, including
  0.084567 for mean age, is not a sampling standard error. The age-30-plus
  share is validation only.
- **PSID wealth:** 2005/2007 weighted net-wealth to gross labor-earnings ratio;
  wealth ages 18–85 and earnings ages 18–65 under the retained sample definitions.
  Reference-person bootstrap standard error is 0.417310. Wealth includes home
  equity. Older-household wealth dispersion is validation only; the retained
  PSID income filter/age support and model pension denominator do not coincide.
- **SCF bequests:** the 2007 mortality-weighted child-directed annual estate flow
  divided by net wealth is 0.007291. The empirical calculation is complete;
  recipient eligibility, timing and model entry mapping remain open. No survey
  standard error was established; a synthetic working uncertainty is not
  empirical precision.
- **AHS rooms:** 2007 occupied households, head ages 18–85, positive weights and
  valid rooms; 37,793 observations, literal public-use room topcode 21. Fay BRR
  standard error is 0.008934. The retained loss weight 128.020702 is not the
  inverse of this variance, and model room exposure remains an approximation.
- **ACS housing/ownership:** accepted national 2005/2006 targets. Ownership and
  recent-parent rows retain a DUE structure restriction absent from the model;
  family rooms uses resident own children, a minor-child screen and rooms capped
  at nine. The old 42-metro bootstrap weights were retained; they do not
  establish national sampling variance. The recent-parent model observer uses
  births into previously empty-dependent homes versus currently empty homes,
  including former parents. It is not the exact ACS oldest-child-age estimator.
- **PSID first-birth rooms:** the authoritative A2h Sun–Abraham contrast is
  calendar years \(+3/+4\) versus \(-3/-2\), baseline reference persons/spouses,
  biological first births and confirmed-childless controls. It uses person/year
  fixed effects, age/education controls and individual clustering: 117,853 fitted
  observations, 9,310 clusters, standard error 0.050270. Tommaso selected the
  rounded 1.465 target for this working calibration. The model destination-period
  contrast does not replicate the empirical panel estimator or cohort selection.
  Historical 0.600/0.720 targets are not current. Contemporaneous-income-controlled
  OLS is a robustness exercise: 1.300 versus approximately 1.461 without income on
  the matched sample. It has not replaced the target.

These limits remain visible. Removing, reweighting or replacing a target needs an
identification argument and an explicit contract change; this refresh makes none.

## Active timing comparison: execution state and acceptance

The deployment receipt identifies remote root
`/scratch/td2248/projects/soft_timing_calibration_20261002_v2`.
Attempt-one smoke 19085105 failed repeat verification in a reused interpreter.
The repaired driver uses a fresh child interpreter within the same absolute
deadline. Attempt-two smoke 19086529 passed both arms, including the full
14-row fit, 31-parameter record, 17-plot packet and selected-point repeat. No
economic change or gate relaxation was introduced by this process repair.
The production array 19086987 was subsequently submitted.

The alternative timing at the common starting coordinate has loss
**57.186333** and price **0.730719**, versus **23.078309** and **0.726639** for
the original reference. This is a held-coordinate comparison after equilibrium
solution, not a ranking of separately recalibrated specifications. Both complete
native smoke tables and parameter records are linked by
[the smoke collection readout](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/smoke_collection/readout.json).

The author subsequently expanded the design to 48 chains, 24 matched starts per
arm. Array 19087556 adds 40 chains (indices 4–23 per arm) in remote v3. The
[expanded plan](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/expanded_start_plan.json)
retains bounds, targets, gates and the passed v2 economic implementation. Starts
combine the original four, seven other historical soft candidates, eight nearby
variations and five broader starts. The 240-chain proposal was superseded.

**Final matched search, verified as of October 3, 2026, 06:00 New York.** Both
arrays are terminal, with zero active Slurm tasks at collection. Of 48 chains,
46 passed the fresh native selected-point and exact-repeat gates. Original and
alternative chain 23 ended with no admissible selected point; neither enters
the winner comparison. The lowest verified original-timing loss is
**18.44530519407432** (v3 chain 15), and the lowest verified experimental
post-interest loss is **13.771131463467462** (v3 chain 13). Both winner packets
have 14 target-fit rows, 31 parameter rows, 17 standard plots in each native
root and exact repeat, and exact repeated tables and plot hashes. See the
[complete 48-chain collection receipt](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/collection.json),
[original fit](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/original_target_fit.csv),
[alternative fit](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/alternative_target_fit.csv),
[original parameters](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/original_parameters.csv), and
[alternative parameters](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/alternative_parameters.csv).
The selected original estimate of \(h_P=2.6\) contacts its upper bound; the
alternative estimate is \(2.59376\), near the same bound. Early fertility is
0.533079 and 0.533805, respectively, against a target of 0.809528. The fit
improvement therefore does not resolve this target, and neither optimizer has
a convergence certificate. The common target fingerprint is
`db60605ef444b747b2ed7f9482b6b8c90ef5e8b1d75c49a3ac1ff6f4275c9ba1`,
the weight fingerprint is
`2391cd2d4a39a6669a405be34ff14116ad314354353173cb619f7fa7c66043b0`,
and the selected soft checkpoint SHA-256 is
`b5e21a8584fa536a3740039b3e54480f3c318a63b96266bd31406712c5da991c`.
These are experimental recalibrations; the existing working reference remains
unchanged.

Each chain had a six-hour budget, at most 250 objective evaluations,
a retained 1,800-second finalization reserve, one CPU and 24 GiB; case lifecycle
limits and checkpoints are defined in the driver plan. The combined maximum
budget was 288 core-hours, with at most 48 concurrent single-core chains.
No failed chain was silently retried, and the selected points remain experimental.

The input contract pins target SHA-256
`db60605ef444b747b2ed7f9482b6b8c90ef5e8b1d75c49a3ac1ff6f4275c9ba1`
and weight SHA-256
`2391cd2d4a39a6669a405be34ff14116ad314354353173cb619f7fa7c66043b0`.
The v2 upload archive SHA-256 is
`f7a8fec4ff370fd3690c0d0068ca595b75a17dd8aaac3bd47f6009ef73ecd68b`;
the manifest verifies the pinned source set and enumerates the driver repair.
The v3 expansion archive SHA-256 is
`604096fe0bc9d8d7f7f3de7b52fbf3a7cad411dd4478f4405e4add68f54371a0`;
its [stage and submission evidence](output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/deployment/expanded/status.json)
pins the parent v2 production/archive and new start table. The parent v2
receipt above identifies the passed two-arm native smoke. The terminal deployment state and
submission receipt supersede older preparation prose/flags.
Use the owning “Compare interest timing and credit” chat and its existing
monitor; do not create a duplicate monitor from this documentation refresh.
Older “prepared/not submitted” prose in design notes is superseded by the
structured deployment receipt.

## Experimental narrower wealth target: completed local comparison

**Verified October 3, 2026, 08:07 New York.** The isolated local post-interest
transaction-timing search ended with ten terminal chains, all ten passing fresh
native selected-point and exact-repeat checks. Chain 2 has the lowest verified
loss **48.170707377609034** under its *experimental* target contract. The
new pooled 2005/2007 PSID aggregate net-wealth/annual gross-labor-earnings
target is **4.45838713455674**, excluding business/farm equity, other real
estate and vehicles while retaining catch-all other assets; the model moment
at chain 2 is **6.561732864006831**. The old numerical weight
**7.595098472533724** was retained for a controlled sensitivity, and no new
standard error has been estimated. Entry wealth and income distributions, the
bequest target and all other targets, model observers and economic inputs remain
unchanged relative to the alternative-timing arm. The model bequest denominator
and the new PSID wealth numerator have not been reconciled.

The [complete new-target collection](output/model/fixed_reference_economics_20260928/alternative_wealth_local_20261003_v1/collection/RESULTS.md)
includes the authoritative 14-row fit and 31 parameter records; its native
`target_fit.csv` uses the **old** wealth target only as a diagnostic, while
`target_fit_new_contract.csv` is authoritative for this experimental arm. The
[three-arm comparison](output/model/fixed_reference_economics_20260928/alternative_wealth_local_20261003_v1/comparison/COMPARISON.md)
links every target, model value, gap, weight, contribution, parameter bound and
source hash. Arithmetic rescoring of the same saved moments yields old-contract
scores **18.445305**, **13.771131**, and **15.580542**, and new-contract scores
**51.765179**, **53.064444**, and **48.170707**, for original timing,
alternative timing and new-wealth chain 2 respectively. These scores are not
new solves. The age-25 children-ever-born stock is **0.535240** at chain 2
against **0.809528**, essentially unresolved across the three points. The ten
positive-weight moments equal the ten free coordinates in count, but informative
rank and optimizer convergence are uncertified. The unequal 48-chain and
10-chain budgets do not establish that the new wealth target is unreachable.
Neither the target nor the chain-2 parameter vector is adopted; the working
soft selected point remains unchanged.

## Interactive inspection and numerical readiness

Tommaso wants direct Python parameter edits, quick fixed-price solves, policy
inspection and whole-population aggregates from the existing engine.
[MODEL_PLAYGROUND.md](code/model/tools/MODEL_PLAYGROUND.md) documents
the preferred three-script entry: [run_model.py](code/model/run_model.py),
[plot_model_policies.py](code/model/plot_model_policies.py), and
[plot_model_aggregates.py](code/model/plot_model_aggregates.py). Internal and
common external parameters are explicit editable literals. Each run writes its
full native arrays and effective inputs to a separate `tmp/model_runs/` folder;
both plotters load those arrays and use editable Matplotlib calls. The older
interactive interface remains optional. A fixed-price household solve is not
the normalized calibration GE; it reports numerical, renewal and housing
diagnostics without claiming a new equilibrium.

The [October 3 workflow verification](output/model/fixed_reference_economics_20260928/model_control_scripts_v1/verification.json)
at the authenticated selected soft point matched all 11 checked baseline arrays
exactly. The household solve took 6.2 seconds; the full run with 17 standard
figures, saving and round-trip validation took 9.9 seconds. Saved loading
reproduced 67 solution arrays and 49 parameter arrays exactly, without reference
initialization or solving. The two plotters produced eight policy and seven
aggregate figures. A failed validation preserves the last complete run.
The verification receipt records source hashes and its verified-as-of time.
The browser explorer requires its server; opening its HTML alone is insufficient.

The explorer distinguishes inherited tenure from the chosen branch, exposes
asset values and grid indices, and reports aggregates independently of slice
selectors. Conditional versus averaged policy comparisons use the same age,
income and family state. Distribution weighting must distinguish the beginning
distribution after fertility from the saved post-decision distribution.

The [asset-grid diagnosis](output/model/fixed_reference_economics_20260928/asset_grid_diagnosis_v1/README.md)
now contains six sequential fixed-price tests with unchanged economics,
prices and entry rules. On the 120-node \([-12,3000]\) grid, 99.8665% of
pre-decision mass lies in \([-6.4,33.6606]\), with zero endpoint mass.
Extending the high-wealth tail to 6000 with 124 nodes changed none of the tested
occupied averaged policies/aggregates. Refining the occupied region to 214 and
402 nodes changed mean assets by approximately +0.763% and +0.962%; all-age
ownership moved from 66.650% to 66.410% and 66.333%. This supports attention to
core resolution, not a claim of full policy, target or equilibrium convergence.
Large pointwise differences can coexist with small weighted average differences.

The older matched 160×15 versus 120×9 timing comparison found a 2.503× speedup
under its retained credit setting. It does not certify the current zero-renter-debt
calibration, dated transitions or an automatic coarse-to-fine production handoff.

## Historical policy and transition results to retain

**Hard/quarter purchase rules.** All 24 fresh-calibration selected postchecks
completed; best hard loss is 88.588403, quarter loss 48.319938. The hard housing
floor reaches its upper bound 2.6; the quarter floor is 2.572098. Full fit,
bounds and diagnostics are in the fresh-calibration readout linked above.
These are accepted selected solutions, not proof of optimizer convergence or
author adoption.

The earlier 97.011220/51.556036 points produced accepted temporary tax-shock
paths at 48/64 periods. First-period birth-flow changes were approximately
−0.542%/−0.539% for hard and −0.431%/−0.429% for quarter, respectively.
These are first-four-year birth-flow changes at those historical calibrations.
Permanent terminal steady states passed, implying population-scale changes
−2.324% hard and −2.105% quarter. All four permanent dated transitions failed
root or terminal gates, including the recovered hard-64 run. There is no
accepted permanent transition impact. See
[mechanism deployment](output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/mechanism_deployment/README.md).
The older quarter fixed-price birth increase, approximately +0.184%, and its
negative equilibrium response support a rent-offset mechanism at that point;
they are not a demonstrated mechanism for the current soft solution.

**Historical successive-surprise contract.** The authorized historical design
fits four successive surprises, each believed permanent until the next arrives,
or one permanent 2007 shock to the final 2020–2023 window. Carry the complete
household state and both entry queues between surprises; do not substitute an
announced four-shock path. The block0506 attempts 18801439/18801451 failed before
fitting shocks. A later pension-correction diagnostic is not an estimated history.
The [historical estimator packet](output/model/fixed_reference_transition_20260928/four_shock_v1/README.md)
and [bounded diagnostic](output/model/fixed_reference_transition_20260928/four_shock_v1/budget_diagnostic_v2/README.md)
preserve the request and launch history; their old draft-plan disabled flags do
not negate subsequent author launch authorization. No accepted four-surprise fit
is established by the reviewed receipts.

**Normalized-v1 one-shock transition.** Historical job 18995772 hit its time
cap after two trials. Its best shock benefit was \(\psi=0.129489443\), within
bounds \([0.00171989,0.34397799]\); the 2020–2023 model statistic 1.703658
missed target 1.645750 by 0.057908. Four-window target fit and all shock
parameters are linked in the transition readout above. Exact replay and
intermediate housing roots passed, but strict terminal-state/horizon acceptance
failed and final diagnostic plots were not produced. Monitors were paused;
there was no approved automatic budget extension. It is not a certified policy
initializer for the current soft point.

**Paper-facing artifacts.** The [Corina source map](latex/corina_progress_20260930/source_map.md)
and [slides](latex/corina_progress_20260930/corina_progress.tex) use the earlier
hard/quarter calibration and policy packets. They do not display the latest
88.588/48.320 results or the current soft 23.078 point. The coordinated draft,
slides and mock sources are listed in [latex/README.md](latex/README.md).
Experimental physical-floor/nonlinear-benefit and timing objects should not
be described as fully synchronized with the author manuscript. This refresh
does not change manuscript wording or adopt a new paper specification.

## Outstanding decisions and checks

1. **Timing and purchase constraint:** collect both recalibrated arms with native
   postchecks and complete tables before deciding whether to adopt the
   alternative timing. The hard-versus-soft economic choice remains distinct;
   prior categorical recommendations were not supported by a completed comparison.
2. **Identification and weights:** count is ten moments for ten free parameters.
   Check informative rank/substitution and the two near-bound fertility/continuation
   noise parameters. National ACS uncertainty and the early-fertility weight remain
   working choices requiring their own decision.
3. **Empirical observers:** resolve or explicitly maintain female versus household
   exposure, Sun–Abraham versus stationary first-birth observation, recent-parent
   residence/age definitions, and old-wealth income/age support. Do not call current
   proxies exact measurement equivalence.
4. **Entry wealth, estates and creditors:** current entry mapping is provisional.
   A retained gate reports positive entrant funding 0.011515 per period against
   provisional available net estates 0.120274, but this is not full donor-recipient,
   spouse, negative-estate/creditor or physical-resource settlement. The SCF
   calculation being complete does not close these mappings.
5. **Demographic and population interpretation:** explicitly state household,
   reproductive-member and resident-person units; retain the implemented queue and
   three-plus adjustment. Do not infer literal genealogical tracking, a factor-two
   correction, total fertility effects or a geographic migration closure.
6. **Numerical readiness:** complete relevant grid/target/GE convergence checks
   before production policy claims. Fixed-price grid tests and refactor parity
   do not certify a dated transition or coarse-to-fine handoff.
7. **Policy acceptance:** historical permanent paths and the one-shock normalized-v1
   path failed required gates. A current soft-policy result needs its own reconciled
   closure, reference bridge and complete visual packet.
8. **Publication contract:** distinguish an exploratory working calibration from
   author-adopted economics, preserve the frozen September 14 reference, and
   reconcile accepted changes with the authorized paper representations.

## How to maintain this note

Replace superseded state within its section. Record the source identity and
verification time for mutable claims; do not paste another chronological
transcript here. Put detailed chronology in daily notes or the experiment README,
durable preferences/gotchas in `memory/AGENT_MEMORY.md`, and old snapshots in
`calibration_archive/`. Advisory size budgets are a review signal, not permission
to truncate unresolved items, discard evidence or silently adopt experiments.
