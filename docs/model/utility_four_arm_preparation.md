# Overnight utility comparison: experimental contract and verification

September 26, 2026. An eight-hour, four-arm recalibration on Torch. Tommaso
delegated the experimental choice for tonight: “i am not sure, do what you
think its best, i'll assess tmrw.” The reviewed reference-rent normalization
below is the chosen assumption for this comparison. This authorization permits
the experiment, not adoption into the paper's model. Objective execution also
requires the separately frozen source, target, budget and launch contract.

## What yesterday already established

The immutable reference is `nightpair_20260925_v1`, under
`/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/`.
Its reviewed launch-lock SHA256 is
`6443195fa3f7de0dce5cc8a4c05e2709b99d586421c96a35dba3c07e93e061a1`.
The reference source and results remain untouched.

Reuse its eight entry/transition helper checks; four native seed-by-tax
preflights; complete 728-objective execution with 4,241 completed stationary
solves; numerical, fiscal, entry, purchase and value checks; original selected
point against both exact repetitions; complete target/parameter reporting;
and unchanged 17-figure diagnostic packets. The existing reviews are evidence
of those checks, not proof of convergence or of any new utility arm.

Yesterday already implemented the half-at-16/half-at-20 birth-entry queue,
division by 2.1 exactly once, unchanged dependent departure probability 2/9,
the 0.99 annual discount-factor cap, first-birth rooms target 1.465,
wealth/earnings target 6.927, bequest-flow target 0.007, externally fixed
bequest shift, adopted housing-cost inputs, and the national 2005–2006 ACS
vintage for the 2007 reference economy. Full precision and definitions remain
in the pinned executable objective and its complete provenance contract.
None of these needs another general reconciliation.

## The actual changes

The common fiscal change is author-adopted: the pension divided by mean gross
working-household earnings is fixed at the CPS 2007 ratio, 0.229. If workers
have mass \(W\), retirees mass \(R\), and average gross earnings \(\bar y\),
equal retiree pensions \(b=\rho\bar y\) imply
\[
\tau W\bar y=Rb,
\qquad \tau=\rho R/W.
\]
The native demographic recursion derives a baseline tax of 8.028%; it is not
the frozen 8.751% tax and is not fitted to utility outcomes. Each completed
objective must independently certify the pension ratio on its solved
distribution. Future transitions hold this baseline tax fixed and balance
equal pensions against the actual tax base. This comparison runs stationary
recalibrations, not policy transitions.

The four experimental arms combine the same household equivalence scale with
either a parenthood housing requirement or child-dependent housing expenditure
shares; each is paired with either \(v(m)=m\) or \(v(m)=m^{0.86}\). Here \(m\)
means children currently at home, not children ever born. The power 0.86 is a
fixed sensitivity restriction, not an external estimate of this model's
curvature. The fertility benefit coefficient \(\psi\) is separately normalized
to completed fertility 2.1 in every objective, using the same tolerance and
normalizer. The first-ever-birth utility cost remains a separate parameter.

The floor arms estimate eight structural coordinates. The share arms replace
the single housing requirement by two share coefficients, so estimate nine.
Seven coordinates are common: housing supply level, discounting, owner-service
premium, first-birth cost, two fertility taste-shock scales and bequest strength.
All twelve positive-weight target rows remain, plus the unscored fertility
normalization. Having twelve rows does not establish local rank or global
identification; curvature is not added as another free parameter.

## Experimental common units for the share alternative

Let \(\alpha_m\) be the nonhousing expenditure share and \(r^*>0\) a fixed
reference rental price in the model's period units. The experiment uses the
same reference rent in every arm: 0.110, taken from the frozen low-tax selected
2007 economy. This is an explicit experimental choice, not a measured welfare
equivalence scale. Keep the historical unnormalized share implementation as
a separate reference.

For an unconstrained renter with discretionary expenditure \(E\), the
optimal unnormalized Cobb–Douglas composite is \(K(\alpha,r^*)E\), where
\[
K(\alpha,r^*)=\alpha^\alpha
\left(\frac{1-\alpha}{r^*}\right)^{1-\alpha}.
\]
Normalize the composite by
\[
A(\alpha_m)=\frac{K(\alpha_0,r^*)}{K(\alpha_m,r^*)},
\qquad Q_m=A(\alpha_m)c^{\alpha_m}s^{1-\alpha_m},
\qquad u_m=U\!\left(\frac{Q_m}{e(m)}\right)+\psi v(m).
\]
At the reference prices, a household needs \(e(m)E_0\) to attain the same
material living standard as a childless household spending \(E_0\). The
floor specification additionally needs \(r^*h_P\) while children are home.
This comparison concerns material utility and excludes the direct benefit of
children. It does not equate total utility across household sizes, eliminate
tenure/purchase constraints, or imply that owner households attain the rental
reference allocation.

Since \(\alpha_m=\alpha_0\) in childless states, \(A(\alpha_0)=1\) exactly.
The floor arms also retain their existing material composite unchanged. The
native CRRA multiplier becomes \(e(m)^{\sigma-1}A^{1-\sigma}\); with
\(\sigma=2\), it is \(e(m)/A\). Conditional spending shares are unchanged by
the multiplicative normalization, but utility comparisons across child states
change. Thus this is an economic restriction to assess in the experimental
results. A bounded independent mathematical review verified the formula,
childless invariance and native array routing; the author delegated its use
for tonight without adopting it as the preferred specification.

## Finite search and critical path

The read-only timing reduction covers every one of yesterday's 728 objective
receipts. The sum of native solve times per normalized objective was 714.569
seconds at the median and 1,025.743 seconds at the 90th percentile, with a
maximum of 1,914.117 seconds. These fields exclude setup, observers, scoring,
checkpoint writing and reporting. The proposal adds an explicitly assumed
120-second allowance, yielding 1,145.743 seconds per objective for planning.
Every smoke and search case records full wall time. The historical estimate
is a planning input; the shared absolute deadline governs even if actual case
times are longer.

Each of four arm jobs has ten single-thread workers. Sampled completed
nightpair workers used about 6.3 GiB; the proposed allocation is 120 GiB per
arm, with 40 simultaneous numerical workers across Torch. This is a memory
allowance, not a new peak-memory measurement for the share arms. Slurm enforces
the allocation. All four controllers must become ready within thirty minutes
before the shared clock begins. Local computation is not required, and the
cluster owns progress, termination and export.

Within each arm, two identical smoke objectives run sequentially through the
same case loop. The initial population contains 40 points: the same two
previously fitted common structural seeds, 14 full-range joint draws, twelve
medium joint draws and twelve local joint draws. Every structural coordinate
can move. Medium and local Gaussian standard deviations are 0.15 and 0.04 in
the declared transformed unit intervals, reflected at the bounds. Full-range
draws cover those complete intervals. The share starting coefficients are
explicit search seeds, not estimates.

Three subsequent differential-evolution generations contain 40 joint trials
each, with mutation factor 0.7 and crossover probability one. Common parameter
draws and donor-index streams match across arms; arm-specific survivor
selection can diverge. The full inherited bounds remain, including the 0.99
discount-factor cap. Selection freezes before two final repetitions. Neither
failed cases nor duplicates create replacement search slots.

The first verified smoke supplies the identical initial-population seed, so
there is no third unchanged evaluation of that point. This is at most 163 new
objectives per arm, or 652 overall, for the same 160-member-per-arm search
history across stages. At the inherited
normalizer's theoretical limit of 23 stationary solves per objective, the
absolute solve-count bound is 14,996. **That is not expected workload.** The
observed median was six solves and maximum fourteen. The objective wall cap
is 3,100 seconds, and the shared eight-hour deadline takes precedence over
case counts.

The conservative planning path is two sequential smoke waves, followed by
sixteen search waves (four per initial/generation barrier), then a reserved
75 minutes for both repetitions scheduled concurrently within each arm and
15 minutes for export. This totals 7.229 hours, leaving 46 minutes of slack.
Four DE generations would exceed the conservative eight-hour plan. Actual
work that reaches a deadline is reported incomplete; no automatic extension,
failed-run retry or gate relaxation is permitted. Heartbeats and latest/best
receipts must continue at least once per minute/case.

## Common unresolved assumptions

The native first-birth room response remains an unmatched proxy for the PSID
estimator. The model bequest observer sums positive estates of all households,
whereas the SCF target is child-directed. The inherited entry-wealth marginal
and income-rank coupling are held fixed; the common-scale candidate is not
adopted. The entry stock and any later simulated inheritance must not be
funded twice from the same estate pool. The estate recipient and mortality
diagnostics remain separate and cannot silently alter these arms.

## Verification and execution

Torch job `18566460` completed in 50 seconds. Nineteen focused tests and all
four native preflights passed, with both structural seeds bound in each arm.
The preflights used zero equilibrium solves. The floor/linear branch preserves
native shared arrays exactly at the same parameters; the other branches
restrict changes to the declared benefit/material arrays. All four fiscal
bindings reproduce the adopted ratio.

The integrated source snapshot `utility_four_arm_preparation_20260925_v2`
passed 54 tests and all four two-seed native preflights in Torch job
`18567413` (52 seconds, zero equilibrium solves). Tests cover subprocess
deadlines and owned process groups, complete finite barriers, explicit
inadmissible proposals, refusal to repeat failed smokes, source/contract
binding, both original-versus-repeat comparisons, and complete report rows.

Source entry points are `code/model/tools/run_e5f_utility_comparison.py`,
`e5f_utility_comparison_runtime.py`, and `e5f_utility_comparison_design.py`.
`code/cluster/check_e5f_utility_comparison.sh` runs the bounded preparation
checks on Torch. The compact packet is indexed at
`output/model/utility_comparison/README.md`.

The full controller is `run_e5f_utility_comparison_search.py`; the matching
collector is `collect_e5f_utility_comparison.py`. The cluster launcher reserves
four jobs with ten CPUs and 120 GiB each. Both exact smoke objectives execute
sequentially through the real subprocess loop; all four arms must pass the
scientific comparison before any starts its search. A failed smoke cannot
trigger selected-point repetitions. Named infeasible parameter proposals
consume their search slots; every other error or timeout stops new dispatch
without a retry. More than half inadmissible proposals at a completed barrier
stops further generations. Every unrun or incomplete slot remains visible.

The final readout compares each original selected point separately against
both repetitions and retains all target/parameter rows and all 17 standard
figures. The collector renders every PDF page for subsequent visual review;
an exported packet is not automatically marked visually reviewed. Report
layout is checked separately using an explicitly labeled historical fixture,
which is not a new model result. The compact output index records the launch
contract, job IDs and reporting validation. No author draft, slides, mock,
frozen result or estate job is changed by this experiment.
