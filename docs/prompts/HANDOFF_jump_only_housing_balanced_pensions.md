# Overnight handoff: parenthood housing requirement and balanced pensions

Prepared for Tommaso to give to the agent already responsible for the overnight calibration and related reruns. This packet records the preference decision reached in this conversation and reconciles it with the pension repair. It does not replace that agent's existing assignment, authorize a duplicate search, or claim that the combined model has been implemented.

Repository: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`.

## 1. Author decision and immediate objective

Combine the Social Security correction with a simpler household-needs specification in the forthcoming recalibration. Retain the current concave equivalence scale inside CRRA, but replace the first-child housing jump plus per-child housing slope with one housing requirement for parenthood. Additional dependent children continue to increase needs through the scale. They no longer each add an unavoidable increment to the housing minimum.

This is the selected preference direction, superseding this conversation's earlier recommendation to retain both housing parameters. Earlier discussion of a linear scale was exploratory and was not adopted. Do not launch a menu of preference alternatives instead of implementing this decision. Preserve the remainder of the overnight agent's authorized objectives and prioritize a verified combined benchmark before interpreting new policy results.

The preference change is a deliberate economic simplification. The pension change repairs a fiscal inconsistency. These are separate reasons for change and should be documented separately even though they will be incorporated into the same rerun.

## 2. Exact proposed utility

Let c be nonhousing consumption, s housing services, and m the number of currently dependent children. m is not completed fertility, total births ever, or a top-bin average. Maintain the model's existing child independence and aging transitions.

The externally fixed scale is

\[
e(m)=\left(\frac{2+0.7m}{2}\right)^{0.7}=(1+0.35m)^{0.7},\qquad e(0)=1.
\]

For a two-adult household, the childless normalization is one. Both 0.7 coefficients remain fixed; do not add a free equivalence-scale parameter.

Replace the old housing requirement

\[
\bar h_{\mathrm{old}}(m)=h_1\mathbf 1\{m>0\}+h_m m
\]

with

\[
\bar h_{\mathrm{new}}(m)=\bar h_P\mathbf 1\{m>0\},\qquad \bar h_P\ge0.
\]

The consumption–housing composite and flow utility are

\[
Q(c,s;m)=c^{\alpha_0}\left[s-\bar h_P\mathbf 1\{m>0\}\right]^{1-\alpha_0},
\]

\[
u_t(c,s;m)=\frac{[Q(c,s;m)/e(m)]^{1-\sigma}}{1-\sigma}+\psi_t m,
\qquad \sigma=2.
\]

Equivalently, at the maintained curvature,

\[
u_t(c,s;m)=-\frac{e(m)}{Q(c,s;m)}+\psi_t m.
\]

Keep alpha constant as in the active floor specification. Keep the existing child utility term, fertility-choice architecture, conception risk, other lifecycle preferences, and bequest block unless the overnight agent has separate explicit author instructions changing them. The formula above describes the material-consumption/children flow block; it does not erase those other model components.

Maintain zero nonhousing consumption floor. Do not add an outside multiplier e(m). Do not restore a child-dependent Cobb–Douglas expenditure share. The domain is c>0 and s>hbar_new(m). The housing requirement is in services; use the existing service-to-physical-housing mapping consistently for owners.

For a comparison that preserves the previous first-child requirement, initialize hbar_P = h1_old + hm_old. Simply zeroing hm while keeping h1 unchanged would also lower the first-child requirement. This initialization is not a new calibrated estimate: hbar_P must subsequently be jointly re-estimated. Review its bounds as a combined requirement rather than blindly inheriting the old jump parameter's interpretation.

## 3. Economic reason and functional-form implications

The interpretation is that starting a family requires suitable housing, while later children share housing and other household resources. Declining marginal child costs are intentional and economically plausible. The author explicitly rejected treating them as a defect merely because they complicate fitting the fertility distribution. Do not restore the housing slope just to counter scale concavity or obtain a preferred result.

The old specification applied both an increasing scale and an increasing housing minimum. Those channels are mathematically coherent and not proven to double count. However, the additional housing slope needs its own economic and empirical justification. A first-birth housing response alone does not justify an unavoidable housing increment for every later child.

The scale at m=0,1,2,3 is approximately 1, 1.2338, 1.4498, 1.6528. Its successive increments are 0.2338, 0.2160, and 0.2030. These are analytical inputs, not estimated expenditure effects. Concavity means additional resources needed to preserve a material living standard rise by progressively smaller amounts.

For an interior renter, let X=c+r*s be current expenditure and r the rental price of a housing-service unit. Define the discretionary-composite price

\[
P_Q(r)=\frac{r^{1-\alpha_0}}{\alpha_0^{\alpha_0}(1-\alpha_0)^{1-\alpha_0}}.
\]

Among parents, conditional utility is

\[
v(X,r,m)=-\frac{P_Q(r)e(m)}{X-r\bar h_P}+\psi_t m,\quad m\ge1,
\]

provided X>r*hbar_P and the renter optimum is interior. Since e''<0, this continuous extension is convex in m at fixed X. Equivalently, the material utility costs of additional children decline. The actual model chooses discrete fertility dynamically; this calculation is not an impossibility theorem about interior family sizes or model fit.

Conditional housing demand is

\[
s^*=\frac{(1-\alpha_0)X}{r}+\alpha_0\bar h_P,\quad m\ge1.
\]

Thus there is no further child-count shift in housing allocation among parents at fixed current expenditure and rental price. Larger families can still choose larger homes through endogenous spending, saving, tenure, sorting, and lifecycle behavior. Housing prices and rents still affect fertility because housing remains in the scaled bundle. Do not translate this conditional identity into “additional children have no housing cost.”

For a fixed effective material standard q=Q/e, required expenditure is r*hbar_P+P_Q*e(m)*q among parents. This remains a housing setup cost plus scaled discretionary expenditure. Removing the slope does not establish that an imported total-needs scale empirically measures the remaining bundle exactly. Joint consumption, housing, and saving evidence remains useful.

Calibration partly adjusts cost levels through child tastes and housing requirements. It cannot make the scale's shape unrestricted. A common linear child benefit changes the attractiveness of children but not this conditional curvature. Retain this distinction when assessing fit.

Literature interpretation: [Scholz, Seshadri, and Khitatrakun (2006)](https://users.ssc.wisc.edu/~aseshadr/Publications/optimality.pdf) supplies the borrowed scale shape, not our complete utility: their aggregation includes an outside household multiplier and exogenous family composition. [Dustmann, Fitzenberger, and Zimmermann (2022)](https://academic.oup.com/ej/article/132/645/1709/6459206) supplies a scale-plus-housing-minimum precedent, but not our jump-only dynamic fertility specification. [de la Croix and Pommeret (2021)](https://perso.uclouvain.be/david.delacroix/pdfpubli/jet21.pdf) supplies an endogenous-fertility precedent for inside consumption scaling with a separate child reward, without our housing sector. These support ingredients; no exact published replication is claimed.

## 4. Pension problem and agreed correction

Read the latest canonical status before implementation: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/CALIBRATION_STATUS.md`.

The old pension formula used a reference worker/retiree age ratio while both payroll tax and pension income were held fixed. The shortcut relied on age masses that ceased to apply after retirement mortality and later population reweighting. Actual payroll revenue and pension outlays therefore did not balance, including in the actual stationary distribution. Existing household-budget, population, and property-tax checks did not test this separate Social Security condition.

The historical review supports the author's recollection that fixed tax and an internally determined balanced pension had previously been communicated. The repair is not an invitation to reopen which fiscal instrument should adjust. Maintain the external payroll tax at 17.9%; pensions adjust to balance actual payroll revenue and retiree exposure in each stationary endpoint and every transition date.

Write B_t for actual taxable payroll in model period units, E_t for actual pension-payment exposure, and b_t for the pension benefit in matching period units. Require

\[
\tau B_t=b_t E_t,\qquad \tau=0.179.
\]

Use the existing actual-household accounting convention. Do not substitute resident-person counts, reference age masses, or unweighted income-state probabilities. Integrate taxable earnings before payroll tax; avoid taxing already net earnings again. In the common-benefit specification, retiree exposure corresponds to eligible household heads. Preserve any explicitly active exposure weighting rather than inventing a new one. Property-tax receipts and rebates remain a separate fiscal contract.

The pension path must enter both backward household valuation and forward policy evaluation consistently with perfect foresight. Updating an accounting spreadsheet or applying pensions only during forward simulation is not the repair. Endogenous pension changes affect saving, wealth, housing and fertility, so old parameters, equilibria and policies need re-evaluation.

Latest inspected repair status: implementation is isolated in `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf`, branch `codex/balanced-social-security`, with source commit a654219c recorded in the repair note. Recheck HEAD and status because the responsible agent may have progressed. Seventy pure tests passed. The compiled smoke retained two failed assertions requiring occupied saving responses to future income. A subsequent diagnostic showed conditional policy responses at unoccupied states and occupied value responses, but did not certify the original suite or a corrected full equilibrium. Preserve that distinction; resolve the test appropriately rather than declaring every occupied household must change its discrete policy.

No corrected jointly cleared production equilibrium or new calibration was certified in the inspected status. The pre-announcement stationary economy, the announced age-reweighted initial state, and the long-run terminal economy are distinct objects requiring separate checks. The terminal stationary boundary is not the 2023 economy.

Detailed source and receipts: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/social_security_repair/README.md`.

## 5. Implementation and calibration contract

Work with the agent's current isolated repair branch; do not overwrite or duplicate its pension work. Preserve old snapshots and outputs. Add a distinct specification identifier and complete source/parameter/target fingerprints for the combined model.

The existing primitive can represent the proposal with `hbar_first_child_jump = hbar_P`, `hbar_child_rooms = 0`, and the child-room-floor specification enabled. Verify this representation consistently in all active household evaluators, feasibility checks, policy replay, stationary initialization, historical transition, terminal solution, and any separately authorized fertility nest. It is acceptable to retain existing parameter field names if their contract is unambiguous: one estimated parenthood requirement, per-child slope fixed exactly to zero. Keep `delta_alpha` and `delta_alpha_jump` at zero as required for constant alpha, including the zero-floor boundary where activation logic may switch branches. Do not leave a redundant slope in search vectors or revive it through an old candidate file, profile, launcher, or checkpoint loader.

Remove one free coordinate if the active search currently estimates both housing parameters. Do not infer the final total dimension from an old note; enumerate the actual active vector, transformations, bounds and external restrictions. Reconcile candidate loading and restart schemas. Keep informative housing targets: a larger-family room gap becomes a stronger test of the remaining mechanisms, not a reason to delete the moment. Report an unreachable target and investigate measurement, implementation and mechanism before proposing any change in targets or weights.

Maintain the agent's currently authorized target contract. This conversation does not adopt the earlier proposed early-target redesign or authorize dropping moments. Preserve the author-selected initial fertility normalization of 2.1, the immediate-announcement/perfect-foresight convention, and the constant post-2023 fertility-preference baseline. Other pre-existing author instructions in the overnight task still apply. Re-solving the initial normalization can change its preference intercept; distinguish this from a fixed-parameter mechanical comparison.

Warm starts from old calibrated parameters are useful. Old value functions, distributions, prices, pension schedules and policy results are not certified solutions of the revised model. Recompute the required stationary endpoints and jointly clear the dated housing and Social Security conditions before interpreting fitted or counterfactual results as equilibria.

### Concrete implementation locations verified for this handoff

All paths in this list are relative to the isolated worktree `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/`. Recheck line locations and branch status before editing.

- `code/model/intergen_eqscale_seq_optimized/solver.py`, around 2230–2290: shared floor, consumption-share and scale construction. Existing inputs already express jump-only preferences; a new utility kernel is not intrinsically required.
- `code/model/intergen_eqscale_seq_optimized/parameters.py`, around 531: nonnegative floor and housing-cap validation. Preserve a feasible total parenthood requirement.
- `code/model/intergen_eqscale_seq_optimized/kernels.py`: scalar and compiled renter/owner evaluations consume shared inputs. Verify agreement rather than modifying one path in isolation.
- `code/model/intergen_eqscale_seq_optimized/e5f_floor_profile.py`, around 18: currently appends a positively bounded, log-transformed `hbar_child_rooms` to the free domain. Remove it from the revised free domain and fix zero; zero cannot be supplied as a free log-transformed coordinate. Update “floor per child” metadata.
- `code/model/intergen_eqscale_seq_optimized/e5f_income_entry_profile.py`: the distinct `first_birth_fixed_cost` is not removed by this housing decision. Preserve it absent separate author instruction.
- `code/model/tools/run_e5f_transition_calibration.py`, `configure_first_child_room_jump`, around 573: appends the jump parameter with the old incremental-jump bounds `[0, 0.5]`. Reconcile total-parenthood-floor bounds explicitly; the initial sum of old parameters may exceed that old upper bound.
- `code/model/tools/e5f_matched_pf_moments.py`, around 18–56: `PARAMETER_NAMES` currently includes both housing terms, and the observer calls `require_identified(11)`. Update searched-coordinate reporting and dimensionality while retaining fixed-zero slope consistency checks.
- `code/model/intergen_eqscale_seq_optimized/e5_profile.py`: inherited twelve-target system. No row deletion is authorized by removing a parameter.
- `code/model/tools/run_e5f_matched_pf_history.py`: the historical/person-tail path must receive the revised shared household parameters and the same dated fiscal paths in both passes.

For the inherited eleven-coordinate transition search, this change yields ten searched coordinates, all else equal. The initial child-preference intercept remains separately normalized to 2.1. A different already-authorized calibration design needs its own count; do not blindly substitute ten everywhere. Identification also requires informative variation and adequate local rank, not only more rows than coordinates.

## 6. Tests and overnight execution

Before the long run, verify directly: zero housing floor when m=0; the same positive floor for every m>=1; increasing scale with the stated values; CRRA algebra and domain; correct treatment of current dependents; preserved first-child initialization; and consistent owner housing-service units. Test both relevant value and policy paths. Verify that metadata and every loaded candidate enforce a zero per-child slope.

For pensions, verify actual payroll and exposure aggregation, time units, endpoint balance, dated balance, anticipation through both passes, and fresh replay of the converged solution. Keep housing-market and Social Security residual gates separate; passing one cannot conceal failure of the other. Preserve existing numerical tolerances and budget/accounting checks.

Use Torch for the overnight computation under the project's long-run rules. Estimate solve count and wall-clock cost using the corrected solver, smoke-test the exact search/collection loop, and set a fixed overnight budget and stopping conditions. Keep checkpoints, a latest-completed summary, a best-so-far summary, and progress at least every case or five minutes. If no checkpoint or heartbeat appears for thirty minutes, investigate. Do not start a second uncoordinated search from this handoff.

If feasibility or equilibrium convergence blocks a production search, return a concrete diagnosis and the best verified checkpoint; do not silently relax gates, substitute an old endpoint, or present a fixed-price diagnostic as a calibrated equilibrium. Perform the already-authorized policy and other reruns after the necessary benchmark checks, within the stated time budget.

## 7. Morning deliverable

Report the exact utility and pension closure actually run; source/specification identifiers; complete free/fixed parameter inventory and bounds; every active target's value, model counterpart, gap, weight and loss contribution; separate housing and Social Security residuals by date; endpoint and horizon checks; and completed versus incomplete objectives. Include the stable diagnostic figure packet. Do not report only a scalar loss or selected successful moments.

Give particular attention to family-size shares and income gradients, housing differences among parents, and consumption/saving behavior. Declining marginal child costs are intentional. Their quantitative implications must be assessed, not “corrected” by silently reinstating the removed slope. Distinguish targeted fit from untargeted validation and diagnostic runs from accepted benchmark/policy results.

No new model code, calibration, policy, target contract, or job was changed in the conversation preparing this handoff. The author-controlled manuscript under `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_draft/` remains read-only.
