# Recent-parent ownership from realized fertility flows

Prepared 2026-09-11. Read-only design review; no observer implementation, model solve, empirical-target change, weight selection, or SMM promotion.

**A useful passive observer is feasible with the existing saved distribution and policies.** Select actual births into households with no pre-birth dependents, carry that selected mass through the actual post-fertility location and tenure choices, and compare ownership with all currently empty-dependent homes, including former parents. This corrects both defects of the old all-dependent-parents versus lifetime-childless comparison. It is a selected-flow ownership contrast, not a causal birth effect.

Under a synchronized four-year observation convention, this construction has an exact interpretation inside the count model: a current birth into an empty home is the only way to have current dependents while having no surviving child from an earlier model birth date. Equating that model group with ACS oldest resident child under four additionally requires a declared timing and residence mapping. The current annual person bridge does not itself establish that mapping. The earlier contract's characterization of the problem as an “irreducible” mismatch is too strong if it is read as ruling out a passive flow observer or passive history bookkeeping. Its narrower conclusion—that the exact ACS group is not recoverable by a static mask of the current count distribution alone—remains correct.

## Empirical object kept fixed

The audited early target is

\[
\Pr(O=1\mid NCHILD>0,\ ELDCH<4)
-\Pr(O=1\mid NCHILD=0)=0.16289550916123285.
\]

It uses ACS household heads ages 30–55, 2005–2006, HHWT, the existing common sample, 42 MET2013 cities, and `UNITSSTR` 3:10 in both groups. `ELDCH` is the **oldest resident own child**: the treatment group requires every resident own child to be under four. The control includes former parents whose own children no longer reside in the household. The audited metro-bootstrap SE is 0.006079524679464666. These numbers and definitions remain those in the [empirical observer contract](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/observer_contract/observer_contract.json). No precision number is adopted as an SMM weight here.

## Existing model objects and exact timing

Let \(g^-\) be `evaluation.g_pre`, \(g^F\) be `evaluation.g_post_fertility`, and \(g^C\) be `evaluation.g_current`. Their axes are liquid wealth, tenure, location, age, income, lifetime parity \(n\), and child/readiness state \(c\). The first array is already the population after the actual price-feasibility gate. Do not substitute the ungated inherited distribution or gate the selected group separately.

The maintained sequential evaluator applies fertility, then current location and tenure/housing transactions. Its current transport does **not** advance age, saving, income, readiness, or child departure. Those operations occur in the separate cohort advance. Therefore the family-state mask commutes with current location/tenure transport. The policies and maps in this evaluation must remain fixed throughout observation. See [evaluate_period](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/tools/run_dynamic_population_transition.py:469), [PeriodEvaluation](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/tools/run_dynamic_population_transition.py:122), and [realize_current_choices](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/intergen_eqscale_seq_optimized/solver.py:3538).

For \(n>0\), \(c\) is the number of dependent children at home under `independent_count`. For \(n=0\), valid readiness states all describe physically empty homes. With readiness disabled the sole childless state is 0; with it enabled, states 0 and 1 mean unready and settled, respectively. Only the settled state can have a current first birth. The parameter default is readiness disabled, but the observer should use the supplied parameter object and readiness helpers, not hard-code that default. Readiness transitions occur during cohort advance, so current unready mass must not be assigned the settled state's current birth probability. See [readiness definitions](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/intergen_eqscale_seq_optimized/parameters.py:827).

Define the physical-empty mask

\[
E(n,c)=\mathbf 1\{n=0,\ c\in\mathcal R_0(P)\}
       +\mathbf 1\{n>0,\ c=0\},
\]

where \(\mathcal R_0(P)\) is `readiness_childless_states(P)`. Invalid family-state cells carrying positive mass must fail validation; they must not be interpreted as additional readiness states.

## Exact passive algorithm for a synchronized snapshot

The implementation can reuse the actual fertility kernel on selected submass; there is no need to duplicate its behavioral equations or construct a forced-birth branch.

1. Require the supplied maintained sequential, independent-count evaluation, with its owned `fert2_probs`, matching dimensions, finite nonnegative mass, and the correct configured model module. Reject alternative joint/nested choice architectures rather than applying the sequential factorization to them. Do not call a solving/configuration routine to repair a missing input.
2. Form \(g_E^-=E g^-\), retaining every wealth, tenure, location, age and income coordinate. Apply `apply_sequential_fertility(g_E^-, policy.fert_probs, P, policy_continuation_birth_probs(policy,P))` with exactly the saved policy. Because this kernel is linear in input mass at fixed probabilities, this is the realized fertility transition of the empty-home subpopulation.
3. From its output \(g_E^F\), select \(B^F=\mathbf 1\{n>0,c=1\}g_E^F\). Every unit of this mass is an actual birth into a pre-birth empty home. All non-birth mass in \(g_E^F\) remains physically empty. A raw `c>0` mask is wrong when childless readiness is active.
4. Let \(K\) be `realize_current_cross_section` with this evaluation's location probabilities, conditional tenure probabilities or choices, and location/tenure transaction maps. Compute \(B^C=K B^F\). Use the exact same current transport as the realized population, with no policy re-solve and no counterfactual first birth.
5. Set the comparison mass \(C^C=E g^C\). Equivalently, \(C^C=K(Eg^F)=K(Eg_E^F)\), up to numerical scatter tolerances. This includes unready never-parents, settled never-parents without a realized birth, and every formerly parented empty-dependent home, including top-coded parity homes. Households with a realized current birth are absent from this control.
6. Use explicitly supplied uniform annual-age overlap \(a_j=|[18+4j,22+4j)\cap[30,56)|/4\). Weights are 1 at cell starts 30, 34, 38, 42, 46 and 50, and 1/2 at 54; all others are zero. For each selected distribution \(X\), compute \(D_X=\sum_s a_{j(s)}X_s\) and \(N_X=\sum_s a_{j(s)}\mathbf1\{tenure(s)>0\}X_s\). Return \(N_B/D_B-N_C/D_C\) only when both denominators are finite and strictly positive.

The model's fertile ages end in the cell beginning at 42 for the maintained configuration. Consequently the full 30–55 control includes ages at which the model cannot have a current birth, just as the empirical definition does not age-match the groups internally. Do not silently standardize both groups to the same age distribution. That would define a different statistic. The declared age overlap uses actual model age mass; it does not import ACS weights.

For independent verification, at any eligible age \(j\), with fecundity \(\pi_j\), the only positive entries of the tagged birth mass before transport are

\[
B^F_{x,1,1}=\pi_j g^-_{x,0,r_*}p^{(1)}_x,
\qquad
B^F_{x,n+1,1}=\pi_j g^-_{x,n,0}p^{(+)}_{x,n,0},
\quad n=1,\ldots,N-2,
\]

where \(x\) retains all other state coordinates, \(r_*\) is the settled readiness state, and \(N=P.n\_parity\). Use the exact kernel's fertility-age restriction. Its continuation-risk pools are captured from the pre-birth distribution, so a new first birth cannot receive another birth in the same period. The literal top-code representative used in population accounting is not an additional timed birth and must not multiply this household-group mass. See [apply_sequential_fertility](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/tools/run_e5f_open_population_transition.py:673) and [owned continuation probabilities](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/tools/run_dynamic_population_transition.py:109).

The location/tenure step matters economically. A selected household's birth changes its family state before its housing choice; movers enter the destination as renters after the location wealth map. Conditional tenure probabilities are then evaluated at that destination wealth, location and post-birth family state, normalized across tenure choices as in the kernel, and transported through the tenure transaction map. An ownership indicator at the origin, or tenure probabilities read at the origin when the household moves, is not the required outcome. Calling the existing transport preserves this selection and sequencing. The selected flow is not the matched forced-birth/no-birth housing-response estimand used elsewhere in the diagnostic packet.

## What this proves, and where timing still matters

**Synchronized model-time equivalence.** Suppose births are assigned to model observation dates four years apart, observation is after current fertility and housing but before the next child-departure transition, and modeled dependents stand for resident own children. Any surviving child already present in \(g^-\) was born at an earlier date and is at least four. A current birth makes all current children under four if and only if no pre-birth child was present. With at most one explicit birth per household per period, this is exactly \(B^C\). Empty controls are exactly \(C^C\). No additional child-age state, past-parity record, or Bellman solve is required for that model-time statistic.

**The actual calendar bridge is different information.** [advance_person_state_block](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/demographic_transition/four_year_bridge.py:30) spreads a decision-date flow into four equal annual person birth flows in years \(t+1,\ldots,t+4\). [The matched birth-path module](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/tools/e5f_matched_pf_birth_path.py:1) states the same convention. These are aggregate person cohorts, not linked child vintages inside particular households. They do not turn `g_current` into an annual-interview resident-child distribution.

At an interior annual observation, a child born late in the preceding block can still be under four, while some of the current block's births have not happened. For example, a prior-block child born at \(t-1\) is two at \(t+1\), but that household has a positive pre-birth dependent count and is absent from the selected current empty-home birth flow. The current flow can also include a birth assigned to \(t+3\), which is not yet a parent at \(t+1\). Thus this is not generally the whole annual ACS under-four stock or even a universally nested subset of it.

Uniform annual births do allow a boundary alignment: at the end of a four-year birth block, that block's children have ages 0–3 while the prior block's have ages 4–7, under the corresponding integer-year convention. But interpreting the evaluator's post-choice population as that endpoint also requires an explicit head-age assignment and a placement of child departures relative to observation. The saved evaluator calls that population current at the decision date and excludes the subsequent maturation step. A block-end relabeling must not silently bypass this timing difference.

The existing [initial fertility observer](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/tools/e5f_initial_fertility_observer.py:1) already treats uniform birth time and constant post-cell parity as alternative diagnostic projections. The housing flow design should retain an equally explicit projection label. The immediate feasible statistic above corresponds to a synchronized/post-fertility snapshot with constant within-cell housing outcomes, not a claim that the annual bridge uniquely implies period-start births.

**A uniform-birth-time sensitivity is possible but is not just a universal one-half weight.** If, as an additional diagnostic assumption, birth time and household observation time are uniform within a four-year cell and the modeled post-birth housing outcome is held constant, the age-selected post-birth exposure for a birth flow is

\[
\omega_j^B=\frac14\int_{u_L}^{u_R}\frac{u}{4}\,du
           =\frac{u_R^2-u_L^2}{32},
\]

where \([u_L,u_R]\subseteq[0,4]\) is the part of that head-age cell in the empirical age interval. It is 1/2 for a complete cell, but 1/8 for the first two years of a cell. This illustrates why annual-age overlap and time-since-birth exposure are different objects. In this specific 30–55 current-birth contrast, all maintained eligible birth cells overlap the sample in full, so the same one-half factor would cancel from the selected-birth ownership ratio under these extra homogeneity assumptions. It still does not supply ownership of future parents before birth, earlier-block recent-child households, or current empty-home exposure. Those components cannot be silently borrowed from the constant post-cell control. This sensitivity therefore needs a separate complete exposure contract before it can replace the proposed snapshot statistic.

## What additional history can and cannot repair

Two households can have identical current \((n,c)=(1,1)\) and the same economic coordinates while their resident children were born one versus five years ago. An annual oldest-child criterion distinguishes them; the current count array does not. Saved current fertility probabilities do not recover that missing vintage. This establishes the limitation of a single static count distribution, not the impossibility of a passive observer.

A richer **passive** cohort calculation could carry child birth-vintage tags alongside the existing household distribution. Policies would continue to depend only on their present model states. Tagged masses would use the same actual birth probabilities, saving, location, tenure, income and survival kernels, and would aggregate back to the original count distribution. No new decision state or reoptimization is intrinsically necessary.

The current independent child-departure kernel makes this extension particularly feasible. Given recent and older counts \((r,o)\), independent identical departures can be applied separately: \(r'\sim Binomial(r,1-\mu)\), \(o'\sim Binomial(o,1-\mu)\), independently. Their sum has precisely the existing \(Binomial(r+o,1-\mu)\) law. Equivalently, conditional on the total number remaining, surviving vintage labels are hypergeometric. This is a policy-preserving lift of the existing count transition, not an assumption that child ages affect choices. See [child-count transition construction](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/intergen_eqscale_seq_optimized/parameters.py:896) and [cohort advance](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/intergen_eqscale_seq_optimized/solver.py:5281).

Such bookkeeping requires initial vintage masses or a consistent passive lifecycle reconstruction from entry. For annual interviews it also needs a specified within-period birth schedule and departure timing. The aggregate person bridge does not supply household links or a unique within-year residence process. These are additional observer inputs/assumptions; they are not reasons to alter the Bellman state automatically.

Even an internally exact tagged observer would preserve the model's residence limitations: dependent exit is stochastic rather than a literal age-18 cutoff; ACS own children can be resident adults; the model lacks adult-child return or residence states; the model allows at most one explicit birth per four-year period while ACS permits multiple under-four children; and the 42-metro DUE structure sample has no exact counterpart in this national housing-product state space. These differences must remain visible. A vintage tag does not manufacture missing empirical residence or structure choices.

## Proposed verification and disposition

A future default-off adapter can implement the snapshot algorithm without changing the current empirical row, followed by small synthetic checks:

- **Readiness and former parents:** with empty masses 0.2 unready, 0.4 settled and 0.1 former-parent, fecundity one and actual birth probabilities 0.25 and 0.5 for the latter two eligible groups, selected birth mass is 0.15 and current empty mass is 0.55. Unready mass stays in the control. A positive-dependent origin never enters the selected birth group.
- **Selection and tenure:** give those two birth origins different actual probabilities and different post-birth conditional ownership policies; the returned parent ownership must weight realized flows, not equalize origins or use a forced birth. Include a moving owner whose location transaction first creates destination renter mass.
- **Kernel identities:** verify \(Eg_E^F=Eg^F\), \(\sum B^F\) equals the empty-subpopulation birth flow, and \(\sum B^C=\sum B^F\). Check \(B^C\leq g^C\) and the equivalent control constructions within the production transport's pruning/scatter tolerances. Fixed-policy transport is linear apart from numerical pruning; synthetic comparisons can disable pruning explicitly.
- **One birth and age exposure:** a first birth cannot immediately receive a continuation birth; a birth at positive pre-birth dependent count is excluded; top-code weights do not scale selected household mass. Check the 30–55 overlap fractions and zero birth flow outside fertility ages.
- **Failure behavior:** missing continuation probabilities, invalid states, nonfinite mass, nonpositive selected denominator, unsupported choice architecture, or missing timing/proxy opt-in returns an explicit failure/unavailable status. There is no denominator floor and no invented zero.

Suggested descriptive name: `ownership_current_birth_from_empty_dependent_home_minus_current_empty_home_30_55`. Require separate explicit consent to the family-residence proxy and a named `synchronized_post_fertility_snapshot` timing convention, in addition to the existing age-projection choice. Save both groups' denominators and owner numerators, first-versus-continuation contributions, readiness setting, age weights, feasibility projection mass, policy provenance, and all mapping warnings. Always label it diagnostic-only, with actual SMM weight and loss contribution unset.

**Recommended next step for the lead:** approve or revise the explicit snapshot convention, then implement and synthetically verify this passive observer in a separate bounded task. The present [initial housing observer](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/tools/e5f_initial_housing_observer.py:240) and the existing fit table remain unchanged and continue to report the empirical recent-parent row unavailable. This note establishes a feasible, more faithful diagnostic path; it does not certify an exact ACS mapper or fill that numerical row without review.

Verification for this note: direct inspection of the named fertility, transport, readiness, child-departure and annual-birth bridge functions; algebraic reconciliation of masks and flow conservation. No numerical model evaluation or new empirical estimation was run.
