# September 5: independent review of the simple allocation argument

Read-only reviewer report. The lead corrected two transcription errors: the
first display now defines old housing value without including the separate
real cost K, consistently with the derivative below; one misspelled local path
was repaired. No substantive reviewer conclusion was changed. This review
preceded the new mixed-tenure equilibrium constructions. Those constructions
are checked separately in `simple_assessment_new_claims_review.md`. Line
references describe the reviewed checkpoint, not subsequent layout edits.

The smallest rigorous result is valid, but it is conditional on an eligible owner–owner pair; the rental cap is not a premise of the local proof.

**Local proposition and proof check — proved.** Fix fertility, tenure, cohort sizes, all other allocations, future occupancy, estates, and external claims. Let a young owner \(i\) receive \(\epsilon\) housing from an old owner \(j\), with \(x_i=c_i-\chi n_i\), \(s_i=h_i-\kappa n_i\), and cost \(C(0)=0\), \(C'(0+)=K\). If

\[
MV_i^Y\equiv \frac{\alpha x_i}{s_i}>MV_j^O+K,\qquad
MV_j^O\equiv\frac{\gamma c_j^2}{h_j^2},
\]

and \(h_i<h_O^{\max}\), \(h_j^2>0\), then a sufficiently small transfer is a Pareto improvement under the permitted direct allocation that relaxes private finance.

Exactly compensate the old household with

\[
D_j(\epsilon)=c_j^2\!\left[\left(\frac{h_j^2}{h_j^2-\epsilon}\right)^\gamma-1\right].
\]

Give \(D_j\) to \(j\), reduce \(i\)’s current consumption by \(D_j+C\), and move \(\epsilon\) housing from \(j\) to \(i\). Current goods and housing exactly balance. The old household is indifferent, while

\[
\Delta U_i=
\log\!\left(1-\frac{D_j(\epsilon)+C(\epsilon)}{x_i}\right)
+\alpha\log\!\left(1+\frac{\epsilon}{s_i}\right),
\quad
\Delta U_i'(0+)=
\frac{MV_i^Y-MV_j^O-K}{x_i}>0.
\]

This derivation is correct and does not require a saving Euler equation. See the restricted allocation problem and proof in [first_best_review.md](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/simplified_olg_amendments/first_best_review.md:65) and [first_best_review.md](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/simplified_olg_amendments/first_best_review.md:103).

Financial feasibility is also correctly separated from real feasibility. The title/claim settlement leaves the young household’s future wealth plus liquidation value unchanged:

\[
\Delta a_i'=\phi P_t\epsilon-qP_{t+1}\epsilon,\qquad
\Delta(a_{i,t+1}+P_{t+1}H_i)=0.
\]

The old seller replaces the estate title with a bond; the consolidated external-asset change is zero. Thus no creditor loses and no foreign resources are created, provided direct current-good transfers and title/financial settlement are admitted. [first_best_review.md](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/simplified_olg_amendments/first_best_review.md:157)–[185]. The actual limitation is institutional, not algebraic: this is not an equilibrium reached while retaining the private down-payment/nonnegative-saving restrictions, nor a complete whole-economy first best.

**Does purchase rationing generate the MRS gap? — yes, conditionally.** In the zero-capital-gains core, the original dated owner problem is [amendment proposal](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex:110)–[122]. Let \(\lambda_i=1/x_i\) be the current-goods multiplier and \(\eta_i>0\) the down-payment multiplier. Its saving and housing FOCs imply

\[
MV_i^Y=\frac{\alpha x_i}{s_i}
=qr_t+\frac{\eta_i}{\lambda_i}(1-\phi_t)P_t>qr_t.
\]

For a slack old owner, \(MV_j^O=\gamma c_j^2/h_j^2=qr_t\). Hence strict purchase rationing yields the needed gap, and supports the proposition whenever \(K<MV_i^Y-qr_t\). This is correctly recorded in [first_best_review.md](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/simplified_olg_amendments/first_best_review.md:125).

Indispensable conditions for that corollary are: a strictly positive down-payment multiplier (not merely equality numerically); positive young saving; slack young physical owner cap; slack young old-age retention and estate constraints; fixed tenure/taste; and slack old-owner retention and estate-composition constraints. Positive consumption and space are also needed. Without them, compare actual \(MV_j^O\) and \(MV_i^Y\) directly; the short corollary no longer follows. The dated original section’s gains-tax terms are intentionally absent from the later core, which explicitly sets the gains tax to zero [amendment proposal](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex:80).

**Rental cap.** It is neither necessary nor sufficient for this local owner-to-owner inefficiency result.

- It is **not necessary**: the proposition never reallocates rental housing, and the owner’s credit wedge alone can generate \(MV_i^Y>MV_j^O\).
- It is **not sufficient**: a binding rental cap does not ensure an eligible young owner, a slack old donor, physical owner-cap slack, or a gap exceeding \(K\).
- It is a reason renting is not an unconstrained substitute for space. That is economically useful for the motivating narrative.

A capped renter’s high housing MRS alone does **not** prove inefficiency if the planner must retain \(h\le h_R^{\max}\): giving that renter more housing in the same tenure is infeasible. Turning the renter into an owner is a tenure/title-conversion problem, not this theorem. This is exactly why the direct proof holds rental occupancy fixed [first_best_review.md](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/simplified_olg_amendments/first_best_review.md:146).

At fixed prices, a renter whose cap is slack chooses housing until its consumption-unit MRS equals \(qr_t\); removing the cap therefore removes this particular renter wedge. But it does not mechanically eliminate an owner-side credit wedge: a sufficiently favorable realized ownership taste \(\xi_i\) can keep a household in ownership. In general equilibrium, removing or loosening the cap changes tenure selection, rents, prices, and the set of donors; it cannot establish that the cap is necessary for, or that its presence guarantees, misallocation.

**What OLG contributes.** It provides the economic story: old owners carry title, young households need liquid down-payment cash before current income/rebates arrive, and title/estate claims have dated settlement. The local welfare proof itself is a two-household reallocation with fixed continuation utilities and resources. It does not show “dynamic inefficiency” merely because agents are young and old.

**Recommended claim.** “When an equilibrium contains a liquidity-constrained young owner and a slack old owner with a consumption-unit housing-value gap exceeding real transfer cost, relaxing private financing in a direct allocation permits a local Pareto-improving reallocation of existing housing toward the young, holding fertility fixed.”

**Do not claim.** “Borrowing limits and rental caps make every equilibrium inefficient,” or that a capped renter’s high MRS proves inefficiency under retained physical rental limits.

**Remaining mathematical task.** If the text wants an aggregate statement rather than an explicitly conditional one, prove existence of a positive-measure set of eligible pairs under a stated equilibrium/tenure-selection configuration. That existence is not implied by age, strict purchase rationing alone, or the rental cap.