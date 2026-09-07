# Assessment of Fable’s review

Fable 5.1 completed the focused review in about 6.5 minutes at high effort. It found no decisive error in the compensated reallocation, conditional fertility derivative, or transition construction. It treated the numerical bounds as supplied evidence; this was not a replay of the certificate.

My assessment: the algebraic checks agree with the note. Its first three suggestions can improve exposition, but none changes the results. The marginal-value gap can be explained without adding a borrowing multiplier to the main text. The payment can be defined by the total differential already in the appendix. The zero-tax primitive condition is separate from the condition used at the positive-tax reference economy.

I would keep the positive-saving assumption pending a deliberate simplification of the settlement argument. Also, an unexpected intermediary capital loss belongs to the intermediary’s residual claimant under the chosen ownership structure; it cannot simply be assigned to external bondholders. That ownership issue matters for a future welfare analysis of the credit reform, which the note expressly does not claim. Finally, when the buyer’s retention bound loosens, its original future allocation remains feasible; its optimized value need not be exactly unchanged if that bound had bound. The direct-allocation proof needs feasibility of the unchanged allocation, so this imprecision in the reviewer’s wording does not undermine the result.

No revisions were made in response to the review. The author can choose the small explanatory changes during discussion.

---

# Fable’s report, verbatim

I'll assess the three items directly from the note, checking the algebra by hand.

**Verdict.** I find no decisive error. The compensated reallocation, the fertility derivative, and the transition construction are all correct as stated. The issues below are gaps in exposition or scope, not mistakes.

**(1) Compensated reallocation, Proposition 1 and Appendix A.1.**

The compensation in \eqref{eq:amend_direct_compensation} holds $\log c_j^2+\gamma\log h_j^2$ fixed with $e$ fixed, so $D_j'(0)=\gamma c_j^2/h_j^2=MV_j^O$. The buyer's gain $\log(1-D_j/x_i)+\alpha\log(1+\epsilon/s_i)$ has derivative $(MV_i^Y-MV_j^O)/x_i$. Goods and housing sums are unchanged. The settlement in \eqref{eq:amend_compensation}–\eqref{eq:amend_settlement} closes: using $u_t=P_t+q\tau^pP_t-qP_{t+1}$, the seller's receipts $P_t\epsilon+L+q\tau^pP_t\epsilon$ minus $D_j$ minus the bond $qP_{t+1}\epsilon$ are exactly zero, and the buyer's borrowing is $qP_{t+1}\epsilon$, repaid at $t+1$ by selling the extra title. The seller's estate receives $P_{t+1}\epsilon$ from the bond in place of the title, so $e$ is unchanged. The buyer's old-age resources $a_{t+1}+P_{t+1}H$ are unchanged, and the buyer's retention bound loosens, so $V^O_{t+1}$ is unchanged. Total tax on the stock and the rebate are unchanged, so future cohorts face identical prices and allocations. This is a genuine Pareto improvement conditional on the assumed margins. The claim in Section 4 that $MV_i^Y>qr_t$ under a binding down payment is correct but only asserted (see item 1 below).

**(2) Fertility derivative, Proposition 2 and Appendix A.2.**

I re-derived \eqref{eq:amend_owner_resources}: from the young owner's budget and $a_{t+1}$, $qW=w-c-uh$ and the Euler equation gives $qW=kx$, so $(1+k)x+\chi n+uh=w$. The FOC \eqref{eq:amend_fertility_foc} follows. Differentiating with $(1+k)\dd x=-\chi\dd n-p\dd h$ gives exactly \eqref{eq:amend_fertility_sign}, and \eqref{eq:amend_dated_fertility} is also correct. For Proposition 3: the share $g(v)$ in \eqref{eq:amend_unconstrained_share} is the correct Cobb–Douglas housing expenditure share, $g(\alpha/(1+k))=\alpha/(1+k+\alpha)$, and $a_h=(1-q)b/[(1-\phi)(y+b)]$ at stationary prices with $\tau^p=0$. The step "$(1+k)\kappa MV^Y>\alpha\chi$ and $p<MV^Y$ imply a positive numerator" checks out: multiply $(1+k)\kappa x>\chi s$ by $\alpha x$ to get $(1+k)\alpha\kappa x^2>\chi s^2MV^Y>\chi s^2p$. The renter version of \eqref{eq:amend_transition_choices} also reproduces from the capped old-renter problem.

**(3) Transition, Section 6 and Appendix B.**

The construction is a standard saddle-path argument. I checked the example: mean fertility $\frac{11}{21}\cdot\frac34+\frac{10}{21}\cdot\frac{9}{40}=\frac12$; both fertility FOCs hold at the stated allocations; $u=9717/18500$; old-owner housing $\beta\gamma x_O/(qu)$ gives housing clearing equal to $\bar H=68104/68019$ exactly; owner and renter budgets both give the same $w$, consistent with the stated $y$ and $T$; saving is positive in both tenures; the old-owner retention and estate bounds are slack ($h^2\approx0.457<1$, $e=0.64>Ph^2$). The eigenvalue count (four stable including zero, two unstable) matches four boundary restrictions in \eqref{eq:amend_actual_boundary} against two jump variables, and the stable/unstable inversion in \eqref{eq:amend_sequence_inverse} is correct in the weighted norm. At $\delta=0$ the boundary restriction on $M$ is satisfied by the baseline tail because $L(A+PH+TB)$ equals $(\beta\gamma/q)Y_{t-1}\int\pi x^O$ under correct expectations, so the identification of the continuation with the baseline is valid. Given the reported signs, impact fertility and endpoint population follow by integrating the derivatives over $[\phi_0,\phi_0+\delta]$, uniformly in $t_p$ because of \eqref{eq:amend_uniform_baseline}. Identity \eqref{eq:amend_cumulative_fertility} follows from \eqref{eq:amend_transition_population} and convergence of both paths. The endpoint ranking reduces to a ranking of $\bar h^{Y*}+\bar h^{O*}$ via \eqref{eq:amend_population_scale}, which is the simple economics behind the result.

**Problems, in priority order.**

1. **Missing derivation, Section 4 after \eqref{eq:amend_oldmrs}.** The inequality $MV_i^Y>qr_t$ is the entire source of the gap but is only asserted. Smallest repair: from the reduced problem \eqref{eq:amend_reduced_problem}, the first-order condition in $s$ gives $\alpha x/s=u_t+\mu x$ with $\mu>0$ the cap multiplier, so $MV_i^Y-MV_j^O=\mu x_i$. One sentence. (Exposition.)

2. **Loose definition of $p$, Proposition 2.** "Net lifetime payment ... including any accompanying transfers" is only made precise by \eqref{eq:amend_dated_fertility}. Repair: define $p\,\dd h\equiv u\,\dd h+h\,\dd u-\dd\widetilde w$ at fixed $\vartheta$, then state Proposition 2. This also makes the rental case $p=u_t+qu_{t+1}$ immediate. (Exposition.)

3. **Proposition 3 does not cover the example economy.** Its hypotheses require $\tau^p=0$, but the reference economy has $\tau^p>0$. The note never says which route verifies the sign there. Repair: state that in the example $\kappa u/\chi=1.75>\alpha/(1+k)\approx0.24$, so the direct sufficient condition preceding Proposition 3 applies; Proposition 3 is a parametric convenience only. (Exposition.)

4. **Proposition 1 hypothesis "uniformly positive young saving" is unused in the direct-allocation proof.** The planner relaxes finance, so saving matters only for the private settlement in Appendix A.1. Repair: move the saving assumption to the settlement remark or drop it. Harmless as written. (Exposition.)

5. **Unexpected-price losses at $t_p$ are unallocated.** With a surprise $P_{t_p}$, the rental intermediary that bought at $t_p-1$ earns a nonzero return, and the note does not say who absorbs it. This does not affect the population or fertility comparisons, which only use $\mathcal I$ and clearing, but it should be stated (e.g. borne by external bondholders) before any welfare statement about the credit reform is added. (Scope, optional extension.)

**Not problems.** The zero stable root, the negative unstable root, and the large modulus $80.9$ are all compatible with the argument. The claim that a temporary reform cannot move the endpoint is correct given local uniqueness of $Z^*(\lambda)$. The rebate treatment, tax reserve wording, and estate bound $e\geq P_{t+1}h^2$ are internally consistent.
