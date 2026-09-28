# Fertility, housing costs and supply: four useful results

**2007 stationary reference — block0506, September 28 verified export**

A short analytical note, September 28, 2026. All saved preferences, including the child-benefit parameter, remain fixed. The only numerical comparison reuses the completed experiment with a permanent 10% increase in the housing asset price and its implied rent. No model was solved for this note. Identities below describe the full sequential model; signs require the stated conditions.

## 1. A birth responds to the value of success relative to waiting

A household enters the period with age, liquid wealth, income, location, tenure, children ever born and children currently at home. Write this inherited state as \(s\), and distinguish children ever born \(n\) from children at home \(m\). The reference has one pooled housing market, stochastic earnings, survival, sequential births and tenure choice. Each period lasts four years.

Let \(W_1(s;x)\) and \(W_0(s;x)\) be values after a successful birth and after no birth, respectively, at a specified housing-price experiment \(x\). Each value includes subsequent housing/tenure choice, consumption, saving, future income and child transitions, future birth opportunities, survival and bequests. Thus these are full optimized values, not current utility or next-period values alone. Let \(F\) be the fixed utility cost of a successful first birth, \(\pi_j\) the age-specific conception probability conditional on trying, and \(\kappa_n\) the attempt taste-shock scale.

The attempt probability \(a\) and realized birth probability \(h\) are

\[
G=W_1-W_0-\mathbf{1}_{n=0}F,\qquad a=\Lambda(\pi_jG/\kappa_n),\qquad h=\pi_j a.
\]

Here \(\Lambda\) is the logistic function. The attempt value is \(\pi_j(W_1-\mathbf{1}_{n=0}F)+(1-\pi_j)W_0\); the wait value is \(W_0\). Conception is realized before current housing and saving choices. The first-birth and subsequent-birth scales are 0.176 and 0.332. There is no simultaneous fertility-tenure nest in this reference.

At a smooth interior state, with conception, preferences and the feasible branches fixed, differentiation gives the exact local identity

\[
\frac{\partial\log h}{\partial\log x}=\frac{(1-a)\pi_j}{\kappa_n}\frac{\partial(W_1-W_0)}{\partial\log x}.
\]

**Conditional sign result:** a higher housing cost reduces birth probability precisely when it reduces the success-versus-wait value gap. Housing needs, saving constraints, continuation values and existing-owner valuation all enter that gap. No global negative sign follows from the model alone. Smaller taste scales amplify a given value-gap change, holding the attempt probability fixed; they do not by themselves rank responses across birth orders. At zero conception, unavailable births, feasibility changes or discrete-policy switches, use levels and one-sided/finite differences rather than this log derivative.

The model allows transitions from zero to one, one to two, and two into a three-or-more bin. It does not separately model fourth births. The completed-fertility observer weights that terminal bin by 3.602; raw birth events and adjusted births used for renewal are therefore different objects. In levels, adjusted births equal raw births plus 0.602 times entry into the top bin (using the unrounded weight in calculations).

<!-- pagebreak -->

## 2. Aggregate responses combine behavior and composition

Let \(\mu_s\) be occupied pre-choice household mass and \(B=\sum_s\mu_s h_s\) raw births per period. For first births, restrict the sum to \(n=0\); for subsequent births, restrict it to the relevant states. At positive masses and probabilities,

\[
\frac{d\log B}{d\log x}=\sum_s\omega_s\frac{d\log h_s}{d\log x}+\sum_s\omega_s\frac{d\log\mu_s}{d\log x},\qquad \omega_s=\frac{\mu_s h_s}{B}.
\]

This is a product-rule identity: occupied masses weight level responses, and birth-contribution shares weight elasticities. Use the level product rule at zero cells. For births per household, subtract the elasticity of total household mass. On impact, inherited pre-choice mass is fixed, so the composition term is zero, even though housing transactions and post-choice wealth change. A recomputed cohort allows its occupied distribution to change. A transition additionally carries realized births into future entry and clears dated markets under an explicit expectations path.

Completed fertility \(C=\sum_n v_n p_n\) is a terminal cohort stock, where \(p_n\) is the terminal child-count distribution and \(v=(0,1,2,3.602)\) are the fixed observer weights. Its derivative follows the entire sequence of cohort transitions. It is not the elasticity of current births and has no one-date impact counterpart.

### Measured finite changes: prescribed price and implied rent +10%

The following log-change elasticity is computed from full-precision saved values: \(e=\log(Y_1/Y_0)/\log(1.1)\). It is a finite-change elasticity, not an infinitesimal derivative. The midpoint arc measure, also retained in the calculation JSON, uses proportional differences relative to the two-point means.

<!-- elasticities:start -->
| Outcome / distribution | Reference | Shock | Change % | Log elasticity |
|---|---:|---:|---:|---:|
| Raw birth flow / impact | 0.115 | 0.111 | -4.186 | -0.449 |
| First-birth flow / impact | 0.050 | 0.046 | -8.388 | -0.919 |
| Second-birth flow / impact | 0.042 | 0.041 | -1.138 | -0.120 |
| Entry into 3+ flow / impact | 0.024 | 0.023 | -0.632 | -0.067 |
| Rooms per household / impact | 5.848 | 5.664 | -3.139 | -0.335 |
| Raw birth flow / cohort | 0.115 | 0.109 | -5.130 | -0.553 |
| Rooms per household / cohort | 5.848 | 5.453 | -6.751 | -0.733 |
| Completed fertility / cohort | 2.100 | 1.988 | -5.328 | -0.574 |
| Renewal-adjusted birth flow / impact | 0.130 | 0.125 | -3.796 | -0.406 |
<!-- elasticities:end -->

Birth flows are events per household per four-year period, not conditional hazards. The first/second/top-bin rows distinguish birth order without equating the top bin to exactly three lifetime children. The immediate policies anticipate permanently higher prescribed prices; these are not myopic responses. The cohort calculation keeps entry endowments and fiscal inputs fixed and uses normalized entry. It does not enforce birth-derived demographic renewal.

The experiment changes only the asset price and its implied rent. Preferences, earnings, interest, taxes, pension, survival, entry endowments, credit rules, housing menus and supply primitives are unchanged. The immediate decline in births is mainly first births. The cohort's replacement shortfall is 5.328%; excess housing supply is 0.545 on impact and 0.757 in the cohort. These residuals rule out interpreting either column as a cleared equilibrium.

<!-- pagebreak -->

## 3. Supply elasticity governs price transmission under a stated closure

Let \(Q\) be the asset price per room, \(r\) rent per room, \(D(Q,z)\) occupied aggregate housing demand, and \(H^S=AQ^\eta\) physical supply. The scale \(A\) shifts supply; \(z\) shifts demand. Hold population and fiscal inputs fixed, and specify the same continuation-price beliefs when differentiating demand and fertility. Define \(\epsilon_D=\partial\log D/\partial\log Q\) and \(\epsilon_{D,z}=\partial\log D/\partial\log z\). At a differentiable clearing point,

\[
(\eta-\epsilon_D)d\log Q=\epsilon_{D,z}d\log z-d\log A.
\]

For a pure supply shift with no other direct household effect,

\[
\frac{d\log B}{d\log A}=-\frac{\epsilon_{B,Q}}{\eta-\epsilon_D}.
\]

These are local implicit-differentiation identities when the denominator is nonzero. If demand slopes down, a supply expansion lowers prices; if fertility also falls with prices, it raises births. Holding demand and fertility derivatives fixed, larger supply elasticity dampens price and fertility responses to a given demand shift or proportional intercept shift. Fixed physical stock means \(\eta=0\) and \(A=\bar H\), not that households cannot change housing. The saved supply elasticity is 0.630. Neither the local housing-demand derivative nor the equilibrium fertility-to-supply derivative has been measured here. Substituting the 10% secants would only give a closure-dependent approximation, so no such numerical supply elasticity is reported.

The reference stationary rent map is \(r=uQ\), with \(u=i+\delta+\tau_H\), where \(i\) is the period interest rate, \(\delta\) depreciation and \(\tau_H\) the property-tax rate. Thus the saved 10% experiment changes asset prices and rents together: it does not identify a rent-only elasticity. The dated perfect-foresight implementation instead uses

\[
r_t=(1+i+\delta+\tau_H)Q_t-Q_{t+1}.
\]

Expected capital gains therefore matter for dated rent; current asset prices alone do not determine renter costs. An owner also has a valuation channel through the inherited house and net sale proceeds. Holding pre-choice states fixed preserves their physical house and liquid position, not their total market-valued wealth. The saved increase in post-transaction liquid assets on impact is consistent with this accounting, but does not isolate a causal wealth effect. A multi-date equilibrium requires derivatives of the entire price and pension paths, not this scalar formula.

## 4. Closed stationary renewal restricts the endpoint

Let \(b(Q,p)\) be adjusted births per normalized household, \(e(Q,p)\) adult entry, \(d(Q,p)\) rooms demanded, \(Y_W(Q,p)\) gross worker earnings and \(N_R(Q,p)\) retiree exposure. The pension is \(p\), the fixed payroll tax is \(\tau\), and total household mass is \(N\). A positive closed stationary endpoint must satisfy

\[
\frac{b(Q,p)}{2.1e(Q,p)}=1,\qquad \tau Y_W(Q,p)=pN_R(Q,p),\qquad Nd(Q,p)=H^S(Q).
\]

The first equation follows from the adopted birth-to-household conversion: half of adjusted births enter after 16 years and half after 20, each divided by 2.1. At stationarity the lags no longer change the entry flow. With fixed survival and entry law, the adjusted lifetime births per entrant must equal replacement. This is an equilibrium restriction, not permission to renormalize fertility preferences after a shock. An unrestricted long-run completed-fertility elasticity across these stationary endpoints is consequently the wrong object; prices, pensions, population levels, transition births and age composition are the useful margins.

<!-- pagebreak -->

### A pure supply-scaling candidate

Under the present one-market closure, earnings and entry endowments are per household, there is no outside entry or fixed aggregate grant, and the provisional estate ledger scales with household mass. Multiplying the supply intercept by \(k\) therefore admits the following algebraic candidate: leave prices, pension, policies and the normalized distribution unchanged; multiply population, birth/entry queues, deaths, housing demand and fiscal flows by \(k\). Both sides of all three stationary equations scale consistently.

Along this candidate family, price, per-household birth-flow and completed-fertility elasticities with respect to the supply intercept are zero; population and aggregate birth-flow elasticities are one. These are conditional scaling identities, not measured policy results or claims of uniqueness, stability or attainability from the inherited population. A fixed-stock change has the same scale logic if only its level changes and all stated homogeneity conditions hold. The reference's small renewal residual persists in relative terms. Native queue/one-step scaling checks remain necessary. A positive short-run birth response and a zero stationary per-household response are compatible because population can adjust over time.

## Evidence, limits and reproduction

The authenticated reference is the primary export under `output/model/fertility_identification_20260928/resume_v1/selected_export/primary/`. Common-primary reference loss is 19.581; this is an identity, not a new calibration result. The complete 14-row target-fit table, 31-row parameter/restriction table and all 17 standard plots are retained in the existing [30-page economic packet](../../pdf/fixed_reference_economics_block0506.pdf), with [full-precision reference tables](../fertility_identification_20260928/resume_v1/selected_export/primary/) and [case receipts/tables](../fixed_reference_economics_20260928/fixed_price_v1/). This analytical note supplements that packet.

The source check uses the frozen Torch solver, not current constructor defaults. Source references: `intergen_eqscale_seq_optimized/solver.py` for sequential fertility and post-birth decisions; `parameters.py` and `child_preferences.py` for saved scales, conception and child utility; `run_e5f_perfect_foresight_transition.py:rents_from_asset_prices` for dated rent. Exact locations and hash receipts are recorded in the calculation receipt. The [transition preparation note](../fixed_reference_transition_20260928/preparation_v1/transition_readiness.md) owns the fiscal, population and estate closures and the remaining implementation work.

Reference manifest SHA256: `147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4`.
Checkpoint SHA256: `b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d`.
Source manifest SHA256: `07d84336a3112b251afe505908113d9c00585b91c34bd50f0dee108435db496d`.
Fixed-price executed driver SHA256: `96d6923a252f57bc4d8c44fd6479b13f48ba217d74edf8ef629d120428b03b44`.

The original experiment's two controls reproduce 113 arrays, all 14 fit rows, all 31 parameter values and all 17 plot hashes. This note rechecks compact receipt/table identities and recomputes finite changes without loading a checkpoint. Existing high-wealth ownership/housing patterns, retirement profiles, empirical observer approximations and provisional estate settlement remain limitations. This note adds no global sign, optimum, welfare or efficiency claim.

**Unmeasured local response:** a useful next calculation would hold this exact contract fixed and evaluate asset-price factors 0.990, 0.995, 1.005 and 1.010, with one exact control, five lifecycle solves maximum. Compare central derivatives at both step sizes on identical inherited states, and separately for normalized cohorts; report branch/feasibility switches. Suggested cap: one Torch worker, 20 minutes, no retries; stop on control mismatch, scientific-gate failure or timeout. This plan is not launched and would still not identify rent-only, supply-equilibrium or transition elasticities.

**Reproduction:** `build_note.py` calculates the table/JSON and renders `theory_note.pdf` from this README on Torch only. `launch.sh` supplies a five-minute, one-CPU, four-GiB allocation with zero model solves. `calculation_receipt.json` retains full precision and authenticated inputs; `qa/` contains rendered pages for review. All work is confined to this folder.
