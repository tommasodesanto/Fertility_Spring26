# Dated pensions implied by the observed historical age bridge

**Verdict: the maintained economic transition determines the pensions for 2007, 2011, 2015, 2019 and the inherited 2023 stock before any Bellman solve.** The first four are historical decisions; 2023 is the first person-tail decision, but its inherited household distribution has already received the observed-age bridge. This does not establish a formula for the post-2023 transition. It also does not imply bitwise equality in the current floating-point forward implementation: the numerical qualifications below require an actual dated marginal/fiscal audit before removing root coordinates.

This is a read-only derivation from the isolated `tmp/e5f_matched_pf` source. No household, equilibrium or diagnostic model solve was run, and no code was changed. The approved initial-state object supplies the actual ACS-reweighted 2007 household distribution, not the pre-announcement stationary age weights.

## Exact economic recursion

Let \(M_t(j,z)\) be household-head mass at model age \(j\) and earnings state \(z\), summing wealth, tenure, location and family states. Let \(A_t(j)=\sum_zM_t(j,z)\), and let \(q_t(j,\cdot)=M_t(j,\cdot)/A_t(j)\) be the conditional earnings distribution. There is one location, fixed wages and payroll tax, exogenous age survival \(s_j\), a common row-stochastic earnings matrix \(\Pi_z\), and entrant earnings weights \(\pi\).

Start with \(M_{2007}\) aggregated from `old.initial_state.g_pre`. For each historical advance, before age reweighting,

\[
\widetilde M_{t+4}(0,\cdot)=E_{t+4}\pi,
\qquad
\widetilde M_{t+4}(j+1,\cdot)=s_jM_t(j,\cdot)\Pi_z.
\]

The observed bridge multiplies the entire age cell by one scalar to match its target \(A^{obs}_{t+4}(j)\). Consequently, for positive source age mass and survival,

\[
M_{t+4}(0,\cdot)=A^{obs}_{t+4}(0)\pi,
\qquad
M_{t+4}(j+1,\cdot)=A^{obs}_{t+4}(j+1)
       \bigl[q_t(j,\cdot)\Pi_z\bigr].
\]

These formulas apply at 2011, 2015, 2019 and 2023. Survival and total entry cancel from the conditional earnings recursion because neither selects income states and the age totals are imposed afterward. They remain necessary for the actual demographic accounting and for positive support; zero pre-bridge age mass cannot be repaired by this formula or by the existing bridge.

The target age mass is the initial household scale times the HH-3 total restricted to model ages relative to 2007, times that year's national ACS head-age share. It is not merely a fixed-total ACS-share normalization after 2007. Total household scale cancels from the pension ratio at each date, but retain it when producing absolute fiscal ledgers and checking distributions.

If the initial conditional earnings distribution equals an invariant \(\pi\), every conditional distribution above remains \(\pi\). The general recursion is preferable for verification because it does not silently impose that stronger identity on an actual saved marginal. With fixed income/survival/age-target primitives, the economic pension prefix is also independent of the nine searched preference, housing and bequest coordinates; caching requires those demographic/income contracts to match.

## Why household choices disappear

Current fertility subtracts and adds the same parent's mass within a fixed age/income cell; a birth does not create a new household head. Current housing choice moves wealth/tenure and leaves age/income unchanged. One market removes location selection. Saving and stochastic child maturation conserve parent mass; the latter uses row-stochastic binomial transitions. The earnings transition applies once when the incumbent advances an age cell. Therefore, economically,

\[
M_t^{pre}=M_t^{after\ feasibility}=M_t^{post\ fertility}=M_t^{current}.
\]

The feasibility projection is not an economic obstacle to this identity: `_censor_entry_dead_mass` moves only wealth, retaining age, income and every other column index. A state with no valid wealth destination is rejected by the existing gate. The projection-size gate must remain; this argument does not authorize additional projection.

The first four queue amounts due are inherited waiting vintages. The current birth choices are appended and are not due during those four advances. Actual entry uses outside flow plus retention times the **adjusted** queue amount. The separate **raw** queue is diagnostic and does not set entry. More generally, positive entry totals cancel from age-0 reweighting because all entrants draw the same \(\pi\); this cancellation does not authorize changing the queue, outside-flow or retention contract.

## PAYGO formula and anticipation

Write \(a_j\) for exogenous age earnings efficiency, \(w\) for the one-market wage, \(J_R\) for the first retired age index, and \(d_z=1+s_z(z-1)\) for the retirement-income multiplier. With period length \(\Delta=4\),

\[
b_t=\frac{\tau\Delta w\sum_{j<J_R,z}a_jzM_t(j,z)}
           {\sum_{j\ge J_R,z}d_zM_t(j,z)},\qquad\tau=0.179.
\]

This is a **four-year base pension**. `fiscal_accounts` multiplies annual gross earnings by four once; `bind_social_security_income` writes \(b_t\) directly into retired base income, and the solver applies \(d_z\) afterward. Property-tax transfers and resident nonheads do not enter the payroll base or benefit exposure. Require finite positive payroll revenue and retiree exposure.

Compute the five benefits from the marginal recursion, assemble the full pension vector, and supply that exact vector to both backward induction and forward replay. Eliminating a fiscal root coordinate does not eliminate its effect on current or earlier household choices: the pension still enters every relevant continuation problem. The preference sequence remains fully announced in 2007, linear through 2023 and constant afterward; changing numerical pension guesses must not alter that sequence.

For a six-date mapping, a permissible starting vector is the five derived values followed by the independently balanced terminal pension as a **2027 numerical guess**. The sixth value is not established by this derivation. The 2023 decision itself can use the fifth derived value: the joined routine transfers the observed-reweighted 2023 household stock into the person branch before any person-demographic advancement.

## Numerical obstacles and safe use

1. `advance_cohort_one_period_markov_income` uses stored tenure probabilities directly; in the sequential solver these are stored as `float32`. Current-choice realization separately renormalizes tenure probabilities. Rounded row sums, branch pruning at `1e-15`, clipping tiny negative fertility mass and scalar cohort mass correction can produce tiny, choice-dependent conditional-income errors. The economic recursion is exact; bitwise identity with current code is not proved.
2. Cohort-total and observed-age gates alone do not bound income-composition errors: opposite errors across earnings states can cancel in a total. Check the actual age-by-income marginal and the actual fiscal residual at every prefix date. Do not reinterpret such a discrepancy as a new economic pension mechanism or relax the gates.
3. The current `solve_social_security_path` requires equal-length price and pension vectors and a square `2*T` Jacobian. Genuine elimination of five pensions needs an explicit reduced adapter: `T` housing residuals plus `T-5` free fiscal residuals, with all `T` fiscal ledgers retained as acceptance checks. Merely reporting an exogenous prefix while allowing the root to change it does not perform this reduction.
4. No extension after the inherited 2023 stock is established here. The person/headship law then owns population advancement and birth-related composition; its pensions remain actual-distribution root objects. No current historical equilibrium, horizon or calibration certification follows from this note.

Practical sequence: first use the five values as starting guesses in the existing full root and audit their marginals on a genuine mapping. If the prefix remains within the unchanged marginal/fiscal gates under a distinct trial, fix it in a reviewed reduced root. If a prefix gate fails, diagnose the forward numerical discrepancy; retain the full root as a diagnostic fallback rather than assert exact numerical exogeneity.

## Four decisive tests

1. **Non-invariant earnings fixture:** start with different income distributions across ages, a nontrivial Markov matrix, positive age survival, heterogeneous wealth/tenure and stochastic child maturation. Compare one actual cohort advance plus observed-age reweight to the derived marginal, including age-0 entrant composition and distinct raw/adjusted queues.
2. **Current-stage invariant:** on saved valid policy states, compare age/income marginals before projection, after projection, after fertility and after current housing choice. Include a state where the existing wealth projection moves positive mass. Retain its existing feasibility gate.
3. **Independent fiscal arithmetic:** hand-compute the four-year payroll base and benefits with a nonzero retirement-income dispersion coefficient; compare both the marginal ledger and `fiscal_accounts` on the full distribution, including common-mass rescaling.
4. **Actual historical prefix replay:** on a full five-date prefix from an already authorized mapping, compare predicted and actual marginals/pensions at all dates; repeat at a distinct trial price/pension path. Confirm identical pension slices reach backward/forward calls and that the inherited 2023 stock is included. This checks the shortcut, not historical equilibrium or the post-2023 formula.

## Source map

All paths below are relative to `tmp/e5f_matched_pf/code/model/`.

- `tools/e5f_approved_initial_state.py:57–74`: ACS2007 reweight, adjusted renewal flow, separate raw queue, four slots and observed future age targets.
- `tools/run_e5f_open_population_transition.py:567–650`: HH-3/ACS age targets and uniform within-age reweighting; `:314–348` renewal law; `:351–368` queue timing; `:678–740` sequential fertility; `:788–825` cohort mass correction; `:828–904` survival/cohort advance and entrant placement.
- `intergen_eqscale_seq_optimized/solver.py:3947–3974`: wealth-only projection; `:3538–3652` current housing realization; `:5281–5426` actual cohort advance with income transition, child aging and numerical pruning; `:2423–2426` sequential tenure-probability storage.
- `tools/run_dynamic_population_transition.py:371–389`: exogenous entrant earnings composition; `:392–421` all-age feasibility projection; `:469–524` ordering of current stages.
- `tools/run_e5f_perfect_foresight_transition.py:508–549`: backward fiscal binding; `:618–697` forward fiscal binding, current observation, queue pop, survival advance, entry and observed reweight, in that order.
- `tools/run_e5f_matched_pf_history.py:80–122`: four historical decisions, observed 2023 inheritance, and the first person-tail decision; `tools/run_e5f_perfect_foresight_person_demography.py:727–799` checks the inherited 2023 head-age identity and evaluates its choices before population advancement.
- `tools/e5f_social_security.py:29–52,78–123`: period-income binding and exact fiscal ledger; `tools/e5f_balanced_history.py:190–238` actual dated accounting and common anticipated pension path; `tools/e5f_social_security_root.py:68–75,103–135,168` current full-vector root dimensions.
