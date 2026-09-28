You are a senior quantitative macroeconomist and economic theorist advising an economics job-market paper. Work autonomously from the evidence included below. Deliver the analysis itself, not a proposal, a list of questions, or an offer to continue.

## Objective

Develop three or four simple, useful analytical results for the FULL implemented lifecycle model of fertility and housing. The main questions are how fertility responds to housing costs, how housing supply affects that response, and why impact, cohort and stationary responses can differ. Aim for material an advisor could understand and a researcher could turn into a short model-results section. I do not want welfare analysis, optimal policy, a general existence/uniqueness theorem, a literature survey or an unrelated two-good toy model.

A four-page preliminary note is included. Independently check it, improve it, and look for one genuinely useful implication beyond mechanical differentiation if the actual model supports one. Do not assume its conclusions are correct. Equally, do not invent objections or novelty. A standard but informative result with a transparent economic condition is more valuable than an elaborate weak proposition.

You have no access to the local repository or Torch cluster. File contents embedded below are the available evidence. You may perform small calculations from supplied scalars and consult a few primary literature sources if essential. Do not claim to have inspected omitted code, run the structural model or obtained new counterfactuals. Where information is missing, give the strongest conditional result available, identify the missing object, and continue the independent parts without asking me to supply more material.

## Non-negotiable reference and scope

Use exactly this label: **2007 stationary reference — block0506, September 28 verified export**.

The baseline approximates the 2007 economy under deliberately imposed replacement stationarity. Calibration normalized the child-benefit parameter to completed fertility 2.1. Economic counterfactuals freeze ALL saved preferences, including that parameter. No recalibration, target revisions, fertility renormalization, replacement of the reference, or silent changes to timing, earnings, entry, credit or fiscal rules are allowed.

Reference manifest SHA256: 147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4.
Checkpoint SHA256: b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d.
Common-primary reference loss: 19.581 (full precision is in the evidence). It is a reference identity, not a new result or a reason to reopen calibration.

The package contains selected frozen source excerpts, not the entire solver. Their parent files were hash-verified on Torch. Use the source excerpts for timing and formulas, the saved parameter table for actual settings, the case-derived calculation receipt for numerical values, and the transition preparation note for the intended closure and outstanding implementation. Current constructor defaults, inactive optional source branches, historical calibrations and old transition fits must not supersede the saved reference.

## Essential model facts

- Households differ in age, liquid financial position, income, inherited housing/tenure and family state. The current reference is one pooled market; do not invent spatial migration responses. It has stochastic earnings, retirement, survival, saving, tenure/housing choice and sequential fertility.
- Each period is four years. Distinguish children ever born n from children currently at home m. At most one modeled birth event occurs per period. States represent zero, one, two and a three-or-more bin. The terminal fertility observer uses weights (0,1,2,3.602); the precise saved top-bin weight is in the supplied data. Fourth births are not separately modeled.
- Timing is: inherit state; choose whether to try for a birth; realize stochastic conception; choose housing/tenure and saving conditional on the outcome; then apply the next-age, income, survival and child-aging transitions. The reference has no simultaneous fertility-tenure nest and no readiness gate.
- A first-birth utility cost is paid on SUCCESS, not on every attempt. The attempt logit compares the conception-weighted success/failure value with waiting. First/subsequent attempt taste scales are distinct. Full optimized post-fertility values include current choices AND continuation, survival and bequests.
- Current child utility and the equivalence scale depend on m. Housing shares change at the first child at home. The compensation factor uses a FIXED reference rent; do not recompute it at the counterfactual rent. Housing floors and the later-child housing-share loading are zero in this reference. Do not replace these preferences with a different fertility aggregator.
- Asset price Q, period rent r, owner carrying costs and sale proceeds are distinct. At a constant price, r=(i+depreciation+property tax)Q. The dated perfect-foresight rent identity subtracts next period's asset price. Initial owners have a valuation channel. Keeping inherited physical housing and liquid assets fixed does not keep market-valued net worth fixed.
- The financed share is phi, so the purchase down-payment threshold uses (1-phi). Borrowing/feasibility constraints, discrete housing choices and rental caps matter. An envelope formula must respect moving constraints and changes in feasible branches.
- Stationary supply is H0[(u Q)/r_bar]^eta, with fixed user-cost rate u and reference rent r_bar; eta=0.630. The notation A Q^eta in the preliminary note absorbs the fixed normalization into A. Fixed PHYSICAL stock means fixing H, not fixing household housing choices or merely retaining an elastic supply curve.

## What to work out

1. **Household fertility sensitivity.** Verify the attempt/birth formula and its derivative using the implemented ordering and conception mixture. Be explicit about which probability is differentiated, what the price experiment changes now and in continuation, and what is held fixed. Then examine whether an envelope representation or a concise sufficient condition can expose more economics than saying “the success-minus-wait value gap falls.” If it can, identify current renter expenditure, initial-owner valuation, continuation and constraint terms. State exactly what is conditional and what has no general sign. Do not assume the marginal child increases housing demand in every state; do not discard continuation or binding constraints to manufacture a sign. If no useful stronger result follows, say so and retain the useful exact identity.

2. **Aggregation and timing.** Give the occupied-state decomposition into behavior and composition. Distinguish conditional birth hazards, first/subsequent event flows, births per household, aggregate births and completed-cohort fertility. Explain why the inherited-state price experiment and the recomputed normalized cohort differ. The latter is not automatically a demographic steady state. An exact finite-change decomposition can be useful if it adds clarity without expanding the note.

3. **Supply transmission.** Derive the small-shock price/fertility mapping under a clearly stated closure. State which housing-demand elasticity belongs in its denominator and why demand, fertility and supply derivatives must use the same population and expectations assumptions. Distinguish a permanent prescribed-price experiment, a one-date market-clearing calculation, and a genuine dated equilibrium. Explain the limits of using a 10% price secant as a local derivative. Do not report an uncomputed numerical supply elasticity. A compact multi-date implicit system is welcome only if it materially clarifies the economics; do not develop a sequence-space solution project.

4. **Closed renewal and supply scaling.** Independently audit the stationary replacement restriction and the pure supply-scaling candidate. The adopted entry rule maps half of adjusted birth vintages into adult households after 16 years and half after 20 years, each divided by 2.1. Preferences stay fixed. Stationary renewal, PAYGO and housing must all hold. Explain what the restriction does and does not imply for measured completed fertility and for finite transitions. Check the homogeneity conditions behind a larger supply intercept being absorbed by population at unchanged prices and per-household behavior. Distinguish an algebraic candidate family from verified equilibria, uniqueness, stability and reachability. If a zero long-run elasticity is a consequence of closure, say so rather than presenting it as an unrestricted fertility prediction.

Use judgment about which three or four results deserve space. Do not fill the response with near-duplicate formulas merely to satisfy this list.

## Existing numerical evidence: reuse, do not overinterpret

One completed experiment permanently increased the asset price and implied rent by 10%, preserving all other primitives, fiscal inputs and preferences. Two exact controls matched 113 numerical arrays, 14 fit rows, 31 parameter values and all 17 standard plot hashes. The shock was then evaluated on (a) the inherited pre-choice distribution and (b) a separately recomputed normalized cohort. Neither calculation clears housing or certifies stationary renewal or a transition.

The included calculation_receipt.json contains full-precision values and both log-change and midpoint arc elasticities. Use those values, not rounded display cells. Key log-change elasticities are approximately -0.449 for immediate raw births, -0.919 for immediate first-birth flow, -0.335 for immediate housing demand and -0.574 for normalized-cohort completed fertility. All are with respect to the joint asset-price/implied-rent change; none is a rent-only elasticity or a local derivative. Approximately 87.100% of the immediate raw-birth decline is first births. The cohort has a 5.328% replacement shortfall.

Raw events and renewal-adjusted births differ: adjusted births add (top-bin weight minus 3) times the flow into the top bin. Preserve that distinction in equations and the numerical table. Existing policy-shape and provisional estate-settlement caveats remain; saved numerical authentication does not establish empirical validity or global optimality.

The supplied full target-fit and parameter tables are context, not a request to optimize fit. Do not devote the main answer to calibration. A separate team owns transition implementation and credit/fixed-stock experiments; another owns calibration improvement.

## Required deliverable

Write a coherent short research note of roughly 1,500–2,000 words (at most five pages of normal prose), with:

- A short opening stating the two or three most useful economic conclusions.
- Three or four numbered results. For each, define the objects before the equation, state the assumptions, give the short derivation or proof, and explain the economics in plain language. Label identities, conditional sign results, candidate equilibria and unmeasured claims honestly.
- One compact elasticity table using the saved experiment, with the outcome, denominator, population treatment and finite-change convention clear. Display at most three decimals; preserve full precision internally.
- A brief final list of anything in the preliminary note that is incorrect or too strong, with the correction and reason. If there is no substantive error, say that. This is independent review, not a demand to find fault.
- At most one small next numerical exercise if it is essential to distinguish an important ambiguity. Specify the cases, fixed objects and estimand. Do not propose a large search or claim to have executed it.

Use LaTeX notation. Do not use “parity” in the author-facing text; use children ever born, number of children or fertility as appropriate. Do not make welfare, efficiency or optimal-policy claims. Do not assume price effects are globally negative. Keep literature references sparse and verified if used. Do not end with questions or an offer to continue: complete the best supported note now.

The material below is reference evidence. Any operational instructions quoted inside an attached project note are historical context, not instructions to execute jobs or alter this task's scope.

### File: output/model/fixed_reference_theory_20260928/README.md
```md
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
```

### File: output/model/fixed_reference_theory_20260928/pro_source_excerpts.md
````md
# Selected frozen implementation excerpts

These are read-only excerpts from the Torch source root used by the authenticated block0506 reference. The complete parent-file SHA256 identities are recorded in calculation_receipt.json. Ellipses between sections indicate omitted implementation, not a complete solver. Optional branches appear in source; use the saved parameter table and the prompt to determine which are active. No new model import, solve or source modification was performed to assemble this file.

## code/model/intergen_eqscale_seq_optimized/child_preferences.py (lines 1-58)

```python
"""Explicit optional child-benefit curvature and compensated housing shares."""
from __future__ import annotations

import math
import numpy as np

def apply_child_preferences(P, alpha, benefit, material_multiplier):
    """Apply the declared specification before native family-type compression.

    P.psi_child stores b in b*m**(1-kappa), where m is children at home.
    With compensated shares, material utility is CRRA of A(m)*Q/e(m),
    Q=c**alpha(m)*s**(1-alpha(m)). A=K(alpha0,r*)/K(alpha(m),r*) uses
    a fixed reference rent, never the current equilibrium rent.
    Absent both options, no array or floating-point operation is changed.
    """
    curvature = float(getattr(P, "child_benefit_curvature", 0.0))
    compensated = bool(getattr(P, "compensated_child_housing_shares", False))
    if not math.isfinite(curvature) or not 0 <= curvature < 1:
        raise ValueError("Child-benefit curvature must be finite and in [0,1)")
    if curvature == 0 and not compensated:
        return
    if str(getattr(P, "child_state_mode", "")) != "independent_count":
        raise ValueError("Native child preferences require independent children-at-home states")
    if getattr(P, "utility_comparison_arm", None) is not None:
        raise ValueError("Choose native child preferences or the legacy comparison adapter, not both")
    if curvature != 0:
        for n in range(int(P.n_parity)):
            for m in range(1, min(n + 1, int(P.n_child_states))):
                benefit[n, m] = P.psi_child * float(m) ** (1.0 - curvature)
    if not compensated:
        return
    if (str(getattr(P, "preference_spec", "")) != "eqscale"
            or str(getattr(P, "eqscale_form", "")) != "power"
            or bool(getattr(P, "child_room_floor", False))
            or float(getattr(P, "hbar_first_child_jump", 0.0)) != 0
            or float(getattr(P, "hbar_child_rooms", 0.0)) != 0
            or float(P.delta_alpha) != 0):
        raise ValueError("Compensated first-child shares require power equivalence scale, no floor and no later loading")
    base_alpha = float(P.alpha_cons)
    loading = float(P.delta_alpha_jump)
    reference_rent = float(getattr(P, "utility_reference_rent", math.nan))
    if (not math.isfinite(reference_rent) or reference_rent <= 0
            or not math.isfinite(base_alpha) or not math.isfinite(loading)
            or loading < 0 or not .05 <= base_alpha - loading <= base_alpha <= .95):
        raise ValueError("Explicit positive reference rent and unclipped interior housing shares required")
    log_k = alpha * np.log(alpha) + (1 - alpha) * np.log((1 - alpha) / reference_rent)
    log_k0 = (base_alpha * math.log(base_alpha)
              + (1 - base_alpha) * math.log((1 - base_alpha) / reference_rent))
    factors = np.where(alpha == base_alpha, 1.0, np.exp(log_k0 - log_k))
    material_multiplier *= factors ** (1.0 - float(P.sigma))
```

## code/model/intergen_eqscale_seq_optimized/solver.py (lines 137-195)

```python
def renter_borrowing_floor(P: SimpleNamespace, b: Any, j: int) -> np.ndarray:
    """Renter floor; all renter debt is unsecured."""

    return debt_rule_at_age(P, b, j)

def native_due_owner_floor(b, collateral_floor, *, death_floor=-np.inf):
    """DUE stayer rule in asset units: service interest, do not grow excess debt.

    The separate net-estate floor preserves no negative estates at possible
    death; it can require repayment following a sufficiently large price fall.
    """
    return np.maximum(np.minimum(np.asarray(b, dtype=float),
                                 np.asarray(collateral_floor, dtype=float)), death_floor)

def native_due_death_floor(P, j, price, house):
    death_possible = (j == int(P.J) - 1 or
        (bool(getattr(P, "use_age_survival", False)) and float(P.survival_probs[j]) < 1.0))
    return -(1.0 - float(P.psi)) * float(price) * float(house) if death_possible else -np.inf

def owner_borrowing_floor(
    P: SimpleNamespace,
    b: Any,
    collateral_floor: Any,
    j: int,
    *,
    stay_on: bool = False,
    stay_orig: bool = False,
    amort: float = 0.0,
) -> np.ndarray:
    """Owner floor after separating secured from unsecured debt.

    Prices and therefore ``collateral_floor`` are those of the current solver
    iterate.  This convention is inert in stationary equilibrium.
    With ``stay_on``, the stayer rule replaces the taper/line rollover: debt
    may not rise (``stay_orig``), and with ``amort`` it must fall by at least
    that share; with no debt the collateral floor applies under
    ``stay_orig``.  Defaults reproduce the legacy floor bit for bit.
    """

    b_arr = np.asarray(b, dtype=float)
    bf_arr = effective_owner_collateral_floor(P, collateral_floor, j)
    if bool(getattr(P, "native_purchase_income", False)):
        if stay_on:
            if bool(getattr(P, "native_due_stayer_credit", False)):
                return native_due_owner_floor(b_arr, bf_arr)
            raise ValueError("native purchase-income floor does not combine stayer rules")
        return np.broadcast_to(bf_arr, np.broadcast_shapes(b_arr.shape, bf_arr.shape)).copy()
    if not stay_on:
        current_unsecured = b_arr - bf_arr
        return bf_arr + debt_rule_at_age(P, current_unsecured, j)
    standard = bf_arr + debt_rule_at_age(P, b_arr - bf_arr, j)
    amort_floor = b_arr * (1.0 - float(amort))
    if stay_orig:
        return np.where(b_arr < 0.0, np.maximum(amort_floor, b_arr), bf_arr)
    return np.where(b_arr < 0.0, np.maximum(amort_floor, standard), standard)

```

## code/model/intergen_eqscale_seq_optimized/solver.py (lines 2549-2618)

```python
    for nn in range(P.n_parity):
        for cs in range(P.n_child_states):
            if independent_child_maturation_active(P):
                nk = cs if cs <= nn else 0
                kp = nk > 0
            else:
                nk = nn
                kp = (cs >= 1) and (cs < csm1)
            if readiness_gate_active(P) and nn == 0 and cs == 1:
                # E6c reuses the otherwise invalid (childless, cs=1) cell for
                # the settled state. It has childless preferences, not the
                # child-at-home consumption and housing adjustments.
                kp = False
            if kp:
                c_bar[nn, cs] = P.c_bar_0 + P.c_bar_n * nk
                if str(P.child_housing_spec).lower() == "linear_only":
                    h_bar[nn, cs] = P.h_bar_0 + P.h_bar_n * nk
                else:
                    h_bar[nn, cs] = P.h_bar_0 + P.h_bar_jump + P.h_bar_n * nk
                psi_v[nn, cs] = P.psi_child * nk
                g_bar[nn, cs] = g0 + gn * nk
                if str(getattr(P, "preference_spec", "stone_geary")).lower() == "eqscale":
                    c_bar[nn, cs] = 0.0
                    if child_room_floor_active:
                        h_bar[nn, cs] = (
                            float(P.hbar_first_child_jump)
                            + float(P.hbar_child_rooms) * nk
                        )
                        alpha_bar[nn, cs] = P.alpha_cons
                    else:
                        h_bar[nn, cs] = 0.0
                        alpha_bar[nn, cs] = np.clip(
                            P.alpha_cons - (P.delta_alpha_jump + P.delta_alpha * nk), 0.05, 0.95
                        )
                    if eqscale_form == "power":
                        # Imposed Scholz-Seshadri-Khitatrakun (2006, JPE 114(4), p.619;
                        # Citro-Michael 1995) scale relative to a childless couple:
                        #   e(n) = ((2 + 0.7 n)/2)**0.7.
                        # Flow utility is multiplied by escale, u = escale * x**(1-sigma)/(1-sigma),
                        # while per-equivalent CRRA utility is u(x/e) = e**(sigma-1) * x**(1-sigma)/(1-sigma),
                        # so the multiplier is e**(sigma-1); at the baseline sigma = 2 the
                        # multiplier equals the scale itself. n is the literal parity state
                        # (under L4, nn = 3 is the top-coded 3+ bin, scaled at n = 3).
                        escale[nn, cs] = (((2.0 + 0.7 * nk) / 2.0) ** 0.7) ** (float(P.sigma) - 1.0)
                    elif eqscale_form == "sqrt":
                        # Declared robustness alternative: square-root household-size scale,
                        # e(n) = sqrt((2 + n)/2) relative to a childless couple.
                        escale[nn, cs] = (((2.0 + nk) / 2.0) ** 0.5) ** (float(P.sigma) - 1.0)
                    else:
                        escale[nn, cs] = 1.0 + P.gamma_e * nk
            else:
                c_bar[nn, cs] = P.c_bar_0
                h_bar[nn, cs] = P.h_bar_0
                g_bar[nn, cs] = g0
                if str(getattr(P, "preference_spec", "stone_geary")).lower() == "eqscale":
                    c_bar[nn, cs] = 0.0
                    h_bar[nn, cs] = 0.0

    apply_child_preferences(P, alpha_bar, psi_v, escale)
    triples = np.column_stack(
        [
            c_bar.reshape(-1, order="F"),
            h_bar.reshape(-1, order="F"),
            psi_v.reshape(-1, order="F"),
        ]
    )
    unique_triples, type_map = np.unique(triples, axis=0, return_inverse=True)
    birth_dp = np.zeros((P.n_parity, P.n_child_states, 1 + P.n_house, 1 + P.n_house), dtype=bool)
    for nn in range(P.n_parity):
        for cs in range(P.n_child_states):
```

## code/model/intergen_eqscale_seq_optimized/solver.py (lines 3013-3130)

```python
def _tenure_location_stage(
    Vd: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    SD: SimpleNamespace,
    ctx: SimpleNamespace,
    dp_choice: np.ndarray,
    Vd_stay: np.ndarray | None = None,
    bmo_purchase: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, np.ndarray, np.ndarray]:
    """Tenure choice + location logit for savings values ``Vd``.

    Returns ``(VH, tcj, prj_or_None, VI, lpj)`` where ``prj_or_None`` is the
    tenure-probability block when the tenure logit is active and ``None``
    otherwise (the caller then leaves ``tenure_probs`` untouched, as before).
    ``Vd_stay`` supplies the stayer (to == tn) values; ``None`` reads the
    standard ``Vd`` there, exactly as before.
    """
    Nb = len(b_grid)
    I = P.I
    npar = P.n_parity
    ncs = P.n_child_states
    nt = ctx.hcost.shape[1]
    purchase_floor = ctx.bmo if bmo_purchase is None else bmo_purchase
    transaction_support = bool(getattr(P, "native_purchase_income", False))
    birth_entry_grant = SD.birth_entry_grant
    tenure_choice_kappa = max(float(getattr(P, "tenure_choice_kappa", 0.0)), 0.0)
    use_tenure_logit = tenure_choice_kappa > 0.0
    if Vd_stay is None:
        Vd_stay = Vd

    if use_tenure_logit and NUMBA_AVAILABLE and bool(getattr(P, "use_tenure_kernel", True)):
        VH, tcj, prj = tenure_logit_kernel(
            Vd, b_grid, ctx.heq, ctx.hcost, dp_choice, purchase_floor, SD.birth_dp, birth_entry_grant, tenure_choice_kappa, Vd_stay, transaction_support
        )
        prj_full: np.ndarray | None = prj
    elif (not use_tenure_logit) and NUMBA_AVAILABLE and bool(getattr(P, "use_tenure_kernel", True)):
        VH, tcj = tenure_choice_kernel(
            Vd, b_grid, ctx.heq, ctx.hcost, dp_choice, purchase_floor, SD.birth_dp, birth_entry_grant, Vd_stay, False, transaction_support
        )
        prj_full = None
    else:
        VH = np.zeros((Nb, nt, I, npar, ncs))
        tcj = np.zeros((Nb, nt, I, npar, ncs), dtype=np.int16)
        prj_full = np.zeros((Nb, nt, I, npar, ncs, nt)) if use_tenure_logit else None
        for id_ in range(I):
            for to in range(nt):
                sp = ctx.heq[id_, to] if to > 0 else 0.0
                Vopt = np.zeros((Nb, npar, ncs, nt))
                if to == 0:
                    Vopt[:, :, :, 0] = Vd[:, 0, id_, :, :]
                else:
                    ba = np.clip(b_grid + sp, b_grid[0], b_grid[-1])
                    Vopt[:, :, :, 0] = interp_on_grid(b_grid, Vd[:, 0, id_, :, :], ba)
                for tn in range(1, nt):
                    hc = ctx.hcost[id_, tn]
                    Vow = Vd_stay[:, tn, id_, :, :] if to == tn else Vd[:, tn, id_, :, :]
                    if to == tn:
                        Vopt[:, :, :, tn] = Vow
                    elif to == 0:
                        bab = b_grid - hc
                        Vb = interp_on_grid(b_grid, Vow, bab)
                        for nn in range(npar):
                            for cs in range(ncs):
                                dpn = dp_choice[id_, tn, nn, cs]
                                bmn = ctx.bmo[id_, tn, nn, cs]
                                if SD.birth_dp[nn, cs, to, tn]:
                                    bag = np.maximum(bab, bmn)
                                    Vb[:, nn, cs] = interp_vector(b_grid, Vow[:, nn, cs], bag)
                                elif birth_entry_grant[id_, tn, nn, cs] > 0:
                                    gfix = birth_entry_grant[id_, tn, nn, cs]
                                    babg = bab + gfix
                                    Vg = interp_vector(b_grid, Vow[:, nn, cs], babg)
                                    inf_m = ((b_grid + gfix) < dpn) | (babg < bmn)
                                    Vg[inf_m] = -1e10
                                    Vb[:, nn, cs] = Vg
                                else:
                                    inf_m = (b_grid < dpn) | (bab < bmn)
                                    Vb[inf_m, nn, cs] = -1e10
                        Vopt[:, :, :, tn] = Vb
                    else:
                        bar = b_grid + sp - hc
                        Vrs = interp_on_grid(b_grid, Vow, bar)
                        for nn in range(npar):
                            for cs in range(ncs):
                                dpn = dp_choice[id_, tn, nn, cs]
                                bmn = ctx.bmo[id_, tn, nn, cs]
                                dpc = dpn - sp
                                inf_m = (b_grid < dpc) | (bar < bmn)
                                Vrs[inf_m, nn, cs] = -1e10
                        Vopt[:, :, :, tn] = Vrs
                tc = np.argmax(Vopt, axis=3)
                if use_tenure_logit:
                    ls, pr = logsumexp(Vopt / tenure_choice_kappa, axis=3)
                    pr[np.max(Vopt, axis=3) <= DEAD_VALUE_CUTOFF, :] = 0.0
                    VH[:, to, id_, :, :] = tenure_choice_kappa * ls
                    assert prj_full is not None
                    prj_full[:, to, id_, :, :, :] = pr.astype(np.float32)
                else:
                    VH[:, to, id_, :, :] = np.max(Vopt, axis=3)
                tcj[:, to, id_, :, :] = tc

    kl = P.kappa_loc
    if NUMBA_AVAILABLE and bool(getattr(P, "use_loc_kernel", True)):
        VI, lpj = location_logit_kernel(VH, ctx.iidx, ctx.iwt, ctx.loc_shift, kl)
    else:
        VI = np.zeros((Nb, nt, I, npar, ncs))
        lpj = np.zeros((Nb, nt, I, I, npar, ncs))
        for io in range(I):
            for to in range(nt):
                Va = np.zeros((Nb, I, npar, ncs))
                Va[:, io, :, :] = VH[:, to, io, :, :]
                idx = ctx.iidx[:, io, to]
                wt = ctx.iwt[:, io, to]
                for id_ in range(I):
                    if id_ == io:
                        continue
                    Vdst = VH[:, 0, id_, :, :]
```

## code/model/intergen_eqscale_seq_optimized/solver.py (lines 3348-3555)

```python
    for j in range(J - 1, -1, -1):
        in_fert = (j + 1 >= P.A_f_start) and (j + 1 <= P.A_f_end)
        s_next = float(P.debt_taper_weights[j + 1])
        D_next = float(P.debt_caps[j + 1])
        renter_floor = np.maximum(renter_borrowing_floor(P, b_grid, j), b_grid[0])
        for zz, z_value in enumerate(z_grid):
            if j == J - 1:
                Vnr = Vbq
            else:
                Vnr = np.zeros((Nb, nt, I, npar, ncs))
                next_values = V if continuation_V is None else continuation_V
                for znext in range(Nz):
                    transition_weight = Pi_z[zz, znext]
                    if transition_weight > 0.0:
                        Vnr += transition_weight * next_values[
                            :, :, :, j + 1, znext, :, :
                        ]
                if bool(getattr(P, "use_age_survival", False)):
                    survival = float(P.survival_probs[j])
                    Vnr = survival * Vnr + (1.0 - survival) * Vbq
            if natural_credit and j < J - 1:
                survival = float(P.survival_probs[j]) if bool(getattr(P, "use_age_survival", False)) else 1.0
                # Income is the leading axis for the strict support operator.
                dated_next = np.moveaxis(next_values[:, :, :, j + 1, :, :, :], 3, 0)
                Vnr = native_solvency_continuation(dated_next, Pi_z[zz], Vbq, survival)
            Vc = apply_child_aging(Vnr, P, Nb, nt, I, npar, ncs, age_index=j)
            if natural_credit:
                child_bad = apply_child_aging((Vnr <= DEAD_VALUE_CUTOFF).astype(float),
                    P, Nb, nt, I, npar, ncs, age_index=j) > 0.0
                Vc[child_bad] = -1e10
            # Parent-age newborn exemption (m-d): continuation with the
            # birth-period child safe.  At fertile ages the housing/saving +
            # tenure/location stages are re-solved under Vc_ex and the success
            # branch reads VI_ex at birth destinations; the wait branch and
            # all non-destination uses keep VI.  Constant mode skips this
            # (Vc_ex is None) bit for bit.
            Vc_ex = (
                apply_child_aging_exempt(Vnr, P, Nb, nt, I, npar, ncs, j)
                if parent_age_maturation_active(P)
                and independent_child_maturation_active(P)
                else None
            )
            Vd, cd, hd, bd = _savings_stage(
                Vc, P, b_grid, SD, ctx, r_hat, j, float(z_value),
                s_next, D_next, renter_floor,
            )
            if stay_active:
                assert bp_pol_stay is not None
                Vd_s, cd_s, _, bd_s = _savings_stage(
                    Vc, P, b_grid, SD, ctx, r_hat, j, float(z_value),
                    s_next, D_next, renter_floor, stay_floor=True,
                )
                bp_pol_stay[:, :, :, j, zz, :, :] = bd_s
                if c_pol_stay is not None:
                    c_pol_stay[:, :, :, j, zz, :, :] = cd_s
            else:
                Vd_s = Vd

            c_pol[:, :, :, j, zz, :, :] = cd
            hR_pol[:, :, :, j, zz, :, :] = hd
            bp_pol[:, :, :, j, zz, :, :] = bd

            dp_choice = ctx.dp_arr
            bmo_purchase = None
            if purchase_income:
                income_for_purchase = np.array([
                    income_at_state(P, i, j, float(z_value)) for i in range(I)
                ], dtype=float).reshape(I, 1, 1, 1) / Rg
                dp_choice = ctx.dp_arr - income_for_purchase
                bmo_purchase = np.maximum(ctx.bmo - income_for_purchase, b_grid[0])
            if natural_credit:
                dp_choice = np.full_like(ctx.dp_arr, -np.inf)
                bmo_purchase = np.full_like(ctx.bmo, b_grid[0])
            if bool(getattr(P, "use_pti_constraint", False)):
                income_j = np.array([income_at_state(P, i, j, float(z_value)) for i in range(I)], dtype=float)
                dp_choice = pti_adjusted_downpayment(ctx.dp_arr, ctx.hcost, income_j, P, b_grid)

            if joint_active:
                joint_result = joint_nested.bellman_block(
                    Vd, (b_grid, ctx.heq, ctx.hcost, dp_choice, ctx.bmo, SD.birth_dp, SD.birth_entry_grant),
                    P, j, fec, tenure_choice_kernel,
                )
                value, joint_prob, product, wait_prob = joint_result[:4]
                if getattr(P, "fertility_nest_choice", False):
                    joint.failure_probabilities[:, :, :, j, zz] = joint_result[4]
                V[:, :, :, j, zz] = value
                joint.probabilities[:, :, :, j, zz] = joint_prob
                joint.products[:, :, :, j, zz] = product
                joint.wait_probabilities[:, :, :, j, zz] = wait_prob
                # This wait-menu kernel is a fallback only. Every forward call
                # replaces it with exact joint-selected mass for its own pool.
                for product_index in range(nt):
                    tenure_probs[:, :, :, j, zz, :, :, product_index] = np.sum(
                        wait_prob * (product == product_index), axis=-1
                    )
                tenure_choice[:, :, :, j, zz] = (np.argmax(wait_prob, axis=-1)
                    if getattr(P, "fertility_nest_choice", False) else product[..., 0])
                loc_probs[:, :, 0, 0, j, zz] = (value[:, :, 0] > DEAD_VALUE_CUTOFF)
                fert_probs[:, :, :, j, zz, :2] = joint_nested.action_marginals(joint_prob[..., 0, 0, :, :])
                fert_value[:, :, :, j, zz] = value[..., 0, 0]
                for nn in range(1, npar - 1):
                    for cs in range(nn + 1):
                        fert2_probs[:, :, :, j, zz, :, nn - 1, cs] = joint_nested.action_marginals(joint_prob[..., nn, cs, :, :])
                continue

            VH, tcj, prj_full, VI, lpj = _tenure_location_stage(
                Vd, P, b_grid, SD, ctx, dp_choice, Vd_s, bmo_purchase,
            )
            if prj_full is not None:
                tenure_probs[:, :, :, j, zz, :, :, :] = prj_full
            tenure_choice[:, :, :, j, zz, :, :] = tcj
            loc_probs[:, :, :, :, j, zz, :, :] = lpj

            # Exact newborn-exempt success values: re-solve the housing/saving
            # + tenure/location stages under Vc_ex and read the success branch
            # off VI_ex at birth-destination states.  Only at fertile ages in
            # parent-age mode, so constant mode stays bitwise identical.
            VI_ex = None
            if (
                in_fert
                and Vc_ex is not None
                and bool(getattr(P, "sequential_births", False))
            ):
                Vd_ex, _, _, _ = _savings_stage(
                    Vc_ex, P, b_grid, SD, ctx, r_hat, j, float(z_value),
                    s_next, D_next, renter_floor,
                )
                if stay_active:
                    Vd_ex_s, _, _, _ = _savings_stage(
                        Vc_ex, P, b_grid, SD, ctx, r_hat, j, float(z_value),
                        s_next, D_next, renter_floor, stay_floor=True,
                    )
                else:
                    Vd_ex_s = Vd_ex
                _, _, _, VI_ex, _ = _tenure_location_stage(
                    Vd_ex, P, b_grid, SD, ctx, dp_choice, Vd_ex_s, bmo_purchase,
                )

            if in_fert:
                pi_j = float(fec[j])
                if bool(getattr(P, "sequential_births", False)):
                    Vfa = np.empty((Nb, nt, I, 2))
                    settled_cs = readiness_settled_state(P)
                    Vfa[:, :, :, 0] = VI[:, :, :, 0, settled_cs]
                    # Exact exempt success value: re-optimized housing/saving +
                    # tenure/location value with the newborn safe next period.
                    if VI_ex is not None:
                        first_dest = VI_ex[:, :, :, 1, 1]
                    else:
                        first_dest = VI[:, :, :, 1, 1]
                    Vfa[:, :, :, 1] = (
                        pi_j
                        * (
                            first_dest
                            - float(P.first_birth_fixed_cost)
                        )
                        + (1.0 - pi_j) * VI[:, :, :, 0, settled_cs]
                    )
                    lf = Vfa / P.kappa_fert
                    ls, pr = logsumexp(lf, axis=3)
                    pr[np.max(Vfa, axis=3) <= DEAD_VALUE_CUTOFF, :] = 0.0
                    fert_probs[:, :, :, j, zz, :2] = pr
                    fert_value[:, :, :, j, zz] = P.kappa_fert * ls
                    if readiness_gate_active(P):
                        # Unsettled households cannot attempt a first birth.
                        # Settled households retain the existing entry logit.
                        V[:, :, :, j, zz, 0, 0] = VI[:, :, :, 0, 0]
                        V[:, :, :, j, zz, 0, 1] = fert_value[:, :, :, j, zz]
                    else:
                        V[:, :, :, j, zz, 0, 0] = fert_value[:, :, :, j, zz]
                    V[:, :, :, j, zz, 1:, :] = VI[:, :, :, 1:, :]
                    childless_copy_start = 2 if readiness_gate_active(P) else 1
                    V[:, :, :, j, zz, 0, childless_copy_start:] = VI[
                        :, :, :, 0, childless_copy_start:
                    ]
                    # Entry (childless wait/try) keeps kappa_fert; upward attempts
                    # at every parity use the continuation scale when set — margin-specific Gumbel scales on the same sequential choice tree.
                    kf_cont_raw = getattr(P, "kappa_fert_continuation", None)
                    kf_cont = float(P.kappa_fert) if kf_cont_raw is None else float(kf_cont_raw)
                    # Parity nn may try for birth nn+1.  Under the historical
                    # shared clock this is only child state 1.  Under the
                    # repaired specification, cs is the current at-home count
                    # and a birth maps (nn, cs) to (nn+1, cs+1).
                    for nn in range(1, npar - 1):
                        child_states = range(0, nn + 1) if independent_child_maturation_active(P) else (1,)
                        for cs in child_states:
                            destination_cs = birth_destination_child_state(P, cs)
                            V2 = np.empty((Nb, nt, I, 2))
                            V2[:, :, :, 0] = VI[:, :, :, nn, cs]
                            if VI_ex is not None:
                                cont_dest = VI_ex[:, :, :, nn + 1, destination_cs]
                            else:
                                cont_dest = VI[:, :, :, nn + 1, destination_cs]
                            V2[:, :, :, 1] = (
                                pi_j * cont_dest
                                + (1.0 - pi_j) * VI[:, :, :, nn, cs]
                            )
                            l2, p2 = logsumexp(V2 / kf_cont, axis=3)
                            p2[np.max(V2, axis=3) <= DEAD_VALUE_CUTOFF, :] = 0.0
                            if independent_child_maturation_active(P):
                                fert2_probs[:, :, :, j, zz, :, nn - 1, cs] = p2
                            else:
                                fert2_probs[:, :, :, j, zz, :, nn - 1] = p2
                            V[:, :, :, j, zz, nn, cs] = kf_cont * l2
                else:
                    Vfa = np.zeros((Nb, nt, I, npar))
                    Vfa[:, :, :, 0] = VI[:, :, :, 0, 0]
                    if pi_j < 1.0:
```

## code/model/intergen_eqscale_seq_optimized/solver.py (lines 4280-4297)

```python
def eval_renter(bp, Rv, Vbar, b_grid, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta, vinterp=None, es=1.0):
    if vinterp is None:
        vinterp = make_value_interp(b_grid, Vbar, "linear")
    surplus = Rv - dc - bp
    ss = np.maximum(surplus, 1e-10)
    if es != 1.0:
        f = es * Kr * ss ** oms / oms + pc + beta * vinterp(bp)
    else:
        f = Kr * ss ** oms / oms + pc + beta * vinterp(bp)
    cm = surplus > cc
    if np.any(cm):
        ct = np.maximum(Rv[cm] - cb_c - ri * hRmax - bp[cm], 1e-10)
        if es != 1.0:
            f[cm] = es * (ct**alpha * ht_cap_c ** (1 - alpha)) ** oms / oms + pc + beta * vinterp(bp[cm])
        else:
            f[cm] = (ct**alpha * ht_cap_c ** (1 - alpha)) ** oms / oms + pc + beta * vinterp(bp[cm])
    f[surplus <= 1e-10] = -1e10
    return f
```

## code/model/intergen_eqscale_seq_optimized/solver.py (lines 4322-4335)

```python
def eval_owner(bp, Rv, Vbar, b_grid, oc, cb_c, pc, Ko_c, alpha, oms, beta, vinterp=None, es=1.0):
    if vinterp is None:
        vinterp = make_value_interp(b_grid, Vbar, "linear")
    ct_raw = Rv - oc - cb_c - bp
    ct = np.maximum(ct_raw, 1e-10)
    if es != 1.0:
        f = es * Ko_c * ct ** (alpha * oms) / oms + pc + beta * vinterp(bp)
    else:
        f = Ko_c * ct ** (alpha * oms) / oms + pc + beta * vinterp(bp)
    f[ct_raw <= 1e-10] = -1e10
    return f

def build_forward_tenure_transition_maps(
```

## code/model/intergen_eqscale_seq_optimized/parameters.py (lines 885-925)

```python
        jr = int(getattr(P, "J_R", P.J))
        work_mean = float(np.mean(income_profile[: max(jr, 1)]))
        if work_mean > 0:
            income_profile = income_profile / work_mean
    return income_profile

def get_fecundity_by_age(P: SimpleNamespace) -> np.ndarray:
    """Per-period conception probability by age index j (length J).

    omega1 == 0 -> all ones (production behavior, including beyond the
    terminal age): this exact rule is the bitwise-nesting guarantee.
    """
    J = int(P.J)
    w1 = float(getattr(P, "fecundity_omega1", 0.0))
    if w1 == 0.0:
        return np.ones(J, dtype=float)
    w2 = float(getattr(P, "fecundity_omega2", 0.0))
    terminal = float(getattr(P, "fecundity_terminal_age", 45.0))
    ages = float(P.age_start) + np.arange(J, dtype=float) * float(P.da)
    pi = 1.0 - w1 * np.exp(w2 * (ages - float(P.age_start)))
    pi = np.clip(pi, 0.0, 1.0)
    terminal_decay = float(getattr(P, "fecundity_terminal_decay", 0.0))
    if terminal_decay < 0.0 or not np.isfinite(terminal_decay):
        raise ValueError("fecundity_terminal_decay must be finite and nonnegative.")
    if terminal_decay > 0.0:
        tail_start = float(getattr(P, "fecundity_tail_start_age", 40.0))
        if not np.isfinite(tail_start):
            raise ValueError("fecundity_tail_start_age must be finite.")
        pi *= np.exp(-terminal_decay * np.maximum(ages - tail_start, 0.0))
    pi[ages >= terminal] = 0.0
    return pi

def fecundity_active(P: SimpleNamespace) -> bool:
    return float(getattr(P, "fecundity_omega1", 0.0)) != 0.0

def readiness_gate_active(P: SimpleNamespace) -> bool:
    """Whether the default-off E6c childless readiness state is active."""
    return bool(getattr(P, "readiness_gate_enabled", False))
```

## code/model/tools/run_e5f_perfect_foresight_transition.py (lines 358-386)

```python
def rents_from_asset_prices(
    prices: Sequence[float], terminal_price: float, P: SimpleNamespace
) -> np.ndarray:
    """Apply the one-period owner/renter no-arbitrage identity.

    Owners earn next period's asset price and pay depreciation and property
    tax.  Hence r_t + p_{t+1} = (R + delta + tau_H) p_t.  At a constant
    price this reduces exactly to the stationary user-cost identity used by
    the existing model.
    """
    current = np.asarray(prices, dtype=float).reshape(-1)
    if current.size < 1 or np.any(~np.isfinite(current)) or np.any(current <= 0.0):
        raise ValueError("Asset prices must be finite and strictly positive.")
    terminal = float(terminal_price)
    if not math.isfinite(terminal) or terminal <= 0.0:
        raise ValueError("The terminal asset price must be finite and positive.")
    next_prices = np.r_[current[1:], terminal]
    carrying_factor = float(P.R_gross) + float(P.delta) + float(P.tau_H)
    user_cost = float(getattr(P, "user_cost_rate", carrying_factor - 1.0))
    if not math.isclose(user_cost, carrying_factor - 1.0, rel_tol=0.0, abs_tol=2e-14):
        raise ValueError("Stationary user cost disagrees with interest, depreciation and tax")
    # The equivalent expression avoids subtracting two price-sized terms.
    # Constant paths reproduce the stationary Bellman's rent bit for bit.
    rents = user_cost * current + (current - next_prices)
    if np.any(~np.isfinite(rents)) or np.any(rents <= 0.0):
        raise ValueError(
            "The candidate asset-price path implies a nonpositive renter price."
        )
    return rents
```

## code/model/tools/run_e5f_open_population_transition.py (lines 767-805)

```python
def calendar_topcode_birth_accounting(
    g_pre: np.ndarray,
    g_post: np.ndarray,
    explicit_births: float,
    P: SimpleNamespace,
) -> dict[str, float]:
    """Translate the explicit 0/1/2/3+ state into measured child units.

    Entry into the last parity state identifies families reaching the 3+ bin.
    The additional children represented by that bin are used only for aggregate
    population renewal; household choices continue to use the existing state.
    """
    top_state = int(P.n_parity) - 1
    if top_state != 3 or str(getattr(P, "fertility_units", "")).lower() != "literal_topcode":
        return {
            "explicit_birth_children": float(explicit_births),
            "top_bin_entry_birth_flow": 0.0,
            "topcode_adjusted_birth_children": float(explicit_births),
        }
    top_weight = float(getattr(P, "tfr_top_bin_weight", top_state))
    top_before = float(np.sum(g_pre[:, :, :, :, :, top_state, :]))
    top_after = float(np.sum(g_post[:, :, :, :, :, top_state, :]))
    top_entry = top_after - top_before
    if top_entry < -2e-12:
        raise RuntimeError(
            f"Top-bin mass fell during the fertility stage: {top_entry:.3e}"
        )
    top_entry = max(top_entry, 0.0)
    adjusted = float(explicit_births) + (top_weight - top_state) * top_entry
    if adjusted + 1e-14 < float(explicit_births):
        raise RuntimeError("Top-code adjustment reduced the birth flow")
    return {
        "explicit_birth_children": float(explicit_births),
        "top_bin_entry_birth_flow": top_entry,
        "topcode_adjusted_birth_children": adjusted,
    }

def children_at_home_units(distribution: np.ndarray, P: SimpleNamespace) -> float:
```

# Transition closure evidence

Source: output/model/fixed_reference_transition_20260928/preparation_v1/transition_readiness.md. This is the preparation team's current closure assessment, not a completed transition result. The selected sections below are reproduced verbatim.

## Closures and their status

| Object | Classification | Contract for this preparation |
|---|---|---|
| Preferences | Estimated or calibration-normalized, then fixed | All saved preferences; no post-shock normalization. Historical preference shocks are a separate, pending specification. |
| Earnings | Externally estimated, fixed | Approved B15 15-state Markov process and saved age profile; use completed measurement audit, not constructor defaults. |
| Initial population | Empirically normalized reference | Exact saved pre-choice distribution, mass one. No age reweighting, reset or rescaling along a path. |
| Adult entry | Author-fixed reduced-form normalization | Half of each birth vintage enters after 16 years and half after 20; adjusted children map to households at 1/2.1. Preserve both raw and adjusted queues. This is not a new estimate of child survival. |
| Entry wealth/income | Empirically normalized; coupling approximation retained | Full saved conditional entrant distribution and grid. No zero-wealth fallback or frontier censoring of inherited population. |
| Survival | Externally fixed | Exact saved age survival and terminal death; no mortality change. |
| Geography and outside entry | National target scope; closed-path assumption explicit | One pooled market; no spatial migration in I=1. No new immigration, retention parameter, age bridge or old quota defaults. National calibration does not itself validate a population forecast. |
| PAYGO pension | Empirical baseline normalization; fixed tax closure | Hold saved payroll tax at 8.028%; solve each date's equal pension from actual worker/retiree masses. Baseline benefit 0.918 is a starting value, not a fixed transition benefit. |
| Property tax | Externally fixed | Saved annual rate 1.060%, period rate 4.239%, zero household rebate. Do not import old equal rebates. |
| Estates/entry funding | Outstanding substantive settlement; provisional reference ledger retained | Net positive estates fund actual next-period positive entrant assets; residual sink; funding shortage and negative estates fail. Donor utility stays unchanged. Lender counterparties, recipient law and physical settlement remain unresolved. |
| Housing, credit experiment | Estimated intercept, fixed external elasticity | Retain saved absolute supply curve (elasticity 0.630). Do not rescale it by population. |
| Housing, fixed-stock experiment | Author-fixed counterfactual | H equals actual reference supply at the reference equilibrium price, not H0. Prices/rents and individual housing/tenure choices may change. Constant gross stock entails replacement of depreciation; it is not zero gross construction. |
| Credit experiment | Author-adopted comparison; implementation outstanding | Replace artificial purchaser and incumbent debt restrictions by lifetime no-default solvency; keep interest, repayment, prices and death settlement. No arbitrary debt-floor substitute. |
| Terminal population | Endogenous equilibrium object; outstanding | Solve renewal, PAYGO and housing jointly at fixed preferences. A normalized stationary price root is not enough. |

The reference has a small measured renewal discrepancy: entry exceeds potential
birth-derived entry by \(4.889\times10^{-8}\) per reference household. Seed
prehistory from actual saved entry, then let actual births enter the queue.
Report resulting no-shock drift rather than adjusting child benefit or queues.

## Natural solvency and the new steady state

A natural borrowing limit is the most debt a household can repay under every
modeled future event with positive probability. Construct feasible sets backward
at each existing decision node, preserving the order of fertility realizations
and subsequent housing choices. Feasibility must be a separate Boolean object;
large negative utility is not itself proof of infeasibility.

For each chosen saving/tenure branch, require the resulting successor state to
be feasible for every reachable income and family-state outcome. Whenever death
has positive probability, require the inherited timing's net liquidation estate
\(b'+(1-\psi_{sell})q_t h\geq0\). The terminal condition is the same repayment
condition. With positive death risk and no default or life insurance, this
condition can itself exclude unsecured debt even when future earnings are
positive; removing an artificial limit does not authorize unpaid death-state
liabilities. The current model values this estate at today's price after saving;
changing its timing would be an additional economic change.

The current `native_solvency_credit` uses value cutoffs, first feasible grid
nodes and the grid's lower end; it is a prototype, not yet a certified natural
limit. It also conflicts with saved `native_due_stayer_credit=True`. Replacing
the incumbent-owner restriction is part of the authorized credit change, but
must be explicit. Setting financed share to one only removes a down payment;
it does not remove the collateral-based debt limit. Verify the new feasible
sets against terminal/worst-income cases and an expanded/refined debt grid,
without altering the reference entry distribution or relaxing occupied gates.

For a closed positive stationary population, let \(B(q,p)\), \(E(q,p)\) and
\(d(q,p)\) denote adjusted births, new household entry and housing demand per
normalized household, with pension \(p\). At fixed preferences the endpoint
must satisfy
\[
B(q,p)/(2.1E(q,p))=1,\qquad
\tau_{pay}Y_W(q,p)=p N_R(q,p),\qquad
N d(q,p)=H^S(q).
\]
Here \(Y_W\) is gross worker earnings and \(N_R\) is retiree exposure in the
normalized distribution. After solving renewal and PAYGO, housing determines
the level \(N=H^S(q)/d(q,p)\); for fixed stock replace \(H^S(q)\) by \(\bar H\).
One-step native distribution and queue reproduction must then pass. If no
positive root is found, report that failure and the search range; do not force
replacement fertility, rescale entrants or claim nonexistence from a timeout.

## Supply scaling: a candidate endpoint, not a transition

For a pure 10% increase in the housing-supply intercept, a candidate endpoint
keeps prices, pension, household policies and the normalized distribution
unchanged and multiplies population, both entry queues, births, deaths, housing
demand and fiscal flows by 1.1. This follows from the level-linear forward
population operator and the per-household earnings and entry laws; both sides
of PAYGO and the provisional estate-funding ledger scale by the same factor.
The absolute supply schedule also scales by 1.1. The reference's small renewal
residual scales in levels and is unchanged in relative terms.

This algebra applies to the present single-market, exogenous-earnings closure
with no outside entry and no fixed aggregate transfer. A fixed immigration
flow, aggregate fiscal grant, population-dependent wages or amenities, a
non-scaling estate allocation, or normalization of population at each date
would break it. The reviewed estate ledger is homogeneous in its economic
flows; its absolute numerical tolerance is not an economic transfer. Native
one-step and scaling checks are still required, and no uniqueness, stability
or transition claim follows. The economic-analysis chat owns this supply
comparison; no additional supply solve was launched here.

## What each calculation can establish
````

### File: output/model/fixed_reference_theory_20260928/calculation_receipt.json
```json
{
  "reference_label": "2007 stationary reference \u2014 block0506, September 28 verified export",
  "slurm_job": "18743282",
  "model_solves": 0,
  "checkpoint_loaded": false,
  "reference_manifest_sha256": "147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4",
  "checkpoint_sha256": "b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d",
  "inputs": {
    "output/model/fixed_reference_economics_20260928/fixed_price_v1/control_1/receipt.json": "4a912ae4384301a0c3d243596dc23a86d67245567e1336702345752426d48414",
    "output/model/fixed_reference_economics_20260928/fixed_price_v1/control_2/receipt.json": "9d7732b31f15a0c22940dadc4f1a95663fb77e5a010af926e22ce42887d66cb6",
    "output/model/fixed_reference_economics_20260928/fixed_price_v1/price_110/receipt.json": "8d425fbfcf81100d7b386ce9d6df31cd200fcc122f9ad94e23eea5e9a08485a6",
    "output/model/fixed_reference_economics_20260928/fixed_price_v1/comparison.json": "5439f92c85577caf23be0e5847094f4c408c026a52cb9a081d696bb2f7371e67"
  },
  "verified_source_hashes": {
    "code/model/intergen_eqscale_seq_optimized/solver.py": "b637a655a9344b63f4461ee0fa4796c04bd98188477c4e6ace2c48ae0fc8aec1",
    "code/model/intergen_eqscale_seq_optimized/parameters.py": "66f86697c2c58ca3864305bf13dd2be71a008905b2beb573f1a4ebafabef5464",
    "code/model/intergen_eqscale_seq_optimized/child_preferences.py": "69544858b0137e09d37ac9e07ee6470ddc2d5849236b2c3e33d270737164f09d",
    "code/model/tools/run_e5f_perfect_foresight_transition.py": "fc5519e89ef8c31906aa76fdad6a08d703ab427b7d7ada775d04bc0e240c8898",
    "code/model/tools/run_e5f_open_population_transition.py": "af52983e0d6d97f16dff74f2491742825d2cbe506240417d47048773d31f249d"
  },
  "source_locations": {
    "attempts": "solver.py:3486-3550",
    "continuation": "solver.py:3348-3453",
    "forward_timing": "solver.py:5533",
    "dated_rent": "run_e5f_perfect_foresight_transition.py:358-386",
    "adjusted_births": "run_e5f_open_population_transition.py:767"
  },
  "asset_price_ratio": 1.1,
  "metrics": [
    {
      "outcome": "Raw birth flow",
      "distribution": "Inherited states",
      "reference": 0.1154188014126581,
      "shock": 0.11058774873106168,
      "percent_change": -4.1856721976551325,
      "log_change_elasticity": -0.4486189487685127,
      "midpoint_arc_elasticity": -0.44889011512724014
    },
    {
      "outcome": "First-birth flow",
      "distribution": "Inherited states",
      "reference": 0.05016405889154707,
      "shock": 0.04595622439855762,
      "percent_change": -8.388145987322593,
      "log_change_elasticity": -0.9192041421173992,
      "midpoint_arc_elasticity": -0.9193119425801278
    },
    {
      "outcome": "Second-birth flow",
      "distribution": "Inherited states",
      "reference": 0.04164532699159618,
      "shock": 0.04117137921240421,
      "percent_change": -1.138057528729708,
      "log_change_elasticity": -0.12009031582813423,
      "midpoint_arc_elasticity": -0.12017989870925892
    },
    {
      "outcome": "Entry into 3+ flow",
      "distribution": "Inherited states",
      "reference": 0.023609415529510624,
      "shock": 0.02346014512009592,
      "percent_change": -0.6322494905819376,
      "log_change_elasticity": -0.06654658018988893,
      "midpoint_arc_elasticity": -0.06659672523913761
    },
    {
      "outcome": "Rooms per household",
      "distribution": "Inherited states",
      "reference": 5.847942329061136,
      "shock": 5.664363821043602,
      "percent_change": -3.139198331441262,
      "log_change_elasticity": -0.334647070678697,
      "midpoint_arc_elasticity": -0.334871972487732
    },
    {
      "outcome": "Raw birth flow",
      "distribution": "Normalized cohort",
      "reference": 0.1154188014126581,
      "shock": 0.10949738012158954,
      "percent_change": -5.1303784293320165,
      "log_change_elasticity": -0.5525814936512565,
      "midpoint_arc_elasticity": -0.5528719466256156
    },
    {
      "outcome": "First-birth flow",
      "distribution": "Normalized cohort",
      "reference": 0.05016405889154707,
      "shock": 0.04812573775866602,
      "percent_change": -4.063309823648497,
      "log_change_elasticity": -0.43522831963926945,
      "midpoint_arc_elasticity": -0.43549529299396705
    },
    {
      "outcome": "Second-birth flow",
      "distribution": "Normalized cohort",
      "reference": 0.04164532699159618,
      "shock": 0.03939778831543378,
      "percent_change": -5.396856834899966,
      "log_change_elasticity": -0.582094008734459,
      "midpoint_arc_elasticity": -0.5823852158274099
    },
    {
      "outcome": "Entry into 3+ flow",
      "distribution": "Normalized cohort",
      "reference": 0.023609415529510624,
      "shock": 0.02197385404748539,
      "percent_change": -6.927581413359684,
      "log_change_elasticity": -0.7532490393730397,
      "midpoint_arc_elasticity": -0.7534955574986499
    },
    {
      "outcome": "Rooms per household",
      "distribution": "Normalized cohort",
      "reference": 5.847942329061136,
      "shock": 5.453161126072182,
      "percent_change": -6.750771139227963,
      "log_change_elasticity": -0.7333361161570704,
      "midpoint_arc_elasticity": -0.7335925465758191
    },
    {
      "outcome": "Completed fertility",
      "distribution": "Normalized cohort",
      "reference": 2.0999983368073374,
      "shock": 1.988120311708367,
      "percent_change": -5.327529224097405,
      "log_change_elasticity": -0.5744079748226275,
      "midpoint_arc_elasticity": -0.5746992025124807
    },
    {
      "outcome": "Renewal-adjusted birth flow",
      "distribution": "Inherited states",
      "reference": 0.12964015530498035,
      "shock": 0.12471918818584646,
      "percent_change": -3.795866417744753,
      "log_change_elasticity": -0.4060202245813809,
      "midpoint_arc_elasticity": -0.40627683687011323
    },
    {
      "outcome": "Renewal-adjusted birth flow",
      "distribution": "Normalized cohort",
      "reference": 0.12964015530498035,
      "shock": 0.1227335381449392,
      "percent_change": -5.32752922409977,
      "log_change_elasticity": -0.5744079748228896,
      "midpoint_arc_elasticity": -0.5746992025127422
    }
  ],
  "first_birth_share_impact_decline": 0.8709974347867161,
  "pdf_pages": 4,
  "pdf_sha256": "8a9a0fa66ede31ea39180d5bbc2204204b6f04389252878e8b654b197f31ce00",
  "readme_sha256": "d9a1df8debea9122558dcd437c2ddf6e53bcdb5c715a2dbc020768c78f9c9450",
  "status": "calculation_and_render_pass_visual_review_pending"
}
```

### File: output/model/fertility_identification_20260928/resume_v1/selected_export/primary/target_fit.csv
```
moment,target,model,gap,weight,loss_contribution,role
initial_normalization,2.1,2.0999983368073374,-1.6631926627042048e-06,,,normalization
cps_childlessness,0.19827875100684264,0.2010621765888024,0.002783425581959764,35532.3042455214,0.2752850337303754,scored
cps_exactly_one,0.21365532522014702,0.2094219325930911,-0.00423339262705591,26952.820824310795,0.4830380277051848,scored
nchs_mean_age,25.976263860992496,25.93278056306455,-0.043483297927945586,139.82806784479274,0.26438651897923604,scored
nchs_share30,0.2492780130410667,0.22364558977321622,-0.02563242326785048,0.0,0.0,validation
wealth_earnings,6.92658379107299,6.326314855155861,-0.6002689359171294,7.595098472533724,2.736687113167318,scored
bequest_wealth,0.007291023472616158,0.007057418371647596,-0.00023360510096856208,5165289.256198346,0.2818767727196904,scored
old_dispersion,3.51593508651872,3.0685912991959516,-0.44734378732276836,0.0,0.0,validation
mean_rooms,5.729434240102641,5.847942329061135,0.11850808895849418,128.02070205233477,1.7979441381033878,scored
ownership_30_55,0.6762604168538028,0.6548177832031233,-0.021442633650679443,2339.3623724673616,1.075607326077079,scored
first_birth_rooms,1.465,1.6221055857475672,0.15710558574756717,137.5652749002964,3.3954088234135944,scored
family_rooms,0.38509964969278165,0.3526264368889658,-0.032473212803815876,0.0,0.0,validation
recent_parent_ownership,0.12760836356692162,0.11954774603336116,-0.008060617533560466,27055.822957508266,1.7579130016044182,scored
early_fertility,0.8095276384290021,0.5354261243026106,-0.2741015141263915,100.0,7.51316400463804,scored
```

### File: output/model/fertility_identification_20260928/resume_v1/selected_export/primary/parameters.csv
```
parameter,estimate,lower,upper,near_bound,status
H0,6.293507689200028,0.2,80.0,False,free in evening DUE calibration
beta_annual,0.9634753808828738,0.94,0.99,False,free in evening DUE calibration
chi,1.0938718913975698,0.1,5.0,False,free in evening DUE calibration
first_birth_fixed_cost,0.6209410745816288,0.0,8.0,False,free in evening DUE calibration
kappa_fert,0.17561183772819416,0.02,50.0,True,free in evening DUE calibration
kappa_fert_continuation,0.33176505374023935,0.02,50.0,True,free in evening DUE calibration
theta0,0.12451612635999697,0.0,8.0,False,free in evening DUE calibration
delta_alpha_jump,0.13487659325900622,0.0,0.25,False,free in evening DUE calibration
child_benefit_curvature,0.061454973811513665,0.0,0.8,False,free in evening DUE calibration
tenure_choice_kappa,0.012087458119671006,0.001,0.1,False,free in evening DUE calibration
psi_child,0.1355551166583114,,,,normalized to completed fertility 2.1
child_benefit_CRRA_coefficient,0.1272245805140582,,,,derived from normalized one-child benefit
theta1,0.008193084126995582,,,,fixed external restriction
sigma,2.0,,,,fixed
alpha_cons,0.733,,,,fixed CEX childless expenditure share
delta_alpha,0.0,,,,fixed zero later-child loading
h_P,0.0,,,,no housing floor
utility_reference_rent,0.11046592704873838,,,,fixed substantive utility normalization
q_annual,0.020000000000000018,,,,author-retained 2% annual real rate
financed_share,0.8,,,,inherited credit contract
housing_supply_elasticity,0.63,,,,fixed provisional external mapping
payroll_tax,0.08028070961950022,,,,derived from adopted pension ratio
pension_period,0.917784047463731,,,,balanced PAYGO
annual_depreciation,0.01416143718381309,,,,adopted
period_depreciation,0.05545379079326218,,,,compounded
annual_property_tax,0.010598360773872594,,,,adopted
period_property_tax,0.042393443095490375,,,,linear period convention
selling_cost,0.06,,,,retained
rental_cap,6.0,,,,retained provisional
wealth_grid_nodes,160.0,,,,retained exact grid
income_states,15.0,,,,retained B15
```
