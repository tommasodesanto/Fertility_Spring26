# Intergenerational Housing and Fertility: Paper Issues Ledger

This ledger records substantive issues identified while reviewing
`latex/intergenerational_housing_fertility_paper_draft.tex`. It distinguishes
decisions that are accepted for the current draft from open model or exposition
issues. Add items here as the paper review proceeds; do not treat an issue as
resolved merely because the paper and code currently agree.

| ID | Issue | Current paper and code | Status | Required follow-up |
|---|---|---|---|---|
| BQ-1 | Terminal housing valuation | The estate includes the gross house value, \(W^B=b+Ph\), with no sale cost at death. | Accepted provisionally | Retain for the current draft. Revisit whether inherited housing should enter gross or net of liquidation costs in the final model. |
| BQ-2 | Bequest normalization | Bequest utility is not recentered at zero estate wealth. With \(\sigma>1\) and \(\theta_n>0\), this both raises the marginal value of estate wealth for larger families and lowers their terminal utility level at a fixed estate. | Open; high priority | Decide whether to restore the zero-estate normalization. If restored in code, re-estimate before reporting quantitative results. Do not describe the current formulation as an innocuous normalization. |
| BQ-3 | Bequest-to-entry wealth link | Entrant wealth is drawn from an exogenous distribution, while estates provide warm-glow utility but are not assigned to descendants or used to generate entrant wealth. | Open | Decide whether to retain the reduced-form open bequest loop or introduce an inheritance-assignment rule that closes the intergenerational wealth distribution. |
| HF-1 | Uniform owner debt floor | Buyers and incumbent owners both face \(b'\ge-\phi PH_k\). Unlike the May write-up, an incumbent receives no exception that preserves pre-existing leverage after a price decline. | Accepted for current draft | Revisit whether forced deleveraging after house-price declines is intended before the final model is estimated. |
| HF-2 | Capital-gains lock-in absent quantitatively | The analytical model includes capital-gains taxation and stepped-up basis, but the quantitative model uses only sale costs, bequests, and a direct retention experiment. | Accepted for current draft; open extension | Add tax basis as a state and recompute retention and policy effects before presenting the capital-gains channel quantitatively. |
| HF-3 | Rental-size cap | The hard cap \(\bar h^R=6\) is motivated by the DUE rental-support rule, discussed extensively in the paper, and supported by the existing cap profile. | Accepted for current draft | Retain six rooms as the benchmark and 6.5 as the stated robustness case. |
| HF-4 | Housing-supply calibration | The active code fixes the supply level at \(\bar H=4\) and the elasticity at \(\eta=1\). The level is not internally calibrated, and the elasticity is below the population-weighted U.S. metro average reported by Saiz. | Open | Calibrate or invert \(\bar H\) to an aggregate housing-quantity or price anchor. Set \(\eta\) from external evidence and report sensitivity. |
| HT-1 | Within-period choice timing | The current paper and code nest fertility outside the contemporaneous housing and saving problem: an eligible household chooses family size, then housing, consumption, and saving conditional on that choice. | Open | Decide whether fertility should precede housing, follow housing, or be chosen jointly with housing; determine the corresponding taste-shock nesting and re-estimate if the model timing changes. |
| INC-1 | Lifecycle labor-income process | Labor income is a deterministic household-earnings profile \(e_a\) times a persistent mean-one component \(z\), discretized with a five-state Rouwenhorst chain. The architecture is standard, but the current five-point age profile and the empirical source for the AR(1) parameters are not documented. | Open; urgent | Choose and cite an external source for the household-earnings profile and earnings-risk parameters. Replace the current age-profile inputs if the chosen source implies different values, and document their normalization and four-year conversion. No solver or state-space change is required. |
| LC-1 | Lifecycle fertility profile | The original figure labeled a conditional family-size index as expected births. The paper figure now allocates reported completed fertility across ages using reconstructed family-size choice flows. | Corrected in paper figure; code diagnostic open | Repair the generic `fert_by_age` statistic for the Markov-income state dimensions and regenerate after the pending recalibration. |
| LC-2 | Sequential fertility and child cohorts | Fertility is a one-time completed-family-size choice. All children arrive in one cohort and share the same child-age stage, so the model cannot represent sequential births or children of different ages. | Open; urgent | Add parity and child-cohort ages to the state, allow birth choices over successive fertile periods, redefine the fertility moments consistently, and re-estimate before treating the lifecycle timing or policy results as final. |

## Review notes

### BQ-1: Gross terminal estate

The May write-up valued an owner's terminal estate at
\(b+(1-\psi)Ph\), as if the house were liquidated and the sale cost paid. The
current model instead allows the house to pass directly to the estate at market
value. This convention is acceptable for the current draft, but remains a
modeling choice to revisit.

### BQ-2: Family-size-dependent normalization

The May write-up subtracted \(\theta_1^{1-\sigma}\), setting bequest utility to
zero at zero estate wealth. The active code omits that term. This omission was
identified in the July 9 audit, but it was not included in the subsequent code
repairs. The paper was changed to describe the code instead. Because family
size multiplies the bequest term, the omitted constant is not behaviorally
irrelevant: it changes fertility values as well as saving incentives.

### BQ-3: Bequest-to-entry wealth link

The analytical and quantitative models take the distribution of entrant wealth
as given. Old households value the estates they leave, but those estates are
not assigned to particular descendants and do not generate the next cohort's
entry-wealth distribution. This open bequest loop is now stated explicitly in
the paper rather than described as anonymous inheritance. Closing the loop
would require an inheritance-assignment rule and a stationary distributional
condition linking estates to entrant wealth.

### HF-1: Incumbent deleveraging

The May write-up allowed an incumbent who kept the same house to satisfy
\(b'\ge\min\{b,-\phi PH_k\}\), avoiding forced deleveraging after a price
decline. The active code and current paper instead apply the origination-style
debt floor \(b'\ge-\phi PH_k\) to every owner. The uniform rule is accepted for
the current draft because it matches the estimated model and keeps the finance
block simple. It remains an outstanding modeling issue because price declines
can force an incumbent to reduce debt even without refinancing or moving.

### HF-2: Capital gains and stepped-up basis

The quantitative model does not track an owner's tax basis and therefore cannot
implement the capital-gains and stepped-up-basis mechanism emphasized in the
analytical section. The current draft states this omission explicitly. For the
present version, old-owner retention is generated through proportional sale
costs, the bequest motive, and a direct retention experiment. A later version
should determine whether adding tax basis materially changes the lifecycle
allocation of family-sized homes and the associated policy conclusions.

### HF-3: Rental-size availability

The current six-room cap is the maintained reduced-form availability limit for
large rental housing. The paper links it to the DUE rental-support rule and the
matched ACS rental distribution, and reports the existing profile over 5.5,
6.0, 6.5, and 7.0 rooms. Six rooms remains the benchmark and 6.5 the robustness
case. This issue is settled for the current draft.

### HF-4: Supply level and elasticity

The one-market code inherited \(\bar H=4\) and a unit-elastic supply curve.
These assumptions matter for how demand changes divide between prices and
quantities. The supply level should be calibrated or inverted to an aggregate
rooms-per-household or price anchor. The supply elasticity should be borrowed
from external evidence and varied transparently in sensitivity checks; Saiz's
population-weighted U.S. metropolitan average is 1.75, while later quantitative
work commonly rounds this benchmark to 2.

### HT-1: Within-period choice timing

The current nesting lets a household's fertility option value incorporate the
housing, consumption, and saving choices made in the same period. It therefore
allows a birth to change current housing immediately, but it treats housing as
chosen conditional on family size rather than as a prior commitment or a joint
choice. This timing is retained for the present exposition only. Before the
model is finalized, decide whether the economically relevant sequence is
fertility then housing, housing then fertility, or a joint fertility-housing
choice, and revise the Bellman nesting, choice shocks, measurement, and
calibration consistently if it changes.

### INC-1: Lifecycle labor-income process

The retained structure separates the common lifecycle path of household labor
earnings, \(e_a\), from persistent earnings differences among households of the
same age, \(z\). The current five-point schedule for \(e_a\) is not traceable to
an external estimate, and the paper reports the persistence and innovation
variance without naming their empirical source. Before the income process is
treated as final, select a literature or data anchor for both objects, update
the input values if necessary, and document how the annual process is mapped to
the model's four-year periods.

### LC-1: Fertility over the lifecycle

The original green series was not an aggregate birth flow: it used the expected
family-size index among selected childless states and did not apply the paper's
completed-fertility mapping. The replacement series reconstructs pre-choice
childless mass from the saved post-choice distribution and the probability of
choosing zero children, computes family-size choice flows by age, and allocates
the solution's reported completed fertility across those ages. The corrected
profile reveals a large entry-age spike as well as late family formation.

### LC-2: Sequential fertility

The current model chooses completed family size in a single period. Every
child implied by that choice enters the same cohort and ages through the same
dependent-child stage. The next fertility extension should instead permit
births in successive periods and track the parity and ages of children already
in the household. This is required to distinguish birth timing from completed
fertility internally rather than reconstructing timing from a one-shot choice.
