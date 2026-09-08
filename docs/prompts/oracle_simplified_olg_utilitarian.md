# Housing allocation and fertility: a utilitarian planner

We are developing a short illustrative economic theory note, alongside a richer quantitative model. Please work as an independent economic theorist. The task is mathematical: derive a small number of useful analytical propositions and check every feasibility and equilibrium step. The author wants something that reads like a clear problem set, not an elaborate research program. Keep the main answer compact (roughly 2,500 words if possible), with essential derivations in a short appendix. Take the reasoning time needed to resolve the difficult steps. Do not stop at a collection of first-order conditions or an outline when an actual proof is possible.

## The change in direction

We previously spent substantial effort proving Pareto improvements. The author now wants a SEPARATE note based on an equally weighted utilitarian objective. The Pareto material can remain supporting appendix material. The desired order is:

1. A direct-allocation planner: establish a general welfare statement and, crucially, when the planner moves housing from old to young households. The planner can bypass household financing and market implementation but respects real resources, the housing stock, physical tenure size limits, and intertemporal feasibility.
2. A planner restricted to taxes and transfers: households still optimize, housing markets clear, and the original mortgage and physical restrictions remain. Establish a useful constrained utilitarian welfare result with an actual market allocation, ideally with the same old-to-young housing direction.
3. Connect the intervention to fertility on the transition path. Separate dated or finite-path results that do not require convergence from conclusions about eventual population that require a convergent path and positive stationary endpoints.

The author is happy to specify a utilitarian criterion. We understand that failure to maximize one equally weighted welfare function is not itself proof of Pareto inefficiency. We want the economic result to reveal the financing and housing mechanism, not merely restate that a redistributive planner dislikes unequal wealth. Explain how to isolate the financing contribution. Do not continue optimizing the previous no-loser Pareto exercise.

## Author preferences and non-negotiable scope

- The attached LaTeX note is the exact candidate household model. Preserve its notation, positive child goods and space costs, heterogeneous income and wealth, ownership choice, two physical size limits, and mortgage timing unless a change is separately justified. It is a proposal, not an adopted author manuscript.
- Current income is available when buying a house. There is known old-age income. Old owners may resize within the owner size cap.
- Keep the conventional mortgage: principal and accumulated interest are repaid on entering old age. Do not silently replace it by an interest-only or amortizing covenant.
- A 20 percent down payment means financed share phi=0.8 and the bound is q a' + phi P h >= 0 in the attached timing convention.
- Welfare comparisons should initially hold individual fertility and cohort sizes fixed. Any subsequent fertility response comes from household choices, not a social objective that values adding people.
- Endowments of successive entering cohorts are exogenous; warm-glow estates do not fund their initial wealth. Do not silently add dynastic altruism or an estate distribution law.
- No numerical reference equilibrium, unspecified open neighborhood, numerical eigenvalue calculation or simulation as the central proof. Give explicit interpretable inequalities, preferably in primitives. High patience must be allowed. The author especially dislikes conditions expressed only using borrowing multipliers.
- Distinguish existence of some improvement, the direction of that improvement, and the allocation of the global planner optimum. A current housing reallocation is not the same statement as changing a household's lifetime housing profile.
- Do not infer a full transition theorem from a stationary comparison or a fixed-price derivative. We have not proved general transition convergence for this candidate model.
- A short, sound result with explicit limitations is preferable to an overstated broad claim. But do attempt to prove the broad result first and identify the exact obstruction if it fails.

## Welfare weights: make the choice explicit

We have not yet settled the age/cohort normalization. For a local dated comparison, one candidate gives equal weight to the remaining lifetime utility of every household currently alive. Current young utility is u^y + beta V_next; current old remaining utility is V. Another gives equal weights to lifetime utilities measured from each cohort's birth, which puts coefficient beta on the current old's remaining utility (past utility is fixed). For a full infinite-horizon planner, a well-defined criterion additionally needs summable cohort weights or another explicit welfare-comparison convention.

Please show which statements depend on this choice. Do not change weights merely to obtain the preferred sign. A locally finite intervention that leaves all later cohorts unchanged would avoid some global welfare aggregation issues, if it can be implemented through the permitted transfers.

## A potentially much simpler direct-allocation argument: verify independently

Write p_t=(1+q tau^p)P_t-qP_(t+1)>0 and L_t=(1-phi_t+q tau^p)P_t>0. Let m_next=1/c_next^2 be marginal old utility of resources, Lambda=1/(c-chi n) the young marginal utility of current resources, and mu>=0 the multiplier on q a'+phi P h>=0. At slack young housing caps, the original first-order conditions appear to imply

Lambda = (beta/q)m_next + mu,

alpha/(h-kappa n) = (beta/q)p_t m_next + mu L_t.

For a current old owner with slack housing cap and estate floor, gamma/h_old=p_t m_current. At a stationary equilibrium, a current young type and its predecessor in the old cohort have m_next=m_current=m. Thus the raw marginal welfare gain from transferring housing from old to young, using equal remaining-utility weights, would be

[(beta/q)-1]p m + mu L.

This is positive if beta>=q and finance is strict. At beta=q the term remaining is mu L, which might isolate the financing distortion. For equal birth-normalized lifetime weights, subtract beta gamma/h_old instead: the gap becomes beta(1/q-1)p m+mu L, positive even without a borrowing wedge. This illustrates why the normalization matters.

Check the algebra, exact planner feasibility, ownership and income matching under heterogeneous F, and all missing conditions. The attached note's appendix already provides fully primitive sufficient conditions for stationary existence and positive mass of owners with strict finance, slack young/old caps and slack estate floor. Can those be reused without adding a new numerical construction? Do not treat exact type matching as a positive mass atom when F is continuous.

The direct planner's housing-only first-order condition under equal current remaining-utility weights is alpha/(h_y-kappa n)=gamma/h_o at slack physical caps. Can the equilibrium wedge above deliver a short proposition about the direction of the planner's housing reallocation, rather than merely restating that an observed marginal utility difference would justify a transfer? What is required to strengthen a local directional improvement to an aggregate comparison with the global utilitarian housing allocation?

## The difficult priority: transfers followed by genuine market choices

At fixed prices, moving current cash from old j to young i gives marginal equal-remaining-utility welfare gain 1/(c_i-chi n_i)-1/c_j^2. That is only the direct effect. A valid constrained result needs all induced price, rebate, inherited-title, and future effects.

In particular, do not mistake making a proposed bundle affordable for implementing it. The previous direct planner compensates an old donor for surrendering housing. At unchanged prices, a positive lump-sum compensation would generally make that old household demand MORE housing, not voluntarily donate it. Thus its budget identities are not a transfers-only proof.

Please formulate a natural transfer authority precisely. Distinguish a balanced one-time old-to-young tax/transfer from richer committed age- and type-specific transfers with later taxes and government borrowing repaid. Transfers may fund the required down payment but the household mortgage rule itself stays unchanged. Do not call unrestricted government lending a genuinely restrictive benchmark without explaining what remains constrained.

Seek an analytical theorem for an actual transfer equilibrium. Possible routes include a small intervention with exact supporting prices and finite future settlement, or a local market-clearing response with explicit, economically meaningful sufficient conditions. Fixed fertility and tenure can be the initial welfare benchmark, followed by a separate tenure/fertility extension. If the full claim cannot be proved generically, show the obstacle and derive the narrowest useful sufficient condition. Avoid replacing the model by a one-period partial-equilibrium problem without clearly identifying that change.

## Fertility and the path

The household fertility condition is

theta/n = chi/(c-chi n)+alpha kappa/(h-kappa n).

One useful exact differential is

dn = [(chi/x^2) dc + (alpha kappa/s^2) dh] /
     [theta/n^2+chi^2/x^2+alpha kappa^2/s^2],

where x=c-chi n and s=h-kappa n. For a finite comparison, evaluating the new consumption/housing bundle at the original n may yield a clean necessary-and-sufficient sign test. We already have mortgage-phi comparative statics in the attached note. A current cash transfer need not have the same sign conditions as a mortgage relaxation; please derive the response to the actual intervention in your constrained theorem. Also account for endogenous tenure if claiming cohort-average fertility.

The intended demographic experiment starts at a steady state, suffers an exogenous persistent decline in the fertility preference theta, and introduces housing policy later during the ensuing transition. We compare the policy path with continuation without policy, from the same inherited state. With Y_(t+1)=nu nbar_t Y_t and O_(t+1)=Y_t, finite population differences follow products of fertility ratios. At positive stationary endpoints both fertility levels equal 1/nu, while total population can differ. The housing identity implies final population is larger only if equilibrium lifetime housing per cohort is smaller. None of these accounting identities by itself proves the desired policy signs.

Please label separately:

- A dated household or equilibrium fertility result valid at any admissible state, independent of eventual convergence.
- A finite population-path comparison, with all future price assumptions stated.
- A long-run population result conditional on existence and convergence, or a genuine convergence theorem if you can supply one analytically.

## Deliverable

Lead with a short verdict and the strongest three statements that the new note can responsibly make, in the requested order. State the exact assumptions and prove the nontrivial steps. Identify precisely where an additional model or planner choice is needed. Then recommend the simplest coherent main-text version and what should remain appendix material. We are trying to close a short theory note, not create a large new theory project.
