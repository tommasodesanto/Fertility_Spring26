# Fertility nests with conception-contingent housing

Author decision, September 7: try simultaneous fertility groups, with the
retained sequential model as fallback. First check the complete menu, then
numerics and the retained-parameter objective; no recalibration is authorized.
The monitoring routine was deleted at the author's request and must not be
recreated. This note is the first specification checkpoint, not model adoption.

## Finding

A simultaneous GEV choice over complete housing plans can represent the current
sequential probabilities exactly. It requires additional outcome-dependent
housing subnests beyond the simple four-alternative illustration. Those extra
correlations have NOT been adopted. They preserve housing-product dispersion
and observable probabilities, not the literal reuse of original housing shocks.
There is no reason to run another calibration merely to verify an algebraically
equivalent representation. A full implementation or empirical result is not
claimed by the arithmetic test below.

## Maintained menu

The retained national profile has one location and six housing products:
renter plus owner sizes 2, 4, 6, 8 and 10. Let q_bh be the optimized value of
product h after biological outcome b=0 (no birth) or b=1 (birth), including
transactions, feasibility, consumption, saving and continuation. All products
remain separate alternatives. Let p be conception probability, k the housing
logit scale, f the first- or later-birth fertility scale, and C the first-birth
fixed cost (zero for later births). Define

    L_b = k log sum_h exp(q_bh/k).

The original solver uses

    I_wait = L_0
    I_attempt = (1-p)L_0 + p L_1 - p C
    V = f log[exp(I_wait/f) + exp(I_attempt/f)].

Complete alternatives are wait(h0) and attempt(h0,h1). The latter has
deterministic value (1-p)q_0h0 + p q_1h1 - p C. It permits renting after
failure and owning after success, or the reverse. There is no tenure commitment.

## Equivalent simultaneous distribution

Use fertility groups at the top, with scale f. Within wait, housing scale is k.
Within attempt, denote the more likely biological outcome by M and the less
likely by m, with probabilities w_M and w_m. Group complete plans by their
housing choice h_M, with scale w_M k; within those groups, choose h_m with
scale w_m k. The valid GEV scale order is

    f >= k > 0, and w_M k >= w_m k > 0.

All plan-level GEV shocks are observed before choosing a complete plan.
The biological outcome remains unresolved and independent. Center each leaf
shock once, using the common marginal root scale f, so that the expected
maximum has no Euler-constant term.

Expanding the logsums gives exactly the I_attempt above. Conditional attempt
plan probabilities factor as logit(q_0/k) times logit(q_1/k). Hence fertility
probabilities and realized product/family probabilities match production state
by state. With the same continuation and transition maps, backward induction
gives the same lifecycle solution. Selected event-study cohorts must also retain
the original origin weights and destination continuation policies; aggregate
fertility or ownership equality alone is insufficient evidence.

At p=0 or p=1, collapse only the zero-probability housing coordinate. Do not
remove an attempt label that production retains. Preserve all existing age,
parity and readiness restrictions. The construction is presently restricted
to the one-location national model.

## Economic qualifications

The additional subnest scales depend on conception probabilities and their
orientation switches with the more likely outcome. This is a special
correlation restriction chosen to preserve the existing recursion. It is not
implied by simultaneous taste revelation or fertility grouping alone.

Leaf errors are attached to entire contingent plans. They do not equal a
fertility shock plus probability-weighted, reused outcome-specific housing
shocks. In particular, additive reused shocks have zero rectangular differences
across the Cartesian housing-plan menu; general positive-scale GEV leaf shocks
do not. Thus this is an observationally equivalent simultaneous interpretation,
not evidence that revealing the original independent shocks sooner is innocuous.

A simple single-scale housing nest over all contingent attempt plans instead
has conditional scales k/(1-p) and k/p. It changes housing dispersion and creates
a menu-size bonus from the nearly irrelevant outcome as its probability tends
to zero. Neither that version nor the deeper-tree representation should be
introduced as a silent implementation detail.

## Verification and next decision

The independent Astra/max review confirmed the derivation, scale ordering,
centering and endpoint treatment. Lead implemented an independent full-plan
enumeration check: 1,320 synthetic menus, six products, masks, conception
probabilities from zero to one, both retained fertility scales, flat-logit
boundaries, and deterministic utility spreads through 1,000. The largest
value error is 1.78e-15; the largest realized product/family probability error
is 1.23e-15. This is pure arithmetic, not a Bellman, equilibrium or objective run.

Reproduce from this isolated branch with an installed NumPy environment:

    python code/model/tools/check_fertility_nest_menu.py --output output/model/fertility_nest_menu_check/summary.json

The baseline scientific source remains revision ac676c2a, which was used by the
retained diagnostic snapshot. No production source was edited. Existing saved
choice captures only contain the renter and maximum-owner values; they cannot
verify the complete six-product menu. No new lifecycle run was made.

Author decision needed before adoption: accept this equivalent simultaneous
contingent-plan interpretation, including its explicit extra correlations, or
choose a different correlation structure and assess a genuinely changed model.
Sequential remains the retained fallback in either case.

Reference: Kenneth Train, Discrete Choice Methods with Simulation, chapter 4,
especially the multilevel GEV scale restrictions:
https://eml.berkeley.edu/books/choice2nd/Ch04_p76-96.pdf
