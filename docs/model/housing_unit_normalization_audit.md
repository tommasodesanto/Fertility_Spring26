# Housing-unit normalization: equation audit

October 1, 2026. User requested the unit-equivalence audit and a transparent advisor note, subsequently expanded to explain the dollar conversion and reference rents. No model code, economic specification, calibration input, manuscript, or deck changed. The separately requested note is [editable LaTeX](../../latex/housing_unit_normalization_note.tex); its [PDF](../../output/model/fixed_reference_economics_20260928/room_unit_audit_v1/housing_unit_normalization_note.pdf) and arithmetic evidence belong to the existing fixed-reference economics output tree.

## Finding and limits

The implemented economic equations admit a coherent change of housing coordinates. At fixed monetary units, let \(H'=\lambda H\). Converting all room quantities by \(\lambda\), all per-room prices and reference rents by \(1/\lambda\), and all utility-valued coefficients by the common material-utility factor preserves budgets, feasible physical allocations, choice-value ratios, demand/supply balance and population accounting. This is an equation-level result, not a certificate of an executed rescaled equilibrium.

The standalone arithmetic record has **1,260 comparisons**, maximum relative discrepancy \(1.911\times10^{-14}\); all are below \(10^{-12}\). These are scalar evaluations of transcribed equations, not model tests, simulations, recomputed policies, or solved equilibria. Six inspected engine source files match the free-psi floor experiment's source pins. Results: [arithmetic.json](../../output/model/fixed_reference_economics_20260928/room_unit_audit_v1/arithmetic.json). The runnable arithmetic recipe is retained in that packet solely to reproduce these calculations; it imports no model modules.

This establishes no independent empirical validation of the linear room-to-services mapping, tiny continuous rental quantities, house values/rents relative to earnings, Cobb--Douglas substitution, or the provisional supply elasticity. It also does not prove that past changes of grid, reference rent, utility coefficients or parameter bounds were consistently converted.

## Specification identities

1. **Adopted 2007 stationary reference, block0506, September 28 verified export.** Actual parameters are in [fixed_reference_manifest.json](../../output/model/fertility_identification_20260928/fixed_reference_manifest.json). The reference uses \(H_{own}=(2,4,6,8,10)\), rental cap 6, baseline \(\alpha_0=0.733\), \(\sigma=2\), compensated child-dependent shares and no active housing floor. The raw serialized legacy `h_bar_0=4` is inactive: `shared.py` zeros housing floors in this branch. Reference price and aggregates come from the [saved diagnostic summary](../../output/model/fertility_identification_20260928/resume_v1/selected_export/primary/standard_diagnostics/summary.json). This is a normalized-population housing-clearing specification. Complete [parameter/bound table](../../output/model/fertility_identification_20260928/resume_v1/selected_export/primary/parameters.csv) and [target-fit table](../../output/model/fertility_identification_20260928/resume_v1/selected_export/primary/target_fit.csv) retain the full calibration record; no recalibration or fit comparison was performed here.
2. **Separate experimental October 1 09:41 NY floor point, chain 7 case0173_nm.** The [31-row parameter table](../../output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/deployment/monitor_snapshot/verified_global_20261001T0941NY_chain7_0173/ROOT/parameters.csv) and [14-row fit table](../../output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/deployment/monitor_snapshot/verified_global_20261001T0941NY_chain7_0173/ROOT/target_fit.csv) identify this point. Physical parent floor is 2.3 rooms, alpha is common, compensation is off, unsecured credit is zero and entry is nonnegative. It is not adopted and is not a D=.25 pilot or D=.53 diagnostic. The [closure](../../output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/deployment/monitor_snapshot/verified_global_20261001T0941NY_chain7_0173/ROOT/closure.json) sets price through birth renewal and population through absolute housing clearing. H0 is fixed and not identified by per-household targets under this closure.

## Exact conversion recipe

For ten-room units, \(\lambda=0.1\); for tenths-of-a-room units, \(\lambda=10\). Monetary quantities and time periods remain unchanged. The recipe is a statement of an equivalent specification, not changes applied to the project.

| Object | Transformation | Economics |
|---|---|---|
| Owner grid, renter policy, cap, physical floor, room thresholds, H0 | multiply by lambda | Same physical space and capacity |
| Asset price p, rent r, supply reference rent r_bar, active utility_reference_rent | divide by lambda | Same cost per physical dwelling |
| Consumption, income, financial wealth, wealth grid, theta1 | unchanged | Same monetary unit and feasible saving choices |
| Alpha, sigma, beta, financing share, rates, child equivalence scale, owner premium chi, fertility curvature | unchanged | Dimensionless preferences, technology and institutions |
| Child-benefit amplitude, fixed birth utility cost, theta0, fertility/tenure/location utility shock scales, additive location amenities/moving disutilities | multiply by K | Same comparisons in a common utility unit |
| Derived child-benefit CRRA coefficient | multiply by K | It is linear in the child-benefit amplitude |
| Room-level moment targets and uncertainty | multiply by lambda | Same empirical measurements in new coordinates |
| Weights on squared room-level gaps | divide by lambda squared | Same objective contribution |
| Search bounds and finite-difference scales | transform with their parameter | Same economic search region |

Here \(K=\lambda^{(1-\alpha_0)(1-\sigma)}\). At \(\lambda=0.1\), \(K=1.849\); at \(\lambda=10\), \(K=0.541\). `chi` is an owner service premium, not the child-benefit coefficient used in the external Fable response. Do not rescale `chi` as an additive utility weight. `theta1` is a monetary bequest shifter and stays fixed; `theta0` is utility-valued and converts.

The reference has one location. Location probabilities are degenerate, while its move disutility is inactive for the realized one-location choices. The general recipe still converts all additive location disutilities and utility shock scales. Inactive legacy fields are not evidence of an active housing requirement.

Optional nonlinear size cost is zero in the inspected reference. If activated, a cost \(\eta p(H-H_{ref})^\nu\) additionally needs \(\eta'=\eta\lambda^{1-\nu}\) and \(H_{ref}'=\lambda H_{ref}\). An active quadratic renter-cost wedge would likewise require coefficient-specific dimensional conversion. These inactive extensions are outside the arithmetic certification.

## Why the household problem is equivalent

### Budgets and financing

Owner asset value is `p_hat * H_own`; sale proceeds are `(1-selling_cost)*p_hat*H_own`; down payments are `(1-financed_share)*asset_value`; running costs are `(depreciation+property_tax)*asset_value`. Hence \(p'H'=pH\). Rent payments satisfy \(r'H'=rH\). Income-aware lending limits and wealth interpolation are functions of monetary quantities, so their economic feasible sets stay unchanged. Renter floors/caps convert together with their quantities; the floor's monetary cost \(r h_P\) stays unchanged.

Sources: [household.py:158](../../code/model/refactor_lab/engine/household.py#L158), owner scalar objective [kernels.py:173](../../code/model/refactor_lab/engine/kernels.py#L173), renter allocation [kernels.py:647](../../code/model/refactor_lab/engine/kernels.py#L647), payment-to-income restriction [household.py:1520](../../code/model/refactor_lab/engine/household.py#L1520). Mortgage/solvency objects are monetary; the audit does not adopt or change the credit contract.

### Material utility and reference compensation

With a common weight, \(u^{mat}=e(m)^{\sigma-1}[c^\alpha s^{1-\alpha}]^{1-\sigma}/(1-\sigma)\). Effective owner services are \(s=\chi(H-h_P)\), renter services are \(s=H-h_P\), and the floor is zero for childless households. Converting both H and h_P gives \(s'=\lambda s\), so \(u^{mat\prime}=K u^{mat}\). Physical parent-owner infeasibility at \(H\le h_P\) is preserved.

The adopted share reference requires an additional step. Let \(\alpha_m\) be the family-specific consumption weight, \(k(a,r)=a^a[(1-a)/r]^{1-a}\), and \(A_m=k(\alpha_0,r_*)/k(\alpha_m,r_*)\). The code uses \(A_m\) to compensate the material composite. Under \(r_*'=r_*/\lambda\),

\[
A_m'=\lambda^{\alpha_m-\alpha_0}A_m,
\qquad A_m'c^{\alpha_m}(\lambda s)^{1-\alpha_m}
=\lambda^{1-\alpha_0}A_mc^{\alpha_m}s^{1-\alpha_m}.
\]

The same K therefore applies to every family state. Holding the numerical reference rent fixed would fail this common-factor property. Source: [child_preferences.py:39](../../code/model/refactor_lab/engine/child_preferences.py#L39), called after family objects are built in [shared.py:240](../../code/model/refactor_lab/engine/shared.py#L240).

### Fertility, bequests and dynamic choices

Child flow benefit is \(\psi_{child}m^{1-\gamma}\). Multiply psi_child by K, leaving gamma unchanged. Bequest utility is linear in theta0, with its monetary arguments and theta1 unchanged, and therefore also multiplies by K. The first-birth utility cost, additive amenities and moving disutilities convert by K. With unchanged beta and transition probabilities, backward induction yields \(V'=KV\) for each corresponding feasible allocation.

The inclusive value obeys \((K\kappa)\log\sum_j\exp(KV_j/(K\kappa))=K\kappa\log\sum_j\exp(V_j/\kappa)\), and the choice probabilities are identical. Thus the birth-attempt and tenure odds, their continuation values and optimal saving/consumption are invariant at the equation level. Arithmetic checks of softmax ratios use explicitly illustrative value contrasts, not saved Bellman values or freshly computed policy probabilities.

Sources: child benefit [child_preferences.py:30](../../code/model/refactor_lab/engine/child_preferences.py#L30); first birth [household.py:975](../../code/model/refactor_lab/engine/household.py#L975); later-birth scales [household.py:1007](../../code/model/refactor_lab/engine/household.py#L1007); tenure [household.py:590](../../code/model/refactor_lab/engine/household.py#L590); bequests [household.py:1480](../../code/model/refactor_lab/engine/household.py#L1480). The one-location case avoids any additional active spatial-choice economics.

## Why supply curvature and population are equivalent

Source: [distribution.py:3055](../../code/model/refactor_lab/engine/distribution.py#L3055), inverse supply [equilibrium.py:320](../../code/model/refactor_lab/engine/equilibrium.py#L320). Rent is user-cost rate times asset price. User-cost rate does not convert.

\[
S'=\lambda H_0\left(\frac{r/\lambda}{\bar r/\lambda}\right)^\xi=\lambda S,
\qquad r'(Q')=\frac{r(Q)}{\lambda}.
\]

The raw inverse-supply derivative satisfies \(dr'/dQ'=(dr/dQ)/\lambda^2\), while \(d\log r'/d\log Q'=1/\xi\). At xi=.63, a 5% quantity expansion along the curve implies an 8.052% price increase. This is not the GE effect of a 5% demand shift. If reference rent is held numerically fixed instead of converted, the equivalent supply intercept would be \(H_0'=\lambda^{1+\xi}H_0\), not simply lambda H0; converting reference rent explicitly is clearer.

Fixed-population housing demand and supply both convert by lambda. In the later renewal-price experiment, \(N=S/\bar H\), so \(N'=\lambda S/(\lambda\bar H)=N\). Per-household fertility and entry probabilities stay unchanged and hence the renewal root is unchanged in physical-price units. Source: [phase_b_pilot.py:93](../../output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/phase_b_pilot.py#L93), especially lines 118--121. Estimating H0 to match a quantity level is a separate economic calibration choice; unit conversion does not endogenize H0. Under the newer population closure, changing H0 alone rescales N and cannot independently repair per-household demand or birth incentives.

## Benchmark population one and the supply coefficient

This is distinct from converting rooms into ten-room housing units. Keep rooms per household and the earnings-based monetary unit fixed. At benchmark rent \(r_b>0\) and household housing demand \(\bar H_b>0\), normalizing benchmark household population to one makes housing clearing

\[
\bar H_b=H_0(r_b/\bar r)^\xi,
\qquad H_0=\bar H_b(\bar r/r_b)^\xi.
\]

Thus a positive coefficient consistent with that benchmark price and demand always exists. This guarantees the aggregate housing identity, not fit to empirical targets or the separate birth-renewal and fiscal conditions. Searching over prices and preferences and deriving this coefficient need not be the same numerical procedure as searching over H0 and solving for prices; matching their admissible sets, roots, bounds and demographic definitions remains a separate question.

For an existing solution with population \(N^*\), a pure change of population units divides population and all aggregate housing quantities by \(N^*\), giving \(H_0^{(1)}=H_0/N^*\). The retained floor illustration has \(N^*=0.9222615666500361\) and \(H_0=6.293507689200028\), so \(H_0^{(1)}=6.82399431655833\). At unchanged period rent \(0.12965123763644737\), this yields supply \(5.977125401038295\), equal to its housing demand per household. These are arithmetic on the [saved closure](../../output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/deployment/monitor_snapshot/verified_global_20261001T0941NY_chain7_0173/ROOT/closure.json), not a rerun, new calibration estimate or adoption of that coefficient. A pure change of population units preserves this solution's price and household choices; it does not repair its price-level gap. An actual recalibration can change choices and prices.

The complete supply curve and all initial aggregate stocks and flows must use the same population units, and those units must be retained during transitions. An externally fixed aggregate that is left in its old units would make this an economic change rather than a consistent reexpression of the same solution. Removing the supply reference rent separately uses \(A=H_0/\bar r^\xi\) and preserves the curve algebraically; this does not alter the population normalization or the meaning of a room.

## Arithmetic coverage and remaining checks

The retained standalone recipe evaluates lambda=.1 and 10 for both identified specifications, consumption levels .1/1/10, children-at-home states 0--3, owner/renter services, physical sizes .35/2/2.3/2.4/4/6/10, all five purchase nodes, financial costs, reference compensation, floor feasibility, child flow utilities, bequests, illustrative choice-value ratios, supply, inverse supply, supply slopes, population accounting and squared room-moment weights. It uses saved JSON/CSV inputs and Python standard-library math; it never imports model packages.

**Established:** dimensional conversion and common-factor equivalence of the inspected economic equations; arithmetic identities; source identity for the six reviewed files.

**Not established:** full solved policy/distribution equivalence, floating-point or grid convergence, absence of historical conversion errors, empirical validity of price levels, services or elasticities. Fixed utility sentinels (`-1e10`, DEAD_VALUE_CUTOFF `-1e9`), housing regularizers (`1e-10`), supply demand floor, very-small taste-scale safeguards and absolute price tolerances/bounds do not automatically convert. Their units and activation would need attention in an executed numerical comparison. No source repair or new model run is implied by this note.

The physical meaning of ROOMS comes from the existing [2007 AHS builder](</Users/tommasodesanto/.codex/worktrees/54a9/Fertility_Spring26/code/data/ahs_supply_snapshot/build_ahs_2007_room_target.py>) and its authoritative Census definitions. The [accepted input reconciliation](accepted_input_reconciliation_20260926.md) records that a joint empirical price--quantity anchor was deferred. That outstanding validation is separate from the factor-of-ten mathematical question.
