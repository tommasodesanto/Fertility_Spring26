# One maturation-to-entry rule: complete accounting and implementation assessment

**September 24, 2026. Derivation and exact toy checks; no model solve or specification/source edit.** The primary construction retains the author's \(1/2.1\) entry transformation, ordinary independent departure probability \(\mu=2/9\), the couple decision-maker, existing entrant wealth and pooled inheritance rules, and the current household decision states. It replaces birth-vintage entry with actual ordinary departures plus release of remaining dependents when the parental household exits. The tagged top-bin ledger below is a **proposed implementation detail**, not an adopted specification.

The author's objection is correct: **applying the same factor to both departure routes is consistent.** Twelve departing units, whether ten ordinary plus two released at death or any other disjoint split, produce \(12/2.1=5.714\) entrant household units. There is no reason to replace that factor with \(1/2\) merely to implement death release. The earlier objection conflated a replacement-normalized household model with a literal census identity. The former is coherent; the latter is a different additional requirement.

## 1. Objects, timing and one-period operator

A household is the existing couple decision unit. The distribution \(g_t(x)\) counts its mass immediately before fertility choices. Its state \(x\) includes adult model age, wealth, housing/tenure, location, earnings, children ever born \(n\), and represented children currently dependent \(m\). Children can remain economically dependent beyond chronological adulthood. The ordinary departure process uses the existing binomial transition and includes the newborn in the current period's draw. [S1–S3]

Initially suppress the top-bin adjustment and define \(H_t=\sum_xg_t(x)\), dependent stock \(C_t=\sum_xm(x)g_t(x)\), and births \(B_t\). Section 3 gives the exact weighted extension. Calendar stocks are dated consistently: births and household choices occur at \(t\); departures and exits during \(t\to t+1\) supply entrants present at the beginning of \(t+1\). New entrants make no choices in period \(t\).

1. **Fertility.** Apply existing fertility choices and realized conception probabilities. This preserves household mass and produces post-fertility \(\bar g_t\), with dependent stock \(\bar C_t=C_t+B_t\). Carry the weighted top-bin moment through this same stage if used.
2. **Household exits.** For each adult age \(j\), let \(s_j\) be the existing survival probability, setting \(s_J=0\) at the final age regardless of the optional earlier-age survival switch. The exiting household mass is \(X_t=\sum_x(1-s_j)\bar g_t(x)\). Release **all** attached dependents of these households: \(R_t=\sum_x(1-s_j)m(x)\bar g_t(x)\). None of this mass enters the ordinary-departure calculation. This is forward-mass ordering: current consumption, housing, saving and estate policies are still evaluated for all households, including those that exit; earlier splitting is valid here because survival depends on adult age alone.
3. **Surviving incumbents.** Advance \(s_j\bar g_t\) through the existing location, tenure, saving and earnings transports. These preserve household and dependent counts before child departure. Apply the existing child kernel
   \[
   K(m'\mid m)=\binom m{m'}(1-\mu)^{m'}\mu^{m-m'}.
   \]
   Its loss of dependent units is \(M_t\), the ordinary departure flow. With constant \(\mu\), \(M_t=\mu\sum_xs_jm(x)\bar g_t(x)\); a before/after stock difference is a useful independent check. The survivor distribution contains \(H_t-X_t\) households and \(C_t+B_t-M_t-R_t\) dependents.
4. **Single entrant flow.** Put \(F_t=M_t+R_t\), choose the retained conversion \(c=1/f\), \(f=2.1\), and insert
   \[
   E_{t+1}=c\rho_t F_t+I_t
   \]
   households **once** at age index zero, \(n=m=0\), with the existing entrant tenure/wealth/earnings distribution. Here \(\rho_t\) is the existing regional retention fraction if that closure uses one, and \(I_t\) is its explicit outside household inflow. Closed operation uses \(\rho=1,I=0\). Both departure routes use the same conversion and entrant distribution; no genealogical matching is required. Current dated code is one-market-only; preserve that scope. [S3–S5]
5. **Optional observed-age bridge.** Only where the existing historical contract requires empirical reweighting, apply it after entrant insertion. Measure its net household and dependent additions as \(J_t^H,J_t^C\), including any top-bin moment. A bridge is an externally imposed distribution adjustment; record its wealth, household, child and fiscal effects separately. It cannot be silently carried into a closed forecast.
6. **Estates and other accounts.** Record estates at the existing death event and pay them through the existing receiver rule, if active. New entry does not itself authorize another payment of those estates. Record exogenous entrant wealth separately. See section 6.

This gives the same stock-flow law in every period, including initial and terminal dates:
\[
\boxed{C_{t+1}=C_t+B_t-M_t-R_t+J_t^C,\quad
H_{t+1}=H_t-X_t+c\rho_t(M_t+R_t)+I_t+J_t^H.}
\]
During a dated transition, retain distribution **levels**. If the implementation stores a normalized density, carry its changing household-mass scalar and use it in every flow and market aggregate; normalizing each date back to one without that scalar erases demographic change. Roundoff correction is not an economic normalization.

No new orphan-care, marriage, parent-age-cap, child-age, child-death or individual-spouse state is needed. Dependents survive the parental exit by construction; this is an explicit assumption of the proposed rule.

## 2. What \(1/2.1\) means, and stationary replacement

Define the **replacement-unit index** \(Q_t=fH_t+C_t\). It is an index assigning \(2.1\) reproductive units to one household; it is not a census of actual adults and children. The above operator implies exactly
\[
\boxed{Q_{t+1}-Q_t=B_t-fX_t-(1-\rho_t)F_t+fI_t+fJ_t^H+J_t^C.}
\]
All retained maturation cancels, irrespective of route. No omitted flow is needed to make this *normalized* accounting internally consistent. Household choices may remain those of a couple; the demographic index's weight does not claim that the household contains 2.1 adults. No new mortality, migration or residual person class is being assumed by this recommendation.

If, instead, a result is advertised as **literal resident people**, with individually counted children and exactly two adults per counted household, then \(P_t=2H_t+C_t\) obeys
\[
\Delta P_t=B_t-2X_t-(1-\rho_t)F_t+(2c-1)\rho_tF_t
 +2I_t+2J_t^H+J_t^C.
\]
For \(c=1/2.1\), \(2c-1=-1/21\). This does not invalidate the normalized model. It identifies the additional assumption needed for a literal census interpretation: either (i) use \(c=1/2\); (ii) explicitly send \(F/21\) people outside that counted population; or (iii) track them in a separate resident class with its own subsequent flows and resource demands. Calling (ii) mortality or emigration would be a **new economic assumption**, not a deduction. Another exact weighting is \(2H+(20/21)C\), a rescaling of \(Q\); it likewise ceases to count every represented child literally. Merely rescaling all stocks by the same constant cannot turn \(1/2.1\) into literal pairing.

Let \(b\) be expected lifetime births in the selected demographic units **per entering household**, including any mortality before fertility. Since every dependent ultimately exits either normally or at terminal parental exit, steady-state \(F=B=bE\). Finite household lifetimes imply \(X=E\). The stationary renewal condition is therefore
\[
E=\rho c bE+I.
\]
In the closed normalized reference, \(c=1/2.1\) and \(b=2.1\) give exact replacement. Under the documented support—fertility ends at 42, mortality starts at 66, and post-fertility survival is independent of fertility state—the weighted completed-fertility observer at ages 46+ can equal this lifetime \(b\). Verify that equality in the selected parameter object and saved observer. If mortality moves into reproductive ages, survivor-conditioned completed fertility alone no longer establishes it. [S6,S7]

With literal \(c=1/2\), the same \(b=2.1\) gives 1.05 descendant households per entrant. Constant size would require, for example, retention \(20/21\); that is another assumption, not an accounting necessity under the retained normalization. The conversion and the reference fertility remain separately named quantities even when their values are chosen to match.

**Two boundary consequences are real but predate this change.** In a closed population, a terminal positive constant-size stationary state requires \(cb=1\). If permanently lower fertility gives \(cb<1\), an arbitrary positive normalized terminal household distribution is not a closed stationary endpoint. An open inflow can support it, or one must solve the demographic decline and its compatible economic tail; a growing/declining normalized composition is not itself a constant-level equilibrium. At exact replacement with no outside inflow, demography does not determine the population level, so retain the model's separately specified scale normalization. For the existing open stationary closure, replace its maturation counter with the unified one in \(S=I/(e_0-\rho c f_0)\), where \(e_0\) is entry and \(f_0\) is departing child units per unit scale; require a positive denominator if \(I>0\). [S2]

## 3. Retain the 3+ birth weights without a larger household problem

The current decision state stops at three children; aggregate renewal weights the top group as \(w_3=3.602\) (full precision in the receipt). Entry into that group currently contributes \(1+\delta\) birth units, where \(\delta=w_3-3\). A blanket multiplier on every dependent in a top-bin household is wrong: it changes the weights of earlier children when the third is born, including children already gone. [S1,S3]

A minimal exact **first-moment ledger** solves this under the existing symmetric departure law. Label the represented third child at its birth, solely for accounting. Let \(z_t(x)\) equal \(\delta\) times the household mass in cell \(x\) whose third child is still dependent. This is a forward distribution statistic, **not a new state on which households optimize**.

- On a transition into \(n=3\), add \(\delta\) times that birth-branch household mass to \(z\), at the resulting \((n,m)\) and existing asset/income/location coordinates.
- Carry \(z\) through all policy transports using exactly the probabilities and interpolation weights applied to its associated \(g\). Choices depend on the existing state, so this propagation is linear.
- Parent death releases \((1-s_j)z\).
- Conditional on \(m\to m'\), an attached tagged child survives with probability \(m'/m\), because all \(m\) children face the same hazard. Thus the exact moment transition is
  \[
  g'_{m'}\mathrel{+}=sK(m'\mid m)g_m,\qquad
  z'_{m'}\mathrel{+}=sK(m'\mid m)\frac{m'}m z_m\quad(m>0).
  \]
  The extra ordinary departure flow on that branch is \(sK(m'\mid m)(1-m'/m)z_m\); set \(z=0\) at \(m=0\).

Use \(C^w=\sum_x[m(x)g(x)+z(x)]\), \(B^w=B^{\rm explicit}+\delta B^{\rm third}\), and the corresponding weighted \(M^w,R^w\) in **all** equations above. Keep the parallel literal-state ledger \(C^0=\sum mg\) for household cost/housing diagnostics. These are different units and must be labeled. The existing top-bin resource approximation remains: decisions use the represented \(m\); the added \(\delta\) demographic weight shares the represented third child's timing. This preserves the source's top-bin birth timing convention and supplies its missing departure ledger; it is not a newly estimated family-size/resource model.

This moment closes exactly because maturation is symmetric and household choices do not depend on the tag. It would need reconsideration under an asymmetric newborn exemption or identity-dependent departure law. An unconditioned scalar extra-child stock would not be enough when state-dependent policies sort households: use the matching forward cells. Dense storage of the \(n=3\) slice costs at most one quarter of a full four-\(n\)-state distribution array; no Bellman array grows. Bounds are \(0\le z(x)\le\delta g(x)\), with \(z=0\) when \(n<3\) or \(m=0\).

The smaller alternative is to use literal births and departures throughout and leave 3.602 only in reported completed fertility. That is coherent, but then the reported target of 2.1 need not imply literal lifetime births of 2.1 or replacement at \(1/2.1\). **Recommendation: keep the existing adjusted renewal units and implement the tagged moment.** This is the only genuinely new small representation choice needed for the primary construction; it avoids silently changing the fertility normalization's demographic meaning.

## 4. Tiny transition checks

All results use \(c=10/21\), no retention loss or outside inflow unless stated, and fractional continuum household masses.

| Situation | Event accounting | Next state and interpretation |
|---|---|---|
| No deaths; two households each have a newborn | \(H=2,C=0,B=2; M=4/9,R=X=0\) | \(C'=14/9\), \(H'=418/189=2.212\). \(Q\) rises exactly by the two birth units. |
| One surviving parent household's one dependent leaves | \(H=1,C=1;M=1,R=X=0\) | \(H'=31/21=1.476,C'=0\); \(Q'=Q=3.1\). This is a conditional realized branch, not an assertion that departure probability is one. |
| One household exits with two dependents | \(H=1,C=2;X=1,R=2,M=0\) | \(H'=20/21=0.952,C'=0\); \(Q\) drops by 2.1 for the exiting household. No dependent also enters the ordinary stream. |
| Half of one household mass with two dependents survives | \(X=1/2,R=1,M=(1/2)2(2/9)=2/9\) | \(F=11/9\), \(C'=7/9\), \(H'=409/378=1.082\). Using \(2/9\) on the full original child stock would double-count departures from dead households. |
| Third birth; only one of the first two children remains | \(m:1\to2\), \(z:0\to\delta\) | Weighted stock rises by exactly \(1+\delta=1.602\). If one of the two slots leaves, its conditional surviving weighted stock is \(1+\delta/2=1.301\). If the parent instead exits, all \(2+\delta=2.602\) units are released. Earlier departures never change weight. |
| Terminal-age parent with two remaining dependents | Set \(s_J=0\) | Same death-release case; no separate ordinary draw and no dependent tail left behind. |

The exact test enumerates all departing subsets for one, two and three represented children, both tagged statuses, and survival probabilities zero, one third and one: 18 cases. It also checks 36 state-mixture/choice-probability combinations, the normalized and literal identities with retention/outside inflow/age-bridge adjustments, one terminally closed child cohort, and estate nonduplication. These are pure arithmetic checks, not tests of an implemented model change. `hand_checks.json` contains exact fractions and unrounded values.

## 5. Initial state, legacy queue, age meaning and terminal state

**Initial state.** Build \(g_0\) and, for weighted accounting, \(z_0\) under the new operator, using the matching **pre-fertility** distribution: the stationary forward sweep mutates its distribution through births, so an unlabeled saved `sol.g` is not automatically the required starting phase. If the old reference is stationary, replay its saved household policies from childless entry through adult ages to recover \(z_0\) and compare the unified implied entry with its existing entry. This is a forward calculation, not a Bellman solve. If it reproduces \(g_0\), the market accounts, and \(E=cF+I\), the reference equilibrium can be reused. There is no need to assert in advance that every old equilibrium is invalid. At the 2.1 replacement reference this equality can hold exactly after top-bin and death release are made consistent.

An arbitrary dated checkpoint's \(g_0\) does not identify \(z_0\): two households with the same \((n=3,m=1)\) may have different surviving-child identities. Recover the moment using saved past transports/birth branches, or replay the historical path from a reconstructed stationary initial state. A uniform-tag imputation is a possible explicit approximation, not exact recovery.

**Legacy queue.** The old birth queue is a second entry clock, not a separate verified resident population. Remove it from entry production and from the new checkpoint schema; preserve old arrays as audit metadata only. Never add its due arrivals to \(c(M+R)\), and never subtract its dollar/household units from \(m\). The old code allows queued adulthood and dependent residence to overlap or leave gaps, so there is no exact person-by-person splice from just the old \((g,\text{queue})\). A new-model transition starts from the explicitly reconstructed state; a consistent historical comparison replays that history. [S3,S4]

**Age.** The retained geometric law implies a first departure after four years and mean \(4/\mu=18\) years absent parental exit. Compulsory release weakly shortens that mean. If a birth occurs at parental age \(a\), its probability of departure on step \(k\ge1\) is
\[
\Pr(T=k\mid a)=\left[\prod_{r=0}^{k-2}s_{a+4r}(1-\mu)\right]
\left[1-s_{a+4(k-1)}(1-\mu)\right].
\]
The bracket splits into death release \((1-s)\) and ordinary departure \(s\mu\). Terminal \(s=0\) makes the probabilities sum to one. No recalibration of \(\mu\) is implied. Under the historical support, children born by parental age 42 are already at least 24 when the parent first faces mortality at 66; a release during 66→70 enters the next date when the youngest would be 28. Calling every entrant model age 18 is therefore an **economic-lifecycle reset**, affecting earnings, fertility opportunities and remaining lifetime. It is already implicit for early geometric departures if they become entrants. Exact chronological ages cannot be recovered from \((n,m)\) alone. This assessment retains the requested state-space restriction and does not propose reopening it. [S6,S7]

**Terminal state.** Replace queue-gap tests with joint \((g,z)\), total-scale, entry, and dependent-stock convergence tests against a terminal state solving the same closure. Distinguish terminal adult age, where all remaining attached children are released, from the finite calendar horizon, where they remain in the endpoint distribution and must not all be forced out. The household terminal value, fiscal closure, and last endogenous price must correspond to the new demographic tail. A prescribed terminal value alone does not prove convergence. [S4]

## 6. Inheritances, wealth and current household values

The source's optional estate pool calculates death-weighted positive net estates and pays eligible existing ages 45–65; entrants draw the existing exogenous wealth distribution. Activation must be read from the actual parameter object. Preserve these rules. The inspected stationary solver explicitly closes the estate pool, whereas the two inspected dated drivers contain no explicit estate-transfer fixed point. If the selected run activates the pool, verify its downstream dated wiring or add the dated generated-versus-paid equation; carrying a stationary lump sum through changing deaths would not establish that identity. This is conditional implementation work, not evidence that every existing run activates or violates the pool. A warm-glow bequest utility term is a preference, not another cash payment. [S2,S5]

Let \(\mathcal E_t\) be net estates, \(T_t^{\rm pool}\) their pooled payments, \(T_t^{\rm new}\) any separately authorized payments to entrants, and \(L_t\) explicitly specified taxes/other disposition. Enforce \(\mathcal E_t=T_t^{\rm pool}+T_t^{\rm new}+L_t\). Under the retained pool, \(T_t^{\rm new}=0\). A household exit with estate 150 and three eligible recipient households pays 50 each, totaling 150. Its two released dependents create \(20/21\) new households; they do **not** each receive 150 again. If entrants' exogenous mean beginning wealth is 7, their separate asset inflow is \((20/21)7=20/3\). Record the existing funding/boundary convention for this inflow; calling it a second inheritance does not fund it. If that convention is absent, it remains a preexisting resource-closure item, not a reason to invent parent-child wealth matching.

At fixed prices, transfers and the existing warm-glow rule, release-at-death does not alter surviving parents' child kernel, utility, budgets or continuation values. It changes who enters the next period. The tagged moment likewise never enters household optimization. Thus household policy arrays should reproduce at identical aggregate inputs. New entrant numbers and timing can subsequently change housing demand, taxes, pension recipient/contributor masses, pooled-transfer denominators, prices and household policies in equilibrium. These are economic effects of the new entry law, despite the small accounting implementation.

## 7. Targeted code map and validation gates

| Source and location | Required scoped change |
|---|---|
| `intergen_eqscale_seq_optimized/parameters.py`, binomial kernels; `e5_maturation_repair.py:17–33` | Preserve \(\mu\), existing child kernel, newborn timing and household state dimensions. Add a checked demographic rule identifier and conversion/unit metadata rather than silently reusing conflicting defaults. |
| `intergen_eqscale_seq_optimized/e5_profile.py:175–184`; `solver.py:4625–4629,5619–5658` | The profile counter currently uses \(1/2\), and counts ordinary departures only. Expose **raw** ordinary and death-release flows plus weighted counterparts; apply the single retained conversion once. Cover both Markov/non-Markov and any enabled compiled paths. Do not mistake the counter for the active stationary closure. |
| `solver.py:924–1056,1070–1082,4595–4610,5689–5725` | Reconcile actual normalized entry and open/closed scale closures with the same unified flows. Add terminal release and tagged forward moment. Scale every stock, flow and auxiliary moment consistently if a stationary distribution is normalized. |
| `tools/run_e5f_open_population_transition.py:740–785,828–921,1214–1267` | Keep birth measurement, add weighted child ledger, record \(R\) before survival deletion, return \(M,R\), and replace both operational birth queues with the unified entry. Preserve explicit retention/outside inflow and measure bridge residuals in both ledgers. |
| `tools/run_e5f_perfect_foresight_transition.py:73–78,346–445,650–700,805–806` | New initial/checkpoint state carries \(g,z\); use same operator as the temporary-equilibrium path; remove queue advancement and queue convergence criterion. Keep death and fiscal timing aligned with each dated policy. Version schemas to reject ambiguous old restarts. |
| `tools/run_dynamic_population_transition.py:371–390` | Reuse existing entrant distribution for both streams. No new inherited-wealth or parental-origin kernel. |
| `solver.py:1129–1188,1263–1306` | Retain estate rules; verify deaths and recipient dates agree with forward accounting and no estate is paid twice. |

Before a production edit, pin actual serialized parameters **and** adapter overrides; constructor defaults are insufficient. Match the retained/frozen branch instead of modifying the September 14 reference or a running bundle. [S7]

Minimal meaningful checks after implementation: (a) the deterministic cases above through the actual shared operator; (b) household, literal-state child and weighted-child residuals at every stage, including terminal age and \(m=0\); (c) nonnegativity and tag bounds; (d) ordinary/death disjointness and exactly one entrant insertion; (e) fixed-input household-policy reproduction; (f) new stationary operator applied once reproduces its stationary state, with \(F=B\) in the same units; (g) dated versus stationary identical-policy agreement, including retention/bridge disabled; (h) estate generated-versus-paid and entrant-wealth boundary accounts; (i) explicit rejection of legacy queue restarts; (j) source/parameter hashes. Only then perform one fixed-policy path comparison, a fresh equilibrium/path where needed, and the standard diagnostics. Recalibration is a subsequent decision based on the actual changed fit.

## 8. Costs, measurable effects, and smallest decisions

**Implementation cost.** With literal units alone, the core is a few flow accumulations and replacement of queue insertion in two drivers. The primary weighted construction adds the tagged first moment and checkpoint/initial/terminal plumbing. This is a medium integration task across existing paths, not a new household model. Rough planning allowance: about 1–2 focused engineering days for a verified experimental implementation and saved-policy replay, contingent on frozen-source availability; numerical equilibrium/calibration time is additional. This is a scope estimate, not a measured runtime.

**Compute and memory.** No new Bellman state, household choice or optimization is introduced. Flow accumulation is \(O(N)\) in occupied forward-grid cells; the tag requires additional linear transport for only the top-state slice, with \(O(N)\) worst-case work and at most \(8N/4\) bytes per dense tagged array when \(n\) has four cells. A transition stores it per date only if restart/diagnostics require that; otherwise two buffers suffice. A generic full-array implementation might add roughly one forward pass; a top-slice implementation is smaller. No runtime multiplier is claimed without profiling. Fresh GE cost is driven by path length, market/fiscal iterations and numerical convergence, not an exponential expansion of household states.

**What saved outputs can establish without a new household solve.** Actual \(M,R,X,B\), age/location incidence, weighted flows after recovering the tag, \(c(M+R)\) versus the old scheduled queue, the first-step entry difference, and a forward path at frozen policies/prices are calculable. The historical September 12 receipt has literal \(R=0.005972\) per unit household mass and four-year step: at unchanged \(c=1/2.1\), release adds \(R/2.1=0.002844\) to an ordinary-only entry counter. This is historical and literal-state evidence, **not** the latest weighted effect or the difference from the old birth queue; that queue already schedules births after parental death. No fresh checkpoint was loaded here. [S6]

The weighted unified law can preserve the old stationary entrant level at the 2.1 reference while materially changing a transition's **timing**: a fertility perturbation can affect entries next period, with its subsequent lag distribution spread geometrically, instead of the existing queue's twenty-year effect lag. Its price, fertility, pension, wealth and welfare effects require a new consistent solve; their sign is not determined by bookkeeping.

**Recommended defaults for making this concrete:** retain \(1/2.1\) as replacement-normalized entry; retain the age-18 economic lifecycle and \(2/9\) ordinary departure; use the tagged top-bin moment; retain current estate/entrant-wealth rules; replace rather than supplement the queue; use one operator everywhere. The only small new author choice is whether to keep adjusted renewal units with that tagged representation or intentionally switch renewal to literal state births. A literal-census population model, cap, orphan-care state, genealogical inheritance or new mortality assumption is not required. The second material decision arises only when defining a production terminal policy closure: a closed low-fertility population cannot be assigned a positive constant-size stationary endpoint. Preserve the existing explicit open closure where applicable, and do not hide that choice in demographic normalization.

Evidence S1–S7, exact source hashes, derivation assumptions, and verification receipts are in `receipt.json`. The accepted household model, source files, canonical status and papers were not changed.
