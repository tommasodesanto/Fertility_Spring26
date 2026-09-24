# How large is the entry-timing change, and what is the narrow alternative?

**September 24, 2026.** Pure timing arithmetic plus reused compact saved evidence; no household/equilibrium solve, checkpoint load, source change or new adopted rule. This assessment qualifies the earlier preference for unifying departure and entry: **the change is large in entry timing, even though the mean delay moves only two years. I would preserve chronological entry for now and make the separate dependency interpretation explicit, rather than adopt full coupling as a population-accounting repair.** That preserves an approximation; it does not establish literal residence/resource consistency.

## Measured kernel change: ordinary departures are the main issue

The inspected dated driver has four waiting slots, pops the oldest value before appending this date's births, and inserts the popped mass at the *next* date. Thus a date-zero birth pulse first produces entrants at date five, twenty years later. The parser help and saved-metadata formula explicitly state this convention. The perfect-foresight driver uses the same helper. This is the reviewed dated implementation, not a claim about every historical model. Source locations and hashes are in `timing_assessment.json`.

Keep the accepted conversion \(1/2.1\) on both alternatives. Normalize a birth pulse to produce one eventual entrant-household unit, so the conversion cancels. With the unchanged departure probability \(\mu=2/9\), coupling departure to entry gives a delay \(T=4K\), with \(K\) geometric on \(1,2,\ldots\). Ignoring parental death to isolate this clock change,
\[
\Pr(T\le4k)=1-(7/9)^k,\qquad
\Pr(T=4k)=(2/9)(7/9)^{k-1}.
\]

| Years since the birth pulse | Coupled entry: cumulative entrants (%) | Existing twenty-year queue (%) |
|---|---:|---:|
| 4 | 22.222 | 0 |
| 8 | 39.506 | 0 |
| 12 | 52.949 | 0 |
| 16 | 63.405 | 0 |
| 20 | 71.537 | 100 |
| 28 | 82.782 | 100 |
| 40 | 91.899 | 100 |

The coupled delay has mean **18 years**, standard deviation **15.875 years**, and median **12 years**; the queue has mean/median twenty and zero dispersion. **63.405% enter before eighteen; 28.463% enter after twenty.** These are analytically implied shares of the illustrative birth pulse, not observed household shares or an estimated transition outcome. Under the reviewed mortality support, death release does not affect the fractions through twenty years.

For a **permanent birth increase** normalized to yield one extra entrant-household unit every four years in the long run, the coupled extra entry flow equals the cumulative pulse column divided by 100: 0.222 after four years, 0.634 after sixteen, 0.715 after twenty. The queue gives zero before twenty and one thereafter. Across all affected birth cohorts, cumulative extra entries by year twenty are **2.496 versus 1.000**. This is a kernel calculation at a fixed birth path, with no descendant-fertility or equilibrium feedback. It shows why similar mean durations do not imply similar short-run housing or population dynamics; the sign and magnitude of those economic effects remain unmeasured.

The first-period 22.222% departure probability already exists. What coupling newly does is assign that departing mass to an adult decision household four years after birth, initialized at model age eighteen. This is a structural timing approximation, not a newly discovered probability error or solver failure.

## Parental death: smaller flow, and no missing queued adulthood under the reviewed support

The saved September 12 stationary diagnostic—not the latest selected checkpoint—reports literal dependent flows per four-year step and unit household mass:
\[
M=0.109302306,\quad R=0.005972242,\quad B=0.115274548.
\]
Here \(M\) is ordinary departure conditional on parental survival and \(R\) is remaining dependents removed with parental exit. **\(R/(M+R)=5.181\%\)** of this historical departure flow; it is also 1.200% of the contemporaneous dependent stock. It is **not** a current cohort probability, a share of all households, or the increase versus the existing birth queue. At fixed conversion it would add 0.002844 to an ordinary-only entrant counter, which is a different comparison.

The decisive timing fact is that the diagnostic's latest birth occurs at parental node 42, while its first positive mortality applies at node 66. Consequently, even the youngest attached offspring is already **24** at the start of the first possible death period and **28** at the next-date entry boundary. Its twenty-year queue entry already occurred when its parent reached 62. Earlier births entered even earlier. The current source applies survival at the current age node and inserts entrants next date, confirming this alignment; terminal parental exit is later still.

Therefore **keeping the existing queue and adding only death-triggered entry for offspring who have not already entered adds exactly zero households under this reviewed age support**. The queue does not discard offspring when their parent dies. Adding entrants for every remaining dependent at parental death would instead create a second entrant representation of that birth vintage. This zero-effect conclusion is conditional on maintaining the documented 42/66 support and twenty-year queue; the latest selected parameter object has not been freshly reopened here.

For scale only, using the actual historical survival schedule with births individually placed at parent ages 18, 30 and 42 gives death-release route probabilities 2.320%, 4.930% and 10.479%, and mean unified entry delays 17.675, 17.310 and 16.533 years. These are **age-specific illustrative calculations without birth-age weights**, not aggregate estimates or population bounds. Every death-triggered release in these examples occurs after the queue's entry date. Forced release mainly shortens the far tail of the geometric process; it does not explain the much earlier mass of coupled entrants.

## Three possibilities, with a recommendation

| Option | What it changes or preserves | Remaining approximation and cost |
|---|---|---|
| **1. Preserve chronological adult entry; separate the demographic ledger from dependency. Recommended narrow default.** | Retain the twenty-year queue, \(1/2.1\), current departure probability, policies and decision states. Record every not-yet-entered birth in an aggregate pre-entry register, regardless of whether it remains attached to its parent. Parent death does not cancel that register. Under reviewed timing, no extra death entry is needed. | \(m\) becomes an economic-dependence/support indicator rather than an exclusive count of people who have not yet formed a household. If it must mean literal co-residence, renaming does not repair residence, housing and resource use. Small accounting/reporting cost; no new household solve merely to expose the registers. |
| **2. Adopt one departure-to-entry clock.** | Replace the queue with ordinary plus death-triggered departures. This is the previously derived construction; it keeps \(\mu\), \(1/2.1\), and decision states. | Accept the wide entry-age distribution quantified above and the adult-age reset. Tagged top-bin forward accounting is sufficient for aggregate units; new transitions/fit must be assessed. No household-state expansion, but a substantive entry-timing change. |
| **3. Impose an actual minimum or narrow range for dependent-child departure ages.** | A split *adult-entry* queue—for example half at sixteen and half at twenty—can retain an eighteen-year mean using only forward registers. But it still leaves the existing dependency clock separate. | Making the household's dependent count itself obey those ages generally requires richer household history: \((n,m)\) does not reveal which children are old enough. A parent-age-dependent hazard is feasible without new states, but cannot enforce a child's age threshold for all birth timings. This is a larger specification change, not the narrow next step. |

Children can remain individually counted under option 1: the register counts each represented child once in the pre-entry population, and \(m\) counts an individual support relationship. What cannot simultaneously be claimed without more structure is that \(m\) is an exclusive co-resident pre-adult population while those same people also generate independent adult households through the queue. The existing top-bin and replacement weights remain normalized quantities as before.

Let \(N_t\) count birth units that have not yet generated adult entries and \(A_t\) the birth units due for entry. The minimal queue ledger is \(N_{t+1}=N_t+B_t-A_t\), with household additions \(A_t/2.1\). Ordinary departures from \(m\) do not delete anyone from \(N\); parental death does not either. This closes a demographic replacement-unit ledger without counting \(m\) a second time. For a literal membership partition, distinguish attached pre-entry individuals, detached pre-entry individuals, and already-entered individuals still attached. A register can count these groups, but supplies no missing housing/consumption behavior for them. This is why the recommendation preserves chronology while explicitly retaining an economic approximation, rather than claiming a complete literal-person model.

A death-plus-queue hybrid can be deduplicated with **forward-only** information if younger mortality is introduced. For each not-yet-due birth cohort \(k\), carry an attached-child moment \(a_k(x)\) over the existing parent cells. Under symmetric thinning, transport it as \(sK(m'\mid m)(m'/m)a_k(x)\); death-triggered early entries remove \((1-s)a_k\) from the matching future queue slot. Ordinary early leavers remain in a detached waiting register until their due date; already-entered offspring never create another entrant. At most five young birth-age buckets plus short aggregate queues suffice for the twenty-year horizon. These moments add linear forward memory/transport, not Bellman decision states, because existing policies ignore child identity/age. They do **not** make age-restricted dependence Markov in \((n,m)\), or fix the resource interpretation. Under the reviewed age support, the early-death channel is identically zero, so this machinery is unnecessary now.

**Narrow next step:** keep the current entry clock provisionally; use the selected saved birth path to compare queued and departure-based entrant paths at fixed policies/prices, before deciding whether the large entry-timing change is economically intended. Recheck only the selected parameter object's birth cutoff, survival and queue length. If the author requires \(m\) to remain literal co-resident children, isolate that requirement as a modeling question; do not claim an accounting label has solved it. The unrelated subreplacement terminal-endpoint issue remains separate and is not a reason to choose either clock.

Reproduce: `python3 output/model/demographic_accounting_review/maturation_derivation_20260924/timing_arithmetic.py`. Exact rational checks validate the geometric moments, pulse/step kernels and queue lag; subset enumeration validates the conditional cohort-moment formula. JSON pins the compact historical receipt and inspected source hashes. No implementation has been adopted.
