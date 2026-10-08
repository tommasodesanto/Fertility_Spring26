# The sale-to-rent screen: what it is and what it does

Base 14.402, fixed price $P=0.7794$. Code is the engine the 14.402 solves import:
`tmp/rental_menu_precaution_20261007/root/code/model/experiments/birth_count_choice/model/engine/`.
Periods are four years. $R=1.0824$ per period, $\varsigma=0.06$ selling cost, rent $r=0.1405$ per room per period,
parent room floor $h_P=2.596$, owner rungs $h\in\{2,4,6,8,10\}$.

## 1. Timing within a period, as coded

An owner enters with financial position $b$ (negative = mortgage debt) and house $h$. Then:

1. Interest accrues on $b$, giving $Rb$ (`household.py:317`). There is one rate for both signs of $b$.
2. Income $y$ arrives. It includes the property-tax rebate.
3. The household picks a tenure and house. A transaction moves the sale proceeds $S=(1-\varsigma)Ph$ and the purchase cost $Ph'$ through the
   same budget.
4. It consumes and chooses end-of-period $b'$.

So the budget for every branch is $c + \text{housing cost} + b' = Rb + y + S\cdot\mathbb 1_{\text{sell}} - Ph'\cdot\mathbb 1_{\text{buy}}$.
The branch code passes $b + (S-Ph')/R$ into the saving step, which applies $R$ and adds $y$. That is identical.

## 2. The owner's options and their constraints

| option | eligibility check before saving | end-of-period floor | where |
|---|---|---|---|
| **stay** | none | $b'\ge\max\{\min(b,-\varphi Ph),\ \text{death floor}\}$ | `kernels.py:976-982` |
| **resize** to another owned $h'$ | $Rb+S+y\ge(1-\varphi)Ph'$ | $b'\ge-\varphi Ph'$ | `kernels.py:369-372`, `household.py:951-956`, floor `kernels.py:959` |
| **sell and rent** | $Rb+S\ge0$ (**no $y$**) | $b'\ge0$ (renters cannot borrow) | screen `kernels.py:341`, renter floor `household.py:47` |
| (renter **buys**, for comparison) | $Rb+y\ge(1-\varphi)Ph'$ | $b'\ge-\varphi Ph'$ | `kernels.py:363-366` |

Parents have one more restriction. A parent cannot live in an owned home smaller than $h_P$, so a 2-room owner is infeasible with a child at home (`kernels.py:923`).

The death floor $b'\ge-(1-\varsigma)Ph$ keeps the net estate non-negative. It applies only at ages where death is possible, so it plays no role at 22–33.

## 3. Why purchases count income and the sale screen doesn't

- The screen is not a general rule. It is switched on only in the "corrected" credit mode, i.e. when a renter borrowing limit $\bar d$ is set (`household.py:562`, `credit.py` `bind_engine_credit`). The production point sets $\bar d=0$. In the older "reference" mode the screen does not exist.
- The rule was written in the Sept 29 credit fixture under the original timing as $b+S\ge0$. Its stated purpose is sale solvency: an owner who sells must not become a renter whose debt exceeds what $\bar d$ allows. When interest moved first it became $Rb+S\ge0$. Income was never added (`credit.py` docstring).
- Purchases count income because of the purchase-income timing, which makes income available before the transaction. The same timing applies to resizes, but it was never extended to the sale screen. I found no written reason for the difference.
- With $\bar d=0$, the renter's own floor already enforces $b'\ge0$ after income. So without the screen, selling into renting still could never leave a renter in debt. The screen only removes moves that are budget-feasible.

The screen fails exactly when the entering LTV exceeds $(1-\varsigma)/R=86.8\%$, whatever the household's income. At 80% financing almost nobody enters above that: 1.6% of childless owners 22–33, an artefact of the asset grid. At 95%, 60.5% do.

## 4. Worked example (LTV 95 solve, age 26, median income)

A childless household bought the **2-room starter** last period at 95% financing. It is typical of the barred group, 88% of whom own 2 rooms.

- $Ph = 0.7794\times2 = 1.559$. It enters at the grid node $b=-1.512$, an LTV of 97%.
- $Rb=-1.636$ and $S=0.94\times1.559=1.465$, so $Rb+S=-0.171<0$: **it may not sell into renting.**
- With income $y=2.888$, $Rb+S+y=2.717$. Renting the parent minimum of 2.6 rooms costs $0.365$. **The move is affordable.**

What it can do instead:

- **Stay in 2 rooms.** This is fine while childless: interest-only, $b'\ge-1.512$. With a baby it is not allowed, because 2 rooms are below $h_P$.
- **Trade up to 6 owned rooms.** $Ph'=4.677$. Eligibility is $2.717\ge0.05\times4.677=0.234$, which passes. After the purchase, resources are $-1.636+1.465-4.677+2.888=-1.960$, upkeep is $0.458$, and the floor is $b'\ge-4.443$. So it can borrow its way into the bigger house.

What the model chooses (tenure probabilities):

| | waiting branch: rent / 2 / 4 / 6 rooms | birth branch: rent / 4 / 6 rooms | first-birth attempt probability |
|---|---|---|---|
| current screen | 0 / 22% / 51% / 27% | 0 / 2% / **98%** | 0.381 |
| screen with income | 46% / 11% / 27% / 15% | **75%** / 1% / 25% | 0.410 (+7.6%) |

So under the current screen, having a child forces this household to trade up to a 6-room owned home at about 95% leverage, with a 6% cost on any later exit. Under the income-counting screen, three quarters of these households would instead rent about 3–4 rooms and keep the option to adjust. Being forced into the bigger, illiquid, highly levered house makes the birth branch less attractive, so fewer try. This interpretation is my inference; the probabilities are measured.

Across all barred childless owners 22–33 at LTV 95, the mean attempt probability is 0.118 under the current screen and 0.142 counting income. That compares the same states, but the two solves have slightly different barred groups, so it is not an exact matched-state comparison. On the birth branch, 0% of the barred group rents under the current screen and 47% would rent with income counted. In the stationary solve, completed fertility under 95% financing goes from −1.35% to −1.06% relative to 80%. Births on impact do not change (+0.42% and +0.41%).

## 5. Two defensible readings

**A. A deliberate pre-income solvency rule.** On this reading, a transaction must be solvent before income arrives. Lenders require the old loan repaid at closing, and closing does not wait for the next paycheck.
- Then **purchases and resizes should obey it too**: $Rb\ge(1-\varphi)Ph'$ and $Rb+S\ge(1-\varphi)Ph'$, with no $y$.
- That tightens every down-payment test, which is the DUE / Sommer–Sullivan–Verbrugge convention noted on Oct 1.
- It would change ownership, the calibration fit and every credit result. It needs a refit.

**B. An omission.** On this reading, income counts for every transaction, as the purchase-income timing already says.
- Then the screen should be $Rb+S+y\ge0$. With $\bar d=0$ it is then redundant with the renter floor and could simply be dropped.
- It is inert at 80% financing (+0.00% here), so the production calibration is unchanged.
- It changes any experiment that pushes entering LTVs above 86.8%:
  - 95% financing: completed fertility −1.06% instead of −1.35%.
  - Transitions with house-price falls above about 8%: an owner at 80% then enters above the cutoff. This would hit the un-rebated property tax (price −19%) but not the rebated one (−0.7%). This is inferred, not checked.

Under either reading the two rules should match. Today purchases follow B and sales follow A.
