# Joint budget: cached implementation check

Reference: accepted soft post-interest fixed-price \(\phi=.8\) versus approved \(.95\), with identical public inputs except financed share. [Illustrative JSON](joint_budget.json) records array/input hashes and four matching production pins; [driver](extract_joint_budget.py) reproduces it. Standard graphs remain linked in [the diagnostic report](../README.md).

For every housing alternative, `_savings_stage` jointly optimizes consumption \(c\) and ending assets \(b'\), anticipating continuation utility. Tenure interpolates these **conditional optimized values**, without prior saving commitment. Fertility uses conditional tenure inclusive values for waiting and successful birth. Linear interpolation approximates off-grid optimized values. Nested taste shocks and birth-outcome contingency remain sequential; this is not a flat simultaneous logit. See [household.py:864–991](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/production/engine/household.py:864) and [exhaustive saving](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/production/engine/kernels.py:505).

For a purchase, \(Q=ph\) is its price, \(S\) the proceeds from selling inherited housing, and \(K=(\delta+\tau_H)Q\) owner carrying costs here. The implemented budget is

\[
Rb+S+y=Q+K+c+b',\qquad b'\ge-\phi Q,\qquad c>\bar c.
\]

The screen \(Rb+S+y\ge(1-\phi)Q\) is necessary, **not sufficient**: consumption, carrying costs, collateral and continuation value still constrain the choice. Purchase timing is [income-adjusted](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/production/engine/household.py:896). Current \(y\) is four-year after-tax income; beginning \(b\) is before interest. This interval accounting does not establish actual closing dates within the interval.

At age 22, income state 5, inherited renting and \(b=S=0\), \(y=2.061181\), \(h=4\), \(Q=3.104228\), \(K=.303740\). Baseline down payment is \(.620846\): screen slack is \(1.440336\), but maximum consumption at the collateral floor is only \(1.136596\), not the whole income. Saved controls are the exact implemented two-node transaction mixture:

| Financed share | Branch | Consumption | Ending assets | Floor binding | Owner-four-room probability |
|---|---|---:|---:|---|---:|
| .80 | WAIT | 1.136596 | −2.483382 | Both nodes | 4.762720% |
| .80 | SUCCESS | 1.136596 | −2.483382 | Both nodes | .012276% |
| .95 | WAIT | 1.272596 | −2.619383 | Neither node | 9.280753% |
| .95 | SUCCESS | 1.141585 | −2.488372 | Neither node | .014526% |

Here \(\bar c=0\). No consumption lower bound binds. At a binding saving floor, its local multiplier is \(\mu=U_c-\beta V'_{c,+}\ge0\). Holding continuation and menus fixed, the direct current borrowing-capacity derivative is \(Q\mu\). Baseline conditional derivatives are **.395485 WAIT and .364198 SUCCESS**, a child-relative difference of **−.031288**; weighting by branch probabilities gives **.018836 and .000044708**, respectively. This establishes positive current borrowing-capacity value for that same house, not a larger child-relative gain. The full permanent-policy successful-birth-minus-wait value change is **−.00170830** here; it also includes future credit, other housing alternatives and probability changes. Neither local derivative identifies that total causal attribution.

Correction to the provisional capture concern: executed inputs omit `child_maturation_mode`; the verified helper defaults to **constant**, so newborn-exempt `VI_ex` is inactive. Saved \((n,m)=(1,1)\) controls represent SUCCESS. Canonical transition normalization changes entries by at most \(1.11\times10^{-16}\), exactly reproducing the saved transition. All 48 illustrative endpoint budgets and saving inequalities pass; consumption surplus remains positive.

The [age-22 renter extension](age22_renter_envelope.json) ([per-state CSV](age22_renter_envelope.csv)) audits all 249 baseline susceptible nodes and five owner products. Weight \(W=G_{\mathrm{pre}}\pi^2a p_{\mathrm{wait}}/\kappa\) measures first-birth-flow sensitivity to successful-birth-minus-wait utility. The current borrowing-only tenure derivative is \(\sum_h q_hQ_h\mu_h\); the rental-floor derivative is zero here. Future continuation, menus and screens are held fixed.

| Earnings group | Baseline W | Weighted current child-relative derivative per unit φ | Permanent .80→.95 child-relative value change |
|---|---:|---:|---:|
| All | .0297661 | −.000590073 | −.000441432 |
| Low | .0004389 | −.000614889 | −.001518351 |
| Mid | .0181506 | −.004094423 | −.000738991 |
| High | .0111766 | +.005101908 | +.000084085 |

This slice covers **21.5401%** of all-age/all-beginning-tenure first-birth W. Every positive-q branch passes budgets and floor/interior/kink inequalities: uncovered positive-qW is zero. Zero-q branches have **unmeasured μ**; their numerical fixed-menu contribution is zero, without excluding menu-opening effects. They touch 52.1511% of WAIT W and 100% of SUCCESS W. Saved probability sums drift by at most \(5.57\times10^{-8}\); retained probabilities reproduce the executed policy, without renormalization. Using both saved birth probabilities gives the same total W as \(a(1-a)\) at reported precision.

The current local first-birth-flow derivative is **−.0000175642 per unit φ**. Selected owner-floor-binding probabilities are 3.1678% WAIT and 1.5961% SUCCESS (10.7423% and 7.7030% conditional on ownership). Permanent gap changes use identical fully supported W states. This is not a permanent-policy causal decomposition or welfare claim. No solve or capture was added; broader attribution remains unresolved.

Regenerate the 17 standard graphs from saved .95 arrays/scalars through unchanged canonical `write_diagnostics`:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python output/model/credit_mechanism_20261004/diagnostics/joint_budget/age22_renter_envelope.py --standard-graphs write
```

[Extension driver](age22_renter_envelope.py) without arguments regenerates CSV/JSON. `--case PATH --output PATH` selects another saved case and an owned graph-replay folder. Cached canonical summary reconstruction passed exactly. The lead ran the command above: all 17 standard graphs regenerated successfully without a model solve.
