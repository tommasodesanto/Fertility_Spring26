# Verified-winner borrowing diagnostic

Reference: verified floor-search chain 7, case 0173_nm, original 14-target loss 31.28400725566496, \(q_0=0.719168368828958\). All 31 recorded effective parameter fields, the 120×9 grid, and entry distribution are fixed; the GE comparison estimates none of them. The calibration contract has 10 free parameters, while the report carries 31 stored parameter fields. Experimental borrowing changes only: disable the DUE stayer rule, clear the fixed unsecured credit limit, and enable the existing native lifetime-solvency rule. The estate/fiscal/entry contracts and original target weights remain fixed. This experiment is **not** a calibration or an adopted model specification.

The six price cells hold \(q\) at \(0.99q_0,q_0,1.01q_0\). Their PRE-impact columns apply each policy to the *same inherited reference distribution*; their stationary-cohort columns allow the cohort distribution to change. These prescribed-price cases are not general-equilibrium roots. All three expanded-credit cells pass the occupied-state support guard, while unoccupied continuation alternatives and full natural support remain unverified.

| Stationary-cohort object | Verified reference | Credit at prescribed \(q_0\) | Credit GE |
|---|---:|---:|---:|
| Owner asset price \(q\) | 0.71916837 | 0.71916837 | 0.65505720 |
| Unit rent / owner user cost | 0.12965124 | 0.12965124 | 0.11809332 |
| Population scale \(N\) | 0.92226157 | Not closed | 0.81992468 |
| Completed fertility | 2.100000 | 2.002188 | 2.100000 |
| Childlessness | 0.202533 | 0.244450 | 0.216278 |
| Mean age at first birth, target measure | 25.9564 | 25.6285 | 25.3029 |
| Children ever born by age 25, capped at 3 | 0.531218 | 0.539120 | 0.584755 |
| Rooms per household | 5.977125 | 5.965737 | 6.339064 |
| Ownership, all ages | 0.659374 | 0.724575 | 0.748113 |
| First-birth rooms response, target measure | 1.223193 | 1.201782 | 1.181390 |

The GE/reference population ratio is 0.8890370, or −11.10%, using the same housing-supply-per-household closure. The credit GE selected point is \(q=0.6550572019731878\), \(N=0.8199246807842695\), and renewal residual \(3.7892859961\times10^{-9}\). Seven new native lifecycle solves were used, plus one prior presolve failed attempt reserved conservatively against the 12-attempt envelope. The root and fresh selected repeat match all 14 target CSV rows, 31 parameter CSV rows, and 17 standard plot SHA-256 hashes. The selected repeat also matched 107 numerical arrays. The job exited 0. Its result is a **support-limited stationary GE diagnostic**, not a full natural-support certificate, grid-convergence certificate, transition, or production adoption.

At \(q_0\), holding the inherited reference PRE distribution fixed, expanded credit changes births per household from 0.11525385 to 0.11859850, first births from 0.04991269 to 0.05273499, ownership from 0.65937360 to 0.70569611, and rooms per household from 5.97712540 to 6.11527112. In the separate prescribed-price stationary cohort, births per household move from 0.11525385 to 0.10974661 and completed fertility from 2.100000 to 2.002188. These are different estimands.

The central price elasticity below is \(\log(y_{1.01}/y_{0.99})/\log(1.01/0.99)\), computed within each arm; it does not trace GE equilibria.

| Outcome | Reference PRE impact | Credit PRE impact | Reference cohort | Credit cohort |
|---|---:|---:|---:|---:|
| Births per household | −0.5502 | −0.5826 | −0.5195 | −0.5368 |
| First births | −1.1162 | −1.1558 | −0.4062 | −0.4187 |
| Ownership | −0.2089 | −0.2875 | −0.3889 | −0.3768 |
| Rooms per household | −0.3343 | −0.3582 | −0.6616 | −0.6875 |
| Nonhousing consumption | +0.1224 | +0.1251 | −0.1148 | −0.1076 |

Evidence: [comparison.csv](comparison.csv), [all six fixed-price cases](responses/), [selected GE target table](ge_replacement_18958548/root_04/target_fit.csv), [selected GE parameter table](ge_replacement_18958548/root_04/parameters.csv), [17 standard plots](ge_replacement_18958548/root_04/standard_diagnostics/), [repeat receipt](ge_replacement_18958548/selected_repeat/receipt.json), and [GE completion receipt](ge_replacement_18958548/completed.json). The prior GE submission [failed before its first lifecycle solve](ge_failed_18958106/ge/lower_95/failure.json) because a driver called a nonexistent serialization API; the reviewed replacement used the original absolute deadline and preserved that failure evidence. No experimental result establishes broad economic infeasibility or optimizer convergence.
