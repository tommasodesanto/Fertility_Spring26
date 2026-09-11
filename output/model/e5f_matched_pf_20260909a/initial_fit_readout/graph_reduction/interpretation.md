# Occupied ownership-probability reversals: CSV interpretation

**The reversals are not confined to empty states, but the largest visible jumps have extremely small mass.** A targeted check of one conditional slice is justified; these CSVs do not establish an optimizer error. The extracted ages are **30 and 42**, not 50. No new solve or checkpoint load was performed.

All six CSV SHA256 hashes match the extraction receipt. Independent reconstruction reproduces every age/type drop count and lower-node mass in the saved screen; permanent-group rates also reproduce from the combined-type owner numerators and mass denominators.

## Definitions and scope

A drop is an adjacent-wealth decrease in owner probability exceeding $10^{-7}$, with both endpoints passing the original graph mask. An occupied drop has lower-node **post-fertility/pre-tenure** mass above $10^{-12}$. This is the appropriate population for the conditional childless tenure decision. Each lower node is counted once. The selected slice pools all 15 earnings states at renter tenure, location 0, parity 0 and child state 0, at the stated age(s). Global shares divide by total household mass, computed from all ages and all 15 combined states; this is approximately one in each case.

These are shares of lower nodes associated with a probability reversal, not losses of owners, changes in equilibrium ownership, or welfare effects. The two displayed age slices are not a comprehensive audit of all ages or family/tenure states.

## Occupied-node results

| Case | Age | All valid drops | Occupied drops | Lower-node mass | Share of selected slice | Share of global households | Largest occupied drop |
|---|---|---:|---:|---:|---:|---:|---:|
| old_old | 30 | 361 | 252 | 0.002537497 | 18.3117% | 0.253750% | 7.1647 pp |
| old_old | 42 | 355 | 302 | 0.001427826 | 26.6439% | 0.142783% | 7.6724 pp |
| old_old | 30 + 42 | 716 | 554 | 0.003965323 | 20.6353% | 0.396532% | 7.6724 pp |
| new_balanced | 30 | 361 | 253 | 0.002198857 | 16.0095% | 0.219886% | 6.1836 pp |
| new_balanced | 42 | 331 | 278 | 0.001249656 | 23.1079% | 0.124966% | 7.8354 pp |
| new_balanced | 30 + 42 | 692 | 531 | 0.003448513 | 18.0148% | 0.344851% | 7.8354 pp |

The two selected slices contain 1.921618% of global households in `old_old` and 1.914265% in `new_balanced`. At least one-percentage-point drops affect lower-node mass 0.000619465 in `old_old` (125 occupied intervals; 0.0619465% of global households) and 0.000464464 in `new_balanced` (109 intervals; 0.0464464% globally). The full counts above include much smaller reversals.

The largest occupied decline in `old_old` is 7.6724 pp, at age 42, earnings index 2, wealth 5.465116 to 5.604651; lower-node mass is only $1.9240\times10^{-9}$. In `new_balanced`, its counterpart maximum is 7.8354 pp at wealth 4.906977 to 5.046512, same age/index, with mass $5.5122\times10^{-9}$. The largest unrestricted graph drops are 27.7436 pp and 25.9278 pp respectively. They are not the economically most exposed states.

## Three concrete states in the revised, balanced case

Indices below are zero-based. All three have renter/location0/parity0/child-state0. The saved inclusive value increases across each pair; that does not certify every conditional branch action.

| Priority | Age; earnings index | Wealth pair | Owner-probability pair | Decline | Lower-node mass | Global mass share |
|---|---|---|---|---:|---:|---:|
| 1 | 30; 8 | 0.860465 → 1.000000 | 27.9660% → 26.6428% | 1.3232 pp | 0.00010059048 | 0.01005905% |
| 2 | 42; 4 | 2.255814 → 2.395349 | 39.5618% → 34.3638% | 5.1980 pp | 4.50596987e-06 | 0.00045060% |
| 3 | 42; 2 | 4.906977 → 5.046512 | 43.4035% → 35.5681% | 7.8354 pp | 5.51220263e-09 | 0.00000055% |

1. **Start with age 30, index 8, wealth indices 52–53.** This state has the largest mass-times-drop score in the revised case; its lower node is 0.732380% of the selected age-30 childless-renter slice. Earnings combine middle permanent income 0.821954 and changing state 3, giving $z=1.143228$. Conditional rental housing is 6 rooms at both endpoints; the deterministic selected tenure is renter at both. The six-room owner probability falls 27.7012% to 26.3760%, accounting for essentially the total reversal. At price 0.633328602 and financed share 0.8, the 6- and 8-room down-payment thresholds are 0.759994 and 1.013326. Both wealth nodes are between them: **this decrease does not cross a down-payment threshold**. The same interval in `old_old` also declines, 30.8571% to 29.4385%, with lower mass 0.0001040367.

2. **Age 42, index 4, wealth indices 62–63:** lower permanent income with the highest changing earnings state, $z=0.584903$. Rental housing rises 4.155082 to 4.211339 while the four-room owner probability falls 39.5585% to 34.3532%; renter remains selected. All five owner products clear their simple down-payment threshold at both wealth nodes. Thus immediate product affordability entry does not explain this 5.1980-pp decline either.

3. **Age 42, index 2, wealth indices 81–82:** maximum occupied decline, but tiny mass. The four-room owner probability falls 43.3571% to 35.3808% while the six-room owner probability rises 0.0461% to 0.1869%. Rental housing rises 4.282884 to 4.327946 and renter remains selected. All products clear their simple down-payment thresholds; this is not a new affordability entry. This point is useful for shape diagnosis but is lower priority than state 1.

The five down-payment thresholds in `new_balanced` are 0.253331, 0.506663, 0.759994, 1.013326 and 1.266657 for 2, 4, 6, 8 and 10 rooms. These statements concern the down-payment restriction alone; the CSV does not certify all budget/continuation feasibility conditions for every branch.

## The permanent-income gradient has late-life exceptions

In both cases ownership strictly rises low → middle → high permanent income at **every age node from 18 through 70**. At ages **74, 78 and 82**, middle permanent income has higher ownership than high permanent income. At age 82 in `old_old`, low also exceeds high. Therefore the earlier positive-gradient statement is valid for pooled prime-age ownership, not universally across the lifecycle.

| Case | Age | Low permanent | Middle permanent | High permanent | Middle − high |
|---|---:|---:|---:|---:|---:|
| old_old | 30 | 1.4075% | 34.9579% | 87.3515% | -52.3936 pp |
| old_old | 42 | 8.6364% | 69.8093% | 97.7884% | -27.9791 pp |
| old_old | 74 | 89.8010% | 99.4667% | 98.3809% | +1.0859 pp |
| old_old | 78 | 97.6050% | 99.7631% | 98.3380% | +1.4251 pp |
| old_old | 82 | 99.5590% | 99.9275% | 98.8105% | +1.1171 pp |
| new_balanced | 30 | 0.8142% | 33.2070% | 87.0038% | -53.7968 pp |
| new_balanced | 42 | 4.4189% | 67.0535% | 97.5852% | -30.5317 pp |
| new_balanced | 74 | 78.0673% | 98.4898% | 98.2913% | +0.1985 pp |
| new_balanced | 78 | 92.8123% | 99.4265% | 98.2548% | +1.1716 pp |
| new_balanced | 82 | 98.5092% | 99.8523% | 98.7406% | +1.1117 pp |

These are mass-weighted three-group rates, not averages of the 15 displayed line heights. They establish an observed late-life crossing, not its economic explanation. The combined income index still mixes permanent type with the changing earnings component; no claim of a globally monotone 15-state ordering is supported.

## Recommended next decision

Run one bounded saved-continuation audit of **state 1 at its two adjacent wealth nodes** in `new_balanced`, without equilibrium reclearing or recalibration. Re-evaluate the renter and all affordable owner-product branch optima under the same saved prices, continuation values, constraints and exhaustive saving specification. Cross-check feasible stored/cross-node actions, and report branch values, saving/consumption/housing, best-minus-stored objective gains, and reconstructed tenure probabilities. This separates a valid relative-value reversal from a suboptimal saved action. Merely repeating the probability range or inclusive-value monotonicity test would not answer that question.

A nonmonotone probability is economically possible with discrete owner products and continuous rental housing; the CSV evidence has not established that explanation. It also has not established an optimizer error. The first audit should remain confined to this one conditional slice; no broad well-behaved-policy claim is justified. No audit job is launched by this interpretation.

Full precise summaries, three states, every age-specific permanent-group rate, threshold values and source hashes are in `interpretation.json`. The raw six CSVs and original `receipt.json` remain unchanged.
