# Saved-case calibration fit diagnosis

This is a read-only check of the **80% financed-share** calibration. It uses the fresh-postchecked hard selected point (Torch continuation chain 11, loss 97.011220) and quarter selected point (local chain 54, loss 51.556036). The full [hard target table](hard_verified_best_chain11_target_fit.csv), [hard 31-parameter table](hard_verified_best_chain11_parameters.csv), [quarter target table](quarter_verified_best_chain54_target_fit.csv), and [quarter 31-parameter table](quarter_verified_best_chain54_parameters.csv) are saved alongside this report. The early-fertility target is children ever born by age 25, capped at three, as implemented in `code/model/tools/e5f_calibration_runtime.py:95`. Neither selected fit has optimizer convergence certification.

| Loss block | Hard contribution | Hard share | Quarter contribution | Quarter share |
|---|---:|---:|---:|---:|
| Mean rooms, ownership at ages 30–55, first-birth rooms, recent-parent ownership | 85.155 | 87.8% | 38.176 | 74.0% |
| Early fertility | 8.128 | 8.4% | 8.315 | 16.1% |
| Other scored moments | 3.728 | 3.8% | 5.065 | 9.8% |
| **Total** | **97.011** | **100%** | **51.556** | **100%** |

| Moment | Target | Hard model | Quarter model |
|---|---:|---:|---:|
| Mean occupied rooms | 5.729 | 6.241 | 6.131 |
| Ownership, ages 30–55 | 0.676 | 0.609 | 0.646 |
| Room increase at first birth | 1.465 | 1.063 | 1.167 |
| Recent-parent ownership difference | 0.128 | 0.101 | 0.117 |
| Early fertility | 0.810 | 0.524 | 0.521 |

For an exploratory persistence check, I read only four saved case lists: hard continuation chain 11 and local chain 51; quarter original Torch chain 35 and local chain 54. A case counts here only if its status is `passed`, its base loss is finite, and it has a 14-row target table. “Near best” means loss no more than 25% above the selected winner within the same purchase rule. These are repeated optimizer evaluations, not independent draws.

| Rule | Near-best cases | Mean rooms model range | First-birth room increase range | Early-fertility model range |
|---|---:|---:|---:|---:|
| Hard | 98 | 6.136–6.451 | 1.005–1.095 | 0.521–0.546 |
| Quarter | 171 | 5.986–6.214 | 1.135–1.197 | 0.521–0.529 |

Every case in these near-best subsets has mean rooms above its 5.729 target, first-birth room increase below 1.465, and early fertility below 0.810. This establishes persistence in the sampled neighborhoods, not that the targets are unreachable. The search varies ten parameters together and these cases do not identify a causal response to any single parameter.

There is no universal first-birth-versus-mean-rooms tradeoff in these receipts. The hard local chain 51 case `0066_nm` moves first-birth rooms from 1.063 to 1.093 and early fertility from 0.524 to 0.540, while mean rooms rises from 6.241 to 6.352; its loss is 98.291. The quarter original chain 35 case `0062_nm` raises first-birth rooms from 1.167 to 1.197 and improves mean rooms from 6.131 to 6.103; early fertility rises only from 0.521 to 0.523, and loss is 51.970. Thus the housing moments can move favorably together, while the early-fertility gap remains large in these examples.

Evidence: hard continuation `search/cases.json` at `/scratch/td2248/projects/purchase_restart_controller_v2/results/chain_11/`; quarter original `search/cases.json` at `/scratch/td2248/projects/purchase_rules_overnight_v1/results/chain_35/`; local case lists at `local_runtime/runs/local10_v1/chain51/search/cases.json` and `chain54/search/cases.json`. No model was rerun.

## Age-25 fertility accounting

The age-25 target is **children ever born, capped at three**, not the share who have become mothers. Let $m=\Pr(N_{25}\geq1)$ and $s=\mathbb E[N_{25}\mid N_{25}\geq1]$. Then $\mathbb E[N_{25}]=m s$. The [measurement audit](../../../fertility_identification_20260928/measurement_audit_v1/early_fertility_decomposition.json) gives data motherhood $m=0.457254$, children per mother $s=1.770410$, and total $0.809528$. The selected [hard observer](../collection/readout/hard/selected_root/observers.json) gives $m=0.437884$, $s=1.197637$, and total $0.524426$; the [quarter observer](../collection/readout/quarter/selected_root/observers.json) gives $m=0.434716$, $s=1.198878$, and total $0.521171$.

A symmetric product decomposition writes the data-minus-model gap as

$$\Delta(ms)=(m_D-m_M)\frac{s_D+s_M}{2}+(s_D-s_M)\frac{m_D+m_M}{2}.$$

| Rule | Total gap | Motherhood part | Children-per-mother part | Share from children per mother |
|---|---:|---:|---:|---:|
| Hard | 0.285101 | 0.028745 | 0.256356 | 89.9% |
| Quarter | 0.288357 | 0.033462 | 0.254895 | 88.4% |

This is **accounting, not a causal decomposition**. It shows that the selected points get motherhood closer than they get the number of children among mothers. The saved observer metadata allows at most one explicit birth transition per four-year cell. Both selected fits put exactly zero age-25 mass in the three-or-more-children state. Timing, spacing, and the effective age-25 support therefore deserve a focused check. The data value of 1.770 children per mother is below two, so these receipts do not establish that the target is mechanically unreachable.
