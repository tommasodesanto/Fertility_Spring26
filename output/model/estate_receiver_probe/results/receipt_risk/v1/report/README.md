# Inheritance uncertainty: completed fixed-price test

Torch job 18571183; three household solves. Experimental, not recalibrated or funded general equilibrium.

![Supplemental lifecycle comparison](receipt_risk_lifecycle.png)

The no-receipt control exactly reproduces the retained net-valuation case. All 12 numerical tests, the ordered-loop smoke and native calendar/budget/purchase/fiscal/value/probability checks pass. No occupied wealth is clipped. The two receipt cases have the same expected transfer at each age and the same fixed house price. Housing and estate funding residuals below are deliberately reported.

All economic changes are experimental: age-pooled published receipt profiles, zero receipts outside supported model nodes26–78, IID receipt risk, and start-period liquid-wealth timing. The certain and lottery cases share this timing. Entry wealth, B15 earnings, preferences, the inherited 8.751% tax and the frozen target/weight contract remain fixed. The target table includes the inherited experimental first-birth rooms target1.465; this does not adopt it as the future empirical/model observation contract.

| Moment | Target | No receipts | Certain mean | Lottery |
|---|---:|---:|---:|---:|
| initial_normalization | 2.100 | 2.100 | 2.196 | 2.139 |
| cps_childlessness | 0.198 | 0.183 | 0.151 | 0.170 |
| cps_exactly_one | 0.214 | 0.231 | 0.225 | 0.229 |
| nchs_mean_age | 25.976 | 26.345 | 26.302 | 26.361 |
| nchs_share30 | 0.249 | 0.247 | 0.245 | 0.248 |
| wealth_earnings | 6.927 | 6.077 | 5.986 | 6.275 |
| bequest_wealth | 0.007 | 0.007 | 0.007 | 0.007 |
| old_dispersion | 3.516 | 2.985 | 2.879 | 2.808 |
| mean_rooms | 5.608 | 6.547 | 6.718 | 6.701 |
| ownership_30_55 | 0.676 | 0.757 | 0.737 | 0.759 |
| first_birth_rooms | 1.465 | 0.791 | 0.736 | 0.760 |
| family_rooms | 0.385 | 0.221 | 0.205 | 0.212 |
| recent_parent_ownership | 0.128 | 0.090 | 0.036 | 0.074 |

| Case | Inherited loss | House-market residual | Paid − generated estates |
|---|---:|---:|---:|
| no_receipt | 281.190 | 0.001% | -0.114690 |
| conditional_mean | 584.392 | 2.846% | -0.001528 |
| receipt_lottery | 388.612 | 2.581% | -0.006904 |

[Full target fits: every gap, weight and loss contribution](full_target_fit.csv)

[Every parameter, estimate, bound, restriction and near-bound flag](parameters.csv)

[Estate accounts and run diagnostics](accounts.csv)

The 17 standard figures for each case are retained unchanged:

- [no_receipt, figures1–9](no_receipt_standard_contact_1.png)
- [no_receipt, figures10–17](no_receipt_standard_contact_2.png)

- [conditional_mean, figures1–9](conditional_mean_standard_contact_1.png)
- [conditional_mean, figures10–17](conditional_mean_standard_contact_2.png)

- [receipt_lottery, figures1–9](receipt_lottery_standard_contact_1.png)
- [receipt_lottery, figures10–17](receipt_lottery_standard_contact_2.png)
