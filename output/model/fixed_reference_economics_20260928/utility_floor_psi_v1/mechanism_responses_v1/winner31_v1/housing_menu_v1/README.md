# Owned-room menu experiment, frozen winner31

Reference: verified chain7 case0173, loss31.2840072557; prescribed q=.719168368828958; Stone–Geary parent housing floor2.3;120 wealth×9 income; financed share.8; zero unsecured credit. Existing owned menu[2,4,6,8,10], continuous renter option unchanged. Saved baseline from[completed mortgage baseline](../purchase_ltv_v1/local_run/retry5/results/baseline_80_80/target_fit.csv), not rerun.

Two authorized experiments change only P.H_own and P.n_house, with shared/context arrays recomputed: add3=[2,3,4,6,8,10], nt7; full literal additions=[1,2,3,4,5,6,7,8,10], nt10. All31 effective scalar parameters, original14 target rows/weights, entry distributions, housing units and tenure-choice shock scale remain unchanged. Adding one-room support and additional logit alternatives is an economic experiment: the unchanged unnormalized logsumexp admits variety value. No implicit utility or weight normalization. No CES work.

[Driver](run_housing_menu.py) uses the existing authenticated fixed-price runtime and unchanged gates. Fixed-PRE impact replay is omitted because the state axes change; these are stationary cohort outcomes. Prices do not clear markets and birth renewal is reported rather than imposed. No GE, recalibration, adoption, or pure-grid causal claim. All two lifecycle calls completed in38.22 seconds total (7.27/10.86 seconds lifecycle time), one core/thread, within120 seconds percell and300 seconds overall. Original0LC missing-optional-attribute initialization receipt retained under results/; resolved using existing engine false-default convention, with no model edit.

| Outcome | Saved baseline | Add3 | Add1,3,5,7 |
|---|---:|---:|---:|
| Completed fertility |2.100000|2.094746|2.091195|
| Childlessness |.202533|.204053|.204894|
| Early fertility |.531218|.528612|.527699|
| First-birth rooms response |1.223193|1.353624|1.293638|
| Mean rooms |5.977125|5.918940|5.792074|
| Ownership30–55 |.656965|.673527|.696148|
| Recent-parent ownership moment |.118374|.076510|.086665|
| Mean first-birth age |25.95638|25.97684|25.97894|
| Original weighted loss |31.2840|89.8826|64.2750|

The added3-room option is used by8.098% of all households in add3, but its current-parent mass is only2.5225e-6 of all households (about.000669% of current parents). Among ages18–42, its parent mass is7.8617e-7. Households without current children include empty nesters and are not necessarily ever-childless. Thus this fixed-price experiment does not support the hypothesis that3 rooms meaningfully relaxes parents' discrete housing jump; fertility instead falls slightly. A3-room house still leaves only.7 rooms above the physical floor for parents. This is a source-grounded interpretation of realized choices, not an isolated welfare or lifetime value-gap decomposition. The recent-parent ownership moment is realized ownership after location/tenure transactions among current births from empty dependent homes minus ownership in current empty homes, ages30–55 with uniform within-cell age projection. Current empty homes include former parents; this selected-birth contrast is not a forced-birth causal response. The observer metadata retains its unresolved diagnostic measurement caveats.

The recent-parent ownership model moment falls to.076510/.086665 against target.12760836, contributing70.645/45.355 weighted loss units, which worsens its fit despite the higher aggregate ownership.

First-birth-conditioned tenure masses and child/no-child choice value gaps were not extracted. Supplemental age diagnostics use model age-cell starts; target clocks remain unchanged.

Complete fit tables, with target/model/gap/weight/loss contribution: [add3](results_run/add3/target_fit.csv), [full](results_run/add1_3_5_7/target_fit.csv), [all cases including saved baseline](results_run/complete_target_comparison.csv). All31 native parameter tables: [saved baseline](../purchase_ltv_v1/local_run/retry5/results/baseline_80_80/parameters.csv), [add3](results_run/add3/parameters.csv), [full](results_run/add1_3_5_7/parameters.csv). New experimental menu/support metadata is in each closure, beyond the31 scalar parameter rows. Added-size use: [add3](results_run/add3/menu_use.json), [full](results_run/add1_3_5_7/menu_use.json). Standard17-panel diagnostics remain in each case's standard_diagnostics/. [Initializer](results_run/initializer.json), [verification](results_run/verification.json), [completion](results_run/completed.json) retain source/target pins, axes, gates and counts. No further solve authorized by this packet.
