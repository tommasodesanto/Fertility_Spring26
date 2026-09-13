# Successive surprise history: recovered final-window candidate

Actual household history is carried through the three fitted windows into 2019. Each shock is unexpected, and households forecast its preference level to persist. The final candidate below passes all finite-path housing, pension, household and replay gates; its fertility error exceeds the 0.005 acceptance tolerance. Search was stopped by the author on September 12 at22:48EDT. The terminal-distance check fails, and no horizon certificate or production transition-policy result is claimed.

| Window | Target | Model | Model − target | Preference | Fertility fit accepted |
|---|---:|---:|---:|---:|---|
| 2007–2011 | 1.974875 | 1.974856 | -0.000019 | 0.147087157 | Yes |
| 2011–2015 | 1.861000 | 1.861069 | +0.000069 | 0.136087157 | Yes |
| 2015–2019 | 1.755375 | 1.755536 | +0.000161 | 0.125555932 | Yes |
| 2019–2023 | 1.645750 | 1.633313 | -0.012437 | 0.111555932 | No |

These are scalar shock-fitting roots, not a new weighted structural SMM calibration: the structural parameters and their original target/weight contract remain unchanged. Every historical window and its preference proposal is shown. The final-window relative miss is -0.7557%.

Reproduction: root_receipt.json records exactly zero market/fiscal/reproduction differences; source/expected_transition.csv and source/fertility.json are the saved native aggregate path and fertility measurements. The prior three accepted windows are in source/realized_fit.json. Six four-year forecast dates are used per surprise.

Cluster job17559194: candidate_path_20260911a/batches/finite_sequences_20260912/recover_saved_2019. The author stopped this and the other no-rebate runs at22:48EDT. No further trial is running; the intended baseline must rebate property-tax revenue. Standard policy diagnostics are retained in each admissible candidate’s accepted_graphs folder.


The continuation now uses the saved 2019 surprise forecast, with preference constant after its final historical change. Model period fertility is1.6262,1.6106,1.5997,1.5993,1.6287 in windows ending2027,2031,2035,2039,2043. These values have no long-horizon certificate; the final upturn must not be interpreted as a robust prediction.

## 2023 Data/Model readout

One fixed-coordinate cluster replay, job17575193, completed in248.5seconds and reproduced the entire saved path with maximum absolute difference0. No search was run. `source/readout_2023/verification.json` and `measurement_verification.json` record exact replay and successful observers. The inherited2019state is pinned by SHA256 `7787c71ae5278ef41bade3dd797b92f07db782b29ffe502277ed3d819c627a5f`; the approved initial checkpoint remains `120ffc45c0fb8756f4182f999c96b7c0236adf315cb938190ec31cd2068c87c2`.

`figures/validation_2023.csv`, `.tex`, `.pdf` and `.png` contain all13 initial-calibration moment families. Rows are untargeted comparisons at the carried2023model state. Empirical sources are CPS2024 (nearest fertility supplement), NCHS2023, ACS2023, explicitly pooled PSID samples, and the external bequest benchmark. The CSV preserves exact definitions, gaps and vintages. It does not silently label pooled wealth or event-study estimates as2023observations.

Completed fertility is children ever born at40–44, not period fertility or the initial2.1normalization. First-birth timing compares model2023–2027flows with annual2023birth counts on the same model age-cell midpoint mapping. First-birth rooms use the dated2019–2023 matched branch. National household/maternal proxies, the42-metro ACS sample, stochastic dependents as proxies for resident children, and model pension income versus PSID family income remain the same measurement limitations documented in the observer JSON; no new targets or weights were adopted.

The cluster replay driver and exact arguments are recorded in `/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/carried_readout_20260913/submission.json`. The full2023state remains there in `source/state_2023.pkl.gz`; only small readout artifacts are copied locally.

Regenerate the figures, full table and two-page PDF (no model solve):

```sh
/opt/anaconda3/bin/python -B code/model/tools/build_e5f_patch_readout.py --base output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/patch_readout_fit --sequence-base output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/recovered_sequence --pdf output/pdf/e5f_current_history_and_2023_fit.pdf
```
