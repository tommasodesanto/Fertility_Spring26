# Successive surprise history: recovered final-window candidate

As of September 12, 22:43 EDT. Actual household history is carried through the three fitted windows into 2019. Each shock is unexpected, and households forecast its preference level to persist. The final candidate below passes all finite-path housing, pension, household and replay gates; its fertility error exceeds the 0.005 acceptance tolerance, so search continues. The terminal-distance check fails, and no horizon certificate or production transition-policy result is claimed.

| Window | Target | Model | Model − target | Preference | Fertility fit accepted |
|---|---:|---:|---:|---:|---|
| 2007–2011 | 1.974875 | 1.974856 | -0.000019 | 0.147087157 | Yes |
| 2011–2015 | 1.861000 | 1.861069 | +0.000069 | 0.136087157 | Yes |
| 2015–2019 | 1.755375 | 1.755536 | +0.000161 | 0.125555932 | Yes |
| 2019–2023 | 1.645750 | 1.633313 | -0.012437 | 0.111555932 | No |

These are scalar shock-fitting roots, not a new weighted structural SMM calibration: the structural parameters and their original target/weight contract remain unchanged. Every historical window and its preference proposal is shown. The final-window relative miss is -0.7557%.

Reproduction: root_receipt.json records exactly zero market/fiscal/reproduction differences; source/expected_transition.csv and source/fertility.json are the saved native aggregate path and fertility measurements. The prior three accepted windows are in source/realized_fit.json. Six four-year forecast dates are used per surprise.

Cluster job17559194: candidate_path_20260911a/batches/finite_sequences_20260912/recover_saved_2019. The author stopped this and the other no-rebate runs at22:48EDT. No further trial is running; the intended baseline must rebate property-tax revenue. Standard policy diagnostics are retained in each admissible candidate’s accepted_graphs folder.


Regenerate the provisional fertility figure (no model solve):

```sh
/opt/anaconda3/bin/python -B code/model/tools/build_e5f_patch_readout.py --base output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/patch_readout_fit --sequence-base output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/recovered_sequence
```
