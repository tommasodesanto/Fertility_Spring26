# Recovered 2023 permanent-shock profile

Collected from Torch recovery job `17707689` (retry2) on 2026-09-13.  This is
the 2023 cross-section of permanent preference-shock iteration 3, evaluated at
the saved fixed-price path; it is not a converged finite equilibrium.

`verification.json` passes exact reproduction of all five saved aggregate rows
(maximum absolute gap `0.0`, tolerance `2e-10`), with 22 backward and five
forward dates and zero root solves.  `model_2023.json` SHA-256 is
`abcb26c19639b51d19c770ddd4429d09f52e44f90369870889b2e59d05be001e`.

Run `/Users/tommasodesanto/miniconda3/bin/python
code/model/tools/build_e5f_patch_readout.py --permanent-profile
PATH/model_2023.json --figures-dir PATH/figures` to rebuild the two 2023
figures and their verification sidecars.  The six large-owner allocation
shares sum to 100 percent; the age-42 completed-fertility value is
`1.8230999869363615`.
