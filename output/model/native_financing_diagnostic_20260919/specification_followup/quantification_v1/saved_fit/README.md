# Saved income-fit tradeoffs

This deterministic audit reads the 96 saved proposals from `overnight/final_search`: 89 completed receipts and 7 rejected proposals. It runs no household or equilibrium solves. The completed receipts agree exactly on the 13 target rows (12 scored plus the separate 2.1 normalization), the nine free-parameter bounds from `plan.remote.json` (with the actual beta upper bound 0.99), and the source-fingerprint dictionary. Every scored contribution is independently recomputed as `actual_weight * gap^2`; each total reproduces the saved objective.

Reproduce with: `code/model/.venv/bin/python code/model/tools/analyze_e5f_saved_income_fit_tradeoffs.py`.

The ranked table sums losses by the explicit blocks `fertility` (four rows), `housing` (three), `ownership` (two), and `wealth` (three). Near-bound flags use a distance of at most 1% of the actual plan-bound span. The room/ownership file is a nondominated set under absolute target-unit gaps among these evaluated cases. It is an observed evaluated-search frontier only; it does not establish an attainable frontier, identification, derivatives, or causal effects.

Five factual findings from the saved search:

1. The overnight selected case 60 has total loss 353.658872914; the retained incumbent from the prior search (case 3) has loss 502.745614113.
2. The selected case ranks 1 of 89 by the unchanged weighted objective.
3. The selected case's largest block is housing, with loss 171.277; this is a descriptive decomposition, not a new score.
4. The observed room/ownership nondominated set contains 4 evaluated cases; its coordinates are raw absolute target-unit deviations.
5. Seven proposals are retained separately as rejected cases and are excluded from fit rankings and frontier construction.

Limits: this is a finite, adaptive proposal sample. It cannot establish global optimization, local identification, an attainable tradeoff curve, or causal comparative statics. The target rows retain the project's existing model-observer and empirical-sample qualifications; no target, weight, parameter bound, or objective was changed.
