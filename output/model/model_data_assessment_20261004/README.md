# Model–data plots and older-age wealth check

Verified case: `20261003T175652812716Z_b1c72f13`, post-interest soft chain 13. The existing `code/model/plot_model_aggregates.py` writes the three-page comparison under that case's `aggregate_plots/model_data_assessment/`. No model, input, target or calibration change was made.

The fertility page also compares mean children ever born (capped at three) across five income groups, pooling ages 18–45. CPS 2024 uses family money income; the model uses current labor earnings. Tied income categories split fractionally across equal-weight groups. The model uses the observed CPS pooled age weights; the conditional age composition within each income group can differ.

The [older-age extraction](oldwealth_model_check/extract.py) reads the pinned case and PSID cache; its [summary](oldwealth_model_check/summary.json) records hashes, samples and quantiles. Reproduce with the pinned model Python, one thread. Current target and parameter tables remain in the pinned case.

The reported older-age wealth-to-income dispersion, 2.904, reproduces exactly using the authenticated observer's age-overlap weights. It is a zero-weight validation statistic against 3.516 in PSID 2003/2005. The model uses pension income while the data use total family income; this is an unresolved comparison limitation. The native whole-cell statistic, 2.740, is a different age operator.

The separate wealth-only diagnostic for ages 76–84 gives median/90th-percentile wealth of 3.005/8.726 in the model and 3.528/19.473 in PSID 2005/2007, each in its own mean working-age earnings units. This is not a proposed replacement target. The plotted mean is sensitive to the empirical upper tail; no arithmetic error was found. Model late-life financial positions fall while gross housing value remains broadly stable. The causal role of saving, borrowing and bequest rules remains untested.

Empirical older-age target provenance: `output/model/e5f_matched_pf_20260909a/design_research/wealth/build_initial_wealth.R` and `old_wealth_results.csv`. Current observer receipt: pinned case `native/phase_b_ge/selected_root/observers.json`; it supersedes historical age-mask notes in the September observer contract.
