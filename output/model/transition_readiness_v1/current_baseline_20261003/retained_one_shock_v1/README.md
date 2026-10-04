# Retained one-shock experimental reference

The one-birth Estate-A numerical fit completed on October 4, 2026. The permanent 2007 preference level is **0.11999694638724082**, compared with baseline **0.17892072066041628**: a **32.93% reduction**. Original bounds are **[0.001789207206604163, 0.35784144132083257]**; the estimate is not near a bound. The final-window gap is within the original ±0.005 tolerance.

| Birth window | Target | Model | Model − target | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| 2008–2011 | 1.974875 | 1.6266185669112152 | −0.34825643308878473 | 0 | 0 |
| 2012–2015 | 1.861000 | 1.6413074510640842 | −0.21969254893591583 | 0 | 0 |
| 2016–2019 | 1.755375 | 1.6416996993988433 | −0.1136753006011566 | 0 | 0 |
| 2020–2023 | 1.645750 | 1.6431340251170146 | −0.0026159748829854834 | 1 | 0.0000068433245884109135 |

The retained fertility statistic is the household-rate analogue used in the original target contract. It is not a newly redefined female-exposure measure.

- [Numerical completion](complete.json) and [scalar fit](fit.json): selected candidate 4 freshly reproduces the fitted result; both 24/32 roots and original accounting/replay/historical-horizon gates pass.
- [Exact 2023 state](state_2023/actual_2023.pkl.gz), with [checkpoint receipt](state_2023/checkpoint_receipt.json): actual inherited distribution, both entry queues, price/pension/preference forecast and continuation values; no reconstruction or rescaling.
- Standard diagnostics: 17 unchanged plot names in each of [date 0](diagnostics/date_000/standard_diagnostics/), [date 16](diagnostics/date_016/standard_diagnostics/) and [date 31](diagnostics/date_031/standard_diagnostics/).
- Slide-layout fertility path: [PNG](slide_plot/fertility_2007_2063.png), [PDF](slide_plot/fertility_2007_2063.pdf), [plotted values](slide_plot/fertility_2007_2063_plotted_values.csv). It was drawn from candidate 1, whose four historical measurements exactly equal the final reproduction; its source pins remain unchanged.
- [All 31 baseline parameters](baseline_parameters.csv), [estimated parameter and bounds](shock_parameters.csv), [full empirical provenance and fit rows](fertility_fit.csv).
- [Retention audit](retention_receipt.json) hashes all copied artifacts and records their original paths. All 87 original source hashes are preserved: 28 copied from the runtime snapshot (including both files subsequently edited elsewhere), plus 59 still-matching working files. No model call was made during retention.

**Collection limitation:** after the numerical fit, exact state export and diagnostics completed, panel finalization detected concurrent edits to `estate_contract.py` and `observer_adapters.py`. The [failed candidate receipt](candidate_receipt.json) and [nonzero process exit](manager_terminal.json) are preserved. This package is a separately audited experimental reference, not an accepted panel receipt or a claim that the later source is equivalent.

**Scientific limitations:** both terminal checks are false. Historical four-window stability passes, but the maximum difference between the 24/32-date projections through 2063 is 0.003201; the long projection is diagnostic. Full 104/128-period certification and estate funding/recipient closure remain outstanding. Experiments must preserve these disclosures and validate their own continuation. The saved `state_experiment_ready` flag remains false.

The author requested this reference remain available while the [separate two-shock experiment](../two_shock_v1/README.md) fits the midpoint and final fertility windows. Do not overwrite this package or silently promote its certification flags.
