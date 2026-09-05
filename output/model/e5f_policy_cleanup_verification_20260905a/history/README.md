# E5F transition-calibration pilot

This packet measures the existing twelve-row E5 target system on the
simulated 2023 distribution. Date 0 preserves conditional states from
an old steady-state benchmark but reweights age masses to the observed
2007 householder-age profile. In panel mode the active structural parameters
and the terminal child-preference change are estimated jointly on the
dated 2023 target system. This remains a search packet, not a
promoted estimate.

The fertility stock targets use the model cohort centered at age 42.
First-birth timing combines old-steady-state prehistory with dated
first-birth hazards during 2007--2023 and labels four-year age cells
by their midpoints. The previous period-TFR
normalization is therefore not treated as the completed-fertility
calibration target.

The 2007--2023 household-count and age paths are matched by an explicit
Census/ACS bridge. This is an externally normalized migration and
household-formation input: the code does not attribute that growth to
post-2007 births, which cannot mature before the terminal date.

Primary files: `summary.json`, `candidate_summary.csv`,
`target_fit_long.csv`, `parameter_table.csv`, and the two diagnostic
figures. Each case folder contains the full transition and cohort timing
ledger.
