# Pre-2007 fertility stocks by age: extra calibration candidates

Author-requested experimental extension, September13,2026. The aim is to improve
the inherited age profile of children ever born in the approximate2007stationary
calibration, then assess the subsequent completed-fertility path. The original
objective remains retained. Point extraction does not activate new targets,
weights, standard errors or a calibrated model.

The extraction reads only the June2004 and June2006 partitions of the existing
IPUMS CPS fixed-width source, verifies both partition hashes and the loader hash,
and keeps femaleSEX2, the stated age band, FREVER0–20 and positive FRSUPPWT.
Weights use the loader's division by10000. The pooled estimates use the same
observed survey weights as the existing early calibration source.

Files: `age_profile_candidates.csv` and `.json` contain five age groups
(20–24 through40–44), separately by year and pooled. They report weighted
0/1/2/3+ shares, childlessness, exactly-one conditional on motherhood, and three
different means: uncapped empirical children ever born, empirical children ever
born capped atfive, and observed child-count groups coded with the maintained
model3+ representative3.602359422009. These means must not be conflated. The
model-coded mean is an encoding of empirical shares, not a direct count mean.
No uncertainty estimate or SMM weight was computed.

Proposed new scored candidates are mean children ever born at25–29,30–34,35–39,
40–44 and childlessness at25–29/35–39. The mean definition and weights must be
pinned explicitly before the experimental objective is run. The existing
40–44 childlessness/exactly-one rows remain in the original objective and are
not duplicated. Full child-count shares are useful diagnostics; their adding-up
constraint must be respected. Ages20–24 are a boundary diagnostic because the
model starts households at18 without prior children.

Verification: all15 profiles have shares summing to one; model-coded means
reproduce their weighted-share formula. The2004,2006 and pooled40–44 mean,
childlessness and one-child observations reproduce the original authoritative
extraction with maximum absolute difference4.22e-15. Positive-weight FREVER999
exclusions are counted explicitly (zero in these selected partitions); such
records never enter the observations. The lead reviewed the byte fields and
sample, repaired missing verification metadata, and reran the builder.

Reproduce without a model solve or full raw-file scan:

```sh
/opt/anaconda3/bin/python -B output/model/e5f_matched_pf_20260909a/design_research/fertility_contract/age_profile/extract_age_profile.py
```

The original partition receipt is
`output/model/e5f_matched_pf_20260909a/parameter_target_audit/fertility/fertility_availability.json`.
The original working target/weight contract remains in
`output/model/e5f_matched_pf_20260909a/initial_calibration_contract/`.
The bounded experimental search and its relationship to the two demographic
versions are specified in `docs/model/e5f_two_closure_overnight_plan.md`, section2a.
