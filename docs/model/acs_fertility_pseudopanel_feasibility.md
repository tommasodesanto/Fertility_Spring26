# ACS fertility instruments and pseudo-panels: feasibility memo

September 17, 2026. Scope: design and metadata review; no new estimation.

**Recommendation: proceed to a small measurement pilot, not a causal pseudo-panel specification.** Existing ACS files contain the necessary household links and housing outcomes. Birth histories remain incomplete, and aggregating observations does not repair instrument validity. The useful sequence is roster validation, pooled microdata estimates, then an independently evaluated cohort extension. The [sample and variable inventory](acs_fertility_pseudopanel_inventory.md) records verified inputs and remaining gaps.

## What this could identify

An instrument is an event that changes fertility while affecting the housing outcome only through the specified fertility treatment. Same-sex first two children can encourage a third child; twins at first birth can induce reaching two children sooner. These are different treatments and populations. The [Angrist–Evans study](https://www.nber.org/papers/w5778) provides a labor-supply benchmark, not evidence that either instrument is excluded from housing demand.

For housing, same-sex siblings may share bedrooms more readily. This possible direct channel was already flagged in the project's corrected PSID audit; it is a threat, not an established estimate. In a constant-effect illustration, let \(\rho\) be the instrument's housing reduced form, \(\pi\) its fertility first stage, \(\beta\) the fertility effect, and \(\gamma\) the direct housing effect:

\[
\rho=\beta\pi+\gamma,\qquad \beta_{IV}=\beta+\gamma/\pi.
\]

For example, a first stage of 0.05 and a direct effect of −0.02 rooms imply −0.40 rooms of IV bias. These are illustrative numbers. Report reduced forms and sensitivity to \(\gamma\); testing only families with exactly two children conditions on a response to the instrument and cannot establish exclusion. Ownership is not automatically immune to a direct housing-demand channel.

Twins also change spacing and simultaneous care needs. Maternal selection into twinning is documented by [Bhalotra–Clarke](https://www.iza.org/de/publications/dp/10405/the-twin-instrument). Age adjustment and measured balance are necessary diagnostics, not proof of random assignment. Neither design directly identifies housing's effect on fertility or the model's aggregate fertility-preference shocks. Even a credible local fertility-to-housing effect needs an explicit model measurement mapping before becoming a calibration moment.

## Why the pseudo-panel is a separate exercise

A pseudo-panel follows means for fixed population groups across repeated cross-sections. Begin nationally with maternal birth-cohort groups and, if precision permits, broad birthplace groups. Do not define membership by current tenure, income, education, or residence. Immigration, mortality, child departures, and changing eligibility still alter composition. Estimated birth year from integer age also needs sensitivity checks.

The key distinction is population variation versus sampling variation. Averaging sibling-sex instruments within large, stable eligible groups drives their shares toward their population probabilities. Cohort and year effects can absorb much of the remaining variation; random annual fluctuations are not a new persistent shock. Grouping by instrument status preserves contrasts, but a time-invariant instrument is absorbed by group effects, so dynamic interactions require their own identification argument. [Imbens–Wooldridge, section 6](https://cemmap.ac.uk/wp-content/legacy/resources/imbens_wooldridge/lecture_34.pdf) explains the population and rank conditions behind pseudo-panel identification.

Post-birth surveys can describe outcomes by approximate child age. They cannot recover the same mother's pre-birth housing, nor assign future sibling sex to women observed before birth. A synthetic event profile is therefore not a longitudinal event study. Full age, period, and birth-cohort effects also cannot be separately identified without restrictions.

## Bounded next-stage protocol

1. **Validate measurement in PSID first.** Using observed household rosters, reconstruct ACS-style child order and same-age twin proxies; compare with PSID fertility histories. Audit absent older children, ambiguous links, and misclassification by maternal/child age. PSID's existing same-birth-year twin proxy is not clinical truth; where exact birth information is unavailable, report agreement with that proxy only. No future roster information may enter the ACS-style reconstruction.
2. **Pilot pooled ACS microdata.** Retain all ages while building links, then select mothers, initially ages 21–35 with oldest linked child under 12; show a broader age sensitivity. This restriction mitigates departures but cannot certify complete histories. Exclude ambiguous multiple-birth proxies from same-sex samples. Analyze twins-at-first-birth and same-sex-first-two separately. Match event ages and treatments to PSID; a five-year endpoint means comparing ACS observations near the fifth event anniversary, not treating current children as future five-year fertility. Do not fabricate baseline changes from cross-sectional levels.
3. **Validate aggregation separately.** On PSID common support, compare full-panel cohort trajectories with repeated cross-section resamples using identical definitions and weights; resample original family clusters and match biennial interview timing. Start with rooms and ownership means, then assess estimator differences, sampling uncertainty, and composition sensitivity. Rare twins may make IV validation uninformative; do not treat an imprecise agreement as success.
4. **Precommit a decision gate.** Report exclusions, effective cell sizes, linkage errors, first stages, reduced forms, and weak-IV-robust intervals. Freeze a substantive equivalence tolerance before examining validation estimates. Stop the causal extension if measurement or exclusion remains unresolved; descriptive cohort trajectories may still be useful. Do not expand to metro cells or search specifications for significance.

This memo completes the authorized feasibility stage. A subsequent pilot should first time one annual metadata/roster pass, set a wall-time cap from that evidence, and write one compact diagnostic packet. No download, full-data regression, target revision, or model run was performed. Luna supplied bounded file discovery; the lead corrected its incomplete raw-data inventory through direct header/codebook checks and reviewed the identification argument.
