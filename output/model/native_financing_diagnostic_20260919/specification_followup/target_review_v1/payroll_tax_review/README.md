# Pension payroll-tax source reviews

- **September 25 focused DUE review:** `fable_due_rate_review_20260925/` contains the completed fresh Fable 5.1 review, prior-version source evidence, and independently retrieved OECD data. DUE's earlier 17.9% describes personal income tax for 2012–2016; its newer 15.6% prose and unchanged table require a version/denominator qualification. The OASI proxy of 8.751% measures a different fiscal object, not the total tax on working income. See the packet's lead review for corrections to overstrong Fable interpretations. No fiscal specification or new run was adopted.

- **September 25 bounded Fable review:** `fable_options_review_20260925/lead_review.json` records the author-shortened fifteen-minute attempt. Fable 5.1 collected sources but returned no final opinion before its research cutoff; one session stopped with exit143, without retry. No specification was adopted. The lead overview separates baseline generosity from the transition financing rule; a closure comparison should preserve the same initial steady state.

- **September 25 literature comparison:** `olg_paygo_source_comparison_20260925.md` compares four U.S. OLG papers, their payroll bases, numerical rates, benefit calibration, and PAYGO closure. It distinguishes externally targeted rates from endogenous financing rates and confirms that pension income is outside the payroll base. This review does not adopt a tax rate or modify a run.

- **September 24 commute run:** `commute_calibration_20260924/` retains the
  cluster-only driver, Slurm chain, frozen experimental target contract,
  preflight receipt, job IDs, and shared one-hour deadline. The completed results and reviewed report are indexed there; the lower PAYGO
  rate remains an experiment, not an adopted parameter.
- **September 24, 2026:** `paygo_rate_review_20260924.md` and its JSON receipt propose a 2007 OASI-based flat payroll rate of 8.75%, pending author choice. PAYGO remains maintained, pensions remain endogenous, and no numerical parameter or source code was changed. The report distinguishes the combined employee–employer incidence convention, taxable earnings cap, coverage, and current income-base mapping.
- Earlier source checks are retained in `lead_review.json` and `sommer_sources.json`; September 24 is the current PAYGO-specific recommendation, not a decision to adopt earlier fiscal-structure suggestions.
