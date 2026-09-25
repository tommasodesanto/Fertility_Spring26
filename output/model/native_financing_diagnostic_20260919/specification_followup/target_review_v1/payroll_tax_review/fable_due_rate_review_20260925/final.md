Both deliverables are saved. The source question is settled, so I'm stopping here.

**Where 17.9% came from.** It is the 2012–2016 average of the OECD Taxing Wages "Average income tax rate" for a single childless worker at 100% of the average wage. The January 2025 AEA copy (p. 26) and the February 2025 NBER copy (p. 28, footnote 29, OECD.Stat accessed Sept 14, 2023) state this in so many words, and the February Table 1 (p. 25) prints the same row as the December table. The December 5, 2025 version moved the calibration window to 2005–2009 and rewrote the prose to 15.6% with a new footnote, but left Table 1 at 0.179. That is a stale-table version inconsistency, not two tax concepts.

**Does OECD verify the numbers?** I pulled the series from the OECD SDMX API and re-read the raw CSV rows.
- 2012–2016 mean of the gross-wage-based rate is 17.98%, which rounds to 18.0%. The printed 17.9% is consistent with a 2023 data vintage or truncation.
- 2005–2009 mean of the same rate is 16.95%, which does not match 15.6%. Dividing each year's income tax by labour costs instead gives 15.58%, matching the printed value. So the December prose almost certainly uses the tax-wedge decomposition convention. This is arithmetic plus inference; the decomposition dataflow itself could not be retrieved after six URL attempts.

**Tax concept.** The OECD measure is personal income tax as a share of gross wage earnings. It excludes employee and employer Social Security contributions, which OECD reports separately. DUE's label "payroll tax" is a misnomer relative to its source.

**DUE fiscal closure.** Equation (2.1) taxes working-age earnings only. The text says revenue goes to pensions and wasteful spending. The equilibrium definition on p. 20 has no government budget constraint, and the pension level has no Table 1 row or calibration sentence in the pages read. DUE does not have a pension-financing closure.

**Recommendation in the file.** Neither DUE number can calibrate a pension-only levy. Using 17.9% requires adopting DUE's structure with a general income tax and a wasteful residual. The 8.751% OASI proxy is the right kind of object for the maintained PAYGO closure, with its already-documented base and incidence caveats unresolved.

**Open items.** Confirming the labour-cost "Income tax" measure in the OECD decomposition dataflow (one UI check), which value the December model was actually solved with, and DUE's pension level.

Files written:
- `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/payroll_tax_review/fable_due_rate_review_20260925/preliminary.md`
- `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/payroll_tax_review/fable_due_rate_review_20260925/final.md`
