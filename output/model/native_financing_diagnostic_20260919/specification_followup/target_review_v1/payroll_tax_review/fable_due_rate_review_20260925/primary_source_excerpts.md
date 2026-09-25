# Primary-source excerpts recovered from the authorized Fable WebFetch payloads

These passages are extracted from the PDF binaries returned and saved by Fable's WebFetch calls, not copied from Fable's report. WebFetch itself returned a short wrapper saying each PDF was unreadable binary; the exact source files were saved under the Claude session's `tool-results/` directory. I used local text extraction on only those fetched PDF payloads; no rendering was performed.

## NBER Working Paper 33512, February 2025

- Source: https://www.nber.org/system/files/working_papers/w33512/w33512.pdf
- WebFetch call: `toolu_01VQYztitsHetkD4XtrdgCxT`; tool-result PDF saved at `/Users/tommasodesanto/.claude/projects/-Users-tommasodesanto-Desktop-Projects-Fertility-Fertility-Spring26/dd7f94a8-c2ae-40a0-b44c-a4e119ccc498/tool-results/webfetch-1790307716258-1nmr2z.pdf` (8,809,966 bytes). The WebFetch text wrapper could not parse the binary, but the saved returned PDF is readable by text extraction.
- Extracted PDF page 29 (printed page 28): “The average tax rate for the period 2012–2016 was 17.9%.” It follows the sentence identifying the U.S. OECD-reported “payroll tax rate.”
- Footnote 29 on the same page points to OECD.Stat's variable “Average income tax rate” in the Taxing Wages section and gives access date September 14, 2023. This source definition is personal income tax, not an OASI contribution series.
- The extracted PDF page also states the adjacent sentence about the U.S. OECD rate and its 2012–16 window. This supports the period/source provenance directly from the version cited.

## AEA/ASSA 2025 program version, January 2025

- Source: https://www.aeaweb.org/conference/2025/program/paper/R5DzakN3
- WebFetch call: `toolu_01RMQoGKFqqcmDW5SZDs1b7V`; returned PDF saved at `/Users/tommasodesanto/.claude/projects/-Users-tommasodesanto-Desktop-Projects-Fertility-Fertility-Spring26/dd7f94a8-c2ae-40a0-b44c-a4e119ccc498/tool-results/webfetch-1790307867630-fs0k7u.pdf` (8.8 MB). The wrapper could not parse the binary; extracted PDF page 26 contains the same “2012–2016 was 17.9%” passage and OECD.Stat footnote 29.

## OECD primary data response

The Fable WebFetch call `toolu_01TKk1fkrdxzPD7S6EawUQJm` returned a table from the OECD SDMX endpoint with metadata: measure “Average income tax rate”; unit “Percentage of gross wage earnings”; household “Single person, no children”; earnings “100% of average wage.” The exact series API request and machine-readable row evidence are independently preserved in the lead's `lead_oecd_income_tax.csv` and `lead_arithmetic.json` in this same packet. The unfiltered OECD component response call `toolu_01Fmdb8sGAQqy6vUrphpyhiQ` also returned separate AV_ITR, employee SSC, employer SSC, average-tax-wedge rows; those do not establish the Dec. 2025 footnote's denominator because the OECD decomposition dataflow fetch failed.

## Limits

The fetched NBER and AEA PDFs establish the earlier-version sentence and footnote. They do not by themselves prove why the Dec. 2025 Table 1 remained at 0.179. “Stale table” is an inference from the cross-version comparison. The OECD arithmetic and the 15.58% labour-cost conversion do not independently establish which denominator the Dec. 2025 authors used for 15.6%.
