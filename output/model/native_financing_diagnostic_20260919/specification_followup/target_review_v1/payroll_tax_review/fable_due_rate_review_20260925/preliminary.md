# Provisional answer: provenance of DUE's 17.9% tax rate (preliminary, unverified beyond local evidence)

Status: PROVISIONAL. Written within the 8-minute window from local evidence only (`lead_review.json`, `paygo_rate_review_20260924.md`, `olg_paygo_source_comparison_20260925.md`, prior dispatch receipts). No new primary-source retrieval yet. Treat every claim below as a hypothesis to be confirmed or corrected in `final.md`.

## Provisional lead answer

Where 17.9% came from, as far as local evidence shows: the prior lead review recorded that the public Dec. 5, 2025 copy of *Dynamic Urban Economics* (Greaney, Parkhomenko, Van Nieuwerburgh) prints 0.179 as the tax rate in Table 1 (printed p. 27) and states 15.6% over 2005–2009 in the prose on p. 31, with footnote 33 citing OECD United States "income tax, single person without children, 100% of average wage, accessed May 15, 2025." The prior check did NOT reproduce either number from OECD. So the printed 17.9% is currently **not verified** against its stated source, and the paper's own two numbers disagree. The most likely explanation (inference, not fact): the OECD "income tax" concept for a single childless worker at 100% of the average wage is the OECD Taxing Wages personal income tax as a share of gross wage earnings, which excludes social security contributions, so DUE's number is a general income tax rate, not a payroll contribution and not a pension-financing rate. A single-year "accessed May 2025" value near 17–18% versus a 2005–09 average near 15.6% would be consistent with an average-versus-point mismatch across paper versions, but this must be checked against the OECD series before being asserted.

## Provisional comparison table

| Object | Value | Definition / tax base and denominator (provisional) | Period / vintage | Source location | What it can finance |
|---|---:|---|---|---|---|
| DUE Table 1 | 17.9% | Unknown pending PDF read; likely OECD personal income tax for single childless worker at 100% average wage, share of gross wage | Unknown (possibly a recent single year, accessed May 15, 2025) | DUE Dec. 5, 2025 copy, Table 1, printed p. 27 | In DUE: pensions plus wasteful public spending (per prior lead review of pp. 10–11) |
| DUE prose | 15.6% | Same concept, average over 2005–2009 (per prior lead review) | 2005–2009 average | DUE Dec. 5, 2025 copy, p. 31, fn. 33 | Same as above |
| OECD Taxing Wages, US, single no children, 100% AW | not yet retrieved | Personal income tax as % of gross wage earnings, or possibly the "net personal average tax rate" including employee SSC; to be verified | 2005–2009 mean and latest year to be computed | To be retrieved (OECD Data Explorer / Taxing Wages) | Not a program-specific levy; a general income-tax measure |
| Our proposal | 8.7510174% | 2007 combined employee+employer OASI statutory rate 10.6% times taxable/uncapped covered earnings 5,268,200/6,381,306 (SSA 2025 vintage, $ millions); applied to uncapped model household gross working earnings | 2007, SSA 2025 statistical vintage | `paygo_rate_review_20260924.md` | Only the common PAYGO pension in our model |

## Arithmetic (verified locally)

Our proxy: 0.106 × (5,268,200 / 6,381,306) = 0.106 × 0.825567... = 0.0875101742... (dimensionless share of gross working earnings; SSA 2025 Supplement Tables 4.B1–4.B2, revised 2007 rows, $ millions).

DUE/OECD averages: not yet reproduced.

## Provisional recommendation

- Neither DUE number should be described as a pension payroll rate until the OECD concept is verified. If the OECD measure is personal income tax (excluding social security contributions), 17.9% and 15.6% are general income-tax rates and cannot calibrate a pension-only levy.
- Using 17.9% in our model would require a general income tax plus an explicit residual for non-pension spending, matching DUE's fiscal statement that revenue finances pensions and wasteful spending.
- Unresolved: exact Table 1 label and note, the 2005–09 vs point-year discrepancy, and the OECD concept.

## Next checks (in progress)

1. Read local `tmp/pdfs/DUE_latest_checked.pdf` pages around printed pp. 10–11, 27, 31: exact row label, footnote attachment, fiscal equations.
2. Retrieve OECD Taxing Wages US series for single person, no children, 100% AW: personal income tax and net personal average tax rate; compute 2005–2009 mean and latest value.
3. Write `final.md`.
