# Published Boar–Gorea–Midrigan unsecured bound: evidence receipt

Read-only research checked September 30, 2026. No downloaded MATLAB code was executed, no data were bulk extracted, and no model run occurred. Large PDFs/ZIP remain outside the repository.

**Established:** the published model imposes \(a'\ge\underline a\), with \(\underline a<0\). Its official baseline replication implements \(\underline a=-0.4\) model units. The authors' official README states average annual model income is 5.4614 units and one model unit is 12,896 2016 USD. Consequently the allowance is \(0.4/5.4614=0.0732413\), or about **7.32% of average annual income**, approximately **5,158 2016 USD** using their rounded dollar conversion. This is a stock bound, even though the model period is quarterly; do not multiply the bound by four.

## Published article verification

Article: [ReStud 89(3), 1120–1154](https://academic.oup.com/restud/article/89/3/1120/6372706), DOI 10.1093/restud/rdab063. The publisher HTML was directly read through the in-app browser; standard web/CURL access failed. Publisher printed page numbers were not independently retrieved, so references below are exact sections/equations, not guessed pages.

Section 3.1, **Liquid assets**, states: “We impose a borrowing limit” followed by \(\underline a<0\). Section 3.2, **Recursive formulation**, renter option \(V^1\), explicitly restricts \(a'\ge\underline a\) alongside budget equation (4). The \(V^2\) purchase-without-mortgage and \(V^3\) purchase-with-mortgage branches also invoke the borrowing constraint on liquid assets, the latter alongside LTV/PTI equations (1)–(2). It is therefore not a renter-only allowance.

Section 4.1 and Table 2 Panel B were checked directly. Table 2 reports the assigned and estimated parameters but does not report \(\underline a\). Section 4.1.2 targets the mean/median liquid asset holdings, hand-to-mouth fractions and homeowner 90th percentile; it does not state that the unsecured bound targets the 10th percentile. Table 3 Panels A/B report 10th-percentile liquid holdings as an additional model implication (all households: data −0.04, model −0.05; homeowners: −0.04/−0.04). These are distribution quantiles, not the bound.

The article's Data Availability Statement explicitly links [DOI 10.5281/zenodo.5112964](https://doi.org/10.5281/zenodo.5112964), establishing that this is its official replication package. The article was published online September 20, 2021; issue date May 2022.

## Official replication evidence

[Zenodo v2](https://zenodo.org/records/5112964), `Replication_ReStud.zip`, 83,896,919 bytes. Independently calculated MD5 is `175959450457a21b65987eedcb4b409c`, exactly the publisher record's checksum.

`Replication_ReStud/README.pdf`, printed page 1, routes baseline Tables 2–4 to `model_replication/Main nu 3/start.m`; it provides the income/unit conversion above. The retained excerpt is [replication_readme_excerpt.txt](replication_readme_excerpt.txt).

`Main nu 3/start.m` lines 94–96 fix `amin=-0.4`, `amax=100` and construct the liquid asset grid. Lines 156–157 solve the owner/renter maximum-consumption boundaries so next-period assets equal `amin`; this makes the lower grid point an enforced economic choice bound, not merely a display axis. `savings.m` lines 7/11 define owner/renter next-period assets, and line 17 compares assets to that bound. The same scalar is set once before the backward age loop, so it does not taper with age or current earnings. By contrast `lmin=-1` (start.m113) is an intermediate liquidity interpolation grid, not the unsecured asset bound. `atmin` (117) is another transformed-state grid and is not a second economic debt allowance.

`Main nu 3/objective.m` line 92 again fixes `amin=-0.4`. Its lines 16–24 map only nine estimated coordinates; `start_calibration.m` line 16 lists those nine, and none is the unsecured bound. Thus the bound is fixed during the supplied estimation routine. Code excerpts with original member paths and line numbers are retained in [replication_code_excerpts.txt](replication_code_excerpts.txt); member hashes are in [source_receipt.json](source_receipt.json).

## Publisher appendix and limitations

The actual [publisher supplementary PDF](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/restud/89/3/10.1093_restud_rdab063/1/rdab063_supplementary_data.pdf) was downloaded from the article's observed supplement link and searched. Its cover says December 2020, but publisher provenance now verifies its role as the published supplement. Its hash differs from the earlier author-hosted appendix file. No `borrowing limit` or `0.036` appears; all ten `lower bound` hits refer to the liquid interest rate. The single `borrowing` hit concerns secondary homes. The choice of −0.4 has **no verified calibration rationale in the checked published text, appendix, README or baseline parameter code**. It may be a fixed numerical/economic convention selected by the authors; do not present a rationale as established.

The April 2017 draft's 0.036 annual-income allowance and 10th-percentile rationale are historical evidence only and are not the published calibration. They differ from the verified official replication value. No manuscript or project parameter was changed or adopted by this receipt.
