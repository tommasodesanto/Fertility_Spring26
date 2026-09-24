# Payroll tax for the maintained PAYGO pension system — September 24, 2026

**Proposal, not an adopted parameter change:** set the flat pension payroll rate to **8.75%** of uncapped gross working earnings. This is an aggregate OASI contribution-rate proxy, obtained by applying the 2007 combined employee–employer OASI rate to the fraction of covered earnings subject to Social Security tax. It is appropriate when the model's earnings endowment is treated as a proxy for covered cash labor earnings, including self-employment. The implementation audit confirms an uncapped household earnings object; it does not establish that the empirical earnings measure equals covered cash earnings.

PAYGO means current workers' contributions finance current pension payments. Maintain the existing rule: the pension payment adjusts to exhaust payroll receipts. No replacement-rate parameter, separate pension target, fiscal residual, model run, or parameter edit is proposed here. The rate recommendation is independent of the DUE interpretation and supersedes no author choice.

## Implemented base (independent worker receipt, lead checked)

The active parameter is currently 0.179. Working income is four-year income, \(4(1-\tau_{\mathrm{pay}})\widehat w_i a_j z\), where \(\widehat w_i a_j z\) is gross modeled annual household earnings. The inspected formulas have no individual-spouse earnings bases, statutory ceiling, or employer-cost gross-up. Payroll revenue uses the full working earnings base; pension and asset income are not taxed by this levy. A common pension is payroll revenue divided by the implemented survival-weighted retiree exposure. Sources: `code/model/intergen_eqscale_seq_optimized/parameters.py:818–859`; `code/model/tools/e5f_social_security.py:78–115`; `code/model/tools/e5f_stationary_paygo.py:51–65`. The rate itself is dimensionless: do not multiply 8.75% by four; the code already multiplies the earnings flow by the period length. These findings support an **aggregate flat approximation**, with the incidence and coverage limitations below, rather than an exact household statutory calculation.

## Official 2007 rate and flat-rate calculation

OASI is Old-Age and Survivors Insurance; OASDI adds Disability Insurance; HI is Medicare Hospital Insurance. In 2007 the employee and employer each paid OASI 5.3%, DI 0.9%, and HI 1.45%. Thus combined OASI was 10.6%, OASDI 12.4%, and HI 2.9%. The OASDI cap was **$97,500 per worker**, whereas HI had no earnings cap. Employee-only OASI was 5.3%; employee-only OASDI was 6.2%. Sources: [SSA, 2008 Supplement, Table 2.A3](https://www.ssa.gov/policy/docs/statcomps/supplement/2008/2a1-2a7.html#table2.a3), [SSA, 2007 Fast Facts, general information](https://www.ssa.gov/policy/docs/chartbooks/fast_facts/2007/fast_facts07.html).

Using the revised historical **2007** row in [SSA's 2025 Supplement, Tables 4.B1–4.B2](https://www.ssa.gov/policy/docs/statcomps/supplement/2025/4b.html#table4.b1), in millions of current dollars:

| Matched earnings category | Uncapped covered earnings | Taxable earnings |
|---|---:|---:|
| Wage and salary | 5,900,235 | 4,973,300 |
| Self-employment | 481,071 | 294,900 |
| Total | 6,381,306 | 5,268,200 |

Both totals include wages and self-employment; this does not divide a wages-plus-self-employment numerator by wages alone. The source's self-employment denominator is reported net earnings. The taxable fraction is 82.557%. Accordingly,

\[
\widehat\tau_{\mathrm{OASI}}
=0.106\frac{5{,}268{,}200}{6{,}381{,}306}
=0.08751017424959717.
\]

This is a statutory-rate-times-effective-base calculation, **not a direct estimate of cash receipts divided by income**. SSA's contribution table is based on reported earnings and does not make the flat rate an exact replica of every individual's liability. All historical amounts above use the same 2025 statistical vintage; the year being calibrated remains 2007. The corresponding OASDI calculation is 10.237%; the wage-only OASI calculation is 8.935%. These are definition checks, not additional recommended parameters.

The **alternative is 10.6%**, if the author deliberately wants the conventional flat statutory OASI abstraction used in some lifecycle models. Applied to uncapped earnings, it raises roughly 21% more revenue than the cap-adjusted proxy at the same aggregate earnings. It should not be described as the observed effective tax on uncapped household earnings.

## Why OASI, and what the incidence assumption means

For a retirement pension, OASI is the closest conventional U.S. program benchmark. Including DI would direct a disability-program levy to old-age pensions; including HI would additionally direct health-insurance revenue to pensions. OASI itself includes survivor benefits, including some paid before retirement. Its use here is therefore a retirement-program proxy, not an exact retired-worker-only institutional match. Conesa–Krueger explicitly use OASI for this reason, excluding disability and Medicare.

The recommendation attributes **both statutory sides of the OASI burden to households**. If the model income endowment is cash wages before employee taxes, subtracting the combined rate is a maintained incidence approximation; it does not reproduce actual paycheck withholding. For example, below the cap, a $100 cash wage implies $5.30 employee OASI and $5.30 employer OASI: the worker retains $94.70 before other taxes, while pension contributions total $10.60. A model deducting $10.60 from $100 leaves $89.40. This cannot be called an exact cash-budget mapping. If model income instead denotes labor cost inclusive of only the employer OASI payment, the consistent below-cap combined wedge is $10.60/$105.30 = 10.066%. With other employer benefits included, use that broader compensation denominator. These accounting translations are not alternative fiscal structures.

A national uncapped earnings base also includes some noncovered employment. The proposed 8.75% preserves an explicit covered-earnings proxy rather than claiming a fully national effective rate. It does not mechanically apply the individual $97,500 cap to household earnings: separate earners have separate caps. The aggregate approach avoids inventing a household cap where the model lacks individual earners.

Self-employment needs its own caution: in the ordinary 2007 calculation, Schedule SE multiplies business profit by 92.35% before applying the self-employment tax and shares the OASDI cap with wage earnings. Thus total business income is not automatically the statutory net-earnings base. [IRS, 2007 Schedule SE, p. 1 line 4 and p. 2 lines 4a, 7–10](https://www.irs.gov/pub/irs-prior/f1040sse--2007.pdf).

## Relevant lifecycle literature

All three papers below explicitly finance public pensions with payroll taxes. Their numerical rates are not interchangeable calibration targets.

| Paper and economy | Pension tax/base | Pension budget closure | Relevance and primary location |
|---|---|---|---|
| Conesa and Krueger (1999), United States | Flat tax on labor earnings, no distinct employer payment in the household budget. Historical OASI benchmark **10.7%**; excludes Medicare and DI. | PAYGO. They choose the replacement rate to match the observed **10.7% OASI tax**; given their demographic dependency ratio, this implies **50%** replacement. Revenues finance pensions. | Closest simple proportional-tax, uniform-pension precedent. Its 10.7% is historical, not the 2007 statutory split. [Published article, pp. 759–762 and 765–766](https://editorialexpress.com/jrust/econ698s/conesa.pdf#page=9). |
| Fuster, Ayşe İmrohoroğlu and Selahattin İmrohoroğlu (2007), United States | Separate flat social-security tax on earnings. Baseline payroll **10.3%**, distinct from **17%** general labor-income tax. Household equation combines the wedges; no explicit employee/employer division. | Payroll rate adjusts each period to finance a progressive pension formula; benchmark average-earner replacement **44%**. | Useful dynastic/heterogeneous-household PAYGO precedent. The 10.3% is endogenous model output, not an external U.S. rate to import. [Author-linked published PDF, pp. 116–118, 120, 123 and Table 2 p. 125](https://drive.google.com/file/d/1vbwO9zf5CA0NM35VKVQwStIptjrw3tvU/view); [author source page](https://sites.google.com/usc.edu/ayse-imrohoroglu/research). |
| Huggett and Ventura (1999), United States | Separate proportional Social Security tax on labor earnings; equation (5) does not cap the tax, while benefit-crediting earnings are capped. No separate employer wedge. | Explicit PAYGO balance; baseline pension schedule determines required payroll rate. | Confirms the distinction between pension tax and general income tax. Not a pure OASI numerical reference: benchmark transfers also represent hospital/medical benefits. [Published article, pp. 504–508 and 512–514](https://editorialexpress.com/jrust/econ698s/hugget.pdf#page=7). |

A useful additional base check is Huggett–Parra's **September 2008 manuscript**, later published in JPE (2010): it explicitly uses the combined **10.6% OASI** rate and caps taxable earnings at 2.42 times average earnings (pp. 9 and 14). Its partial-equilibrium tax-transfer setup is not used here as evidence of a separate annual PAYGO budget closure. [Primary manuscript](https://economics.ucr.edu/wp-content/uploads/2019/11/Huggett.pdf#page=9). The earlier Sommer source checks remain in `lead_review.json`; they do not identify a financed-pension tax rate.

## Interpretation of the maintained closure

The recommended scalar controls the contribution pool; the existing PAYGO equation determines pensions from that pool and the worker/retiree distribution. It does not separately promise a replacement rate. In actual 2007 OASI, net contributions were **$560.9 billion** and benefit payments **$489.1 billion**; other income and expenditures also existed. Therefore a pure PAYGO model calibrated to the contribution burden will not simultaneously reproduce that year's benefit outlays without additional calibration changes, which are outside this recommendation. [SSA, 2008 Trustees Report, Table II.B1](https://www.ssa.gov/oact/TR/TR08/II_cyoper.html).

Verification: source text and table definitions checked; arithmetic recomputed using the displayed integers. No empirical reconstruction, model solve, GE run, or code change. Recommendation awaits the author's choice; the independent income-base receipt was incorporated and lead checked.
