# Income Risk, Subsistence, and Feasibility: Decision Memo

**Date:** 2026-07-17

**Applies to:** the post-M5 intergenerational housing--fertility model

**Decision:** use realistic income risk only with an empirically measured,
asset-tested, budget-balanced residual transfer. Keep its policy guarantee
separate from the Stone--Geary preference minimum and use feasibility as a
model test. Do not make low income risk, an unexplained income truncation, or
unsecured credit the paper baseline.

## Executive decision

The literature does not contain a device that lets a household permanently
consume a Stone--Geary minimum exceeding its cash resources without assigning
the shortfall to somebody. The verified choices are to (i) pay a residual
transfer, (ii) impose a support restriction that rules out the state, or
(iii) add a resource-providing margin such as labor supply or family support.
Standard housing papers largely avoid the issue because utility is homothetic,
rental housing is continuous or its minimum is calibrated to be affordable,
and consumption has no positive subsistence intercept.

The recommended M6 specification distinguishes the required bundle from the
policy guarantee:

\[
  B_j(n)=\bar c(n)+r_j\bar h(n), \qquad
  T_j(x,n)=\max\{0,G_j(n)-x\}.
\]

where \(x\) is current cash-on-hand after the ordinary tax schedule and
realizable assets and \(G_j(n)\) is a measured safety-net guarantee, not a
function of the preference parameter \(\bar c(n)\). A recipient must choose
\(a'=0\) and cannot use unsecured debt. A proportional levy \(\tau_T\) on
labor income closes

\[
  \tau_T\int y\,d\mu=\int T_j(x,n)\,d\mu .
\]

The tax is an equilibrium residual, not a calibrated parameter. The transfer
must be disciplined by statutory generosity plus receipt and outlay moments by
age and parent status. Feasibility requires \(G_j(n)\geq B_j(n)\) for states
that can reach zero other resources. If that restriction or the transfer
moments fail while realistic income risk and the existing
wealth/fertility/tenure targets are retained, the joint preference/support
specification is rejected. We should then add family co-residence/support or
lower and re-identify \(\bar c_0\), not quietly shrink income risk or set
\(G_j(n)=B_j(n)\) by construction.

## Two distinctions that change the diagnosis

First, a high benefit per recipient is not the same as mass participation. At
the verified age-22 probe, a dead mass of \(0.0016\) and an illustrative
resource gap of \(0.11\) mean annual income imply a mechanical transfer cost
of only \(0.000176\) of mean income, or 0.0176%, for that age cross-section
before behavioral responses. The solved receipt rate and outlay—not the size
of the guarantee alone—determine whether the institution is a narrow safety
net or a broad program.

Second, the Hubbard--Skinner--Zeldes (HSZ) \$7,000 precedent is not a universal
young-adult floor. The detailed construction uses a "female-headed family
with two children" for the nonelderly and an elderly program bundle for retirees
([HSZ 1995, pp. 29--31](https://www.nber.org/system/files/working_papers/w4884/w4884.pdf)).
It therefore validates the possible *size of a targeted guarantee*, not the
interpretation of \(\bar c_0\simeq 0.30\) as a universal childless-household
subsistence preference.

## 1. Floor levels in verified structural models

Income shares use the contemporaneous Census median household income when the
paper's dollar year is verified. They are descriptive scale comparisons, not
harmonized benefit-replacement rates.

| Model | Verified device and value | Share of median household income | Funding in the displayed model |
|---|---|---:|---|
| Hubbard, Skinner, and Zeldes (1994) | Eq. (5), p. 176: \(TR_s=\max[0,(\bar C+M_s)-(A_{s-1}(1+r)+E_s)]\). The same page sets \(\bar C=\$7{,}000\) in 1984 dollars ([paper](https://business.columbia.edu/sites/default/files-efs/imce-uploads/szeldes/pdfs/expanding.pdf)). | 31.2% of the 1984 median \$22,420 ([Census, Table A](https://www2.census.gov/prod2/popscan/p60-151.pdf)). | No government budget equation was found in the verified household problem, pp. 174--177; the transfer enters Eq. (5). **Reading inference: not fiscally closed.** |
| Hubbard, Skinner, and Zeldes (1995) | Eq. (6), p. 19 repeats the residual-transfer formula. Pages 29--31 construct approximately \$6,937 for a nonelderly mother with two children and round it to \$7,000; the elderly calculation is approximately \$6,893 ([working-paper version](https://www.nber.org/system/files/working_papers/w4884/w4884.pdf)). | 31.2% of the same 1984 median. | No government budget equation was found in the verified model section, pp. 17--20. **Reading inference: not fiscally closed.** |
| French (2005) | The core problem has utility in consumption and hours (Eq. (1), p. 5) and the constraint \(A_{t+1}\geq0\) (Eq. (5), p. 9); no residual-transfer or positive consumption-floor equation was found ([working-paper version](https://fraser.stlouisfed.org/files/docs/historical/frbchi/workingpapers/frbchi_workingpaper_2000-02.pdf)). | Not applicable. | Not applicable. Final-publication equivalence of these page references: **NOT VERIFIED**. |
| French and Jones (2011) | Eq. (6), p. 7: \(tr_t=\max\{0,\underline C-(A_t+Y_t+ss_t)\}\). Table 5, p. 26 sets \(\underline C=\$4{,}380\) in 1998 dollars ([paper](https://www.chicagofed.org/-/media/publications/working-papers/2001/wp2001-19-pdf.pdf)). | 11.3% of the 1998 median \$38,885 ([Census](https://www2.census.gov/library/publications/1999/demographics/p60-206.pdf)). | No government budget equation was found in the verified household problem, pp. 5--9. **Reading inference: not fiscally closed.** |
| Kaplan and Violante (2014) | Eq. (4), p. 1209 requires \(c\geq0,s\geq0\), not a positive floor. Eq. (8), p. 1210 permits unsecured liquid borrowing; p. 1220 calibrates its limit to 74% of quarterly income ([published paper](https://cpb-us-w2.wpmucdn.com/voices.uchicago.edu/dist/2/981/files/2018/07/kaplan_violante_ecmtra_2014-1sxn3k8.pdf)). | Not applicable. | Not applicable. |
| Ameriks, Briggs, Caplin, Shapiro, and Tonetti (2020) | The government-care option sets wealth to zero and consumption to \(\omega_G\) (Eq. (7), p. 5). Table 4, p. 25 estimates \(\omega_G=29.45\) thousand with standard error 36.67; p. 26 describes the estimate as imprecise ([NBER version](https://www.nber.org/system/files/working_papers/w20973/revisions/w20973.rev1.pdf)). | **NOT VERIFIED:** the relevant dollar base year was not established, so no defensible income share is reported. | No government budget equation was found around Eq. (7). **Reading inference: an external option, not fiscally closed.** |
| Scholz, Seshadri, and Khitatrakun (2006) | P. 620: \(T(e_j,a_j,n_j)=\max\{0,\bar c n_j/g(1,2)-[e_j+(1+r)a_j]\}\). The same page reports \$8,159 for 1992, in 1992 dollars ([paper](https://users.ssc.wisc.edu/~aseshadr/Publications/optimality.pdf)). | 26.5% of the 1992 median \$30,786 ([Census, Table A](https://www2.census.gov/prod2/popscan/p60-186rd.pdf)). | No government budget equation was found in the verified transfer and household-budget discussion, pp. 619--621. **Reading inference: not fiscally closed.** |
| De Nardi, French, and Jones (2010) | Eq. (10), p. 44: \(b_t=\max\{0,\underline c+m_t-[a_t+y_n(ra_t+y_t,t)]\}\); receipt implies \(c_t=\underline c\) and \(a_{t+1}=0\). P. 58 estimates \(\underline c=\$2{,}700\) in 1998 dollars ([paper](https://users.nber.org/~denardim/research/De_Nardi_French_Jones_JPE_2010.pdf)). | 6.9% of the 1998 median. | No government budget equation was found in the verified model section. **Reading inference: not fiscally closed.** |

The verified floors span roughly 7%--31% of median household income. The top
of that range is real, but it is attached to program eligibility and family
composition. None of the verified papers establishes a universal
Stone--Geary *preference* intercept near 30% for young childless households.

## 2. Stone--Geary or child minima with persistent income risk

The closest fertility model does not solve our problem with a transfer.
Sommer (2016) imposes \(q\geq\bar q\) only for parents in Eqs. (2)--(3), p. 31,
and calibrates \(\rho=0.95\), persistent-shock standard deviation 0.21, and
transitory-shock standard deviation 0.17 in Table 2, p. 31. Footnote 12, p. 32
says the case would require "an exogenous transfer from an unmodeled
government" and reports that "no such cases are detected"
([paper](https://www.kamilasommer.net/Fertility.pdf)). Thus feasibility is an
ex post numerical property of her calibration, not a funded institution; the
minimum is child quality, not childless nonhousing consumption.

A direct Stone--Geary portfolio model instead restricts the support. Achury,
Hubar, and Koulovatianos impose \(k_0>\chi/r_f\) in Assumption 1, p. 7, which
keeps initial resources above the subsistence annuity; pp. 7--10 contain no
labor-income process ([paper](https://www.nottingham.ac.uk/cfcm/documents/papers/10-01.pdf)).
This is a mathematical feasibility restriction, not an empirical income-risk
solution, and the paper supplies no calibrated \(\chi\)-to-income ratio.

Doepke and Kindermann's fertility environment is not a positive-consumption
floor precedent: their displayed instantaneous payoff is linear in
consumption, \(u(c_g,d_g,v_g,b)=c_g-d_g+v_gb\) (Eq. (7), p. 25)
([paper](https://mdoepke.github.io/research/Doepke_Kindermann_1118.pdf)).
Verified quantitative Barro--Becker and Daruich implementations that combine a
universal Stone--Geary floor of our magnitude with persistent risk were not
located: **NOT VERIFIED**.

**Conclusion.** In the verified set, a positive minimum is handled either by
an HSZ residual transfer, an ex ante support restriction, or a calibration
for which the constraint does not bind in simulation. There is no verified
paper that supplies a fourth accounting solution.

## 3. Minimum housing with income risk

| Model | Verified poorest-household treatment | What it implies for us |
|---|---|---|
| Sommer and Sullivan (2018) | Eq. (20), p. 253 is CRRA over a Cobb--Douglas composite of \(c\) and housing services, with no Stone--Geary intercept. P. 246 says renters may choose a unit "smaller than the minimum house size that is available for purchase." Table 1, p. 254 sets annual \(\rho=0.90\) and innovation standard deviation 0.20 ([paper](https://kamilasommer.net/Taxes.pdf)). | A low rental product is enough only because consumption has no positive intercept. The paper does not provide a transfer, default, co-residence, or homelessness state in the verified model pages. |
| Kaplan, Mitman, and Violante (2020) | Eqs. (1)--(2), p. 7 require positive \(c\) and housing services but no positive intercept; renters cannot borrow. P. 18 sets annual earnings \(\rho=0.97\), innovation standard deviation 0.20, and says the minimum rental size is calibrated; Table 1, p. 20 reports rental sizes \(\{1.17,1.50,1.92\}\) ([paper](https://www.nber.org/system/files/working_papers/w23694/w23694.pdf)). | The verified pages report no separate feasibility program. Positive endowments plus a calibrated rental menu avoid an explicit subsistence gap. |
| Boar, Gorea, and Midrigan (2020) | The household problem on pp. 12--14 uses homothetic consumption--housing utility, continuous rental services, liquid borrowing, and home production; no positive nonhousing intercept appears. Table 2 reports annualized \(\rho=0.964\), persistent-shock standard deviation 0.150, and transitory standard deviation 0.327 ([paper](https://www.nber.org/system/files/working_papers/w23345/w23345.pdf)). | Their feasibility is produced jointly by flexible housing, credit, home production, and no Stone--Geary floor. Importing only the shock process is not a like-for-like exercise. |
| Greaney, Parkhomenko, and Van Nieuwerburgh (2025) | Eq. (8), p. 9 is Cobb--Douglas in consumption and housing. P. 30 makes the rental set continuous with lower bound zero ([current author version](https://www.andrii-parkhomenko.com/files/Dynamic_Urban_Economics.pdf)). | Letting rental services approach zero removes a discrete housing floor; it does not fund positive \(\bar c_0\). Publication status and final pagination: **NOT VERIFIED**. |
| Favilukis, Ludvigson, and Van Nieuwerburgh | The model's period utility is CES over nonhousing consumption and housing services and housing is continuous (Eqs. (1)--(3), pp. 7--8); the paper states there is no explicit rental market on p. 7 ([NYU archive version](https://archive.nyu.edu/handle/2451/29876)). | It is not a discrete-renter-floor precedent. Exact final-publication pagination: **NOT VERIFIED**. |

No verified paper in this housing set contains doubling up, moving in with
family, or homelessness as the bottom housing state. Kaplan (2012) is the
direct structural precedent for adding family support: p. 447 gives parents
"both monetary support, through explicit financial transfers" and support
through shared residence, with a "reduction in per capita direct housing costs"
([paper](https://gregkaplan.me/s/kaplan_jpe_2012.pdf)). Co-residence alone,
however, does not close our nonhousing subsistence gap; parental resources or
an in-kind consumption flow must also enter the budget.

## 4. What the imported income processes measure

| Source | Verified income concept and process | Importability |
|---|---|---|
| Sommer and Sullivan (2018) | P. 254 says evidence from Card, HSZ, and Heathcote--Storesletten--Violante places \(\rho\) in 0.88--0.96 and innovation standard deviation in 0.12--0.25; the authors then set 0.90 and 0.20. Their model state is labor productivity before the progressive tax schedule (budget equation (7), p. 248) ([paper](https://kamilasommer.net/Taxes.pdf)). | The 0.90/0.20 pair is a literature calibration, not an estimate on a Sommer--Sullivan PSID sample. It is pre-tax productivity and does not include a transfer floor. |
| HSZ (1995) | Table 2, p. 28 reports AR coefficients 0.955/0.954/0.959 and innovation variances 0.033/0.026/0.020 by education. P. 27 defines net earnings to include unemployment insurance and subtract income taxes, while means-tested support is represented by the separate floor on pp. 29--31 ([working-paper version](https://www.nber.org/system/files/working_papers/w4884/w4884.pdf)). | Closest verified persistent-risk estimate paired with an explicit floor, but the floor's demographic construction must also be imported or rebuilt. |
| Boar, Gorea, and Midrigan (2020) | P. 18 builds taxable income from "wages (net of pension contributions), social security income, pension income, unemployment compensation and other transfers," then subtracts federal and state income taxes. Table 2 reports annualized \(\rho=0.964\), persistent standard deviation 0.150, and transitory standard deviation 0.327 ([paper](https://www.nber.org/system/files/working_papers/w23345/w23345.pdf)). | Best verified post-tax-and-transfer process in this survey. It belongs to a model with home production, flexible rent, and borrowing, so its shocks should be used as an empirical target or sensitivity—not copied alone as a structural package. |
| Blundell, Pistaferri, and Preston (2008) | P. 10 defines disposable income as labor income plus transfers, net of taxes; Data Appendix p. 46 gives the PSID construction. The process separates permanent and transitory innovations rather than reporting a single stationary AR(1) pair ([author final](https://web.stanford.edu/~pista/aer_final.pdf)). | Useful empirical moments, not a direct five-state \((\rho,\sigma)\) import. The baseline continuously married sample is described on pp. 7--8 and excludes most low-income SEO households, so bottom-tail use requires care. |

Krueger--Perri and a directly transferable Heathcote--Storesletten--Violante
post-transfer \((\rho,\sigma)\) pair were not verified in this pass:
**NOT VERIFIED**. The implementable empirical strategy is to estimate the
model's exact income concept on the project's chosen PSID sample and separately
model the safety-net schedule. Mixing a post-transfer process with an explicit
floor without matching definitions would double count insurance.

## 5. Funding verdict

For HSZ (1994, Eq. (5), p. 176), HSZ (1995, Eq. (6), p. 19), French--Jones
(2011, Eq. (6), p. 7), Scholz--Seshadri--Khitatrakun (2006, p. 620), De
Nardi--French--Jones (2010, Eq. (10), p. 44), and Ameriks et al. (2020, Eq.
(7), p. 5), the verified floor appears in the household problem but no
government budget equation was found. These are therefore classified as
unfunded/partial-equilibrium institutions **by reading inference**. Sommer
(2016, footnote 12, p. 32) explicitly describes a possible government
transfer outside the modeled equilibrium. No verified paper in this floor set
funds the guarantee with an endogenous payroll or income tax.

That absence is not a reason for us to omit funding. The present project solves
market-clearing and distributional outcomes and may make welfare statements;
the fiscal incidence is part of the mechanism.

## 6. The \(\bar c_0\) tension

There is precedent for a transfer guarantee near 30% of median household
income: HSZ's \$7,000 is 31.2%, and Scholz--Seshadri--Khitatrakun's family-size
adjusted \$8,159 is 26.5%. There is no verified precedent in this survey for a
universal, young-childless Stone--Geary preference intercept near 30% of mean
household income. The closest directly verified Stone--Geary paper enforces
the support restriction \(k_0>\chi/r_f\) (Achury et al., Assumption 1, p. 7)
and reports no comparable calibrated ratio. Thus our \(\bar c_0\) is an
outlier in interpretation, even if a benefit bundle of similar scale is not
an outlier for selected families.

This matters for identification. Once \(\bar c_0\) determines eligibility and
benefit amounts, it is no longer disciplined only by consumption/saving and
fertility behavior; program receipt and spending provide direct identifying
variation. If those implications are wrong, the response is to reject or
re-estimate \(\bar c_0\), not to conceal its resource implications with a
numerical sentinel.

## Ranked implementation menu

### 1. Empirically measured, funded HSZ residual transfer — recommended

**Mechanism.** Use the displayed \(T_j(x,n)\), enforce \(a'=0\) for recipients,
and solve \(\tau_T\) from the aggregate budget. Construct \(G_j(n)\) from the
relevant cash, food, and housing-assistance rules; do not define it using
\(\bar c(n)\). This applies the asset-test logic of De Nardi--French--Jones Eq.
(10), p. 44, while turning \(G_j(n)\geq B_j(n)\) into a refutable support
restriction. If the housing component varies with local rent, report the
induced subsidy to expensive locations and its effect on spatial sorting.

**Expected effects.** It removes all empty budget sets only when the empirical
guarantee satisfies the support restriction. The 100% phase-out and asset test
reduce saving near eligibility, may delay ownership, insure bad income draws,
and—if \(G_j(n)\) rises with children—subsidize fertility at the margin. A
location-varying housing allowance can reduce migration away from high-rent
cities. The financing levy slightly lowers resources for all earners.

**Identification.** Add at least (i) receipt rates and (ii) transfer outlays as
a share of income, each split by young childless versus parent households. The
exact SIPP/CPS/administrative target definitions, program aggregation, and
sampling uncertainty are **NOT VERIFIED and must be constructed before
calibration**. If \(G_j(n)\) is fixed from rules and the tax is an equilibrium
residual, the four split moments would move the current system from 15/14 to
19 moments/14 free parameters. They overidentify the bottom income process and
the compatibility of \(\bar c(n)\) with observed support. Do not count the tax
rate as a free SMM parameter.

**Implementation cost.** Medium: household budget/choice logic, distribution
aggregation, one fiscal residual, new empirical targets, and diagnostic plots.

### 2. Re-estimate or externally restrict the hard preference floor

**Mechanism.** If an empirically measured \(G_j(n)\) lies below \(B_j(n)\),
lower \(\bar c(n)\) and re-estimate it with direct low-resource consumption or
budget-share moments, or replace the hard Stone--Geary domain with a finite-
utility nonhomothetic specification. The latter functional form and its
literature precedent are **NOT VERIFIED in this survey** and require a separate
specification memo.

**Expected effects.** A lower or soft minimum reduces precautionary saving and
the effective resource cost of children, so fertility and ownership can move
materially. It restores a clean separation between preferences and public
policy but changes the existing identification block.

**Implementation cost.** Low for an externally restricted lower floor; medium
for re-estimation; high for a new preference specification and full
recalibration.

### 3. Re-estimate disposable income and the transfer schedule jointly

**Mechanism.** Estimate pre-program resources and safety-net receipts
separately on the project sample; use Boar--Gorea--Midrigan's disposable-income
definition (p. 18) as a comparison and HSZ's split between net earnings and a
separate floor (pp. 27--31) as the structural template.

**Expected effects.** Better bottom-tail persistence can reduce or increase
receipt incidence. A post-transfer process plus a separate floor risks double
counting. It does not itself fix the accounting problem created by a positive
\(\bar c_0\).

**Implementation cost.** Medium/high data work; low model-code cost after the
transfer institution exists.

### 4. Add parental co-residence and family transfers

**Mechanism.** Add a parental-home option with lower housing cost plus an
explicit transfer or household-public-good flow, following Kaplan (2012,
p. 447). Link the resources to the parent state so support is not free.

**Expected effects.** Strongly changes young precautionary saving, mobility,
tenure entry, and the fertility value of living near parents. It is an
economically rich treatment of the actual binding group but requires
intergenerational matching or a reduced-form parent-resource state.

**Implementation cost.** High.

### 5. Add endogenous labor supply or home production

**Mechanism.** Let low-resource households increase hours, as in French's
consumption-hours problem (Eq. (1), p. 5), or substitute home production, as
in Boar--Gorea--Midrigan's household problem (pp. 12--14).

**Expected effects.** Raises resources in some states but cannot guarantee
feasibility after sufficiently bad persistent shocks or at an hours bound. It
changes child costs, labor supply after births, and housing demand, creating a
large new parameter and target block.

**Implementation cost.** High.

### 6. Permit limited unsecured credit or delinquency

**Mechanism.** Add a liquid borrowing limit as in Kaplan--Violante Eq. (8),
p. 1210, with default/delinquency if the debt can become unpayable.

**Expected effects.** Smooths transitory shortfalls but does not fund a
persistent gap between income and subsistence. It would add young debt and
could delay ownership. The verified binding households are renters without
debt, so this is not the primary missing margin.

**Implementation cost.** Medium with exogenous limits; high with endogenous
default.

### 7. Truncate the income process

**Mechanism.** Bound the lowest income realization or restrict the state space
so cash resources always exceed \(B_j(n)\), analogous mathematically to
Achury et al.'s \(k_0>\chi/r_f\) restriction (Assumption 1, p. 7).

**Expected effects.** Cheapest and numerically safe, but it suppresses the
very bottom-tail risk the income upgrade is intended to add.

**Implementation cost.** Low. **Reject as the paper baseline unless the bound
is measured, economically interpreted, and the affected preference parameter
is re-identified.**

## Ex ante M6 gates

Before any broad search, run a small exact-loop diagnostic with the proposed
income process and transfer block. Promotion requires:

1. zero empty-budget mass at every reachable state and age;
2. fiscal residual within the model's equilibrium tolerance;
3. the externally measured \(G_j(n)\) satisfies \(G_j(n)\geq B_j(n)\) wherever
   zero-resource states are reachable;
4. receipt and outlay moments within their assigned empirical uncertainty,
   including the young-childless/parent split;
5. \(\bar c_0\) interior to its search bounds and locally identified once the
   new moments are included;
6. no material deterioration in the full active wealth, fertility, ownership,
   and housing-response target system;
7. stable policy functions and the standard boundary-state diagnostic packet.

If gates 2--5 fail, do not launch more optimizer chains. The failure means the
current income support and preference floor cannot jointly describe the data.

## Paper-facing language

> We model a means-tested residual transfer whose statutory guarantee is
> measured independently of the Stone--Geary preference parameters, with
> recipients exhausting liquid assets. A proportional labor-income levy is
> determined in equilibrium to balance aggregate transfer outlays. We target
> receipt and spending by age and parenthood and reject parameterizations in
> which observed support cannot finance the admissible minimum bundle.

## Bottom line for tomorrow

Do not alter the M5 production run. For M6, first construct the statutory
guarantee and receipt/outlay targets; only then implement the funded residual-
transfer block on the fixed M5 reference point. The mechanism earns a
calibration run only if the small diagnostic passes the seven gates above.
