clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/dp_at_first_purchase_v1"
capture mkdir "`outdir'"
capture log close _all
log using "`outdir'/dp_at_first_purchase_v1.log", replace text

* ==================================================================
* DP-AT-FIRST-PURCHASE (V1)
*
* Object: density of down-payment share (1 - mortgage / value) at the
* year of first owner-tenure spell. Mass at <= 5/10/20 percent =
* direct evidence that buyers bunch at the collateral frontier.
*
* Sample: PSID heads observed as renter at least once and then
* observed as owner; first owner-year per ID; YEARMOVED_ within the
* same biennial wave window so we are close to origination.
* ==================================================================

use ID year AGEREP SEX HOMEOWN HOMEVALUER HOMEMORTOTR HOMEEQUITYR ///
    HOMEMOR1MR HOMEMOR2MR YEARMOVED_ IW ///
    using "`dta'", clear

* Recode tenure: owner = 1, renter = 0, other = missing
gen byte own = .
replace own = 1 if HOMEOWN == 1
replace own = 0 if HOMEOWN == 2
* HOMEOWN == 3 (neither own nor rent) -> missing

keep if inrange(year, 1984, 2019)
keep if !missing(AGEREP) & inrange(AGEREP, 18, 65)

sort ID year

* ------------------------------------------------------------------
* Identify each ID's first observed owner-year
* and require a prior renter-year in the panel
* ------------------------------------------------------------------
by ID: gen any_renter_before = sum(own == 0)
by ID: gen first_own_year = year if own == 1 & any_renter_before[_n-1] >= 1 & ///
    sum(own == 1) == 1

* Flag the row that is the first owner-spell start (preceded by renter)
gen byte first_purchase = !missing(first_own_year)

* ------------------------------------------------------------------
* Require move-in year close to first owner-year, so reported value
* and mortgage are at/near origination, not 4+ years post-purchase.
* PSID is biennial post-1997; allow YEARMOVED_ in [year-2, year].
* ------------------------------------------------------------------
gen byte recent_move = inrange(YEARMOVED_, year - 2, year) if !missing(YEARMOVED_)

* ------------------------------------------------------------------
* Build sample and DP measure
* ------------------------------------------------------------------
keep if first_purchase == 1

di _n "=== First-purchase rows before screens ==="
count

keep if recent_move == 1

di _n "=== After recent-move screen (YEARMOVED_ in [t-2, t]) ==="
count

keep if HOMEVALUER > 0 & !missing(HOMEVALUER) & !missing(HOMEMORTOTR)

* Construct DP share at first observation
gen dp_share = 1 - HOMEMORTOTR / HOMEVALUER
gen byte cash_buyer = HOMEMORTOTR == 0 & HOMEVALUER > 0

di _n "=== After value > 0 and mortgage non-missing ==="
count
sum dp_share, detail

* Trim implausible (DP > 1 means negative reported mortgage; DP < 0 means
* mortgage > value at origination, possible with reporting noise but
* rare and obscures the bunching plot)
gen byte in_range = inrange(dp_share, 0, 1.5)

di _n "=== Distribution of DP-share, in [0, 1.5] ==="
sum dp_share if in_range, detail

* ------------------------------------------------------------------
* Headline stats: bunching at collateral floors
* ------------------------------------------------------------------
preserve
    keep if in_range
    gen byte at_or_below_035 = dp_share <= 0.035
    gen byte at_or_below_050 = dp_share <= 0.05
    gen byte at_or_below_100 = dp_share <= 0.10
    gen byte at_or_below_200 = dp_share <= 0.20

    di _n "=== Share of first-time buyers at/below collateral floors ==="
    sum at_or_below_035 at_or_below_050 at_or_below_100 at_or_below_200 ///
        cash_buyer [aweight=IW]

    collapse (mean) share_le_035 = at_or_below_035 ///
             (mean) share_le_050 = at_or_below_050 ///
             (mean) share_le_100 = at_or_below_100 ///
             (mean) share_le_200 = at_or_below_200 ///
             (mean) share_cash   = cash_buyer ///
             (p50) median_dp_share = dp_share ///
             (count) n_first_purchases = ID ///
        [aweight=IW]
    list
    export delimited using "`outdir'/dp_share_bunching_summary_v1.csv", replace
restore

* ------------------------------------------------------------------
* FIGURE 1: kernel density of DP-share at first purchase
* Vertical lines at common collateral floors
* ------------------------------------------------------------------
preserve
    keep if in_range
    twoway (kdensity dp_share [aweight=IW], lwidth(thick) lcolor(navy) ///
            bwidth(0.025)) ///
        , xline(0.035, lcolor(gs9) lpattern(dot)) ///
          xline(0.05,  lcolor(gs9) lpattern(dot)) ///
          xline(0.10,  lcolor(gs9) lpattern(dot)) ///
          xline(0.20,  lcolor(cranberry) lpattern(dash)) ///
          xtitle("Down-payment share at first purchase (1 - mortgage / value)") ///
          ytitle("Density") ///
          xlabel(0(0.1)1) ///
          title("Bunching at the collateral frontier") ///
          subtitle("PSID first-time owners, 1984-2019, ages 18-65, YEARMOVED in [t-2, t]") ///
          note("Dotted lines: 3.5% (FHA), 5%, 10% conventional floors. Dashed: 20% GSE.") ///
          graphregion(color(white)) scheme(s1color) legend(off)
    graph export "`outdir'/dp_share_density_v1.png", replace width(1600)
restore

* ------------------------------------------------------------------
* FIGURE 2: histogram with bins around the policy thresholds
* ------------------------------------------------------------------
preserve
    keep if in_range
    twoway (histogram dp_share [fweight=int(IW)], width(0.025) start(0) ///
            fcolor(navy%40) lcolor(navy)) ///
        , xline(0.20, lcolor(cranberry) lpattern(dash)) ///
          xline(0.10, lcolor(gs9) lpattern(dot)) ///
          xline(0.05, lcolor(gs9) lpattern(dot)) ///
          xtitle("Down-payment share at first purchase") ///
          ytitle("Density") ///
          xlabel(0(0.1)1) ///
          title("Down-payment share at first purchase") ///
          subtitle("PSID first-time owners, 1984-2019") ///
          graphregion(color(white)) scheme(s1color) legend(off)
    graph export "`outdir'/dp_share_histogram_v1.png", replace width(1600)
restore

* ------------------------------------------------------------------
* FIGURE 3: kernel density EXCLUDING cash buyers
* Cash buyers (DP = 100%) are a separate group and the 1.0 spike
* visually distracts. Conditional on a mortgage, the collateral-floor
* story is sharper.
* ------------------------------------------------------------------
* Block 3a: stats + CSV (no nested preserve)
preserve
    keep if in_range
    keep if cash_buyer == 0
    di _n "=== Mortgaged first-time buyers (no cash) ==="
    count

    gen byte mort_le_035 = dp_share <= 0.035
    gen byte mort_le_050 = dp_share <= 0.05
    gen byte mort_le_100 = dp_share <= 0.10
    gen byte mort_le_200 = dp_share <= 0.20

    collapse (mean) share_le_035 = mort_le_035 ///
             (mean) share_le_050 = mort_le_050 ///
             (mean) share_le_100 = mort_le_100 ///
             (mean) share_le_200 = mort_le_200 ///
             (p50)  median_dp = dp_share ///
             (count) n_mortgaged = ID [aweight=IW]
    list
    export delimited using ///
        "`outdir'/dp_share_bunching_summary_nocash_v1.csv", replace
restore

* Block 3b: figure (separate top-level preserve)
preserve
    keep if in_range
    keep if cash_buyer == 0
    twoway (kdensity dp_share [aweight=IW], lwidth(thick) lcolor(navy) ///
            bwidth(0.025)) ///
        , xline(0.035, lcolor(gs9) lpattern(dot)) ///
          xline(0.05,  lcolor(gs9) lpattern(dot)) ///
          xline(0.10,  lcolor(gs9) lpattern(dot)) ///
          xline(0.20,  lcolor(cranberry) lpattern(dash)) ///
          xtitle("Down-payment share at first purchase (mortgaged buyers)") ///
          ytitle("Density") ///
          xlabel(0(0.1)0.9) ///
          title("Bunching at the collateral frontier (mortgaged buyers)") ///
          subtitle("PSID first-time buyers with a mortgage, 1984-2019") ///
          note("Dotted lines: 3.5% (FHA), 5%, 10% conventional floors. Dashed: 20% GSE." ///
               "Cash buyers (DP = 100%) excluded; reported separately.") ///
          graphregion(color(white)) scheme(s1color) legend(off)
    graph export "`outdir'/dp_share_density_nocash_v1.png", replace width(1600)
restore

* ------------------------------------------------------------------
* TABLE: by year of first purchase, to show the bunching is not just
* a 2007-2010 crisis-era artifact.
* ------------------------------------------------------------------
preserve
    keep if in_range
    gen year_bin = .
    replace year_bin = 1 if inrange(year, 1984, 1994)
    replace year_bin = 2 if inrange(year, 1995, 2004)
    replace year_bin = 3 if inrange(year, 2005, 2010)
    replace year_bin = 4 if inrange(year, 2011, 2019)
    label define yb 1 "1984-1994" 2 "1995-2004" 3 "2005-2010" 4 "2011-2019"
    label values year_bin yb

    gen byte at_or_below_200 = dp_share <= 0.20
    gen byte at_or_below_100 = dp_share <= 0.10

    collapse (mean) share_le_200 = at_or_below_200 ///
             (mean) share_le_100 = at_or_below_100 ///
             (p50) median_dp = dp_share ///
             (count) n = ID [aweight=IW], by(year_bin)
    list
    export delimited using "`outdir'/dp_share_by_year_bin_v1.csv", replace
restore

log close
