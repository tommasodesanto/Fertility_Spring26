clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/fertility_wealth_v1"
capture log close _all
log using "`outdir'/fertility_wealth_v10b.log", replace text

* ==================================================================
* V10b: TIMING-BASED DOSE-RESPONSE + COMPOSITE FIGURE
*
* Key insight from v10: parenthood by 45 is ~80% at all wealth levels.
* But TIMING differs. Try parenthood by 35 as the fertility outcome.
* Also: composite figure combining ownership gradient + within-wealth gap.
* ==================================================================

use ID year AGEREP SEX DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR NETWORTHR NETWORTH2R HOMEOWN IW using "`dta'", clear

drop if year > DEATHYEAR
drop if missing(AGEREP) | AGEREP < 18
keep if inrange(year, 1984, 2019)

replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3

gen first_birth_year = RELCHI1BYEAR
gen pre_first_birth = (missing(first_birth_year) | year < first_birth_year)
gen own = HOMEOWN

* --- Entry characteristics at 25-30 ---
preserve
    keep if inrange(AGEREP, 25, 30)
    keep if !missing(own) & !missing(NETWORTHR)
    collapse (mean) own_rate_2530 = own ///
        (mean) nw_dollars = NETWORTHR ///
        (mean) inc_dollars = INCFAMR ///
        (mean) weight_young = IW, by(ID)
    gen renter_2530 = (own_rate_2530 < 0.5)
    gen nw_k = nw_dollars / 1000
    tempfile entry
    save `entry', replace
restore

* --- Tenure at 31-35 ---
preserve
    keep if inrange(AGEREP, 31, 35) & !missing(own)
    collapse (mean) own_rate_3135 = own, by(ID)
    gen owner_by_35 = (own_rate_3135 >= 0.5) if !missing(own_rate_3135)
    keep ID owner_by_35
    tempfile t3135
    save `t3135', replace
restore

* --- Fertility ---
preserve
    collapse (min) min_year = year (min) min_age = AGEREP ///
        (first) first_birth_year_b = first_birth_year, by(ID)
    gen birth_year = min_year - min_age
    gen age_at_first_birth = first_birth_year_b - birth_year if !missing(first_birth_year_b)
    keep ID birth_year age_at_first_birth
    tempfile fert
    save `fert', replace
restore

use `entry', clear
merge 1:1 ID using `t3135', keep(match master) nogen
merge 1:1 ID using `fert', keep(match) nogen

gen observed_by_45 = (2019 - birth_year) >= 45
gen observed_by_40 = (2019 - birth_year) >= 40
gen observed_by_35 = (2019 - birth_year) >= 35

gen parent_by_45 = (!missing(age_at_first_birth) & age_at_first_birth <= 45)
gen parent_by_40 = (!missing(age_at_first_birth) & age_at_first_birth <= 40)
gen parent_by_35 = (!missing(age_at_first_birth) & age_at_first_birth <= 35)
gen childless_by_45 = 1 - parent_by_45

* Trim
_pctile nw_dollars [pw=weight_young], p(1 99)
local lo = r(r1)
local hi = r(r2)

* ==================================================================
* 1. TIMING: Parenthood by 35 as outcome (more variation than by 45)
*    Sample: all young adults observed by 35
* ==================================================================
di _n "============================================="
di    "1. DOLLAR-BIN: Ownership by 35 + Parenthood by 35"
di    "============================================="

preserve
    keep if observed_by_35 == 1 & !missing(owner_by_35)
    drop if nw_dollars < `lo' | nw_dollars > `hi'

    di "Sample (all, observed by 35):"
    count

    gen nw_bin = .
    replace nw_bin = 1 if nw_dollars < 0
    replace nw_bin = 2 if nw_dollars >= 0     & nw_dollars < 10000
    replace nw_bin = 3 if nw_dollars >= 10000  & nw_dollars < 25000
    replace nw_bin = 4 if nw_dollars >= 25000  & nw_dollars < 50000
    replace nw_bin = 5 if nw_dollars >= 50000  & nw_dollars < 100000
    replace nw_bin = 6 if nw_dollars >= 100000 & !missing(nw_dollars)

    tab nw_bin

    collapse (mean) owner_by_35 parent_by_35 ///
        (p50) median_nw_k = nw_k ///
        (count) n_obs = ID [pweight=weight_young], by(nw_bin)
    list

    export delimited using "`outdir'/dose_response_by35_v10b.csv", replace

    twoway (bar owner_by_35 nw_bin, barw(0.6) fcolor(navy%50) lcolor(navy)) ///
           (connected parent_by_35 nw_bin, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(D) msize(vlarge) lwidth(thick)) ///
        , legend(order(1 "Homeowner by 35" 2 "Parent by 35") ///
                 position(6) rows(1) size(medium)) ///
          xtitle("Net worth at 25-30 (2019 $)", size(medium)) ///
          ytitle("Share", size(medium)) ///
          xlabel(1 `""< $0""' 2 `""$0-10k""' 3 `""$10-25k""' 4 `""$25-50k""' ///
                 5 `""$50-100k""' 6 `""$100k+""', angle(20)) ///
          ylabel(0(0.1)0.9, angle(h)) ///
          title("Wealth predicts both ownership and early parenthood") ///
          subtitle("All young adults at 25-30 (PSID 1984-2019)") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/dose_response_by35_v10b.png", replace width(1600)
restore

* ==================================================================
* 2. TIMING: Renters only, parenthood by 35
* ==================================================================
di _n "============================================="
di    "2. RENTERS: Ownership transition + Parenthood by 35"
di    "============================================="

preserve
    keep if renter_2530 == 1 & observed_by_35 == 1 & !missing(owner_by_35)
    drop if nw_dollars < `lo' | nw_dollars > `hi'

    di "Sample (renters, observed by 35):"
    count

    gen nw_bin = .
    replace nw_bin = 1 if nw_dollars < 0
    replace nw_bin = 2 if nw_dollars >= 0     & nw_dollars < 10000
    replace nw_bin = 3 if nw_dollars >= 10000  & nw_dollars < 25000
    replace nw_bin = 4 if nw_dollars >= 25000  & nw_dollars < 50000
    replace nw_bin = 5 if nw_dollars >= 50000  & nw_dollars < 100000
    replace nw_bin = 6 if nw_dollars >= 100000 & !missing(nw_dollars)

    tab nw_bin

    collapse (mean) owner_by_35 parent_by_35 ///
        (p50) median_nw_k = nw_k ///
        (count) n_obs = ID [pweight=weight_young], by(nw_bin)
    list

    export delimited using "`outdir'/dose_response_renters_by35_v10b.csv", replace

    twoway (bar owner_by_35 nw_bin, barw(0.6) fcolor(navy%50) lcolor(navy)) ///
           (connected parent_by_35 nw_bin, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(D) msize(vlarge) lwidth(thick)) ///
        , legend(order(1 "Became owner by 35" 2 "Became parent by 35") ///
                 position(6) rows(1) size(medium)) ///
          xtitle("Net worth at 25-30 (2019 $)", size(medium)) ///
          ytitle("Share", size(medium)) ///
          xlabel(1 `""< $0""' 2 `""$0-10k""' 3 `""$10-25k""' 4 `""$25-50k""' ///
                 5 `""$50-100k""' 6 `""$100k+""', angle(20)) ///
          ylabel(0(0.1)0.9, angle(h)) ///
          title("The down-payment margin: wealth, ownership, and fertility") ///
          subtitle("Among renters at 25-30, outcomes by age 35") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/dose_response_renters_by35_v10b.png", replace width(1600)
restore

* ==================================================================
* 3. LPOLY: All young adults, parenthood by 35 (timing outcome)
* ==================================================================
di _n "============================================="
di    "3. LPOLY: Ownership + parenthood by 35"
di    "============================================="

preserve
    keep if observed_by_35 == 1 & !missing(owner_by_35)
    drop if nw_dollars < `lo' | nw_dollars > `hi'

    di "Sample for lpoly (all, by 35):"
    count

    gen nw_plot = min(nw_dollars/1000, 200)

    twoway (lpoly owner_by_35 nw_plot [aw=weight_young], ///
            lcolor(navy) lwidth(thick) bwidth(25) degree(1)) ///
           (lpoly parent_by_35 nw_plot [aw=weight_young], ///
            lcolor(cranberry) lwidth(thick) lpattern(dash) bwidth(25) degree(1)) ///
        , legend(order(1 "Pr(own by 35)" 2 "Pr(parent by 35)") ///
                 position(6) rows(1) size(medium)) ///
          xtitle("Net worth at 25-30 ($1000s, 2019 real)", size(medium)) ///
          ytitle("Probability", size(medium)) ///
          ylabel(0(0.1)1.0, angle(h)) ///
          xline(0, lcolor(gs10) lpattern(dot)) ///
          title("Wealth, ownership, and early parenthood") ///
          subtitle("Local polynomial smooth, all young adults at 25-30") ///
          note("Epanechnikov kernel, bandwidth = $25k.") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/lpoly_by35_v10b.png", replace width(1600)
restore

* ==================================================================
* 4. COMPOSITE: Ownership transition + childlessness by transition,
*    within wealth bins. Three series in one figure.
* ==================================================================
di _n "============================================="
di    "4. COMPOSITE: 3 series by wealth bin"
di    "============================================="

preserve
    keep if renter_2530 == 1 & observed_by_45 == 1 & !missing(owner_by_35)
    drop if nw_dollars < `lo' | nw_dollars > `hi'

    gen nw_bin = .
    replace nw_bin = 1 if nw_dollars < 0
    replace nw_bin = 2 if nw_dollars >= 0     & nw_dollars < 10000
    replace nw_bin = 3 if nw_dollars >= 10000  & nw_dollars < 25000
    replace nw_bin = 4 if nw_dollars >= 25000  & nw_dollars < 50000
    replace nw_bin = 5 if nw_dollars >= 50000  & nw_dollars < 100000
    replace nw_bin = 6 if nw_dollars >= 100000 & !missing(nw_dollars)

    * Childlessness by transition status within bin
    collapse (mean) childless_by_45 owner_by_35 ///
        (count) n_obs = ID [pweight=weight_young], by(nw_bin owner_by_35)

    * Also need overall ownership rate per bin
    bysort nw_bin: egen own_rate = total(n_obs * owner_by_35)
    bysort nw_bin: egen total_n = total(n_obs)
    replace own_rate = own_rate / total_n

    list, sep(2)
    export delimited using "`outdir'/composite_v10b.csv", replace

    * Connected plot: childlessness for owners vs renters by wealth bin
    twoway (connected childless_by_45 nw_bin if owner_by_35==1, ///
            mcolor(navy) lcolor(navy) msymbol(O) msize(large) lwidth(medthick)) ///
           (connected childless_by_45 nw_bin if owner_by_35==0, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(S) msize(large) lwidth(medthick)) ///
        , legend(order(1 "Became owner by 35" 2 "Stayed renter by 35") ///
                 position(6) rows(1) size(medium)) ///
          xtitle("Net worth at 25-30 (2019 $)", size(medium)) ///
          ytitle("Share childless by 45", size(medium)) ///
          xlabel(1 `""< $0""' 2 `""$0-10k""' 3 `""$10-25k""' 4 `""$25-50k""' ///
                 5 `""$50-100k""' 6 `""$100k+""', angle(20)) ///
          ylabel(0(0.1)0.8, angle(h)) ///
          title("Ownership gap in childlessness by entry wealth") ///
          subtitle("Among renters at 25-30, PSID 1984-2019") ///
          note("At every wealth level, renters who transition to ownership" ///
               "are 20-35pp less likely to remain childless.") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/composite_v10b.png", replace width(1600)
restore

* ==================================================================
* 5. CONSTRAINT MARGIN: Where does ownership bite most?
*    Show the GAP in childlessness (renter - owner) by wealth bin,
*    alongside the ownership transition rate.
* ==================================================================
di _n "============================================="
di    "5. CONSTRAINT MARGIN: gap + ownership rate"
di    "============================================="

preserve
    keep if renter_2530 == 1 & observed_by_45 == 1 & !missing(owner_by_35)
    drop if nw_dollars < `lo' | nw_dollars > `hi'

    gen nw_bin = .
    replace nw_bin = 1 if nw_dollars < 0
    replace nw_bin = 2 if nw_dollars >= 0     & nw_dollars < 10000
    replace nw_bin = 3 if nw_dollars >= 10000  & nw_dollars < 25000
    replace nw_bin = 4 if nw_dollars >= 25000  & nw_dollars < 50000
    replace nw_bin = 5 if nw_dollars >= 50000  & nw_dollars < 100000
    replace nw_bin = 6 if nw_dollars >= 100000 & !missing(nw_dollars)

    collapse (mean) childless_by_45 owner_by_35 ///
        (count) n_obs = ID [pweight=weight_young], by(nw_bin owner_by_35)

    reshape wide childless_by_45 n_obs, i(nw_bin) j(owner_by_35)
    gen gap = childless_by_451 - childless_by_450
    rename childless_by_450 cless_renter
    rename childless_by_451 cless_owner

    * Need ownership rate per bin
    gen own_rate = n_obs1 / (n_obs0 + n_obs1)

    list

    export delimited using "`outdir'/margin_gap_v10b.csv", replace

    twoway (bar own_rate nw_bin, barw(0.6) fcolor(navy%40) lcolor(navy)) ///
           (connected gap nw_bin, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(D) msize(vlarge) lwidth(thick) ///
            yaxis(2)) ///
        , legend(order(1 "Ownership transition rate" 2 "Owner-renter childlessness gap") ///
                 position(6) rows(1) size(small)) ///
          xtitle("Net worth at 25-30 (2019 $)", size(medium)) ///
          ytitle("Ownership rate", size(medium) axis(1)) ///
          ytitle("Childlessness gap (owner - renter)", size(medium) axis(2)) ///
          xlabel(1 `""< $0""' 2 `""$0-10k""' 3 `""$10-25k""' 4 `""$25-50k""' ///
                 5 `""$50-100k""' 6 `""$100k+""', angle(20)) ///
          ylabel(0(0.1)0.7, angle(h) axis(1)) ///
          ylabel(-0.5(0.1)0, angle(h) axis(2)) ///
          title("The ownership-fertility constraint margin") ///
          subtitle("Among renters at 25-30") ///
          note("Gap = owner childlessness - renter childlessness (negative = owners less childless)." ///
               "Wealth increases ownership rate; at every level, owners are less childless.") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/margin_gap_v10b.png", replace width(1600)
restore

log close
