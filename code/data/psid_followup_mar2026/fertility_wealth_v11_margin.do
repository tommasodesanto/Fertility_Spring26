clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/fertility_wealth_v1"
capture log close _all
log using "`outdir'/fertility_wealth_v11.log", replace text

* ==================================================================
* V11: TWO MOTIVATIONAL FIGURES
*
* Figure A: Aggregate decomposition — where does childlessness live?
*   4 cells: owner×parent, owner×childless, renter×parent, renter×childless
*   + counterfactual: if stayer-renters had transitioner fertility
*
* Figure B: "Stuck" population by wealth bin
*   Share who are BOTH non-owners AND childless, by entry NW
*   Overlay with ownership transition rate
*
* Use expanded samples where possible:
*   - Parenthood by 40 (not 45) to gain ~25% more obs
*   - All young adults (not just renters) for Figure B variant
* ==================================================================

use ID year AGEREP SEX DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR NETWORTHR HOMEOWN IW using "`dta'", clear

drop if year > DEATHYEAR
drop if missing(AGEREP) | AGEREP < 18
keep if inrange(year, 1984, 2019)

replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3

gen first_birth_year = RELCHI1BYEAR
gen pre_first_birth = (missing(first_birth_year) | year < first_birth_year)
gen own = HOMEOWN

* ==================================================================
* BUILD PERSON-LEVEL DATA
* ==================================================================

* --- Entry characteristics at 25-30 ---
preserve
    keep if inrange(AGEREP, 25, 30)
    keep if !missing(own) & !missing(NETWORTHR)
    collapse (mean) own_rate_2530 = own ///
        (mean) nw_dollars = NETWORTHR ///
        (mean) inc_dollars = INCFAMR ///
        (mean) weight_young = IW, by(ID)
    gen renter_2530 = (own_rate_2530 < 0.5)
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

* --- Tenure at 31-40 ---
preserve
    keep if inrange(AGEREP, 31, 40) & !missing(own)
    collapse (mean) own_rate_3140 = own, by(ID)
    gen owner_by_40 = (own_rate_3140 >= 0.5) if !missing(own_rate_3140)
    keep ID owner_by_40
    tempfile t3140
    save `t3140', replace
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

* --- Merge ---
use `entry', clear
merge 1:1 ID using `t3135', keep(match master) nogen
merge 1:1 ID using `t3140', keep(match master) nogen
merge 1:1 ID using `fert', keep(match) nogen

gen observed_by_45 = (2019 - birth_year) >= 45
gen observed_by_40 = (2019 - birth_year) >= 40
gen observed_by_35 = (2019 - birth_year) >= 35

gen parent_by_45 = (!missing(age_at_first_birth) & age_at_first_birth <= 45)
gen parent_by_40 = (!missing(age_at_first_birth) & age_at_first_birth <= 40)
gen parent_by_35 = (!missing(age_at_first_birth) & age_at_first_birth <= 35)

gen nw_k = nw_dollars / 1000

di _n "=== SAMPLE SIZES ==="
di "All at 25-30 with wealth:"
count
di "Renters at 25-30:"
count if renter_2530 == 1
di "Renters, observed by 45, with tenure at 31-35:"
count if renter_2530 == 1 & observed_by_45 == 1 & !missing(owner_by_35)
di "Renters, observed by 40, with tenure at 31-35:"
count if renter_2530 == 1 & observed_by_40 == 1 & !missing(owner_by_35)
di "Renters, observed by 40, with tenure at 31-40:"
count if renter_2530 == 1 & observed_by_40 == 1 & !missing(owner_by_40)
di "All adults, observed by 45, with tenure at 31-35:"
count if observed_by_45 == 1 & !missing(owner_by_35)
di "All adults, observed by 40, with tenure at 31-35:"
count if observed_by_40 == 1 & !missing(owner_by_35)

* ==================================================================
* FIGURE A: AGGREGATE DECOMPOSITION
* Where does childlessness live among initial renters?
*
* Use: renters at 25-30, observed by 40, tenure measured at 31-35
* (larger sample than by-45)
* ==================================================================
di _n "============================================="
di    "FIGURE A: Aggregate decomposition"
di    "============================================="

preserve
    keep if renter_2530 == 1 & observed_by_40 == 1 & !missing(owner_by_35)

    * Trim 1-99 on NW
    _pctile nw_dollars [pw=weight_young], p(1 99)
    drop if nw_dollars < r(r1) | nw_dollars > r(r2)

    di _n "Sample size:"
    count

    * Four cells
    gen cell = .
    replace cell = 1 if owner_by_35 == 1 & parent_by_40 == 1
    replace cell = 2 if owner_by_35 == 1 & parent_by_40 == 0
    replace cell = 3 if owner_by_35 == 0 & parent_by_40 == 1
    replace cell = 4 if owner_by_35 == 0 & parent_by_40 == 0

    label define cell_lbl 1 "Owner, parent" 2 "Owner, childless" ///
        3 "Renter, parent" 4 "Renter, childless"
    label values cell cell_lbl

    tab cell [aw=weight_young]

    * Compute weighted shares
    collapse (count) n_obs = ID [pweight=weight_young], by(cell)
    egen total = sum(n_obs)
    gen share = n_obs / total
    list

    * Key numbers for annotation
    di _n "=== DECOMPOSITION ==="
    forval c = 1/4 {
        sum share if cell == `c'
        local sh`c' = r(mean)
        di "Cell `c': " %5.3f `sh`c''
    }

    * Childlessness rate among transitioners
    local cless_owners = `sh2' / (`sh1' + `sh2')
    * Childlessness rate among stayer-renters
    local cless_renters = `sh4' / (`sh3' + `sh4')
    * Share who stay renters
    local share_renters = `sh3' + `sh4'
    * Share who become owners
    local share_owners = `sh1' + `sh2'

    di _n "Childless rate, transitioners: " %5.3f `cless_owners'
    di "Childless rate, stayer-renters: " %5.3f `cless_renters'
    di "Share who stay renters: " %5.3f `share_renters'
    di "Share who become owners: " %5.3f `share_owners'

    * Counterfactual: if stayer-renters had transitioner fertility
    local cf_cless_renters = `share_renters' * `cless_owners'
    local observed_cless = `sh2' + `sh4'
    local cf_cless = `sh2' + `cf_cless_renters'
    local excess = `observed_cless' - `cf_cless'
    di _n "Observed childlessness: " %5.3f `observed_cless'
    di "Counterfactual childlessness: " %5.3f `cf_cless'
    di "Excess childlessness from ownership barrier: " %5.3f `excess'
    di "  = " %4.1f `=`excess'*100' " pp"

    * --- Bar chart: observed vs counterfactual ---
    clear
    set obs 4
    gen bar = _n
    gen height = .
    gen group = .

    * Bar 1: observed childlessness among owner-transitioners
    replace height = `sh2' in 1
    replace group = 1 in 1
    * Bar 2: observed childlessness among stayer-renters
    replace height = `sh4' in 2
    replace group = 1 in 2
    * Bar 3: counterfactual childlessness among stayer-renters
    replace height = `cf_cless_renters' in 3
    replace group = 2 in 3
    * Bar 4: same owner-childless (unchanged)
    replace height = `sh2' in 4
    replace group = 2 in 4

    * Actually, let me do this differently — paired bars

    clear
    set obs 2
    gen scenario = _n
    gen owner_childless = .
    gen renter_childless = .
    replace owner_childless = `sh2' in 1
    replace renter_childless = `sh4' in 1
    replace owner_childless = `sh2' in 2
    replace renter_childless = `cf_cless_renters' in 2

    gen total_childless = owner_childless + renter_childless

    list

    * Save for export
    export delimited using "`outdir'/decomposition_v11.csv", replace

    * Stacked bar
    graph bar owner_childless renter_childless, ///
        over(scenario, relabel(1 "Observed" 2 `""Counterfactual:" "renters get" "owner fertility""')) ///
        stack ///
        bar(1, fcolor(navy%60) lcolor(navy)) ///
        bar(2, fcolor(cranberry%60) lcolor(cranberry)) ///
        legend(order(2 "Childless stayer-renters" 1 "Childless owner-transitioners") ///
               position(6) rows(1) size(small)) ///
        ytitle("Share of all initial renters at 25-30") ///
        ylabel(0(0.05)0.45, angle(h) format(%4.2f)) ///
        title("Childlessness attributable to the ownership barrier") ///
        subtitle("PSID: renters at 25-30, parenthood by 40") ///
        note("Counterfactual assigns stayer-renters the owner-transitioner parenthood rate." ///
             "Excess childlessness = " %4.1f `=`excess'*100' " pp.") ///
        graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/decomposition_v11.png", replace width(1600)
restore

* ==================================================================
* FIGURE A2: Simpler version — just the two bars
* Childlessness rate by tenure transition, with N annotations
* ==================================================================
di _n "============================================="
di    "FIGURE A2: Simple two-bar childlessness gap"
di    "============================================="

preserve
    keep if renter_2530 == 1 & observed_by_40 == 1 & !missing(owner_by_35)
    _pctile nw_dollars [pw=weight_young], p(1 99)
    drop if nw_dollars < r(r1) | nw_dollars > r(r2)

    gen childless_by_40 = 1 - parent_by_40

    collapse (mean) childless_rate = childless_by_40 ///
        (count) n_obs = ID [pweight=weight_young], by(owner_by_35)
    list

    graph bar childless_rate, ///
        over(owner_by_35, relabel(1 "Stayed renter" 2 "Became owner") ///
             label(labsize(medium))) ///
        bar(1, fcolor(cranberry%70) lcolor(cranberry)) ///
        ytitle("Share childless by 40", size(medium)) ///
        ylabel(0(0.1)0.7, angle(h)) ///
        title("The ownership–fertility gap among initial renters") ///
        subtitle("PSID: renters at 25-30, followed to 40") ///
        blabel(bar, format(%4.2f) size(large)) ///
        graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/simple_gap_v11.png", replace width(1600)
restore

* ==================================================================
* FIGURE B: "STUCK" POPULATION BY WEALTH
* Share who are both non-owners AND childless, by entry NW bin
* This is the policy-movable margin
*
* Use: all young adults at 25-30, observed by 40 (largest sample)
* ==================================================================
di _n "============================================="
di    "FIGURE B: Stuck population by wealth"
di    "============================================="

preserve
    keep if observed_by_40 == 1 & !missing(owner_by_35)
    _pctile nw_dollars [pw=weight_young], p(1 99)
    drop if nw_dollars < r(r1) | nw_dollars > r(r2)

    di "Sample (all young adults, observed by 40):"
    count

    gen childless_by_40 = 1 - parent_by_40

    * "Stuck" = non-owner AND childless
    gen stuck = (owner_by_35 == 0 & childless_by_40 == 1)
    * Renter-parent (had kids without buying)
    gen renter_parent = (owner_by_35 == 0 & parent_by_40 == 1)
    * Owner-childless
    gen owner_childless = (owner_by_35 == 1 & childless_by_40 == 1)
    * Owner-parent
    gen owner_parent = (owner_by_35 == 1 & parent_by_40 == 1)

    * Dollar bins
    gen nw_bin = .
    replace nw_bin = 1 if nw_dollars < 0
    replace nw_bin = 2 if nw_dollars >= 0     & nw_dollars < 10000
    replace nw_bin = 3 if nw_dollars >= 10000  & nw_dollars < 25000
    replace nw_bin = 4 if nw_dollars >= 25000  & nw_dollars < 50000
    replace nw_bin = 5 if nw_dollars >= 50000  & nw_dollars < 100000
    replace nw_bin = 6 if nw_dollars >= 100000 & !missing(nw_dollars)

    tab nw_bin

    collapse (mean) stuck renter_parent owner_childless owner_parent ///
        (mean) owner_by_35 childless_by_40 ///
        (count) n_obs = ID [pweight=weight_young], by(nw_bin)

    * The "stuck" share should decline with wealth
    list

    export delimited using "`outdir'/stuck_by_wealth_v11.csv", replace

    * --- Version 1: Stacked bar showing all four cells ---
    graph bar owner_parent renter_parent owner_childless stuck, ///
        over(nw_bin, relabel(1 `""< $0""' 2 `""$0-" "10k""' 3 `""$10-" "25k""' ///
             4 `""$25-" "50k""' 5 `""$50-" "100k""' 6 `""$100k+""')) ///
        stack ///
        bar(1, fcolor(navy%40) lcolor(navy%60)) ///
        bar(2, fcolor(eltgreen%60) lcolor(eltgreen%80)) ///
        bar(3, fcolor(orange%50) lcolor(orange%70)) ///
        bar(4, fcolor(cranberry%70) lcolor(cranberry)) ///
        legend(order(4 "Non-owner, childless" 3 "Owner, childless" ///
                     2 "Non-owner, parent" 1 "Owner, parent") ///
               position(6) rows(1) size(vsmall)) ///
        ytitle("Share of young adults in each cell") ///
        ylabel(0(0.1)1, angle(h)) ///
        title("Where childlessness lives, by entry wealth") ///
        subtitle("All young adults at 25-30 (PSID), outcomes by 40") ///
        note("Non-owner childless (dark red) = policy-movable margin." ///
             "Sample: all young adults with wealth data at 25-30.") ///
        graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/stuck_stacked_v11.png", replace width(1600)

    * --- Version 2: Two-series — "stuck" share + ownership rate ---
    twoway (bar stuck nw_bin, barw(0.6) fcolor(cranberry%60) lcolor(cranberry)) ///
           (connected owner_by_35 nw_bin, ///
            mcolor(navy) lcolor(navy) msymbol(O) msize(vlarge) lwidth(thick) ///
            yaxis(1)) ///
        , legend(order(1 `""Non-owner &" "childless by 40""' ///
                       2 "Ownership rate by 35") ///
                 position(6) rows(1) size(medium)) ///
          xtitle("Net worth at 25-30 (2019 $)", size(medium)) ///
          ytitle("Share", size(medium)) ///
          xlabel(1 `""< $0""' 2 `""$0-10k""' 3 `""$10-25k""' 4 `""$25-50k""' ///
                 5 `""$50-100k""' 6 `""$100k+""', angle(20)) ///
          ylabel(0(0.1)0.9, angle(h)) ///
          title("The ownership barrier and childlessness") ///
          subtitle("All young adults at 25-30 (PSID 1984-2019)") ///
          note("Bars: share who remain non-owners AND childless by 40." ///
               "Line: ownership rate by 35. As ownership rises, the stuck share falls.") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/stuck_vs_ownership_v11.png", replace width(1600)
restore

* ==================================================================
* FIGURE B2: Same as B but RENTERS ONLY (initial renters at 25-30)
* More directly comparable to the model's constrained population
* ==================================================================
di _n "============================================="
di    "FIGURE B2: Stuck population by wealth (renters only)"
di    "============================================="

preserve
    keep if renter_2530 == 1 & observed_by_40 == 1 & !missing(owner_by_35)
    _pctile nw_dollars [pw=weight_young], p(1 99)
    drop if nw_dollars < r(r1) | nw_dollars > r(r2)

    di "Sample (renters, observed by 40):"
    count

    gen childless_by_40 = 1 - parent_by_40
    gen stuck = (owner_by_35 == 0 & childless_by_40 == 1)

    gen nw_bin = .
    replace nw_bin = 1 if nw_dollars < 0
    replace nw_bin = 2 if nw_dollars >= 0     & nw_dollars < 10000
    replace nw_bin = 3 if nw_dollars >= 10000  & nw_dollars < 25000
    replace nw_bin = 4 if nw_dollars >= 25000  & nw_dollars < 50000
    replace nw_bin = 5 if nw_dollars >= 50000  & nw_dollars < 100000
    replace nw_bin = 6 if nw_dollars >= 100000 & !missing(nw_dollars)

    collapse (mean) stuck ///
        (mean) owner_by_35 childless_by_40 ///
        (count) n_obs = ID [pweight=weight_young], by(nw_bin)

    list

    export delimited using "`outdir'/stuck_renters_v11.csv", replace

    twoway (bar stuck nw_bin, barw(0.6) fcolor(cranberry%60) lcolor(cranberry)) ///
           (connected owner_by_35 nw_bin, ///
            mcolor(navy) lcolor(navy) msymbol(O) msize(vlarge) lwidth(thick)) ///
        , legend(order(1 `""Non-owner &" "childless by 40""' ///
                       2 `""Ownership" "transition rate""') ///
                 position(6) rows(1) size(medium)) ///
          xtitle("Net worth at 25-30 (2019 $)", size(medium)) ///
          ytitle("Share", size(medium)) ///
          xlabel(1 `""< $0""' 2 `""$0-10k""' 3 `""$10-25k""' 4 `""$25-50k""' ///
                 5 `""$50-100k""' 6 `""$100k+""', angle(20)) ///
          ylabel(0(0.1)0.7, angle(h)) ///
          title("The ownership barrier and childlessness") ///
          subtitle("Among renters at 25-30 (PSID 1984-2019)") ///
          note("Bars: share who remain non-owners AND childless by 40." ///
               "As ownership rises with wealth, the stuck share falls.") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/stuck_renters_vs_own_v11.png", replace width(1600)
restore

* ==================================================================
* FIGURE C: LPOLY — dual smooth, ALL young adults, observed by 40
* Largest sample for maximum smoothness
* ==================================================================
di _n "============================================="
di    "FIGURE C: Local polynomial, all adults, by 40"
di    "============================================="

preserve
    keep if observed_by_40 == 1 & !missing(owner_by_35)
    _pctile nw_dollars [pw=weight_young], p(2 98)
    drop if nw_dollars < r(r1) | nw_dollars > r(r2)

    di "Sample for lpoly (all, by 40):"
    count

    gen childless_by_40 = 1 - parent_by_40
    gen stuck = (owner_by_35 == 0 & childless_by_40 == 1)
    gen nw_plot = min(nw_dollars/1000, 200)

    * Dual smooth: ownership and parenthood
    twoway (lpoly owner_by_35 nw_plot [aw=weight_young], ///
            lcolor(navy) lwidth(thick) bwidth(25) degree(1)) ///
           (lpoly parent_by_40 nw_plot [aw=weight_young], ///
            lcolor(cranberry) lwidth(thick) lpattern(dash) bwidth(25) degree(1)) ///
           (lpoly stuck nw_plot [aw=weight_young], ///
            lcolor(gs6) lwidth(medthick) lpattern(shortdash) bwidth(25) degree(1)) ///
        , legend(order(1 "Pr(own by 35)" 2 "Pr(parent by 40)" ///
                       3 `""Pr(non-owner" "& childless)""') ///
                 position(6) rows(1) size(small)) ///
          xtitle("Net worth at 25-30 ($1000s, 2019 real)", size(medium)) ///
          ytitle("Probability", size(medium)) ///
          ylabel(0(0.1)1.0, angle(h)) ///
          xline(0, lcolor(gs10) lpattern(dot)) ///
          title("Ownership, parenthood, and the stuck margin") ///
          subtitle("Local polynomial smooth, all young adults at 25-30") ///
          note("Local linear, Epanechnikov kernel, bw = $25k. Trimmed p2/p98.") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/lpoly_three_v11.png", replace width(1600)
restore

log close
