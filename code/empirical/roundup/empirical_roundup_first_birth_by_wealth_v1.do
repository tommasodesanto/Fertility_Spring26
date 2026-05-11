clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "`root'/Fertility_Spring26/april_26_discrete_time/output"
local outdir "`out_root'/empirical_roundup_first_birth_by_wealth_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/empirical_roundup_first_birth_by_wealth_v1.log", replace text

global weight IW
global ref_sa 2
global bg "250 250 250"
global graphregion_slides graphregion(fcolor("$bg") lcolor("$bg")) plotregion(fcolor("$bg") lcolor("$bg"))

cap which eventstudyinteract
if _rc {
    di as error "eventstudyinteract is not installed. Aborting."
    exit 198
}

cap which svmat2
if _rc {
    di as error "svmat2 is not installed. Aborting."
    exit 198
}

use ID year AGEREP EDUYEAR SEX DEATHYEAR RELCHI1BYEAR RELCHI2BYEAR ///
    ACTUALROOMS_ HOMEOWN WHYMOVED1_ NETWORTH2R ${weight} using "`dta'", clear

replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3

drop if year > DEATHYEAR
drop if AGEREP < 18
keep if inrange(year, 1984, 2019)

bysort ID: egen first_child_year = min(RELCHI1BYEAR)
bysort ID: egen second_child_year = min(RELCHI2BYEAR)
drop if missing(first_child_year)

bysort ID: egen year_entry = min(year)
drop if first_child_year < year_entry
drop year_entry

drop if !missing(second_child_year) & second_child_year == first_child_year

rename ACTUALROOMS_ rooms
rename HOMEOWN own

gen moved_for_size = 1 if WHYMOVED1_ == 3 & WHYMOVED1_ != . & year > 1974
replace moved_for_size = 0 if WHYMOVED1_ != 3 & WHYMOVED1_ != . & year > 1974

xtset ID year
gen K = year - first_child_year

preserve
    keep if K == -2
    keep ID first_child_year NETWORTH2R ${weight}
    rename NETWORTH2R networth2r_pre
    rename ${weight} iw_pre
    keep if !missing(networth2r_pre, iw_pre)

    _pctile networth2r_pre [pw=iw_pre], p(1 99)
    keep if networth2r_pre >= r(r1) & networth2r_pre <= r(r2)

    xtile wealth_q = networth2r_pre [pw=iw_pre], n(5)

    tempfile baseline
    save `baseline', replace

    collapse (p50) median_nw_pre = networth2r_pre ///
        (count) n_obs = ID [pweight=iw_pre], by(wealth_q)
    export delimited using "`outdir'/wealth_quintile_summary_v1.csv", replace
restore

merge m:1 ID first_child_year using `baseline', keep(match) nogen

tempfile analysis
save `analysis', replace

tempname posth
postfile `posth' ///
    str18 outcome byte wealth_q ///
    double pre_event_mean coef_p3 se_p3 coef_p5 se_p5 sample_obs sample_ids rc ///
    using "`outdir'/first_birth_by_wealth_summary_v1.dta", replace

local outcomes "rooms own moved_for_size"

foreach outcome in `outcomes' {
    di as text "Running outcome: `outcome'"

    forvalues q = 1/5 {
        di as text "  Wealth quintile `q'"
        use `analysis', clear
        keep if wealth_q == `q'

        sum first_child_year
        gen lastcohort = first_child_year == r(max)

        cap drop L*event F*event
        forvalues l = 0/10 {
            gen L`l'event = K == `l'
        }
        gen L11event = (K > 10 & K != .)
        forvalues l = 1/5 {
            gen F`l'event = K == -`l'
        }
        gen F7event = (K < -6 & K != .)
        drop F${ref_sa}event

        capture noisily eventstudyinteract `outcome' L*event F*event, ///
            vce(cluster ID) absorb(year) cohort(first_child_year) control_cohort(lastcohort) ///
            covariates(i.AGEREP i.EDUYEAR)

        local rc_run = _rc
        if `rc_run' != 0 {
            post `posth' ("`outcome'") (`q') (.) (.) (.) (.) (.) (.) (.) (`rc_run')
            continue
        }

        quietly summarize `outcome' if e(sample) & K == -${ref_sa}
        local pre_event_mean = r(mean)

        quietly count if e(sample)
        local sample_obs = r(N)

        egen tag_id = tag(ID) if e(sample)
        quietly count if tag_id == 1
        local sample_ids = r(N)
        drop tag_id

        matrix b = e(b_iw)
        matrix var = e(V_iw)
        matrix vardiag = vecdiag(var)
        matrix combine = b \ vardiag
        matrix rownames combine = b variance
        matrix results = combine'
        matrix drop combine

        preserve
            clear
            svmat2 results, names(col) rnames(coeff)
            gen outcome = "`outcome'"
            gen wealth_q = `q'
            gen pre_event_mean = `pre_event_mean'
            gen sample_obs = `sample_obs'
            gen sample_ids = `sample_ids'
            gen se = sqrt(variance)
            drop variance

            replace b = . if b == 0 & se == 0
            replace se = . if b == . & se == 0

            drop if coeff == "_cons"
            gen relative_time = subinstr(coeff, "event", "", .)
            replace relative_time = subinstr(relative_time, "L", "", .)
            replace relative_time = subinstr(relative_time, "F", "-", .)
            replace relative_time = subinstr(relative_time, "o.", "", .)
            destring relative_time, replace

            assert relative_time != -${ref_sa}
            assert relative_time != .

            set obs `=_N+1'
            replace relative_time = -${ref_sa} if missing(relative_time)
            replace b = 0 if relative_time == -${ref_sa}
            replace se = 0 if relative_time == -${ref_sa}

            gen ci_lo = b - 1.96*se
            gen ci_hi = b + 1.96*se

            sort relative_time
            save "`outdir'/`outcome'_q`q'_estimates_v1.dta", replace

            drop if missing(b)
            twoway ///
                (rarea ci_lo ci_hi relative_time, fcolor(dkgreen%10) lcolor(dkgreen%10)) ///
                (connected b relative_time, mcolor(dkgreen) msymbol(circle) lcolor(dkgreen) lwidth(medthick) lpattern(dash)), ///
                ytitle("`outcome'") xtitle("Years relative to first birth") ///
                xline(-1.5, lwidth(thin) lpattern(dash) lcolor(gray)) ///
                yline(0, lwidth(vthin) lpattern(solid) lcolor(gray)) ///
                title("`outcome', pre-birth wealth quintile `q'") ///
                $graphregion_slides legend(off)
            graph export "`outdir'/`outcome'_q`q'_v1.png", replace

            quietly summarize b if relative_time == 3
            local coef_p3 = r(mean)
            quietly summarize se if relative_time == 3
            local se_p3 = r(mean)
            quietly summarize b if relative_time == 5
            local coef_p5 = r(mean)
            quietly summarize se if relative_time == 5
            local se_p5 = r(mean)
        restore

        post `posth' ("`outcome'") (`q') ///
            (`pre_event_mean') (`coef_p3') (`se_p3') (`coef_p5') (`se_p5') ///
            (`sample_obs') (`sample_ids') (0)
    }

    preserve
        clear
        local any_files 0
        forvalues q = 1/5 {
            capture append using "`outdir'/`outcome'_q`q'_estimates_v1.dta"
            if _rc == 0 {
                local any_files 1
            }
        }

        if `any_files' == 1 {
            keep wealth_q relative_time b
            drop if missing(b)

            twoway ///
                (connected b relative_time if wealth_q == 1, lcolor(cranberry) mcolor(cranberry) msymbol(circle)) ///
                (connected b relative_time if wealth_q == 2, lcolor(navy) mcolor(navy) msymbol(triangle)) ///
                (connected b relative_time if wealth_q == 3, lcolor(dkgreen) mcolor(dkgreen) msymbol(square)) ///
                (connected b relative_time if wealth_q == 4, lcolor(orange) mcolor(orange) msymbol(diamond)) ///
                (connected b relative_time if wealth_q == 5, lcolor(black) mcolor(black) msymbol(plus)), ///
                ytitle("`outcome'") xtitle("Years relative to first birth") ///
                xline(-1.5, lwidth(thin) lpattern(dash) lcolor(gray)) ///
                yline(0, lwidth(vthin) lpattern(solid) lcolor(gray)) ///
                legend(order(1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4" 5 "Q5")) ///
                title("First-birth event study by pre-birth wealth: `outcome'") ///
                $graphregion_slides
            graph export "`outdir'/`outcome'_overlay_v1.png", replace
        }
    restore
}

postclose `posth'

use "`outdir'/first_birth_by_wealth_summary_v1.dta", clear
sort outcome wealth_q
export delimited using "`outdir'/first_birth_by_wealth_summary_v1.csv", replace

di as result "First-birth by wealth event studies complete."
di as result "Summary CSV: `outdir'/first_birth_by_wealth_summary_v1.csv"

log close _all
