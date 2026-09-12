clear all
set more off
cd "/Users/tommasodesanto/Desktop/Projects/Fertility" 

*************************************************************************
* 1. DEFINING GLOBAL MACROS
*************************************************************************
global path "/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs"
global tables "$path/Tables/tests"
global graphs "$path/Graphs/tests"
global weight "IW"  // Original weight variable

* Define all possible outcomes for later expansion
global all_outcomes log_INCFAMR log_EARNINDR log_INCFAMFR log_EARNINDFR log_EARNINDRRC log_EARNINDFRRC
global other_outcomes mv_n mv_s rooms own

* Start with just EARNINDR for initial testing
global outcomes log_e  

* LATER: To run for more/all outcomes, uncomment this line:
* global outcomes $all_outcomes


*Other globals
global ref_sa 2 // this means that we normalize K=-2 to zero
global groups mvr n_mvr
global ei_vars f_c_y //f_m_y


*************************************************************************
* 2. DATA PREPARATION - OPTIMIZE BY KEEPING ONLY NEEDED VARIABLES
*************************************************************************
use "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta", clear


* Keep only essential variables for our analysis
**first some cleaning
	*  recode owning: equal 0 if rent
replace HOMEOWN=0 if HOMEOWN==2 
* random code if neither
replace HOMEOWN=10 if HOMEOWN==3 
	**for this exercise, drop people for which we have ambiguity
	replace HOMEOWN=. if HOMEOWN==10
	
	**drop dead
	drop if year>DEATHYEAR


 gen moved_for_size=1 if WHYMOVED1_==3 & WHYMOVED1_!=. & year>1974
 replace moved_for_size=0  if WHYMOVED1_!=3 & WHYMOVED1_!=.
  gen moved_for_neigh=1 if WHYMOVED1_==6 & WHYMOVED1_!=. & year>1974
 replace moved_for_neigh=0  if WHYMOVED1_!=6 & WHYMOVED1_!=.

keep ID year AGEREP EDUYEAR SEX RELCHI1BYEAR MOVEDFREF_ DEATHYEAR moved_for_size  moved_for_neigh  ACTUALROOMS_ HOMEOWN ///
     EARNINDR EDUYEAR ${weight}
     
* Add other income variables only if needed
foreach var in INCFAMR INCFAMFR EARNINDFR EARNINDRRC EARNINDFRRC {
    * Check if we need this variable based on outcomes selected
    if strpos("$outcomes", "`var'") {
        qui merge 1:1 ID year using "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta", ///
            keepusing(`var') keep(master match) nogen
    }
}


* Basic cleaning and sample selection
drop if AGEREP < 18
drop if year > DEATHYEAR

* Set up panel structure
xtset ID year

* Clean mobility variables
rename MOVEDFREF_ movedthisyear

replace movedthisyear = . if movedthisyear == 8 | movedthisyear == 9 
replace movedthisyear = 0 if movedthisyear == 5

xtset ID year
gen change_own= (HOMEOWN==1 & L.HOMEOWN==0)
gen moved_to_own=(movedthisyear==1 & change_own==1)

* Generate gender variable
gen female = (SEX == 2)
gen male = 1-female
drop SEX



* Generate logged income variables (only for variables we need)
foreach var in EARNINDR {
    gen log_`var' = ln(`var')
    replace log_`var' = 0 if `var' == 0
    label var log_`var' "Log `var'"
}
* Add logging for other income variables if needed
foreach var in INCFAMR INCFAMFR EARNINDFR EARNINDRRC EARNINDFRRC {
    if strpos("$outcomes", "log_`var'") {
        gen log_`var' = ln(`var')
        replace log_`var' = 0 if `var' == 0
        label var log_`var' "Log `var'"
    }
}

* Create compact fertility variables
gen hadfirstchild = 0
replace hadfirstchild = 1 if year == RELCHI1BYEAR
bysort ID: egen onechild = max(hadfirstchild)

* Drop individuals with children before sample entry
bysort ID: egen year_entry = min(year)
gen hadchildbeforesample = 0
replace hadchildbeforesample = 1 if RELCHI1BYEAR < year_entry
drop if hadchildbeforesample == 1
drop hadchildbeforesample year_entry

*************************************************************************
* 3. CREATING WEIGHTS & EVENT STRUCTURE (OPTIMIZED)
*************************************************************************
* Convert analytical weights to frequency weights
sum ${weight}, detail
gen wt_norm = ${weight}/r(mean)
gen fw_wt = round(wt_norm*10)
replace fw_wt = 1 if fw_wt < 1 & !missing(fw_wt)
drop wt_norm ${weight}  // Drop original weight to save memory

* Set up event study structure
gen move_event = movedthisyear
gen move_year = year if move_event == 1
bysort ID: egen first_move_year = min(move_year)
drop move_year move_event  // Drop intermediate variables

*************************************************************************
* 3b. Creating event variables
*************************************************************************

* Generate time relative to move and treatment indicator
gen time_to_move = year - first_move_year
bysort ID: egen ever_moved = max(movedthisyear)


* Post-move fertility indicator
gen birth_post_move = (RELCHI1BYEAR >= first_move_year & RELCHI1BYEAR <= first_move_year + 3)
bysort ID: egen had_child_post_move = max(birth_post_move)

rename RELCHI1BYEAR first_child_year

gen time_to_event = year - first_child_year

gen post_birth_window = (time_to_event > -2 & time_to_event <= 1)
gen moved_post_birth = movedthisyear & post_birth_window


by ID: egen mvr=max(moved_post_birth)
gen n_mvr=1-mvr



rename  moved_for_size mv_s
rename moved_for_neigh mv_n
rename ACTUALROOMS_ rooms 
rename HOMEOWN own
* Compress data to save memory
compress

/*
*************************************************************************
* 4. CALLAWAY & SANT'ANNA IMPLEMENTATION WITH REDUCED COMPLEXITY
*************************************************************************
* Install required packages if needed
foreach package in csdid drdid event_plot {
    cap which `package'
    if _rc {
        ssc install `package', replace
    }
}

* Start with just EARNINDR
foreach outcome in $outcomes {
    local shortname = subinstr("`outcome'", "log_", "", .)
    
    * Full sample with reduced bootstraps
    di "Estimating `outcome' for full sample"
    csdid `outcome', ivar(ID) time(year) gvar(first_move_year) ///
          agg(event) method(dripw) notyet wboot(reps(99)) ///
          [fweight=fw_wt] rseed(1234)
    
    csdid_stats
    estat event
    matrix b = r(table)
    matrix V = r(V)
    
    * Save results
    preserve
        clear
        svmat b
        svmat V
        
        * Add metadata
        gen outcome = "`shortname'"
        gen sample = "all"
        gen event_time = _n - 11
        
        * Generate confidence intervals
        gen estimate = b1
        gen stderr = sqrt(V1)
        gen lb = estimate - 1.96*stderr
        gen ub = estimate + 1.96*stderr
        
        * Keep only essential data
        keep outcome sample event_time estimate stderr lb ub
        save "$csresults/cs_`shortname'_all.dta", replace
    restore
    
    * Gender sub-analysis
    di "Estimating `outcome' for females"
    cap csdid `outcome' if female == 1, ivar(ID) time(year) ///
          gvar(first_move_year) agg(event) method(dripw) ///
          notyet wboot(reps(99)) [fweight=fw_wt] rseed(1234)
          
    if _rc == 0 {
        csdid_stats
        estat event
        matrix b = r(table)
        matrix V = r(V)
        
        preserve
            clear
            svmat b
            svmat V
            
            gen outcome = "`shortname'"
            gen sample = "female"
            gen event_time = _n - 11
            
            gen estimate = b1
            gen stderr = sqrt(V1)
            gen lb = estimate - 1.96*stderr
            gen ub = estimate + 1.96*stderr
            
            keep outcome sample event_time estimate stderr lb ub
            save "$csresults/cs_`shortname'_female.dta", replace
        restore
    }
    else {
        di "WARNING: Female subsample estimation failed, continuing with next step"
    }
    
    di "Estimating `outcome' for males"
    cap csdid `outcome' if female == 0, ivar(ID) time(year) ///
          gvar(first_move_year) agg(event) method(dripw) ///
          notyet wboot(reps(99)) [fweight=fw_wt] rseed(1234)
          
    if _rc == 0 {
        csdid_stats
        estat event
        matrix b = r(table)
        matrix V = r(V)
        
        preserve
            clear
            svmat b
            svmat V
            
            gen outcome = "`shortname'"
            gen sample = "male"
            gen event_time = _n - 11
            
            gen estimate = b1
            gen stderr = sqrt(V1)
            gen lb = estimate - 1.96*stderr
            gen ub = estimate + 1.96*stderr
            
            keep outcome sample event_time estimate stderr lb ub
            save "$csresults/cs_`shortname'_male.dta", replace
        restore
    }
    else {
        di "WARNING: Male subsample estimation failed, continuing with next step"
    }
}
*/

*************************************************************************
* 4b. SUN & ABRAHAM IMPLEMENTATION 
*************************************************************************
rename log_EARNINDR log_e
rename first_move_year f_m_y
rename first_child_year f_c_y

program do_estimations
    syntax, outcome_vars(str) ei(str) name(str) [ controls(str)] 
	foreach var in `outcome_vars' {	
	estimations_sa, y(`var') i(ID) t(year)  ///
			 ei(`ei')  controls(`controls')  naming(`name')
			  	 
	}

end		


	program estimations_sa,
		syntax, y(varlist) i(varlist) t(varlist) ei(varlist) naming(str) [controls(str)] 
		preserve
		
		capture gen K = `t'-`ei' 		
		sum `ei'
		capture gen lastcohort = `ei'==r(max) // dummy for the latest- or never-treated cohort
		cap drop L*event F*event
		forvalues l = 0/10 {
			cap gen L`l'event = K==`l'
		}
		cap gen L11event = (K>10 & K!=.)	
		forvalues l = 1/5 {
			cap gen F`l'event = K==-`l'
		}
		cap gen F7event = (K<-6 & K!=.)		
		drop F${ref_sa}event // normalize K=-${ref_sa} (and also K=-6) to zero
				
		eventstudyinteract `y' L*event F*event , ///
				vce(cluster `i') absorb(`i' `t') cohort(`ei') control_cohort(lastcohort) ///
				 covariates(i.AGEREP i.EDUYEAR) 
			
			summ `y' if e(sample) & K==-${ref_sa}
			save_estimates_sa, name(`y'_`ei'_`naming') prebirthmean(`r(mean)')
		
		
		restore		
	end
	
		program save_estimates_sa
			syntax , name(str)	prebirthmean(real) 
			
			matrix b = e(b_iw)
			matrix var = e(V_iw)
			matrix vardiag=vecdiag(var)
			matrix combine = b \ vardiag
			matrix rownames combine = b variance
			matrix `name' = combine'
			matrix drop combine	
			preserve
			clear
			svmat2 `name', names(col) rnames(coeff)
			gen prebirth_mean = `prebirthmean'
			save "$tables/`name'_temp.dta", replace	
			restore			
		end	

program construct_estimates_datasets
    syntax, outcome_vars(str) name(str) [controls(str)]
	
		foreach var in `outcome_vars' {								
			construct_dataset_sa, y(`var') naming(`name') ref_lead($ref_sa)
		}
end			

	program construct_dataset_sa
		syntax, y(str) ref_lead(int) naming(str)

		use "$tables/`y'_`naming'_temp.dta", clear
		
		*drop *11
		*drop *12
		
		*rename *1 b 
		*rename *2 variance 
		*keep b variance coef prebirth_mean
			gen se = sqrt(variance)
		
			drop variance			
			
			replace b = . if b == 0 & se == 0
			replace se = . if b == . & se == 0
			
			drop if coeff == "_cons"	
			gen relative_time = subinstr(coeff,"event","",.)
			replace relative_time = subinstr(relative_time,"L","",.)
			replace relative_time = subinstr(relative_time,"F","-",.)	
			replace relative_time = subinstr(relative_time,"o.","",.)				
			destring relative_time, replace

			* set the coeff and CI to zero in the relative year
			assert relative_time != -`ref_lead' 
			assert relative_time != .			
			set obs `=_N+1'
			replace relative_time = -`ref_lead' if relative_time == .
			replace b = 0 if relative_time == -`ref_lead'
			replace se = 0 if relative_time == -`ref_lead'
			
			gen ci_lo = b -1.96*se
			gen ci_hi = b +1.96*se				
			
			sort relative_time
			order coeff relative_time b se ci_lo ci_hi
			
			save "$tables/`y'_`naming'_estimates.dta", replace
	
	end
	

	program graph_eventstudies
	syntax, outcome_vars(str) name(str) [gen(str)]
	
	global color2005 dkgreen
	
	foreach y in `outcome_vars' {		
		
		use "$tables/`y'_`name'_estimates.dta", clear	
	    drop if b == . 
		local xlabel 
		if coeff == "F4event" local xlabel  xlabel(-4 0 7)
	
			
			twoway (rarea ci_lo ci_hi relative_time, fcolor(${color2005}%10) lcolor(${color2005}%10)) ///
				   (connected b relative_time, mcolor(${color2005}) msymbol(circle) lcolor(${color2005}) lwidth(medthick) lpattern(dash))  , ///
				   ytitle("${`y'title}") `titlegap' xtitle(Years relative to childbirth) ///
				   xline(-1.5, lwidth(thin) lpattern(dash) lcolor(gray)) `xlabel' ///
				   yline(0, lwidth(vthin) lpattern(solid) lcolor(gray)) ///			   
				   `graphregion'
				   
			graph export "$graphs/`y'_`name'.png", replace 
		
	}
end	

program graph_both_groups
	syntax, outcome_vars(str) [name(str)] [controls(str)]  w(str) m(str)
	
	local naming `name'
	
	global color_women   pink
	global color_men  ebblue
	
	foreach y in `outcome_vars' {	
	
		use "$tables/`y'_`naming'_`w'_estimates.dta", clear
		renvars b* ci* prebirth_mean*, postfix(_f)
		merge 1:1 relative_time using "$tables/`y'_`naming'_`m'_estimates.dta"
		drop if b == . & b`p'_f == .			
		local xlabel 
		if coeff == "F4event" local xlabel  xlabel(-4 0 7)
		
		
			local rounding 0.01
			if inlist("`y'", "leis", "leisi") local rounding 1
			local pbm ⁠ =round( ⁠=prebirth_mean`p'',`rounding')'	
			local pbf ⁠ =round( ⁠=prebirth_mean`p'_f',`rounding')'	
			
			twoway (rarea ci_lo ci_hi relative_time, fcolor(${color_men}%10) lcolor(${color_men}%10)) ///
				   (connected b relative_time, mcolor(${color_men}) msymbol(square) lcolor(${color_men}) lwidth(medthick)) ///
				   (rarea ci_lo_f ci_hi_f relative_time, fcolor(${color_women}%10) lcolor(${color_women}%10)) ///
				   (connected b_f relative_time, mcolor(${color_women}) msymbol(circle) lcolor(${color_women}) lwidth(medthick)), ///
				   ytitle("${`y'title}") `titlegap' xtitle(Years relative to childbirth) ///
				   xline(-1.5, lwidth(thin) lpattern(dash) lcolor(gray)) `xlabel' ///
				   yline(0, lwidth(vthin) lpattern(solid) lcolor(gray)) ///			   
				   legend(order(2 "Non-Movers" 4 "Movers" ) rows(1) region(fcolor(none))) ///
				    note("Pre-birth Mean Non-Movers: `pbm'; Pre-birth Mean Movers: `pbf'.")
				   
			graph export "$graphs/`y'`w'_vs`m'_`naming'.png", replace 
		}	

	
end


foreach ei in $ei_vars {
	preserve
			do_estimations, outcome_vars($other_outcomes)  name(all) ei(`ei')
		construct_estimates_datasets, outcome_vars($other_outcomes) name(`ei'_all)
		graph_eventstudies, outcome_vars($other_outcomes) name(`ei'_all)	
		restore
	foreach group in $groups {
	    preserve
		keep if `group' == 1
		do_estimations, outcome_vars($outcomes)  name(`group') ei(`ei')
		construct_estimates_datasets, outcome_vars($outcomes) name(`ei'_`group')
		graph_eventstudies, outcome_vars($outcomes) 	name(`ei'_`group')			

		restore
	}
		do_estimations, outcome_vars($outcomes)  name(all) ei(`ei')
		construct_estimates_datasets, outcome_vars($outcomes) name(`ei'_all)
		graph_eventstudies, outcome_vars($outcomes) name(`ei'_all)	
				graph_both_groups, outcome_vars($outcomes) name(`ei') w(mvr) m(n_mvr)
}




	
/*
*************************************************************************
* 5. CREATE EVENT STUDY GRAPHS FROM SAVED RESULTS
*************************************************************************
* Define outcome titles
local title_INCFAMR "Family Income (Real)"
local title_EARNINDR "Individual Earnings (Real)"
local title_INCFAMFR "Family Income (Family Size Adjusted)"
local title_EARNINDFR "Individual Earnings (Family Size Adjusted)"
local title_EARNINDRRC "Couple Combined Earnings (Real)"
local title_EARNINDFRRC "Couple Combined Earnings (Family Size Adjusted)"

* Process each outcome that has been analyzed
foreach outcome in $outcomes {
    local shortname = subinstr("`outcome'", "log_", "", .)
    local outcome_title = "`title_`shortname''"
    
    * Check if results exist before attempting to graph
    cap confirm file "$csresults/cs_`shortname'_all.dta"
    if _rc == 0 {
        * Create overall effect graph
        use "$csresults/cs_`shortname'_all.dta", clear
        
        twoway (connected estimate event_time, lcolor(blue) mcolor(blue)) ///
               (rcap lb ub event_time, lcolor(blue*.5)), ///
               xline(0, lpattern(dash) lcolor(red)) yline(0, lpattern(dash)) ///
               ylabel(, angle(horizontal)) ///
               xlabel(-10(2)10, labsize(small)) ///
               xtitle("Years Relative to Move", size(small)) ///
               ytitle("Effect on `outcome_title'", size(small)) ///
               title("Mobility Impact on `outcome_title'", size(medium)) ///
               subtitle("Callaway & Sant'Anna Estimator", size(small)) ///
               legend(off) scheme(s1mono)
        graph export "$graphs/cs_`shortname'_all_event_study.pdf", replace
        
        * Check if gender-specific results exist
        cap confirm file "$csresults/cs_`shortname'_female.dta"
        cap confirm file "$csresults/cs_`shortname'_male.dta"
        if _rc == 0 {
            * Create gender comparison if both exist
            use "$csresults/cs_`shortname'_female.dta", clear
            rename estimate f_est
            rename lb f_lb
            rename ub f_ub
            keep event_time f_est f_lb f_ub
            
            merge 1:1 event_time using "$csresults/cs_`shortname'_male.dta", ///
                  keepusing(estimate lb ub) nogen
            rename estimate m_est
            rename lb m_lb
            rename ub m_ub
            
            twoway (connected f_est event_time, lcolor(red) mcolor(red) lpattern(solid)) ///
                   (rcap f_lb f_ub event_time, lcolor(red*.5)) ///
                   (connected m_est event_time, lcolor(blue) mcolor(blue) lpattern(dash)) ///
                   (rcap m_lb m_ub event_time, lcolor(blue*.5)), ///
                   xline(0, lpattern(dash) lcolor(black)) yline(0, lpattern(dash)) ///
                   ylabel(, angle(horizontal)) ///
                   xlabel(-10(2)10, labsize(small)) ///
                   xtitle("Years Relative to Move", size(small)) ///
                   ytitle("Effect on `outcome_title'", size(small)) ///
                   title("Gender Differences in Mobility Effects", size(medium)) ///
                   subtitle("`outcome_title' (Callaway & Sant'Anna)", size(small)) ///
                   legend(order(1 "Females" 3 "Males") pos(6) rows(1)) ///
                   scheme(s1mono)
            graph export "$graphs/cs_`shortname'_gender_event_study.pdf", replace
        }
    }
}

*************************************************************************
* 6. SIMPLE SUMMARY STATISTICS
*************************************************************************
use "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta", clear

* Keep only essential variables
keep ID year AGEREP MOVEDFREF_ DEATHYEAR ${weight}

* Add variables needed for analysis
foreach var in $outcomes {
    local shortname = subinstr("`var'", "log_", "", .)
    cap merge 1:1 ID year using "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta", ///
        keepusing(`shortname') keep(master match) nogen
}

* Basic cleaning
drop if AGEREP < 18
drop if year > DEATHYEAR
rename MOVEDFREF_ movedthisyear
replace movedthisyear = . if movedthisyear == 8 | movedthisyear == 9 
replace movedthisyear = 0 if movedthisyear == 5

* Simplified event structure
gen move_year = year if movedthisyear == 1
bysort ID: egen first_move_year = min(move_year)
gen time_to_move = year - first_move_year
bysort ID: egen ever_moved = max(movedthisyear)

* Create pre-post windows
gen pre_move = (time_to_move >= -3 & time_to_move < 0)
gen post_move = (time_to_move >= 0 & time_to_move < 3)

* Process each outcome
foreach var in $outcomes {
    local shortname = subinstr("`var'", "log_", "", .)
    
    * Calculate pre and post averages
    gen pre_`shortname' = `shortname' if pre_move
    gen post_`shortname' = `shortname' if post_move
}

* Collapse to individual level
preserve
    collapse (mean) pre_* post_* (first) ever_moved [aw=${weight}], by(ID)
    
    * Calculate percent changes
    foreach var in $outcomes {
        local shortname = subinstr("`var'", "log_", "", .)
        gen pct_`shortname' = ((post_`shortname' - pre_`shortname')/pre_`shortname')*100 ///
                             if !missing(pre_`shortname') & !missing(post_`shortname')
    }
    
    * Create summary statistics table for each outcome
    foreach var in $outcomes {
        local shortname = subinstr("`var'", "log_", "", .)
        local title = "`title_`shortname''"
        
        estpost tabstat pre_`shortname' post_`shortname' pct_`shortname', ///
                by(ever_moved) statistics(mean sd n) columns(statistics)
        esttab . using "$tables/`shortname'_by_mobility.tex", replace ///
             cells("mean(fmt(%9.0fc)) sd(fmt(%9.0fc)) count") ///
             noobs title("`title' Before and After Mobility") ///
             label mtitles("Non-Movers" "Movers")
    }
restore
*/

*************************************************************************
* END OF CODE
*************************************************************************
