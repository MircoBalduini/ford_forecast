capture program drop _all
program define ford
*
* selection of the variable Y and take logs (if y>0)
g y = `1'
g x = y
*
* defines the forecast horizon
local H = `2' 
* sets primary S(easonality) on the basis of dataset periodicity  
****************************MUST BE EXTENDED TO OTHER PERIODICITIES!
local S = 1
qui tsset 
if r(unit) == "monthly" {
	local S = 12
	gen month = month(dofm(t))
}
else if r(unit) == "quarterly" {
	local S = 4
	gen quarter = quarter(dofq(t))
}
*
* counts TO(tal) obs and defines a max sample of [maxyears] years
* to be used in easch forecast iteration
qui summ y
local TO = r(N)
local maxyears = 40 
local screensample "(uncut data source)"
if `TO'> `maxyears'*`S' {
	local cutsample = `TO'-`maxyears'*`S'
	drop if _n <= `cutsample'
	local screensample "(data source cut to `maxyears' years)"
}
*
* logs or levels? (use the whole sample to prevent log problems
qui summ y
local min = r(min)
local screenlog "untrasformed levels [x=y]"
if `min'>0 {
	drop x
	g y_begin = y
	g x = log(y)
	replace y = x
	local screenlog "transformed in logs [x=log(y)]"
}
*
* T counts the training obs 
qui tsset
local tmx = r(tmaxs)
local tmn = r(tmins)
summ x
local T = r(N)
*
* compare variances to choose between x or d.x 
qui summ x , d 
local varlev = r(Var)
qui summ d.x , d 
local vardif = r(Var)
local Irho = 0
local screeniro "is modelled in levels"
if `vardif' <= 1.2*`varlev' {
	g ytemp = d.y
	drop y
	rename ytemp y
	local screeniro "is modelled (by Delta) in differences"
	local Irho = 1
	qui summ y
***********************************************	local T = r(N)
}
*
	local IA = 0
	local screenia "not present"
	local IR = 0
	local screenir "not present"

if `S'>1 & `T'>3*`S'+`Irho'{
*	ANOVA-test for additive (IA) and autoregressive (IR) seasonality
	qui reg y b`S'.month  
 	qui testparm b`S'.month
	if r(p) < 0.1 {
		local IA = 1
		local screenia "present"
	}
	qui reg y l`S'.y
	local df = e(N)-1
	local ti = abs(_b[l`S'.y]/_se[l`S'.y])
	if t(`df', `ti') > .95 {
		local IR = 1
		local screenir "present"
	}
*	
}
di ""
di ""
di "PRELIMINARY STAGE OUTCOMES "
di " all data from `tmn' to `tmx'"   
di " primay seasonality, S = `S'"
di " variable to be forecast y = `1'"
di " data `screenlog'"
di " Irho = `Irho' [x `screeniro']"
di " available obs for modelling x = `T' `screensample'"
di " additive seasonality `screenia', i.e. IA=`IA'"
di " autoregressive seas. `screenir', i.e. IR=`IR'"
di ""
di ""
* save locals in scalars to use them outside the program
scalar primary_S = `S'
scalar I_rho = `Irho'
scalar estimation_T = `T'
scalar I_A= `IA'
scalar I_R= `IR'
scalar sample_H = `H'
*
********************************************************************************
***************************** DELTA ESTIMATION *********************************
********************************************************************************
*
if `Irho' == 1 {
*
*********************** Irho == 1 (variable in differences)*********************
*
	if `IA' == 1 {
		*
		******* (seasonality) seasonal table and 3x5 means
		*
		gen year =year(dofm(t))
		sort year month
		save storico, replace
		keep year month y
		reshape wide y, i(year) j(month)
		tsset year
		tsappend, add(3)
		foreach j of num 1/`S' {
			carryforward y`j', replace
			g Y`j' = (1/15)*f3.y`j' + (2/15)*f2.y`j' + (3/15)*f1.y`j' +  ///
			(3/15)*y`j' + (3/15)*l1.y`j' + (2/15)*l2.y`j' + (1/15)*l3.y`j'
		*	the denominator of the ratio above is 15 because it is a MA 3*5 		
		*	carryforward is used not to lose last obs
		}
		*
		******** merge new and old df
		*
		keep year Y* 
		reshape long Y, i(year) j(month)
************************************************                                            
*		drop if _n > `T' 
		sort year month
		save mediemobili, replace
		merge 1:1 year month using storico
		keep if _m == 3
		*
		gen rev_trend = _N-_n+1
		gen rev_year = ceil(rev_trend/`S')
		sort rev_year
		*
		*
		by rev_year: egen mmss = mean(Y)
		g s = Y - mmss
		by rev_year: egen z = mean(y)
		*
		*
		tsset t
		*
		gen y_temp = x
		*
		*
		****** generate dr
		*
		local dr = (z[_N]+z[_N-`S']+z[_N-2*`S']+z[_N-3*`S']+ z[_N-4*`S'] + ///
		z[_N-5*`S'])/6
		*
		****** generate dm and delta_s
		*
		gen Dy = y_temp-l`S'.y_temp
		sum Dy , meanonly
		local delta_s = r(mean)/`S'
		*
		sum z if Dy != . , meanonly
		local dm = r(mean)
		*
		*
	}
	*
    ******************************** case IA == 0 ******************************
	*
	else {
		tsset t
		gen y_temp = x
		*
		rename y z
		gen year =year(dofm(t))
		gen rev_trend = _N-_n+1
		*
		******** generate dr
		*
		local dr = (z[_N]+z[_N-1]+z[_N-2]+z[_N-3]+z[_N-4]+z[_N-5])/6
		*
		******** generate dm and delta_s
		*
		gen Dy = y_temp-l`S'.y_temp
		sum Dy , meanonly
		local delta_s = r(mean)/`S'
		*
		sum z if Dy != . , meanonly
		local dm = r(mean)
		*
		g s = 0
	}
	*
	******* compute d1, d2 (cases Irho == 1 together) 
	*
	save storico2, replace
	*
	if `T' <= 6 {
		local d1 = `dm'
		local d2 = `dm'
	}
	else {
		*	
		******** compute d1
		*
		qui summ z
		local max = r(max)
		local min = r(min)
		*
		if abs(`min') < abs(`max') {
			sort z
			drop if z >= z[_N]
			qui summ z
			local d1 = r(mean)
		}
		else {
			gsort -z
			drop if z <= z[_N]
			qui summ z
			local d1 = r(mean)
		}
		*
		******** compute d2
		*
		foreach i in 1 2 {
			qui summ z
			local max = r(max)
			local min = r(min)
			if abs(`min') < abs(`max') {
				sort z
				drop if z >= z[_N]
			}
			else {
				gsort -z
				drop if z <= z[_N]
			}
		}
		*
		qui summ z
		local d2 = r(mean)
		*
		*
		drop _all
		use storico2, replace
		tsset t month
		*
	}
	********************* compute amin's (cases Irho == 1 together) ************
	*
	*
	if `T' > 1+2*`S' {
		if `dm'*`delta_s' <= 0 {
			local dm = 0
		}
		else if abs(`dm') <= abs(`delta_s') {
			local dm = `dm'
		}	
		else {
			local dm = `delta_s'
		}
	}
	*
	*
	if `dr'*`dm' <= 0 {
		local dr_star = 0
	}
	else if abs(`dr') <= abs(`dm') {
		local dr_star = `dr'
	}	
	else {
		local dr_star = `dm'
	}
	*
	******* compute amin for forecast T+1
	*
	if `dr_star'*`d1' <= 0 {
		local amin1 = 0
	}
	else if abs(`dr_star') <= abs(`d1') {
		local amin1 = `dr_star'
	}	
	else {
		local amin1 = `d1'
	}
	*
	****** compute amin for forecast T+h
	*
	if `dr_star'*`d2' <= 0 {
		local aminh = 0
	}
	else if abs(`dr_star') <= abs(`d2') {
		local aminh = `dr_star'
	}	
	else {
		local aminh = `d2'
	}
	***************************** compute forecast ***************************** 
	*
	******* forecast begins in "start"
	tsset t
	qui summ t
	local start = r(max)+1
	************************************ GOLI:
	local obs = r(max)-r(min)+1
	
	tsappend, add(`H') 
	replace month = month(dofm(t))
	replace year =year(dofm(t))
	*
	gen trend = _n
	*
	gen exog_esse = 0
	local counter = 1
	while `counter' <= `H'{
		forval i = 1/`S' {
	************************************
			replace exog_esse = s[`obs'-`S'+`i'] if trend == `obs'+`counter'
* 			replace exog_esse = s[`T'-`S'+`i'] if trend == `T'+`counter'
			local counter = `counter'+1
		}
	}	
	*
	*********** compute exog_amin 
	*
	g exog_amin=0
	replace exog_amin = `amin1' if trend ==`T'+1 
	replace exog_amin = `aminh' if trend >`T'+1 
	*
	forecast create pippo, replace
	forecast identity x = l.x + exog_amin + exog_esse
	forecast exogenous exog_amin exog_esse
	qui forecast solve , prefix(delta_) begin(`start') periods(`H')
	*
	*
}
*
******************************* cases Irho == 0 ********************************
*
else { 
	if `IA' == 1 {
		*
		************************ case IA == 1 **********************************
		*
		gen year =year(dofm(t))
		sort year month
		save storico, replace
		keep year month y
		reshape wide y, i(year) j(month)
		tsset year
		tsappend, add(5)
		foreach j of num 1/`S' {
			carryforward y`j', replace
			g Y`j' = (1/35)*f5.y`j' + (2/35)*f4.y`j' + (3/35)*f3.y`j' + ///
			(4/35)*f2.y`j' + (5/35)*f1.y`j' +  (5/35)*y`j' + (5/35)*l1.y`j' ///
			+ (4/35)*l2.y`j' + (3/35)*l3.y`j' + (2/35)*l4.y`j' + (1/35)*l5.y`j'
		*		the denominator of the ratio above is 35 because it is a MA 7*5 		
		*		carryforward is used not to lose last obs
		}
		*
		******** merge new and old df 
		*
		keep year Y* 
		reshape long Y, i(year) j(month)
		sort year month
		save mediemobili, replace
		merge 1:1 year month using storico
		keep if _m == 3
		*
		*
		gen rev_trend = _N-_n+1
		gen rev_year = ceil(rev_trend/`S')
		sort rev_year
		*
		********** compute s and z
		*
		by rev_year: egen mmss = mean(Y)
		g s = Y - mmss
		by rev_year: egen z = mean(y)
		*
		sort t
		local S_tilde = 1
		gen trend = _n
		*		
	}
	else { /***************** case IA = 0, Irho = 0 ***************************/
		tsset t
		gen year = year(dofm(t))
		rename y z
		local S_tilde = max(2,`S')
		gen trend = _n
		gen rev_trend = _N-_n+1
	}	
	*
	************************* compute FORECAST ********************************* 
	*
	tsset t
	gen delta_x = x
	local begin = `T'+1
	gen s_hat = s
	*
	*
	if `S' >= `T' {
		forval k = 1/`H' { /* exception if t_series is short */
			qui sum delta_x
			local m_z = r(mean)
			tsappend, add(1)
			replace trend = _n
			replace delta_x = `m_z' if trend == `T'+`k'
		}
	}
	else {
		*
		qui sum z if rev_trend <= `S_tilde'
		local m_S_tilde = r(mean)
		*
		******** forecast 1st step ahead 
		*
		tsappend, add(1)
		replace trend = _n
		replace rev_trend = _N-_n+1
		if `IA' == 1 {
			replace s_hat = s[`T'-`S'+`IA'] if trend == `begin'
		}
		replace delta_x = `m_S_tilde'+s_hat if trend == `begin'
		*
		********* forecast H steps ahead 
		*
		local counter = 2
		while `counter' <= `H'{
			*
			************ compute step-by-step m(r)
			*
			if `IA' == 0 {
				qui sum delta_x if rev_trend <= 6*`S_tilde'
				local m_6S_tilde = r(mean)
				qui sum delta_x if rev_trend <= `S_tilde'
				local m_S_tilde = r(mean)
				tsappend, add(1)
				replace trend = _n
				replace rev_trend = _N-_n+1
			}
			else {
				replace rev_year = ceil(rev_trend/`S')
				sort rev_year
				by rev_year: egen k`counter' = mean(delta_x)
				sort t
				tsappend, add(1) 
				replace trend = _n
				replace rev_trend = _N-_n+1
				qui sum k`counter' if rev_year <= 6*`S_tilde'
				local m_6S_tilde = r(mean)
				qui sum k`counter' if rev_year <= `S_tilde'
				local m_S_tilde = r(mean)
			}	
			*
			************** compute seasonal component
			*
			local begin = `begin'+1
			if `IA' == 0 {
				replace s_hat = 0
			}
			else {
				replace s_hat = s_hat[`T'-`S'+`counter'] if trend == `begin'
			}
			replace delta_x =((`m_S_tilde'+`m_6S_tilde')/2)+s_hat if trend == `begin'
			drop k`counter'
			local counter = `counter'+1
		}
		*
		replace year =year(dofm(t))
		replace month = month(dofm(t))
		*
	}
}
*
********************************************************************************
****************************** RHO ESTIMATION **********************************
********************************************************************************
*
local Ir = `Irho' 
local Idelta = 0
local Itau = 0
*
di `IR'
di `T'
di `Ir'
*
if `IA' == 1 {
	*
	****************** case IA == 1
	*
	* compute centered seasonal dummies
	*
	forval i = 2/12 {
		gen seas`i' = (month == `i')
		replace seas`i' = seas`i'-1/`S'
	}
*
	if `Ir' == 1 & `IR' == 1 & `Idelta' == 0 & `Itau' == 0 { 
		*
		******************* case Ir == 1 & IR == 1 *******************
		*
		di ""
		di "******* case IR =`IR' **********"
		di "******* case Ir =`Ir' **********"
		di "******* case IA =`IA' **********"
		di ""
		qui reg x l.x l`S'.x seas*
		est store rho_6
		if _b[l.x]> .5 & _b[l.x] + 2*_se[l.x] > .9 {
			local Idelta = 1
			constraint 1 l.x = 1
			qui cnsreg x l.x l`S'.x seas*, constraint(1)
			est store rho_2_1 
		}	
		else if _b[l.x]< 0 {
			local Ir = 0
			qui reg x seas*
			est store rho_2_2
			if `Idelta' == 0 & `T'-`IR'-`Ir'-`IA' > 10  {
				predict eps, resid
				g cum_eps = sum(eps) if eps!=.
				qui ttest cum_eps == 0
				local pval = r(p)
				if `pval' < 0.01 {						
					gen trend_a = _n/`S' 
					local Itau = 1
					qui reg x l.x l`S'.x trend_a seas*
					if _b[l.x]< -0.5 {
						local Itau = 0
						est restore rho_6
					}
					else {
						est store rho_2_2_3
					}
				}
				else {
					est restore rho_2_2
				}
			}
			else restore rho_2_2
		}
		else {
			if `Idelta' == 0 & `T'-`IR'-`Ir'-`IA' > 10  {
				predict eps, resid
				g cum_eps = sum(eps) if eps!=.
				qui ttest cum_eps == 0
				local pval = r(p)
				if `pval' < 0.01 {						
					gen trend_a = _n/`S' 
					local Itau = 1
					qui reg x l.x l`S'.x trend_a seas*
					if _b[l.x]< -0.5 {
						est restore rho_6
					}
					else {
						est store rho_3
					}
				}
				else {
					est restore rho_6
				}
			}
			else {
				est restore rho_6
			}
		}
	}
	else if `Ir' == 1 & `IR' == 0 & `Idelta' == 0 & `Itau' == 0 {
		*
		******************* case Ir == 1 & IR == 0 *******************
		*
		di ""
		di "******* case IR =`IR' **********"
		di "******* case Ir =`Ir' **********"
		di "******* case IA =`IA' **********"
		di ""
		qui reg x l.x seas*
		est store rho_6
		if _b[l.x]> .5 & _b[l.x] + 2*_se[l.x] > .9 {
			local Idelta = 1
			constraint 1 l.x = 1
			qui cnsreg x l.x seas*, constraint(1)
			est store rho_2_1 
		}	
		else if _b[l.x]< 0 {
			local Ir = 0
			qui reg x seas*
			est store rho_2_2
			if `Idelta' == 0 & `T'-`Ir'-`IA' > 10  {
				predict eps, resid
				qui ttest eps == 0
				local pval = r(p)
				if `pval' < 0.01 {						
					gen trend_a = _n/`S' 
					local Itau = 1
					qui reg x l.x trend_a seas*
					if _b[l.x]< -0.5 {
						local Itau = 0
						est restore rho_6
					}
					else {
						est store rho_2_2_3
					}
				}
				else {
					est restore rho_2_2
				}
			}
			else {
				est restore rho_2_2
			}
		}
		else {
			if `Idelta' == 0 & `T'-`IR'-`Ir'-`IA' > 10  {
				predict eps, resid
				qui ttest eps == 0
				local pval = r(p)
				if `pval' < 0.01 {						
					gen trend_a = _n/`S' 
					local Itau = 1
					qui reg x l.x trend_a seas*
					if _b[l.x]< -0.5 {
						est restore rho_6
					}
					else {
						est store rho_3
					}
				}
				else {
					est restore rho_6
				}
			}
			else {
				est restore rho_6
			}
		}
	}
	else {
		*
		******************* case Ir == 0 & IR == 0 *******************
		*
		di ""
		di "******* case IR =`IR' **********"
		di "******* case Ir =`Ir' **********"
		di "******* case IA =`IA' **********"
		di ""
		qui reg x seas*
		est store rho_6
		if `Idelta' == 0 & `T'-`IR'-`Ir'-`IA' > 10  {
			predict eps, resid
			qui ttest eps == 0
			local pval = r(p)
			if `pval' < 0.01 {						
				gen trend_a = _n/`S' 
				local Itau = 1
				qui reg x trend_a seas*
				est store rho_3
			}
			else {
				est restore rho_6
			}
		}
		else {
			est restore rho_6
		}
	}
}
else {
	*
	****************** case IA == 0
	*
	if `Ir' == 1 & `IR' == 1 & `Idelta' == 0 & `Itau' == 0 { 
		*
		******************* case Ir == 1 & IR == 1 *******************
		*
		di ""
		di "******* case IR =`IR' **********"
		di "******* case Ir =`Ir' **********"
		di "******* case IA =`IA' **********"
		di ""
		qui reg x l.x l`S'.x 
		if _b[l.x]> .5 & _b[l.x] + 2*_se[l.x] > .9 {
			local Idelta = 1
			constraint 1 l.x = 1
			qui cnsreg x l.x l`S'.x, constraint(1)
			est store rho_2_1 
		}	
		else if _b[l.x]< 0 {
			local Ir = 0
			qui reg x
			est store rho_2_2
			if `Idelta' == 0 & `T'-`IR'-`Ir'-`IA' > 10  {
				predict eps, resid
				g cum_eps = sum(eps) if eps!=.
				qui ttest cum_eps == 0
				local pval = r(p)
				if `pval' < 0.01 {						
					gen trend_a = _n/`S' 
					local Itau = 1
					qui reg x l.x l`S'.x trend_a 
					if _b[l.x]< -0.5 {
						local Itau = 0
						est restore rho_6
					}
					else {
						est store rho_2_2_3
					}
				}
				else {
					est restore rho_2_2
				}
			}
			else restore rho_2_2
		}
		else {
			if `Idelta' == 0 & `T'-`IR'-`Ir'-`IA' > 10  {
				predict eps, resid
				g cum_eps = sum(eps) if eps!=.
				qui ttest cum_eps == 0
				local pval = r(p)
				if `pval' < 0.01 {						
					gen trend_a = _n/`S' 
					local Itau = 1
					qui reg x l.x l`S'.x trend_a
					if _b[l.x]< -0.5 {
						est restore rho_6
					}
					else {
						est store rho_3
					}
				}
				else {
					est restore rho_6
				}
			}
			else {
				est restore rho_6
			}
		}
	}
	else if `Ir' == 1 & `IR' == 0 & `Idelta' == 0 & `Itau' == 0 {
		*
		******************* case Ir == 1 & IR == 0 *******************
		*
		di ""
		di "******* case IR =`IR' **********"
		di "******* case Ir =`Ir' **********"
		di "******* case IA =`IA' **********"
		di ""
		qui reg x l.x 
		est store rho_6
		if _b[l.x]> .5 & _b[l.x] + 2*_se[l.x] > .9 {
			local Idelta = 1
			constraint 1 l.x = 1
			qui cnsreg x l.x, constraint(1)
			est store rho_2_1 
		}	
		else if _b[l.x]< 0 {
			local Ir = 0
			qui reg x
			est store rho_2_2
			if `Idelta' == 0 & `T'-`Ir'-`IA' > 10  {
				predict eps, resid
				qui ttest eps == 0
				local pval = r(p)
				if `pval' < 0.01 {						
					gen trend_a = _n/`S' 
					local Itau = 1
					qui reg x l.x trend_a
					if _b[l.x]< -0.5 {
						local Itau = 0
						est restore rho_6
					}
					else {
						est store rho_2_2_3
					}
				}
				else {
					est restore rho_2_2
				}
			}
			else {
				est restore rho_2_2
			}
		}
		else {
			if `Idelta' == 0 & `T'-`IR'-`Ir'-`IA' > 10  {
				predict eps, resid
				qui ttest eps == 0
				local pval = r(p)
				if `pval' < 0.01 {						
					gen trend_a = _n/`S' 
					local Itau = 1
					qui reg x l.x trend_a 
					if _b[l.x]< -0.5 {
						est restore rho_6
					}
					else {
						est store rho_3
					}
				}
				else {
					est restore rho_6
				}
			}
			else {
				est restore rho_6
			}
		}
	}
	else {
		*
		******************* case Ir == 0 & IR == 0 *******************
		*
		di ""
		di "******* case IR =`IR' **********"
		di "******* case Ir =`Ir' **********"
		di "******* case IA =`IA' **********"
		di ""
		qui reg x 
		est store rho_6
		if `Idelta' == 0 & `T'-`IR'-`Ir'-`IA' > 10  {
			predict eps, resid
			qui ttest eps == 0
			local pval = r(p)
			if `pval' < 0.01 {						
				gen trend_a = _n/`S' 
				local Itau = 1
				qui reg x trend_a 
				est store rho_3
			}
			else {
				est restore rho_6
			}
		}
		else {
			est restore rho_6
		}
	}
}
*
******************** compute rho forecast
*
local mod_name = e(_estimates_name)
local sigmahat = e(rmse)
local esse = 1.645*`sigmahat'
*
if `IA' == 1 {
	*
	if "`mod_name'" == "rho_6" | "`mod_name'" == "rho_2_1"  | "`mod_name'" == "rho_2_2" {
		*
		*
		forecast create pluto, replace
		forecast estimates `mod_name'
		forecast exogenous seas* 
		forecast solve, prefix(rho_) begin(`inizio') periods(`H')
		*
	}
	else {
		*
		forecast create pluto, replace
		forecast estimates `mod_name'
		forecast exogenous trend_a seas* 
		forecast solve, prefix(rho_) begin(`inizio') periods(`H')
		*
	}
	*
	*
	if `Idelta' == 1 {
		gen upci = rho_x + `esse'
		gen loci = rho_x - `esse'
		*
		di "adjustment of the random walk forecsst"
		*
		gen indice = _n
		gen x_tilde = rho_x
		local TH = `T'+`H'
		replace x_tilde = loci if indice > `TH' & rho_x > 0 & loci > 0
		replace x_tilde = upci if indice > `TH' & rho_x <= 0 & upci < 0
		drop rho_x
		rename x_tilde rho_x
	}
}
else {
	*
	if "`mod_name'" == "rho_6" | "`mod_name'" == "rho_2_1"  | "`mod_name'" == "rho_2_2" {
		*
		*
		forecast create pluto, replace
		forecast estimates `mod_name'
		forecast solve, prefix(rho_) begin(`inizio') periods(`H')
		*
	}
	else {
		*
		forecast create pluto, replace
		forecast estimates `mod_name'
		forecast exogenous trend_a
		forecast solve, prefix(rho_) begin(`inizio') periods(`H')
		*
	}
	*
	*
	if `Idelta' == 1 {
		gen upci = rho_x + `esse'
		gen loci = rho_x - `esse'
		*
		di "adjustment of the random walk forecsst"
		*
		gen indice = _n
		gen x_tilde = rho_x
		local TH = `T'+`H'
		replace x_tilde = loci if indice > `TH' & rho_x > 0 & loci > 0
		replace x_tilde = upci if indice > `TH' & rho_x <= 0 & upci < 0
		drop rho_x
		rename x_tilde rho_x
	}
}
*
***************************************************************
************************ Averaging ****************************
***************************************************************
*
gen x_tilde = (delta_x + rho_x)/2
*
*
***************************************************************
************************ Calibration **************************
***************************************************************
*
keep x t Time ipi ipi_sa x_tilde delta_x rho_x year month seas*
qui tsset t 
*
******************* generate indicator functions and variables
*
local T_c = `T'+`H'
*
*
if `T' > 4*`S' {
	local I4 = 1
}
else {
	local I4 = 0
}
*
*
if `IR' * `Irho' * `I4' > 0 {
	local Iplus = 2
}
else {
	local Iplus = 0
}
*
*
if `IA' == 0 {
	local I_notA = 2
}
else {
	local I_notA = 0
}
*	
*
if `T_c'-`Irho'-`Iplus'-`I_notA'-`IA'-1 > 10 & `S' != 24 & `T' > 3*`S' {
	local I6 = 1
}
else {
	local I6 = 0
}
*
*
if `Irho' == 1 & `S' == 4 {
	local I5 = 1
}
else if `Irho' == 1 & `S' == 12 {
	local I5 = 1
}
else if `Irho' == 1 & `S' == 13 {
	local I5 = 1
}
else {
	local I5 = 0
}
*
*
local S_minus = `S'+1
*
*
******************* generate variables
*
*
gen trend = _n
gen d_t = 1 if trend < `T' - min(2*`S', `T_c'/2) 
replace d_t = 0 if d_t != 1
gen D_t = trend*d_t
*
*
gen S_1t = sin(2*_pi*trend/`S') 
gen C_1t = cos(2*_pi*trend/`S') 
*
*
*
***************** calibration tree
*
if `Irho' == 1 & `Iplus' == 2 {  
	if `IA' == 1 {
		if `I6' == 1 {
			if `I5' == 1 {
				qui reg x_tilde l.x_tilde l`S'.x_tilde l`S_minus'.x_tilde ///
				 seas* d_t D_t
			}
			else {		/* case I5 == 0  */
				qui reg x_tilde l.x_tilde l`S'.x_tilde l`S_minus'.x_tilde ///
				seas* d_t
			}
		}
		else {		/* case I6 == 0  */
			qui reg x_tilde l.x_tilde l`S'.x_tilde l`S_minus'.x_tilde seas*
		}
	}
	else {		/* case IA == 0 */
		if `I6' == 1 {
			if `I5' == 1 {
				qui reg x_tilde l.x_tilde l`S'.x_tilde l`S_minus'.x_tilde ///
				S_1t C_1t  d_t D_t
			}
			else {		/* case I5 == 0 */
				qui reg x_tilde l.x_tilde l`S'.x_tilde l`S_minus'.x_tilde ///
				S_1t C_1t d_t
			}
		}
		else {		/* case I6 == 0  */
			qui reg x_tilde l.x_tilde l`S'.x_tilde l`S_minus'.x_tilde S_1t C_1t 
		}
	}
}
else if `Irho' == 1 & `Iplus' == 0 { 
	if `IA' == 1 {
		if `I6' == 1 {
			if `I5' == 1 {
				qui reg x_tilde l.x_tilde seas* d_t D_t
			}
			else {		/* case I5 == 0 */
				qui reg x_tilde l.x_tilde seas* d_t
			}
		}
		else {		/* case I6 == 0 */
			qui reg x_tilde l.x_tilde seas*
		}
	}
	else {		/* case IA == 0 */
		if `I6' == 1 {
			if `I5' == 1 {
				qui reg x_tilde l.x_tilde S_1t C_1t d_t D_t
			}
			else {		/* case I5 == 0 */
				qui reg x_tilde l.x_tilde S_1t C_1t d_t
			}
		}
		else {		/* case I6 == 0 */
			qui reg x_tilde l.x_tilde S_1t C_1t 
		}
	}
}
else {		/* case without autoregressive component */
	if `IA' == 1 {
		if `I6' == 1 {
			qui reg x_tilde seas* d_t
		}
		else {		/* case I6 == 0 */
			qui reg x_tilde seas*
		}
	}
	else {		/* case IA == 0 */
		if `I6' == 1 {
			qui reg x_tilde S_1t C_1t d_t
		}
		else {		/* case I6 == 0 */
			qui reg x_tilde S_1t C_1t 
		}
	}
}
*
*
predict ford_x
predict dev_ford, stdp
g lb = ford_x - invnormal(.975)*dev_ford if trend > `T'
g ub = ford_x + invnormal(.975)*dev_ford if trend > `T'
*
*
replace ford_x = x if trend <= `T'
*
*
tsline ford_x /*ub lb*/ x_tilde rho_x delta_x  if trend >= `T'
*
*
************************************
* recover the levels of the forecast
************************************
*
g delta_`1' =exp(delta_x)
g rho_`1'   =exp(rho_x)
g avg_`1'   =exp(x_tilde)
g ford_`1'  =exp(ford_x)
*
keep t Time `1' delta_`1' rho_`1' avg_`1' ford_`1'
*
*
*
end
*
*
******************************** MAIN PROGRAM **********************************
******************************** MAIN PROGRAM **********************************
******************************** MAIN PROGRAM **********************************
*
*
* FAI ATTENZIONE A CAMBIARE LE DIRECTORY!!
*
cd C:\Users\Admin\OneDrive\Desktop\Tesi\STATA
*cd g:\data\tesi
*
capture log close
log using tesi_Mirco.log, replace
*
import excel ipi.xlsx, sheet("stata") cellrange(a1:c790) firstrow clear
g      t = tm(1955m1) + _n -1
format t %tm
tsset  t
*
compress
order t Time ipi ipi_sa
summ
* salvo la banca dati completa prima del calcolo delle
* previsioni ricorsive (oppure rolling)
qui save whole_dataset, replace
summarize
set more off
******************** run the program to make the section 2.1 
******************** of Doornik et al (IJF, 2020)
ford ipi 18
*
local S = primary_S
local Irho = I_rho
local T = estimation_T
local IA = I_A
local IR = I_R
local H = sample_H
*
summ
*
*
log close
