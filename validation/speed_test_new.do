*===============================================================================
* Speed benchmarks
*===============================================================================

* Warm up ctools to load plugin before benchmarks
sysuse auto, clear
csort price

* List of obs
local n_obs_list 1e5 1e6 1e7

* Setup: matrix to store results and track current row
matrix A = J(4,5,.)
local current = 1

* Iterate over obs counts
foreach N in `n_obs_list' {
	
	* Create dataset
	clear
	set obs `=`N''
	forval i=1/10 {
		gen float`i' = rnormal()
	}
	forval i=1/5 {
		gen int`i' = ceil(runiform()*1000) 
	}
	tempfile t1
	save `t1'
	
	* Sort times
	timer on 1
	sort float1
	timer off 1
	
	timer on 2
	sort int1
	timer off 2
	
	* Csort times
	timer on 3
	noi csort float2, verbose
	timer off 3
	
	timer on 4
	noi csort int2, verbose
	timer off 4
	
	* Save timers to a row
	qui timer list
	mat A[`current',1] = `N'
	mat A[`current',2] = r(t1)
	mat A[`current',3] = r(t2)
	mat A[`current',4] = r(t3)
	mat A[`current',5] = r(t4)
	timer clear
	local current = `current' + 1
	
}

clear
svmat A
rename (A1 A2 A3 A4 A5) (n_obs time_sort_float time_sort_int time_csort_float time_csort_int)
gen relative_time_float = time_sort_float/time_csort_float
gen relative_time_int = time_sort_int/time_csort_int

twoway (scatter relative_time_float n_obs)
gen n_obs_2 = _n
twoway (scatter relative_time_float n_obs_2), xlabel(1 "100k" 2 "1m" 3 "10m" 4 "100m") xtitle("")
