*===============================================================================
* Speed benchmarks
*===============================================================================

* Warm up ctools to load plugin before benchmarks
sysuse auto, clear
csort price

* List of obs
local n_obs_list 1e5 1e6 1e7 1e8
local n_obs_list 1e5 1e6 1e7 1e8

* Setup: matrix to store results and track current row
matrix A = J(4,3,.)
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
	
	* Csort times
	timer on 2
	noi csort float2, verbose
	timer off 2
	
	
	* Save timers to a row
	qui timer list
	mat A[`current',1] = `N'
	mat A[`current',2] = r(t1)
	mat A[`current',3] = r(t2)
	timer clear
	local current = `current' + 1
	
}

clear
svmat A
rename (A1 A2 A3) (n_obs time_sort_float time_csort_float)
gen rel_time = time_sort_float/time_csort_float

twoway (scatter rel_time n_obs)
gen n_obs_2 = _n
twoway (scatter rel_time n_obs_2), xlabel(1 "100k" 2 "1m" 3 "10m" 4 "100m") xtitle("")


/*

-------------------------------------------------------
csort timing breakdown:
-------------------------------------------------------
  C plugin internals:
    Data load:                0.9943 sec
    Sort (compute order):     1.1066 sec
    Sort (apply permute):     3.2922 sec
    Data store:               2.1853 sec
    Memory cleanup:           0.1238 sec
  -----------------------------------------------------
    C plugin total:           7.7022 sec
  -----------------------------------------------------
  Stata overhead:
    Pre-plugin parsing:       0.0670 sec
    Plugin call overhead:     0.0018 sec
    Post-plugin sort:         1.1020 sec
  -----------------------------------------------------
    Stata overhead total:     1.1708 sec
-------------------------------------------------------
    Wall clock total:         8.8730 sec
-------------------------------------------------------
r; t=80.25 12:00:39

*/
