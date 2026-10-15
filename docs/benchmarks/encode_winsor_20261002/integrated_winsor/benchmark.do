* One frozen plugin per process. Timers include the entire public command.
args buildpath label outdir cases reps threadlist rows refreps
clear all
set more off
set linesize 255
set seed 928143
adopath ++ "`buildpath'"
ctools, version
about
which rangestat
which gstats
program perf_clock, plugin using("`outdir'/clock.plugin")
if "`reps'" == "" local reps 3
if "`refreps'" == "" local refreps `reps'
if "`threadlist'" == "" local threadlist 1 8
if "`rows'" == "" local rows 1000000
if "`cases'" == "" local cases rs20 rs200 rs1000 rs5000 rs_single rs_ties rs_median rs_skew w20 w200 w1000 w_single w_wide
foreach case of local cases {
    clear
    set seed 928143
    quietly set obs `rows'
    local size 200
    if inlist("`case'", "rs20", "w20") local size 20
    if inlist("`case'", "rs1000", "w1000") local size 1000
    if "`case'" == "rs5000" local size 5000
    if inlist("`case'", "rs_single", "w_single") local size `rows'
    quietly gen long id = _n
    quietly gen long g = floor((_n-1)/`size')
    quietly gen double key = mod(_n-1,`size')
    if "`case'" == "rs_ties" quietly replace key = floor(key/20)
    if "`case'" == "rs_median" quietly replace key = floor(key/10)
    if "`case'" == "rs_skew" quietly replace g = 0 if _n <= .9*`rows'
    quietly gen double x = 100 + 5*rnormal()
    quietly replace x = . if mod(id,37)==0
    quietly replace x = .a if mod(id,107)==0
    local vars x
    if "`case'" == "w_wide" {
        forvalues k=2/8 {
            quietly gen double x`k' = 100 + 5*rnormal()
            local vars `vars' x`k'
        }
    }
    * Randomize row order to include sorting and restoration costs.
    quietly gen double shuffle = runiform()
    sort shuffle
    drop shuffle
    local stat mean
    if "`case'" == "rs_median" local stat median
    local isrange = (substr("`case'",1,2)=="rs")
    if `isrange' {
        local command crangestat (`stat') result=x, interval(key -5 5) by(g)
        local reference rangestat (`stat') result=x, interval(key -5 5) by(g)
        local results result
    }
    else {
        local command cwinsor `vars', by(g) suffix(_w) cuts(1 99)
        local reference gstats winsor `vars', by(g) suffix(_w) cuts(1 99)
        local results
        foreach v of local vars {
            local results `results' `v'_w
        }
    }
    tempfile base expected
    quietly save `base'
    foreach threads of local threadlist {
        forvalues rep=0/`reps' {
            quietly use `base', clear
            plugin call perf_clock, "start"
            quietly `command' threads(`threads')
            plugin call perf_clock, "stop"
            di "PERF,`label',`case',ctools,`rows',`threads',`rep'," %12.7f scalar(__ctools_perf_seconds)
            if `isrange' {
                di "PHASE,`label',`case',`threads',`rep'," %12.7f scalar(_crangestat_time_load) "," %12.7f scalar(_crangestat_time_sort) "," %12.7f scalar(_crangestat_time_compute) "," %12.7f scalar(_crangestat_time_store)
            }
            else {
                di "PHASE,`label',`case',`threads',`rep'," %12.7f scalar(_cwinsor_time_load) "," %12.7f scalar(_cwinsor_time_sort) "," %12.7f scalar(_cwinsor_time_winsor) "," %12.7f scalar(_cwinsor_time_store)
            }
            if `rep' == 0 {
                if "`label'" == "b1" quietly save "`outdir'/expected_`case'_`threads'.dta", replace
                else cf _all using "`outdir'/expected_`case'_`threads'.dta"
                quietly keep id `results'
                quietly save `expected', replace
                quietly use `base', clear
                quietly `reference'
                foreach v of local results {
                    rename `v' ref_`v'
                }
                quietly merge 1:1 id using `expected', assert(match) nogenerate
                foreach v of local results {
                    assert missing(`v') == missing(ref_`v')
                    assert `v' == ref_`v' | abs(`v'-ref_`v')/max(abs(`v'),abs(ref_`v')) < 1e-7 if !missing(`v')
                }
                di "CORRECT,`label',`case',`threads',all_rows"
            }
        }
    }
    * Reference timings use the same saved input, outside comparison/merge work.
    forvalues rep=0/`refreps' {
        quietly use `base', clear
        plugin call perf_clock, "start"
        quietly `reference'
        plugin call perf_clock, "stop"
        di "PERF,`label',`case',reference,`rows',0,`rep'," %12.7f scalar(__ctools_perf_seconds)
    }
}
di "GROUP_BENCHMARK_COMPLETE"
