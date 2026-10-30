* Named benchmark workloads. All original sizes, seeds and timing checks are retained.
if "`__workload'" == "" {
    gettoken __workload 0 : 0
    macro shift
}

if "`__workload'" == "encode_performance" {
* Frozen public commands; setup, comparisons and disk reads are outside timers.
args buildpath label outdir cases reps threadlist rows refreps datadir
if `"`datadir'"'=="" local datadir `"`outdir'"'
clear all
set more off
set linesize 255
adopath ++ "`buildpath'"
ctools, version
about
which encode
program perf_clock, plugin using("`outdir'/clock.plugin")
if "`cases'" == "" local cases e10 e1000 e10000 e50000 e_long e_unicode e_filtered e_existing
foreach case of local cases {
    clear
    set seed 928143
    local n = `rows'
    if "`case'" == "e_long" local n = min(`rows',100000)
    quietly set obs `n'
    local categories 1000
    if "`case'" == "e10" local categories 10
    if "`case'" == "e10000" local categories 10000
    if "`case'" == "e50000" local categories 50000
    quietly gen long id = _n
    quietly gen str12 x = "value_" + string(mod(_n-1,`categories'), "%05.0f")
    if "`case'" == "e_long" {
        quietly recast str2045 x
        mata: st_sstore(., "x", st_sdata(., "x") :+ invtokens(J(1,2034,"a"), ""))
    }
    if "`case'" == "e_unicode" {
        quietly recast str40 x
        quietly replace x = "é東京" + x
    }
    quietly replace x = "" if mod(id,37)==0
    quietly gen double shuffle = runiform()
    sort shuffle
    drop shuffle
    local subset
    if "`case'" == "e_filtered" local subset if mod(id,3)==0 in 10/`n'
    local labelopt
    if "`case'" == "e_existing" {
        quietly encode x, gen(seed) label(existing)
        * Keep seed attached so save/use retains the existing label definition.
        local labelopt label(existing)
    }
    local command cencode x `subset', gen(y) `labelopt'
    local reference encode x `subset', gen(y) `labelopt'
    tempfile base expected
    quietly save `base'
    foreach threads of local threadlist {
        forvalues rep=0/`reps' {
            quietly use `base', clear
            plugin call perf_clock, "start"
            quietly `command' threads(`threads')
            plugin call perf_clock, "stop"
            di "PERF,`label',`case',ctools,`n',`threads',`rep'," %12.7f scalar(__ctools_perf_seconds)
            foreach phase in load collect sort encode labels total {
                di "PHASE,`label',`case',`threads',`rep',`phase'," %12.7f scalar(_cencode_time_`phase')
            }
            if `rep' == 0 {
                quietly decode y, gen(text)
                quietly keep id y text
                if "`label'" == "b1" quietly save "`datadir'/expected_`case'_`threads'.dta", replace
                else cf _all using "`datadir'/expected_`case'_`threads'.dta"
                quietly save `expected', replace
                quietly use `base', clear
                quietly `reference'
                quietly decode y, gen(text)
                cf id y text using `expected'
                di "CORRECT,`label',`case',`threads',codes_and_label_bytes_all_rows"
            }
        }
    }
    forvalues rep=0/`refreps' {
        quietly use `base', clear
        plugin call perf_clock, "start"
        quietly `reference'
        plugin call perf_clock, "stop"
        di "PERF,`label',`case',reference,`n',0,`rep'," %12.7f scalar(__ctools_perf_seconds)
    }
}
di "ENCODE_BENCHMARK_COMPLETE"
exit
}

if "`__workload'" == "group_performance" {
* One frozen plugin per process. Timers include the entire public command.
args buildpath label outdir cases reps threadlist rows refreps datadir
if `"`datadir'"'=="" local datadir `"`outdir'"'
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
                if "`label'" == "b1" quietly save "`datadir'/expected_`case'_`threads'.dta", replace
                else cf _all using "`datadir'/expected_`case'_`threads'.dta"
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
exit
}

if "`__workload'" == "cpsmatch_paths" {
/* One workload per process. Reloads, warmups and comparisons are untimed. */
args buildpath clockpath case N teams reps label seed reference
clear all
set more off
set linesize 255
set processors 4
adopath ++ "`buildpath'"
program path_clock, plugin using("`clockpath'")
quietly set obs `N'
set seed `seed'
set sortseed `seed'
quietly gen long row_id=_n
quietly gen double x1=rnormal()
quietly gen double x2=rnormal()
quietly gen double ps=invlogit(.6*x1+.4*x2)
quietly gen byte treat=runiform()<ps
quietly gen double y=2*treat+x1+x2+rnormal()
local opts "pscore(ps) outcome(y)"
local restriction ""
if substr("`case'",1,2)=="nr" local opts "`opts' noreplacement"
if "`case'"=="nr_equal" quietly replace ps=.5
if "`case'"=="nr_desc" local opts "`opts' descending"
if "`case'"=="nr_common" local opts "`opts' common caliper(.01)"
if substr("`case'",1,6)=="radius" local opts "`opts' radius caliper(.05)"
if "`case'"=="radius_narrow" local opts "pscore(ps) outcome(y) radius caliper(.00001)"
if "`case'"=="radius_default" local opts "pscore(ps) outcome(y) radius"
if "`case'"=="radius_gaps" {
    quietly replace ps=cond(treat,.2,.8) if mod(row_id,2)==0
    quietly replace ps=round(ps,.1)
    local opts "pscore(ps) outcome(y) radius caliper(.01) common"
}
if "`case'"=="radius_bounds" {
    quietly replace ps=mod(row_id,8)/8
    local opts "pscore(ps) outcome(y) radius caliper(.125)"
}
if "`case'"=="restricted" {
    local restriction "if mod(row_id,3)!=0 in 11/`=`N'-10'"
    quietly replace y=.a if mod(row_id,17)==0
    quietly replace ps=.b if mod(row_id,19)==0
    local opts "pscore(ps) outcome(y) radius caliper(.05)"
}
if "`case'"=="kernel" local opts "`opts' kernel bwidth(.06)"
tempfile base
quietly save `base'
foreach team of local teams {
    forvalues rep=0/`reps' {
        quietly use `base', clear
        plugin call path_clock, "start"
        quietly cpsmatch treat `restriction', `opts' threads(`team')
        plugin call path_clock, "stop"
        local wall : display %21.9f scalar(__ctools_perf_seconds)
        local results "RESULT,`label',`case',`N',`team',`rep'"
        foreach value in N n_treated n_controls n_matched n_off_support common_min common_max att att_se att_t {
            local value : display %24.16g r(`value')
            local results "`results',`=strtrim("`value'")'"
        }
        local sample "PATH,`label',`case',`N',`team',`rep',`=strtrim("`wall'")'"
        foreach stage in load setup sort match stats store total {
            capture confirm scalar _cpsmatch_time_`stage'
            if _rc local value "0"
            else local value : display %21.9f scalar(_cpsmatch_time_`stage')
            local sample "`sample',`=strtrim("`value'")'"
        }
        assert _N==`N'
        if `rep'>0 {
            di "`sample'"
            di "`results'"
        }
        if "`reference'"!="" & `rep'==1 {
            local path "`reference'/`case'_`N'_`team'.dta"
            if "`label'"=="b1" quietly save "`path'", replace
            else {
                ds _weight, not
                quietly cf `r(varlist)' using "`path'"
                preserve
                rename _weight candidate_weight
                quietly merge 1:1 row_id using "`path'", keepusing(_weight) nogen
                assert candidate_weight==_weight if _weight==0 | missing(_weight)
                assert candidate_weight>0 if _weight>0 & !missing(_weight)
                assert abs(candidate_weight-_weight)/max(abs(candidate_weight),abs(_weight))<=10^(-7) | candidate_weight==_weight
                restore
                di "PARITY,`label',`case',`N',`team',PASS"
            }
        }
    }
}
di "PATH_COMPLETE,`label',`case'"
exit
}

if "`__workload'" == "winsor_references" {
* End-to-end comparisons with both winsor2 and gstats, one frozen ctools build.
args buildpath clockpath rows reps
clear all
set more off
set linesize 255
adopath ++ "`buildpath'"
ctools, version
which winsor2
which gstats
program perf_clock, plugin using("`clockpath'")
foreach case in grouped single wide noby {
    clear
    set seed 928143
    quietly set obs `rows'
    quietly gen long id=_n
    quietly gen long g=floor((_n-1)/200)
    if "`case'"=="single" quietly replace g=0
    quietly gen double x=100+5*rnormal()
    quietly replace x=. if mod(id,37)==0
    quietly replace x=.a if mod(id,107)==0
    local vars x
    if "`case'"=="wide" {
        forvalues k=2/8 {
            quietly gen double x`k'=100+5*rnormal()
            local vars `vars' x`k'
        }
    }
    quietly gen double shuffle=runiform()
    sort shuffle
    drop shuffle
    local byopt by(g)
    if "`case'"=="noby" local byopt
    tempfile base expected
    quietly save `base'
    foreach operation in ctools gstats winsor2 {
        local command cwinsor `vars', `byopt' cuts(1 99) suffix(_w) threads(8)
        if "`operation'"=="gstats" local command gstats winsor `vars', `byopt' cuts(1 99) suffix(_w)
        if "`operation'"=="winsor2" local command winsor2 `vars', `byopt' cuts(1 99) suffix(_w)
        forvalues rep=0/`reps' {
            quietly use `base', clear
            plugin call perf_clock, "start"
            quietly `command'
            plugin call perf_clock, "stop"
            di "WINSOR_REF,`case',`operation',`rows',`rep'," %12.7f scalar(__ctools_perf_seconds)
            if `rep'==0 {
                foreach v of local vars {
                    * winsor2 suffix() turns extended missings into ordinary .;
                    * compare all missing locations and all finite values.
                    quietly replace `v'_w=. if missing(`v'_w)
                }
                if "`operation'"=="ctools" quietly save `expected'
                else {
                    foreach v of local vars {
                        rename `v'_w `v'_ref
                    }
                    quietly merge 1:1 id using `expected', assert(match) nogenerate
                    foreach v of local vars {
                        assert missing(`v'_w)==missing(`v'_ref)
                        assert `v'_w==`v'_ref | abs(`v'_w-`v'_ref)/max(abs(`v'_w),abs(`v'_ref))<1e-7 if !missing(`v'_w)
                    }
                }
            }
        }
    }
}
di "WINSOR_REFERENCES_COMPLETE"
exit
}

di as error "Unknown workload: `__workload'"
exit 198
