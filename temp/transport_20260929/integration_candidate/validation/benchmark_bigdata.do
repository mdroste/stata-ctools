* Reproducible end-to-end sort/merge profiles. Run through the stata shell alias.
* args: plugin directory, output directory, label, N, K, repetitions, cases,
*       csort options (quoted), threads
args plugin_dir out label n k reps cases sortopts threads
clear all
set more off
set linesize 255
set seed 24092026
set sortseed 24092026
if "`threads'" == "" local threads 8
adopath ++ "`plugin_dir'"
capture log close _all
log using "`out'/`label'_`n'_`k'.log", text replace
local master "`out'/master_`n'_`k'.dta"
local small "`out'/using_group.dta"
local large "`out'/using_id_`n'.dta"
capture confirm file "`master'"
if _rc {
    set obs `n'
    gen long id = _n
    gen int group = mod(id * 17, 10000)
    gen str6 skey = "g" + string(group, "%05.0f")
    gen double x = rnormal()
    forvalues j = 5/`k' {
        local type = mod(`j', 5)
        if `type' == 0 gen str12 v`j' = "v`j'_" + string(mod(id, 100000), "%05.0f")
        if `type' == 1 gen byte v`j' = mod(id + `j', 100)
        if `type' == 2 gen int v`j' = mod(id + `j', 30000)
        if `type' == 3 gen float v`j' = mod(id + `j', 100000) / 17
        if `type' == 4 gen double v`j' = id / (`j' + 1)
    }
    gen double shuffle = runiform()
    sort shuffle
    drop shuffle
    save "`master'", replace
}
capture confirm file "`small'"
if _rc {
    clear
    set obs 10000
    gen int group = _n - 1
    gen double u = group / 7
    gen str9 us = "u" + string(group, "%05.0f")
    save "`small'", replace
}
capture confirm file "`large'"
if _rc {
    clear
    set obs `n'
    gen long id = _n + floor(`n' / 10)
    gen double u = id / 7
    gen str12 us = "u" + string(id, "%09.0f")
    gen double shuffle = runiform()
    sort shuffle
    drop shuffle
    save "`large'", replace
}
* Warm up the plugin outside measured work.
sysuse auto, clear
csort price, threads(`threads')
forvalues rep = 1/`reps' {
    foreach case of local cases {
        use "`master'", clear
        timer clear 1
        timer on 1
        if "`case'" == "sort_id" csort id, verbose threads(`threads') `sortopts'
        if "`case'" == "sort_group" csort group, verbose threads(`threads') `sortopts'
        if "`case'" == "sort_string" csort skey x, verbose threads(`threads') `sortopts'
        if "`case'" == "sort_float" csort x, verbose threads(`threads') `sortopts'
        if "`case'" == "merge_group" cmerge m:1 group using "`small'", verbose threads(`threads')
        if "`case'" == "merge_id" cmerge 1:1 id using "`large'", verbose threads(`threads')
        timer off 1
        quietly timer list 1
        local elapsed = r(t1)
        di "BENCH,`label',`n',`k',`rep',`case',`threads',`elapsed'"
        if substr("`case'",1,4) == "sort" {
            di "PHASE,load," _csort_time_load ",sort," _csort_time_sort ",permute," _csort_time_permute ",store," _csort_time_store ",stream," _csort_time_stream ",total," _csort_time_total
            assert _N == `n'
            if "`case'" == "sort_id" assert id == _n
            if "`case'" == "sort_group" assert group >= group[_n-1] if _n > 1
            if "`case'" == "sort_string" assert skey >= skey[_n-1] if _n > 1
            if "`case'" == "sort_float" assert x >= x[_n-1] if _n > 1
        }
        else {
            di "PHASE,load," _cmerge_p2_load_master ",sort," _cmerge_p2_sort_master ",join," _cmerge_p2_merge_join ",permute," _cmerge_p2_permute ",store," _cmerge_p2_store ",meta," _cmerge_p2_write_meta ",total," _cmerge_p2_total
            if "`case'" == "merge_group" {
                assert _N == `n'
                assert _merge == 3
                assert u == group / 7
                assert us == "u" + string(group, "%05.0f")
            }
            else {
                assert _N == `n' + floor(`n' / 10)
                assert _merge == cond(id <= floor(`n'/10), 1, cond(id > `n', 2, 3))
                assert u == id / 7 if _merge != 1
                assert us == "u" + string(id, "%09.0f") if _merge != 1
            }
        }
        * Check every payload cell once per case; repeated timings use the same
        * saved input and retain the ordering/match checks above on every run.
        if `rep' == 1 {
            assert group == mod(id * 17, 10000) if id <= `n'
            assert skey == "g" + string(group, "%05.0f") if id <= `n'
            forvalues j = 5/`k' {
                local type = mod(`j', 5)
                if `type' == 0 assert v`j' == "v`j'_" + string(mod(id, 100000), "%05.0f") if id <= `n'
                if `type' == 1 assert v`j' == mod(id + `j', 100) if id <= `n'
                if `type' == 2 assert v`j' == mod(id + `j', 30000) if id <= `n'
                if `type' == 3 assert v`j' == float(mod(id + `j', 100000) / 17) if id <= `n'
                if `type' == 4 assert v`j' == id / (`j' + 1) if id <= `n'
            }
        }
        di "CHECK,`case',passed"
    }
}
di "BIGDATA_COMPLETE"
log close
exit, clear
