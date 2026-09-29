/* End-to-end transfer workloads with wide mixed payloads.
   args: plugin directory, output directory, label, rows, repetitions, threads
   Run through the stata shell alias. Timings exclude input generation and verification. */
args plugin_dir out label n reps threads
clear all
set more off
set linesize 255
set seed 26092026
set sortseed 26092026
if "`threads'" == "" local threads 8
adopath ++ "`plugin_dir'"
capture log close _all
log using "`out'/wide_`label'_`n'.log", text replace
program transport_clock, plugin using("`out'/clock.plugin")
local master "`out'/wide_master_`n'.dta"
local using "`out'/wide_using_`n'.dta"
capture confirm file "`master'"
if _rc {
    quietly set obs `n'
    quietly gen long id = _n
    quietly gen int group = mod(id*17,10000)
    quietly gen str32 skey = "g" + string(group,"%05.0f") + "xxxxxxxxxxxxxxxxxxxxxxxxxx"
    quietly gen double x = id/7
    forvalues j = 5/20 {
        if mod(`j',2)==0 {
            quietly gen str32 v`j' = substr("v`j'_" + string(id,"%09.0f") + "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",1,32)
            quietly replace v`j' = "" if mod(id,17)==0
        }
        else {
            local kind = mod(`j',5)
            if `kind' == 0 quietly gen byte v`j' = mod(id+`j',101)-50
            if `kind' == 1 quietly gen int v`j' = mod(id+`j',30000)
            if `kind' == 2 quietly gen long v`j' = id+`j'
            if `kind' == 3 quietly gen float v`j' = id/`j'
            if `kind' == 4 quietly gen double v`j' = id/`j'
            quietly replace v`j' = .z if mod(id,103)==0
        }
    }
    quietly gen double shuffle = runiform()
    quietly sort shuffle
    drop shuffle
    quietly gen long origin = _n
    quietly save "`master'", replace
}
capture confirm file "`using'"
if _rc {
    quietly use "`master'", clear
    keep group
    quietly duplicates drop
    quietly gen double u = group/7
    quietly gen str32 us = "u" + string(group,"%05.0f") + "xxxxxxxxxxxxxxxxxxxxxxxxxx"
    quietly save "`using'", replace
}
sysuse auto, clear
csort price, threads(`threads')
forvalues rep = 1/`reps' {
    foreach case in sort_id sort_string sort_stream merge_group {
        use "`master'", clear
        plugin call transport_clock, "start"
        if "`case'" == "sort_id" csort id, verbose threads(`threads') nostream
        if "`case'" == "sort_string" csort skey x, verbose threads(`threads') nostream
        if "`case'" == "sort_stream" csort skey x, verbose threads(`threads') stream(4)
        if "`case'" == "merge_group" cmerge m:1 group using "`using'", verbose threads(`threads')
        plugin call transport_clock, "stop"
        local elapsed = scalar(__ctools_perf_seconds)
        di "COMMAND,`label',`n',`rep',`case',`threads',`elapsed'"
        if "`case'" == "merge_group" {
            di "PHASE,load," _cmerge_p2_load_master ",store," _cmerge_p2_store ",total," _cmerge_p2_total
            assert _merge == 3
            assert u == group/7
            assert us == "u" + string(group,"%05.0f") + "xxxxxxxxxxxxxxxxxxxxxxxxxx"
        }
        else {
            di "PHASE,load," _csort_time_load ",store," _csort_time_store ",stream," _csort_time_stream ",total," _csort_time_total
            if "`case'" == "sort_id" assert id == _n
            if inlist("`case'", "sort_string", "sort_stream") {
                assert skey >= skey[_n-1] if _n > 1
                assert x >= x[_n-1] if _n > 1 & skey == skey[_n-1]
            }
        }
        assert _N == `n'
        assert group == mod(id*17,10000)
        assert x == id/7
        assert skey == "g" + string(group,"%05.0f") + "xxxxxxxxxxxxxxxxxxxxxxxxxx"
        if `rep' == 1 {
            forvalues j = 5/20 {
                if mod(`j',2)==0 {
                    assert v`j' == cond(mod(id,17)==0,"",substr("v`j'_" + string(id,"%09.0f") + "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",1,32))
                }
                else {
                    assert v`j' == .z if mod(id,103)==0
                    local kind = mod(`j',5)
                    if `kind' == 0 assert v`j' == mod(id+`j',101)-50 if mod(id,103)!=0
                    if `kind' == 1 assert v`j' == mod(id+`j',30000) if mod(id,103)!=0
                    if `kind' == 2 assert v`j' == id+`j' if mod(id,103)!=0
                    if `kind' == 3 assert v`j' == float(id/`j') if mod(id,103)!=0
                    if `kind' == 4 assert v`j' == id/`j' if mod(id,103)!=0
                }
            }
        }
        di "CHECK,`case',passed"
    }
}
di "TRANSPORT_COMMANDS_COMPLETE"
log close
exit, clear
