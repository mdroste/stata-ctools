* Named benchmark workloads. All original sizes, seeds and timing checks are retained.
if "`__workload'" == "" {
    gettoken __workload 0 : 0
    macro shift
}

if "`__workload'" == "io_performance" {
* Frozen-build I/O comparison; prepare_io_performance.py creates ABBA drivers.
* threads=0 in the output denotes cimport excel's automatic worker count.
* rep=0 warms each path and verifies candidate imports against the baseline.
* Export artifacts from rep=0 support exact CSV/ZIP-part comparisons.
args buildpath label outdir cases reps datadir
if `"`datadir'"'=="" local datadir `"`outdir'"'
clear all
set more off
set linesize 255
adopath ++ "`buildpath'"
ctools, version
program perf_clock, plugin using("`outdir'/clock.plugin")
if "`cases'"=="" local cases numeric mixed wide strings
if "`reps'"=="" local reps 3
foreach case of local cases {
    clear
    set seed 9182
    local N 200000
    local K 10
    if "`case'"=="large" local N 1000000
    if "`case'"=="small" local N 1000
    if "`case'"=="wide" {
        local N 30000
        local K 100
    }
    if "`case'"=="strings" local N 50000
    quietly set obs `N'
    forvalues j=1/`K' {
        quietly gen double v`j'=rnormal()*1000
        quietly replace v`j'=. if mod(_n,101)==0
    }
    if "`case'"=="mixed" | "`case'"=="strings" {
        forvalues j=1/5 {
            quietly gen str80 s`j'="row "+string(_n)+", quoted "+char(34)+"text"+char(34)+" & <xml>"
        }
    }
    if "`case'"=="strings" {
        forvalues j=1/5 {
            quietly replace s`j'=s`j'+s`j'
        }
    }
    tempfile base
    quietly save `base'
    if substr("`label'",1,1)=="b" {
        quietly cexport delimited using "`datadir'/`case'.csv", replace
        quietly cexport excel using "`datadir'/`case'.xlsx", firstrow(variables) replace
    }
    foreach threads in 1 8 {
        foreach op in export_csv import_csv export_xlsx import_xlsx {
            forvalues rep=0/`reps' {
                quietly use `base', clear
                if "`op'"=="import_xlsx" & `threads'!=8 continue
                local verbose
                if `rep'==0 local verbose verbose
                di "PROFILE_BEGIN,`label',`case',`op',`threads',`rep'"
                plugin call perf_clock, "start"
                if "`op'"=="export_csv" cexport delimited using "`datadir'/output.csv", replace threads(`threads') `verbose'
                if "`op'"=="import_csv" cimport delimited using "`datadir'/`case'.csv", clear threads(`threads') `verbose'
                if "`op'"=="export_xlsx" cexport excel using "`datadir'/output.xlsx", replace firstrow(variables) threads(`threads') `verbose'
                if "`op'"=="import_xlsx" cimport excel using "`datadir'/`case'.xlsx", clear firstrow
                plugin call perf_clock, "stop"
                local reported_threads `threads'
                if "`op'"=="import_xlsx" local reported_threads 0
                di "PERF,`label',`case',`op',`reported_threads',`rep'," %12.6f scalar(__ctools_perf_seconds)
                if `rep'==0 & substr("`op'",1,6)=="export" {
                    local ext csv
                    if "`op'"=="export_xlsx" local ext xlsx
                    quietly copy "`datadir'/output.`ext'" "`datadir'/`label'_`case'_`threads'.`ext'", replace
                }
                if `rep'==0 & substr("`op'",1,6)=="import" {
                    if substr("`label'",1,1)=="b" quietly save "`datadir'/expected_`case'_`op'.dta", replace
                    else cf _all using "`datadir'/expected_`case'_`op'.dta"
                }
            }
        }
    }
}
exit
}

if "`__workload'" == "transport_commands" {
/* Whole-command transport acceptance benchmark.
   Args: frozen plugin/ado directory, output directory, variant label,
         repetitions (default 7), threads (default 8).
   The output directory must contain clock.plugin, built from
   validation/performance/benchmark_clock.c. Run only through stata.
   Rep 0 is a warmup; reps 1..7 are timed. Fixture generation, loading,
   native reference sorting and every value assertion are outside timing. */
args plugin_dir out label reps threads
clear all
set more off
set linesize 255
set seed 29092026
set sortseed 29092026
if "`reps'" == "" local reps 7
if "`threads'" == "" local threads 8
adopath ++ "`plugin_dir'"
capture log close workload
log using "`out'/commands_`label'.log", text replace name(workload)
which csort
program transport_clock, plugin using("`out'/clock.plugin")

program define adaptive_commands
    version 16
    args label reps threads
    foreach case in tiny_numeric narrow_numeric wide_numeric str2045_full {
        if "`case'" == "tiny_numeric" {
            local n 100
            local k 20
        }
        if "`case'" == "narrow_numeric" {
            local n 100000
            local k 20
        }
        if "`case'" == "wide_numeric" {
            local n 200000
            local k 128
        }
        if "`case'" == "str2045_full" {
            local n 100000
            local k 4
        }
        clear
        quietly set obs `n'
        /* 65537 is coprime to every fixture row count. This yields unique,
           shuffled IDs without timing random-number generation. */
        quietly gen long id = mod((_n-1)*65537,`n') + 1
        if "`case'" == "str2045_full" {
            local pad "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
            forvalues p = 1/6 {
                local pad "`pad'`pad'"
            }
            forvalues j = 2/`k' {
                quietly gen str2045 v`j' = substr("é" + string(id,"%09.0f") + "_`j'_" + "`pad'",1,2045)
                quietly replace v`j' = "" if mod(id+`j',17)==0
            }
        }
        else {
            forvalues j = 2/`k' {
                local kind = mod(`j'-2,5)
                if `kind' == 0 quietly gen byte v`j' = mod(id+`j',101)-50
                if `kind' == 1 quietly gen int v`j' = mod(id+7*`j',60001)-30000
                if `kind' == 2 quietly gen long v`j' = id*17+`j'
                if `kind' == 3 quietly gen float v`j' = id/`j'
                if `kind' == 4 quietly gen double v`j' = id/7+`j'/16
                quietly replace v`j' = . if mod(id,97)==0
                quietly replace v`j' = .a if mod(id,101)==0
                quietly replace v`j' = .z if mod(id,103)==0
            }
        }
        tempfile fixture
        quietly save "`fixture'", replace
        /* Native sorting supplies an independent whole-dataset checksum. */
        quietly sort id
        quietly datasignature
        local expected_signature "`r(datasignature)'"
        forvalues rep = 0/`reps' {
            quietly use "`fixture'", clear
            plugin call transport_clock, "start"
            quietly csort id, verbose threads(`threads') nostream
            plugin call transport_clock, "stop"
            local elapsed : display %21.9f scalar(__ctools_perf_seconds)
            local elapsed = strtrim("`elapsed'")
            local load : display %21.9f scalar(_csort_time_load)
            local load = strtrim("`load'")
            local store : display %21.9f scalar(_csort_time_store)
            local store = strtrim("`store'")

            assert _N == `n'
            unab result_vars : *
            local observed_k : word count `result_vars'
            assert `observed_k' == `k'
            assert id == _n
            local sortedby : sortedby
            assert "`sortedby'" == "id"
            if "`case'" == "str2045_full" {
                forvalues j = 2/`k' {
                    assert v`j' == cond(mod(id+`j',17)==0,"",substr("é" + string(id,"%09.0f") + "_`j'_" + "`pad'",1,2045))
                    assert strlen(v`j') == 2045 if mod(id+`j',17)!=0
                    assert v`j' == "" if mod(id+`j',17)==0
                }
            }
            else {
                forvalues j = 2/`k' {
                    local kind = mod(`j'-2,5)
                    if `kind' == 0 local expression "mod(id+`j',101)-50"
                    if `kind' == 1 local expression "mod(id+7*`j',60001)-30000"
                    if `kind' == 2 local expression "id*17+`j'"
                    if `kind' == 3 local expression "float(id/`j')"
                    if `kind' == 4 local expression "id/7+`j'/16"
                    assert v`j' == cond(mod(id,103)==0,.z,cond(mod(id,101)==0,.a,cond(mod(id,97)==0,.,`expression')))
                }
            }
            quietly datasignature
            assert "`r(datasignature)'" == "`expected_signature'"
            /* Emit only after all assertions pass. Keep rep 0 distinguishable. */
            di "COMMANDCSV,`label',`case',`n',`k',`threads',`rep',`elapsed',`load',`store'"
        }
    }
end

di "COMMANDCSV_HEADER,variant,case,rows,columns,threads,rep,elapsed_seconds,load_seconds,store_seconds"
capture noisily adaptive_commands "`label'" `reps' `threads'
local rc = _rc
di "ADAPTIVE_COMMANDS_COMPLETE RC=`rc'"
log close workload
exit `rc'
}

if "`__workload'" == "bigdata" {
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
capture log close workload
log using "`out'/`label'_`n'_`k'.log", text replace name(workload)
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
log close workload
exit
}

di as error "Unknown workload: `__workload'"
exit 198
