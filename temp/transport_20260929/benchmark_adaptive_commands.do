/* Whole-command transport acceptance benchmark.
   Args: frozen plugin/ado directory, output directory, variant label,
         repetitions (default 7), threads (default 8).
   The output directory must contain clock.plugin, built from
   validation/benchmark_clock.c. Run only through oldstata.
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
capture log close _all
log using "`out'/commands_`label'.log", text replace
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
log close
exit, clear
