* Frozen public commands; setup, comparisons and disk reads are outside timers.
args buildpath label outdir cases reps threadlist rows refreps
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
    if "`case'" == "e_long" quietly replace x = x
    if "`case'" == "e_long" {
        quietly recast str2045 x
        quietly replace x = x + dup(2034,"a")
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
        drop seed
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
                if "`label'" == "b1" quietly save "`outdir'/expected_`case'_`threads'.dta", replace
                else cf _all using "`outdir'/expected_`case'_`threads'.dta"
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
