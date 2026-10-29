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
