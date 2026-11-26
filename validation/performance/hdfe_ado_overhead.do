* Time one frozen ado/plugin directory per Stata process. Run again in a new
* process for the other directory, using the same plugin and arguments.
* Usage: do validation/performance/hdfe_ado_overhead.do directory label [groups] [reps] [commands]
* Defaults reproduce the 10M-row, 1M-individual balanced panel.
args directory variant groups reps commands
if `"`directory'"' == "" | "`variant'" == "" error 198
if "`groups'" == "" local groups 1000000
if "`reps'" == "" local reps 5
if "`commands'" == "" local commands creghdfe civreghdfe cppmlhdfe
clear all
set more off
adopath ++ `"`directory'"'
set seed 20261125
set obs `groups'
generate i = _n
expand 10
bysort i: generate t = _n
generate x1 = rnormal()
generate x2 = rnormal()
generate x3 = rnormal()
generate y = x1 + x2 + x3 + i/`groups' + t/10 + rnormal()

foreach cmd of local commands {
    local model y x1 x2 x3
    if "`cmd'" == "civreghdfe" {
        set seed 20261126
        generate z = x1 + rnormal()
        local model y x2 x3 (x1 = z)
    }
    if "`cmd'" == "cppmlhdfe" {
        set seed 20261127
        replace y = rpoisson(exp(.1*x1 + .1*x2 + .1*x3 + i/`groups' + t/10))
    }
    * Warm the estimator and its display routines before measuring.
    `cmd' `model', absorb(i t)
    forvalues pair = 1/`reps' {
        di "ADO_BENCH_BEGIN command=`cmd' variant=`variant' pair=`pair' rows=" _N
        `cmd' `model', absorb(i t) verbose
        di "ADO_BENCH_END command=`cmd' variant=`variant' pair=`pair'"
    }
}
di "ADO_BENCH_COMPLETE"
