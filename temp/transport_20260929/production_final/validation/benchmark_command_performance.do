* Repeated end-to-end benchmarks of one frozen build per Stata session.
* Use prepare_command_performance.py; baseline/candidate sessions run ABBA.
* A benchmark-only monotonic clock measures elapsed time.
args buildpath label outdir cases
clear all
set more off
set linesize 255
adopath ++ "`buildpath'"
ctools, version
program perf_clock, plugin using("`outdir'/clock.plugin")
set seed 95183
capture program drop perf_run
program define perf_run
    args case threads usingfile verbose
    if substr("`case'",1,2)=="ip" {
        local byopt by(g)
        if "`case'"=="ip_single" local byopt
        cipolate y x, gen(result) `byopt' epolate threads(`threads') `verbose'
    }
    else if substr("`case'",1,2)=="sp" {
        if "`case'"=="sp_multi" csplit text, gen(part) parse("::" ":" "|" ";" "--") threads(`threads') `verbose'
        else csplit text, gen(part) parse("|") threads(`threads') `verbose'
    }
    else crangejoin key -1 1 using "`usingfile'", by(g) threads(`threads') `verbose'
end
if "`cases'"=="" local cases ip_sorted ip_shuffled ip_single ip_skew ip_string sp_short sp_long sp_multi rj_sorted rj_string rj_dense rj_wide
foreach case of local cases {
    clear
    set seed 95183
    local N 1000000
    if "`case'"=="sp_long" local N 200000
    if substr("`case'",1,2)=="rj" local N 300000
    if "`case'"=="rj_dense" local N 100000
    quietly set obs `N'
    quietly gen long id=_n
    if substr("`case'",1,2)=="ip" {
        quietly gen long g=floor((_n-1)/1000)
        quietly gen double x=mod(_n-1,1000)
        if "`case'"=="ip_shuffled" | "`case'"=="ip_string" {
            quietly replace g=floor(runiform()*1000)
            quietly replace x=runiform()*1000
        }
        if "`case'"=="ip_single" quietly replace x=runiform()*1000000
        if "`case'"=="ip_skew" {
            quietly replace g=0 if _n<=950000
            quietly replace x=runiform()*1000000
        }
        if "`case'"=="ip_string" {
            rename g gn
            quietly gen str12 g="group"+string(gn,"%05.0f")
            drop gn
        }
        quietly gen double y=rnormal()
        quietly replace y=. if runiform()<.6
    }
    else if substr("`case'",1,2)=="sp" {
        if "`case'"=="sp_short" quietly gen str80 text="NY|"+string(mod(_n,1000))+"|business|active|2024-01-31"
        if "`case'"=="sp_long" {
            local longtext
            forvalues k=1/16 {
                local longtext "`longtext'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuv|"
            }
            quietly gen str1800 text="`longtext'end"
        }
        if "`case'"=="sp_multi" quietly gen str180 text="NY::"+string(mod(_n,1000))+"|business--active;2024::01|31;end"
    }
    else {
        quietly gen long g=floor((_n-1)/200)
        quietly gen double key=mod(_n-1,200)
        if "`case'"=="rj_string" {
            quietly replace g=floor(runiform()*1500)
            quietly replace key=floor(runiform()*200)
            rename g gn
            quietly gen str12 g="group"+string(gn,"%05.0f")
            drop gn
        }
        if "`case'"=="rj_dense" quietly replace key=key/10
        quietly gen str20 text="row"+string(_n)
        if "`case'"=="rj_wide" {
            forvalues k=1/12 {
                quietly gen double v`k'=runiform()
            }
        }
    }
    tempfile base usingfile
    if substr("`case'",1,2)=="rj" {
        quietly save `usingfile'
        quietly replace key=key+.5
    }
    quietly save `base'
    quietly perf_run `case' 12 "`usingfile'" ""
    if substr("`label'",1,1)=="b" quietly save "`outdir'/expected_`case'.dta", replace
    else cf _all using "`outdir'/expected_`case'.dta"
    quietly use `base', clear
    di "PROFILE_BEGIN label=`label' case=`case'"
    perf_run `case' 12 "`usingfile'" verbose
    di "PROFILE_END label=`label' case=`case'"
    foreach threads in 1 4 12 {
        forvalues rep=0/6 {
            quietly use `base', clear
            plugin call perf_clock, "start"
            quietly perf_run `case' `threads' "`usingfile'" ""
            plugin call perf_clock, "stop"
            if `rep'>0 di "PERF,`label',`case',`threads',`rep'," %12.6f scalar(__ctools_perf_seconds) "," _N
        }
    }
}
di "PERFORMANCE_COMPLETE label=`label'"
