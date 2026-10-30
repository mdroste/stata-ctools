* Named benchmark workloads. All original sizes, seeds and timing checks are retained.
if "`__workload'" == "" {
    gettoken __workload 0 : 0
    macro shift
}

if "`__workload'" == "command_performance" {
* Repeated end-to-end benchmarks of one frozen build per Stata session.
* Use prepare_command_performance.py; baseline/candidate sessions run ABBA.
* A benchmark-only monotonic clock measures elapsed time.
args buildpath label outdir cases datadir
if `"`datadir'"'=="" local datadir `"`outdir'"'
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
    if substr("`label'",1,1)=="b" quietly save "`datadir'/expected_`case'.dta", replace
    else cf _all using "`datadir'/expected_`case'.dta"
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
exit
}

if "`__workload'" == "newcommands" {
* Run via a Stata driver that captures errors and exits cleanly.
* Reference rangejoin must be on adopath. Timings are end-to-end, one run each.
clear all
set more off
adopath ++ "build"
ctools, version
creturn list
which ipolate
which split
which rangejoin
local N = 1000000
set seed 95183
clear
set obs `N'
gen long id=_n
gen int g=mod(_n,1000)
gen double x=floor((_n-1)/1000)
gen double y=rnormal()
replace y=. if runiform()<.6
tempfile base
save `base'
timer clear
timer on 1
ipolate y x, gen(ref) by(g) epolate
timer off 1
sort id
tempfile interpolated
save `interpolated'
use `base', clear
timer on 2
cipolate y x, gen(ref) by(g) epolate verbose
timer off 2
cf _all using `interpolated'
clear
set obs `N'
gen str80 text="NY|123|business|active|2024-01-31"
timer on 3
split text, gen(ref) parse("|")
timer off 3
timer on 4
csplit text, gen(got) parse("|") verbose
timer off 4
forvalues j=1/5 {
    assert ref`j'==got`j'
}
clear
set obs 200000
gen long id=_n
gen long g=floor((_n-1)/200)
gen double key=mod(_n-1,200)
gen str20 text="row"+string(_n)
tempfile using master joined
save `using'
replace key=key+.5
save `master'
timer on 5
rangejoin key -1 1 using `using', by(g)
timer off 5
save `joined'
use `master', clear
timer on 6
crangejoin key -1 1 using `using', by(g) verbose
timer off 6
cf _all using `joined'
quietly timer list 1
local ipolate=r(t1)
quietly timer list 2
local cipolate=r(t2)
quietly timer list 3
local split=r(t3)
quietly timer list 4
local csplit=r(t4)
quietly timer list 5
local rangejoin=r(t5)
quietly timer list 6
local crangejoin=r(t6)
di "BENCH ipolate rows=`N' native=`ipolate' ctools=`cipolate'"
di "BENCH split rows=`N' native=`split' ctools=`csplit'"
di "BENCH rangejoin rows=200000+200000 output=`=_N' native=`rangejoin' ctools=`crangejoin'"
exit
}

if "`__workload'" == "audit" {
* Reproducible end-to-end audit: all commands on 40 mixed-type columns.
* validation/bench.py audit freezes the plugin/ado files and validates completion.
args buildpath outdir label sizes threads reps cases
clear all
set more off
set linesize 255
set graphics off
adopath ++ "`buildpath'"
ctools, version
program audit_clock, plugin using("`outdir'/clock.plugin")
di "AUDIT_SETTINGS label=`label' sizes=`sizes' columns=40 threads=`threads' reps=`reps' stata=`c(stata_version)' processors=`c(processors)'"
if `"`cases'"'=="" local cases export import sort merge sample bsample encode decode destring ipolate split rangejoin winsor rangestat binscatter reghdfe ivreghdfe ppmlhdfe psmatch2 qreg group_numeric group_mixed group_count group_label winsor_by reghdfe_strcluster sample_by bsample_by
foreach N of local sizes {
        clear
        set seed 95183
        set sortseed 95183
        quietly set obs `N'
        * Base columns: 16 numeric (all five storage types), four strings.
        quietly gen long id = _n
        quietly gen long firm = mod(id-1,`N'/200)+1
        quietly gen int month = floor((id-1)/(`N'/200))
        quietly gen int group = floor(runiform()*100)+1
        quietly gen double z = rnormal()
        quietly gen double x1 = .7*z + rnormal()
        quietly gen float x2 = rnormal()
        quietly gen float x3 = rnormal()
        quietly gen double y = 1 + x1 + .5*x2 - .2*x3 + .1*sin(firm) + rnormal()
        quietly gen double yi = cond(runiform()<.3,.,y)
        quietly gen long count = rpoisson(exp(.3+.1*x1+.1*x2+.1*sin(firm)))
        quietly gen byte treat = runiform()<invlogit(.2*x1-.1*x2)
        quietly gen byte category = mod(id,4)+1
        label define category 1 "North" 2 "South" 3 "East" 4 "West", replace
        label values category category
        quietly gen str8 region = cond(category==1,"North",cond(category==2,"South",cond(category==3,"East","West")))
        quietly gen str8 sector = "S" + string(mod(id*73,1000),"%04.0f")
        quietly gen str16 number = string(mod(id*37,1000000)/100,"%12.2f")
        quietly gen str40 code = region + "|" + sector + "|" + strtrim(number)
        quietly gen float noise = rnormal()
        quietly gen int year = 2000 + mod(month,20)
        quietly gen byte flag = mod(id,2)
        forvalues j=1/4 {
            quietly gen double extra_d`j' = rnormal()
            quietly gen float extra_f`j' = rnormal()
            quietly gen int extra_i`j' = mod(id*`j',20000)
            quietly gen byte extra_b`j' = mod(id,100)
            quietly gen str24 extra_s`j' = "payload" + string(mod(id*`j',10000),"%05.0f")
        }
        assert c(k)==40
        tempfile base lookup offers inputcsv outputcsv
        quietly save `base'
        preserve
        quietly keep if month==0
        keep firm noise region
        rename noise lookup_value
        rename region lookup_region
        quietly save `lookup'
        restore
        * Forty-variable using input; a narrow interval yields one match per row.
        preserve
        quietly recast double month
        quietly replace month = month+.25
        quietly save `offers'
        restore

        local c_export `"cexport delimited using "`outputcsv'", replace nolabel"'
        local c_import `"cimport delimited using "`inputcsv'", clear"'
        local c_sort "csort x1 id,"
        local c_merge `"cmerge m:1 firm using "`lookup'", nogen"'
        local c_sample "csample 10,"
        local c_bsample "cbsample,"
        local c_encode "cencode region, gen(result)"
        local c_decode "cdecode category, gen(result)"
        local c_destring "cdestring number, gen(result)"
        local c_ipolate "cipolate yi month, by(firm) gen(result) epolate"
        local c_split `"csplit code, gen(part) parse("|")"'
        local c_rangejoin `"crangejoin month -.5 .5 using "`offers'", by(firm) keepusing(x2 region)"'
        local c_qreg "cqreg y x1 x2 x3,"
        local c_winsor "cwinsor x1 x2 x3, cuts(1 99) replace"
        local c_rangestat "crangestat (mean) y, interval(month -2 2) by(firm)"
        local c_psmatch2 "cpsmatch treat x1 x2, outcome(y) logit neighbor(1)"
        local c_binscatter "cbinscatter y x1, nquantiles(20)"
        local c_reghdfe "creghdfe y x1 x2 x3, absorb(firm group) tolerance(1e-8)"
        local c_ivreghdfe "civreghdfe y x2 x3 (x1=z), absorb(firm group) tolerance(1e-8)"
        local c_ppmlhdfe "cppmlhdfe count x1 x2 x3, absorb(firm group) tolerance(1e-8)"
        local c_group_numeric "_ctools_group firm, generate(result)"
        local c_group_mixed "_ctools_group sector group, generate(result)"
        local c_group_count "_ctools_group sector group, count"
        local c_group_label "_ctools_group region, generate(result) label"
        local c_winsor_by "cwinsor x1 x2 x3 extra_d*, cuts(1 99) by(sector) suffix(_w)"
        local c_reghdfe_strcluster "creghdfe y x1 x2 x3, absorb(firm group) vce(cluster sector)"
        local c_sample_by "csample 10, by(sector)"
        local c_bsample_by "cbsample, strata(sector)"

        * Input generation/export is outside the timed region.
        quietly cexport delimited using "`inputcsv'", replace nolabel
        foreach cmd of local cases {
            forvalues rep=0/`reps' {
                quietly use `base', clear
                set seed 95183
                set sortseed 95183
                local command `"`c_`cmd''"'
                if substr("`cmd'",1,6)!="group_" local command `"`command' threads(`threads')"'
                if `rep'==0 & substr("`cmd'",1,6)!="group_" local command `"`command' verbose"'
                di "AUDIT_START,`label',`N',`cmd',`threads',`rep'"
                plugin call audit_clock, "start"
                capture noisily `command'
                local rc = _rc
                local groups .
                if substr("`cmd'",1,6)=="group_" & !`rc' local groups = r(N_groups)
                plugin call audit_clock, "stop"
                local elapsed = scalar(__ctools_perf_seconds)
                di "AUDIT,`label',`N',`cmd',`threads',`rep'," %12.6f `elapsed' ",`rc'," _N "," c(k)
                if `rc' exit `rc'
                if `rep'==0 {
                    local expected = `N'
                    if inlist("`cmd'","sample","sample_by") local expected=`N'/10
                    assert _N==`expected'
                    if inlist("`cmd'","qreg","reghdfe","reghdfe_strcluster","ivreghdfe","ppmlhdfe") {
                        assert e(N)==`N'
                        di "AUDIT_ESTIMATE,`label',`N',`cmd',`threads'"
                        matrix list e(b), format(%21x)
                        matrix list e(V), format(%21x)
                    }
                    quietly datasignature
                    di "AUDIT_CHECK,`label',`N',`cmd',`threads',`r(datasignature)',`groups'"
                }
            }
        }
        erase `base'
        erase `lookup'
        erase `offers'
        capture erase `inputcsv'
        capture erase `outputcsv'
}
di "AUDIT_COMPLETE label=`label'"
exit
}

di as error "Unknown workload: `__workload'"
exit 198
