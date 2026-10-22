* Reproducible end-to-end audit: all commands on 40 mixed-type columns.
* run_performance_audit.py freezes the plugin/ado files and validates completion.
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
        local c_ppmlhdfe "cpplmhdfe count x1 x2 x3, absorb(firm group) tolerance(1e-8)"
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
