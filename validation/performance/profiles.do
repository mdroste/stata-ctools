* Named benchmark workloads. All original sizes, seeds and timing checks are retained.
if "`__workload'" == "" {
    gettoken __workload 0 : 0
    macro shift
}

if "`__workload'" == "profile_readme" {
* One-shot README profiler. Run from the repository root; see docs/performance/README.md.
* Arguments: output directory, quoted observation counts, Stata/MP processors,
* quoted subset of pairs (default all). ctools commands take no threads()
* option, so they use the ctools default.
args outdir sizes threads pairs
clear all
set more off
set linesize 255
set graphics off
if `"`outdir'"' == "" local outdir "docs/benchmarks/readme_profile"
if `"`sizes'"' == "" local sizes "1000000 10000000 25000000"
if "`threads'" == "" local threads = c(processors)
local allpairs "export import sort merge sample bsample encode decode destring ipolate split rangejoin winsor rangestat binscatter reghdfe ivreghdfe ppmlhdfe psmatch2 qreg"
if `"`pairs'"' == "" local pairs "`allpairs'"
capture mkdir "`outdir'"
log using "`outdir'/profile.log", text replace name(profile)

program define profile_readme
    args outdir sizes threads pairs allpairs
    set processors `threads'
    if `"`: list pairs - allpairs'"' != "" | `"`: list dups pairs'"' != "" error 198
    adopath ++ "build"
    local candidate : environment CTOOLS_TEST_PLUGIN_DIR
    if `"`candidate'"' != "" adopath ++ `"`candidate'"'
    about
    ctools, version
    di "PROFILE_SETTINGS sizes=`sizes' variables=20 processors=`threads' seed=95183 pairs=`pairs'"
    foreach cmd in cimport cexport csort cmerge csample cbsample cencode cdecode ///
        cdestring cipolate csplit crangejoin cqreg cwinsor crangestat cpsmatch ///
        cbinscatter creghdfe civreghdfe cppmlhdfe rangejoin rangestat gstats ///
        psmatch2 binscatter reghdfe ivreghdfe ppmlhdfe ftools ivreg2 {
        capture noisily which `cmd'
    }
    tempname results
    file open `results' using "`outdir'/timings.csv", write text replace
    file write `results' "command,observations,method,seconds,rc,check_rc" _n
    foreach N of local sizes {
        if `N' < 1000 | mod(`N',200) != 0 error 198
        clear
        set seed 95183
        set sortseed 95183
        quietly set obs `N'
        * Exactly 20 variables: 16 numeric (all five storage types), four strings.
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
        assert c(k)==20
        tempfile base lookup offers inputcsv outputcsv
        quietly save `base'
        preserve
        quietly keep if month==0
        keep firm noise region
        rename noise lookup_value
        rename region lookup_region
        quietly save `lookup'
        restore
        * Twenty-variable using input; a narrow interval yields one match per row.
        preserve
        quietly recast double month
        quietly replace month = month+.25
        quietly save `offers'
        restore

        local ref_export `"export delimited using "`inputcsv'", replace nolabel"'
        local c_export `"cexport delimited using "`outputcsv'", replace nolabel"'
        local ref_import `"import delimited using "`inputcsv'", clear"'
        local c_import `"cimport delimited using "`inputcsv'", clear"'
        local ref_sort "sort x1 id, stable"
        local c_sort "csort x1 id"
        local ref_merge `"merge m:1 firm using "`lookup'", nogen"'
        local c_merge `"cmerge m:1 firm using "`lookup'", nogen"'
        local ref_sample "sample 10"
        local c_sample "csample 10"
        local ref_bsample "bsample"
        local c_bsample "cbsample"
        local ref_encode "encode region, gen(result)"
        local c_encode "cencode region, gen(result)"
        local ref_decode "decode category, gen(result)"
        local c_decode "cdecode category, gen(result)"
        local ref_destring "destring number, gen(result)"
        local c_destring "cdestring number, gen(result)"
        local ref_ipolate "ipolate yi month, by(firm) gen(result) epolate"
        local c_ipolate "cipolate yi month, by(firm) gen(result) epolate"
        local ref_split `"split code, gen(part) parse("|")"'
        local c_split `"csplit code, gen(part) parse("|")"'
        local ref_rangejoin `"rangejoin month -.5 .5 using "`offers'", by(firm) keepusing(x2 region)"'
        local c_rangejoin `"crangejoin month -.5 .5 using "`offers'", by(firm) keepusing(x2 region)"'
        local ref_qreg "qreg y x1 x2 x3"
        local c_qreg "cqreg y x1 x2 x3"
        local ref_winsor "gstats winsor x1 x2 x3, cuts(1 99) replace"
        local c_winsor "cwinsor x1 x2 x3, cuts(1 99) replace"
        local ref_rangestat "rangestat (mean) y, interval(month -2 2) by(firm)"
        local c_rangestat "crangestat (mean) y, interval(month -2 2) by(firm)"
        local ref_psmatch2 "psmatch2 treat x1 x2, outcome(y) logit neighbor(1)"
        local c_psmatch2 "cpsmatch treat x1 x2, outcome(y) logit neighbor(1)"
        local ref_binscatter "binscatter y x1, nquantiles(20)"
        local c_binscatter "cbinscatter y x1, nquantiles(20)"
        local ref_reghdfe "reghdfe y x1 x2, absorb(firm group) tolerance(1e-8)"
        local c_reghdfe "creghdfe y x1 x2, absorb(firm group) tolerance(1e-8)"
        local ref_ivreghdfe "ivreghdfe y x2 (x1=z), absorb(firm group) tolerance(1e-8)"
        local c_ivreghdfe "civreghdfe y x2 (x1=z), absorb(firm group) tolerance(1e-8)"
        local ref_ppmlhdfe "ppmlhdfe count x1 x2, absorb(firm group) tolerance(1e-8)"
        local c_ppmlhdfe "cppmlhdfe count x1 x2, absorb(firm group) tolerance(1e-8)"

        foreach cmd of local pairs {
            foreach method in ref c {
                quietly use `base', clear
                set seed 95183
                set sortseed 95183
                local command `"``method'_`cmd''"'
                di `"PROFILE_START N=`N' method=`method' command=`command'"'
                timer clear 99
                timer on 99
                capture noisily `command'
                local rc = _rc
                timer off 99
                quietly timer list 99
                local elapsed = r(t99)
                local checkrc = 0
                if `rc'==0 {
                    local expected = `N'
                    if "`cmd'"=="sample" local expected = `N'/10
                    capture assert _N==`expected'
                    local checkrc = _rc
                    if inlist("`cmd'","qreg","reghdfe","ivreghdfe","ppmlhdfe") {
                        capture assert e(N)==`N'
                        if _rc local checkrc = _rc
                    }
                    if "`cmd'"=="import" {
                        capture assert c(k)==20
                        if _rc local checkrc = _rc
                    }
                }
                file write `results' "`cmd',`N',`method',`elapsed',`rc',`checkrc'" _n
                di "PROFILE_RESULT `cmd',`N',`method',`elapsed',`rc',`checkrc'"
            }
        }
        * Release disk space before generating the next scale.
        erase `base'
        erase `lookup'
        erase `offers'
        capture erase `inputcsv'
        capture erase `outputcsv'
    }
    file close `results'
    di "PROFILE_COMPLETE"
end
capture noisily profile_readme `"`outdir'"' `"`sizes'"' `threads' `"`pairs'"' `"`allpairs'"'
local rc = _rc
di "PROFILE_EXIT rc=`rc'"
log close profile
exit `rc'
}

if "`__workload'" == "readme_commands" {
* Native/SSC comparisons for the README. Prepared by prepare_readme_benchmarks.py.
* Generation, input reloads, and exact comparisons are outside the timed region.
args buildpath refpath outdir N threads reps
clear all
set more off
set linesize 255
set processors `threads'
set seed 95183
set sortseed 95183
adopath ++ "`refpath'"
adopath ++ "`buildpath'"
ctools, version
about
which ipolate
which split
which rangejoin
which rangestat
program perf_clock, plugin using("`outdir'/clock.plugin")
program define readme_run
    args case method threads usingfile
    local cap
    if "`method'"=="ctools" local cap threads(`threads')
    if "`case'"=="ipolate" {
        local command ipolate
        if "`method'"=="ctools" local command cipolate
        `command' sales month, generate(sales_i) by(firm) epolate `cap'
    }
    else if "`case'"=="split" {
        local command split
        if "`method'"=="ctools" local command csplit
        `command' code, generate(part) parse("|") `cap'
    }
    else {
        local command rangejoin
        if "`method'"=="ctools" local command crangejoin
        `command' date -1 1 using "`usingfile'", by(firm) keepusing(offer_id price venue) `cap'
    }
end

foreach case in ipolate split rangejoin {
    clear
    set seed 95183
    quietly set obs `N'
    quietly gen long row_id=_n
    local expected_rows `N'
    tempfile base offers expected
    if "`case'"=="ipolate" {
        * 200 monthly extracts, interleaved firms, 35% missing sales.
        quietly gen long firm=mod(row_id-1,`N'/200)+1
        quietly gen int month=ym(2000,1)+floor((row_id-1)/(`N'/200))
        format month %tm
        quietly gen byte sector=1+mod(firm,20)
        quietly gen int employment=20+mod(firm*73,1980)
        quietly gen str2 region=substr("NYSFTXFL",1+2*mod(firm,4),2)
        quietly gen double sales=exp(8+ln(employment)+.002*(month-ym(2000,1))+.15*rnormal())
        quietly replace sales=. if runiform()<.35
    }
    else if "`case'"=="split" {
        * Six heterogeneous fields in transaction records, not a constant string.
        quietly gen str64 code=substr("NYSFTXFL",1+2*mod(row_id,4),2)+"|"+ ///
            cond(mod(row_id,3)==0,"wholesale",cond(mod(row_id,3)==1,"retail","online"))+"|SKU"+ ///
            string(mod(row_id*73,1000000),"%06.0f")+"|"+string(mod(row_id*37,100000000),"%08.0f")+ ///
            "|"+cond(mod(row_id,2)==0,"USD","EUR")+"|"+string(td(01jan2020)+mod(row_id,1461),"%tdCCYY-NN-DD")
        quietly gen int quantity=1+mod(row_id,50)
        quietly gen double amount=quantity*(10+mod(row_id*17,10000)/100)
    }
    else {
        * N offers: 500 daily quotes per firm, numeric and string payloads.
        rename row_id offer_id
        quietly gen long firm=floor((offer_id-1)/500)+1
        quietly gen int date=td(01jan2020)+mod(offer_id-1,500)
        format date %td
        quietly gen double price=100+mod(firm*13,10000)/100+.01*mod(offer_id-1,500)
        quietly gen str3 venue=cond(mod(offer_id,2)==0,"NYC","CHI")
        quietly save `offers'
        * N inquiries: two adjacent daily offers (one at each firm's last date).
        clear
        quietly set obs `N'
        quietly gen long row_id=_n
        quietly gen long firm=floor((row_id-1)/500)+1
        quietly gen double date=td(01jan2020)+mod(row_id-1,500)+.5
        format date %td
        quietly gen int quantity=1+mod(row_id*17,50)
        quietly gen str2 region=substr("NYSFTXFL",1+2*mod(firm,4),2)
        local expected_rows=2*`N'-`N'/500
    }
    quietly save `base'
    * Full-size warmups also verify every output cell against the reference.
    foreach method in reference ctools {
        quietly use `base', clear
        quietly readme_run `case' `method' `threads' "`offers'"
        assert _N==`expected_rows'
        if "`case'"=="rangejoin" sort row_id offer_id
        else sort row_id
        if "`method'"=="reference" quietly save `expected'
        else cf _all using `expected'
    }
    di "README_CHECK,`case',`N',PASS,`expected_rows'"
    erase `expected'
    * Four repeats by default: AB, BA, AB, BA to balance ordering.
    forvalues rep=1/`reps' {
        local methods reference ctools
        if mod(`rep',2)==0 local methods ctools reference
        foreach method of local methods {
            quietly use `base', clear
            plugin call perf_clock, "start"
            quietly readme_run `case' `method' `threads' "`offers'"
            plugin call perf_clock, "stop"
            assert _N==`expected_rows'
            di "README_BENCH,`case',`N',`method',`threads',`rep'," %12.6f scalar(__ctools_perf_seconds) "," _N
        }
    }
    erase `base'
    if "`case'"=="rangejoin" erase `offers'
}
di "README_BENCHMARK_COMPLETE N=`N'"
exit
}

di as error "Unknown workload: `__workload'"
exit 198
