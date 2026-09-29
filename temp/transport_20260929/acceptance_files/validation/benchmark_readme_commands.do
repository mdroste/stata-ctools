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
