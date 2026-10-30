* Named benchmark workloads. All original sizes, seeds and timing checks are retained.
if "`__workload'" == "" {
    gettoken __workload 0 : 0
    macro shift
}

if "`__workload'" == "csort_performance" {
/* Focused csort benchmark; prepare_csort_performance.py supplies frozen builds.
   Each case starts from deterministic, unmarked data. Native stable sorting,
   reloads, warmups, and exact whole-dataset comparisons are outside timing. */
args buildpath label outdir cases threads reps datadir keepdata
clear all
set more off
set linesize 255
adopath ++ "`buildpath'"
which csort
ctools, version
di "CSORT_SETTINGS label=`label' stata=`c(stata_version)' processors=`c(processors)'"
program csort_clock, plugin using("`outdir'/clock.plugin")

foreach spec of local cases {
    local colon = strpos("`spec'", ":")
    local case = substr("`spec'", 1, `colon'-1)
    local N = real(substr("`spec'", `colon'+1, .))
    clear
    set seed 95183
    set sortseed 95183
    quietly set obs `N'
    quietly gen long row_id = _n
    local keys "key"
    local K 4
    if "`case'" == "wide" local K 32
    if "`case'" == "verywide" local K 256
    if inlist("`case'", "ties", "constant_ties") local K 5

    if inlist("`case'", "ties", "constant_ties") {
        /* Four equal-valued key columns; row_id witnesses stability within ties. */
        forvalues j = 1/4 {
            quietly gen byte key`j' = floor((row_id-1)/ceil(`N'/8))
            if "`case'" == "constant_ties" quietly replace key`j' = 0
        }
        local keys "key1 key2 key3 key4"
    }
    else if "`case'" == "panel" {
        quietly gen long firm = floor(runiform()*10000)
        quietly gen int period = floor(runiform()*100)
        local keys "firm period"
    }
    else if "`case'" == "mixed" {
        quietly gen str8 region = "G" + string(floor(runiform()*10000), "%07.0f")
        quietly gen double key = rnormal()
        local keys "region key"
    }
    else if inlist("`case'", "long_string", "long_string_large", "long_prefix") {
        quietly gen str64 key = ""
        local chunks 8
        if "`case'" == "long_prefix" {
            quietly replace key = "0123456789012345678901234567890123456789"
            local chunks 3
        }
        forvalues j = 1/`chunks' {
            quietly replace key = key + string(floor(runiform()*100000000), "%08.0f")
        }
    }
    else if inlist("`case'", "nearly_sorted", "nearly_sorted_middle") {
        quietly gen int key = floor((_n-1)/max(1,ceil(`N'/1000)))
        /* One adjacent inversion, at the end or away from regular partitions. */
        local inversion = `N'
        if "`case'" == "nearly_sorted_middle" local inversion = max(2,floor(`N'*.61803398875))
        if `N'>1 quietly replace key = key[`inversion'-1]-1 in `inversion'
    }
    else if "`case'" == "integer" quietly gen int key = floor(runiform()*1000)
    else if "`case'" == "sorted" quietly gen long key = floor((_n-1)/10)
    else if "`case'" == "reverse" quietly gen long key = floor((`N'-_n)/10)
    else quietly gen double key = rnormal()

    local first = c(k)+1
    if `first'<=`K' {
        forvalues j = `first'/`K' {
            quietly gen double value`j' = row_id/7 + `j'/16
            quietly replace value`j' = .a if mod(row_id+`j',101)==0
        }
    }
    if "`case'" == "verywide" order key, last
    assert c(k)==`K'
    local marked : sortedby
    assert "`marked'" == ""
    local base "`datadir'/`case'_base.dta"
    local expected "`datadir'/`case'_expected.dta"
    quietly save "`base'", replace
    quietly sort `keys', stable
    quietly save "`expected'", replace

    foreach team of local threads {
        /* Validate the plugin's permutation without the ado's final Stata sort. */
        quietly use "`base'", clear
        quietly csort `keys', threads(`team') nosortedby
        assert _N==`N' & c(k)==`K'
        quietly cf _all using "`expected'"
        di "CSORT_CHECK,`label',`case',`team',nosortedby,PASS"

        /* Rep 0 warms the complete public command and is never reported as a sample. */
        forvalues rep = 0/`reps' {
            quietly use "`base'", clear
            plugin call csort_clock, "start"
            quietly csort `keys', threads(`team')
            plugin call csort_clock, "stop"
            local wall : display %21.9f scalar(__ctools_perf_seconds)
            local wall = strtrim("`wall'")
            local record "CSORT_PERF,`label',`case',`N',`K',`team',`rep',`wall'"
            foreach stage in load sort permute store stream cleanup total {
                local value : display %21.9f scalar(_csort_time_`stage')
                local value = strtrim("`value'")
                local record "`record',`value'"
            }
            local record "`record',`=scalar(_csort_stream)',`=scalar(_csort_threads_max)'"
            assert _N==`N' & c(k)==`K'
            local marked : sortedby
            assert "`marked'" == "`keys'"
            quietly cf _all using "`expected'"
            if `rep'>0 di "`record'"
        }
    }
    if !`keepdata' {
        erase "`base'"
        erase "`expected'"
    }
}
di "CSORT_PERFORMANCE_COMPLETE label=`label'"
exit
}

di as error "Unknown workload: `__workload'"
exit 198
