* Exact counts and end-to-end timings; run in separate baseline/candidate sessions.
if "`__workload'" == "" {
    gettoken __workload 0 : 0
    macro shift
}
if "`__workload'" == "cdistinct_performance" {
args buildpath clockpath label rows repetitions
if "`rows'" == "" local rows 1000000
if "`repetitions'" == "" local repetitions 5
adopath ++ "`buildpath'"
program distinct_clock, plugin using("`clockpath'")
clear
set seed 71932
quietly set obs `rows'
quietly gen long id = _n
quietly gen double low = floor(runiform()*100)
quietly gen double high = runiform()
quietly gen str16 textlow = "key" + string(low, "%03.0f")
quietly gen str24 texthigh = "key" + string(id, "%012.0f")
* Shuffle the unique strings without changing the other workloads.
quietly replace texthigh = "key" + string(mod(id*999983, `rows'), "%012.0f")
quietly gen byte selected = mod(id,100)==0
foreach case in low high textlow texthigh joint sorted selected multi {
    local keys `case'
    local options
    local restriction
    if "`case'" == "joint" {
        local keys low textlow
        local options joint
    }
    if "`case'" == "sorted" local keys id
    if "`case'" == "selected" {
        local keys low
        local restriction if selected
    }
    if "`case'" == "multi" local keys low high textlow texthigh
    foreach threads in 1 4 {
        quietly cdistinct `keys' `restriction', `options' threads(`threads')
        local expected_n = r(N)
        local expected_d = r(ndistinct)
        forvalues rep = 1/`repetitions' {
            plugin call distinct_clock, "start"
            quietly cdistinct `keys' `restriction', `options' threads(`threads')
            local n = r(N)
            local d = r(ndistinct)
            plugin call distinct_clock, "stop"
            local seconds = scalar(__ctools_perf_seconds)
            assert `n' == `expected_n' & `d' == `expected_d'
            assert id == _n
            di "DISTINCT_PERF,`label',`case',`threads',`rep',`seconds',`n',`d'"
        }
        cdistinct `keys' `restriction', `options' threads(`threads') verbose
    }
}
di "DISTINCT_PERFORMANCE_COMPLETE"
}
