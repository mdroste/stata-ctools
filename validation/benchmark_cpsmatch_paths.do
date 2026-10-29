/* One workload per process. Reloads, warmups and comparisons are untimed. */
args buildpath clockpath case N teams reps label seed reference
clear all
set more off
set linesize 255
set processors 4
adopath ++ "`buildpath'"
program path_clock, plugin using("`clockpath'")
quietly set obs `N'
set seed `seed'
set sortseed `seed'
quietly gen long row_id=_n
quietly gen double x1=rnormal()
quietly gen double x2=rnormal()
quietly gen double ps=invlogit(.6*x1+.4*x2)
quietly gen byte treat=runiform()<ps
quietly gen double y=2*treat+x1+x2+rnormal()
local opts "pscore(ps) outcome(y)"
local restriction ""
if substr("`case'",1,2)=="nr" local opts "`opts' noreplacement"
if "`case'"=="nr_equal" quietly replace ps=.5
if "`case'"=="nr_desc" local opts "`opts' descending"
if "`case'"=="nr_common" local opts "`opts' common caliper(.01)"
if substr("`case'",1,6)=="radius" local opts "`opts' radius caliper(.05)"
if "`case'"=="radius_narrow" local opts "pscore(ps) outcome(y) radius caliper(.00001)"
if "`case'"=="radius_default" local opts "pscore(ps) outcome(y) radius"
if "`case'"=="radius_gaps" {
    quietly replace ps=cond(treat,.2,.8) if mod(row_id,2)==0
    quietly replace ps=round(ps,.1)
    local opts "pscore(ps) outcome(y) radius caliper(.01) common"
}
if "`case'"=="radius_bounds" {
    quietly replace ps=mod(row_id,8)/8
    local opts "pscore(ps) outcome(y) radius caliper(.125)"
}
if "`case'"=="restricted" {
    local restriction "if mod(row_id,3)!=0 in 11/`=`N'-10'"
    quietly replace y=.a if mod(row_id,17)==0
    quietly replace ps=.b if mod(row_id,19)==0
    local opts "pscore(ps) outcome(y) radius caliper(.05)"
}
if "`case'"=="kernel" local opts "`opts' kernel bwidth(.06)"
tempfile base
quietly save `base'
foreach team of local teams {
    forvalues rep=0/`reps' {
        quietly use `base', clear
        plugin call path_clock, "start"
        quietly cpsmatch treat `restriction', `opts' threads(`team')
        plugin call path_clock, "stop"
        local wall : display %21.9f scalar(__ctools_perf_seconds)
        local results "RESULT,`label',`case',`N',`team',`rep'"
        foreach value in N n_treated n_controls n_matched n_off_support common_min common_max att att_se att_t {
            local value : display %24.16g r(`value')
            local results "`results',`=strtrim("`value'")'"
        }
        local sample "PATH,`label',`case',`N',`team',`rep',`=strtrim("`wall'")'"
        foreach stage in load setup sort match stats store total {
            capture confirm scalar _cpsmatch_time_`stage'
            if _rc local value "0"
            else local value : display %21.9f scalar(_cpsmatch_time_`stage')
            local sample "`sample',`=strtrim("`value'")'"
        }
        assert _N==`N'
        if `rep'>0 {
            di "`sample'"
            di "`results'"
        }
        if "`reference'"!="" & `rep'==1 {
            local path "`reference'/`case'_`N'_`team'.dta"
            if "`label'"=="b1" quietly save "`path'", replace
            else {
                ds _weight, not
                quietly cf `r(varlist)' using "`path'"
                preserve
                rename _weight candidate_weight
                quietly merge 1:1 row_id using "`path'", keepusing(_weight) nogen
                assert candidate_weight==_weight if _weight==0 | missing(_weight)
                assert candidate_weight>0 if _weight>0 & !missing(_weight)
                assert abs(candidate_weight-_weight)/max(abs(candidate_weight),abs(_weight))<=10^(-7) | candidate_weight==_weight
                restore
                di "PARITY,`label',`case',`N',`team',PASS"
            }
        }
    }
}
di "PATH_COMPLETE,`label',`case'"
