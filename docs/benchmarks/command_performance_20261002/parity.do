args buildpath label out
clear all
set more off
adopath ++ "`buildpath'"
local stats count mean sum min max sd variance median iqr first last firstnm lastnm p1 p5 p10 p25 p75 p90 p95 p99 skewness kurtosis
foreach size in 200 1000 {
    foreach exclude in 0 1 {
        foreach threads in 1 8 {
            clear
            set obs `=8*`size''
            gen long id=_n
            gen long g=floor((_n-1)/`size')
            gen double key=floor(mod(_n-1,`size')/10)
            gen double x=mod(_n*7919,997)/4
            replace key=.a if mod(id,131)==0
            replace x=.a if mod(id,17)==0 | g==2
            gen double shuffle=mod(id*7919,8*`size'+1)
            sort shuffle id
            drop shuffle
            local opts
            if `exclude' local opts excludeself
            local spec
            foreach stat of local stats {
                local spec `spec' (`stat') v_`stat'=x
            }
            crangestat `spec', interval(key -2 2) by(g) threads(`threads') `opts'
            if "`label'"=="baseline" save "`out'/allstats_`size'_`exclude'_`threads'.dta", replace
            else cf _all using "`out'/allstats_`size'_`exclude'_`threads'.dta"
            di "ALLSTATS_PASS,`label',`size',`exclude',`threads'"
        }
    }
}
