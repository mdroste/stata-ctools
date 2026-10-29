* End-to-end comparisons with both winsor2 and gstats, one frozen ctools build.
args buildpath clockpath rows reps
clear all
set more off
set linesize 255
adopath ++ "`buildpath'"
ctools, version
which winsor2
which gstats
program perf_clock, plugin using("`clockpath'")
foreach case in grouped single wide noby {
    clear
    set seed 928143
    quietly set obs `rows'
    quietly gen long id=_n
    quietly gen long g=floor((_n-1)/200)
    if "`case'"=="single" quietly replace g=0
    quietly gen double x=100+5*rnormal()
    quietly replace x=. if mod(id,37)==0
    quietly replace x=.a if mod(id,107)==0
    local vars x
    if "`case'"=="wide" {
        forvalues k=2/8 {
            quietly gen double x`k'=100+5*rnormal()
            local vars `vars' x`k'
        }
    }
    quietly gen double shuffle=runiform()
    sort shuffle
    drop shuffle
    local byopt by(g)
    if "`case'"=="noby" local byopt
    tempfile base expected
    quietly save `base'
    foreach operation in ctools gstats winsor2 {
        local command cwinsor `vars', `byopt' cuts(1 99) suffix(_w) threads(8)
        if "`operation'"=="gstats" local command gstats winsor `vars', `byopt' cuts(1 99) suffix(_w)
        if "`operation'"=="winsor2" local command winsor2 `vars', `byopt' cuts(1 99) suffix(_w)
        forvalues rep=0/`reps' {
            quietly use `base', clear
            plugin call perf_clock, "start"
            quietly `command'
            plugin call perf_clock, "stop"
            di "WINSOR_REF,`case',`operation',`rows',`rep'," %12.7f scalar(__ctools_perf_seconds)
            if `rep'==0 {
                foreach v of local vars {
                    * winsor2 suffix() turns extended missings into ordinary .;
                    * compare all missing locations and all finite values.
                    quietly replace `v'_w=. if missing(`v'_w)
                }
                if "`operation'"=="ctools" quietly save `expected'
                else {
                    foreach v of local vars {
                        rename `v'_w `v'_ref
                    }
                    quietly merge 1:1 id using `expected', assert(match) nogenerate
                    foreach v of local vars {
                        assert missing(`v'_w)==missing(`v'_ref)
                        assert `v'_w==`v'_ref | abs(`v'_w-`v'_ref)/max(abs(`v'_w),abs(`v'_ref))<1e-7 if !missing(`v'_w)
                    }
                }
            }
        }
    }
}
di "WINSOR_REFERENCES_COMPLETE"
