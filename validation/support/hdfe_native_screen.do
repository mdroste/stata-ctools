* Compare native screening to Stata marksample, including parallel-sized
* inputs. Use identical new plugins; only the candidate ado defers marking.
args baseline candidate
if `"`baseline'"' == "" | `"`candidate'"' == "" error 198
clear all
set more off
set seed 20261129
set obs 131103
generate long g = ceil(_n/9)
generate byte t = mod(_n-1,9)+1
generate double x1 = rnormal()
generate double x2 = rnormal()
generate double z = x1 + rnormal()
generate double y = x1+x2+g/15000+t/9+rnormal()
generate double yp = rpoisson(exp(.1*x1+.1*x2+g/15000+t/9))
generate byte expected = 0
generate str20 cluster_string = ""
foreach v in g t x1 x2 z y yp {
    clonevar original_`v' = `v'
}
local checked = 0
foreach pattern in complete sparse dense {
    foreach v in g t x1 x2 z y yp {
        quietly replace `v' = original_`v'
    }
    if "`pattern'" != "complete" {
        quietly replace x1 = .a in 4095/4097
        quietly replace x2 = .z in 8192
        quietly replace g = . in 16384
        quietly replace t = .b in 32768
        quietly replace y = . in 65536
        quietly replace yp = .c in 65537
        quietly replace z = .d in 98304
        quietly replace g = 20000 in 100000
    }
    if "`pattern'" == "dense" {
        quietly replace x2 = .a if mod(_n,7)==0
        quietly replace g = .z if mod(_n,11)==0
        quietly replace yp = 0 if g == 2
    }
    quietly replace cluster_string = string(original_g)
    quietly replace cluster_string = "excluded_" + string(_n) if missing(x1,x2,g,t)
    foreach cmd in creghdfe civreghdfe cppmlhdfe {
        local model y x1 x2
        if "`cmd'" == "civreghdfe" local model y x2 (x1=z)
        if "`cmd'" == "cppmlhdfe" local model yp x1 x2
        foreach spec in "absorb(g t)" "absorb(g t) vce(cluster g) keepsingletons" "absorb(g t) vce(cluster cluster_string)" {
            foreach qualifier in "" "if t > 2 in 2/130000" {
                foreach variant in baseline candidate {
                    quietly capture adopath - `"`baseline'"'
                    quietly capture adopath - `"`candidate'"'
                    quietly adopath ++ `"``variant''"'
                    discard
                    quietly `cmd' `model' `qualifier', `spec'
                    if "`variant'" == "baseline" {
                        matrix reference_b = e(b)
                        matrix reference_V = e(V)
                        local reference_N = e(N)
                        local reference_N_full = e(N_full)
                        local reference_singletons = e(num_singletons)
                        quietly replace expected = e(sample)
                    }
                    else {
                        assert mreldif(e(b), reference_b) < 1e-12
                        assert mreldif(e(V), reference_V) < 1e-12
                        assert e(N) == `reference_N' & e(N_full) == `reference_N_full'
                        assert e(num_singletons) == `reference_singletons'
                        assert expected == e(sample)
                        * Explicit one-thread execution must retain the same rows.
                        quietly `cmd' `model' `qualifier', `spec' threads(1)
                        assert expected == e(sample)
                        assert mreldif(e(b), reference_b) < 1e-10
                        assert mreldif(e(V), reference_V) < 1e-10
                        local ++checked
                    }
                }
            }
        }
    }
}

* Missing output buffers must not participate in input screening. Automatic
* FE targets retain the plain-input path in OLS and PPML.
foreach v in out_r __hdfe1__ __hdfe2__ {
    generate double reference_`v' = .
}
foreach cmd in creghdfe civreghdfe cppmlhdfe {
    local model y x1 x2
    local resopt residuals(out_r)
    if "`cmd'" == "civreghdfe" local model y x2 (x1=z)
    if "`cmd'" == "cppmlhdfe" {
        local model yp x1 x2
        local resopt d(out_r)
    }
    foreach qualifier in "" "if t > 2 in 2/130000" {
        foreach variant in baseline candidate {
            quietly capture adopath - `"`baseline'"'
            quietly capture adopath - `"`candidate'"'
            quietly adopath ++ `"``variant''"'
            discard
            quietly `cmd' `model' `qualifier', absorb(g t, savefe) `resopt'
            if "`variant'" == "baseline" quietly replace expected = e(sample)
            else assert expected == e(sample)
            foreach v in out_r __hdfe1__ __hdfe2__ {
                if "`variant'" == "baseline" quietly replace reference_`v' = `v'
                else assert `v' == reference_`v' | reldif(`v', reference_`v') < 1e-12
                drop `v'
            }
            if "`variant'" == "candidate" local ++checked
        }
    }
}

* All-missing and no-observation errors, including a PPML negative outcome
* whose row must be excluded for a missing regressor before y validation.
quietly replace x1 = .a
quietly replace yp = -1
foreach cmd in creghdfe civreghdfe cppmlhdfe {
    local model y x1 x2
    if "`cmd'" == "civreghdfe" local model y x2 (x1=z)
    if "`cmd'" == "cppmlhdfe" local model yp x1 x2
    foreach variant in baseline candidate {
        quietly capture adopath - `"`baseline'"'
        quietly capture adopath - `"`candidate'"'
        quietly adopath ++ `"``variant''"'
        discard
        capture quietly `cmd' `model', absorb(g t)
        local rc = _rc
        if "`variant'" == "baseline" local expected_rc = `rc'
        else {
            assert `rc' == `expected_rc'
            assert `rc' == cond("`cmd'"=="civreghdfe",2001,2000)
            local ++checked
        }
    }
}
di "HDFE_NATIVE_SCREEN_PASS comparisons=`checked'"
