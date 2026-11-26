* Compare temporary/output initialization against a frozen ado directory.
* Both directories must contain the same compatible plugin binary.
args baseline candidate
if `"`baseline'"' == "" | `"`candidate'"' == "" error 198
clear all
set more off
set seed 20261128
set obs 1200
generate long g = ceil(_n/12)
generate byte t = mod(_n-1,12)+1
generate double x1 = rnormal()
generate double x2 = rnormal()
generate double z = x1 + rnormal()
generate double y = x1 + x2 + g/100 + t/12 + rnormal()
generate double yp = rpoisson(exp(.2*x1 + .1*x2 + g/100 + t/12))
generate double w = mod(_n,3)+1
replace x1 = .a in 2
replace g = .z in 4
replace z = .b in 6
replace g = 101 in 9
replace yp = 0 if g == 2
generate byte expected = 0
foreach v in out_r out_fe out_t out_g {
    generate double reference_`v' = .
}
local checked = 0

foreach cmd in creghdfe civreghdfe cppmlhdfe {
    local model y x1 x2
    local resopt residuals(out_r)
    if "`cmd'" == "civreghdfe" local model y x2 (x1 = z)
    if "`cmd'" == "cppmlhdfe" {
        local model yp x1 x2
        local resopt d(out_r)
    }
    foreach scenario in constant residuals effects {
        local opts
        local outputs
        if "`scenario'" == "residuals" {
            local opts absorb(g t) `resopt'
            local outputs out_r
        }
        if "`scenario'" == "effects" {
            local opts absorb(out_fe=g out_t=t) `resopt'
            local outputs out_r out_fe out_t
            if "`cmd'" == "creghdfe" {
                local opts `opts' groupvar(out_g)
                local outputs `outputs' out_g
            }
        }
        foreach qualifier in "" "if t > 2 in 2/1100" "[fw=w]" {
            foreach variant in baseline candidate {
                capture adopath - `"`baseline'"'
                capture adopath - `"`candidate'"'
                adopath ++ `"``variant''"'
                discard
                quietly `cmd' `model' `qualifier', `opts'
                if "`variant'" == "baseline" {
                    matrix reference_b = e(b)
                    matrix reference_V = e(V)
                    local reference_N = e(N)
                    if "`cmd'" == "civreghdfe" local reference_N_full = e(N_full)
                    quietly replace expected = e(sample)
                    foreach v of local outputs {
                        quietly replace reference_`v' = `v'
                        local type_`v' : type `v'
                        local format_`v' : format `v'
                    }
                }
                else {
                    assert mreldif(e(b), reference_b) < 1e-12
                    assert mreldif(e(V), reference_V) < 1e-12
                    assert e(N) == `reference_N'
                    if "`cmd'" == "civreghdfe" assert e(N_full) == `reference_N_full'
                    assert expected == e(sample)
                    foreach v of local outputs {
                        assert `v' == reference_`v' | reldif(`v', reference_`v') < 1e-12
                        assert "`: type `v''" == "`type_`v''"
                        assert "`: format `v''" == "`format_`v''"
                    }
                    local ++checked
                }
                foreach v of local outputs {
                    drop `v'
                }
                * st_addvar must respect tempvar lifetime.
                capture ds __*
                assert _rc == 111
            }
        }
    }
}

* by: without absorb() must still create a constant on the active group,
* including groups whose observations do not start at row one.
foreach variant in baseline candidate {
    capture adopath - `"`baseline'"'
    capture adopath - `"`candidate'"'
    adopath ++ `"``variant''"'
    discard
    quietly bysort t: civreghdfe y x2 (x1=z), robust
    if "`variant'" == "baseline" {
        matrix reference_b = e(b)
        matrix reference_V = e(V)
        quietly replace expected = e(sample)
    }
    else {
        assert mreldif(e(b), reference_b) < 1e-12
        assert mreldif(e(V), reference_V) < 1e-12
        assert expected == e(sample)
        local ++checked
    }
}

* Diagnostic-only PPML also publishes correctly initialized tag/certificate
* variables when rows were excluded before the plugin runs.
foreach variant in baseline candidate {
    capture adopath - `"`baseline'"'
    capture adopath - `"`candidate'"'
    adopath ++ `"``variant''"'
    discard
    quietly cppmlhdfe yp x1 x2 if t > 2, absorb(g t) tagsep(out_r) zvarname(out_fe)
    foreach v in out_r out_fe {
        if "`variant'" == "baseline" quietly replace reference_`v' = `v'
        else assert `v' == reference_`v' | reldif(`v', reference_`v') < 1e-12
        drop `v'
    }
    if "`variant'" == "candidate" local ++checked
}

* Clock first use, warm use and re-registration preserve incoming r().
discard
quietly summarize y
local n = r(N)
local mean = r(mean)
_ctools_clock start
assert r(N) == `n' & r(mean) == `mean'
_ctools_lap elapsed start `start'
assert `elapsed' >= 0
assert r(N) == `n' & r(mean) == `mean'
discard
quietly summarize y
_ctools_clock start
assert r(N) == `n' & r(mean) == `mean'

* Cleanup removes both namespaces and both object types, preserving others.
foreach name in __ado_buffer_a __ado_buffer_b _ado_buffer_c {
    scalar `name' = 1
    matrix `name' = J(2,2,1)
}
scalar unrelated_scalar = 7
matrix unrelated_matrix = J(2,2,7)
_ctools_cleanup ado_buffer
_ctools_cleanup ado_buffer
foreach name in __ado_buffer_a __ado_buffer_b _ado_buffer_c {
    capture confirm scalar `name'
    assert _rc == 111
    capture confirm matrix `name'
    assert _rc == 111
}
assert scalar(unrelated_scalar) == 7
assert el(unrelated_matrix,1,1) == 7
di "HDFE_ADO_BUFFERS_PASS comparisons=`checked'"
