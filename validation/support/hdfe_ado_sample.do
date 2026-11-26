* Compare sample preparation with a frozen pre-change ado directory.
* Both directories must contain the same compatible plugin binary.
* Usage: do validation/support/hdfe_ado_sample.do baseline candidate
args baseline candidate
if `"`baseline'"' == "" | `"`candidate'"' == "" error 198
clear all
set more off
set seed 20261125
set obs 1200
generate long g = ceil(_n/12)
generate byte t = mod(_n-1,12)+1
generate double x1 = rnormal()
generate double x2 = rnormal()
generate double z = x1 + rnormal()
generate double y = x1 + x2 + g/100 + t/12 + rnormal()
generate double yp = rpoisson(exp(.2*x1 + .1*x2 + g/100 + t/12))
generate double w = mod(_n,3)+1
replace w = 0 in 1
replace x1 = .a in 2
replace x2 = . in 3
replace g = .z in 4
replace t = . in 5
replace z = .b in 6
replace y = . in 7
replace yp = . in 8
replace g = 101 in 9
generate str5 gs = string(g)
replace gs = "" in 4
generate byte expected = 0
generate byte prior = mod(_n,4) != 0
local checked = 0
local successful = 0

foreach cmd in creghdfe civreghdfe cppmlhdfe {
    local model y x1 x2
    if "`cmd'" == "civreghdfe" local model y x2 (x1 = z)
    if "`cmd'" == "cppmlhdfe" local model yp x1 x2
    foreach spec in "absorb(g t)" "absorb(g t) vce(robust)" "absorb(g t) vce(cluster g)" "absorb(g t) keepsingletons" "absorb(g#t) keepsingletons" "absorb(gs t)" "absorb(g##c.x2 t)" {
        foreach qualifier in "" "if prior in 2/1100" "[fw=w]" {
            foreach variant in baseline candidate {
                capture adopath - `"`baseline'"'
                capture adopath - `"`candidate'"'
                adopath ++ `"``variant''"'
                discard
                capture quietly `cmd' `model' `qualifier', `spec'
                local rc = _rc
                if "`variant'" == "baseline" {
                    local reference_rc = `rc'
                    if "`spec'" == "absorb(g t)" assert !`rc'
                    if !`rc' {
                        matrix reference_b = e(b)
                        matrix reference_V = e(V)
                        local reference_N = e(N)
                        if "`cmd'" == "civreghdfe" local reference_N_full = e(N_full)
                        local reference_names : colfullnames e(b)
                        quietly replace expected = e(sample)
                    }
                }
                else {
                    assert `rc' == `reference_rc'
                    if !`rc' {
                        assert mreldif(e(b), reference_b) < 1e-12
                        assert mreldif(e(V), reference_V) < 1e-12
                        assert e(N) == `reference_N'
                        if "`cmd'" == "civreghdfe" assert e(N_full) == `reference_N_full'
                        assert `"`: colfullnames e(b)'"' == `"`reference_names'"'
                        assert expected == e(sample)
                        local ++successful
                    }
                    local ++checked
                }
            }
        }
    }
    * The fast path must not consume or clear previous estimates early.
    foreach variant in baseline candidate {
        capture adopath - `"`baseline'"'
        capture adopath - `"`candidate'"'
        adopath ++ `"``variant''"'
        discard
        quietly regress y x1 x2 if prior
        quietly `cmd' `model' if e(sample) [aw=w+e(N)/1000], absorb(g t)
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
    * Batching must preserve weight validation, even on missing-FE rows.
    foreach badweight in -1 1.5 {
        quietly replace w = `badweight' in 4
        foreach variant in baseline candidate {
            capture adopath - `"`baseline'"'
            capture adopath - `"`candidate'"'
            adopath ++ `"``variant''"'
            discard
            capture quietly `cmd' `model' [fw=w], absorb(g t)
            local rc = _rc
            if "`variant'" == "baseline" local reference_rc = `rc'
            else {
                assert `rc' == `reference_rc'
                local ++checked
            }
        }
    }
    quietly replace w = 2 in 4
}
* PPML must still reject a negative outcome before dropping its singleton,
* leave named outputs unpublished, and ignore a negative excluded outcome.
replace yp = -1 in 9
foreach variant in baseline candidate {
    capture adopath - `"`baseline'"'
    capture adopath - `"`candidate'"'
    adopath ++ `"``variant''"'
    discard
    capture quietly cppmlhdfe yp x1 x2, absorb(g t) d(fe_output)
    assert _rc == 198
    capture confirm variable fe_output
    assert _rc == 111
    quietly cppmlhdfe yp x1 x2 if _n != 9, absorb(g t)
    capture quietly cppmlhdfe yp x1 x2 if 0, absorb(g t)
    assert _rc == 2000
}
assert `successful' >= 40
di "HDFE_ADO_SAMPLE_PASS comparisons=`checked' successful=`successful'"
