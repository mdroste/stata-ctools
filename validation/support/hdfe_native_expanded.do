* Weighted/factor sample preparation against frozen wrappers, same plugin.
args baseline candidate
if `"`baseline'"' == "" | `"`candidate'"' == "" error 198
clear all
set more off
set seed 20261130
set obs 65539
gen long g = ceil(_n/13)
gen byte t = mod(_n-1,13)+1
gen double x = rnormal()
gen double x2 = rnormal()
gen double z = x+rnormal()
gen double endog = .5*z+rnormal()
gen double y = x+x2+.5*endog+g/5000+t/13+rnormal()
gen double yp = rpoisson(exp(.1*x+.1*x2+g/5000+t/13))
gen byte f = mod(_n,4)
gen double w = 1+mod(_n,3)
gen double iw = w+.5
gen str16 gs = string(g)
gen double off = t/20
gen double exposure = 1+mod(_n,3)
replace x = .a in 4095/4097
replace f = .z in 8192
replace g = .a in 16384
replace gs = "" in 16384
replace y = . in 32768
replace yp = . in 32768
replace f = -1 in 32768
replace exposure = -1 in 32768
replace z = .b in 50000
replace w = 0 in 11
replace w = .a in 12
replace f = 9 in 11/12
replace f = 8 in 4095/4097
replace iw = w if missing(w) | w==0
gen byte expected = 0
local checked = 0
foreach cmd in creghdfe civreghdfe cppmlhdfe {
    local weights none fw aw pw
    foreach shape in plain factor interaction continuous base alllevels selected {
        local rhs x x2
        if "`shape'" == "factor" local rhs x i.f
        if "`shape'" == "interaction" local rhs c.x##i.f
        if "`shape'" == "continuous" local rhs c.x#c.x x2
        if "`shape'" == "base" local rhs x ib2.f
        if "`shape'" == "alllevels" local rhs x ibn.f
        if "`shape'" == "selected" local rhs x i(1/3).f
        local model y `rhs'
        if "`cmd'" == "cppmlhdfe" local model yp `rhs'
        if "`cmd'" == "civreghdfe" local model y `rhs' (endog=z)
        foreach wt of local weights {
            local weight
            local covariance vce(robust)
            if "`wt'" != "none" local weight [`wt'=w]
            if "`wt'" == "iw" {
                local weight [iw=iw]
                local covariance
            }
            foreach variant in baseline candidate {
                quietly capture adopath - `"`baseline'"'
                quietly capture adopath - `"`candidate'"'
                quietly adopath ++ `"``variant''"'
                discard
                capture quietly `cmd' `model' if t>1 in 2/65000 `weight', absorb(g t) `covariance'
                local rc = _rc
                if `rc' {
                    di as error "EXPANDED_CASE command=`cmd' shape=`shape' weight=`wt' variant=`variant' rc=`rc'"
                    error `rc'
                }
                if "`variant'" == "baseline" {
                    matrix ref_b = e(b)
                    matrix ref_V = e(V)
                    local ref_N = e(N)
                    local ref_N_full = e(N_full)
                    local ref_df_a = e(df_a)
                    local ref_sing = e(num_singletons)
                    local ref_names : colfullnames e(b)
                    quietly replace expected = e(sample)
                }
                else {
                    assert mreldif(e(b),ref_b)<1e-12 & mreldif(e(V),ref_V)<1e-12
                    assert e(N)==`ref_N' & e(N_full)==`ref_N_full'
                    assert e(df_a)==`ref_df_a' & e(num_singletons)==`ref_sing'
                    assert `"`: colfullnames e(b)'"' == `"`ref_names'"'
                    assert e(sample)==expected
                    if "`shape'" == "interaction" {
                        quietly `cmd' `model' if t>1 in 2/65000 `weight', absorb(g t) `covariance' threads(1)
                        assert mreldif(e(b),ref_b)<1e-10 & mreldif(e(V),ref_V)<1e-10
                        assert e(sample)==expected
                    }
                    local ++checked
                }
            }
        }
    }
    * Expressions referring to prior estimates must see those estimates at
    * both marksample's evaluation and the later weight-expression staging.
    local model y x i.f
    if "`cmd'" == "cppmlhdfe" local model yp x i.f
    if "`cmd'" == "civreghdfe" local model y x i.f (endog=z)
    foreach variant in baseline candidate {
        quietly capture adopath - `"`baseline'"'
        quietly capture adopath - `"`candidate'"'
        quietly adopath ++ `"``variant''"'
        discard
        quietly regress y x x2 if t>1
        quietly `cmd' `model' if e(sample) [aw=w*(1+e(N)/100000)], absorb(g t)
        if "`variant'" == "baseline" {
            matrix ref_b=e(b)
            matrix ref_V=e(V)
            quietly replace expected=e(sample)
        }
        else {
            assert mreldif(e(b),ref_b)<1e-12 & mreldif(e(V),ref_V)<1e-12
            assert e(sample)==expected
            local ++checked
        }
    }
    * Invalid weights still fail before model/FE missingness is screened.
    foreach bad in -1 1.5 {
        quietly replace w=`bad' in 16384
        foreach variant in baseline candidate {
            quietly capture adopath - `"`baseline'"'
            quietly capture adopath - `"`candidate'"'
            quietly adopath ++ `"``variant''"'
            discard
            capture quietly `cmd' `model' [fw=w], absorb(g t)
            local rc=_rc
            if "`variant'" == "baseline" local ref_rc=`rc'
            else {
                assert `rc'==`ref_rc' & `rc'==cond(`bad'<0,402,401)
                local ++checked
            }
        }
    }
    quietly replace w=2 in 16384
    foreach variant in baseline candidate {
        quietly capture adopath - `"`baseline'"'
        quietly capture adopath - `"`candidate'"'
        quietly adopath ++ `"``variant''"'
        discard
        capture quietly `cmd' `model' [iw=iw], absorb(g t)
        assert _rc==101
    }
}

* Early screening also precedes complex-FE parsing, optional offset/exposure
* validation, and weighted models without absorbed effects.
foreach cmd in creghdfe civreghdfe cppmlhdfe {
    local model y x i.f
    if "`cmd'" == "civreghdfe" local model y x i.f (endog=z)
    if "`cmd'" == "cppmlhdfe" local model yp x i.f
    local specs `""absorb(gs t)" "absorb(g##c.x2 t)" """'
    if "`cmd'" == "cppmlhdfe" local specs `"`specs' "absorb(g t) offset(off)" "absorb(g t) exposure(exposure)""'
    foreach spec of local specs {
        foreach variant in baseline candidate {
            quietly capture adopath - `"`baseline'"'
            quietly capture adopath - `"`candidate'"'
            quietly adopath ++ `"``variant''"'
            discard
            quietly `cmd' `model' [aw=w+1], `spec'
            if "`variant'" == "baseline" {
                matrix ref_b=e(b)
                matrix ref_V=e(V)
                local ref_N=e(N)
                local ref_N_full=e(N_full)
                local ref_names : colfullnames e(b)
                quietly replace expected=e(sample)
            }
            else {
                assert mreldif(e(b),ref_b)<1e-12 & mreldif(e(V),ref_V)<1e-12
                assert e(sample)==expected & e(N)==`ref_N' & e(N_full)==`ref_N_full'
                assert `"`: colfullnames e(b)'"' == `"`ref_names'"'
                local ++checked
            }
        }
    }
}

* Planning does not expand levels, inspect excluded invalid categories, or
* disturb r() used by qualifiers/weight expressions. TS stays on its old path.
quietly summarize x
tempname mean
scalar `mean' = r(mean)
foreach term in "i.f##c.x" "f#c.x" "ib(first).f" "ib(last).f" "ib(freq).f" "i(1/3).f" {
    _ctools_hdfe_screen plan `term'
    assert `__screen_ok' == 1 & r(mean)==scalar(`mean')
}
foreach term in "L.x" "F2.x" "D.x" "S2.x" "cL2.x#i.f" "L(1/2).x" {
    _ctools_hdfe_screen plan `term'
    assert `__screen_ok' == 0 & r(mean)==scalar(`mean')
}
di "HDFE_NATIVE_EXPANDED_PASS comparisons=`checked'"
