* s2_hdfe: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "s2_hdfe"

if "`__case'" == "s2_hdfe" {
do validation/support/validate_setup.do
set seed 20261104

* creghdfe against reghdfe on the same model: e(b), e(V) (relative
* difference at most tol()), the counts e(N), e(df_r), e(df_a),
* e(df_a_nested), e(N_clust), e(N_clust#), e(dkraay_bw), e(N_hdfe) and
* e(N_hdfe_extended), the absorbed degrees of freedom e(dof_table)
* exactly, and the names in e(absvars), e(extended_absvars), e(clustvar#),
* e(vce), e(vcetype), e(title3) and e(title4)
capture program drop _s2h_cmp
program define _s2h_cmp
    syntax, REF(string) CT(string) NAME(string) [TOL(real 1e-7) FIT]
    capture noisily {
        ereturn clear
        quietly `ref'
        local fitstats r2 r2_a r2_within r2_a_within rss tss mss rmse ll ll_0 F
        if "`fit'" != "" {
            foreach s of local fitstats {
                local f_`s' = e(`s')
            }
            local f_report_constant = e(report_constant)
            local f_names : colnames e(b)
        }
        tempname b V dof
        matrix `b' = e(b)
        matrix `V' = e(V)
        matrix `dof' = e(dof_table)
        local scalars N df_r df_a df_a_nested df_m N_clust N_clustervars dkraay_bw N_hdfe N_hdfe_extended num_singletons
        foreach s of local scalars {
            local r_`s' = e(`s')
        }
        local ncl = cond(missing(e(N_clustervars)), 0, e(N_clustervars))
        forvalues d = 1/`ncl' {
            local r_N_clust`d' = e(N_clust`d')
            local r_clustvar`d' `e(clustvar`d')'
        }
        local macros absvars extended_absvars clustvar vce vcetype title3 title4
        foreach m of local macros {
            local r_`m' `"`e(`m')'"'
        }
        ereturn clear
        quietly `ct'
        local db = mreldif(`b', e(b))
        local dV = mreldif(`V', e(V))
        if !(`db' <= `tol') | !(`dV' <= `tol') {
            di as error "`name': mreldif e(b) `db', e(V) `dV'"
            exit 9
        }
        if "`fit'" != "" {
            foreach s of local fitstats {
                if !((missing(`f_`s'') & missing(e(`s'))) | reldif(`f_`s'', e(`s')) <= `tol') {
                    di as error "`name': e(`s') `f_`s'' vs " e(`s')
                    exit 9
                }
            }
            assert `f_report_constant' == e(report_constant)
            assert "`f_names'" == "`: colnames e(b)'"
        }
        tempname dof2
        matrix `dof2' = e(dof_table)
        if rowsof(`dof') != rowsof(`dof2') | mreldif(`dof', `dof2') != 0 {
            matrix list `dof'
            matrix list `dof2'
            exit 9
        }
        foreach s of local scalars {
            if !((missing(`r_`s'') & missing(e(`s'))) | `r_`s'' == e(`s')) {
                di as error "`name': e(`s') `r_`s'' vs " e(`s')
                exit 9
            }
        }
        forvalues d = 1/`ncl' {
            assert `r_N_clust`d'' == e(N_clust`d')
            assert "`r_clustvar`d''" == "`e(clustvar`d')'"
        }
        foreach m of local macros {
            if `"`r_`m''"' != `"`e(`m')'"' {
                di as error `"`name': e(`m') "`r_`m''" vs "`e(`m')'""'
                exit 9
            }
        }
    }
    if _rc test_fail "`name'" "r(`=_rc')"
    else test_pass "`name'"
end

* Saved variables of creghdfe and reghdfe agree: max reldif at most tol()
capture program drop _s2h_vars
program define _s2h_vars
    syntax, PAIRS(string) NAME(string) [TOL(real 1e-7)]
    capture noisily {
        while `"`pairs'"' != "" {
            gettoken a pairs : pairs
            gettoken b pairs : pairs
            confirm numeric variable `a' `b'
            tempvar d
            quietly generate double `d' = reldif(`a', `b')
            quietly count if missing(`a') != missing(`b')
            assert r(N) == 0
            quietly summarize `d'
            if r(N) == 0 | r(max) > `tol' {
                di as error "`name': `a' vs `b' max reldif " r(max)
                exit 9
            }
            drop `d'
        }
    }
    if _rc test_fail "`name'" "r(`=_rc')"
    else test_pass "`name'"
end

* =========================================================================
* HDFE-01: multi-way clustering
* =========================================================================
ctools_fixture nlswork, clear
keep if idcode < 600
generate byte w = 1 + mod(idcode, 3)
generate sgrade = string(grade)

_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode year) vce(cluster idcode year)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode year) vce(cluster idcode year)) ///
    name("HDFE-01 two-way clusters equal to the absorbed effects")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode year) vce(cluster idcode year ind_code)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode year) vce(cluster idcode year ind_code)) ///
    name("HDFE-01 three-way clusters")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode year) cluster(year idcode)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode year) cluster(year idcode)) ///
    name("HDFE-01 cluster() with two variables, reversed order")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode) vce(cluster race#year occ_code)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode) vce(cluster race#year occ_code)) ///
    name("HDFE-01 an interaction as a cluster variable")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode year) vce(cluster idcode#year year)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode year) vce(cluster idcode#year year)) ///
    name("HDFE-01 cluster on an interaction of absorbed effects")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode occ_code) vce(cluster birth_yr year)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode occ_code) vce(cluster birth_yr year)) ///
    name("HDFE-01 absorbed effect nested in one cluster variable")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode year) vce(cluster sgrade year)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode year) vce(cluster sgrade year)) ///
    name("HDFE-01 string cluster variable")
_s2h_cmp, ref(reghdfe ln_wage age tenure [fw=w], absorb(idcode year) vce(cluster idcode year ind_code)) ///
    ct(creghdfe ln_wage age tenure [fw=w], absorb(idcode year) vce(cluster idcode year ind_code)) ///
    name("HDFE-01 three-way clusters, fweights")
_s2h_cmp, ref(reghdfe ln_wage age tenure [pw=w], absorb(idcode year) vce(cluster idcode year)) ///
    ct(creghdfe ln_wage age tenure [pw=w], absorb(idcode year) vce(cluster idcode year)) ///
    name("HDFE-01 two-way clusters, pweights")
_s2h_cmp, ref(reghdfe ln_wage age tenure [aw=w], absorb(idcode occ_code) vce(cluster ind_code year)) ///
    ct(creghdfe ln_wage age tenure [aw=w], absorb(idcode occ_code) vce(cluster ind_code year)) ///
    name("HDFE-01 two-way clusters, aweights")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode year) vce(cluster idcode year) dofadjustments(none)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode year) vce(cluster idcode year) dofadjustments(none)) ///
    name("HDFE-01 two-way clusters, dofadjustments(none)")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode year occ_code) vce(cluster ind_code year)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode year occ_code) vce(cluster ind_code year)) ///
    name("HDFE-01 two-way clusters, three absorbed effects")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode##c.ttl_exp year) vce(cluster idcode year)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode##c.ttl_exp year) vce(cluster idcode year)) ///
    name("HDFE-01 two-way clusters with a heterogeneous slope")
* The covariance correction is announced as in reghdfe
capture noisily {
    quietly creghdfe ln_wage age tenure, absorb(idcode year) vce(cluster idcode year ind_code)
    assert e(N_clustervars) == 3 & "`e(clustvar)'" == "idcode year ind_code"
    assert "`e(title3)'" == "Statistics robust to heteroskedasticity"
    creghdfe
}
if _rc test_fail "HDFE-01 replay of a multi-way clustered fit" "r(`=_rc')"
else test_pass "HDFE-01 replay of a multi-way clustered fit"
capture creghdfe ln_wage age, absorb(idcode) vce(cluster idcode year) wildbootstrap(reps(9) seed(1))
if _rc == 198 test_pass "HDFE-01 wildbootstrap() refuses multi-way clusters"
else test_fail "HDFE-01 wildbootstrap() refuses multi-way clusters" "r(`=_rc')"

* HDFE-27: with an absorbed effect nested in one cluster dimension and 2
* clusters in the other, reghdfe's PSD correction (the full matrix clipped,
* then its coefficient block replaced by that block clipped alone) is
* indefinite, and Stata's ereturn post zeroes it or keeps it depending on
* rounding: reghdfe posts zero with vce(cluster foreign turn). creghdfe posts
* a PSD e(V) whose coefficient block is reghdfe's, the e(V) reghdfe posts with
* noconstant, in either order.
preserve
sysuse auto, clear
tempname h27_Vnc h27_b h27_V1
quietly reghdfe price weight length, absorb(turn) vce(cluster turn foreign) noconstant
matrix `h27_Vnc' = e(V)
quietly reghdfe price weight length, absorb(turn) vce(cluster turn foreign)
matrix `h27_b' = e(b)
foreach order in "turn foreign" "foreign turn" {
    capture noisily {
        quietly creghdfe price weight length, absorb(turn) vce(cluster `order')
        assert mreldif(e(b), `h27_b') < 1e-10
        mata: V = st_matrix("e(V)"); assert(!allof(V, 0)); s = sqrt(diagonal(V)); ///
            assert(min(symeigenvalues(V :/ (s * s'))) > -1e-10); ///
            assert(mreldif(V[1..2, 1..2], st_matrix("`h27_Vnc'")) < 1e-10)
        if "`order'" == "turn foreign" matrix `h27_V1' = e(V)
        else assert mreldif(e(V), `h27_V1') < 1e-12
    }
    if _rc test_fail "HDFE-27 two-way clusters, indefinite PSD correction (`order')" "r(`=_rc')"
    else test_pass "HDFE-27 two-way clusters, indefinite PSD correction (`order')"
}
* An indefinite correction that Stata keeps is posted as reghdfe posts it
_s2h_cmp, ref(reghdfe price weight length mpg, absorb(turn) vce(cluster turn foreign)) ///
    ct(creghdfe price weight length mpg, absorb(turn) vce(cluster turn foreign)) ///
    name("HDFE-27 indefinite correction kept as reghdfe's")
generate byte hi = mpg > 20
_s2h_cmp, ref(reghdfe price weight length, absorb(turn) vce(cluster turn foreign hi)) ///
    ct(creghdfe price weight length, absorb(turn) vce(cluster turn foreign hi)) ///
    name("HDFE-27 three-way clusters with 2-cluster dimensions")
restore

* =========================================================================
* HDFE-04: vce(dkraay [#])
* =========================================================================
tsset idcode year
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode) vce(dkraay 3)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode) vce(dkraay 3)) ///
    name("HDFE-04 vce(dkraay 3)")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode year) vce(dkraay)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode year) vce(dkraay)) ///
    name("HDFE-04 default bandwidth, time effect nested in the periods")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode) vce(dkraay 40)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode) vce(dkraay 40)) ///
    name("HDFE-04 bandwidth beyond the time periods")
_s2h_cmp, ref(reghdfe ln_wage age tenure [aw=w], absorb(idcode) vce(dkraay 2)) ///
    ct(creghdfe ln_wage age tenure [aw=w], absorb(idcode) vce(dkraay 2)) ///
    name("HDFE-04 aweights")
_s2h_cmp, ref(reghdfe ln_wage age tenure [fw=w], absorb(idcode) vce(dkraay 4)) ///
    ct(creghdfe ln_wage age tenure [fw=w], absorb(idcode) vce(dkraay 4)) ///
    name("HDFE-04 fweights")
_s2h_cmp, ref(reghdfe ln_wage age tenure [pw=w], absorb(idcode##c.ttl_exp) vce(dkraay 1)) ///
    ct(creghdfe ln_wage age tenure [pw=w], absorb(idcode##c.ttl_exp) vce(dkraay 1)) ///
    name("HDFE-04 pweights and a slope, bandwidth 1")
_s2h_cmp, ref(reghdfe ln_wage age tenure, noabsorb vce(dkraay 2)) ///
    ct(creghdfe ln_wage age tenure, vce(dkraay 2)) ///
    name("HDFE-04 without absorbed effects")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode year) vce(dkr 2) dofadjustments(none)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode year) vce(dkr 2) dofadjustments(none)) ///
    name("HDFE-04 abbreviation dkr, dofadjustments(none)")
* Header: reghdfe's titles, the period count and the bandwidth
capture noisily {
    quietly creghdfe ln_wage age tenure, absorb(idcode) vce(dkraay 3)
    assert "`e(vce)'" == "dkraay" & "`e(vcetype)'" == "Robust" & e(dkraay_bw) == 3
    assert "`e(title4)'" == "and autocorrelation (Driscoll-Kraay)"
    tempfile hlog
    quietly log using "`hlog'", text replace name(s2h)
    creghdfe
    quietly log close s2h
    local found_dk = 0
    local found_cl = 0
    tempname fh
    file open `fh' using "`hlog'", read text
    file read `fh' line
    while r(eof) == 0 {
        if strpos(`"`line'"', "DK Bandwidth") & strpos(`"`line'"', "= ") local found_dk = 1
        if strpos(`"`line'"', "Number of clusters (year)") local found_cl = 1
        file read `fh' line
    }
    file close `fh'
    assert `found_dk' & `found_cl'
}
if _rc test_fail "HDFE-04 header shows the periods and the DK bandwidth" "r(`=_rc')"
else test_pass "HDFE-04 header shows the periods and the DK bandwidth"
* reghdfe's errors: r(9) for a bad bandwidth, options, or no time variable
local dk_ok = 1
foreach spec in "dkraay 0" "dkraay x" "dkraay 2.5" "dkraay, small" {
    capture reghdfe ln_wage age, absorb(idcode) vce(`spec')
    local rc_ref = _rc
    capture creghdfe ln_wage age, absorb(idcode) vce(`spec')
    if _rc != `rc_ref' | `rc_ref' == 0 local dk_ok = 0
}
tsset, clear
capture reghdfe ln_wage age, absorb(idcode) vce(dkraay)
local rc_ref = _rc
capture creghdfe ln_wage age, absorb(idcode) vce(dkraay)
if _rc != `rc_ref' | `rc_ref' == 0 local dk_ok = 0
if `dk_ok' test_pass "HDFE-04 vce(dkraay) errors as reghdfe"
else test_fail "HDFE-04 vce(dkraay) errors as reghdfe" "return codes differ"

* =========================================================================
* HDFE-02: heterogeneous slopes
* =========================================================================
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode##c.ttl_exp year)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode##c.ttl_exp year)) ///
    name("HDFE-02 id##c.x with a second effect")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode##c.(ttl_exp hours) year#c.age)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode##c.(ttl_exp hours) year#c.age)) ///
    name("HDFE-02 id##c.(x z) and a slope-only effect")
* reghdfe fails on this sample with id##(c.x c.z) (r(3498), missing values
* after partialling out); it is compared on the auto data below
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp year)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp year)) ///
    name("HDFE-02 slope without intercept, id#c.x")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp idcode)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp idcode)) ///
    name("HDFE-02 slope-only effect before its intercept")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(i.idcode#i.race##c.ttl_exp year)) ///
    ct(creghdfe ln_wage age tenure, absorb(i.idcode#i.race##c.ttl_exp year)) ///
    name("HDFE-02 interaction with a slope")
_s2h_cmp, ref(reghdfe ln_wage age tenure [aw=w], absorb(idcode#c.ttl_exp year)) ///
    ct(creghdfe ln_wage age tenure [aw=w], absorb(idcode#c.ttl_exp year)) ///
    name("HDFE-02 slopes, aweights")
_s2h_cmp, ref(reghdfe ln_wage age tenure [fw=w], absorb(idcode##c.ttl_exp year) vce(robust)) ///
    ct(creghdfe ln_wage age tenure [fw=w], absorb(idcode##c.ttl_exp year) vce(robust)) ///
    name("HDFE-02 slopes, fweights, robust")
_s2h_cmp, ref(reghdfe ln_wage age tenure [pw=w], absorb(idcode##c.(ttl_exp hours) year)) ///
    ct(creghdfe ln_wage age tenure [pw=w], absorb(idcode##c.(ttl_exp hours) year)) ///
    name("HDFE-02 several slopes, pweights")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode##c.ttl_exp year) vce(cluster idcode)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode##c.ttl_exp year) vce(cluster idcode)) ///
    name("HDFE-02 slopes, intercept nested in the cluster")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp year) vce(cluster idcode)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp year) vce(cluster idcode)) ///
    name("HDFE-02 a slope is never nested in the cluster")
foreach dof in none pairwise firstpair "clusters continuous" {
    _s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode##c.ttl_exp year occ_code) dofadjustments(`dof')) ///
        ct(creghdfe ln_wage age tenure, absorb(idcode##c.ttl_exp year occ_code) dofadjustments(`dof')) ///
        name("HDFE-02 slopes, dofadjustments(`dof')")
}
* Levels where a slope is constant (with an intercept) or zero (without):
* redundant, and dropped from the projection as reghdfe's invsym() drops them
preserve
clear
set obs 300
generate id = ceil(_n / 6)
generate grp = mod(_n, 5)
generate x = 0.1 * runiform() + id / 7
replace x = 0.1 if id == 3
replace x = 0.3 if id == 5
generate x0 = x
replace x0 = 0 if id == 4
generate z = rnormal()
generate y = z + id / 10 * x + rnormal()
_s2h_cmp, ref(reghdfe y z, absorb(id##c.x grp)) ct(creghdfe y z, absorb(id##c.x grp)) ///
    name("HDFE-02 slopes constant within levels")
_s2h_cmp, ref(reghdfe y z, absorb(id#c.x0 grp)) ct(creghdfe y z, absorb(id#c.x0 grp)) ///
    name("HDFE-02 slopes zero within a level")
restore

preserve
sysuse auto, clear
_s2h_cmp, ref(reghdfe price length, absorb(rep78##(c.mpg c.weight) foreign)) ///
    ct(creghdfe price length, absorb(rep78##(c.mpg c.weight) foreign)) ///
    name("HDFE-02 id##(c.x c.z)")
restore

* savefe and predict, d with slopes
capture noisily {
    quietly reghdfe ln_wage age tenure, absorb(idcode##c.(ttl_exp hours) year, savefe) resid(r_ref) tol(1e-12)
    quietly predict double d_ref, d
    rename (__hdfe1__ __hdfe1__Slope1 __hdfe1__Slope2 __hdfe2__) (a_ref s1_ref s2_ref y_ref)
    quietly creghdfe ln_wage age tenure, absorb(idcode##c.(ttl_exp hours) year, savefe) resid(r_ct) tol(1e-12)
    quietly predict double d_ct, d
    confirm numeric variable __hdfe1__ __hdfe1__Slope1 __hdfe1__Slope2 __hdfe2__
}
if _rc test_fail "HDFE-02 savefe with slopes" "r(`=_rc')"
else _s2h_vars, pairs(a_ref __hdfe1__ s1_ref __hdfe1__Slope1 s2_ref __hdfe1__Slope2 y_ref __hdfe2__ d_ref d_ct r_ref r_ct) ///
    name("HDFE-02 savefe with slopes and predict, d")
capture drop a_ref s1_ref s2_ref y_ref d_ref d_ct r_ref r_ct __hdfe*
capture noisily {
    quietly reghdfe ln_wage age tenure [aw=w], absorb(idcode#c.ttl_exp year, savefe) tol(1e-12)
    rename (__hdfe1__Slope1 __hdfe2__) (s_ref y_ref)
    quietly creghdfe ln_wage age tenure [aw=w], absorb(idcode#c.ttl_exp year, savefe) tol(1e-12)
}
if _rc test_fail "HDFE-02 savefe with a slope-only effect, aweights" "r(`=_rc')"
else _s2h_vars, pairs(s_ref __hdfe1__Slope1 y_ref __hdfe2__) name("HDFE-02 savefe with a slope-only effect, aweights")
capture drop s_ref y_ref __hdfe*

* groupvar(): the mobility groups of the first two intercepts, skipping slopes
capture noisily {
    quietly reghdfe ln_wage age tenure, absorb(idcode##c.ttl_exp year) groupvar(g_ref)
    quietly creghdfe ln_wage age tenure, absorb(idcode##c.ttl_exp year) groupvar(g_ct)
    quietly count if g_ref != g_ct
    assert r(N) == 0
    assert "`: variable label g_ref'" == "`: variable label g_ct'"
    drop g_ref g_ct
}
if _rc test_fail "HDFE-02 groupvar() with slopes" "r(`=_rc')"
else test_pass "HDFE-02 groupvar() with slopes"

* Slopes alone: nothing absorbs the constant, and reghdfe fits without one
* (no _cons, uncentered e(tss), e(report_constant) 0) under every VCE
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp)) fit ///
    name("HDFE-02 slopes only, no constant, unadjusted")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp) vce(robust)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp) vce(robust)) fit ///
    name("HDFE-02 slopes only, no constant, robust")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp) vce(cluster idcode)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp) vce(cluster idcode)) fit ///
    name("HDFE-02 slopes only, no constant, one-way cluster")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp) vce(cluster idcode year)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp) vce(cluster idcode year)) fit ///
    name("HDFE-02 slopes only, no constant, two-way cluster")
_s2h_cmp, ref(reghdfe ln_wage age tenure [aw=w], absorb(idcode#c.(ttl_exp hours) year#c.grade)) ///
    ct(creghdfe ln_wage age tenure [aw=w], absorb(idcode#c.(ttl_exp hours) year#c.grade)) fit ///
    name("HDFE-02 slopes only, no constant, several slopes, aweights")
_s2h_cmp, ref(reghdfe ln_wage age tenure [fw=w], absorb(idcode#c.ttl_exp) vce(robust)) ///
    ct(creghdfe ln_wage age tenure [fw=w], absorb(idcode#c.ttl_exp) vce(robust)) fit ///
    name("HDFE-02 slopes only, no constant, fweights")
_s2h_cmp, ref(reghdfe ln_wage age tenure [pw=w], absorb(idcode#c.ttl_exp)) ///
    ct(creghdfe ln_wage age tenure [pw=w], absorb(idcode#c.ttl_exp)) fit ///
    name("HDFE-02 slopes only, no constant, pweights")
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp) dofadjustments(none)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp) dofadjustments(none)) fit ///
    name("HDFE-02 slopes only, no constant, dofadjustments(none)")
tsset idcode year
_s2h_cmp, ref(reghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp) vce(dkraay 2)) ///
    ct(creghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp) vce(dkraay 2)) fit ///
    name("HDFE-02 slopes only, no constant, dkraay")
tsset, clear
capture noisily {
    quietly reghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp, savefe) resid(r_ref) tol(1e-12)
    quietly predict double d_ref, d
    rename __hdfe1__Slope1 s_ref
    quietly creghdfe ln_wage age tenure, absorb(idcode#c.ttl_exp, savefe) resid(r_ct) tol(1e-12)
    quietly predict double d_ct, d
    * The footnote lists the slope with its inexact marker
    tempfile flog
    quietly log using "`flog'", text replace name(s2f)
    creghdfe
    quietly log close s2f
    local found = 0
    tempname fh
    file open `fh' using "`flog'", read text
    file read `fh' line
    while r(eof) == 0 {
        if strpos(`"`line'"', "idcode#c.ttl_exp |") & strpos(`"`line'"', "?") local found = 1
        file read `fh' line
    }
    file close `fh'
    assert `found'
}
if _rc test_fail "HDFE-02 slopes only, no constant: savefe and predict, d" "r(`=_rc')"
else _s2h_vars, pairs(s_ref __hdfe1__Slope1 d_ref d_ct r_ref r_ct) name("HDFE-02 slopes only, no constant: savefe and predict, d")
capture drop s_ref d_ref d_ct r_ref r_ct __hdfe*

* Refused as documented
capture creghdfe ln_wage, absorb(idcode#c.ttl_exp)
local ok = (_rc == 198)
capture confirm variable __hdfe1__
local ok = `ok' & (_rc == 111)
if `ok' test_pass "HDFE-02 slopes only without regressors exit r(198)"
else test_fail "HDFE-02 slopes only without regressors exit r(198)" "unexpected return code"

* =========================================================================
* HDFE-03: newvar= targets
* =========================================================================
capture noisily {
    quietly reghdfe ln_wage age tenure, absorb(FR=idcode YR=year) tol(1e-12)
    quietly creghdfe ln_wage age tenure, absorb(FC=idcode YC=year) tol(1e-12)
}
if _rc test_fail "HDFE-03 newvar= targets" "r(`=_rc')"
else _s2h_vars, pairs(FR FC YR YC) name("HDFE-03 newvar= targets")
capture noisily {
    quietly reghdfe ln_wage age tenure, absorb(SR = idcode##c.ttl_exp year) tol(1e-12)
    quietly creghdfe ln_wage age tenure, absorb(SC = idcode##c.ttl_exp year) tol(1e-12)
    confirm numeric variable SC SCSlope1
    capture confirm variable __hdfe2__
    assert _rc == 111
}
if _rc test_fail "HDFE-03 target with slopes: name and nameSlope1" "r(`=_rc')"
else _s2h_vars, pairs(SR SC SRSlope1 SCSlope1) name("HDFE-03 target with slopes: name and nameSlope1")
capture noisily {
    quietly reghdfe ln_wage age tenure, absorb(TR=idcode year, savefe) tol(1e-12)
    rename __hdfe2__ y_ref
    quietly creghdfe ln_wage age tenure, absorb(TC=idcode year, savefe) tol(1e-12)
    capture confirm variable __hdfe1__
    assert _rc == 111
}
if _rc test_fail "HDFE-03 a target with savefe" "r(`=_rc')"
else _s2h_vars, pairs(TR TC y_ref __hdfe2__) name("HDFE-03 a target with savefe")
capture drop y_ref __hdfe*
capture reghdfe ln_wage age, absorb(age=idcode)
local rc_ref = _rc
capture creghdfe ln_wage age, absorb(age=idcode)
local ok = (_rc == `rc_ref') & `rc_ref' == 110
capture creghdfe ln_wage age, absorb(FF=idcode) fevce(bootstrap, reps(5) seed(1))
local ok = `ok' & (_rc == 198)
capture confirm variable FF
local ok = `ok' & (_rc == 111)
if `ok' test_pass "HDFE-03 existing target r(110); fevce() with targets r(198)"
else test_fail "HDFE-03 existing target r(110); fevce() with targets r(198)" "unexpected return code"

* civreghdfe absorbs slopes since IV-02 (compared with ivreghdfe in
* s2_iv) and saves newvar= targets since the bootstrap session's Phase 4a
* (IV-03), checked in its own cases
capture civreghdfe ln_wage (tenure = union), absorb(idcode##c.ttl_exp)
local ok = (_rc == 0)
capture drop F
capture civreghdfe ln_wage (tenure = union), absorb(F=idcode)
local ok = `ok' & (_rc == 0)
capture confirm numeric variable F
local ok = `ok' & (_rc == 0)
capture drop F
if `ok' test_pass "HDFE-02 civreghdfe absorbs slopes, saves targets"
else test_fail "HDFE-02 civreghdfe absorbs slopes, saves targets" "unexpected return code"

* =========================================================================
* PPML-04: cppmlhdfe absorb(id##(c.x c.z))
* =========================================================================
sysuse auto, clear
capture noisily {
    quietly ppmlhdfe price length, absorb(rep78##(c.mpg))
    tempname b V
    matrix `b' = e(b)
    matrix `V' = e(V)
    local absvars `e(absvars)'
    quietly cppmlhdfe price length, absorb(rep78##(c.mpg))
    assert mreldif(`b', e(b)) < 1e-7 & mreldif(`V', e(V)) < 1e-7
    assert "`absvars'" == "`e(absvars)'"
}
if _rc test_fail "PPML-04 cppmlhdfe absorb(id##(c.x)) vs ppmlhdfe" "r(`=_rc')"
else test_pass "PPML-04 cppmlhdfe absorb(id##(c.x)) vs ppmlhdfe"
* ppmlhdfe 2.3 fails (r(3201)) with two slopes in one effect; the
* alternative syntax must give what id##c.(x z) gives
capture noisily {
    quietly cppmlhdfe price length, absorb(rep78##c.(mpg weight) foreign)
    tempname b V
    matrix `b' = e(b)
    matrix `V' = e(V)
    local dfa = e(df_a)
    quietly cppmlhdfe price length, absorb(rep78##(c.mpg c.weight) foreign)
    assert mreldif(`b', e(b)) == 0 & mreldif(`V', e(V)) == 0 & e(df_a) == `dfa'
    assert "`e(absvars)'" == "rep78##(c.mpg c.weight) foreign"
    quietly cppmlhdfe price length, absorb(rep78#c.(mpg weight) foreign)
    matrix `b' = e(b)
    quietly cppmlhdfe price length, absorb(PF=rep78#(c.mpg c.weight) foreign)
    assert mreldif(`b', e(b)) == 0
    confirm numeric variable PFSlope1 PFSlope2
}
if _rc test_fail "PPML-04 absorb(id##(c.x c.z)) equals id##c.(x z)" "r(`=_rc')"
else test_pass "PPML-04 absorb(id##(c.x c.z)) equals id##c.(x z)"
capture cppmlhdfe price length, absorb(rep78##(mpg c.weight))
if _rc == 198 test_pass "PPML-04 a non-continuous term in id##( ) is r(198)"
else test_fail "PPML-04 a non-continuous term in id##( ) is r(198)" "r(`=_rc')"

/* HDFE-05: group(), individual() and aggregation(), against reghdfe 6.13.1:
   one regression row per group (patent), and an absorbed individual
   (inventor) effect aggregated within the groups. reghdfe's singleton
   pass with individual effects needs moremata, so the comparisons use
   keepsingletons; the singleton rule is checked against a loop in Stata
   that applies it. With aweights or pweights reghdfe 6.13.1 does not
   minimize the weighted sum of squares, so the reference there is the
   regression on explicit individual-share variables. Robust and cluster
   F statistics are compared at tolerance(1e-12): reghdfe's at its default
   tolerance is off by about 2e-7 (creghdfe solves 100 times tighter). */
clear
set seed 4242
set obs 600
gen patent = _n
gen year = 1 + mod(_n, 12)
gen firm = 1 + mod(_n * 7, 30)
gen k = 1 + floor(runiform() * 4)
gen x = rnormal()
gen x2 = rnormal()
gen w = .5 + runiform()
expand k
gen inventor = 1 + floor(runiform() * 300)
bysort patent inventor: keep if _n == 1
bysort patent: egen double mab = mean(inventor / 100)
bysort patent: egen double sab = total(inventor / 100)
gen double eps = rnormal()
bysort patent: replace eps = eps[1]
gen double y = 1 + 2*x - x2 + year/10 + mab + eps
gen double ys = 1 + 2*x - x2 + year/10 + sab + eps
local T group(patent) individual(inventor)
local K keepsingletons
foreach c in "y|absorb(year inventor)" "ys|absorb(year inventor) aggregation(sum)" ///
    "y|absorb(inventor)" "ys|absorb(inventor) aggregation(sum)" "y|absorb(inventor year)" ///
    "y|absorb(year inventor firm)" "y|absorb(year inventor) vce(cluster firm) tolerance(1e-12)" ///
    "y|absorb(year inventor) vce(cluster year) tolerance(1e-12)" "y|absorb(year inventor) vce(robust) tolerance(1e-12)" ///
    "y|absorb(year inventor) noconstant" "y|absorb(year inventor) tolerance(1e-10)" {
    gettoken dv o : c, parse("|")
    local o : subinstr local o "|" ""
    _s2h_cmp, ref(reghdfe `dv' x x2, `o' `T' `K') ct(creghdfe `dv' x x2, `o' `T' `K') ///
        name("HDFE-05 `dv' `o'") fit
}
_s2h_cmp, ref(reghdfe y x x2 if year != 3, absorb(year inventor) `T' `K') ///
    ct(creghdfe y x x2 if year != 3, absorb(year inventor) `T' `K') name("HDFE-05 if") fit
_s2h_cmp, ref(reghdfe y x i.firm, absorb(year inventor) `T' `K') ///
    ct(creghdfe y x i.firm, absorb(year inventor) `T' `K') name("HDFE-05 factor regressor") fit
foreach o in "absorb(year)" "absorb(year firm) vce(cluster firm)" {
    _s2h_cmp, ref(reghdfe y x x2, `o' group(patent)) ct(creghdfe y x x2, `o' group(patent)) ///
        name("HDFE-05 group() alone `o'") fit
}
* residuals: reghdfe's are accurate to its tolerance
capture noisily {
    quietly reghdfe y x x2, absorb(year inventor) `T' `K' resid(r_ref)
    quietly creghdfe y x x2, absorb(year inventor) `T' `K' resid(r_ct)
    assert e(N) == 600
    count if e(sample)
    assert r(N) == 600
    tempvar dr
    gen double `dr' = abs(r_ct - r_ref)
    summarize `dr', meanonly
    assert r(max) < 1e-4 & r(N) == 600
    drop r_ref r_ct
}
if _rc test_fail "HDFE-05 residuals and e(sample)" "r(`=_rc')"
else test_pass "HDFE-05 residuals and e(sample)"
* aweights and pweights: the weighted fit on explicit individual shares
capture noisily {
    bysort patent: replace w = w[1]
    tempname bc bp
    quietly creghdfe y x x2 [aw=w], absorb(year inventor) `T' `K'
    matrix `bc' = e(b)
    local rss = e(rss)
    quietly creghdfe y x x2 [pw=w], absorb(year inventor) `T' `K'
    matrix `bp' = e(b)
    assert mreldif(`bp', `bc') < 1e-12
    preserve
    egen inum = group(inventor)
    summarize inum, meanonly
    local NI = r(max)
    bysort patent: gen n = _N
    forvalues i = 1/`NI' {
        quietly gen double s`i' = (inum == `i') / n
    }
    collapse (first) y x x2 year w (sum) s1-s`NI', by(patent)
    quietly regress y x x2 i.year s1-s`NI' [aw=w]
    assert reldif(_b[x], `bc'[1, 1]) < 1e-8 & reldif(_b[x2], `bc'[1, 2]) < 1e-8
    assert reldif(e(rss), `rss') < 1e-8
    restore
}
if _rc test_fail "HDFE-05 aweights and pweights are the weighted fit" "r(`=_rc')"
else test_pass "HDFE-05 aweights and pweights are the weighted fit"
* Singletons: a group whose individual is in no other group, and a level of
* another effect with one group, dropped to the fixed point
capture noisily {
    tempvar keep tag pi cnt dropg
    gen byte `keep' = 1
    local more = 1
    while `more' {
        local more = 0
        * an inventor in one kept group: that group goes
        bysort inventor patent: gen byte `pi' = (_n == 1) & `keep'
        bysort inventor: egen `cnt' = total(`pi')
        bysort patent: egen byte `dropg' = max(`cnt' == 1 & `keep')
        quietly count if `dropg' & `keep'
        if r(N) {
            quietly replace `keep' = 0 if `dropg'
            local more = 1
        }
        drop `pi' `cnt' `dropg'
        * a year with one kept group: that group goes
        bysort patent: gen byte `tag' = (_n == 1) & `keep'
        bysort year: egen `cnt' = total(`tag')
        bysort patent: egen byte `dropg' = max(`cnt' == 1 & `keep')
        quietly count if `dropg' & `keep'
        if r(N) {
            quietly replace `keep' = 0 if `dropg'
            local more = 1
        }
        drop `tag' `cnt' `dropg'
    }
    bysort patent: gen byte `tag' = _n == 1
    quietly count if `tag' & !`keep'
    local ns = r(N)
    assert `ns' > 0
    tempname b V
    quietly creghdfe y x x2 if `keep', absorb(year inventor) `T' `K'
    matrix `b' = e(b)
    matrix `V' = e(V)
    local N = e(N)
    local dfa = e(df_a)
    quietly creghdfe y x x2, absorb(year inventor) `T'
    assert e(num_singletons) == `ns' & e(N) == `N' & e(df_a) == `dfa'
    assert mreldif(e(b), `b') < 1e-10 & mreldif(e(V), `V') < 1e-10
    drop `tag'
}
if _rc test_fail "HDFE-05 singletons with individual effects" "r(`=_rc')"
else test_pass "HDFE-05 singletons with individual effects"
* Errors as reghdfe gives them, and what creghdfe refuses
local fails
foreach c in "198|absorb(year inventor) individual(inventor)" "100|absorb(year) `T'" ///
    "9|absorb(year inventor) `T' aggregation(median)" "198|absorb(year fe=inventor) `T'" ///
    "198|absorb(year inventor) `T' vce(bootstrap, reps(5) seed(1))" ///
    "198|absorb(year##c.x2 inventor) `T'" "111|absorb(year) group(nosuchvar)" {
    gettoken r o : c, parse("|")
    local o : subinstr local o "|" ""
    capture creghdfe y x x2, `o' `K'
    if _rc != `r' local fails `"`fails' [`o': `=_rc']"'
}
capture creghdfe y x x2 [fw=k], absorb(year inventor) `T'
if _rc != 198 local fails `"`fails' [fweight: `=_rc']"'
gen double xv = x + rnormal()
capture creghdfe y xv x2, absorb(year inventor) `T' `K'
if _rc != 498 local fails `"`fails' [not constant: `=_rc']"'
expand 2 in 1
capture creghdfe y x x2, absorb(year inventor) `T' `K'
if _rc != 459 local fails `"`fails' [duplicate pair: `=_rc']"'
if `"`fails'"' == "" test_pass "HDFE-05 errors and refusals"
else test_fail "HDFE-05 errors and refusals" `"`fails'"'

* =========================================================================
* S3/S4 items (HDFE-06 to HDFE-26): reghdfe 6.13.1's options, abbreviations,
* errors, output and stored results
* =========================================================================
set seed 20261109
* The display of a command without the lines that name it or count its
* iterations: caller local s2h_lines, the lines joined by "|"
capture program drop _s2h_log
program define _s2h_log
    gettoken name 0 : 0
    gettoken colon 0 : 0, parse(":")
    tempfile f
    tempname h
    quietly log using `"`f'"', text replace name(_s2hlog)
    noisily `0'
    quietly log close _s2hlog
    file open `h' using `"`f'"', read text
    file read `h' line
    local out
    while r(eof) == 0 {
        local keep = !inlist(substr(`"`line'"', 1, 2), ". ", "> ") & ///
            !strpos(`"`line'"', "converged in") & !strpos(`"`line'"', "log:") & ///
            !strpos(`"`line'"', "log type:") & !strpos(`"`line'"', "opened on:") & ///
            !strpos(`"`line'"', "closed on:") & !strpos(`"`line'"', "name:")
        if `keep' local out `"`out'|`line'"'
        file read `h' line
    }
    file close `h'
    c_local `name' `"`out'"'
end

* HDFE-09 and HDFE-08: residuals, res, resi, resid without a name replace
* _reghdfe_resid; residuals(newvar), res(), resid(), resid2() name it
sysuse auto, clear
capture noisily {
    quietly reghdfe price weight, absorb(turn) resid
    rename _reghdfe_resid ref_r
    foreach o in res resi resid residuals {
        quietly creghdfe price weight, absorb(turn) `o'
        assert "`e(resid)'" == "_reghdfe_resid"
        quietly creghdfe price weight, absorb(turn) `o'
        assert "`e(resid)'" == "_reghdfe_resid"
    }
    assert reldif(_reghdfe_resid, ref_r) < 1e-9 | (missing(_reghdfe_resid) & missing(ref_r))
    foreach o in "res(r_a)" "resid(r_b)" "residuals(r_c)" "resid2(r_d)" {
        quietly creghdfe price weight, absorb(turn) `o'
        local v = substr("`o'", strpos("`o'", "(") + 1, 3)
        assert "`e(resid)'" == "`v'"
        assert reldif(`v', ref_r) < 1e-9 | (missing(`v') & missing(ref_r))
    }
    capture creghdfe price weight, absorb(turn) resid residuals(zz)
    assert _rc == 198
    capture reghdfe price weight, absorb(turn) resid residuals(zz)
    assert _rc == 198
    capture creghdfe price weight, absorb(turn) residuals(r_a)
    assert _rc == 110
    * bootstrap: reruns its command with the residuals each time
    quietly bootstrap, reps(3) seed(1) nodots: creghdfe price weight, absorb(turn) resid
    drop r_a r_b r_c r_d _reghdfe_resid ref_r
}
if _rc test_fail "HDFE-08/09 residuals spellings and reruns" "r(`=_rc')"
else test_pass "HDFE-08/09 residuals spellings and reruns"

* HDFE-06 and HDFE-07: absorb() suboptions; reruns drop __hdfe*
capture noisily {
    sysuse auto, clear
    quietly reghdfe price weight, absorb(turn foreign, savefe)
    rename (__hdfe1__ __hdfe2__) (ref1 ref2)
    quietly creghdfe price weight, absorb(turn foreign, savefe)
    quietly creghdfe price weight, absorb(turn foreign, savefe)
    assert reldif(__hdfe1__ + __hdfe2__, ref1 + ref2) < 1e-7 | missing(ref1)
    generate byte __hdfe_other = 1
    quietly creghdfe price weight, absorb(turn)
    capture describe __hdfe*
    assert _rc == 111
    foreach o in generate gen g {
        quietly creghdfe price weight, absorb(turn foreign, `o')
        assert reldif(__hdfe1__ + __hdfe2__, ref1 + ref2) < 1e-7 | missing(ref1)
    }
    quietly creghdfe price weight, absorb(turn, noisily)
    quietly creghdfe price weight, absorb(fe1=turn)
    capture creghdfe price weight, absorb(fe1=turn)
    assert _rc == 110
    replace fe1 = .
    quietly creghdfe price weight, absorb(fe1=turn foreign, savefe replace)
    confirm variable fe1 __hdfe2__
    assert fe1 < . if e(sample)
    * as in reghdfe: rep is not replace, and it is checked after the variables
    capture creghdfe price weight, absorb(fe1=turn, rep)
    assert _rc == 110
    capture creghdfe price weight, absorb(turn, rep)
    assert _rc == 198
    capture reghdfe price weight, absorb(turn, rep)
    assert _rc == 198
    generate fs = 1
    generate fsSlope1 = 1
    capture creghdfe price weight, absorb(fs=turn##c.length, replace)
    assert _rc == 110
    quietly bootstrap, reps(3) seed(1) nodots: creghdfe price weight, absorb(turn, savefe)
}
if _rc test_fail "HDFE-06/07 absorb() suboptions and savefe reruns" "r(`=_rc')"
else test_pass "HDFE-06/07 absorb() suboptions and savefe reruns"

* HDFE-11: vce() types as reghdfe reads them
sysuse auto, clear
local fails
foreach v in rob r robust un unad ols olsx conventional rubbish "cl turn" "clu turn" {
    capture noisily {
        quietly reghdfe price weight, absorb(turn) vce(`v')
        tempname V
        matrix `V' = e(V)
        local vce `e(vce)'
        quietly creghdfe price weight, absorb(turn) vce(`v')
        assert mreldif(e(V), `V') < 1e-9 & "`e(vce)'" == "`vce'"
    }
    if _rc local fails `"`fails' [`v']"'
}
foreach v in conv Robust {
    capture creghdfe price weight, absorb(turn) vce(`v')
    local rc = _rc
    capture reghdfe price weight, absorb(turn) vce(`v')
    if _rc != `rc' | `rc' != 9 local fails `"`fails' [`v': `rc']"'
}
if `"`fails'"' == "" test_pass "HDFE-11 vce() abbreviations"
else test_fail "HDFE-11 vce() abbreviations" `"`fails'"'

* HDFE-12, HDFE-13, HDFE-14, HDFE-21: nosample, timeit, version(), nowarn,
* nohe
capture noisily {
    quietly creghdfe price weight, absorb(turn) nosample
    quietly count if e(sample)
    assert r(N) == 0 & e(N) == 70
    tempvar xb
    quietly predict double `xb', xb
    quietly creghdfe price weight, absorb(turn) nosamp
    quietly creghdfe price weight, absorb(turn) timeit
    matrix t = r(timer)
    assert rowsof(t) == 4 & colsof(t) == 2
    assert "`: colnames t'" == "Time (% Total)"
    capture creghdfe price weight, absorb(turn) version(4)
    assert _rc == 9
    capture creghdfe price weight, absorb(turn) version(5)
    assert _rc == 198
    _s2h_log w0 : creghdfe price weight, absorb(turn) keepsingletons nowarn notable nofootnote
    assert !strpos(`"`w0'"', "WARNING")
    _s2h_log w1 : creghdfe price weight, absorb(turn) keepsingletons notable nofootnote
    assert strpos(`"`w1'"', "WARNING: Singleton observations not dropped")
    quietly creghdfe price weight, absorb(turn) nohe
    quietly creghdfe price weight, absorb(turn) nohea
}
if _rc test_fail "HDFE-12/13/14/21 nosample, timeit, version(), nowarn, nohe" "r(`=_rc')"
else test_pass "HDFE-12/13/14/21 nosample, timeit, version(), nowarn, nohe"

* HDFE-15: solver options, as reghdfe parses them; transform() and
* acceleration() run reghdfe's MAP variants
local fails
foreach o in "technique(map)" "technique(lsmr)" "tech(lsqr)" "transform(cim)" "transform(kac)" ///
    "tra(sym)" "acceleration(sd)" "accel(none)" "acceleration(off)" "acceleration(aitken)" ///
    "transform(cim) accel(none)" "preconditioner(diag)" "precond(none)" "prune" "poolsize(5)" ///
    "pool(0)" "compact" "fastregress" "parallel(1)" "transform(kac) acceleration(cg)" {
    capture noisily {
        quietly reghdfe price weight length, absorb(turn foreign) `o' tolerance(1e-12)
        tempname b
        matrix `b' = e(b)
        quietly creghdfe price weight length, absorb(turn foreign) `o' tolerance(1e-12)
        assert mreldif(e(b), `b') < 1e-8 & e(converged) == 1
    }
    if _rc local fails `"`fails' [`o']"'
}
foreach c in "3498|technique(cg)" "9|technique(foo)" "9|transform(foo)" "3498|transform(rand)" ///
    "9|acceleration(foo)" "9|preconditioner(foo)" {
    gettoken r o : c, parse("|")
    local o : subinstr local o "|" ""
    capture creghdfe price weight, absorb(turn foreign) `o'
    local rc = _rc
    capture reghdfe price weight, absorb(turn foreign) `o'
    if `rc' != `r' | _rc != `r' local fails `"`fails' [`o': `rc']"'
}
_s2h_log w : creghdfe price weight, absorb(turn foreign) transform(kac) notable nofootnote
if !strpos(`"`w'"', "(WARNING: convergence is unlikely with transform=kaczmarz and accel=CG)") ///
    local fails `"`fails' [kaczmarz warning]"'
if `"`fails'"' == "" test_pass "HDFE-15 solver options"
else test_fail "HDFE-15 solver options" `"`fails'"'

* HDFE-16: abort(0) keeps the projection after iterate() iterations
clear
set obs 3000
generate id1 = runiformint(1, 300)
generate id2 = runiformint(1, 200)
generate x = rnormal() + id1/100
generate y = x + id1/50 + id2/30 + rnormal()
capture noisily {
    quietly reghdfe y x, absorb(id1 id2) iterate(2) abort(0)
    local bref = _b[x]
    assert e(converged) == 0
    quietly creghdfe y x, absorb(id1 id2) iterate(2) abort(0)
    assert e(converged) == 0 & e(ic) == 0
    assert reldif(_b[x], `bref') < 1e-6
    capture creghdfe y x, absorb(id1 id2) iterate(2) abort(1)
    assert _rc == 430
    capture creghdfe y x, absorb(id1 id2) iterate(2)
    assert _rc == 430
    quietly creghdfe y x, absorb(id1 id2) miniter(50) min_ok(2) accel_start(3) accel_freq(3)
    assert e(converged) == 1
    foreach o in "foo(1)" "foo" {
        capture creghdfe y x, absorb(id1 id2) `o'
        local rc = _rc
        capture reghdfe y x, absorb(id1 id2) `o'
        assert `rc' == 198 & _rc == 198
    }
}
if _rc test_fail "HDFE-16 abort(0) and the other solver properties" "r(`=_rc')"
else test_pass "HDFE-16 abort(0) and the other solver properties"

* HDFE-18: reghdfe's programming options are refused
local fails
foreach o in nopartialout noregress keepmata "nopartialout varlist_is_touse" {
    capture creghdfe y x, absorb(id1 id2) `o'
    if _rc != 198 local fails `"`fails' [`o': `=_rc']"'
}
if `"`fails'"' == "" test_pass "HDFE-18 programming options refused"
else test_fail "HDFE-18 programming options refused" `"`fails'"'

* HDFE-19 and HDFE-20: a factor-variable depvar, ib#. and ibn. in absorb()
sysuse auto, clear
_s2h_cmp, ref(reghdfe 1.foreign weight, absorb(turn)) ct(creghdfe 1.foreign weight, absorb(turn)) ///
    name("HDFE-19 factor-variable depvar") fit
capture noisily {
    quietly creghdfe 1.foreign weight, absorb(turn)
    assert "`e(depvar)'" == "1.foreign"
    capture creghdfe i.foreign weight, absorb(turn)
    assert _rc == 9
}
if _rc test_fail "HDFE-19 factor-variable depvar names" "r(`=_rc')"
else test_pass "HDFE-19 factor-variable depvar names"
_s2h_cmp, ref(reghdfe price weight, absorb(ib2.rep78)) ct(creghdfe price weight, absorb(ib2.rep78)) ///
    name("HDFE-20 ib#. in absorb()") fit
_s2h_cmp, ref(reghdfe price weight, absorb(ibn.rep78 turn)) ct(creghdfe price weight, absorb(ibn.rep78 turn)) ///
    name("HDFE-20 ibn. in absorb()") fit
_s2h_cmp, ref(reghdfe price weight, absorb(ib2.rep78#foreign)) ct(creghdfe price weight, absorb(ib2.rep78#foreign)) ///
    name("HDFE-20 ib#. in an absorbed interaction") fit

* HDFE-22 and HDFE-23: the header and the absorbed degrees of freedom as
* reghdfe shows them, and e(dof_table) with reghdfe's names
sysuse auto, clear
generate t = _n
local fails
foreach m in "price weight length, absorb(turn)" "price weight length, absorb(turn) vce(robust)" ///
    "price weight length, absorb(turn) vce(cluster turn)" "price weight, absorb(turn) vce(cluster turn foreign)" ///
    "price weight length, absorb(turn##c.length foreign#rep78)" ///
    "price weight length, absorb(turn foreign rep78) vce(cluster rep78)" "price weight length," ///
    "price weight length, absorb(turn) dkraay" {
    if strpos("`m'", "dkraay") {
        quietly tsset t
        local m price weight length, absorb(turn) vce(dkraay 2)
    }
    _s2h_log ref : reghdfe `m' notable
    _s2h_log got : creghdfe `m' notable
    if `"`ref'"' != `"`got'"' {
        di as text `"reghdfe:  `ref'"'
        di as text `"creghdfe: `got'"'
        local fails `"`fails' [`m']"'
    }
    quietly reghdfe `m'
    local rr : rowfullnames e(dof_table)
    local rc : colfullnames e(dof_table)
    quietly creghdfe `m'
    if `"`rr'"' != `"`: rowfullnames e(dof_table)'"' | `"`rc'"' != `"`: colfullnames e(dof_table)'"' {
        local fails `"`fails' [names `m']"'
    }
    capture tsset, clear
}
if `"`fails'"' == "" test_pass "HDFE-22/23 header, footnote and e(dof_table) names"
else test_fail "HDFE-22/23 header, footnote and e(dof_table) names" `"`fails'"'

* HDFE-24, HDFE-25: estat (e(estat_cmd)), e(footnote), e(N_full), e(title3),
* e(title4), e(indepvars)
capture noisily {
    sysuse auto, clear
    foreach m in "absorb(turn)" "absorb(turn) vce(cluster turn)" {
        quietly reghdfe price weight length, `m'
        quietly estat ic
        tempname S
        matrix `S' = r(S)
        local nfull = e(N_full)
        local t3 `"`e(title3)'"'
        quietly creghdfe price weight length, `m'
        assert "`e(estat_cmd)'" == "_creghdfe_estat" & "`e(footnote)'" == "_creghdfe_footnote"
        assert e(N_full) == `nfull' & `"`e(title3)'"' == `"`t3'"'
        quietly estat ic
        assert mreldif(r(S), `S') < 1e-10
        quietly estat summarize
        quietly estat vce
        capture estat foo
        assert _rc == 321
    }
    quietly reghdfe price weight i.rep78 i.foreign#c.length, absorb(turn)
    local ind `e(indepvars)'
    quietly creghdfe price weight i.rep78 i.foreign#c.length, absorb(turn)
    assert "`e(indepvars)'" == "`ind'"
    quietly reghdfe price, absorb(turn)
    local ind `e(indepvars)'
    quietly creghdfe price, absorb(turn)
    assert "`e(indepvars)'" == "`ind'"
    generate t = _n
    tsset t
    quietly reghdfe price weight, absorb(turn) vce(dkraay 2)
    local t4 `"`e(title4)'"'
    quietly creghdfe price weight, absorb(turn) vce(dkraay 2)
    assert `"`e(title4)'"' == `"`t4'"'
    _s2h_log foot : _creghdfe_footnote
    assert strpos(`"`foot'"', "Absorbed degrees of freedom")
    quietly bootstrap, reps(3) seed(2) nodots notable: creghdfe price weight, absorb(turn)
}
if _rc test_fail "HDFE-24/25 estat and stored results" "r(`=_rc')"
else test_pass "HDFE-24/25 estat and stored results"

* HDFE-26: predict labels without absorb()
capture noisily {
    sysuse auto, clear
    quietly reghdfe price weight, resid
    predict double d_ref, d
    predict double xbd_ref, xbd
    drop _reghdfe_resid
    quietly creghdfe price weight, resid
    predict double d_got, d
    predict double xbd_got, xbd
    assert `"`: variable label d_got'"' == `"`: variable label d_ref'"'
    assert `"`: variable label xbd_got'"' == `"`: variable label xbd_ref'"'
    assert `"`: variable label d_got'"' == "d[_cons]"
}
if _rc test_fail "HDFE-26 predict labels without absorb()" "r(`=_rc')"
else test_pass "HDFE-26 predict labels without absorb()"

print_summary s2_hdfe
if $TESTS_FAILED > 0 exit 9
global CTOOLS_COMPONENT_COMPLETE "s2_hdfe"
    global CTOOLS_CASE_COMPLETE "s2_hdfe"
    exit
}

di as error "Unknown s2_hdfe case: `__case'"
exit 198
