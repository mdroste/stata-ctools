* results_scope: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "if_rresults"

if "`__case'" == "if_rresults" {
* An if qualifier that names r() is evaluated with the caller's r(), as the
* native command evaluates it, and wrappers that the native command leaves
* r() alone for do not clear it.
do "validation/support/validate_setup.do"
quietly {
set seed 12
clear
set obs 2000
gen double w = runiform()
gen t = _n
gen double y = cond(runiform() < .3, ., rnormal())
gen str5 s = "k" + string(mod(_n, 7))
gen str9 s2 = "a b " + string(mod(_n, 5))
gen byte d = runiform() < .4
gen byte d2 = runiform() < .5
label define rr_dl 0 "z" 1 "o"
label values d rr_dl
label values d2 rr_dl
gen double x = rnormal()
gen double z = rnormal()
gen str4 n = string(mod(_n, 13))
tempfile base fn fc
save `base'

capture program drop rr_same
program define rr_same
    args a b label
    capture assert (`a' == `b') | (missing(`a') & missing(`b'))
    if _rc test_fail "`label'" "values differ from the native command"
    else test_pass "`label'"
end

foreach op in < > {
    foreach in in "" "in 1/1500" {
        use `base', clear
        summarize w, meanonly
        encode s if w`op'r(mean) `in', gen(ref)
        summarize w, meanonly
        cencode s if w`op'r(mean) `in', gen(got) label(rr_got)
        rr_same ref got "cencode if w`op'r(mean) `in'"

        use `base', clear
        summarize w, meanonly
        decode d if w`op'r(mean) `in', gen(ref)
        summarize w, meanonly
        cdecode d if w`op'r(mean) `in', gen(got)
        rr_same ref got "cdecode if w`op'r(mean) `in'"
        summarize w, meanonly
        cdecode d d2 if w`op'r(mean) `in', gen(got1 got2)
        rr_same ref got1 "cdecode two variables if w`op'r(mean) `in'"
    }

    use `base', clear
    summarize w, meanonly
    ipolate y t if w`op'r(mean), gen(ref)
    summarize w, meanonly
    cinterpolate y t if w`op'r(mean), gen(got)
    rr_same ref got "cinterpolate if w`op'r(mean)"

    use `base', clear
    summarize w, meanonly
    split s2 if w`op'r(mean), gen(ref)
    summarize w, meanonly
    csplit s2 if w`op'r(mean), gen(got)
    rr_same ref3 got3 "csplit if w`op'r(mean)"

    use `base', clear
    summarize w, meanonly
    local m = r(mean)
    cdestring n if w`op'`m', gen(ref)
    summarize w, meanonly
    cdestring n if w`op'r(mean), gen(got)
    rr_same ref got "cdestring if w`op'r(mean)"

    use `base', clear
    summarize w, meanonly
    export delimited t x using "`fn'.csv" if w`op'r(mean), replace
    summarize w, meanonly
    cexport delimited t x using "`fc'.csv" if w`op'r(mean), replace
    mata: st_local("same", strofreal(cat(st_local("fn") + ".csv") == cat(st_local("fc") + ".csv")))
    if `same' test_pass "cexport delimited if w`op'r(mean)"
    else test_fail "cexport delimited if w`op'r(mean)" "file differs from export delimited"

    use `base', clear
    summarize w, meanonly
    count if w`op'r(mean)
    local want = r(N)
    summarize w, meanonly
    cexport excel t x using "`fc'.xlsx" if w`op'r(mean), replace
    import excel using "`fc'.xlsx", clear
    if _N == `want' test_pass "cexport excel if w`op'r(mean)"
    else test_fail "cexport excel if w`op'r(mean)" "rows `=_N', expected `want'"

    capture which psmatch2
    if _rc test_skip "cpsmatch if w`op'r(mean)" "psmatch2 not installed"
    else {
        use `base', clear
        summarize w, meanonly
        psmatch2 d x z if w`op'r(mean), outcome(y)
        keep w t y s s2 d d2 x z n _support _pscore
        rename (_support _pscore) (ref_support ref_pscore)
        summarize w, meanonly
        cpsmatch d x z if w`op'r(mean), outcome(y)
        rr_same ref_support _support "cpsmatch if w`op'r(mean): support"
        rr_same ref_pscore _pscore "cpsmatch if w`op'r(mean): score"
    }
}

* Commands that leave r() as the native command does.
foreach cmd in "_ctools_load" "cbsample" "cbsample 100, strata(d)" {
    use `base', clear
    summarize x
    scalar rr_before = r(mean)
    `cmd'
    if r(mean) == scalar(rr_before) test_pass "`cmd' keeps the caller's r()"
    else test_fail "`cmd' keeps the caller's r()" "r(mean) is `r(mean)'"
}
_ctools_load
if "`__ctools_plugin_version'" != "" & "`__ctools_plugin_revision'" != "" test_pass "_ctools_load returns identity through locals"
else test_fail "_ctools_load returns identity through locals" "empty identity"
}
print_summary "if_rresults"
if $TESTS_FAILED exit 9
global CTOOLS_COMPONENT_COMPLETE "if_rresults"
    global CTOOLS_CASE_COMPLETE "if_rresults"
    exit
}

if "`__case'" == "estimator_margins" {
* margins after the estimators uses the posted weights, and refuses margins
* the data cannot identify, as after the reference commands.
do "validation/support/validate_setup.do"
quietly {
set seed 4
clear
set obs 3000
gen id = ceil(_n/10)
gen cat = 1 + mod(_n*7, 4)
gen r = 1 + mod(id, 3)
gen double x = rnormal() + cat/3
gen double z = rnormal()
gen double d = .5*z + .3*x + rnormal()
gen w = 1 + mod(_n*13, 5)
gen double y = 1 + .5*x + .2*x^2 + .3*cat + .4*d + rnormal() + mod(id,7)/3
gen cnt = rpoisson(exp(.2*x))

capture program drop em_hold
program define em_hold
    args stub
    matrix `stub'_b = r(b)
    matrix `stub'_V = r(V)
    matrix `stub'_t = r(table)
end
capture program drop em_compare
program define em_compare
    args label
    assert_matrix_equal got_b ref_b $DEFAULT_SIGFIGS "`label': margins"
    assert_matrix_equal got_V ref_V $DEFAULT_SIGFIGS "`label': variance"
    local same = 1
    forvalues j = 1/`=colsof(ref_t)' {
        if missing(got_t[1,`j']) != missing(ref_t[1,`j']) local same = 0
    }
    if `same' test_pass "`label': estimable pattern"
    else test_fail "`label': estimable pattern" "not-estimable margins differ"
end

foreach wt in "" "[aw=w]" "[fw=w]" "[pw=w]" {
    reghdfe y c.x##c.x i.cat `wt', absorb(id)
    margins, dydx(x)
    em_hold ref
    creghdfe y c.x##c.x i.cat `wt', absorb(id)
    if "`wt'" != "" {
        local wtype = word("aweight fweight pweight", strpos("afp", substr("`wt'", 2, 1)))
        if "`e(wtype)'" == "`wtype'" & "`e(wexp)'" == "= w" test_pass "creghdfe `wt' posts the weights"
        else test_fail "creghdfe `wt' posts the weights" "wtype=`e(wtype)' wexp=`e(wexp)'"
    }
    margins, dydx(x)
    em_hold got
    em_compare "creghdfe dydx `wt'"
    reghdfe y c.x##c.x i.cat `wt', absorb(id)
    margins cat
    em_hold ref
    creghdfe y c.x##c.x i.cat `wt', absorb(id)
    margins cat
    em_hold got
    em_compare "creghdfe margins cat `wt'"

    ivreghdfe y c.x##i.cat (d=z) `wt', absorb(id)
    margins, dydx(x)
    em_hold ref
    local ref_ok "`e(marginsok)'|`e(marginsnotok)'|`e(wtype)'|`e(wexp)'"
    civreghdfe y c.x##i.cat (d=z) `wt', absorb(id)
    if "`e(marginsok)'|`e(marginsnotok)'|`e(wtype)'|`e(wexp)'" == "`ref_ok'" test_pass "civreghdfe `wt' posts weights and margins permissions"
    else test_fail "civreghdfe `wt' posts weights and margins permissions" "`e(marginsok)'|`e(marginsnotok)'|`e(wtype)'|`e(wexp)' vs `ref_ok'"
    margins, dydx(x)
    em_hold got
    em_compare "civreghdfe dydx `wt'"
    ivreghdfe y c.x##i.cat (d=z) `wt', absorb(id)
    margins cat
    em_hold ref
    civreghdfe y c.x##i.cat (d=z) `wt', absorb(id)
    margins cat
    em_hold got
    em_compare "civreghdfe margins cat `wt'"
}

* r is constant within id, so its margins cannot be identified.
reghdfe y x i.r, absorb(id)
margins r
em_hold ref
creghdfe y x i.r, absorb(id)
margins r
em_hold got
em_compare "creghdfe non-estimable margins"
ppmlhdfe cnt x i.r, absorb(id)
margins r, predict(xb)
em_hold ref
cppmlhdfe cnt x i.r, absorb(id)
margins r, predict(xb)
em_hold got
em_compare "cppmlhdfe non-estimable margins"
* ivreghdfe partials out the constant and posts no estimability
* information, so its margins of r are numbers; civreghdfe matches it.
ivreghdfe y x i.r (d=z), absorb(id)
margins r
em_hold ref
civreghdfe y x i.r (d=z), absorb(id)
margins r
em_hold got
em_compare "civreghdfe margins as ivreghdfe"
}
print_summary "estimator_margins"
if $TESTS_FAILED exit 9
global CTOOLS_COMPONENT_COMPLETE "estimator_margins"
    global CTOOLS_CASE_COMPLETE "estimator_margins"
    exit
}

di as error "Unknown results_scope case: `__case'"
exit 198
