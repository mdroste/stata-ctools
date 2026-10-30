* oracle_tests: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "oracle_tests"

if "`__case'" == "coverage_helpers" {
/* Mutations must be rejected even at small scales or when both calls error. */
do validation/support/validate_setup.do
do validation/support/coverage_helpers.do
coverage_fixture 1
mata: assert(coverage_isnumeric_equal((1e-16,1),(9e-16,1))==0)
mata: assert(coverage_isnumeric_equal((.a,.),(.b,.))==0)
mata: assert(coverage_isnumeric_equal((1,2),(1\2))==0)
mata: assert(coverage_isnumeric_equal((0,1),(0,1))==1)
mata: assert(coverage_isnumeric_equal((1e12),(1e12+1))==1)

program coverage_bad_data
    args mutation
    if "`mutation'"=="value" replace y=0 in 1
    if "`mutation'"=="missing" replace y=.a in 1
    if "`mutation'"=="dropvar" drop y
    if "`mutation'"=="dropobs" drop in 1
    if "`mutation'"=="name" rename y changed_y
    if "`mutation'"=="order" order y
    if "`mutation'"=="label" label var y "wrong label"
    if "`mutation'"=="valuelabel" label define groups 1 "wrong value label", modify
    if "`mutation'"=="note" notes: incorrect extra note
    if "`mutation'"=="format" format y %21.6f
    if "`mutation'"=="roworder" gsort -id
    if "`mutation'"=="error" exit 198
end
foreach mutation in value missing dropvar dropobs name order label valuelabel note format roworder error {
    local before=$TESTS_FAILED
    coverage_data, reference(coverage_bad_data correct) candidate(coverage_bad_data `mutation') name("mutation `mutation'")
    assert $TESTS_FAILED==`before'+1
}
local before=$TESTS_FAILED
coverage_data, reference(coverage_bad_data error) candidate(coverage_bad_data error) name("two failing commands")
assert $TESTS_FAILED==`before'+1
local before=$TESTS_FAILED
coverage_data, reference(coverage_bad_data correct) candidate(coverage_bad_data correct) name("identity")
assert $TESTS_FAILED==`before'

// Empty label mappings and empty datasets must pass identity comparisons.
foreach state in no_labels no_observations {
    preserve
    label values group
    if "`state'"=="no_observations" keep if 0
    local before=$TESTS_FAILED
    coverage_data, reference(coverage_bad_data correct) candidate(coverage_bad_data correct) name("identity `state'")
    assert $TESTS_FAILED==`before'
    restore
}
local before=$TESTS_FAILED
coverage_data, reference(coverage_export_readback, ctools(0) kind(delimited)) candidate(coverage_export_readback, ctools(0) kind(delimited)) name("native CSV identity")
assert $TESTS_FAILED==`before'

program coverage_bad_estimate, eclass
    args mutation
    quietly regress y x z
    tempname b V
    matrix `b'=e(b)
    matrix `V'=e(V)
    if "`mutation'"=="coefficient" matrix `b'[1,1]=0
    if "`mutation'"=="covariance" matrix `V'[1,1]=10*`V'[1,1]
    if "`mutation'"=="names" matrix colnames `b'=wrong z _cons
    if "`mutation'"=="sample" {
        tempvar sample
        gen byte `sample'=e(sample)
        replace `sample'=0 in 1
        ereturn repost b=`b' V=`V', esample(`sample') rename
    }
    else ereturn repost b=`b' V=`V', rename
    if "`mutation'"=="count" ereturn scalar N=e(N)-1
end
foreach mutation in coefficient covariance names sample count {
    local before=$TESTS_FAILED
    coverage_estimate, reference(regress y x z) candidate(coverage_bad_estimate `mutation') name("estimation mutation `mutation'")
    assert $TESTS_FAILED>`before'
}
local before=$TESTS_FAILED
coverage_estimate, reference(regress y x z) candidate(regress y x z) name("required absent diagnostic") required("absent_statistic")
assert $TESTS_FAILED>`before'
local before=$TESTS_FAILED
coverage_display, reference(regress y x z, level(90)) candidate(regress y x z, level(90)) name("native display identity") has("90%") table
assert $TESTS_FAILED==`before'
coverage_display, reference(regress y x z, level(90)) candidate(regress y x z, level(95)) name("wrong confidence level") has("90%") table
assert $TESTS_FAILED>`before'
// The precision adapter must reproduce native results and CSV bytes exactly.
coverage_prepare_bins_reference
forvalues variant=1/3 {
    coverage_fixture `variant'
    tempfile native_csv exposed_csv
    binscatter y x, controls(z) by(group) nquantiles(10) nodraw savedata("`native_csv'") replace
    local native_N=e(N)
    matrix native_coefs=e(y1_coefs)
    coverage_binscatter y x, controls(z) by(group) nquantiles(10) nodraw savedata("`exposed_csv'") replace
    assert e(N)==`native_N'
    matrix exposed_coefs=e(y1_coefs)
    mata: assert(st_matrix("native_coefs")==st_matrix("exposed_coefs"))
    mata: fh=fopen(st_local("native_csv")+".csv","r"); native_bytes=fread(fh,1000000); fclose(fh)
    mata: fh=fopen(st_local("exposed_csv")+".csv","r"); exposed_bytes=fread(fh,1000000); fclose(fh)
    mata: assert(native_bytes==exposed_bytes)
    erase "`native_csv'.csv"
    erase "`native_csv'.do"
    erase "`exposed_csv'.csv"
    erase "`exposed_csv'.do"
}
di "COVERAGE_HELPERS_MUTATIONS_PASS"
    global CTOOLS_CASE_COMPLETE "coverage_helpers"
    exit
}

if "`__case'" == "option_contracts" {
/* Public options without a one-command native equivalent: explicit contracts. */
do validation/support/validate_setup.do
do validation/support/coverage_helpers.do
quietly {
/* Finite sampling distributions: the bound is a predeclared chi-square
   significance level, not a numerical tolerance or a same-seed draw comparison. */
foreach cmd in sample csample bsample cbsample {
    mata: inclusion=J(20,1,0)
    forvalues rep=1/300 {
        clear
        set obs 20
        gen id=_n
        set seed `=84302+`rep''
        if "`cmd'"=="sample" sample 5, count
        else if "`cmd'"=="csample" csample, count(5)
        else `cmd' 5
        assert _N==5
        assert inrange(id,1,20) & id==floor(id)
        if inlist("`cmd'","sample","csample") isid id
        mata: ids=st_data(.,"id"); for(j=1;j<=rows(ids);j++) inclusion[ids[j]]=inclusion[ids[j]]+1
    }
    mata: st_numscalar("coverage_chisq",sum((inclusion:-75):^2:/75))
    assert scalar(coverage_chisq)<invchi2(19,.999999)
    test_pass "`cmd' repeated-seed inclusion distribution and exact counts"
}

/* These parsed options explicitly reject unsupported interfaces. Such a test
   proves rejection, not output equivalence or feature implementation. */
coverage_fixture 1
foreach command in "cqreg y x, absorb(g)" ///
    "civreghdfe y x (endog=z z2), center" ///
    "cbinscatter y x, genxq(unimplemented) nograph" {
    tempfile before
    save `before', replace
    capture `command'
    assert _rc==198
    cf _all using `before'
    test_pass "unsupported interface rejected: `command'"
}
cppmlhdfe count_y x z, absorb(g) d(saved_d)
foreach option in cooksd hat likelihood {
    capture predict double unsupported, `option'
    assert _rc==198
    capture confirm variable unsupported
    assert _rc==111
    test_pass "PPML predict `option' rejects without output"
}

/* Diagnostic output must be meaningful, and saved bins must equal the
   computed bins including their counts, not merely exist on disk. */
forvalues variant=1/3 {
    coverage_fixture `variant'
    tempfile bins
    cbinscatter y x, by(group) nquantiles(10) savedata("`bins'") nograph
    matrix computed=e(bindata)
    preserve
    use `bins', clear
    mkmat by_group bin_id x_mean y_mean n_obs, matrix(saved)
    mata: expected=st_matrix("computed"); expected=select(expected,expected[,2]:<.); st_matrix("computed",expected)
    matrix_min_sigfigs computed saved
    assert r(min_sigfigs)>=$DEFAULT_SIGFIGS
    restore
    test_pass "savedata all bin values and counts variant `variant'"
}

/* Check actual rendered marker shapes/colors/text, not acceptance of options. */
forvalues variant=1/3 {
    coverage_fixture `variant'
    coverage_graph, reference(binscatter y x, nquantiles(10) title(COVERAGE_TITLE) xtitle(COVERAGE_X) ytitle(COVERAGE_Y) colors(red blue)) ///
        candidate(cbinscatter y x, nquantiles(10) title(COVERAGE_TITLE) xtitle(COVERAGE_X) ytitle(COVERAGE_Y) colors(red blue)) ///
        name("graph titles colors `variant'") text("COVERAGE_TITLE COVERAGE_X COVERAGE_Y #FF0000 #0000FF")
    coverage_graph, reference(binscatter y x, nquantiles(10) msymbols(D)) ///
        candidate(cbinscatter y x, nquantiles(10) msymbols(D)) name("graph marker symbols `variant'")
    coverage_graph, reference(binscatter y x, nquantiles(10) by(group) legend(label(1 COVERAGE_LEGEND))) ///
        candidate(cbinscatter y x, nquantiles(10) by(group) legend(label(1 COVERAGE_LEGEND))) ///
        name("graph legend `variant'") text("COVERAGE_LEGEND")
    tempfile label_svg
    cbinscatter y x, nquantiles(10) mlabels(COVERAGE_MARKER_LABEL)
    graph export "`label_svg'.svg", as(svg) replace
    mata: fh=fopen(st_local("label_svg")+".svg","r"); svg=fread(fh,10000000); fclose(fh)
    mata: st_local("label_present",strofreal(strpos(svg,"COVERAGE_MARKER_LABEL")>0))
    if !`label_present' test_fail "mlabels rendered text `variant'" "requested marker label absent from SVG"
    else test_pass "mlabels rendered text `variant'"
    erase "`label_svg'.svg"
    foreach line in lfit qfit {
        local regressors x
        if "`line'"=="qfit" local regressors c.x##c.x
        regress y `regressors'
        local intercept=_b[_cons]
        local slope=_b[x]
        local fit_n=e(N)
        local fit_r2=e(r2)
        if "`line'"=="qfit" local quadratic=_b[c.x#c.x]
        cbinscatter y x, linetype(`line') reportreg nograph
        matrix fit=e(coefs)
        matrix stats=e(fit_stats)
        sigfigs `intercept' `=fit[1,1]'
        assert r(sigfigs)>=$DEFAULT_SIGFIGS
        sigfigs `slope' `=fit[1,2]'
        assert r(sigfigs)>=$DEFAULT_SIGFIGS
        if "`line'"=="qfit" {
            sigfigs `quadratic' `=fit[1,3]'
            assert r(sigfigs)>=$DEFAULT_SIGFIGS
        }
        assert stats[1,2]==`fit_n'
        sigfigs `fit_r2' `=stats[1,1]'
        assert r(sigfigs)>=$DEFAULT_SIGFIGS
        test_pass "`line' reportreg coefficients and fit `variant'"
    }
}

}

/* Presentation controls must change the actual displayed output. */
forvalues variant=1/3 {
    coverage_fixture `variant'
    foreach opt in "level(90)" "eform(Effect)" "title(COVERAGE_TITLE)" "subtitle(COVERAGE_SUBTITLE)" "depname(OUTCOME)" {
        local expected "90%"
        if strpos("`opt'","eform") local expected "Effect"
        if strpos("`opt'","title(") local expected "COVERAGE_TITLE"
        if strpos("`opt'","subtitle(") local expected "COVERAGE_SUBTITLE"
        if strpos("`opt'","depname") local expected "OUTCOME"
        local display_check has("`expected'")
        if strpos("`opt'","eform") | strpos("`opt'","depname") local display_check table
        local reference ivreghdfe y x (endog=z z2), absorb(g) `opt'
        if strpos("`opt'","eform") local reference coverage_iv_replay, estimate(ivreghdfe y x (endog=z z2), absorb(g)) display(`opt')
        coverage_display, reference(`reference') ///
            candidate(civreghdfe y x (endog=z z2), absorb(g) `opt') name("IV display `opt' `variant'") `display_check'
    }
    foreach opt in noheader nofooter nooutput noid {
        local absent "number of obs"
        if "`opt'"=="nofooter" local absent "instrumented:"
        if "`opt'"=="nooutput" local absent "conf. interval"
        if "`opt'"=="noid" local absent "underidentification test"
        coverage_display, reference(ivreghdfe y x (endog=z z2), absorb(g) `opt') ///
            candidate(civreghdfe y x (endog=z z2), absorb(g) `opt') name("IV display `opt' `variant'") absent("`absent'")
    }
    coverage_display, reference(ivreghdfe y x (endog=z z2), absorb(g) first) ///
        candidate(civreghdfe y x (endog=z z2), absorb(g) first) name("IV first-stage display `variant'") has("first.stage regression")
    coverage_display, reference(ivreghdfe y x (endog=z z2), absorb(g) rf) ///
        candidate(civreghdfe y x (endog=z z2), absorb(g) rf) name("IV reduced-form display `variant'") has("reduced.form regression")
    // Factor interactions, a collinear regressor, and an empty cell make these
    // formatting switches observable, instead of testing irrelevant OLS data.
    // h is drawn independently of the absorbed g: with the fixture's h, group
    // is a function of (g,h), and ivreghdfe 1.1.4 then misreports the
    // interaction coefficients (zeros with standard errors, missing standard
    // errors) where ivreg2 with partialled g dummies agrees with civreghdfe.
    gen double duplicate_x=x
    replace h=floor(15*runiform())
    drop if group==2 & h==1
    foreach opt in "plus" "noomitted" "omitted" "vsquish" "noemptycells" "baselevels" "allbaselevels" {
        coverage_display, reference(coverage_iv_replay, estimate(ivreghdfe y x duplicate_x i.group##i.h (endog=z z2), absorb(g)) display(`opt')) ///
            candidate(civreghdfe y x duplicate_x i.group##i.h (endog=z z2), absorb(g) `opt') ///
            name("IV factor display `opt' `variant'") table
    }
    graph drop _all
    cbinscatter y x, nograph
    capture graph describe
    if !_rc test_fail "nograph `variant'" "nograph unexpectedly creates a graph"
    else test_pass "nograph `variant'"
}

quietly {
/* Format switches must affect file bytes, not just equivalent read-back. */
forvalues variant=1/3 {
    coverage_fixture `variant'
    replace y=id-300
    replace s="two,fields" if mod(id,11)==0
    replace s=char(34)+"quoted"+char(34) if mod(id,13)==0
    format y %12.2f
    foreach opt in "quote" "novarnames" "nolabel" "delimiter(;)" "datafmt" {
        coverage_export_bytes, name("CSV formatting `opt' `variant'") opts(`opt')
    }
    coverage_export_bytes, name("CSV CRLF `variant'") crlf
    // The documented noquoteif extension deliberately emits raw string fields.
    preserve
    clear
    set obs 1
    gen str20 raw="has,comma"
    tempfile raw_csv
    cexport delimited raw using "`raw_csv'", novarnames noquoteif replace
    mata: fh=fopen(st_local("raw_csv"),"r"); raw_bytes=fread(fh,1000); fclose(fh)
    mata: st_local("raw_correct",strofreal(raw_bytes=="has,comma"+char(10)))
    if !`raw_correct' test_fail "noquoteif raw CSV `variant'" "field was quoted or escaped"
    else test_pass "noquoteif raw CSV `variant'"
    restore
}

/* Version and environment results must agree with the package and plugin. */
ctools, version
local version "`r(version)'"
local plugin "`r(plugin_version)'"
local revision "`r(plugin_revision)'"
assert "`version'"!="" & "`plugin'"!="" & "`revision'"!=""
ctools, environment_check
assert "`r(version)'"=="`version'" & "`r(plugin_version)'"=="`plugin'" & "`r(plugin_revision)'"=="`revision'"
coverage_data, reference(ctools) candidate(ctools, verbose) name("ctools verbose preserves data")
test_pass "ctools version and environment identity"

/* Exercise update dispatch without network or installation. Verify that the
   local net stub intercepts dispatch before invoking the update branch. */
tempfile update_copy
tempname source_file target_file
file open `source_file' using "build/ctools.ado", read text
file open `target_file' using "`update_copy'", write text replace
file read `source_file' line
local substitutions=0
while r(eof)==0 {
    if strpos(`"`macval(line)'"',"program define ctools,") {
        local line : subinstr local line "program define ctools," "program define coverage_ctools,"
        local ++substitutions
    }
    if strpos(`"`macval(line)'"',"net install ctools,") {
        local line : subinstr local line "net install ctools," "coverage_net install ctools,"
        local ++substitutions
    }
    file write `target_file' `"`macval(line)'"' _n
    file read `source_file' line
}
file close `source_file'
file close `target_file'
assert `substitutions'==2
do "`update_copy'"
program define coverage_net
    if `"`0'"'=="__coverage_probe" {
        global COVERAGE_NET_PROBED 1
        exit
    }
    assert strpos(`"`0'"',"install ctools,")==1
    assert strpos(`"`0'"',"https://raw.githubusercontent.com/mdroste/stata-ctools/main/build")>0
    assert strpos(`"`0'"',"replace")>0
    global COVERAGE_NET_CALLED 1
    exit $COVERAGE_NET_RC
end
global COVERAGE_NET_PROBED 0
coverage_net __coverage_probe
assert $COVERAGE_NET_PROBED==1
foreach rc in 0 631 198 {
    global COVERAGE_NET_CALLED 0
    global COVERAGE_NET_RC `rc'
    coverage_ctools, update
    assert $COVERAGE_NET_CALLED==1
    test_pass "ctools update dispatch mocked status `rc'"
}
program drop coverage_net coverage_ctools
macro drop COVERAGE_NET_PROBED COVERAGE_NET_CALLED COVERAGE_NET_RC

/* No-return is a timing extension: it must not publish a new estimation result. */
coverage_fixture 1
regress y z
matrix before_b=e(b)
local before_cmd "`e(cmd)'"
civreghdfe y x (endog=z z2), absorb(g) noreturn
assert "`e(cmd)'"=="`before_cmd'"
matrix_min_sigfigs before_b e(b)
assert r(min_sigfigs)==15
test_pass "IV noreturn preserves previous results"
}
print_summary "option_contracts"
global CTOOLS_COMPONENT_COMPLETE "option_contracts"
    global CTOOLS_CASE_COMPLETE "option_contracts"
    exit
}

if "`__case'" == "validation_helpers" {
/* Mutation checks for the validation harness. Called first by audit_p1, after
   validate_setup.do. Each bad implementation must record exactly one failure.
   Drop the stubs before returning so subsequent tests load the real ado files. */

capture program drop helper_expect
program define helper_expect
    args failures passes
    assert $TESTS_FAILED == `failures'
    assert $TESTS_PASSED == `passes'
    assert $TESTS_TOTAL == `failures' + `passes'
    reset_counters
end

assert_scalar_equal 1e-16 9e-16 7 "tiny values must retain relative precision"
helper_expect 1 0
assert_scalar_equal 1e-16 1e-16 7 "identical tiny values"
helper_expect 0 1
sigfigs 1e-16 1.000001e-16
assert abs(r(sigfigs) - 6) < 1e-5
matrix helper_a = (1e-16, 1)
matrix helper_b = (9e-16, 1)
assert_matrix_equal helper_a helper_b 7 "tiny matrix element differs"
helper_expect 1 0

/* Specialized cancellation checks must tolerate only machine roundoff.
   A small covariance in small units is still tested on its own scale. */
sigfigs_cancellation 0 1e-16 1
assert r(sigfigs) == 15
sigfigs_cancellation 0 1e-12 1
assert r(sigfigs) < 7
sigfigs_cancellation 1e-24 9e-24 1e-20
assert r(sigfigs) < 7
sigfigs_cancellation 0 . 1
assert r(sigfigs) < 7
sigfigs_cancellation 0 1e-16 0
assert r(sigfigs) < 7

/* Inference comparison rejects wrong or missing results from cpsmatch,
   including cases where the reference SE is zero or undefined. */
foreach se in 0 . 30 {
    assert_psmatch_inference `se' -.6666666666666667 3 -2 7 "corrupt matching SE"
    helper_expect 1 0
}
assert_psmatch_inference 3 0 3 -2 7 "corrupt matching t"
helper_expect 1 0
assert_psmatch_inference 3 -.6666666666666667 3 -2 7 "correct matching inference"
helper_expect 0 1
assert_psmatch_inference . . . -2 7 "singleton undefined SE"
helper_expect 0 1
assert_psmatch_inference 0 . . -2 7 "singleton false zero SE"
helper_expect 1 0
assert_psmatch_inference 0 . 0 0 7 "zero variance undefined t"
helper_expect 0 1
assert_psmatch_inference 0 0 0 0 7 "zero variance false t"
helper_expect 1 0
capture program drop csort
program define csort, rclass
    syntax varlist [, ALGorithm(string) STReam(integer 0) VERBOSE]
    assert `stream' == 2
    if "$HELPER_MUTATION" == "identity" {
        return scalar stream = 1
        exit
    }
    sort `varlist', stable
    if "$HELPER_MUTATION" == "corrupt" replace witness = -witness in 1
    if "$HELPER_MUTATION" == "dropvar" drop witness
    return scalar stream = ("$HELPER_MUTATION" != "ignorestream")
end
clear
set obs 3
gen key = 4 - _n
gen witness = _n
foreach mutation in identity corrupt dropvar ignorestream correct {
    global HELPER_MUTATION "`mutation'"
    benchmark_sort key, stream(2) testname("sort mutation: `mutation'")
    if "`mutation'" == "correct" helper_expect 0 1
    else helper_expect 1 0
}
program drop csort

capture program drop cimport
program define cimport
    gettoken format rest : 0
    import `format' `rest'
    if "$HELPER_MUTATION" == "reorder" gsort -id
    if "$HELPER_MUTATION" == "corrupt" replace value = -value in 1
end
tempfile csv
tempname fh
file open `fh' using "`csv'", write text replace
file write `fh' "id,value" _n "1,100" _n "2,200" _n "3,300" _n
file close `fh'
foreach mutation in reorder corrupt correct {
    global HELPER_MUTATION "`mutation'"
    benchmark_import using "`csv'", testname("import mutation: `mutation'")
    if "`mutation'" == "correct" helper_expect 0 1
    else helper_expect 1 0
}
program drop cimport

/* Use identical coefficients, covariance, and samples, varying only the one
   result being tested. A missing scalar or a tiny coefficient must not escape
   comparison just because all the other output agrees. */
capture program drop helper_ppml_result
program define helper_ppml_result, eclass
    tempname b V
    tempvar sample
    matrix `b' = (1e-16, 2)
    matrix colnames `b' = x _cons
    matrix `V' = (1, 1e-16 \ 1e-16, 1)
    matrix rownames `V' = x _cons
    matrix colnames `V' = x _cons
    if "$HELPER_MUTATION" == "coefficient" matrix `b'[1,1] = 9e-16
    if "$HELPER_MUTATION" == "covariance" {
        matrix `V'[1,2] = 9e-16
        matrix `V'[2,1] = 9e-16
    }
    if "$HELPER_MUTATION" == "covariance_above_roundoff" {
        matrix `V'[1,2] = 1e-12
        matrix `V'[2,1] = 1e-12
    }
    gen byte `sample' = 1
    ereturn post `b' `V', esample(`sample') obs(3)
    ereturn scalar ll = -10
    ereturn scalar r2_p = 0
    ereturn scalar N_clust = 3
    if "$HELPER_MUTATION" == "ll" ereturn scalar ll = .
    if "$HELPER_MUTATION" == "r2_p" ereturn scalar r2_p = .
    if "$HELPER_MUTATION" == "N_clust" ereturn scalar N_clust = .
    if "$HELPER_MUTATION" == "zero" ereturn scalar r2_p = .5
    if "$HELPER_MUTATION" == "r2_small" ereturn scalar r2_p = 1e-12
    if "$HELPER_MUTATION" == "r2_roundoff" ereturn scalar r2_p = -2e-15
end
capture program drop ppmlhdfe
program define ppmlhdfe, eclass
    local mutation "$HELPER_MUTATION"
    global HELPER_MUTATION "correct"
    helper_ppml_result
    global HELPER_MUTATION "`mutation'"
end
capture program drop cppmlhdfe
program define cppmlhdfe, eclass
    helper_ppml_result
end
clear
set obs 3
gen x = _n
gen y = _n
gen g = _n
foreach mutation in ll r2_p N_clust zero r2_small r2_roundoff coefficient covariance correct {
    global HELPER_MUTATION "`mutation'"
    benchmark_ppmlhdfe y x, absorb(g) testname("PPML mutation: `mutation'")
    if inlist("`mutation'", "correct", "r2_roundoff") helper_expect 0 1
    else helper_expect 1 0
}
/* Exercise the quantile helper itself: covariance roundoff may pass, but
   a 1e-12 covariance error or a tiny coefficient error must still fail. */
capture program drop qreg
program define qreg, eclass
    ppmlhdfe `0'
end
capture program drop cqreg
program define cqreg, eclass
    cppmlhdfe `0'
end
foreach mutation in covariance covariance_above_roundoff coefficient correct {
    global HELPER_MUTATION "`mutation'"
    benchmark_qreg y x, testname("quantile mutation: `mutation'")
    if inlist("`mutation'", "covariance", "correct") helper_expect 0 1
    else helper_expect 1 0
}
program drop qreg cqreg ppmlhdfe cppmlhdfe helper_ppml_result

/* A failing reference must not be turned into a nonmissing-value smoke test. */
capture program drop crangestat
program define crangestat
    if strpos("`0'", "__brp_ctools") gen __brp_ctools = 1
    else gen __brs_ctools = 1
end
capture program drop rangestat
program define rangestat
    exit 198
end
benchmark_rangestat mean x, interval(g -1 1) testname("failed range reference")
helper_expect 1 0
benchmark_range_pctile p25 x, interval(g -1 1) testname("incorrect window percentiles")
helper_expect 1 0
program drop crangestat rangestat helper_expect
macro drop HELPER_MUTATION
matrix drop helper_a helper_b
    global CTOOLS_CASE_COMPLETE "validation_helpers"
    exit
}

di as error "Unknown oracle_tests case: `__case'"
exit 198
