/* Public options without a one-command native equivalent: explicit contracts. */
do validation/validate_setup.do
do validation/coverage_helpers.do
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
