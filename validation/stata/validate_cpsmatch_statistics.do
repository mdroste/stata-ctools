/* Compare the new sufficient statistics to Stata on the actual matching
   outputs. Matching weights are held fixed, so dynamic kernel reductions
   cannot disguise a statistics discrepancy. Existing precision policy applies. */
capture do "validation/support/validate_setup.do"
if _rc do "validate_setup.do"

program define statistics_equal
    args actual expected name
    if `actual'==`expected' test_pass "`name'"
    else assert_scalar_equal `actual' `expected' $DEFAULT_SIGFIGS "`name'"
end

program define check_cpsmatch_statistics
    args name
    quietly summarize y if _treated==1 & _support==1
    local expected_mean=r(mean)
    local expected_n=r(N)
    local expected_var=r(Var)
    statistics_equal `=scalar(_cpsmatch_y1)' `expected_mean' "`name': supported treated mean"
    statistics_equal `=scalar(_cpsmatch_n1)' `expected_n' "`name': supported treated N"
    statistics_equal `=scalar(_cpsmatch_var1)' `expected_var' "`name': supported treated variance"
    quietly summarize y if _treated==0 & _support==1 [aw=_weight]
    local expected_mean=r(mean)
    local expected_n=r(N)
    local expected_var=r(Var)
    statistics_equal `=scalar(_cpsmatch_y0)' `expected_mean' "`name': weighted control mean"
    statistics_equal `=scalar(_cpsmatch_n0)' `expected_n' "`name': weighted control N"
    quietly summarize y if _treated==0 & _weight>0 & _weight<.
    statistics_equal `=scalar(_cpsmatch_var0)' `=r(Var)' "`name': unweighted used-control variance"
    tempvar w2
    quietly generate double `w2' = _weight^2 if _treated==0
    quietly summarize `w2', meanonly
    statistics_equal `=scalar(_cpsmatch_w2)' `=r(sum)' "`name': squared reuse weights"
    forvalues arm=0/1 {
        quietly summarize y if _treated==`arm'
        local expected_mean=r(mean)
        local expected_n=r(N)
        local expected_var=r(Var)
        statistics_equal `=scalar(_cpsmatch_y`arm'_u)' `expected_mean' "`name': unmatched mean `arm'"
        statistics_equal `=scalar(_cpsmatch_n`arm'_u)' `expected_n' "`name': unmatched N `arm'"
        statistics_equal `=scalar(_cpsmatch_var`arm'_u)' `expected_var' "`name': unmatched variance `arm'"
        quietly count if _treated==`arm' & _support==1
        statistics_equal `=scalar(_cpsmatch_n`arm'_on)' `=r(N)' "`name': support on `arm'"
        quietly count if _treated==`arm' & _support==0
        statistics_equal `=scalar(_cpsmatch_n`arm'_off)' `=r(N)' "`name': support off `arm'"
    }
end

foreach sample in ordinary shifted tiny large constant missing restricted ties noreplace common empty fractional_codes {
    clear
    set seed 881823
    quietly set obs 2000
    quietly gen long id = _n
    quietly gen double ps = runiform()
    quietly gen double treatment = mod(_n,2)
    quietly gen double y = rnormal()+2*treatment
    local options "neighbor(5)"
    local restriction ""
    if "`sample'"=="shifted" quietly replace y = 1e12 + y/1000
    if "`sample'"=="tiny" quietly replace y = y*1e-140
    if "`sample'"=="large" quietly replace y = y*1e140
    if "`sample'"=="constant" quietly replace y = 3
    if "`sample'"=="missing" {
        quietly replace y = .a if mod(id,17)==0
        quietly replace ps = . if mod(id,19)==0
    }
    if "`sample'"=="restricted" local restriction "if mod(id,3)!=0 in 11/1900"
    if "`sample'"=="ties" {
        quietly replace ps = round(ps,.05)
        local options "ties neighbor(3)"
    }
    if "`sample'"=="noreplace" local options "noreplacement descending"
    if "`sample'"=="common" {
        quietly replace ps = ps/2 if treatment==0
        local options "common caliper(.01)"
    }
    if "`sample'"=="empty" {
        quietly replace ps = treatment
        local options "caliper(.01)"
    }
    if "`sample'"=="fractional_codes" quietly replace treatment = .4 if treatment==1
    foreach team in 1 4 {
        quietly cpsmatch treatment `restriction', pscore(ps) outcome(y) `options' threads(`team')
        check_cpsmatch_statistics "`sample'-`team'"
    }
}

foreach kernel in epan normal biweight uniform tricube {
    clear
    set seed 18731
    quietly set obs 1000
    quietly gen double ps=runiform()
    quietly gen byte treatment=mod(_n,2)
    quietly gen double y=rnormal()+treatment
    quietly cpsmatch treatment, pscore(ps) outcome(y) kernel kerneltype(`kernel') threads(4)
    check_cpsmatch_statistics "kernel-`kernel'"
}
quietly cpsmatch treatment, pscore(ps) outcome(y) radius caliper(.1) threads(4)
check_cpsmatch_statistics "radius"

/* A call without outcome must refresh support counts and clear moment results. */
quietly cpsmatch treatment, pscore(ps) common threads(4)
assert missing(scalar(_cpsmatch_y1)) & missing(scalar(_cpsmatch_y0))
assert scalar(_cpsmatch_n1)==0 & scalar(_cpsmatch_n0)==0
forvalues arm=0/1 {
    quietly count if _treated==`arm' & _support==1
    assert scalar(_cpsmatch_n`arm'_on)==r(N)
}
test_pass "no-outcome statistics reset"
print_summary "cpsmatch_statistics"
di "CPSMATCH_STATISTICS_COMPLETE"
