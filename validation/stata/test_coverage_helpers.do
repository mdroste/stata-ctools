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
