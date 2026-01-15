/*******************************************************************************
 * validate_cimport.do
 *
 * Comprehensive validation tests for cimport delimited vs import delimited
 * Tests all import options: delimiters, varnames, case, bindquotes, stripquotes
 *
 * VERIFICATION: All tests use cf _all to verify byte-for-byte identical data
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

di as text ""
di as text "======================================================================"
di as text "              CIMPORT DELIMITED VALIDATION TEST SUITE"
di as text "======================================================================"

capture mkdir "temp"

/*******************************************************************************
 * Helper: compare_imports - Import with both methods and compare using cf _all
 * Returns 0 if identical, nonzero otherwise
 ******************************************************************************/
capture program drop compare_imports
program define compare_imports, rclass
    syntax using/, [STATAopts(string) CIMPORTopts(string) testname(string)]

    * Import with Stata's import delimited
    import delimited `using', `stataopts' clear
    tempfile stata_data
    quietly save `stata_data', replace
    local stata_n = _N
    local stata_k = c(k)

    * Import with cimport delimited
    cimport delimited `using', `cimportopts' clear
    tempfile cimport_data
    quietly save `cimport_data', replace
    local cimport_n = _N
    local cimport_k = c(k)

    * Check dimensions first
    if `stata_n' != `cimport_n' | `stata_k' != `cimport_k' {
        return local result "fail"
        return local reason "dimensions differ: Stata N=`stata_n' K=`stata_k', cimport N=`cimport_n' K=`cimport_k'"
        exit
    }

    * Compare data using cf _all
    use `stata_data', clear
    capture cf _all using `cimport_data'
    if _rc == 0 {
        return local result "pass"
        return local reason ""
    }
    else {
        * Get more detail on what differs
        return local result "fail"
        return local reason "cf _all comparison failed - data not identical"
    }
end

/*******************************************************************************
 * Create test CSV files
 ******************************************************************************/
noi print_section "Creating Test Data Files"

capture file close _all

* Basic CSV with headers
file open fh using "temp/basic.csv", write replace
file write fh "id,name,value,category" _n
file write fh "1,Alpha,100.5,A" _n
file write fh "2,Beta,200.25,B" _n
file write fh "3,Gamma,300.75,A" _n
file write fh "4,Delta,400.1,C" _n
file write fh "5,Epsilon,500.9,B" _n
file close fh

* Tab-delimited file - use char(9) to ensure actual tab characters
local tab = char(9)
file open fh using "temp/tabfile.tsv", write replace
file write fh "id`tab'name`tab'score`tab'group" _n
file write fh "1`tab'John`tab'85.5`tab'X" _n
file write fh "2`tab'Jane`tab'92.3`tab'Y" _n
file write fh "3`tab'Bob`tab'78.1`tab'X" _n
file write fh "4`tab'Alice`tab'95.7`tab'Z" _n
file close fh

* Semicolon-delimited file
file open fh using "temp/semicolon.csv", write replace
file write fh "code;description;amount" _n
file write fh "A001;Widget;1234.56" _n
file write fh "B002;Gadget;2345.67" _n
file write fh "C003;Gizmo;3456.78" _n
file close fh

* No header file
file open fh using "temp/noheader.csv", write replace
file write fh "1,Alpha,100" _n
file write fh "2,Beta,200" _n
file write fh "3,Gamma,300" _n
file close fh

* File with missing values
file open fh using "temp/missing.csv", write replace
file write fh "id,x,y,z" _n
file write fh "1,10,20,30" _n
file write fh "2,,22," _n
file write fh "3,30,,33" _n
file write fh "4,,,44" _n
file write fh "5,50,52,55" _n
file close fh

* File with quoted fields
file open fh using "temp/quoted.csv", write replace
file write fh `"id,text,value"' _n
file write fh `"1,"Hello, World",100"' _n
file write fh `"2,"Embedded ""quotes""",200"' _n
file write fh `"3,Normal text,300"' _n
file close fh

* Large file
file open fh using "temp/large.csv", write replace
file write fh "id,group,x,y,label" _n
forvalues i = 1/5000 {
    local group = mod(`i'-1, 100) + 1
    local x = runiform() * 100
    local y = runiformint(1, 1000)
    file write fh "`i',`group',`x',`y',Item`i'" _n
}
file close fh

* Case test file
file open fh using "temp/casetest.csv", write replace
file write fh "FirstName,LastName,AGE,mixed_Case" _n
file write fh "John,Doe,25,abc" _n
file write fh "Jane,Smith,30,def" _n
file close fh

noi test_pass "Test files created"

/*******************************************************************************
 * SECTION 1: Plugin check
 ******************************************************************************/
noi print_section "Plugin Check"

capture cimport delimited using "temp/basic.csv", clear
if _rc != 0 {
    noi test_fail "cimport plugin load" "returned error `=_rc'"
    noi print_summary "cimport"
    exit 1
}
noi test_pass "cimport plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Basic comma-delimited import (cf _all verification)
 ******************************************************************************/
noi print_section "Basic Comma-Delimited Import"

compare_imports using "temp/basic.csv"
if "`r(result)'" == "pass" {
    noi test_pass "basic CSV data identical (cf _all)"
}
else {
    noi test_fail "basic CSV" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 3: Tab-delimited import (cf _all verification)
 ******************************************************************************/
noi print_section "Tab-Delimited Import"

compare_imports using "temp/tabfile.tsv", stataopts(delimiters(tab)) cimportopts(delimiters(tab))
if "`r(result)'" == "pass" {
    noi test_pass "tab-delimited data identical (cf _all)"
}
else {
    noi test_fail "tab-delimited" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 4: Semicolon-delimited import (cf _all verification)
 ******************************************************************************/
noi print_section "Semicolon-Delimited Import"

compare_imports using "temp/semicolon.csv", stataopts(delimiters(";")) cimportopts(delimiters(";"))
if "`r(result)'" == "pass" {
    noi test_pass "semicolon-delimited data identical (cf _all)"
}
else {
    noi test_fail "semicolon-delimited" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 5: varnames() option (cf _all verification)
 ******************************************************************************/
noi print_section "varnames() Option"

* varnames(1)
compare_imports using "temp/basic.csv", stataopts(varnames(1)) cimportopts(varnames(1))
if "`r(result)'" == "pass" {
    noi test_pass "varnames(1) data identical (cf _all)"
}
else {
    noi test_fail "varnames(1)" "`r(reason)'"
}

* varnames(nonames)
compare_imports using "temp/noheader.csv", stataopts(varnames(nonames)) cimportopts(varnames(nonames))
if "`r(result)'" == "pass" {
    noi test_pass "varnames(nonames) data identical (cf _all)"
}
else {
    noi test_fail "varnames(nonames)" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 6: case() option
 ******************************************************************************/
noi print_section "case() Option"

* case(preserve) - compare with Stata's case(preserve)
compare_imports using "temp/casetest.csv", stataopts(case(preserve)) cimportopts(case(preserve))
if "`r(result)'" == "pass" {
    noi test_pass "case(preserve) data identical (cf _all)"
}
else {
    noi test_fail "case(preserve)" "`r(reason)'"
}

* case(lower) - compare with Stata's case(lower)
compare_imports using "temp/casetest.csv", stataopts(case(lower)) cimportopts(case(lower))
if "`r(result)'" == "pass" {
    noi test_pass "case(lower) data identical (cf _all)"
}
else {
    noi test_fail "case(lower)" "`r(reason)'"
}

* case(upper) - compare with Stata's case(upper)
compare_imports using "temp/casetest.csv", stataopts(case(upper)) cimportopts(case(upper))
if "`r(result)'" == "pass" {
    noi test_pass "case(upper) data identical (cf _all)"
}
else {
    noi test_fail "case(upper)" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 7: Missing values handling (cf _all verification)
 ******************************************************************************/
noi print_section "Missing Values"

compare_imports using "temp/missing.csv"
if "`r(result)'" == "pass" {
    noi test_pass "missing values data identical (cf _all)"
}
else {
    noi test_fail "missing values" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 8: Quoted fields (cf _all verification)
 ******************************************************************************/
noi print_section "Quoted Fields"

compare_imports using "temp/quoted.csv"
if "`r(result)'" == "pass" {
    noi test_pass "quoted fields data identical (cf _all)"
}
else {
    noi test_fail "quoted fields" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 9: bindquotes() option (cf _all verification)
 ******************************************************************************/
noi print_section "bindquotes() Option"

* bindquotes(strict)
compare_imports using "temp/quoted.csv", stataopts(bindquotes(strict)) cimportopts(bindquotes(strict))
if "`r(result)'" == "pass" {
    noi test_pass "bindquotes(strict) data identical (cf _all)"
}
else {
    noi test_fail "bindquotes(strict)" "`r(reason)'"
}

* bindquotes(loose)
compare_imports using "temp/quoted.csv", stataopts(bindquotes(loose)) cimportopts(bindquotes(loose))
if "`r(result)'" == "pass" {
    noi test_pass "bindquotes(loose) data identical (cf _all)"
}
else {
    noi test_fail "bindquotes(loose)" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 10: stripquotes option (cf _all verification)
 ******************************************************************************/
noi print_section "stripquotes Option"

* Note: Stata's import delimited doesn't have stripquotes, so just verify it works
capture cimport delimited using "temp/quoted.csv", stripquotes clear
if _rc == 0 {
    noi test_pass "stripquotes option accepted"
}
else {
    noi test_fail "stripquotes" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 11: Large file import (cf _all verification)
 ******************************************************************************/
noi print_section "Large File Import"

compare_imports using "temp/large.csv"
if "`r(result)'" == "pass" {
    noi test_pass "large file (5000 rows) data identical (cf _all)"
}
else {
    noi test_fail "large file" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 12: verbose option
 ******************************************************************************/
noi print_section "verbose Option"

capture cimport delimited using "temp/basic.csv", verbose clear
if _rc == 0 {
    noi test_pass "verbose option accepted"
}
else {
    noi test_fail "verbose option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 13: clear option behavior
 ******************************************************************************/
noi print_section "clear Option"

sysuse auto, clear
capture cimport delimited using "temp/basic.csv"
if _rc != 0 {
    noi test_pass "no clear: fails with data in memory"
}
else {
    noi test_fail "no clear" "should fail with data in memory"
}

capture cimport delimited using "temp/basic.csv", clear
if _rc == 0 {
    noi test_pass "clear option works"
}
else {
    noi test_fail "clear option" "failed"
}

/*******************************************************************************
 * SECTION 14: Round-trip tests (cf _all verification)
 ******************************************************************************/
noi print_section "Round-Trip Tests"

* auto dataset round-trip
sysuse auto, clear
export delimited using "temp/auto_rt.csv", replace

compare_imports using "temp/auto_rt.csv"
if "`r(result)'" == "pass" {
    noi test_pass "auto dataset round-trip identical (cf _all)"
}
else {
    noi test_fail "auto round-trip" "`r(reason)'"
}

* census dataset round-trip
sysuse census, clear
export delimited using "temp/census_rt.csv", replace

compare_imports using "temp/census_rt.csv"
if "`r(result)'" == "pass" {
    noi test_pass "census dataset round-trip identical (cf _all)"
}
else {
    noi test_fail "census round-trip" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 15: Additional data types
 ******************************************************************************/
noi print_section "Additional Data Types"

* Create file with various numeric formats
file open fh using "temp/numtypes.csv", write replace
file write fh "int_val,float_val,neg_val,zero_val,large_val" _n
file write fh "1,1.5,-10,0,1000000" _n
file write fh "2,2.75,-20,0,2000000" _n
file write fh "3,3.125,-30,0,3000000" _n
file close fh

compare_imports using "temp/numtypes.csv"
if "`r(result)'" == "pass" {
    noi test_pass "numeric types data identical (cf _all)"
}
else {
    noi test_fail "numeric types" "`r(reason)'"
}

* Create file with long strings
file open fh using "temp/longstr.csv", write replace
file write fh "id,short_str,long_str" _n
file write fh "1,abc,This is a longer string with spaces and stuff" _n
file write fh "2,def,Another longer string value here for testing" _n
file close fh

compare_imports using "temp/longstr.csv"
if "`r(result)'" == "pass" {
    noi test_pass "string types data identical (cf _all)"
}
else {
    noi test_fail "string types" "`r(reason)'"
}

/*******************************************************************************
 * Cleanup and summary
 ******************************************************************************/

local files : dir "temp" files "*.csv"
foreach f of local files {
    capture erase "temp/`f'"
}
local files : dir "temp" files "*.tsv"
foreach f of local files {
    capture erase "temp/`f'"
}

noi print_summary "cimport"

if $TESTS_FAILED > 0 {
    exit 1
}

}
