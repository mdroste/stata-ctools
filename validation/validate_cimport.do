/*******************************************************************************
 * validate_cimport.do
 *
 * Comprehensive validation tests for cimport vs import delimited
 * Tests CSV import functionality across various scenarios
 ******************************************************************************/

do "validate_setup.do"

di as text ""
di as text "======================================================================"
di as text "              CIMPORT VALIDATION TEST SUITE"
di as text "======================================================================"

* Create temp directory for test files
capture mkdir "temp"

/*******************************************************************************
 * Create test CSV files
 ******************************************************************************/
print_section "Creating test files"

* Test file 1: Basic CSV with mixed types
capture file close _all
file open csvfile using "temp/basic.csv", write replace
file write csvfile "id,name,value,price" _n
file write csvfile "1,Alpha,100,10.5" _n
file write csvfile "2,Beta,200,20.75" _n
file write csvfile "3,Gamma,300,30.125" _n
file write csvfile "4,Delta,400,40.0" _n
file write csvfile "5,Epsilon,500,50.99" _n
file close csvfile

* Test file 2: Tab-delimited
file open tsvfile using "temp/basic.tsv", write replace
file write tsvfile "id" _tab "name" _tab "value" _n
file write tsvfile "1" _tab "Apple" _tab "1.5" _n
file write tsvfile "2" _tab "Banana" _tab "2.5" _n
file write tsvfile "3" _tab "Cherry" _tab "3.5" _n
file close tsvfile

* Test file 3: CSV with missing values
file open csvfile using "temp/missing.csv", write replace
file write csvfile "id,name,value,score" _n
file write csvfile "1,A,100,95" _n
file write csvfile "2,B,,85" _n
file write csvfile "3,,300," _n
file write csvfile "4,D,400,75" _n
file write csvfile "5,E,," _n
file close csvfile

* Test file 4: CSV with quoted fields
file open csvfile using "temp/quoted.csv", write replace
file write csvfile `"id,name,description"' _n
file write csvfile `"1,"Smith, John","A person with comma""' _n
file write csvfile `"2,"Doe, Jane","Another person""' _n
file write csvfile `"3,Simple,No quotes needed"' _n
file close csvfile

* Test file 5: CSV without header
file open csvfile using "temp/noheader.csv", write replace
file write csvfile "1,Alpha,100" _n
file write csvfile "2,Beta,200" _n
file write csvfile "3,Gamma,300" _n
file close csvfile

* Test file 6: CSV with various numeric formats
file open csvfile using "temp/numeric.csv", write replace
file write csvfile "integer,decimal,scientific,negative" _n
file write csvfile "1,1.5,1e10,-100" _n
file write csvfile "100,0.001,2.5e-5,-0.5" _n
file write csvfile "999999,123.456789,1e-10,-999999" _n
file close csvfile

* Test file 7: Large CSV
file open csvfile using "temp/large.csv", write replace
file write csvfile "id,group,value1,value2,label" _n
forvalues i = 1/5000 {
    local g = mod(`i', 100) + 1
    local v1 = `i' * 1.5
    local v2 = runiform() * 1000
    file write csvfile "`i',`g',`v1',`v2',item`i'" _n
}
file close csvfile

* Test file 8: CSV with UPPER case headers
file open csvfile using "temp/uppercase.csv", write replace
file write csvfile "ID,NAME,VALUE" _n
file write csvfile "1,Test,100" _n
file write csvfile "2,Data,200" _n
file close csvfile

di as result "[DONE] Test files created"

/*******************************************************************************
 * Basic import tests
 ******************************************************************************/
print_section "Basic import tests"

benchmark_import using "temp/basic.csv", testname("basic CSV")

benchmark_import using "temp/basic.tsv", delimiters(tab) testname("tab-delimited")

benchmark_import using "temp/missing.csv", testname("with missing values")

benchmark_import using "temp/quoted.csv", testname("quoted fields")

benchmark_import using "temp/noheader.csv", varnames(nonames) testname("no header")

benchmark_import using "temp/uppercase.csv", case(lower) testname("case(lower)")

benchmark_import using "temp/numeric.csv", testname("numeric formats")

benchmark_import using "temp/large.csv", testname("large file (5K rows)")

/*******************************************************************************
 * case() option tests
 ******************************************************************************/
print_section "case() option tests"

benchmark_import using "temp/uppercase.csv", case(upper) testname("case(upper)")

benchmark_import using "temp/uppercase.csv", case(preserve) testname("case(preserve)")

/*******************************************************************************
 * NOTE: rowrange() option is parsed but not yet implemented in cimport
 * Tests skipped until feature is implemented
 ******************************************************************************/

/*******************************************************************************
 * bindquotes() option tests
 ******************************************************************************/
print_section "bindquotes() option tests"

* Create a test file with embedded quotes
capture file close _all
file open csvfile using "temp/bindquotes_test.csv", write replace
file write csvfile `"id,name,value"' _n
file write csvfile `"1,"Simple Name",100"' _n
file write csvfile `"2,"Name with ""quotes""",200"' _n
file write csvfile `"3,NoQuotes,300"' _n
file close csvfile

benchmark_import using "temp/bindquotes_test.csv", testname("bindquotes default (strict)")

/*******************************************************************************
 * Round-trip tests (export then reimport)
 ******************************************************************************/
print_section "Round-trip tests"

sysuse auto, clear
export delimited using "temp/auto_export.csv", replace
benchmark_import using "temp/auto_export.csv", testname("auto round-trip")

sysuse census, clear
export delimited using "temp/census_export.csv", replace
benchmark_import using "temp/census_export.csv", testname("census round-trip")

/*******************************************************************************
 * Detailed comparison tests
 ******************************************************************************/
print_section "Detailed comparison tests"

* String variable types
import delimited using "temp/basic.csv", clear
local stata_type : type name

cimport delimited using "temp/basic.csv", clear
local cimport_type : type name

global TESTS_TOTAL = $TESTS_TOTAL + 1
if substr("`stata_type'", 1, 3) == "str" & substr("`cimport_type'", 1, 3) == "str" {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] string type inference"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] string type inference"
}

* Numeric type inference
import delimited using "temp/basic.csv", clear
local stata_id_type : type id
local stata_value_type : type value

cimport delimited using "temp/basic.csv", clear
local cimport_id_type : type id
local cimport_value_type : type value

local id_numeric = ("`stata_id_type'" != "str" & "`cimport_id_type'" != "str")
local value_numeric = ("`stata_value_type'" != "str" & "`cimport_value_type'" != "str")

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `id_numeric' & `value_numeric' {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] numeric type inference"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] numeric type inference"
}

* Missing value counts
import delimited using "temp/missing.csv", clear
quietly count if missing(value)
local stata_miss = r(N)

cimport delimited using "temp/missing.csv", clear
quietly count if missing(value)
local cimport_miss = r(N)

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `stata_miss' == `cimport_miss' {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] missing value handling"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] missing value handling"
}

* Large file statistics
import delimited using "temp/large.csv", clear
quietly summarize value1
local stata_mean = r(mean)
local stata_n = _N

cimport delimited using "temp/large.csv", clear
quietly summarize value1
local cimport_mean = r(mean)
local cimport_n = _N

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `stata_n' == `cimport_n' & abs(`stata_mean' - `cimport_mean') < 1e-6 {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] large file statistics"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] large file statistics"
}

/*******************************************************************************
 * Cleanup temp files
 ******************************************************************************/
capture erase "temp/basic.csv"
capture erase "temp/basic.tsv"
capture erase "temp/missing.csv"
capture erase "temp/quoted.csv"
capture erase "temp/noheader.csv"
capture erase "temp/numeric.csv"
capture erase "temp/large.csv"
capture erase "temp/uppercase.csv"
capture erase "temp/auto_export.csv"
capture erase "temp/census_export.csv"
capture erase "temp/bindquotes_test.csv"

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
print_summary "cimport"

* Return error code if any tests failed
if $TESTS_FAILED > 0 {
    exit 1
}
