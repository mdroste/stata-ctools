/*******************************************************************************
 * validate_cimport.do
 *
 * Comprehensive validation tests for cimport vs import delimited
 * Tests CSV import functionality across various scenarios
 ******************************************************************************/

do "validation/validate_setup.do"

di as text ""
di as text "======================================================================"
di as text "              CIMPORT VALIDATION TEST SUITE"
di as text "======================================================================"

* Create temp directory for test files
capture mkdir "validation/temp"

/*******************************************************************************
 * Create test CSV files
 ******************************************************************************/
print_section "Creating test files"

* Test file 1: Basic CSV with mixed types
capture file close _all
file open csvfile using "validation/temp/basic.csv", write replace
file write csvfile "id,name,value,price" _n
file write csvfile "1,Alpha,100,10.5" _n
file write csvfile "2,Beta,200,20.75" _n
file write csvfile "3,Gamma,300,30.125" _n
file write csvfile "4,Delta,400,40.0" _n
file write csvfile "5,Epsilon,500,50.99" _n
file close csvfile

* Test file 2: Tab-delimited
file open tsvfile using "validation/temp/basic.tsv", write replace
file write tsvfile "id" _tab "name" _tab "value" _n
file write tsvfile "1" _tab "Apple" _tab "1.5" _n
file write tsvfile "2" _tab "Banana" _tab "2.5" _n
file write tsvfile "3" _tab "Cherry" _tab "3.5" _n
file close tsvfile

* Test file 3: CSV with missing values
file open csvfile using "validation/temp/missing.csv", write replace
file write csvfile "id,name,value,score" _n
file write csvfile "1,A,100,95" _n
file write csvfile "2,B,,85" _n
file write csvfile "3,,300," _n
file write csvfile "4,D,400,75" _n
file write csvfile "5,E,," _n
file close csvfile

* Test file 4: CSV with quoted fields
file open csvfile using "validation/temp/quoted.csv", write replace
file write csvfile `"id,name,description"' _n
file write csvfile `"1,"Smith, John","A person with comma""' _n
file write csvfile `"2,"Doe, Jane","Another person""' _n
file write csvfile `"3,Simple,No quotes needed"' _n
file close csvfile

* Test file 5: CSV without header
file open csvfile using "validation/temp/noheader.csv", write replace
file write csvfile "1,Alpha,100" _n
file write csvfile "2,Beta,200" _n
file write csvfile "3,Gamma,300" _n
file close csvfile

* Test file 6: CSV with various numeric formats
file open csvfile using "validation/temp/numeric.csv", write replace
file write csvfile "integer,decimal,scientific,negative" _n
file write csvfile "1,1.5,1e10,-100" _n
file write csvfile "100,0.001,2.5e-5,-0.5" _n
file write csvfile "999999,123.456789,1e-10,-999999" _n
file close csvfile

* Test file 7: Large CSV
file open csvfile using "validation/temp/large.csv", write replace
file write csvfile "id,group,value1,value2,label" _n
forvalues i = 1/5000 {
    local g = mod(`i', 100) + 1
    local v1 = `i' * 1.5
    local v2 = runiform() * 1000
    file write csvfile "`i',`g',`v1',`v2',item`i'" _n
}
file close csvfile

* Test file 8: CSV with UPPER case headers
file open csvfile using "validation/temp/uppercase.csv", write replace
file write csvfile "ID,NAME,VALUE" _n
file write csvfile "1,Test,100" _n
file write csvfile "2,Data,200" _n
file close csvfile

di as result "[DONE] Test files created"

/*******************************************************************************
 * Basic import tests
 ******************************************************************************/
print_section "Basic import tests"

benchmark_import using "validation/temp/basic.csv", testname("basic CSV")

benchmark_import using "validation/temp/basic.tsv", delimiters(tab) testname("tab-delimited")

benchmark_import using "validation/temp/missing.csv", testname("with missing values")

benchmark_import using "validation/temp/quoted.csv", testname("quoted fields")

benchmark_import using "validation/temp/noheader.csv", varnames(nonames) testname("no header")

benchmark_import using "validation/temp/uppercase.csv", case(lower) testname("case(lower)")

benchmark_import using "validation/temp/numeric.csv", testname("numeric formats")

benchmark_import using "validation/temp/large.csv", testname("large file (5K rows)")

/*******************************************************************************
 * Round-trip tests (export then reimport)
 ******************************************************************************/
print_section "Round-trip tests"

sysuse auto, clear
export delimited using "validation/temp/auto_export.csv", replace
benchmark_import using "validation/temp/auto_export.csv", testname("auto round-trip")

sysuse census, clear
export delimited using "validation/temp/census_export.csv", replace
benchmark_import using "validation/temp/census_export.csv", testname("census round-trip")

/*******************************************************************************
 * Detailed comparison tests
 ******************************************************************************/
print_section "Detailed comparison tests"

* String variable types
import delimited using "validation/temp/basic.csv", clear
local stata_type : type name

cimport delimited using "validation/temp/basic.csv", clear
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
import delimited using "validation/temp/basic.csv", clear
local stata_id_type : type id
local stata_value_type : type value

cimport delimited using "validation/temp/basic.csv", clear
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
import delimited using "validation/temp/missing.csv", clear
quietly count if missing(value)
local stata_miss = r(N)

cimport delimited using "validation/temp/missing.csv", clear
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
import delimited using "validation/temp/large.csv", clear
quietly summarize value1
local stata_mean = r(mean)
local stata_n = _N

cimport delimited using "validation/temp/large.csv", clear
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
capture erase "validation/temp/basic.csv"
capture erase "validation/temp/basic.tsv"
capture erase "validation/temp/missing.csv"
capture erase "validation/temp/quoted.csv"
capture erase "validation/temp/noheader.csv"
capture erase "validation/temp/numeric.csv"
capture erase "validation/temp/large.csv"
capture erase "validation/temp/uppercase.csv"
capture erase "validation/temp/auto_export.csv"
capture erase "validation/temp/census_export.csv"

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
print_summary "cimport"

* Return error code if any tests failed
if $TESTS_FAILED > 0 {
    exit 1
}
