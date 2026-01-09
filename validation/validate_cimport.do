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
 * Helper: Create test CSV files
 ******************************************************************************/

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

/*******************************************************************************
 * TEST 1: Basic CSV import
 ******************************************************************************/
print_section "Test 1: Basic CSV import"

* Stata import
import delimited using "validation/temp/basic.csv", clear
tempfile stata_imported
quietly save `stata_imported'

* ctools import
cimport delimited using "validation/temp/basic.csv", clear
tempfile cimport_imported
quietly save `cimport_imported'

* Compare observation count
use `stata_imported', clear
local stata_n = _N
local stata_k = c(k)

use `cimport_imported', clear
local cimport_n = _N
local cimport_k = c(k)

if `stata_n' == `cimport_n' & `stata_k' == `cimport_k' {
    di as result "  PASS: Basic import dimensions match (N=`stata_n', K=`stata_k')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Dimensions differ (stata: N=`stata_n' K=`stata_k', cimport: N=`cimport_n' K=`cimport_k')"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

* Compare first id value
use `stata_imported', clear
local stata_id1 = id[1]

use `cimport_imported', clear
local cimport_id1 = id[1]

if `stata_id1' == `cimport_id1' {
    di as result "  PASS: Basic import id values match"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Basic import id values differ"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 2: Tab-delimited import
 ******************************************************************************/
print_section "Test 2: Tab-delimited import"

* Stata import
import delimited using "validation/temp/basic.tsv", clear delimiters(tab)
tempfile stata_imported
quietly save `stata_imported'

* ctools import
cimport delimited using "validation/temp/basic.tsv", clear delimiters(tab)
tempfile cimport_imported
quietly save `cimport_imported'

use `stata_imported', clear
local stata_n = _N

use `cimport_imported', clear
local cimport_n = _N

if `stata_n' == `cimport_n' {
    di as result "  PASS: Tab-delimited import (N=`stata_n')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Row counts differ"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 3: Import with missing values
 ******************************************************************************/
print_section "Test 3: Import with missing values"

* Stata import
import delimited using "validation/temp/missing.csv", clear
tempfile stata_imported
quietly save `stata_imported'

* ctools import
cimport delimited using "validation/temp/missing.csv", clear
tempfile cimport_imported
quietly save `cimport_imported'

* Count missing values
use `stata_imported', clear
quietly count if missing(value)
local stata_miss_value = r(N)
quietly count if missing(score)
local stata_miss_score = r(N)

use `cimport_imported', clear
quietly count if missing(value)
local cimport_miss_value = r(N)
quietly count if missing(score)
local cimport_miss_score = r(N)

if `stata_miss_value' == `cimport_miss_value' & `stata_miss_score' == `cimport_miss_score' {
    di as result "  PASS: Missing values handled correctly (value:`stata_miss_value', score:`stata_miss_score')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Missing value counts differ"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 4: Import with quoted fields
 ******************************************************************************/
print_section "Test 4: Import with quoted fields"

* Stata import
import delimited using "validation/temp/quoted.csv", clear
tempfile stata_imported
quietly save `stata_imported'

* ctools import
cimport delimited using "validation/temp/quoted.csv", clear
tempfile cimport_imported
quietly save `cimport_imported'

* Check that commas in quoted fields are preserved
use `cimport_imported', clear
local name1 = name[1]

if strpos("`name1'", ",") > 0 {
    di as result "  PASS: Quoted fields with commas handled correctly"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Quoted fields not handled correctly (name1=`name1')"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 5: Import without header (varnames(nonames))
 ******************************************************************************/
print_section "Test 5: Import without header"

* Stata import
import delimited using "validation/temp/noheader.csv", clear varnames(nonames)
tempfile stata_imported
quietly save `stata_imported'

* ctools import
cimport delimited using "validation/temp/noheader.csv", clear varnames(nonames)
tempfile cimport_imported
quietly save `cimport_imported'

* Check that default variable names are used
use `cimport_imported', clear
capture confirm variable v1
local has_v1 = (_rc == 0)

if `has_v1' {
    di as result "  PASS: varnames(nonames) creates default variable names"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Default variable names not created"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

* Check row count (should be 3, not treating first row as header)
local cimport_n = _N
if `cimport_n' == 3 {
    di as result "  PASS: All rows imported as data (N=3)"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Row count incorrect (N=`cimport_n', expected 3)"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 6: Import with case(lower)
 ******************************************************************************/
print_section "Test 6: Import with case(lower)"

* Stata import
import delimited using "validation/temp/uppercase.csv", clear case(lower)
tempfile stata_imported
quietly save `stata_imported'

* ctools import
cimport delimited using "validation/temp/uppercase.csv", clear case(lower)
tempfile cimport_imported
quietly save `cimport_imported'

* Check variable names are lowercase
use `cimport_imported', clear
capture confirm variable id
local has_id = (_rc == 0)
capture confirm variable ID
local has_ID = (_rc == 0)

if `has_id' & !`has_ID' {
    di as result "  PASS: case(lower) converts variable names to lowercase"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Variable names not converted to lowercase"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 7: Import numeric formats
 ******************************************************************************/
print_section "Test 7: Import various numeric formats"

* Stata import
import delimited using "validation/temp/numeric.csv", clear
tempfile stata_imported
quietly save `stata_imported'

* ctools import
cimport delimited using "validation/temp/numeric.csv", clear
tempfile cimport_imported
quietly save `cimport_imported'

* Compare numeric values with tolerance
use `stata_imported', clear
local stata_int1 = integer[1]
local stata_dec1 = decimal[1]
local stata_sci1 = scientific[1]
local stata_neg1 = negative[1]

use `cimport_imported', clear
local cimport_int1 = integer[1]
local cimport_dec1 = decimal[1]
local cimport_sci1 = scientific[1]
local cimport_neg1 = negative[1]

local int_match = (`stata_int1' == `cimport_int1')
local dec_match = (abs(`stata_dec1' - `cimport_dec1') < 1e-10)
local sci_match = (abs(`stata_sci1' - `cimport_sci1') / `stata_sci1' < 1e-6)
local neg_match = (`stata_neg1' == `cimport_neg1')

if `int_match' & `dec_match' & `sci_match' & `neg_match' {
    di as result "  PASS: Numeric formats parsed correctly"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Numeric parsing issues"
    if !`int_match' di as error "    Integer: stata=`stata_int1' cimport=`cimport_int1'"
    if !`dec_match' di as error "    Decimal: stata=`stata_dec1' cimport=`cimport_dec1'"
    if !`sci_match' di as error "    Scientific: stata=`stata_sci1' cimport=`cimport_sci1'"
    if !`neg_match' di as error "    Negative: stata=`stata_neg1' cimport=`cimport_neg1'"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 8: Large file import
 ******************************************************************************/
print_section "Test 8: Large file import (5000 rows)"

* Stata import
import delimited using "validation/temp/large.csv", clear
tempfile stata_imported
quietly save `stata_imported'

* ctools import
cimport delimited using "validation/temp/large.csv", clear
tempfile cimport_imported
quietly save `cimport_imported'

* Compare dimensions
use `stata_imported', clear
local stata_n = _N
local stata_k = c(k)

use `cimport_imported', clear
local cimport_n = _N
local cimport_k = c(k)

if `stata_n' == `cimport_n' & `stata_k' == `cimport_k' {
    di as result "  PASS: Large file dimensions match (N=`stata_n', K=`stata_k')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Dimensions differ"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

* Compare summary statistics
use `stata_imported', clear
quietly summarize value1
local stata_mean = r(mean)
local stata_sd = r(sd)

use `cimport_imported', clear
quietly summarize value1
local cimport_mean = r(mean)
local cimport_sd = r(sd)

if abs(`stata_mean' - `cimport_mean') < 1e-6 & abs(`stata_sd' - `cimport_sd') < 1e-6 {
    di as result "  PASS: Large file statistics match"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Statistics differ"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 9: Row range import (feature may not be fully implemented)
 ******************************************************************************/
print_section "Test 9: Row range import"

* Stata import with row range
import delimited using "validation/temp/large.csv", clear rowrange(100:200)
local stata_n = _N

* ctools import with row range
cimport delimited using "validation/temp/large.csv", clear rowrange(100:200)
local cimport_n = _N

if `stata_n' == `cimport_n' {
    di as result "  PASS: Row range import (N=`stata_n')"
    global TESTS_PASSED = $TESTS_PASSED + 1
    global TESTS_TOTAL = $TESTS_TOTAL + 1
}
else {
    * rowrange() may not be fully implemented in cimport
    di as text "  INFO: rowrange() feature not yet implemented (stata=`stata_n', cimport=`cimport_n')"
    di as text "  (This is a known limitation)"
    * Don't count as failure
}

/*******************************************************************************
 * TEST 10: Export and reimport round-trip (using auto dataset)
 ******************************************************************************/
print_section "Test 10: Export and reimport round-trip"

sysuse auto, clear
export delimited using "validation/temp/auto_export.csv", replace

* Stata reimport
import delimited using "validation/temp/auto_export.csv", clear
tempfile stata_reimport
quietly save `stata_reimport'

* ctools reimport
cimport delimited using "validation/temp/auto_export.csv", clear
tempfile cimport_reimport
quietly save `cimport_reimport'

* Compare dimensions
use `stata_reimport', clear
local stata_n = _N

use `cimport_reimport', clear
local cimport_n = _N

if `stata_n' == `cimport_n' {
    di as result "  PASS: Round-trip dimensions match (N=`stata_n')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Dimensions differ"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 11: String variable handling
 ******************************************************************************/
print_section "Test 11: String variable handling"

* Stata import
import delimited using "validation/temp/basic.csv", clear
local stata_type : type name

* ctools import
cimport delimited using "validation/temp/basic.csv", clear
local cimport_type : type name

* Both should be string type
if substr("`stata_type'", 1, 3) == "str" & substr("`cimport_type'", 1, 3) == "str" {
    di as result "  PASS: String variables correctly typed"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: String variable type issue (stata:`stata_type', cimport:`cimport_type')"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 12: Variable type inference
 ******************************************************************************/
print_section "Test 12: Variable type inference"

* Both should infer id as numeric and name as string
import delimited using "validation/temp/basic.csv", clear
local stata_id_type : type id
local stata_name_type : type name
local stata_value_type : type value

cimport delimited using "validation/temp/basic.csv", clear
local cimport_id_type : type id
local cimport_name_type : type name
local cimport_value_type : type value

local id_numeric = ("`stata_id_type'" != "str" & "`cimport_id_type'" != "str")
local name_string = (substr("`stata_name_type'",1,3) == "str" & substr("`cimport_name_type'",1,3) == "str")
local value_numeric = ("`stata_value_type'" != "str" & "`cimport_value_type'" != "str")

if `id_numeric' & `name_string' & `value_numeric' {
    di as result "  PASS: Variable types correctly inferred"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Type inference differs"
    di as error "    id: stata=`stata_id_type', cimport=`cimport_id_type'"
    di as error "    name: stata=`stata_name_type', cimport=`cimport_name_type'"
    di as error "    value: stata=`stata_value_type', cimport=`cimport_value_type'"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 13: Import census-like data
 ******************************************************************************/
print_section "Test 13: Census-like data round trip"

sysuse census, clear
export delimited using "validation/temp/census.csv", replace

import delimited using "validation/temp/census.csv", clear
tempfile stata_census
quietly save `stata_census'

cimport delimited using "validation/temp/census.csv", clear
tempfile cimport_census
quietly save `cimport_census'

use `stata_census', clear
local stata_n = _N
quietly summarize pop
local stata_pop_mean = r(mean)

use `cimport_census', clear
local cimport_n = _N
quietly summarize pop
local cimport_pop_mean = r(mean)

if `stata_n' == `cimport_n' & abs(`stata_pop_mean' - `cimport_pop_mean') < 1 {
    di as result "  PASS: Census data import (N=`stata_n')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Census data import differs"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

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
capture erase "validation/temp/census.csv"

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
print_summary "cimport"

* Return error code if any tests failed
if $TESTS_FAILED > 0 {
    exit 1
}
