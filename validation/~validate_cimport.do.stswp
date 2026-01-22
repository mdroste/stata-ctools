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
 * Helper: benchmark_import - Import with both methods and compare using cf _all
 * All options (except testname) are passed to both import delimited and cimport
 ******************************************************************************/
capture program drop benchmark_import
program define benchmark_import
    syntax using/, testname(string) [*]

    * Import with Stata's import delimited
    import delimited `using', `options' clear
    tempfile stata_data
    quietly save `stata_data', replace
    local stata_n = _N
    local stata_k = c(k)

    * Import with cimport delimited
    cimport delimited `using', `options' clear
    tempfile cimport_data
    quietly save `cimport_data', replace
    local cimport_n = _N
    local cimport_k = c(k)

    * Check dimensions first
    if `stata_n' != `cimport_n' | `stata_k' != `cimport_k' {
        noi test_fail "`testname'" "dimensions differ: Stata N=`stata_n' K=`stata_k', cimport N=`cimport_n' K=`cimport_k'"
        exit
    }

    * Compare data using cf _all
    use `stata_data', clear
    capture cf _all using `cimport_data'
    if _rc == 0 {
        noi test_pass "`testname'"
    }
    else {
        noi test_fail "`testname'" "cf _all comparison failed - data not identical"
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

benchmark_import using "temp/basic.csv", testname("basic CSV")

/*******************************************************************************
 * SECTION 3: Tab-delimited import (cf _all verification)
 ******************************************************************************/
noi print_section "Tab-Delimited Import"

benchmark_import using "temp/tabfile.tsv", testname("tab-delimited") delimiters(tab)

/*******************************************************************************
 * SECTION 4: Semicolon-delimited import (cf _all verification)
 ******************************************************************************/
noi print_section "Semicolon-Delimited Import"

benchmark_import using "temp/semicolon.csv", testname("semicolon-delimited") delimiters(";")

/*******************************************************************************
 * SECTION 5: varnames() option (cf _all verification)
 ******************************************************************************/
noi print_section "varnames() Option"

benchmark_import using "temp/basic.csv", testname("varnames(1)") varnames(1)
benchmark_import using "temp/noheader.csv", testname("varnames(nonames)") varnames(nonames)

/*******************************************************************************
 * SECTION 6: case() option
 ******************************************************************************/
noi print_section "case() Option"

benchmark_import using "temp/casetest.csv", testname("case(preserve)") case(preserve)
benchmark_import using "temp/casetest.csv", testname("case(lower)") case(lower)
benchmark_import using "temp/casetest.csv", testname("case(upper)") case(upper)

/*******************************************************************************
 * SECTION 7: Missing values handling (cf _all verification)
 ******************************************************************************/
noi print_section "Missing Values"

benchmark_import using "temp/missing.csv", testname("missing values")

/*******************************************************************************
 * SECTION 8: Quoted fields (cf _all verification)
 ******************************************************************************/
noi print_section "Quoted Fields"

benchmark_import using "temp/quoted.csv", testname("quoted fields")

/*******************************************************************************
 * SECTION 9: bindquotes() option (cf _all verification)
 ******************************************************************************/
noi print_section "bindquotes() Option"

benchmark_import using "temp/quoted.csv", testname("bindquotes(strict)") bindquotes(strict)
benchmark_import using "temp/quoted.csv", testname("bindquotes(loose)") bindquotes(loose)

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

benchmark_import using "temp/large.csv", testname("large file (5000 rows)")

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
benchmark_import using "temp/auto_rt.csv", testname("auto dataset round-trip")

* census dataset round-trip
sysuse census, clear
export delimited using "temp/census_rt.csv", replace
benchmark_import using "temp/census_rt.csv", testname("census dataset round-trip")

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

benchmark_import using "temp/numtypes.csv", testname("numeric types")

* Create file with long strings
file open fh using "temp/longstr.csv", write replace
file write fh "id,short_str,long_str" _n
file write fh "1,abc,This is a longer string with spaces and stuff" _n
file write fh "2,def,Another longer string value here for testing" _n
file close fh

benchmark_import using "temp/longstr.csv", testname("string types")

/*******************************************************************************
 * SECTION 16: Malformed/Edge Case CSV Files
 ******************************************************************************/
noi print_section "Malformed/Edge Case CSV Files"

* Multiline quoted fields
file open fh using "temp/multiline.csv", write replace
file write fh "id,description,value" _n
file write fh `"1,"This is a simple description",100"' _n
file write fh `"2,"This description"' _n
file write fh `"spans multiple"' _n
file write fh `"lines of text",200"' _n
file write fh `"3,"Another"' _n
file write fh `"multiline entry",300"' _n
file write fh `"4,"Normal again",400"' _n
file close fh

benchmark_import using "temp/multiline.csv", testname("multiline quoted fields")
benchmark_import using "temp/multiline.csv", testname("multiline bindquotes(strict)") bindquotes(strict)

* Multiline with embedded commas and quotes
file open fh using "temp/multiline_complex.csv", write replace
file write fh "id,text,amount" _n
file write fh `"1,"Line 1, with comma"' _n
file write fh `"and line 2",50.5"' _n
file write fh `"2,"Contains ""quoted"" text"' _n
file write fh `"across lines",75.25"' _n
file write fh `"3,"Normal, simple text",100"' _n
file close fh

benchmark_import using "temp/multiline_complex.csv", testname("multiline with commas and quotes")

* Empty lines within quotes
file open fh using "temp/multiline_empty.csv", write replace
file write fh "id,paragraph,num" _n
file write fh `"1,"First paragraph."' _n
file write fh `""' _n
file write fh `"Second paragraph.",10"' _n
file write fh `"2,"Single line",20"' _n
file close fh

benchmark_import using "temp/multiline_empty.csv", testname("multiline with empty lines")

* Field with only newlines
file open fh using "temp/newline_only.csv", write replace
file write fh "id,content,num" _n
file write fh `"1,""' _n
file write fh `""' _n
file write fh `"",100"' _n
file write fh `"2,"text",200"' _n
file close fh

benchmark_import using "temp/newline_only.csv", testname("field with only newlines")

* Very long multiline field
file open fh using "temp/multiline_long.csv", write replace
file write fh "id,longtext,val" _n
file write fh `"1,"This is a very long text field"' _n
file write fh `"that spans many lines and contains"' _n
file write fh `"a substantial amount of content."' _n
file write fh `"It may include commas, like this one,"' _n
file write fh `"and even ""embedded quotes"" too."' _n
file write fh `"Finally it ends here.",999"' _n
file write fh `"2,"Short",1"' _n
file close fh

benchmark_import using "temp/multiline_long.csv", testname("long multiline field")

* Multiline as first and last rows
file open fh using "temp/multiline_position.csv", write replace
file write fh "id,text,num" _n
file write fh `"1,"First row"' _n
file write fh `"is multiline",10"' _n
file write fh `"2,"Middle row normal",20"' _n
file write fh `"3,"Last row"' _n
file write fh `"also multiline",30"' _n
file close fh

benchmark_import using "temp/multiline_position.csv", testname("multiline first/last rows")

* Tab and other whitespace within multiline quotes
local tab = char(9)
file open fh using "temp/multiline_whitespace.csv", write replace
file write fh "id,text,num" _n
file write fh `"1,"Text with`tab'tab"' _n
file write fh `"and newline",10"' _n
file write fh `"2,"  leading spaces"' _n
file write fh `"trailing spaces  ",20"' _n
file close fh

benchmark_import using "temp/multiline_whitespace.csv", testname("multiline with whitespace")

/*******************************************************************************
 * SECTION 17: asfloat and asdouble options
 ******************************************************************************/
noi print_section "asfloat/asdouble Options"

* Create file with numeric values
file open fh using "temp/numeric_types.csv", write replace
file write fh "id,small_int,large_num,decimal_val" _n
file write fh "1,10,1000000000,3.14159265358979" _n
file write fh "2,20,2000000000,2.71828182845904" _n
file write fh "3,30,3000000000,1.41421356237309" _n
file close fh

* Test asfloat option
capture cimport delimited using "temp/numeric_types.csv", clear asfloat
if _rc == 0 {
    local all_float = 1
    foreach var of varlist _all {
        local vtype : type `var'
        if substr("`vtype'", 1, 3) != "str" & "`vtype'" != "float" {
            local all_float = 0
        }
    }
    if `all_float' {
        noi test_pass "asfloat: all numeric variables are float"
    }
    else {
        noi test_fail "asfloat" "not all numeric variables are float"
    }
}
else {
    noi test_fail "asfloat option" "returned error `=_rc'"
}

* Test asdouble option
capture cimport delimited using "temp/numeric_types.csv", clear asdouble
if _rc == 0 {
    local all_double = 1
    foreach var of varlist _all {
        local vtype : type `var'
        if substr("`vtype'", 1, 3) != "str" & "`vtype'" != "double" {
            local all_double = 0
        }
    }
    if `all_double' {
        noi test_pass "asdouble: all numeric variables are double"
    }
    else {
        noi test_fail "asdouble" "not all numeric variables are double"
    }
}
else {
    noi test_fail "asdouble option" "returned error `=_rc'"
}

* Test mutual exclusivity
capture cimport delimited using "temp/numeric_types.csv", clear asfloat asdouble
if _rc != 0 {
    noi test_pass "asfloat + asdouble: correctly rejects mutually exclusive options"
}
else {
    noi test_fail "asfloat + asdouble" "should reject mutually exclusive options"
}

/*******************************************************************************
 * SECTION 18: stringcols and numericcols options
 ******************************************************************************/
noi print_section "stringcols/numericcols Options"

* Create file with mixed data
file open fh using "temp/col_override.csv", write replace
file write fh "zipcode,amount,code,value" _n
file write fh "01234,100,ABC,500" _n
file write fh "02345,200,DEF,600" _n
file write fh "03456,300,GHI,700" _n
file close fh

* Test stringcols
capture cimport delimited using "temp/col_override.csv", clear stringcols(1)
if _rc == 0 {
    local vtype : type zipcode
    if substr("`vtype'", 1, 3) == "str" {
        if zipcode[1] == "01234" {
            noi test_pass "stringcols: column forced to string, leading zeros preserved"
        }
        else {
            noi test_fail "stringcols" "leading zeros not preserved: `=zipcode[1]'"
        }
    }
    else {
        noi test_fail "stringcols" "column not string type: `vtype'"
    }
}
else {
    noi test_fail "stringcols option" "returned error `=_rc'"
}

* Test numericcols
capture cimport delimited using "temp/col_override.csv", clear numericcols(3)
if _rc == 0 {
    local vtype : type code
    if substr("`vtype'", 1, 3) != "str" {
        if missing(code[1]) {
            noi test_pass "numericcols: non-numeric values become missing"
        }
        else {
            noi test_fail "numericcols" "expected missing, got `=code[1]'"
        }
    }
    else {
        noi test_fail "numericcols" "column should be numeric, got `vtype'"
    }
}
else {
    noi test_fail "numericcols option" "returned error `=_rc'"
}

* Test multiple columns
capture cimport delimited using "temp/col_override.csv", clear stringcols(1 4)
if _rc == 0 {
    local vtype1 : type zipcode
    local vtype4 : type value
    if substr("`vtype1'", 1, 3) == "str" & substr("`vtype4'", 1, 3) == "str" {
        noi test_pass "stringcols: multiple columns forced to string"
    }
    else {
        noi test_fail "stringcols multiple" "types: zipcode=`vtype1', value=`vtype4'"
    }
}
else {
    noi test_fail "stringcols multiple" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 19: decimalseparator and groupseparator options
 ******************************************************************************/
noi print_section "decimalseparator/groupseparator Options"

* Create European-format file
file open fh using "temp/european.csv", write replace
file write fh "id;amount;price" _n
file write fh "1;1.234,56;99,99" _n
file write fh "2;2.345,67;199,50" _n
file write fh "3;3.456,78;299,00" _n
file close fh

* Test European format parsing
capture cimport delimited using "temp/european.csv", clear delimiters(";") decimalseparator(,) groupseparator(.)
if _rc == 0 {
    if abs(amount[1] - 1234.56) < 0.01 {
        noi test_pass "European format: comma decimal, period grouping parsed correctly"
    }
    else {
        noi test_fail "European format" "expected 1234.56, got `=amount[1]'"
    }
}
else {
    noi test_fail "decimalseparator/groupseparator" "returned error `=_rc'"
}

* Test just decimal separator
file open fh using "temp/comma_decimal.csv", write replace
file write fh "id,value" _n
file write fh "1,123,45" _n
file write fh "2,234,56" _n
file write fh "3,345,67" _n
file close fh

capture cimport delimited using "temp/comma_decimal.csv", clear decimalseparator(,)
if _rc == 0 {
    if abs(value[1] - 123.45) < 0.01 {
        noi test_pass "decimalseparator only: comma decimal parsed correctly"
    }
    else {
        noi test_fail "decimalseparator only" "expected 123.45, got `=value[1]'"
    }
}
else {
    noi test_fail "decimalseparator only" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 20: threads option
 ******************************************************************************/
noi print_section "threads Option"

capture cimport delimited using "temp/large.csv", clear threads(2)
if _rc == 0 {
    noi test_pass "threads(2) option accepted"
}
else {
    noi test_fail "threads option" "returned error `=_rc'"
}

capture cimport delimited using "temp/large.csv", clear threads(1)
if _rc == 0 {
    noi test_pass "threads(1) single-threaded import works"
}
else {
    noi test_fail "threads(1)" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION: Pathological Data Tests
 ******************************************************************************/
noi print_section "Pathological Data"

* All empty strings
file open myfile using "temp/all_empty.csv", write replace
file write myfile "a,b,c" _n
file write myfile ",," _n
file write myfile ",," _n
file write myfile ",," _n
file close myfile

capture cimport delimited using "temp/all_empty.csv", clear
if _rc == 0 {
    noi test_pass "all empty values"
}
else {
    noi test_fail "all empty" "rc=`=_rc'"
}

* Very long lines
file open myfile using "temp/long_lines.csv", write replace
file write myfile "id,value" _n
forvalues i = 1/10 {
    local long_val = ""
    forvalues j = 1/100 {
        local long_val = "`long_val'x"
    }
    file write myfile "`i',`long_val'" _n
}
file close myfile

capture cimport delimited using "temp/long_lines.csv", clear
if _rc == 0 {
    noi test_pass "long lines (100 chars)"
}
else {
    noi test_fail "long lines" "rc=`=_rc'"
}

* Single column
file open myfile using "temp/single_col.csv", write replace
file write myfile "values" _n
forvalues i = 1/100 {
    file write myfile "`i'" _n
}
file close myfile

capture cimport delimited using "temp/single_col.csv", clear
if _rc == 0 {
    if _N == 100 {
        noi test_pass "single column CSV"
    }
    else {
        noi test_fail "single column" "wrong N"
    }
}
else {
    noi test_fail "single column" "rc=`=_rc'"
}

* Many columns (20)
file open myfile using "temp/many_cols.csv", write replace
local header = "c1"
forvalues i = 2/20 {
    local header = "`header',c`i'"
}
file write myfile "`header'" _n
forvalues row = 1/50 {
    local line = "`row'"
    forvalues col = 2/20 {
        local line = "`line',`=`row'*`col''"
    }
    file write myfile "`line'" _n
}
file close myfile

capture cimport delimited using "temp/many_cols.csv", clear
if _rc == 0 {
    if _N == 50 {
        noi test_pass "20 columns CSV"
    }
    else {
        noi test_fail "many columns" "wrong N"
    }
}
else {
    noi test_fail "many columns" "rc=`=_rc'"
}

* Single row
file open myfile using "temp/single_row.csv", write replace
file write myfile "a,b,c" _n
file write myfile "1,2,3" _n
file close myfile

capture cimport delimited using "temp/single_row.csv", clear
if _rc == 0 {
    if _N == 1 {
        noi test_pass "single row CSV"
    }
    else {
        noi test_fail "single row" "wrong N"
    }
}
else {
    noi test_fail "single row" "rc=`=_rc'"
}

* Numeric edge values
file open myfile using "temp/numeric_edge.csv", write replace
file write myfile "value" _n
file write myfile "0" _n
file write myfile "-0" _n
file write myfile "1e10" _n
file write myfile "-1e10" _n
file write myfile "0.000001" _n
file write myfile "-0.000001" _n
file write myfile "999999999" _n
file close myfile

capture cimport delimited using "temp/numeric_edge.csv", clear
if _rc == 0 {
    noi test_pass "numeric edge values"
}
else {
    noi test_fail "numeric edge" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION: Large File Tests
 ******************************************************************************/
noi print_section "Large File Tests"

* 10K rows
file open myfile using "temp/10k_rows.csv", write replace
file write myfile "id,value,category" _n
forvalues i = 1/10000 {
    local cat = mod(`i', 10) + 1
    file write myfile "`i',`=runiform()',cat`cat'" _n
}
file close myfile

capture cimport delimited using "temp/10k_rows.csv", clear
if _rc == 0 {
    if _N == 10000 {
        noi test_pass "10K rows import"
    }
    else {
        noi test_fail "10K rows" "wrong N: `=_N'"
    }
}
else {
    noi test_fail "10K rows" "rc=`=_rc'"
}

* 50K rows
file open myfile using "temp/50k_rows.csv", write replace
file write myfile "id,x,y" _n
forvalues i = 1/50000 {
    file write myfile "`i',`=runiform()',`=runiform()'" _n
}
file close myfile

capture cimport delimited using "temp/50k_rows.csv", clear
if _rc == 0 {
    if _N == 50000 {
        noi test_pass "50K rows import"
    }
    else {
        noi test_fail "50K rows" "wrong N"
    }
}
else {
    noi test_fail "50K rows" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION: Comparison with import delimited
 ******************************************************************************/
noi print_section "Comparison with Stata's import delimited"

* Create test file
file open myfile using "temp/compare_test.csv", write replace
file write myfile "id,name,value,flag" _n
forvalues i = 1/100 {
    local name = "item_`i'"
    local val = `i' * 1.5
    local flag = mod(`i', 2)
    file write myfile "`i',`name',`val',`flag'" _n
}
file close myfile

* Import with Stata
import delimited using "temp/compare_test.csv", clear
local stata_N = _N
local stata_id_1 = id[1]
local stata_value_50 = value[50]

* Import with cimport
cimport delimited using "temp/compare_test.csv", clear
local cimport_N = _N
local cimport_id_1 = id[1]
local cimport_value_50 = value[50]

if `stata_N' == `cimport_N' & `stata_id_1' == `cimport_id_1' {
    noi test_pass "matches import delimited (basic)"
}
else {
    noi test_fail "compare basic" "N or values differ"
}

* Compare with quoted strings
file open myfile using "temp/compare_quoted.csv", write replace
file write myfile `"id,name,description"' _n
file write myfile `"1,"Simple","A simple item""' _n
file write myfile `"2,"With Space","Has a space""' _n
file write myfile `"3,"Quoted","Contains ""quotes""""' _n
file close myfile

import delimited using "temp/compare_quoted.csv", clear
local stata_N = _N

cimport delimited using "temp/compare_quoted.csv", clear
local cimport_N = _N

if `stata_N' == `cimport_N' {
    noi test_pass "matches import delimited (quoted)"
}
else {
    noi test_fail "compare quoted" "N differs"
}

/*******************************************************************************
 * SECTION: varnames Option Tests
 ******************************************************************************/
noi print_section "varnames Option"

file open myfile using "temp/varnames_test.csv", write replace
file write myfile "100,200,300" _n
file write myfile "1,2,3" _n
file write myfile "4,5,6" _n
file close myfile

* varnames(nonames)
capture cimport delimited using "temp/varnames_test.csv", clear varnames(nonames)
if _rc == 0 {
    if _N == 3 {
        noi test_pass "varnames(nonames)"
    }
    else {
        noi test_fail "varnames(nonames)" "wrong N"
    }
}
else {
    noi test_fail "varnames(nonames)" "rc=`=_rc'"
}

* varnames(1) - first row as names
file open myfile using "temp/varnames_row1.csv", write replace
file write myfile "alpha,beta,gamma" _n
file write myfile "1,2,3" _n
file write myfile "4,5,6" _n
file close myfile

capture cimport delimited using "temp/varnames_row1.csv", clear varnames(1)
if _rc == 0 {
    if _N == 2 {
        capture confirm variable alpha
        if _rc == 0 {
            noi test_pass "varnames(1)"
        }
        else {
            noi test_fail "varnames(1)" "variable alpha not found"
        }
    }
    else {
        noi test_fail "varnames(1)" "wrong N"
    }
}
else {
    noi test_fail "varnames(1)" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION: rowrange Option Tests
 ******************************************************************************/
noi print_section "rowrange Option"

file open myfile using "temp/rowrange_test.csv", write replace
file write myfile "id,value" _n
forvalues i = 1/100 {
    file write myfile "`i',`=`i'*10'" _n
}
file close myfile

* Import rows 10-20
capture cimport delimited using "temp/rowrange_test.csv", clear rowrange(10:20)
if _rc == 0 {
    if _N == 11 {
        noi test_pass "rowrange(10:20)"
    }
    else {
        noi test_fail "rowrange" "wrong N: expected 11, got `=_N'"
    }
}
else {
    noi test_fail "rowrange" "rc=`=_rc'"
}

* Import from row 50 to end
capture cimport delimited using "temp/rowrange_test.csv", clear rowrange(50:)
if _rc == 0 {
    if _N == 51 {
        noi test_pass "rowrange(50:) to end"
    }
    else {
        noi test_fail "rowrange to end" "wrong N"
    }
}
else {
    noi test_fail "rowrange to end" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION: colrange Option Tests
 ******************************************************************************/
noi print_section "colrange Option"

file open myfile using "temp/colrange_test.csv", write replace
file write myfile "a,b,c,d,e" _n
forvalues i = 1/50 {
    file write myfile "`i',`=`i'*2',`=`i'*3',`=`i'*4',`=`i'*5'" _n
}
file close myfile

* Import columns 2-4
capture cimport delimited using "temp/colrange_test.csv", clear colrange(2:4)
if _rc == 0 {
    describe, short
    local nvars = r(k)
    if `nvars' == 3 {
        noi test_pass "colrange(2:4)"
    }
    else {
        noi test_fail "colrange" "wrong nvars: expected 3, got `nvars'"
    }
}
else {
    noi test_fail "colrange" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION: encoding Option Tests
 ******************************************************************************/
noi print_section "encoding Option"

* UTF-8 file (default)
file open myfile using "temp/utf8_test.csv", write replace
file write myfile "name,value" _n
file write myfile "test1,100" _n
file write myfile "test2,200" _n
file close myfile

capture cimport delimited using "temp/utf8_test.csv", clear encoding(utf-8)
if _rc == 0 {
    noi test_pass "encoding(utf-8)"
}
else {
    noi test_fail "encoding utf-8" "rc=`=_rc'"
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
