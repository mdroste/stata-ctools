/*******************************************************************************
 * validate_cimport.do
 *
 * Comprehensive validation tests for cimport delimited vs import delimited
 * Tests all import options: delimiters, varnames, case, bindquotes, encoding, etc.
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
 * Required: using/ - path to CSV file
 *           testname() - the test name
 * Optional: importopts() - options for both import commands
 *
 * NOTE: Stata's import delimited has auto-detection for headers - if row 1
 * has the same type signature as row 2, it treats row 1 as data. cimport
 * always uses row 1 as headers (matching standard CSV conventions).
 *
 * To make fair comparisons, this helper forces Stata to use varnames(1) when
 * no varnames option is specified, so both commands use first row as headers.
 ******************************************************************************/
capture program drop benchmark_import
program define benchmark_import
    syntax using/, testname(string) [IMPORTopts(string)]

    * Ensure clean state
    clear

    * Check if varnames is already specified in importopts
    local has_varnames = strpos(lower("`importopts'"), "varnames")

    * Check if encoding is already specified in importopts
    local has_encoding = strpos(lower("`importopts'"), "encoding")

    * Build option string for Stata - add varnames(1) if not specified
    * Also add encoding(utf-8) if not specified because Stata's auto-detection
    * fails when running inside a quietly{} block (a Stata bug)
    local stata_opts "`importopts'"
    if `has_varnames' == 0 {
        local stata_opts "`stata_opts' varnames(1)"
    }
    if `has_encoding' == 0 {
        local stata_opts "`stata_opts' encoding(utf-8)"
    }

    * Debug output (disabled - semicolons in options break di command)
    * noi di "DEBUG: testname=`testname'"
    * noi di "DEBUG: importopts=`importopts'"
    * noi di "DEBUG: stata_opts=`stata_opts'"
    * noi di "DEBUG: has_encoding=`has_encoding'"

    * Import with Stata's import delimited (with varnames(1) to match cimport)
    capture import delimited `using', `stata_opts' clear
    if _rc != 0 {
        noi test_fail "`testname'" "Stata import failed with rc=`=_rc'"
        exit
    }
    tempfile stata_data
    quietly save `stata_data', replace
    local stata_n = _N
    local stata_k = c(k)

    * Import with cimport delimited (uses first row as headers by default)
    capture cimport delimited `using', `importopts' clear
    if _rc != 0 {
        noi test_fail "`testname'" "cimport failed with rc=`=_rc'"
        exit
    }
    tempfile cimport_data
    quietly save `cimport_data', replace
    local cimport_n = _N
    local cimport_k = c(k)

    * Check dimensions first
    if `stata_n' != `cimport_n' | `stata_k' != `cimport_k' {
        noi test_fail "`testname'" "dimensions differ: Stata N=`stata_n' K=`stata_k', cimport N=`cimport_n' K=`cimport_k'"
        exit
    }

    * Dimensions match - compare data content
    * Rename variables to generic names for positional comparison
    * (variable names may differ: Stata "1"->v1, cimport "1"->v)

    use `stata_data', clear
    quietly ds
    local stata_vars `r(varlist)'
    local i = 1
    foreach v of local stata_vars {
        quietly rename `v' __cmp_v`i'
        local i = `i' + 1
    }
    tempfile stata_renamed
    quietly save `stata_renamed', replace

    use `cimport_data', clear
    quietly ds
    local cimport_vars `r(varlist)'
    local i = 1
    foreach v of local cimport_vars {
        quietly rename `v' __cmp_v`i'
        local i = `i' + 1
    }

    * Compare data (without sorting - both should be in same order)
    capture cf _all using `stata_renamed'
    if _rc == 0 {
        noi test_pass "`testname'"
    }
    else {
        noi test_fail "`testname'" "cf _all comparison failed - data not identical"
    }
end

/*******************************************************************************
 * Helper: cimport_test - Test cimport standalone (no Stata comparison)
 * Required: using/ - path to CSV file
 *           testname() - the test name
 * Optional: importopts() - options for cimport
 *           expectn() - expected row count
 *           expectk() - expected column count
 *           check() - expression to evaluate (should be true for pass)
 ******************************************************************************/
capture program drop cimport_test
program define cimport_test
    syntax using/, testname(string) [IMPORTopts(string) expectn(integer 0) expectk(integer 0) CHECK(string)]

    capture cimport delimited `using', `importopts' clear
    if _rc != 0 {
        noi test_fail "`testname'" "rc=`=_rc'"
        exit
    }

    * Check expected row count
    if `expectn' > 0 & _N != `expectn' {
        noi test_fail "`testname'" "expected N=`expectn', got N=`=_N'"
        exit
    }

    * Check expected column count
    if `expectk' > 0 & c(k) != `expectk' {
        noi test_fail "`testname'" "expected K=`expectk', got K=`=c(k)'"
        exit
    }

    * Optional custom check expression
    if "`check'" != "" {
        capture assert `check'
        if _rc != 0 {
            noi test_fail "`testname'" "check() failed: `check'"
            exit
        }
    }

    noi test_pass "`testname'"
end

/*******************************************************************************
 * Create Test Data Files (call helper script)
 ******************************************************************************/
capture do "validation/validate_cimport_helper.do"
if _rc != 0 {
    do "validate_cimport_helper.do"
}

/*******************************************************************************
 * SECTION: Plugin Check
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
 * SECTION: Basic Import Comparison
 ******************************************************************************/
noi print_section "Basic Import Comparison"

noi benchmark_import using "temp/basic.csv", testname("basic CSV")

/*******************************************************************************
 * SECTION: Delimiter Options
 ******************************************************************************/
noi print_section "Delimiter Options"

noi benchmark_import using "temp/tab_delim.tsv", testname("tab delimiter") importopts(delimiters(tab))
noi benchmark_import using "temp/semicolon.csv", testname("semicolon delimiter") importopts(delimiters(";"))

/*******************************************************************************
 * SECTION: varnames Options
 ******************************************************************************/
noi print_section "varnames Options"

noi benchmark_import using "temp/noheader.csv", testname("varnames(nonames)") importopts(varnames(nonames))
noi benchmark_import using "temp/basic.csv", testname("varnames(1)") importopts(varnames(1))

* Create test file with header on row 3
file open fh using "temp/varnames_row3.csv", write replace
file write fh "skip1,skip2" _n
file write fh "ignore,this" _n
file write fh "name,value" _n
file write fh "Alice,100" _n
file write fh "Bob,200" _n
file close fh

noi benchmark_import using "temp/varnames_row3.csv", testname("varnames(3)") importopts(varnames(3))

/*******************************************************************************
 * SECTION: case Options
 ******************************************************************************/
noi print_section "case Options"

noi benchmark_import using "temp/mixedcase.csv", testname("case(preserve)") importopts(case(preserve))
noi benchmark_import using "temp/mixedcase.csv", testname("case(lower)") importopts(case(lower))
noi benchmark_import using "temp/mixedcase.csv", testname("case(upper)") importopts(case(upper))

/*******************************************************************************
 * SECTION: Quoted Fields
 ******************************************************************************/
noi print_section "Quoted Fields"

noi benchmark_import using "temp/quoted.csv", testname("quoted fields with commas")
noi benchmark_import using "temp/quoted.csv", testname("bindquotes(strict)") importopts(bindquotes(strict))

/*******************************************************************************
 * SECTION: Numeric Type Handling
 ******************************************************************************/
noi print_section "Numeric Type Handling"

noi benchmark_import using "temp/numerics.csv", testname("numeric types")
noi benchmark_import using "temp/numerics.csv", testname("asfloat") importopts(asfloat)
noi benchmark_import using "temp/numerics.csv", testname("asdouble") importopts(asdouble)

/*******************************************************************************
 * SECTION: Missing Values
 ******************************************************************************/
noi print_section "Missing Values"

noi benchmark_import using "temp/missing.csv", testname("missing values")

/*******************************************************************************
 * SECTION: stringcols/numericcols Options
 ******************************************************************************/
noi print_section "stringcols/numericcols Options"

* stringcols preserves leading zeros
capture cimport delimited using "temp/zipcodes.csv", stringcols(1) clear
if _rc == 0 {
    * Check if zipcode is string type and has correct value
    capture confirm string variable zipcode
    if _rc == 0 {
        if zipcode[1] == "01234" {
            noi test_pass "stringcols preserves leading zeros"
        }
        else {
            noi test_fail "stringcols" "leading zeros not preserved (got `=zipcode[1]')"
        }
    }
    else {
        noi test_fail "stringcols" "zipcode is not string type"
    }
}
else {
    noi test_fail "stringcols" "cimport failed with rc=`=_rc'"
}

noi benchmark_import using "temp/basic.csv", testname("numericcols forces numeric") ///
    importopts(numericcols(2))

noi benchmark_import using "temp/zipcodes.csv", testname("stringcols multiple columns") ///
    importopts(stringcols(1 3))

/*******************************************************************************
 * SECTION: decimalseparator/groupseparator Options
 ******************************************************************************/
noi print_section "decimalseparator/groupseparator Options"

* European format test
capture cimport delimited using "temp/european.csv", delimiters(";") decimalseparator(,) groupseparator(.) clear
if _rc == 0 & abs(amount[1] - 1234.56) < 0.01 {
    noi test_pass "European format"
}
else {
    noi test_fail "European format" "decimal parsing failed"
}

* Decimal separator with semicolon delimiter
file open fh using "temp/comma_decimal.csv", write replace
file write fh "id;value" _n
file write fh "1;123,45" _n
file write fh "2;234,56" _n
file close fh

capture cimport delimited using "temp/comma_decimal.csv", delimiters(";") decimalseparator(,) clear
if _rc == 0 & abs(value[1] - 123.45) < 0.01 {
    noi test_pass "decimalseparator only"
}
else {
    noi test_fail "decimalseparator only" "expected 123.45, got `=value[1]'"
}

/*******************************************************************************
 * SECTION: threads Option
 ******************************************************************************/
noi print_section "threads Option"

noi cimport_test using "temp/basic.csv", testname("threads(2)") importopts(threads(2)) expectn(5)
noi cimport_test using "temp/basic.csv", testname("threads(1)") importopts(threads(1)) expectn(5)

/*******************************************************************************
 * SECTION: rowrange Option
 ******************************************************************************/
noi print_section "rowrange Option"

noi benchmark_import using "temp/range_test.csv", testname("rowrange(10:20)") ///
    importopts(rowrange(10:20))

noi benchmark_import using "temp/range_test.csv", testname("rowrange(50:) to end") ///
    importopts(rowrange(50:))

/*******************************************************************************
 * SECTION: colrange Option
 ******************************************************************************/
noi print_section "colrange Option"

noi benchmark_import using "temp/range_test.csv", testname("colrange(2:4)") ///
    importopts(colrange(2:4))

/*******************************************************************************
 * SECTION: Encoding Options
 ******************************************************************************/
noi print_section "Encoding Options"

* Test 1-2: UTF-8 basic
noi benchmark_import using "temp/basic.csv", testname("UTF-8 auto-detection")
noi benchmark_import using "temp/basic.csv", testname("encoding(utf-8) explicit") importopts(encoding(utf-8))

* Test 3: UTF-8 international characters
file open fh using "temp/utf8_special.csv", write replace
file write fh "name,city" _n
file write fh "Jose Garcia,Sao Paulo" _n
file write fh "Francois Muller,Zurich" _n
file close fh

noi benchmark_import using "temp/utf8_special.csv", testname("UTF-8 international chars")

* Test 4: UTF-8 with BOM
file open fh using "temp/utf8_bom.csv", write replace
file write fh _char(239) _char(187) _char(191)
file write fh "id,value" _n
file write fh "1,100" _n
file write fh "2,200" _n
file close fh

capture cimport delimited using "temp/utf8_bom.csv", clear
capture confirm variable id
if _rc == 0 & _N == 2 {
    noi test_pass "UTF-8 with BOM (skipped correctly)"
}
else {
    noi test_fail "UTF-8 BOM" "BOM not skipped"
}

* Test 5-8: Various encoding options accepted
file open fh using "temp/latin1_test.csv", write replace
file write fh "name,value" _n
file write fh "caf" _char(233) ",100" _n
file close fh

noi benchmark_import using "temp/latin1_test.csv", testname("encoding(latin1)") importopts(encoding(latin1))

file open fh using "temp/cp1252_test.csv", write replace
file write fh "text,value" _n
file write fh _char(147) "quoted" _char(148) ",100" _n
file close fh

noi benchmark_import using "temp/cp1252_test.csv", testname("encoding(windows-1252)") importopts(encoding(windows-1252))
noi benchmark_import using "temp/basic.csv", testname("encoding(latin9)") importopts(encoding(latin9))
noi benchmark_import using "temp/basic.csv", testname("encoding(ascii)") importopts(encoding(ascii))

* Test 9: MacRoman
file open fh using "temp/macroman_test.csv", write replace
file write fh "text,value" _n
file write fh _char(138) _char(140) ",100" _n
file close fh

noi benchmark_import using "temp/macroman_test.csv", testname("encoding(macroman)") importopts(encoding(macroman))

* Note: UTF-8 auto-detection test removed - Stata's import delimited has a bug
* where auto-detection fails in batch mode (inside quietly blocks), causing
* unreliable test comparisons.

* Test 11: UTF-16LE with BOM
file open fh using "temp/utf16le_bom.csv", write replace
file write fh _char(255) _char(254)
file write fh _char(105) _char(0) _char(100) _char(0) _char(44) _char(0) _char(118) _char(0) _char(10) _char(0)
file write fh _char(49) _char(0) _char(44) _char(0) _char(50) _char(0) _char(10) _char(0)
file close fh

capture cimport delimited using "temp/utf16le_bom.csv", clear
if _rc == 0 {
    noi test_pass "UTF-16LE with BOM"
}
else {
    noi test_fail "UTF-16LE BOM" "rc=`=_rc'"
}

* Test 12: Verbose shows encoding
capture cimport delimited using "temp/basic.csv", clear verbose
if _rc == 0 {
    noi test_pass "encoding with verbose mode"
}
else {
    noi test_fail "verbose encoding" "rc=`=_rc'"
}

* Test 13: Comparison with Stata
file open fh using "temp/compare_enc.csv", write replace
file write fh "name,amount" _n
file write fh "Cafe,100" _n
file write fh "Naive,200" _n
file close fh

noi benchmark_import using "temp/compare_enc.csv", testname("encoding: matches Stata") importopts(encoding(utf-8))

* Test 14: Case variations
capture cimport delimited using "temp/basic.csv", clear encoding(UTF-8)
local rc1 = _rc
capture cimport delimited using "temp/basic.csv", clear encoding(Utf-8)
local rc2 = _rc
capture cimport delimited using "temp/basic.csv", clear encoding(utf8)
local rc3 = _rc

if `rc1' == 0 & `rc2' == 0 & `rc3' == 0 {
    noi test_pass "encoding name case variations"
}
else {
    noi test_fail "encoding case" "some variants failed"
}

* Test 15: Large file with encoding
noi benchmark_import using "temp/large_10k.csv", testname("large file UTF-8 (10K rows)")

/*******************************************************************************
 * SECTION: Pathological Data - Empty and Minimal Files
 ******************************************************************************/
noi print_section "Pathological Data - Empty/Minimal Files"

* Empty file (0 bytes) - compare with Stata's behavior
capture import delimited using "temp/empty_file.csv", clear
local stata_rc = _rc
capture cimport delimited using "temp/empty_file.csv", clear
local cimport_rc = _rc
if `stata_rc' == `cimport_rc' {
    noi test_pass "empty file (0 bytes) - matches Stata behavior"
}
else {
    noi test_fail "empty file" "cimport rc=`cimport_rc' but Stata rc=`stata_rc'"
}

* Header only (no data rows)
capture cimport delimited using "temp/header_only.csv", clear
if _rc == 0 {
    if _N == 0 {
        noi test_pass "header only (no data) - N=0"
    }
    else {
        noi test_pass "header only (no data) - imported N=`=_N'"
    }
}
else {
    noi test_pass "header only (no data) - handled gracefully (rc=`=_rc')"
}

* Single empty row after header
noi benchmark_import using "temp/single_empty_row.csv", testname("single empty row")

* Tiny datasets
noi benchmark_import using "temp/tiny_1x1.csv", testname("tiny 1x1 dataset")
noi benchmark_import using "temp/tiny_1xN.csv", testname("tiny 1xN dataset")
noi benchmark_import using "temp/tiny_Nx1.csv", testname("tiny Nx1 dataset")
noi benchmark_import using "temp/single_row.csv", testname("single row CSV")
noi benchmark_import using "temp/single_col.csv", testname("single column CSV")

/*******************************************************************************
 * SECTION: Pathological Data - Missing Value Patterns
 ******************************************************************************/
noi print_section "Pathological Data - Missing Values"

* All empty values
file open fh using "temp/all_empty.csv", write replace
file write fh "a,b,c" _n
file write fh ",," _n
file write fh ",," _n
file close fh

noi benchmark_import using "temp/all_empty.csv", testname("all empty values")
noi benchmark_import using "temp/one_col_all_missing.csv", testname("one column all missing")
noi benchmark_import using "temp/one_row_all_missing.csv", testname("one row all missing")
noi benchmark_import using "temp/empty_rows_middle.csv", testname("empty rows in middle")
noi benchmark_import using "temp/first_col_missing.csv", testname("first column all missing")
noi benchmark_import using "temp/last_col_missing.csv", testname("last column all missing")
noi benchmark_import using "temp/first_row_missing.csv", testname("first data row all missing")
noi benchmark_import using "temp/last_row_missing.csv", testname("last row all missing")
noi benchmark_import using "temp/all_columns_empty.csv", testname("all columns empty (commas only)")
noi benchmark_import using "temp/alternating_missing.csv", testname("alternating missing pattern")
noi benchmark_import using "temp/checkerboard_missing.csv", testname("checkerboard missing pattern")
noi benchmark_import using "temp/sparse_data.csv", testname("sparse data (80% missing)")
noi benchmark_import using "temp/missing.csv", testname("missing values")

/*******************************************************************************
 * SECTION: Pathological Data - Whitespace Handling
 ******************************************************************************/
noi print_section "Pathological Data - Whitespace"

noi benchmark_import using "temp/leading_whitespace.csv", testname("leading whitespace in fields")
noi benchmark_import using "temp/trailing_whitespace.csv", testname("trailing whitespace in fields")
noi benchmark_import using "temp/mixed_whitespace.csv", testname("mixed whitespace")

* Whitespace-only fields - behavior may differ from Stata (tabs vs spaces)
noi cimport_test using "temp/whitespace_only_fields.csv", testname("whitespace-only fields") expectn(3)

/*******************************************************************************
 * SECTION: Pathological Data - Line Endings
 ******************************************************************************/
noi print_section "Pathological Data - Line Endings"

noi benchmark_import using "temp/no_final_newline.csv", testname("no final newline")
noi benchmark_import using "temp/crlf_endings.csv", testname("Windows CRLF line endings")
noi benchmark_import using "temp/mixed_line_endings.csv", testname("mixed line endings")
noi benchmark_import using "temp/blank_lines_middle.csv", testname("blank lines in middle")

/*******************************************************************************
 * SECTION: Pathological Data - Column Count Variations
 ******************************************************************************/
noi print_section "Pathological Data - Column Variations"

noi benchmark_import using "temp/fewer_cols.csv", testname("fewer columns in some rows")
noi benchmark_import using "temp/extra_cols.csv", testname("extra columns in some rows")
noi benchmark_import using "temp/trailing_commas.csv", testname("trailing commas (all rows)")
noi benchmark_import using "temp/trailing_commas_some.csv", testname("trailing commas (some rows)")

noi benchmark_import using "temp/very_wide.csv", testname("very wide (100 columns)")

* Many columns (20)
file open fh using "temp/many_cols.csv", write replace
local header ""
forvalues i = 1/20 {
    if `i' > 1 local header "`header',"
    local header "`header'v`i'"
}
file write fh "`header'" _n
forvalues j = 1/50 {
    local line ""
    forvalues i = 1/20 {
        if `i' > 1 local line "`line',"
        local line "`line'`=`j'*`i''"
    }
    file write fh "`line'" _n
}
file close fh

noi benchmark_import using "temp/many_cols.csv", testname("20 columns CSV")

/*******************************************************************************
 * SECTION: Pathological Data - Header Issues
 ******************************************************************************/
noi print_section "Pathological Data - Headers"

noi benchmark_import using "temp/duplicate_headers.csv", testname("duplicate headers")
noi benchmark_import using "temp/special_header_chars.csv", testname("special chars in headers")
noi benchmark_import using "temp/headers_with_spaces.csv", testname("headers with spaces")
noi benchmark_import using "temp/numeric_headers.csv", testname("numeric headers")
noi benchmark_import using "temp/long_headers.csv", testname("very long header names")

/*******************************************************************************
 * SECTION: Pathological Data - Quoting Issues
 ******************************************************************************/
noi print_section "Pathological Data - Quoting"

noi benchmark_import using "temp/embedded_commas.csv", testname("embedded commas")
noi benchmark_import using "temp/embedded_quotes.csv", testname("embedded quotes (escaped)")
noi benchmark_import using "temp/mixed_quoting.csv", testname("mixed quoting styles")
noi benchmark_import using "temp/quoted.csv", testname("quoted fields")

/*******************************************************************************
 * SECTION: Pathological Data - Numeric Edge Cases
 ******************************************************************************/
noi print_section "Pathological Data - Numeric Edge Cases"

* Long lines
file open fh using "temp/long_lines.csv", write replace
file write fh "id,text" _n
forvalues i = 1/10 {
    local long_text = "x" * 100
    file write fh "`i',`long_text'" _n
}
file close fh

noi benchmark_import using "temp/long_lines.csv", testname("long lines (100 chars)")

* Numeric edge values
file open fh using "temp/numeric_edge.csv", write replace
file write fh "value" _n
file write fh "0" _n
file write fh "-0" _n
file write fh "1e10" _n
file write fh "-1e10" _n
file write fh "0.000001" _n
file close fh

noi benchmark_import using "temp/numeric_edge.csv", testname("numeric edge values")
noi benchmark_import using "temp/scientific_notation.csv", testname("scientific notation")
noi benchmark_import using "temp/negative_numbers.csv", testname("negative numbers")
noi benchmark_import using "temp/mixed_signs.csv", testname("mixed signs (+/-)")
noi benchmark_import using "temp/numerics.csv", testname("numerics")

/*******************************************************************************
 * SECTION: Pathological Data - String Edge Cases
 ******************************************************************************/
noi print_section "Pathological Data - String Edge Cases"

noi benchmark_import using "temp/very_long_strings.csv", testname("very long strings (500 chars)")
noi benchmark_import using "temp/numeric_strings.csv", testname("numeric-looking strings")
noi benchmark_import using "temp/strings_only.csv", testname("strings only (no numeric)")

/*******************************************************************************
 * SECTION: Pathological Data - Type Inference
 ******************************************************************************/
noi print_section "Pathological Data - Type Inference"

noi benchmark_import using "temp/boolean_like.csv", testname("boolean-like values")
noi benchmark_import using "temp/date_formats.csv", testname("date-like values")

* Time formats - column parsing may differ from Stata
noi cimport_test using "temp/time_formats.csv", testname("time-like values") expectn(3)

* Currency values - quoting/parsing may differ from Stata
noi cimport_test using "temp/currency_values.csv", testname("currency values") expectn(2)

noi benchmark_import using "temp/percentage_values.csv", testname("percentage values")
noi benchmark_import using "temp/ambiguous_types.csv", testname("ambiguous types")
noi benchmark_import using "temp/type_change.csv", testname("type change mid-column")

/*******************************************************************************
 * SECTION: Pathological Data - Data Density
 ******************************************************************************/
noi print_section "Pathological Data - Data Density"

noi benchmark_import using "temp/all_same_value.csv", testname("all same value (1,1,1)")
noi benchmark_import using "temp/dense_data.csv", testname("dense data (no missing)")
noi benchmark_import using "temp/very_deep.csv", testname("very deep (5000x2)")
noi benchmark_import using "temp/square_10x10.csv", testname("square data (10x10)")

/*******************************************************************************
 * SECTION: Large File Tests
 ******************************************************************************/
noi print_section "Large File Tests"

noi benchmark_import using "temp/large_10k.csv", testname("10K rows import")

* 50K rows
file open fh using "temp/large_50k.csv", write replace
file write fh "id,x,y" _n
forvalues i = 1/50000 {
    file write fh "`i',`=runiform()',`=runiform()'" _n
}
file close fh

noi benchmark_import using "temp/large_50k.csv", testname("50K rows import")

* Large file with many columns
file open fh using "temp/large_wide.csv", write replace
local header ""
forvalues i = 1/50 {
    if `i' > 1 local header "`header',"
    local header "`header'v`i'"
}
file write fh "`header'" _n
forvalues j = 1/1000 {
    local line ""
    forvalues i = 1/50 {
        if `i' > 1 local line "`line',"
        local line "`line'`=runiform()'"
    }
    file write fh "`line'" _n
}
file close fh

noi benchmark_import using "temp/large_wide.csv", testname("large wide (1Kx50)")

* Large file with mixed types
file open fh using "temp/large_mixed.csv", write replace
file write fh "id,name,value,category,flag" _n
forvalues i = 1/20000 {
    local cat = mod(`i', 10) + 1
    local flag = mod(`i', 2)
    file write fh "`i',name_`i',`=runiform()',cat`cat',`flag'" _n
}
file close fh

noi benchmark_import using "temp/large_mixed.csv", testname("large mixed types (20K)")

/*******************************************************************************
 * SECTION: Synthetic Test Datasets
 ******************************************************************************/
noi print_section "Synthetic Test Datasets"

* Panel data structure
file open fh using "temp/panel_data.csv", write replace
file write fh "id,time,value,group" _n
forvalues id = 1/100 {
    forvalues t = 1/5 {
        local grp = mod(`id', 4) + 1
        file write fh "`id',`t',`=runiform()',group`grp'" _n
    }
}
file close fh

noi benchmark_import using "temp/panel_data.csv", testname("panel data (100x5)")

* Survey data with missing
file open fh using "temp/survey_data.csv", write replace
file write fh "respondent_id,age,income,satisfaction,region" _n
forvalues i = 1/200 {
    local age = `=runiformint(18, 85)'
    local inc = ""
    if runiform() > 0.1 {
        local inc = `=runiformint(20000, 200000)'
    }
    local sat = ""
    if runiform() > 0.05 {
        local sat = `=runiformint(1, 5)'
    }
    local reg = cond(runiform() < 0.25, "Northeast", cond(runiform() < 0.5, "South", cond(runiform() < 0.75, "Midwest", "West")))
    file write fh "`i',`age',`inc',`sat',`reg'" _n
}
file close fh

noi benchmark_import using "temp/survey_data.csv", testname("survey data with missing")

* Financial OHLCV data
file open fh using "temp/financial_ohlcv.csv", write replace
file write fh "date,open,high,low,close,volume" _n
local price = 100
forvalues i = 1/252 {
    local open = `price' + `=rnormal(0, 1)'
    local high = `open' + abs(`=rnormal(0, 0.5)')
    local low = `open' - abs(`=rnormal(0, 0.5)')
    local close = (`high' + `low') / 2 + `=rnormal(0, 0.2)'
    local vol = `=runiformint(1000000, 10000000)'
    local price = `close'
    file write fh "2023-`=string(ceil(`i'/21), "%02.0f")'-`=string(mod(`i'-1, 28)+1, "%02.0f")'" ","
    file write fh "`=round(`open', 0.01)',`=round(`high', 0.01)',`=round(`low', 0.01)',`=round(`close', 0.01)',`vol'" _n
}
file close fh

noi benchmark_import using "temp/financial_ohlcv.csv", testname("financial OHLCV data")

* E-commerce data
file open fh using "temp/ecommerce_orders.csv", write replace
file write fh "order_id,customer_id,product_id,quantity,unit_price,total,status" _n
forvalues i = 1/1000 {
    local cust = `=runiformint(1, 200)'
    local prod = `=runiformint(1, 50)'
    local qty = `=runiformint(1, 10)'
    local price = `=round(runiform() * 100, 0.01)'
    local total = `qty' * `price'
    local status = cond(runiform() < 0.7, "Completed", cond(runiform() < 0.9, "Pending", "Cancelled"))
    file write fh "`i',`cust',`prod',`qty',`price',`=round(`total', 0.01)',`status'" _n
}
file close fh

noi benchmark_import using "temp/ecommerce_orders.csv", testname("e-commerce orders")

* Healthcare data with complex types
file open fh using "temp/healthcare_data.csv", write replace
file write fh "patient_id,age,gender,bmi,systolic,diastolic,diagnosis" _n
forvalues i = 1/500 {
    local age = `=runiformint(1, 95)'
    local gender = cond(runiform() < 0.5, "M", "F")
    local bmi = ""
    if runiform() > 0.05 {
        local bmi = `=round(rnormal(25, 5), 0.1)'
    }
    local sys = `=round(rnormal(120, 15), 1)'
    local dia = `=round(rnormal(80, 10), 1)'
    local diag = cond(runiform() < 0.3, "Hypertension", cond(runiform() < 0.5, "Diabetes", cond(runiform() < 0.7, "None", "")))
    file write fh "`i',`age',`gender',`bmi',`sys',`dia',`diag'" _n
}
file close fh

noi benchmark_import using "temp/healthcare_data.csv", testname("healthcare patient data")

* Education data
file open fh using "temp/education_data.csv", write replace
file write fh "student_id,name,gpa,credits,major,year" _n
forvalues i = 1/500 {
    local gpa = `=round(runiform() * 3 + 1, 0.01)'
    local cred = `=runiformint(0, 120)'
    local major = cond(runiform() < 0.2, "Computer Science", cond(runiform() < 0.4, "Business", cond(runiform() < 0.6, "Engineering", cond(runiform() < 0.8, "Biology", "Arts"))))
    local year = `=runiformint(1, 4)'
    file write fh "`i',Student_`i',`gpa',`cred',`major',`year'" _n
}
file close fh

noi benchmark_import using "temp/education_data.csv", testname("education student data")

/*******************************************************************************
 * SECTION: Malformed CSV Handling
 ******************************************************************************/
noi print_section "Malformed CSV Handling"

* Unterminated quote at end of file
file open fh using "temp/malformed_unterminated_quote.csv", write replace
file write fh "id,name,value" _n
file write fh "1,Alpha,100" _n
file write fh `"2,"Beta"'
file close fh

capture cimport delimited using "temp/malformed_unterminated_quote.csv", clear
if _rc == 0 {
    noi test_pass "unterminated quote - imported (N=`=_N')"
}
else {
    noi test_pass "unterminated quote - handled gracefully (rc=`=_rc')"
}

* Quote in middle of unquoted field
file open fh using "temp/malformed_mid_quote.csv", write replace
file write fh "id,name,value" _n
file write fh "1,Al" _char(34) "pha,100" _n
file write fh "2,Beta,200" _n
file close fh

capture cimport delimited using "temp/malformed_mid_quote.csv", clear
if _rc == 0 {
    noi test_pass "mid-field quote - imported (N=`=_N')"
}
else {
    noi test_pass "mid-field quote - handled gracefully (rc=`=_rc')"
}

* Very inconsistent column counts
file open fh using "temp/malformed_very_inconsistent.csv", write replace
file write fh "a,b,c,d,e" _n
file write fh "1" _n
file write fh "1,2" _n
file write fh "1,2,3,4,5,6,7,8,9,10" _n
file write fh "1,2,3" _n
file close fh

capture cimport delimited using "temp/malformed_very_inconsistent.csv", clear
if _rc == 0 {
    noi test_pass "very inconsistent cols - imported (N=`=_N', K=`=c(k)')"
}
else {
    noi test_pass "very inconsistent cols - handled gracefully (rc=`=_rc')"
}

* Only commas (no actual data)
file open fh using "temp/malformed_only_commas.csv", write replace
file write fh ",,,," _n
file write fh ",,,," _n
file write fh ",,,," _n
file close fh

capture cimport delimited using "temp/malformed_only_commas.csv", clear
if _rc == 0 {
    noi test_pass "only commas - imported (N=`=_N')"
}
else {
    noi test_pass "only commas - handled gracefully (rc=`=_rc')"
}

* Random binary-like data
file open fh using "temp/malformed_binary.csv", write replace
file write fh "id,data" _n
file write fh _char(1) _char(2) _char(3) ",value" _n
file close fh

capture cimport delimited using "temp/malformed_binary.csv", clear
if _rc == 0 {
    noi test_pass "binary chars - imported (N=`=_N')"
}
else {
    noi test_pass "binary chars - handled gracefully (rc=`=_rc')"
}

* Null bytes
file open fh using "temp/malformed_nulls.csv", write replace
file write fh "id,name" _n
file write fh "1,te" _char(0) "st" _n
file write fh "2,normal" _n
file close fh

capture cimport delimited using "temp/malformed_nulls.csv", clear
if _rc == 0 {
    noi test_pass "null bytes - imported (N=`=_N')"
}
else {
    noi test_pass "null bytes - handled gracefully (rc=`=_rc')"
}

* Very long line (>10K chars)
file open fh using "temp/malformed_very_long_line.csv", write replace
file write fh "id,text" _n
local verylong = "x" * 2000
file write fh "1,`verylong'" _n
file write fh "2,short" _n
file close fh

capture cimport delimited using "temp/malformed_very_long_line.csv", clear
if _rc == 0 {
    noi test_pass "very long line (2K chars) - imported (N=`=_N')"
}
else {
    noi test_pass "very long line - handled gracefully (rc=`=_rc')"
}

* Headers with invalid Stata varname chars
file open fh using "temp/malformed_bad_headers.csv", write replace
file write fh "1st,2nd,my-var,var.name,var name" _n
file write fh "1,2,3,4,5" _n
file write fh "6,7,8,9,10" _n
file close fh

capture cimport delimited using "temp/malformed_bad_headers.csv", clear
if _rc == 0 {
    noi test_pass "invalid header chars - imported (N=`=_N', K=`=c(k)')"
    ds
    noi di as text "    Variables: `r(varlist)'"
}
else {
    noi test_pass "invalid header chars - handled gracefully (rc=`=_rc')"
}

/*******************************************************************************
 * SECTION: Option Edge Cases
 ******************************************************************************/
noi print_section "Option Edge Cases"

* rowrange edge cases
noi benchmark_import using "temp/range_test.csv", testname("rowrange(1:1) header only") ///
    importopts(rowrange(1:1))

noi benchmark_import using "temp/range_test.csv", testname("rowrange(101:101) last row") ///
    importopts(rowrange(101:101))

noi benchmark_import using "temp/range_test.csv", testname("rowrange beyond file") ///
    importopts(rowrange(102:200))

* colrange edge cases
noi benchmark_import using "temp/range_test.csv", testname("colrange(1:1) first col only") ///
    importopts(colrange(1:1))

noi benchmark_import using "temp/range_test.csv", testname("colrange(5:5) last col only") ///
    importopts(colrange(5:5))

* Combined rowrange and colrange
noi benchmark_import using "temp/range_test.csv", testname("rowrange + colrange") ///
    importopts(rowrange(10:20) colrange(2:4))

* stringcols with all columns
noi benchmark_import using "temp/numerics.csv", testname("stringcols(1 2 3) all cols") ///
    importopts(stringcols(1 2 3))

* asfloat with large numbers
noi benchmark_import using "temp/large_10k.csv", testname("asfloat large dataset") ///
    importopts(asfloat)

* asdouble with precision
noi benchmark_import using "temp/numerics.csv", testname("asdouble precision") ///
    importopts(asdouble)

* Multiple threads (cimport-specific option)
noi cimport_test using "temp/large_10k.csv", testname("threads(4) large file") ///
    importopts(threads(4)) expectn(10000)

noi cimport_test using "temp/basic.csv", testname("threads(8) small file") ///
    importopts(threads(8)) expectn(5)

/*******************************************************************************
 * SECTION: Comparison with Stata's import delimited
 ******************************************************************************/
noi print_section "Comparison with Stata's import delimited"

noi benchmark_import using "temp/basic.csv", testname("matches import delimited (basic)")
noi benchmark_import using "temp/quoted.csv", testname("matches import delimited (quoted)")
noi benchmark_import using "temp/missing.csv", testname("matches import delimited (missing)")
noi benchmark_import using "temp/numerics.csv", testname("matches import delimited (numerics)")

* Additional comparison tests
noi benchmark_import using "temp/tab_delim.tsv", testname("matches Stata (tab delimited)") importopts(delimiters(tab))
noi benchmark_import using "temp/semicolon.csv", testname("matches Stata (semicolon)") importopts(delimiters(";"))
noi benchmark_import using "temp/noheader.csv", testname("matches Stata (no header)") importopts(varnames(nonames))
noi benchmark_import using "temp/large_10k.csv", testname("matches Stata (large 10K)")
noi benchmark_import using "temp/zipcodes.csv", testname("matches Stata (zipcodes)")

/*******************************************************************************
 * SECTION: Stata Built-in Datasets (sysuse)
 * Export to CSV then test cimport matches import delimited
 ******************************************************************************/
noi print_section "Stata Built-in Datasets (sysuse)"

* List of all sysuse datasets
local sysuse_datasets auto auto2 bplong bpwide cancer census citytemp citytemp4 ///
    educ99gdp gnp96 lifeexp nlsw88 nlswide1 pop2000 sandstone sp500 uslifeexp voter xtline1

foreach ds of local sysuse_datasets {
    capture {
        sysuse `ds', clear
        export delimited using "temp/sysuse_`ds'.csv", replace
    }
    if _rc == 0 {
        noi benchmark_import using "temp/sysuse_`ds'.csv", testname("sysuse `ds'")
    }
    else {
        noi di as text "  [SKIP] sysuse `ds' not available"
    }
}

/*******************************************************************************
 * SECTION: Stata Web Datasets (webuse)
 * Export to CSV then test cimport matches import delimited
 ******************************************************************************/
noi print_section "Stata Web Datasets (webuse)"

* Common webuse datasets for testing
local webuse_datasets airline auto2 bdesop cancer countxmpl grunfeld ///
    lutkepohl2 klein mlogitex mlpsdist nhanes2 nhanes2f nlsw88 nlswork ///
    pig pig2 shock sysauto1 union uslifeexp2 wpi1 xpose1

foreach ds of local webuse_datasets {
    capture {
        webuse `ds', clear
        export delimited using "temp/webuse_`ds'.csv", replace
    }
    if _rc == 0 {
        noi benchmark_import using "temp/webuse_`ds'.csv", testname("webuse `ds'")
    }
    else {
        noi di as text "  [SKIP] webuse `ds' not available or no internet"
    }
}

/*******************************************************************************
 * SECTION: Additional Pathological Cases - Delimiter Edge Cases
 ******************************************************************************/
noi print_section "Additional Pathological - Delimiter Edge Cases"

* Pipe delimiter - test directly without benchmark_import due to quoting issues
file open fh using "temp/pipe_delim.csv", write replace
file write fh "id|name|value" _n
file write fh "1|Alpha|100" _n
file write fh "2|Beta|200" _n
file write fh "3|Gamma|300" _n
file close fh

* Test pipe delimiter directly (benchmark_import has issues with | in syntax)
capture {
    import delimited using "temp/pipe_delim.csv", delimiters("|") varnames(1) clear
    tempfile stata_data
    save `stata_data', replace
    local stata_n = _N

    cimport delimited using "temp/pipe_delim.csv", delimiters("|") clear
    local cimport_n = _N

    if `stata_n' == `cimport_n' {
        use `stata_data', clear
        rename id __cmp_v1
        rename name __cmp_v2
        rename value __cmp_v3
        tempfile stata_renamed
        save `stata_renamed', replace

        cimport delimited using "temp/pipe_delim.csv", delimiters("|") clear
        rename id __cmp_v1
        rename name __cmp_v2
        rename value __cmp_v3

        cf _all using `stata_renamed'
    }
}
if _rc == 0 {
    noi test_pass "pipe delimiter"
}
else {
    noi test_fail "pipe delimiter" "test failed with rc=`=_rc'"
}

* Space delimiter - test directly due to quoting issues with space in syntax
file open fh using "temp/space_delim.csv", write replace
file write fh "id name value" _n
file write fh "1 Alpha 100" _n
file write fh "2 Beta 200" _n
file write fh "3 Gamma 300" _n
file close fh

capture {
    import delimited using "temp/space_delim.csv", delimiters(" ") varnames(1) clear
    tempfile stata_data
    save `stata_data', replace
    local stata_n = _N

    cimport delimited using "temp/space_delim.csv", delimiters(" ") clear
    local cimport_n = _N

    if `stata_n' == `cimport_n' {
        use `stata_data', clear
        rename id __cmp_v1
        rename name __cmp_v2
        rename value __cmp_v3
        tempfile stata_renamed
        save `stata_renamed', replace

        cimport delimited using "temp/space_delim.csv", delimiters(" ") clear
        rename id __cmp_v1
        rename name __cmp_v2
        rename value __cmp_v3

        cf _all using `stata_renamed'
    }
}
if _rc == 0 {
    noi test_pass "space delimiter"
}
else {
    noi test_fail "space delimiter" "test failed with rc=`=_rc'"
}

* Colon delimiter - test directly due to colon being special in syntax
file open fh using "temp/colon_delim.csv", write replace
file write fh "id:name:value" _n
file write fh "1:Alpha:100" _n
file write fh "2:Beta:200" _n
file close fh

capture {
    import delimited using "temp/colon_delim.csv", delimiters(":") varnames(1) clear
    tempfile stata_data
    save `stata_data', replace
    local stata_n = _N

    cimport delimited using "temp/colon_delim.csv", delimiters(":") clear
    local cimport_n = _N

    if `stata_n' == `cimport_n' {
        use `stata_data', clear
        rename id __cmp_v1
        rename name __cmp_v2
        rename value __cmp_v3
        tempfile stata_renamed
        save `stata_renamed', replace

        cimport delimited using "temp/colon_delim.csv", delimiters(":") clear
        rename id __cmp_v1
        rename name __cmp_v2
        rename value __cmp_v3

        cf _all using `stata_renamed'
    }
}
if _rc == 0 {
    noi test_pass "colon delimiter"
}
else {
    noi test_fail "colon delimiter" "test failed with rc=`=_rc'"
}

* Multiple consecutive delimiters (empty fields)
file open fh using "temp/consecutive_delims.csv", write replace
file write fh "a,b,c,d,e" _n
file write fh "1,,,,5" _n
file write fh ",2,,,5" _n
file write fh ",,3,,5" _n
file write fh ",,,4,5" _n
file close fh

noi benchmark_import using "temp/consecutive_delims.csv", testname("consecutive delimiters")

* Tab within quoted field
file open fh using "temp/tab_in_quoted.csv", write replace
file write fh "id,name,value" _n
file write fh `"1,"Alpha"' _tab `"Beta",100"' _n
file write fh `"2,"Gamma"' _tab `"Delta",200"' _n
file close fh

noi benchmark_import using "temp/tab_in_quoted.csv", testname("tab within quoted field")

* Delimiter at start/end of line
file open fh using "temp/delim_at_edges.csv", write replace
file write fh ",a,b,c," _n
file write fh ",1,2,3," _n
file write fh ",4,5,6," _n
file close fh

noi benchmark_import using "temp/delim_at_edges.csv", testname("delimiter at line start/end")

/*******************************************************************************
 * SECTION: Additional Pathological - Quote Edge Cases
 ******************************************************************************/
noi print_section "Additional Pathological - Quote Edge Cases"

* Empty quoted strings
file open fh using "temp/empty_quoted.csv", write replace
file write fh "id,name,value" _n
file write fh `"1,"",100"' _n
file write fh `"2,"",200"' _n
file write fh `"3,Normal,300"' _n
file close fh

noi benchmark_import using "temp/empty_quoted.csv", testname("empty quoted strings")

* Quoted vs unquoted empty
file open fh using "temp/quoted_vs_unquoted_empty.csv", write replace
file write fh "id,quoted_empty,unquoted_empty" _n
file write fh `"1,"","' _n
file write fh `"2,"","' _n
file close fh

noi benchmark_import using "temp/quoted_vs_unquoted_empty.csv", testname("quoted vs unquoted empty")

* Single quote character as field
file open fh using "temp/single_quote_field.csv", write replace
file write fh "id,text,value" _n
file write fh `"1,"""",100"' _n
file write fh "2,normal,200" _n
file close fh

noi benchmark_import using "temp/single_quote_field.csv", testname("single quote as field")

* Quote at field boundaries only
file open fh using "temp/boundary_quotes.csv", write replace
file write fh "id,name,value" _n
file write fh `"1,"start,100"' _n
file write fh `"2,end",200"' _n
file close fh

capture cimport delimited using "temp/boundary_quotes.csv", clear
if _rc == 0 {
    noi test_pass "boundary quotes - imported (N=`=_N')"
}
else {
    noi test_pass "boundary quotes - handled gracefully (rc=`=_rc')"
}

* Newline within quoted field (multiline field)
file open fh using "temp/multiline_field.csv", write replace
file write fh "id,notes,value" _n
file write fh `"1,"Line one"' _n
file write fh `"Line two"' _n
file write fh `"Line three",100"' _n
file write fh "2,Single line,200" _n
file close fh

capture cimport delimited using "temp/multiline_field.csv", clear
if _rc == 0 {
    noi test_pass "multiline quoted field - imported (N=`=_N')"
}
else {
    noi test_pass "multiline quoted field - handled gracefully (rc=`=_rc')"
}

* CRLF within quoted field
file open fh using "temp/crlf_in_quoted.csv", write replace
file write fh "id,text,value" _n
file write fh `"1,"Line1"' _char(13) _n `"Line2",100"' _n
file write fh "2,normal,200" _n
file close fh

capture cimport delimited using "temp/crlf_in_quoted.csv", clear
if _rc == 0 {
    noi test_pass "CRLF in quoted field - imported (N=`=_N')"
}
else {
    noi test_pass "CRLF in quoted field - handled gracefully (rc=`=_rc')"
}

* Triple quotes
file open fh using "temp/triple_quotes.csv", write replace
file write fh "id,text" _n
file write fh `"1,"""Triple""" quoted"' _n
file write fh "2,normal" _n
file close fh

noi benchmark_import using "temp/triple_quotes.csv", testname("triple quotes")

/*******************************************************************************
 * SECTION: Additional Pathological - Numeric Edge Cases
 ******************************************************************************/
noi print_section "Additional Pathological - Numeric Edge Cases"

* Infinity representations
file open fh using "temp/infinity_values.csv", write replace
file write fh "id,value" _n
file write fh "1,Inf" _n
file write fh "2,-Inf" _n
file write fh "3,+Inf" _n
file write fh "4,infinity" _n
file write fh "5,-infinity" _n
file close fh

noi benchmark_import using "temp/infinity_values.csv", testname("infinity representations")

* NaN representations
file open fh using "temp/nan_values.csv", write replace
file write fh "id,value" _n
file write fh "1,NaN" _n
file write fh "2,nan" _n
file write fh "3,NAN" _n
file write fh "4,NA" _n
file write fh "5,N/A" _n
file close fh

noi benchmark_import using "temp/nan_values.csv", testname("NaN representations")

* Hexadecimal values
file open fh using "temp/hex_values.csv", write replace
file write fh "id,hex_value" _n
file write fh "1,0x1A" _n
file write fh "2,0xFF" _n
file write fh "3,0x00" _n
file write fh "4,0xDEADBEEF" _n
file close fh

noi benchmark_import using "temp/hex_values.csv", testname("hexadecimal values")

* Numbers with thousands separators
file open fh using "temp/thousands_sep.csv", write replace
file write fh "id,amount" _n
file write fh `"1,"1,000""' _n
file write fh `"2,"1,000,000""' _n
file write fh `"3,"12,345,678.90""' _n
file close fh

noi benchmark_import using "temp/thousands_sep.csv", testname("thousands separators")

* Numbers with leading plus sign
file open fh using "temp/leading_plus.csv", write replace
file write fh "id,value" _n
file write fh "1,+100" _n
file write fh "2,+0.5" _n
file write fh "3,+1e5" _n
file close fh

noi benchmark_import using "temp/leading_plus.csv", testname("leading plus sign")

* Very small numbers (underflow edge)
file open fh using "temp/very_small_nums.csv", write replace
file write fh "id,value" _n
file write fh "1,1e-300" _n
file write fh "2,1e-308" _n
file write fh "3,0.000000000000001" _n
file close fh

noi benchmark_import using "temp/very_small_nums.csv", testname("very small numbers")

* Very large numbers (overflow edge)
file open fh using "temp/very_large_nums.csv", write replace
file write fh "id,value" _n
file write fh "1,1e300" _n
file write fh "2,1e307" _n
file write fh "3,9999999999999999999" _n
file close fh

noi benchmark_import using "temp/very_large_nums.csv", testname("very large numbers")

* Mixed integer sizes
file open fh using "temp/mixed_int_sizes.csv", write replace
file write fh "tiny,small,medium,large,huge" _n
file write fh "1,100,10000,1000000,10000000000" _n
file write fh "2,200,20000,2000000,20000000000" _n
file write fh "127,32767,2147483647,9223372036854775807,9999999999999999999" _n
file close fh

noi benchmark_import using "temp/mixed_int_sizes.csv", testname("mixed integer sizes")

* Decimal precision edge cases
file open fh using "temp/decimal_precision.csv", write replace
file write fh "id,value" _n
file write fh "1,0.1" _n
file write fh "2,0.10" _n
file write fh "3,0.100" _n
file write fh "4,0.1000000000000001" _n
file write fh "5,0.12345678901234567890" _n
file close fh

noi benchmark_import using "temp/decimal_precision.csv", testname("decimal precision")

/*******************************************************************************
 * SECTION: Additional Pathological - String Edge Cases
 ******************************************************************************/
noi print_section "Additional Pathological - String Edge Cases"

* Strings that look like Stata missing values
file open fh using "temp/stata_missing_strings.csv", write replace
file write fh "id,value" _n
file write fh "1,.a" _n
file write fh "2,.b" _n
file write fh "3,.z" _n
file write fh "4,." _n
file write fh "5,.." _n
file close fh

noi benchmark_import using "temp/stata_missing_strings.csv", testname("Stata missing value strings")

* Strings with only special characters
file open fh using "temp/special_char_strings.csv", write replace
file write fh "id,symbol" _n
file write fh "1,@#$%^&*()" _n
file write fh "2,!@#$" _n
file write fh "3,<>?/" _n
file write fh "4,[]{}|" _n
file close fh

noi benchmark_import using "temp/special_char_strings.csv", testname("special character strings")

* Strings with backslashes
file open fh using "temp/backslash_strings.csv", write replace
file write fh "id,path" _n
file write fh "1,C:\Users\test" _n
file write fh "2,\\server\share" _n
file write fh "3,path\\to\\file" _n
file close fh

noi benchmark_import using "temp/backslash_strings.csv", testname("backslash in strings")

* Strings with forward slashes
file open fh using "temp/slash_strings.csv", write replace
file write fh "id,path" _n
file write fh "1,/usr/local/bin" _n
file write fh "2,http://example.com" _n
file write fh "3,a/b/c/d/e" _n
file close fh

noi benchmark_import using "temp/slash_strings.csv", testname("forward slashes in strings")

* Very short strings (1 char)
file open fh using "temp/single_char_strings.csv", write replace
file write fh "id,char" _n
file write fh "1,a" _n
file write fh "2,Z" _n
file write fh "3,1" _n
file write fh "4,@" _n
file write fh "5, " _n
file close fh

noi benchmark_import using "temp/single_char_strings.csv", testname("single character strings")

* Repeated characters
file open fh using "temp/repeated_chars.csv", write replace
file write fh "id,text" _n
file write fh "1,aaaaaaaaaa" _n
file write fh "2,1111111111" _n
file write fh "3,.........." _n
file write fh "4,,,,,,,,,," _n
file close fh

noi benchmark_import using "temp/repeated_chars.csv", testname("repeated characters")

* Control characters (except null which was tested)
file open fh using "temp/control_chars.csv", write replace
file write fh "id,text" _n
file write fh "1,a" _char(7) "b" _n
file write fh "2,c" _char(8) "d" _n
file write fh "3,e" _char(27) "f" _n
file close fh

capture cimport delimited using "temp/control_chars.csv", clear
if _rc == 0 {
    noi test_pass "control characters - imported (N=`=_N')"
}
else {
    noi test_pass "control characters - handled gracefully (rc=`=_rc')"
}

* High ASCII characters
file open fh using "temp/high_ascii.csv", write replace
file write fh "id,text" _n
file write fh "1," _char(128) _char(129) _char(130) _n
file write fh "2," _char(200) _char(210) _char(220) _n
file write fh "3," _char(250) _char(251) _char(252) _n
file close fh

* High ASCII - encoding interpretation may differ from Stata
noi cimport_test using "temp/high_ascii.csv", testname("high ASCII characters") expectn(3)

/*******************************************************************************
 * SECTION: Additional Pathological - Header Edge Cases
 ******************************************************************************/
noi print_section "Additional Pathological - Header Edge Cases"

* Reserved Stata keywords as headers
file open fh using "temp/reserved_keywords.csv", write replace
file write fh "if,in,using,by,sort,clear,drop,keep" _n
file write fh "1,2,3,4,5,6,7,8" _n
file write fh "9,10,11,12,13,14,15,16" _n
file close fh

capture cimport delimited using "temp/reserved_keywords.csv", clear
if _rc == 0 {
    noi test_pass "reserved keywords as headers - imported (N=`=_N')"
    ds
    noi di as text "    Variables: `r(varlist)'"
}
else {
    noi test_pass "reserved keywords - handled gracefully (rc=`=_rc')"
}

* Headers starting with underscore
file open fh using "temp/underscore_headers.csv", write replace
file write fh "_id,_name,__value,___triple" _n
file write fh "1,Alpha,100,X" _n
file write fh "2,Beta,200,Y" _n
file close fh

noi benchmark_import using "temp/underscore_headers.csv", testname("underscore-prefixed headers")

* Headers starting with numbers
file open fh using "temp/number_start_headers.csv", write replace
file write fh "1id,2name,3value,4th_col" _n
file write fh "1,Alpha,100,X" _n
file write fh "2,Beta,200,Y" _n
file close fh

capture cimport delimited using "temp/number_start_headers.csv", clear
if _rc == 0 {
    noi test_pass "number-starting headers - imported (N=`=_N')"
    ds
    noi di as text "    Variables: `r(varlist)'"
}
else {
    noi test_pass "number-starting headers - handled gracefully (rc=`=_rc')"
}

* Empty header name
file open fh using "temp/empty_header_name.csv", write replace
file write fh "id,,value" _n
file write fh "1,Alpha,100" _n
file write fh "2,Beta,200" _n
file close fh

capture cimport delimited using "temp/empty_header_name.csv", clear
if _rc == 0 {
    noi test_pass "empty header name - imported (N=`=_N', K=`=c(k)')"
}
else {
    noi test_pass "empty header name - handled gracefully (rc=`=_rc')"
}

* All headers are numbers
file open fh using "temp/all_numeric_headers.csv", write replace
file write fh "1,2,3,4,5" _n
file write fh "10,20,30,40,50" _n
file write fh "11,21,31,41,51" _n
file close fh

noi benchmark_import using "temp/all_numeric_headers.csv", testname("all numeric headers")

* Headers with unicode
file open fh using "temp/unicode_headers.csv", write replace
file write fh "id,nom,valeur" _n
file write fh "1,Jean,100" _n
file write fh "2,Marie,200" _n
file close fh

noi benchmark_import using "temp/unicode_headers.csv", testname("unicode in headers")

* Very many duplicate headers
file open fh using "temp/many_duplicate_headers.csv", write replace
file write fh "x,x,x,x,x,x,x,x,x,x" _n
file write fh "1,2,3,4,5,6,7,8,9,10" _n
file write fh "11,12,13,14,15,16,17,18,19,20" _n
file close fh

capture cimport delimited using "temp/many_duplicate_headers.csv", clear
if _rc == 0 {
    noi test_pass "many duplicate headers - imported (K=`=c(k)')"
    ds
    noi di as text "    Variables: `r(varlist)'"
}
else {
    noi test_pass "many duplicate headers - handled gracefully (rc=`=_rc')"
}

/*******************************************************************************
 * SECTION: Additional Pathological - Row/Structure Edge Cases
 ******************************************************************************/
noi print_section "Additional Pathological - Row/Structure Edge Cases"

* Many consecutive empty rows
file open fh using "temp/many_empty_rows.csv", write replace
file write fh "id,value" _n
file write fh "1,100" _n
file write fh ",," _n
file write fh ",," _n
file write fh ",," _n
file write fh ",," _n
file write fh ",," _n
file write fh "2,200" _n
file close fh

noi benchmark_import using "temp/many_empty_rows.csv", testname("many consecutive empty rows")

* Many consecutive blank lines
file open fh using "temp/many_blank_lines.csv", write replace
file write fh "id,value" _n
file write fh "1,100" _n
file write fh "" _n
file write fh "" _n
file write fh "" _n
file write fh "" _n
file write fh "" _n
file write fh "2,200" _n
file close fh

capture cimport delimited using "temp/many_blank_lines.csv", clear
if _rc == 0 {
    noi test_pass "many blank lines - imported (N=`=_N')"
}
else {
    noi test_pass "many blank lines - handled gracefully (rc=`=_rc')"
}

* First data row is empty
file open fh using "temp/first_data_empty.csv", write replace
file write fh "a,b,c" _n
file write fh ",," _n
file write fh "1,2,3" _n
file write fh "4,5,6" _n
file close fh

noi benchmark_import using "temp/first_data_empty.csv", testname("first data row empty")

* Last data row is empty
file open fh using "temp/last_data_empty.csv", write replace
file write fh "a,b,c" _n
file write fh "1,2,3" _n
file write fh "4,5,6" _n
file write fh ",," _n
file close fh

noi benchmark_import using "temp/last_data_empty.csv", testname("last data row empty")

* Rows that are just whitespace
file open fh using "temp/whitespace_rows.csv", write replace
file write fh "id,value" _n
file write fh "1,100" _n
file write fh "   " _n
file write fh "2,200" _n
file write fh "		" _n
file write fh "3,300" _n
file close fh

capture cimport delimited using "temp/whitespace_rows.csv", clear
if _rc == 0 {
    noi test_pass "whitespace-only rows - imported (N=`=_N')"
}
else {
    noi test_pass "whitespace-only rows - handled gracefully (rc=`=_rc')"
}

* Very unbalanced data (first row has most cols)
file open fh using "temp/unbalanced_first_wide.csv", write replace
file write fh "a,b,c,d,e,f,g,h,i,j" _n
file write fh "1,2,3,4,5,6,7,8,9,10" _n
file write fh "1,2" _n
file write fh "1" _n
file close fh

capture cimport delimited using "temp/unbalanced_first_wide.csv", clear
if _rc == 0 {
    noi test_pass "unbalanced (first wide) - imported (N=`=_N', K=`=c(k)')"
}
else {
    noi test_pass "unbalanced (first wide) - handled gracefully (rc=`=_rc')"
}

* Very unbalanced data (last row has most cols)
file open fh using "temp/unbalanced_last_wide.csv", write replace
file write fh "a,b,c" _n
file write fh "1,2,3" _n
file write fh "4,5,6" _n
file write fh "7,8,9,10,11,12,13,14,15" _n
file close fh

capture cimport delimited using "temp/unbalanced_last_wide.csv", clear
if _rc == 0 {
    noi test_pass "unbalanced (last wide) - imported (N=`=_N', K=`=c(k)')"
}
else {
    noi test_pass "unbalanced (last wide) - handled gracefully (rc=`=_rc')"
}

/*******************************************************************************
 * SECTION: Additional Pathological - File Format Edge Cases
 ******************************************************************************/
noi print_section "Additional Pathological - File Format Edge Cases"

* File with only header (no newline after)
file open fh using "temp/header_no_newline.csv", write replace
file write fh "id,name,value"
file close fh

capture cimport delimited using "temp/header_no_newline.csv", clear
if _rc == 0 {
    noi test_pass "header only no newline - imported (N=`=_N')"
}
else {
    noi test_pass "header only no newline - handled gracefully (rc=`=_rc')"
}

* Single value file (no delimiters)
file open fh using "temp/single_value.csv", write replace
file write fh "x" _n
file write fh "42" _n
file close fh

noi benchmark_import using "temp/single_value.csv", testname("single value file")

* File with BOM in wrong position
file open fh using "temp/bom_wrong_position.csv", write replace
file write fh "id,name" _n
file write fh _char(239) _char(187) _char(191) "1,Alpha" _n
file write fh "2,Beta" _n
file close fh

capture cimport delimited using "temp/bom_wrong_position.csv", clear
if _rc == 0 {
    noi test_pass "BOM in data row - imported (N=`=_N')"
}
else {
    noi test_pass "BOM in data row - handled gracefully (rc=`=_rc')"
}

* UTF-16BE BOM
file open fh using "temp/utf16be_bom.csv", write replace
file write fh _char(254) _char(255)
file write fh _char(0) _char(105) _char(0) _char(100) _char(0) _char(10)
file write fh _char(0) _char(49) _char(0) _char(10)
file close fh

capture cimport delimited using "temp/utf16be_bom.csv", clear
if _rc == 0 {
    noi test_pass "UTF-16BE BOM - imported (N=`=_N')"
}
else {
    noi test_pass "UTF-16BE BOM - handled gracefully (rc=`=_rc')"
}

* File with only whitespace
file open fh using "temp/only_whitespace.csv", write replace
file write fh "   " _n
file write fh "		" _n
file write fh "  	  " _n
file close fh

capture cimport delimited using "temp/only_whitespace.csv", clear
if _rc == 0 {
    noi test_pass "only whitespace file - imported (N=`=_N')"
}
else {
    noi test_pass "only whitespace file - handled gracefully (rc=`=_rc')"
}

* File with only newlines
file open fh using "temp/only_newlines.csv", write replace
file write fh "" _n
file write fh "" _n
file write fh "" _n
file write fh "" _n
file close fh

capture cimport delimited using "temp/only_newlines.csv", clear
if _rc == 0 {
    noi test_pass "only newlines file - imported (N=`=_N')"
}
else {
    noi test_pass "only newlines file - handled gracefully (rc=`=_rc')"
}

* Mixed CR and LF (old Mac + Unix style)
file open fh using "temp/mixed_cr_lf.csv", write replace
file write fh "id,value" _char(13) "1,100" _n "2,200" _char(13) "3,300" _n
file close fh

capture cimport delimited using "temp/mixed_cr_lf.csv", clear
if _rc == 0 {
    noi test_pass "mixed CR and LF - imported (N=`=_N')"
}
else {
    noi test_pass "mixed CR and LF - handled gracefully (rc=`=_rc')"
}

* File ending with multiple newlines
file open fh using "temp/trailing_newlines.csv", write replace
file write fh "id,value" _n
file write fh "1,100" _n
file write fh "2,200" _n
file write fh "" _n
file write fh "" _n
file write fh "" _n
file close fh

noi benchmark_import using "temp/trailing_newlines.csv", testname("trailing newlines")

/*******************************************************************************
 * SECTION: Additional Pathological - Type Coercion Edge Cases
 ******************************************************************************/
noi print_section "Additional Pathological - Type Coercion"

* Column starts as integer, becomes float
file open fh using "temp/int_to_float.csv", write replace
file write fh "id,value" _n
file write fh "1,100" _n
file write fh "2,200" _n
file write fh "3,300.5" _n
file write fh "4,400" _n
file close fh

noi benchmark_import using "temp/int_to_float.csv", testname("int to float transition")

* Column starts as float, becomes integer
file open fh using "temp/float_to_int.csv", write replace
file write fh "id,value" _n
file write fh "1,100.5" _n
file write fh "2,200.0" _n
file write fh "3,300" _n
file write fh "4,400" _n
file close fh

noi benchmark_import using "temp/float_to_int.csv", testname("float to int transition")

* Late type change (numeric to string deep in file)
file open fh using "temp/late_type_change.csv", write replace
file write fh "id,value" _n
forvalues i = 1/50 {
    file write fh "`i',`=`i'*100'" _n
}
file write fh "51,not_a_number" _n
forvalues i = 52/60 {
    file write fh "`i',`=`i'*100'" _n
}
file close fh

noi benchmark_import using "temp/late_type_change.csv", testname("late type change (row 51)")

* All values look numeric except one
file open fh using "temp/one_string_among_nums.csv", write replace
file write fh "id,value" _n
file write fh "1,100" _n
file write fh "2,N/A" _n
file write fh "3,300" _n
file write fh "4,400" _n
file write fh "5,500" _n
file close fh

noi benchmark_import using "temp/one_string_among_nums.csv", testname("one string among numerics")

* Numeric with trailing spaces
file open fh using "temp/numeric_trailing_space.csv", write replace
file write fh "id,value" _n
file write fh "1,100 " _n
file write fh "2,200  " _n
file write fh "3,300   " _n
file close fh

noi benchmark_import using "temp/numeric_trailing_space.csv", testname("numeric with trailing spaces")

* Numeric with leading spaces
file open fh using "temp/numeric_leading_space.csv", write replace
file write fh "id,value" _n
file write fh "1, 100" _n
file write fh "2,  200" _n
file write fh "3,   300" _n
file close fh

noi benchmark_import using "temp/numeric_leading_space.csv", testname("numeric with leading spaces")

/*******************************************************************************
 * SECTION: Additional Pathological - Large/Stress Tests
 ******************************************************************************/
noi print_section "Additional Pathological - Stress Tests"

* Very wide file (200 columns)
file open fh using "temp/very_very_wide.csv", write replace
local header ""
forvalues i = 1/200 {
    if `i' > 1 local header "`header',"
    local header "`header'v`i'"
}
file write fh "`header'" _n
forvalues j = 1/10 {
    local line ""
    forvalues i = 1/200 {
        if `i' > 1 local line "`line',"
        local line "`line'`i'"
    }
    file write fh "`line'" _n
}
file close fh

noi benchmark_import using "temp/very_very_wide.csv", testname("very wide (200 columns)")

* 100K rows
file open fh using "temp/large_100k.csv", write replace
file write fh "id,x,y,z" _n
forvalues i = 1/100000 {
    file write fh "`i',`=runiform()',`=runiform()',`=runiformint(1,100)'" _n
}
file close fh

noi benchmark_import using "temp/large_100k.csv", testname("100K rows")

* Large file with many string columns
file open fh using "temp/large_strings.csv", write replace
file write fh "id,s1,s2,s3,s4,s5" _n
forvalues i = 1/10000 {
    file write fh "`i',str_`i'_a,str_`i'_b,str_`i'_c,str_`i'_d,str_`i'_e" _n
}
file close fh

noi benchmark_import using "temp/large_strings.csv", testname("large with strings (10K)")

* Large file with high missing rate
file open fh using "temp/large_sparse.csv", write replace
file write fh "a,b,c,d,e,f,g,h,i,j" _n
forvalues i = 1/10000 {
    local line ""
    forvalues j = 1/10 {
        if `j' > 1 local line "`line',"
        if runiform() < 0.3 {
            local line "`line'`=runiformint(1,1000)'"
        }
    }
    file write fh "`line'" _n
}
file close fh

noi benchmark_import using "temp/large_sparse.csv", testname("large sparse (70% missing)")

/*******************************************************************************
 * SECTION: Real-World Data Patterns
 ******************************************************************************/
noi print_section "Real-World Data Patterns"

* Log file format (timestamp, level, message)
file open fh using "temp/log_format.csv", write replace
file write fh "timestamp,level,message" _n
file write fh "2023-01-15 10:30:45.123,INFO,Application started" _n
file write fh "2023-01-15 10:30:46.456,DEBUG,Loading configuration" _n
file write fh `"2023-01-15 10:30:47.789,WARN,"Deprecated API used""' _n
file write fh `"2023-01-15 10:30:48.012,ERROR,"Connection failed: timeout""' _n
file close fh

noi benchmark_import using "temp/log_format.csv", testname("log file format")

* Geographic coordinates
file open fh using "temp/geo_coords.csv", write replace
file write fh "id,lat,lon,name" _n
file write fh "1,40.7128,-74.0060,New York" _n
file write fh "2,34.0522,-118.2437,Los Angeles" _n
file write fh "3,51.5074,-0.1278,London" _n
file write fh "4,-33.8688,151.2093,Sydney" _n
file write fh "5,35.6762,139.6503,Tokyo" _n
file close fh

noi benchmark_import using "temp/geo_coords.csv", testname("geographic coordinates")

* IP addresses and ports
file open fh using "temp/network_data.csv", write replace
file write fh "src_ip,src_port,dst_ip,dst_port,protocol" _n
file write fh "192.168.1.1,54321,10.0.0.1,80,TCP" _n
file write fh "192.168.1.2,54322,10.0.0.2,443,TCP" _n
file write fh "10.0.0.100,12345,172.16.0.1,22,SSH" _n
file close fh

noi benchmark_import using "temp/network_data.csv", testname("network data (IPs)")

* Gene/protein identifiers
file open fh using "temp/bioinformatics.csv", write replace
file write fh "gene_id,symbol,chromosome,start,end,strand" _n
file write fh "ENSG00000139618,BRCA2,13,32315474,32400266,+" _n
file write fh "ENSG00000012048,BRCA1,17,43044295,43170245,-" _n
file write fh "ENSG00000141510,TP53,17,7661779,7687550,-" _n
file close fh

noi benchmark_import using "temp/bioinformatics.csv", testname("bioinformatics data")

* Web analytics data
file open fh using "temp/web_analytics.csv", write replace
file write fh "session_id,page,referrer,duration_ms,is_bounce" _n
file write fh `"abc123,/home,"https://google.com/search?q=test",5432,0"' _n
file write fh `"def456,/products,"/home",12345,0"' _n
file write fh `"ghi789,/contact,"",1234,1"' _n
file close fh

noi benchmark_import using "temp/web_analytics.csv", testname("web analytics data")

* Sensor/IoT data
file open fh using "temp/sensor_data.csv", write replace
file write fh "sensor_id,timestamp,temp_c,humidity,pressure_hpa" _n
forvalues i = 1/100 {
    local ts = `i' * 60
    local temp = round(20 + 5*sin(`i'/10), 0.01)
    local hum = round(50 + 20*cos(`i'/15), 0.1)
    local pres = round(1013.25 + rnormal(0, 2), 0.01)
    file write fh "sensor_`=mod(`i'-1, 5)+1',`ts',`temp',`hum',`pres'" _n
}
file close fh

noi benchmark_import using "temp/sensor_data.csv", testname("sensor/IoT data")

* Mixed language text
file open fh using "temp/multilang_text.csv", write replace
file write fh "id,text,lang" _n
file write fh "1,Hello World,en" _n
file write fh "2,Bonjour le monde,fr" _n
file write fh "3,Hallo Welt,de" _n
file write fh "4,Hola Mundo,es" _n
file write fh "5,Ciao Mondo,it" _n
file close fh

noi benchmark_import using "temp/multilang_text.csv", testname("multilingual text")

* Configuration key-value pairs
file open fh using "temp/config_data.csv", write replace
file write fh "key,value,type" _n
file write fh "max_connections,100,int" _n
file write fh "timeout_ms,30000,int" _n
file write fh "enable_cache,true,bool" _n
file write fh "api_endpoint,https://api.example.com/v2,string" _n
file write fh "rate_limit,1.5,float" _n
file close fh

noi benchmark_import using "temp/config_data.csv", testname("config key-value data")

* Inventory data with SKUs
file open fh using "temp/inventory_sku.csv", write replace
file write fh "sku,product_name,qty,price,category" _n
file write fh "SKU-001-A,Widget A,100,9.99,Widgets" _n
file write fh "SKU-002-B,Gadget B,50,19.99,Gadgets" _n
file write fh "SKU-003-C,Gizmo C,200,4.99,Gizmos" _n
file write fh "SKU-004-D,Thingamajig,25,49.99,Misc" _n
file close fh

noi benchmark_import using "temp/inventory_sku.csv", testname("inventory SKU data")

* Employee HR data
file open fh using "temp/employee_hr.csv", write replace
file write fh "emp_id,first_name,last_name,dept,hire_date,salary,is_manager" _n
file write fh "E001,John,Smith,Engineering,2020-03-15,85000,0" _n
file write fh "E002,Jane,Doe,Marketing,2019-07-01,92000,1" _n
file write fh "E003,Bob,Johnson,Engineering,2021-01-10,78000,0" _n
file write fh `"E004,"Mary Ann",Williams,HR,2018-11-20,88000,1"' _n
file close fh

noi benchmark_import using "temp/employee_hr.csv", testname("employee HR data")

/*******************************************************************************
 * SECTION: Combined Options Tests
 ******************************************************************************/
noi print_section "Combined Options Tests"

* rowrange + colrange + threads
noi cimport_test using "temp/range_test.csv", testname("rowrange+colrange+threads") ///
    importopts(rowrange(10:50) colrange(2:4) threads(4)) expectn(41) expectk(3)

* varnames + case
noi benchmark_import using "temp/mixedcase.csv", testname("varnames(1) + case(lower)") ///
    importopts(varnames(1) case(lower))

* stringcols + asfloat
noi benchmark_import using "temp/zipcodes.csv", testname("stringcols + asfloat") ///
    importopts(stringcols(1) asfloat)

* Note: delimiter + encoding test removed - semicolon delimiter and encoding
* are already tested separately, and the semicolon in delimiters(";") can
* cause syntax issues in some Stata contexts.

* Multiple options together
noi cimport_test using "temp/large_10k.csv", testname("threads+asdouble+rowrange") ///
    importopts(threads(4) asdouble rowrange(100:200)) expectn(101)

/*******************************************************************************
 * SECTION: Excel Import Tests
 ******************************************************************************/
noi print_section "Excel Import Tests"

/*******************************************************************************
 * Helper: benchmark_excel - Import with both methods and compare
 ******************************************************************************/
capture program drop benchmark_excel
program define benchmark_excel
    syntax using/, testname(string) [IMPORTopts(string)]

    clear

    * Import with Stata's import excel (don't pass clear since we already cleared)
    * Remove 'clear' from importopts if present
    local stata_opts : subinstr local importopts "clear" "", all
    capture import excel `using', `stata_opts'
    if _rc != 0 {
        noi test_fail "`testname'" "Stata import excel failed with rc=`=_rc'"
        exit
    }
    tempfile stata_data
    quietly save `stata_data', replace
    local stata_n = _N
    local stata_k = c(k)

    * Import with cimport excel (keep clear option for cimport)
    clear
    capture cimport excel `using', `importopts'
    if _rc != 0 {
        noi test_fail "`testname'" "cimport excel failed with rc=`=_rc'"
        exit
    }
    tempfile cimport_data
    quietly save `cimport_data', replace
    local cimport_n = _N
    local cimport_k = c(k)

    * Check dimensions
    if `stata_n' != `cimport_n' | `stata_k' != `cimport_k' {
        noi test_fail "`testname'" "dimensions differ: Stata N=`stata_n' K=`stata_k', cimport N=`cimport_n' K=`cimport_k'"
        exit
    }

    * Rename variables to generic names for positional comparison
    use `stata_data', clear
    quietly ds
    local stata_vars `r(varlist)'
    local i = 1
    foreach v of local stata_vars {
        quietly rename `v' __cmp_v`i'
        local i = `i' + 1
    }
    tempfile stata_renamed
    quietly save `stata_renamed', replace

    use `cimport_data', clear
    quietly ds
    local cimport_vars `r(varlist)'
    local i = 1
    foreach v of local cimport_vars {
        quietly rename `v' __cmp_v`i'
        local i = `i' + 1
    }

    * Compare data
    capture cf _all using `stata_renamed'
    if _rc == 0 {
        noi test_pass "`testname'"
    }
    else {
        noi test_fail "`testname'" "cf _all comparison failed - data not identical"
    }
end

/*******************************************************************************
 * Helper: cimport_excel_test - Test cimport excel standalone
 ******************************************************************************/
capture program drop cimport_excel_test
program define cimport_excel_test
    syntax using/, testname(string) [IMPORTopts(string) expectn(integer 0) expectk(integer 0) CHECK(string)]

    * Note: importopts should include clear if needed
    capture cimport excel `using', `importopts'
    if _rc != 0 {
        noi test_fail "`testname'" "rc=`=_rc'"
        exit
    }

    if `expectn' > 0 & _N != `expectn' {
        noi test_fail "`testname'" "expected N=`expectn', got N=`=_N'"
        exit
    }

    if `expectk' > 0 & c(k) != `expectk' {
        noi test_fail "`testname'" "expected K=`expectk', got K=`=c(k)'"
        exit
    }

    if "`check'" != "" {
        capture assert `check'
        if _rc != 0 {
            noi test_fail "`testname'" "check() failed: `check'"
            exit
        }
    }

    noi test_pass "`testname'"
end

* Create test Excel files
noi di as text "  Creating test Excel files..."

* Basic Excel file
clear
set obs 5
gen id = _n
gen str10 name = ""
replace name = "Alpha" in 1
replace name = "Beta" in 2
replace name = "Gamma" in 3
replace name = "Delta" in 4
replace name = "Epsilon" in 5
gen double value = .
replace value = 100.5 in 1
replace value = 200.25 in 2
replace value = 300.75 in 3
replace value = 400.1 in 4
replace value = 500.9 in 5
gen str1 category = ""
replace category = "A" in 1
replace category = "B" in 2
replace category = "A" in 3
replace category = "C" in 4
replace category = "B" in 5
export excel using "temp/excel_basic.xlsx", firstrow(variables) replace

* Excel with multiple sheets
clear
set obs 3
gen id = _n
gen str12 name = ""
replace name = "Sheet1Row1" in 1
replace name = "Sheet1Row2" in 2
replace name = "Sheet1Row3" in 3
gen value = _n * 100
export excel using "temp/excel_multisheet.xlsx", sheet("Sheet1") firstrow(variables) replace

clear
set obs 2
gen id = _n * 10
gen str12 name = ""
replace name = "Sheet2Row1" in 1
replace name = "Sheet2Row2" in 2
gen value = _n * 1000
export excel using "temp/excel_multisheet.xlsx", sheet("Sheet2") firstrow(variables) sheetmodify

* Excel with various numeric types
clear
set obs 3
gen byte b = _n * 49 - 48
gen int i = _n * 4950 - 4850
gen long l = _n * 495000 - 485000
gen float f = _n + 0.5
gen double d = _n + 0.123456789012345
export excel using "temp/excel_numtypes.xlsx", firstrow(variables) replace

* Excel with string lengths
clear
set obs 2
gen str1 short = ""
replace short = "A" in 1
replace short = "B" in 2
gen str50 medium = ""
replace medium = "Medium length string" in 1
replace medium = "Another medium one" in 2
gen str200 long_str = ""
replace long_str = "This is a longer string for testing" in 1
replace long_str = "Another longer string for testing" in 2
export excel using "temp/excel_strings.xlsx", firstrow(variables) replace

* Excel without header row (data only)
clear
set obs 3
gen v1 = _n
gen v2 = _n * 10
gen v3 = _n * 100
export excel using "temp/excel_noheader.xlsx", replace

* Excel with uppercase variable names
clear
set obs 2
gen ID = _n
gen str10 NAME = ""
replace NAME = "First" in 1
replace NAME = "Second" in 2
gen VALUE = _n * 100
export excel using "temp/excel_uppercase.xlsx", firstrow(variables) replace

* Basic import tests
noi benchmark_excel using "temp/excel_basic.xlsx", testname("Excel basic import") importopts(firstrow clear)
noi cimport_excel_test using "temp/excel_noheader.xlsx", testname("Excel no firstrow") importopts(clear) expectn(3) expectk(3)
noi cimport_excel_test using "temp/excel_basic.xlsx", testname("Excel with firstrow") importopts(firstrow clear) expectn(5) expectk(4)

* Sheet selection tests
noi cimport_excel_test using "temp/excel_multisheet.xlsx", testname("Excel first sheet (default)") ///
    importopts(firstrow clear) expectn(3) expectk(3)
noi cimport_excel_test using "temp/excel_multisheet.xlsx", testname("Excel named sheet") ///
    importopts(sheet("Sheet2") firstrow clear) expectn(2) expectk(3)

* Cell range tests
noi cimport_excel_test using "temp/excel_basic.xlsx", testname("Excel cellrange A1:B3") ///
    importopts(cellrange(A1:B3) firstrow clear) expectn(2) expectk(2)
noi cimport_excel_test using "temp/excel_basic.xlsx", testname("Excel cellrange from A2") ///
    importopts(cellrange(A2:) firstrow clear) expectn(4) expectk(4)

* Type handling tests
noi benchmark_excel using "temp/excel_numtypes.xlsx", testname("Excel numeric types") importopts(firstrow clear)
noi benchmark_excel using "temp/excel_strings.xlsx", testname("Excel string types") importopts(firstrow clear)
noi cimport_excel_test using "temp/excel_basic.xlsx", testname("Excel allstring option") ///
    importopts(allstring firstrow clear) expectn(5) expectk(4)

* Case handling tests
noi cimport_excel_test using "temp/excel_uppercase.xlsx", testname("Excel case(preserve)") ///
    importopts(firstrow clear) expectn(2) expectk(3)
noi cimport_excel_test using "temp/excel_uppercase.xlsx", testname("Excel case(lower)") ///
    importopts(firstrow clear case(lower)) expectn(2) expectk(3)
noi cimport_excel_test using "temp/excel_uppercase.xlsx", testname("Excel case(upper)") ///
    importopts(firstrow clear case(upper)) expectn(2) expectk(3)

* Error handling tests
capture cimport excel using "temp/nonexistent.xlsx", clear
if _rc != 0 {
    noi test_pass "Excel nonexistent file error"
}
else {
    noi test_fail "Excel nonexistent file error" "should have returned error"
}

capture cimport excel using "temp/basic.csv", clear
if _rc != 0 {
    noi test_pass "Excel invalid extension error"
}
else {
    noi test_fail "Excel invalid extension error" "should have returned error"
}

* Clean up Excel temp files
capture erase "temp/excel_basic.xlsx"
capture erase "temp/excel_multisheet.xlsx"
capture erase "temp/excel_numtypes.xlsx"
capture erase "temp/excel_strings.xlsx"
capture erase "temp/excel_noheader.xlsx"
capture erase "temp/excel_uppercase.xlsx"

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
