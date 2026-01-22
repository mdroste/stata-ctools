/*******************************************************************************
 * validate_cexport.do
 *
 * Comprehensive validation tests for cexport delimited vs export delimited
 * Tests all export options: delimiter, novarnames, quote, nolabel, if/in
 *
 * VERIFICATION: Exports both CSVs, imports them back into Stata, and compares:
 *   - First attempts cf _all for byte-for-byte comparison
 *   - Falls back to tolerance-based comparison (1e-6) for numeric variables
 *     to handle floating point precision differences in CSV formatting
 *   - String variables must match exactly
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

di as text ""
di as text "======================================================================"
di as text "              CEXPORT DELIMITED VALIDATION TEST SUITE"
di as text "======================================================================"

capture mkdir "temp"

/*******************************************************************************
 * Helper: benchmark_export - Export with both methods, import back, compare
 * Required: testname() - the test criterion/name
 * Optional: exportopts() - options for export commands
 *           importopts() - options for re-importing CSVs
 *           ifcond() - if condition
 *           incond() - in condition
 ******************************************************************************/
capture program drop benchmark_export
program define benchmark_export
    syntax [varlist], testname(string) [EXPORTopts(string) IMPORTopts(string) IFcond(string) INcond(string)]

    * Build if/in conditions
    local ifin ""
    if "`ifcond'" != "" local ifin "`ifin' if `ifcond'"
    if "`incond'" != "" local ifin "`ifin' in `incond'"

    tempfile stata_csv cexport_csv

    * Export with Stata's export delimited
    if "`varlist'" != "" {
        export delimited `varlist' using "`stata_csv'" `ifin', `exportopts' replace
    }
    else {
        export delimited using "`stata_csv'" `ifin', `exportopts' replace
    }

    * Export with cexport delimited
    if "`varlist'" != "" {
        cexport delimited `varlist' using "`cexport_csv'" `ifin', `exportopts' replace
    }
    else {
        cexport delimited using "`cexport_csv'" `ifin', `exportopts' replace
    }

    * Import both CSVs back
    preserve

    import delimited using "`stata_csv'", `importopts' clear
    tempfile stata_data
    quietly save `stata_data', replace
    local stata_n = _N
    local stata_k = c(k)

    import delimited using "`cexport_csv'", `importopts' clear
    tempfile cexport_data
    quietly save `cexport_data', replace
    local cexport_n = _N
    local cexport_k = c(k)

    * Check dimensions first
    if `stata_n' != `cexport_n' | `stata_k' != `cexport_k' {
        restore
        noi test_fail "`testname'" "dimensions differ: Stata N=`stata_n' K=`stata_k', cexport N=`cexport_n' K=`cexport_k'"
        exit
    }

    * Compare data using cf _all
    use `stata_data', clear
    capture cf _all using `cexport_data'
    local cfrc = _rc

    if `cfrc' == 0 {
        restore
        noi test_pass "`testname'"
        exit
    }

    * cf _all failed - try tolerance-based comparison for numeric variables
    use `stata_data', clear
    ds
    local allvars `r(varlist)'

    * Rename variables from stata data with _stata suffix
    foreach v of local allvars {
        rename `v' `v'_stata
    }
    gen long _row = _n
    tempfile stata_renamed
    quietly save `stata_renamed', replace

    * Load cexport data and rename with _cexport suffix
    use `cexport_data', clear
    foreach v of local allvars {
        capture confirm variable `v'
        if _rc != 0 {
            restore
            noi test_fail "`testname'" "variable `v' missing from cexport data"
            exit
        }
        rename `v' `v'_cexport
    }
    gen long _row = _n
    merge 1:1 _row using `stata_renamed', nogen

    local all_match = 1
    local fail_reason ""
    local tol = 1e-6

    foreach v of local allvars {
        * Check if string or numeric
        capture confirm string variable `v'_stata
        if _rc == 0 {
            * String variable - must be exactly equal
            quietly count if `v'_stata != `v'_cexport
            if r(N) > 0 {
                local all_match = 0
                local fail_reason "string variable `v' has `r(N)' mismatches"
                continue, break
            }
        }
        else {
            * Numeric variable - allow tolerance
            quietly count if abs(`v'_stata - `v'_cexport) > `tol' & !missing(`v'_stata) & !missing(`v'_cexport)
            if r(N) > 0 {
                local all_match = 0
                local fail_reason "numeric variable `v' has `r(N)' values differing by >`tol'"
                continue, break
            }
            * Check missing values match
            quietly count if missing(`v'_stata) != missing(`v'_cexport)
            if r(N) > 0 {
                local all_match = 0
                local fail_reason "variable `v' has `r(N)' missing value mismatches"
                continue, break
            }
        }
    }

    restore

    if `all_match' == 1 {
        noi test_pass "`testname'"
    }
    else {
        noi test_fail "`testname'" "`fail_reason'"
    }
end

/*******************************************************************************
 * SECTION 1: Plugin check
 ******************************************************************************/
noi print_section "Plugin Check"

sysuse auto, clear
capture cexport delimited using "temp/test.csv", replace
if _rc != 0 {
    noi test_fail "cexport plugin load" "returned error `=_rc'"
    noi print_summary "cexport"
    exit 1
}
noi test_pass "cexport plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Basic comma-delimited export
 ******************************************************************************/
noi print_section "Basic Comma-Delimited Export"

sysuse auto, clear
benchmark_export, testname("basic export")

/*******************************************************************************
 * SECTION 3: Tab delimiter
 ******************************************************************************/
noi print_section "Tab Delimiter"

sysuse auto, clear
benchmark_export, testname("tab delimiter") exportopts(delimiter(tab)) importopts(delimiters(tab))

/*******************************************************************************
 * SECTION 4: Semicolon delimiter
 ******************************************************************************/
noi print_section "Semicolon Delimiter"

sysuse auto, clear
benchmark_export, testname("semicolon delimiter") exportopts(delimiter(";")) importopts(delimiters(";"))

/*******************************************************************************
 * SECTION 5: novarnames option
 ******************************************************************************/
noi print_section "novarnames Option"

sysuse auto, clear
benchmark_export, testname("novarnames") exportopts(novarnames) importopts(varnames(nonames))

/*******************************************************************************
 * SECTION 6: quote option
 ******************************************************************************/
noi print_section "quote Option"

sysuse auto, clear
benchmark_export, testname("quote option")  exportopts(quote)

/*******************************************************************************
 * SECTION 7: nolabel option
 ******************************************************************************/
noi print_section "nolabel Option"

sysuse auto, clear
benchmark_export, testname("with labels")

sysuse auto, clear
benchmark_export, testname("nolabel") exportopts(nolabel)

/*******************************************************************************
 * SECTION 8: Variable selection
 ******************************************************************************/
noi print_section "Variable Selection"

sysuse auto, clear
benchmark_export make price mpg, testname("variable selection (3 vars)")

/*******************************************************************************
 * SECTION 9: if condition
 ******************************************************************************/
noi print_section "if Condition"

sysuse auto, clear
benchmark_export, testname("if foreign==1") ifcond(foreign == 1)

sysuse auto, clear
benchmark_export, testname("if price>10000") ifcond(price > 10000)

/*******************************************************************************
 * SECTION 10: in condition
 ******************************************************************************/
noi print_section "in Condition"

sysuse auto, clear
benchmark_export, testname("in 1/20") incond(1/20)

sysuse auto, clear
benchmark_export, testname("in 30/50") incond(30/50)

/*******************************************************************************
 * SECTION 11: Combined if and in
 ******************************************************************************/
noi print_section "Combined if and in"

sysuse auto, clear
benchmark_export, testname("if and in combined") ifcond(price > 5000) incond(1/50)

/*******************************************************************************
 * SECTION 12: Census dataset
 ******************************************************************************/
noi print_section "Census Dataset"

sysuse census, clear
benchmark_export, testname("census")

/*******************************************************************************
 * SECTION 13: Large dataset export
 ******************************************************************************/
noi print_section "Large Dataset Export"

clear
set seed 12345
set obs 50000
gen id = _n
gen group = runiformint(1, 100)
gen x = runiform()
gen y = rnormal()
gen str20 label = "item" + string(runiformint(1, 1000))

benchmark_export, testname("large dataset (50K)")

/*******************************************************************************
 * SECTION 14: Panel data (nlswork)
 ******************************************************************************/
noi print_section "Panel Data (nlswork)"

webuse nlswork, clear
keep in 1/10000

benchmark_export, testname("panel data")

/*******************************************************************************
 * SECTION 15: verbose/timeit options
 ******************************************************************************/
noi print_section "verbose/timeit Options"

sysuse auto, clear

capture cexport delimited using "temp/test.csv", verbose replace
if _rc == 0 {
    noi test_pass "verbose option accepted"
}
else {
    noi test_fail "verbose option" "returned error `=_rc'"
}

capture cexport delimited using "temp/test.csv", timeit replace
if _rc == 0 {
    noi test_pass "timeit option accepted"
}
else {
    noi test_fail "timeit option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 16: replace behavior
 ******************************************************************************/
noi print_section "replace Behavior"

sysuse auto, clear
cexport delimited using "temp/repl.csv", replace

capture cexport delimited using "temp/repl.csv"
if _rc != 0 {
    noi test_pass "without replace: fails for existing file"
}
else {
    noi test_fail "without replace" "should fail"
}

capture cexport delimited using "temp/repl.csv", replace
if _rc == 0 {
    noi test_pass "with replace: overwrites file"
}
else {
    noi test_fail "with replace" "failed"
}

/*******************************************************************************
 * SECTION 17: Edge cases
 ******************************************************************************/
noi print_section "Edge Cases"

* Single observation
clear
set obs 1
gen x = 42
gen str5 s = "test"
benchmark_export, testname("single observation")

* Single variable
clear
set obs 100
gen x = runiform()
benchmark_export x, testname("single variable")

* String with commas
clear
input id str40 name
1 "Smith, John"
2 "Doe, Jane"
3 "Regular Name"
end
benchmark_export, testname("strings with commas")

/*******************************************************************************
 * SECTION 18: Missing values
 ******************************************************************************/
noi print_section "Missing Values"

clear
set obs 10
gen x = _n
replace x = . if mod(_n, 3) == 0
gen str10 s = "val" + string(_n)
replace s = "" if mod(_n, 4) == 0

benchmark_export, testname("missing values")

/*******************************************************************************
 * SECTION 19: Numeric precision
 ******************************************************************************/
noi print_section "Numeric Precision"

clear
set obs 100
gen double precise = runiform() * 1e10
gen float approx = runiform()
gen long bigint = runiformint(1, 2147483647)
gen byte tinyint = runiformint(-127, 127)

benchmark_export, testname("numeric precision")

/*******************************************************************************
 * SECTION 20: Special characters in strings
 ******************************************************************************/
noi print_section "Special Characters"

clear
input id str50 text
1 "Normal text"
2 "Text with tabs	here"
3 "Text with newlines here"
4 "Text with ""quotes"""
end

benchmark_export, testname("special characters")

/*******************************************************************************
 * SECTION: Large Dataset Tests
 ******************************************************************************/
noi print_section "Large Datasets"

* 10K rows export
clear
set obs 10000
gen id = _n
gen value = runiform()
gen str20 name = "item_" + string(_n)

capture cexport delimited using "temp/large_10k.csv", replace
if _rc == 0 {
    noi test_pass "10K rows export"
}
else {
    noi test_fail "10K rows" "rc=`=_rc'"
}

* Compare with export delimited
export delimited using "temp/large_10k_stata.csv", replace
cexport delimited using "temp/large_10k_cexport.csv", replace

file open f1 using "temp/large_10k_stata.csv", read
file read f1 line1_stata
file close f1

file open f2 using "temp/large_10k_cexport.csv", read
file read f2 line1_cexport
file close f2

if "`line1_stata'" == "`line1_cexport'" {
    noi test_pass "large export matches Stata header"
}
else {
    noi test_fail "large match" "headers differ"
}

* 50K rows
clear
set obs 50000
gen id = _n
gen x = runiform()
gen y = runiform()

capture cexport delimited using "temp/large_50k.csv", replace
if _rc == 0 {
    noi test_pass "50K rows export"
}
else {
    noi test_fail "50K rows" "rc=`=_rc'"
}

* Many columns (20)
clear
set obs 1000
forvalues i = 1/20 {
    gen var`i' = runiform()
}

capture cexport delimited using "temp/many_cols.csv", replace
if _rc == 0 {
    noi test_pass "20 columns export"
}
else {
    noi test_fail "20 columns" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION: Pathological Data
 ******************************************************************************/
noi print_section "Pathological Data"

* All missing values
clear
set obs 100
gen x = .
gen y = .
gen z = .

capture cexport delimited using "temp/all_missing.csv", replace
if _rc == 0 {
    noi test_pass "all missing values"
}
else {
    noi test_fail "all missing" "rc=`=_rc'"
}

* Empty strings
clear
set obs 50
gen str20 name = ""
gen value = _n

capture cexport delimited using "temp/empty_strings.csv", replace
if _rc == 0 {
    noi test_pass "empty strings"
}
else {
    noi test_fail "empty strings" "rc=`=_rc'"
}

* Single column
clear
set obs 100
gen only_col = runiform()

capture cexport delimited using "temp/single_col.csv", replace
if _rc == 0 {
    noi test_pass "single column"
}
else {
    noi test_fail "single column" "rc=`=_rc'"
}

* Single row
clear
set obs 1
gen a = 1
gen b = 2
gen c = 3

capture cexport delimited using "temp/single_row.csv", replace
if _rc == 0 {
    noi test_pass "single row"
}
else {
    noi test_fail "single row" "rc=`=_rc'"
}

* Very long strings
clear
set obs 10
gen str244 long_text = "a" * 200

capture cexport delimited using "temp/long_strings.csv", replace
if _rc == 0 {
    noi test_pass "long strings (200 chars)"
}
else {
    noi test_fail "long strings" "rc=`=_rc'"
}

* Extreme numeric values
clear
set obs 5
gen double extreme = .
replace extreme = 1e308 in 1
replace extreme = -1e308 in 2
replace extreme = 1e-308 in 3
replace extreme = 0 in 4
replace extreme = . in 5

capture cexport delimited using "temp/extreme_nums.csv", replace
if _rc == 0 {
    noi test_pass "extreme numeric values"
}
else {
    noi test_fail "extreme nums" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION: Real-World Datasets
 ******************************************************************************/
noi print_section "Real-World Datasets"

* auto dataset
sysuse auto, clear
benchmark_export, testname("auto dataset full")

* census dataset
sysuse census, clear
benchmark_export, testname("census dataset")

* nlswork dataset
webuse nlswork, clear
capture cexport delimited using "temp/nlswork.csv", replace
if _rc == 0 {
    noi test_pass "nlswork dataset export"
}
else {
    noi test_fail "nlswork" "rc=`=_rc'"
}

* lifeexp dataset
webuse lifeexp, clear
benchmark_export, testname("lifeexp dataset")

/*******************************************************************************
 * SECTION: Comparison Tests
 ******************************************************************************/
noi print_section "Comparison with export delimited"

* Basic comparison
sysuse auto, clear
export delimited using "temp/compare_stata.csv", replace
cexport delimited using "temp/compare_cexport.csv", replace

* Compare file sizes (should be similar)
capture file open f1 using "temp/compare_stata.csv", read
capture file open f2 using "temp/compare_cexport.csv", read
* Just check both files exist and can be opened
if _rc == 0 {
    noi test_pass "comparison files created"
}
else {
    noi test_fail "comparison" "file error"
}
capture file close f1
capture file close f2

* Round-trip test: export then import
clear
set obs 100
gen id = _n
gen value = runiform()
gen str20 name = "test_" + string(_n)
tempfile original
save `original'

cexport delimited using "temp/roundtrip.csv", replace
import delimited using "temp/roundtrip.csv", clear
local imported_N = _N

if `imported_N' == 100 {
    noi test_pass "round-trip export/import"
}
else {
    noi test_fail "round-trip" "wrong N"
}

/*******************************************************************************
 * SECTION: Edge Cases
 ******************************************************************************/
noi print_section "Edge Cases"

* Variable names with underscores
clear
set obs 10
gen my_var_name = _n
gen another_var = runiform()

capture cexport delimited using "temp/underscores.csv", replace
if _rc == 0 {
    noi test_pass "variable names with underscores"
}
else {
    noi test_fail "underscores" "rc=`=_rc'"
}

* Numeric-looking strings
clear
set obs 5
gen str10 numstr = string(_n * 100)

capture cexport delimited using "temp/numstr.csv", replace
if _rc == 0 {
    noi test_pass "numeric-looking strings"
}
else {
    noi test_fail "numstr" "rc=`=_rc'"
}

* Mixed types
clear
set obs 20
gen byte b = mod(_n, 128)
gen int i = _n * 100
gen long l = _n * 10000
gen float f = runiform()
gen double d = runiform() * 1e10
gen str20 s = "text_" + string(_n)

capture cexport delimited using "temp/mixed_types.csv", replace
if _rc == 0 {
    noi test_pass "mixed variable types"
}
else {
    noi test_fail "mixed types" "rc=`=_rc'"
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

noi print_summary "cexport"

if $TESTS_FAILED > 0 {
    exit 1
}

}
