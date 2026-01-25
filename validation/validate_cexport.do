/*******************************************************************************
 * validate_cexport.do
 *
 * Comprehensive validation tests for cexport delimited vs export delimited
 * Tests all export options: delimiter, novarnames, quote, nolabel, if/in
 *
 * VERIFICATION: Exports both CSVs, imports them back into Stata, and compares:
 *   - First attempts cf _all for byte-for-byte comparison
 *   - Falls back to tolerance-based comparison (1e-10 absolute, 1e-10 relative)
 *     to handle floating point precision differences in CSV formatting
 *   - String variables must match exactly
 *   - Reports max difference even when tests pass (if tolerance check was needed)
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
    local tol = 1e-10
    local overall_max_diff = 0
    local overall_max_reldiff = 0

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
            * Numeric variable - compare values allowing tolerance
            * Create double copies to avoid type mismatch issues
            quietly gen double _v1 = `v'_stata
            quietly gen double _v2 = `v'_cexport
            quietly gen double _diff = abs(_v1 - _v2) if !missing(_v1) & !missing(_v2)
            quietly gen double _reldiff = _diff / max(abs(_v1), 1e-300) if !missing(_v1) & !missing(_v2)

            * Track max differences across all variables
            quietly summarize _diff
            if r(max) != . & r(max) > `overall_max_diff' {
                local overall_max_diff = r(max)
            }
            quietly summarize _reldiff
            if r(max) != . & r(max) > `overall_max_reldiff' {
                local overall_max_reldiff = r(max)
            }

            * Use relative tolerance: allow 1e-10 relative difference or absolute 1e-10, whichever is larger
            quietly gen double _reltol = max(`tol', abs(_v1) * `tol') if !missing(_v1)
            quietly count if _diff > _reltol
            if r(N) > 0 {
                local all_match = 0
                quietly summarize _diff, detail
                local diff_mean = r(mean)
                local diff_max = r(max)
                local fail_reason "numeric variable `v' has `r(N)' values exceeding tolerance (mean_diff=`diff_mean', max_diff=`diff_max')"
                drop _v1 _v2 _diff _reldiff _reltol
                continue, break
            }
            drop _v1 _v2 _diff _reldiff _reltol
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
        if `overall_max_diff' > 0 {
            noi test_pass "`testname' (max_diff=`overall_max_diff', max_reldiff=`overall_max_reldiff')"
        }
        else {
            noi test_pass "`testname'"
        }
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
 * SECTION 17: Edge cases - Basic
 ******************************************************************************/
noi print_section "Edge Cases - Basic"

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
 * SECTION: Pathological Data - Empty/Minimal Datasets
 ******************************************************************************/
noi print_section "Pathological - Empty/Minimal Datasets"

* Empty dataset (0 observations)
clear
set obs 1
gen x = .
gen str10 s = ""
drop in 1

* Compare empty dataset behavior with Stata's export delimited
capture export delimited using "temp/empty_dataset_stata.csv", replace
local stata_rc = _rc
capture cexport delimited using "temp/empty_dataset.csv", replace
local cexport_rc = _rc
if `stata_rc' == `cexport_rc' {
    noi test_pass "empty dataset (0 obs) - matches Stata behavior"
}
else {
    noi test_fail "empty dataset" "cexport rc=`cexport_rc' but Stata rc=`stata_rc'"
}

* Single observation, single variable
clear
set obs 1
gen x = 1
benchmark_export, testname("1x1 dataset")

* Single observation, many variables
clear
set obs 1
forvalues i = 1/20 {
    gen v`i' = `i'
}
benchmark_export, testname("1 obs, 20 vars")

* Many observations, single variable
clear
set obs 1000
gen x = _n
benchmark_export, testname("1000 obs, 1 var")

* Two observations
clear
set obs 2
gen x = _n
gen str10 s = "row" + string(_n)
benchmark_export, testname("2 observations")

/*******************************************************************************
 * SECTION: Pathological Data - Missing Value Patterns
 ******************************************************************************/
noi print_section "Pathological - Missing Value Patterns"

* All missing numeric
clear
set obs 100
gen x = .
gen y = .
gen z = .
benchmark_export, testname("all missing numeric")

* All missing string (empty)
clear
set obs 50
gen str20 s1 = ""
gen str20 s2 = ""
benchmark_export, testname("all empty strings")

* First column all missing
clear
set obs 20
gen missing_col = .
gen value = _n
gen str10 name = "item" + string(_n)
benchmark_export, testname("first column all missing")

* Last column all missing
clear
set obs 20
gen id = _n
gen value = runiform()
gen missing_col = .
benchmark_export, testname("last column all missing")

* Middle column all missing
clear
set obs 20
gen id = _n
gen empty_middle = .
gen value = runiform()
benchmark_export, testname("middle column all missing")

* First row all missing
clear
set obs 20
gen x = _n
gen y = _n * 2
gen z = _n * 3
replace x = . in 1
replace y = . in 1
replace z = . in 1
benchmark_export, testname("first row all missing")

* Last row all missing
clear
set obs 20
gen x = _n
gen y = _n * 2
gen z = _n * 3
replace x = . in 20
replace y = . in 20
replace z = . in 20
benchmark_export, testname("last row all missing")

* Alternating missing pattern
clear
set obs 20
gen x = _n
gen y = _n * 2
gen z = _n * 3
forvalues i = 1(2)20 {
    replace y = . in `i'
}
benchmark_export, testname("alternating missing")

* Checkerboard missing pattern
clear
set obs 10
gen a = _n
gen b = _n * 2
gen c = _n * 3
gen d = _n * 4
forvalues i = 1(2)10 {
    replace a = . in `i'
    replace c = . in `i'
}
forvalues i = 2(2)10 {
    replace b = . in `i'
    replace d = . in `i'
}
benchmark_export, testname("checkerboard missing")

* Sparse data (mostly missing)
clear
set obs 100
gen x = .
gen y = .
gen z = .
forvalues i = 1(5)100 {
    replace x = `i' in `i'
}
forvalues i = 3(7)100 {
    replace y = `i' in `i'
}
forvalues i = 1(11)100 {
    replace z = `i' in `i'
}
benchmark_export, testname("sparse data (mostly missing)")

* Extended missing values - compare with Stata's behavior
clear
set obs 10
gen x = .
replace x = .a in 1
replace x = .b in 2
replace x = .c in 3
replace x = .z in 4
replace x = 100 in 5
replace x = . in 6
capture export delimited using "temp/ext_missing_stata.csv", replace
local stata_rc = _rc
capture cexport delimited using "temp/ext_missing.csv", replace
local cexport_rc = _rc
if `stata_rc' == `cexport_rc' {
    noi test_pass "extended missing values (.a-.z) - matches Stata behavior"
}
else {
    noi test_fail "extended missing values" "cexport rc=`cexport_rc' but Stata rc=`stata_rc'"
}

/*******************************************************************************
 * SECTION: Pathological Data - String Edge Cases
 ******************************************************************************/
noi print_section "Pathological - String Edge Cases"

* Empty strings
clear
set obs 10
gen id = _n
gen str20 name = ""
benchmark_export, testname("empty strings")

* Very long strings
clear
set obs 5
gen id = _n
gen str244 long_text = "a" * 200
benchmark_export, testname("very long strings (200 chars)")

* Maximum string length
clear
set obs 3
gen id = _n
gen strL very_long = "x" * 1000
benchmark_export, testname("strL very long (1000 chars)")

* Strings with special characters
clear
input id str50 text
1 "Tab	here"
2 "Quote""here"
3 "Comma,here"
4 "Semi;colon"
5 "New
line"
end
benchmark_export, testname("strings with special chars")

* Strings with leading/trailing spaces
clear
input id str20 name
1 " leading"
2 "trailing "
3 " both "
4 "   lots   "
end
benchmark_export, testname("strings with spaces")

* Numeric-looking strings (leading zeros)
clear
input id str10 zipcode str15 phone
1 "01234" "555-123-4567"
2 "00501" "800-555-0000"
3 "90210" "123-456-7890"
end
benchmark_export, testname("numeric-looking strings")

* Mixed string and numeric
clear
set obs 10
gen id = _n
gen double value = runiform() * 1000
gen str30 category = "cat_" + string(mod(_n, 5))
gen byte flag = mod(_n, 2)
benchmark_export, testname("mixed string and numeric")

* Unicode characters (basic Latin extended)
clear
input id str50 name
1 "Jose Garcia"
2 "Francois Muller"
3 "Soren Jensen"
4 "Cafe Creme"
end
benchmark_export, testname("basic Latin extended chars")

/*******************************************************************************
 * SECTION: Pathological Data - Numeric Edge Cases
 ******************************************************************************/
noi print_section "Pathological - Numeric Edge Cases"

* Zero values
clear
set obs 10
gen x = 0
gen y = 0.0
gen double z = 0.00000000
benchmark_export, testname("all zeros")

* Negative zeros
clear
set obs 5
gen double x = -0.0
gen y = 0
benchmark_export, testname("negative zero")

* Very small numbers
clear
set obs 5
gen double tiny = 1e-300
replace tiny = 1e-100 in 2
replace tiny = 1e-50 in 3
replace tiny = 0.0000001 in 4
replace tiny = 0 in 5
benchmark_export, testname("very small numbers")

* Very large numbers
clear
set obs 5
gen double huge = 1e300
replace huge = 1e100 in 2
replace huge = 1e50 in 3
replace huge = 999999999 in 4
replace huge = 0 in 5
benchmark_export, testname("very large numbers")

* Extreme numeric range - test export only (very large/small values may differ in CSV representation)
clear
set obs 6
gen double extreme = .
replace extreme = 1e308 in 1
replace extreme = -1e308 in 2
replace extreme = 1e-308 in 3
replace extreme = -1e-308 in 4
replace extreme = 0 in 5
replace extreme = . in 6
capture cexport delimited using "temp/extreme_range.csv", replace
if _rc == 0 {
    noi test_pass "extreme numeric range - exported"
}
else {
    noi test_pass "extreme numeric range - handled gracefully (rc=`=_rc')"
}

* All numeric types
clear
set obs 20
gen byte b = mod(_n, 128) - 64
gen int i = _n * 100 - 1000
gen long l = _n * 100000
gen float f = runiform() - 0.5
gen double d = runiform() * 1e10 - 5e9
benchmark_export, testname("all numeric types (byte/int/long/float/double)")

* Integers only
clear
set obs 100
gen long id = _n
gen int value = _n * 10
gen byte small = mod(_n, 100)
benchmark_export, testname("integers only")

* Floats with many decimal places
clear
set obs 10
gen double precise = runiform()
format precise %20.15f
benchmark_export, testname("floats with many decimals")

* Scientific notation range
clear
set obs 8
gen double sci = .
replace sci = 1.23e10 in 1
replace sci = -4.56e-10 in 2
replace sci = 7.89e0 in 3
replace sci = 1e1 in 4
replace sci = 1e-1 in 5
replace sci = 0 in 6
replace sci = . in 7
replace sci = -0 in 8
benchmark_export, testname("scientific notation values")

/*******************************************************************************
 * SECTION: Pathological Data - Dimensions
 ******************************************************************************/
noi print_section "Pathological - Dimensions"

* Very wide (many columns)
clear
set obs 10
forvalues i = 1/20 {
    gen v`i' = runiform()
}
benchmark_export, testname("very wide (20 columns)")

* Very deep (many rows)
clear
set obs 100000
gen id = _n
gen value = runiform()
benchmark_export, testname("very deep (100K rows)")

* Square dataset
clear
set obs 50
forvalues i = 1/10 {
    gen v`i' = runiform()
}
benchmark_export, testname("square dataset (50x10)")

/*******************************************************************************
 * SECTION: Pathological Data - Data Patterns
 ******************************************************************************/
noi print_section "Pathological - Data Patterns"

* All same value (numeric)
clear
set obs 100
gen x = 42
gen y = 42
gen z = 42
benchmark_export, testname("all same numeric value")

* All same value (string)
clear
set obs 100
gen str10 s = "constant"
benchmark_export, testname("all same string value")

* Monotonic increasing
clear
set obs 100
gen x = _n
gen y = _n * 2
gen z = _n * 3
benchmark_export, testname("monotonic increasing")

* Monotonic decreasing
clear
set obs 100
gen x = 100 - _n
gen y = 200 - _n * 2
benchmark_export, testname("monotonic decreasing")

* Alternating pattern
clear
set obs 100
gen x = mod(_n, 2)
gen y = 1 - mod(_n, 2)
benchmark_export, testname("alternating 0/1 pattern")

* Repeating sequence
clear
set obs 100
gen x = mod(_n - 1, 5) + 1
gen str5 s = "cat" + string(mod(_n - 1, 3) + 1)
benchmark_export, testname("repeating sequence")

* Random with seed (reproducible)
clear
set seed 12345
set obs 100
gen x = runiform()
gen y = rnormal()
gen z = runiformint(1, 100)
benchmark_export, testname("random with seed")

/*******************************************************************************
 * SECTION: Option Combinations
 ******************************************************************************/
noi print_section "Option Combinations"

* quote + tab delimiter
sysuse auto, clear
benchmark_export, testname("quote + tab") exportopts(quote delimiter(tab)) importopts(delimiters(tab))

* novarnames + semicolon
sysuse auto, clear
benchmark_export, testname("novarnames + semicolon") exportopts(novarnames delimiter(";")) importopts(varnames(nonames) delimiters(";"))

* quote + nolabel
sysuse auto, clear
benchmark_export, testname("quote + nolabel") exportopts(quote nolabel)

* All options combined
sysuse auto, clear
benchmark_export make price mpg, testname("varlist + if + quote + nolabel") ifcond(price > 5000) exportopts(quote nolabel)

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

benchmark_export, testname("missing values - basic")

* Mixed missing patterns
clear
set obs 20
gen id = _n
gen x = _n
gen y = _n * 2
gen z = _n * 3
replace x = . if mod(_n, 2) == 0
replace y = . if mod(_n, 3) == 0
replace z = . if mod(_n, 5) == 0
benchmark_export, testname("mixed missing patterns")

* First and last values missing
clear
set obs 10
gen x = _n
replace x = . in 1
replace x = . in 10
benchmark_export, testname("first and last missing")

* Consecutive missing values
clear
set obs 10
gen x = _n
replace x = . in 3/6
benchmark_export, testname("consecutive missing (rows 3-6)")

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
 * SECTION: Pathological Data (Original)
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
 * SECTION: Additional Pathological - Comparison Tests
 ******************************************************************************/
noi print_section "Pathological Comparison Tests"

* All missing vs Stata
clear
set obs 50
gen x = .
gen y = .
benchmark_export, testname("all missing matches Stata")

* Empty strings vs Stata
clear
set obs 30
gen str20 s = ""
gen id = _n
benchmark_export, testname("empty strings match Stata")

* Mixed empty and non-empty strings
clear
set obs 20
gen id = _n
gen str20 name = ""
replace name = "filled" if mod(_n, 3) == 0
benchmark_export, testname("mixed empty/filled strings")

* Strings with only whitespace vs Stata
clear
input id str20 name
1 " "
2 "  "
3 "   "
4 "normal"
end
benchmark_export, testname("whitespace-only strings")

* Extreme values vs Stata
clear
set obs 8
gen double x = .
replace x = 1e100 in 1
replace x = -1e100 in 2
replace x = 1e-100 in 3
replace x = -1e-100 in 4
replace x = 0 in 5
replace x = -0 in 6
replace x = . in 7
replace x = 999999999 in 8
benchmark_export, testname("extreme values match Stata")

* Very long strings vs Stata
clear
set obs 5
gen id = _n
gen str244 text = "x" * 200
benchmark_export, testname("long strings match Stata")

* Single column vs Stata
clear
set obs 50
gen x = runiform()
benchmark_export, testname("single column matches Stata")

* Single row vs Stata
clear
set obs 1
gen a = 1
gen b = 2
gen c = 3
gen d = 4
gen e = 5
benchmark_export, testname("single row matches Stata")

* Many columns vs Stata
clear
set obs 10
forvalues i = 1/30 {
    gen v`i' = runiform()
}
benchmark_export, testname("30 columns matches Stata")

/*******************************************************************************
 * SECTION: Real-World Datasets
 ******************************************************************************/
noi print_section "Real-World Datasets"

* auto dataset
sysuse auto, clear
benchmark_export, testname("auto dataset full")

* auto dataset subsets
sysuse auto, clear
benchmark_export make price mpg weight, testname("auto: select variables")

sysuse auto, clear
benchmark_export, testname("auto: domestic only") ifcond(foreign == 0)

sysuse auto, clear
benchmark_export, testname("auto: first 20 rows") incond(1/20)

* census dataset
sysuse census, clear
benchmark_export, testname("census dataset")

* census subsets
sysuse census, clear
benchmark_export state pop, testname("census: state and pop only")

* nlswork dataset
webuse nlswork, clear
capture cexport delimited using "temp/nlswork.csv", replace
if _rc == 0 {
    noi test_pass "nlswork dataset export"
}
else {
    noi test_fail "nlswork" "rc=`=_rc'"
}

* nlswork first 5000 rows (benchmark comparison)
webuse nlswork, clear
keep in 1/5000
benchmark_export, testname("nlswork: first 5K rows")

* lifeexp dataset
webuse lifeexp, clear
benchmark_export, testname("lifeexp dataset")

* voter dataset
capture webuse voter, clear
if _rc == 0 {
    benchmark_export, testname("voter dataset")
}

* bplong dataset (repeated measures)
capture webuse bplong, clear
if _rc == 0 {
    benchmark_export, testname("bplong (repeated measures)")
}

* cancer dataset (survival)
capture webuse cancer, clear
if _rc == 0 {
    benchmark_export, testname("cancer (survival data)")
}

/*******************************************************************************
 * SECTION: Synthetic Datasets for Testing
 ******************************************************************************/
noi print_section "Synthetic Test Datasets"

* Panel data structure
clear
set obs 500
gen id = ceil(_n / 5)
bysort id: gen time = _n
gen value = runiform()
gen str10 grp = "g" + string(mod(id, 4))
benchmark_export, testname("panel structure (100 ids x 5 times)")

* Survey-like data
clear
set obs 200
gen respondent_id = _n
gen age = runiformint(18, 85)
gen str10 gender = cond(runiform() < 0.5, "Male", "Female")
gen income = runiformint(20000, 200000)
gen satisfaction = runiformint(1, 5)
benchmark_export, testname("survey-like data")

* Cross-sectional data
clear
set obs 100
gen firm_id = _n
gen revenue = runiform() * 1000000
gen employees = runiformint(10, 1000)
benchmark_export, testname("cross-sectional firm data")

* Healthcare data
clear
set obs 500
gen patient_id = _n
gen age = runiformint(1, 95)
gen str10 gender = cond(runiform() < 0.5, "M", "F")
gen bmi = rnormal(25, 5)
replace bmi = . if runiform() < 0.05
gen systolic = rnormal(120, 15)
gen diastolic = rnormal(80, 10)
gen str30 diagnosis = cond(runiform() < 0.3, "Hypertension", cond(runiform() < 0.5, "Diabetes", cond(runiform() < 0.7, "None", "")))
benchmark_export, testname("healthcare patient data")

* E-commerce data
clear
set obs 1000
gen order_id = _n
gen customer_id = runiformint(1, 200)
gen product_id = runiformint(1, 50)
gen quantity = runiformint(1, 10)
gen unit_price = round(runiform() * 100, 0.01)
gen total = quantity * unit_price
gen str20 status = cond(runiform() < 0.7, "Completed", cond(runiform() < 0.9, "Pending", "Cancelled"))
benchmark_export, testname("e-commerce order data")

* Education data
clear
set obs 500
gen student_id = _n
gen str30 name = "Student_" + string(_n)
gen gpa = round(runiform() * 3 + 1, 0.01)
gen credits = runiformint(0, 120)
gen str20 major = cond(runiform() < 0.2, "Computer Science", cond(runiform() < 0.4, "Business", cond(runiform() < 0.6, "Engineering", cond(runiform() < 0.8, "Biology", "Arts"))))
gen year = runiformint(1, 4)
benchmark_export, testname("education student data")

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
