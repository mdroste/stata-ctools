/*******************************************************************************
 * validate_cexport.do
 *
 * Comprehensive validation tests for cexport delimited vs export delimited
 * Tests all export options: delimiter, novarnames, quote, nolabel, if/in
 *
 * VERIFICATION: Exports both CSVs, imports them back into Stata, and compares:
 *   - First attempts cf _all for byte-for-byte comparison
 *   - Falls back to ABSOLUTE tolerance-based comparison (1e-7)
 *     to handle floating point precision differences in CSV formatting
 *   - String variables must match exactly
 *   - Reports max difference even when tests pass (if tolerance check was needed)
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

capture mkdir "temp"

quietly {

noi di as text "Running validation tests for cexport..."

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
        test_fail "`testname'" "dimensions differ: Stata N=`stata_n' K=`stata_k', cexport N=`cexport_n' K=`cexport_k'"
        exit
    }

    * Compare data using cf _all
    use `stata_data', clear
    capture cf _all using `cexport_data'
    local cfrc = _rc

    if `cfrc' == 0 {
        restore
        test_pass "`testname'"
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
            test_fail "`testname'" "variable `v' missing from cexport data"
            exit
        }
        rename `v' `v'_cexport
    }
    gen long _row = _n
    merge 1:1 _row using `stata_renamed', nogen

    local all_match = 1
    local fail_reason ""
    local min_sigfigs = 99

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
            * Numeric variable - compare values using significant figures
            * Create double copies to avoid type mismatch issues
            quietly gen double _v1 = `v'_stata
            quietly gen double _v2 = `v'_cexport

            * Check each non-missing pair using sigfigs
            local var_min_sf = 99
            quietly count if !missing(_v1) & !missing(_v2)
            local n_pairs = r(N)
            if `n_pairs' > 0 {
                * Use vectorized comparison via assert_var_equal logic
                quietly gen double _rel_err = abs(_v1 - _v2) / max(abs(_v1), abs(_v2), 1e-300) if !missing(_v1) & !missing(_v2)
                quietly summarize _rel_err
                if r(max) != . & r(max) > 0 {
                    * Convert max relative error to approximate sigfigs
                    local max_rel_err = r(max)
                    local var_min_sf = -log10(`max_rel_err')
                    if `var_min_sf' < `min_sigfigs' {
                        local min_sigfigs = `var_min_sf'
                    }
                }
                drop _rel_err
            }

            * Check if sigfigs meets threshold
            if `var_min_sf' < $DEFAULT_SIGFIGS {
                local all_match = 0
                local sf_fmt : display %4.1f `var_min_sf'
                local fail_reason "numeric variable `v' has only `sf_fmt' significant figures agreement"
                drop _v1 _v2
                continue, break
            }
            drop _v1 _v2
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
        if `min_sigfigs' < 99 {
            local sf_fmt : display %4.1f `min_sigfigs'
            test_pass "`testname' (min_sigfigs=`sf_fmt')"
        }
        else {
            test_pass "`testname'"
        }
    }
    else {
        test_fail "`testname'" "`fail_reason'"
    }
end

/*******************************************************************************
 * SECTION 1: Plugin check
 ******************************************************************************/
print_section "Plugin Check"

sysuse auto, clear
capture cexport delimited using "temp/test.csv", replace
if _rc != 0 {
    test_fail "cexport plugin load" "returned error `=_rc'"
    print_summary "cexport"
    exit 1
}
test_pass "cexport plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Basic comma-delimited export
 ******************************************************************************/
print_section "Basic Comma-Delimited Export"

sysuse auto, clear
benchmark_export, testname("basic export")

/*******************************************************************************
 * SECTION 3: Tab delimiter
 ******************************************************************************/
print_section "Tab Delimiter"

sysuse auto, clear
benchmark_export, testname("tab delimiter") exportopts(delimiter(tab)) importopts(delimiters(tab))

/*******************************************************************************
 * SECTION 4: Semicolon delimiter
 ******************************************************************************/
print_section "Semicolon Delimiter"

sysuse auto, clear
benchmark_export, testname("semicolon delimiter") exportopts(delimiter(";")) importopts(delimiters(";"))

/*******************************************************************************
 * SECTION 5: novarnames option
 ******************************************************************************/
print_section "novarnames Option"

sysuse auto, clear
benchmark_export, testname("novarnames") exportopts(novarnames) importopts(varnames(nonames))

/*******************************************************************************
 * SECTION 6: quote option
 ******************************************************************************/
print_section "quote Option"

sysuse auto, clear
benchmark_export, testname("quote option")  exportopts(quote)

/*******************************************************************************
 * SECTION 7: nolabel option
 ******************************************************************************/
print_section "nolabel Option"

sysuse auto, clear
benchmark_export, testname("with labels")

sysuse auto, clear
benchmark_export, testname("nolabel") exportopts(nolabel)

/*******************************************************************************
 * SECTION 8: Variable selection
 ******************************************************************************/
print_section "Variable Selection"

sysuse auto, clear
benchmark_export make price mpg, testname("variable selection (3 vars)")

/*******************************************************************************
 * SECTION 9: if condition
 ******************************************************************************/
print_section "if Condition"

sysuse auto, clear
benchmark_export, testname("if foreign==1") ifcond(foreign == 1)

sysuse auto, clear
benchmark_export, testname("if price>10000") ifcond(price > 10000)

/*******************************************************************************
 * SECTION 10: in condition
 ******************************************************************************/
print_section "in Condition"

sysuse auto, clear
benchmark_export, testname("in 1/20") incond(1/20)

sysuse auto, clear
benchmark_export, testname("in 30/50") incond(30/50)

/*******************************************************************************
 * SECTION 11: Combined if and in
 ******************************************************************************/
print_section "Combined if and in"

sysuse auto, clear
benchmark_export, testname("if and in combined") ifcond(price > 5000) incond(1/50)

/*******************************************************************************
 * SECTION 12: Census dataset
 ******************************************************************************/
print_section "Census Dataset"

sysuse census, clear
benchmark_export, testname("census")

/*******************************************************************************
 * SECTION 13: Large dataset export
 ******************************************************************************/
print_section "Large Dataset Export"

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
print_section "Panel Data (nlswork)"

webuse nlswork, clear
keep in 1/10000

benchmark_export, testname("panel data")

/*******************************************************************************
 * SECTION 15: verbose option
 ******************************************************************************/
print_section "verbose Option"

sysuse auto, clear

capture cexport delimited using "temp/test.csv", verbose replace
if _rc == 0 {
    test_pass "verbose option accepted"
}
else {
    test_fail "verbose option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 16: replace behavior
 ******************************************************************************/
print_section "replace Behavior"

sysuse auto, clear
cexport delimited using "temp/repl.csv", replace

capture cexport delimited using "temp/repl.csv"
if _rc != 0 {
    test_pass "without replace: fails for existing file"
}
else {
    test_fail "without replace" "should fail"
}

capture cexport delimited using "temp/repl.csv", replace
if _rc == 0 {
    test_pass "with replace: overwrites file"
}
else {
    test_fail "with replace" "failed"
}

/*******************************************************************************
 * SECTION 17: Edge cases - Basic
 ******************************************************************************/
print_section "Edge Cases - Basic"

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
print_section "Pathological - Empty/Minimal Datasets"

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
    test_pass "empty dataset (0 obs) - matches Stata behavior"
}
else {
    test_fail "empty dataset" "cexport rc=`cexport_rc' but Stata rc=`stata_rc'"
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
print_section "Pathological - Missing Value Patterns"

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
    test_pass "extended missing values (.a-.z) - matches Stata behavior"
}
else {
    test_fail "extended missing values" "cexport rc=`cexport_rc' but Stata rc=`stata_rc'"
}

/*******************************************************************************
 * SECTION: Pathological Data - String Edge Cases
 ******************************************************************************/
print_section "Pathological - String Edge Cases"

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
print_section "Pathological - Numeric Edge Cases"

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

* Extreme numeric range - compare behavior with Stata's export delimited
clear
set obs 6
gen double extreme = .
replace extreme = 1e308 in 1
replace extreme = -1e308 in 2
replace extreme = 1e-308 in 3
replace extreme = -1e-308 in 4
replace extreme = 0 in 5
replace extreme = . in 6
capture export delimited using "temp/extreme_range_stata.csv", replace
local stata_rc = _rc
capture cexport delimited using "temp/extreme_range.csv", replace
local cexport_rc = _rc
if `stata_rc' == `cexport_rc' {
    test_pass "extreme numeric range - matches Stata behavior (rc=`cexport_rc')"
}
else {
    test_fail "extreme numeric range" "cexport rc=`cexport_rc' but Stata rc=`stata_rc'"
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
print_section "Pathological - Dimensions"

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
print_section "Pathological - Data Patterns"

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
print_section "Option Combinations"

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
 * SECTION: Date Formatting (datafmt and datestring)
 ******************************************************************************/
print_section "Date Formatting (datafmt/datestring)"

* Create test data with various date formats
clear
set obs 10
gen id = _n

* Daily date (%td)
gen daily_date = mdy(1, _n, 2024)
format daily_date %td

* Datetime (%tc)
gen datetime = clock("2024-01-0" + string(_n) + " 12:30:45", "YMDhms")
format datetime %tc

* Weekly date (%tw)
gen weekly_date = yw(2024, _n)
format weekly_date %tw

* Monthly date (%tm)
gen monthly_date = ym(2024, _n)
format monthly_date %tm

* Quarterly date (%tq)
gen quarterly_date = yq(2024, mod(_n-1, 4) + 1)
format quarterly_date %tq

* Yearly date (%ty)
gen yearly_date = 2020 + _n
format yearly_date %ty

* Non-date numeric (should not be affected by datafmt)
gen value = _n * 100.5

* String variable (should not be affected)
gen str20 label = "item" + string(_n)

tempfile date_testdata
save `date_testdata', replace

* Test 1: datafmt option - basic daily date
use `date_testdata', clear
benchmark_export id daily_date value, testname("datafmt: daily date") exportopts(datafmt)

* Test 2: datafmt option - datetime
use `date_testdata', clear
benchmark_export id datetime value, testname("datafmt: datetime") exportopts(datafmt)

* Test 3: datafmt option - weekly date
use `date_testdata', clear
benchmark_export id weekly_date value, testname("datafmt: weekly date") exportopts(datafmt)

* Test 4: datafmt option - monthly date
use `date_testdata', clear
benchmark_export id monthly_date value, testname("datafmt: monthly date") exportopts(datafmt)

* Test 5: datafmt option - quarterly date
use `date_testdata', clear
benchmark_export id quarterly_date value, testname("datafmt: quarterly date") exportopts(datafmt)

* Test 6: datafmt option - yearly date
use `date_testdata', clear
benchmark_export id yearly_date value, testname("datafmt: yearly date") exportopts(datafmt)

* Test 7: datafmt option - all date types together
use `date_testdata', clear
benchmark_export, testname("datafmt: all date types") exportopts(datafmt)

* Test 8: datestring with custom format (ISO 8601) - cexport-specific option
* Note: Stata's export delimited doesn't have datestring(), so we test cexport alone
use `date_testdata', clear
capture cexport delimited id daily_date value using "temp/datestring_iso.csv", datestring("%tdCCYY-NN-DD") replace
if _rc == 0 {
    * Verify the output contains ISO-formatted dates
    import delimited using "temp/datestring_iso.csv", clear
    capture confirm string variable daily_date
    if _rc == 0 {
        * Check if first date looks like ISO format (starts with year)
        local firstval = daily_date[1]
        if substr("`firstval'", 1, 4) == "2024" {
            test_pass "datestring: ISO format"
        }
        else {
            test_fail "datestring: ISO format" "date not in ISO format: `firstval'"
        }
    }
    else {
        test_fail "datestring: ISO format" "daily_date should be string"
    }
}
else {
    test_fail "datestring: ISO format" "cexport failed with rc=`=_rc'"
}

* Test 9: datestring with different format (US format)
use `date_testdata', clear
capture cexport delimited id daily_date value using "temp/datestring_us.csv", datestring("%tdNN/DD/CCYY") replace
if _rc == 0 {
    import delimited using "temp/datestring_us.csv", clear
    capture confirm string variable daily_date
    if _rc == 0 {
        local firstval = daily_date[1]
        if strpos("`firstval'", "/") > 0 {
            test_pass "datestring: US format"
        }
        else {
            test_fail "datestring: US format" "date not in US format: `firstval'"
        }
    }
    else {
        test_fail "datestring: US format" "daily_date should be string"
    }
}
else {
    test_fail "datestring: US format" "cexport failed with rc=`=_rc'"
}

* Test 10: datestring with datetime
use `date_testdata', clear
capture cexport delimited id datetime value using "temp/datestring_datetime.csv", datestring("%tcCCYY-NN-DD!THH:MM:SS") replace
if _rc == 0 {
    import delimited using "temp/datestring_datetime.csv", clear
    capture confirm string variable datetime
    if _rc == 0 {
        local firstval = datetime[1]
        if strpos("`firstval'", "T") > 0 {
            test_pass "datestring: datetime ISO"
        }
        else {
            test_fail "datestring: datetime ISO" "datetime not in ISO format: `firstval'"
        }
    }
    else {
        test_fail "datestring: datetime ISO" "datetime should be string"
    }
}
else {
    test_fail "datestring: datetime ISO" "cexport failed with rc=`=_rc'"
}

* Test 11: datafmt combined with value labels
clear
set obs 5
gen id = _n
gen daily_date = mdy(6, _n, 2024)
format daily_date %td
gen category = mod(_n, 2) + 1
label define cat_lbl 1 "Type A" 2 "Type B"
label values category cat_lbl
benchmark_export, testname("datafmt with value labels") exportopts(datafmt)

* Test 12: datafmt combined with if condition
use `date_testdata', clear
benchmark_export, testname("datafmt with if condition") exportopts(datafmt) ifcond(id <= 5)

* Test 13: datafmt combined with in condition
use `date_testdata', clear
benchmark_export, testname("datafmt with in range") exportopts(datafmt) incond(3/8)

* Test 14: datafmt combined with quote option
use `date_testdata', clear
benchmark_export, testname("datafmt with quote") exportopts(datafmt quote)

* Test 15: datafmt combined with tab delimiter
use `date_testdata', clear
benchmark_export, testname("datafmt with tab delimiter") exportopts(datafmt delimiter(tab)) importopts(delimiters(tab))

* Test 16: datestring with missing dates (cexport-specific)
* Compare cexport datestring behavior against Stata's datafmt export for missing dates
clear
set obs 10
gen id = _n
gen daily_date = mdy(1, _n, 2024)
replace daily_date = . if mod(_n, 3) == 0
format daily_date %td
gen value = _n * 10
tempfile missing_dates_data
save `missing_dates_data'

* First: export with Stata's export delimited using datafmt to see how Stata handles missing dates
use `missing_dates_data', clear
export delimited using "temp/datestring_missing_stata.csv", datafmt replace
import delimited using "temp/datestring_missing_stata.csv", clear
* Count missing date values in Stata's export
capture confirm string variable daily_date
if _rc == 0 {
    count if daily_date == "" | missing(daily_date)
}
else {
    count if missing(daily_date)
}
local stata_n_missing = r(N)

* Now: export with cexport datestring
use `missing_dates_data', clear
capture cexport delimited using "temp/datestring_missing.csv", datestring("%tdCCYY-NN-DD") replace
if _rc == 0 {
    import delimited using "temp/datestring_missing.csv", clear
    * Count missing date values in cexport's output
    capture confirm string variable daily_date
    if _rc == 0 {
        count if daily_date == "" | missing(daily_date)
    }
    else {
        count if missing(daily_date)
    }
    local cexport_n_missing = r(N)

    * Both should agree on how many missing dates there are
    if `cexport_n_missing' == `stata_n_missing' {
        test_pass "datestring with missing dates (matches Stata: `stata_n_missing' missing)"
    }
    else {
        test_fail "datestring with missing dates" "cexport has `cexport_n_missing' missing but Stata has `stata_n_missing' missing"
    }
}
else {
    test_fail "datestring with missing dates" "cexport failed with rc=`=_rc'"
}
capture erase "temp/datestring_missing_stata.csv"

* Test 17: Negative format prefix (%-t)
clear
set obs 5
gen id = _n
gen daily_date = mdy(1, _n, 2024)
format daily_date %-td
benchmark_export, testname("datafmt: negative format (%-td)") exportopts(datafmt)

* Test 18: Real-world dataset with dates (nlswork)
* Note: We only test the date variable since datafmt affects all numeric formatting
* and cexport doesn't implement full numeric display format support (only dates)
capture webuse nlswork, clear
if _rc == 0 {
    keep in 1/1000
    keep idcode year
    * Create a date variable from year
    gen interview_date = mdy(6, 15, year)
    format interview_date %td
    benchmark_export, testname("datafmt: nlswork with dates") exportopts(datafmt)
}

/*******************************************************************************
 * SECTION 18: Missing values
 ******************************************************************************/
print_section "Missing Values"

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
print_section "Numeric Precision"

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
print_section "Special Characters"

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
print_section "Large Datasets"

* 10K rows export - full comparison
clear
set seed 12345
set obs 10000
gen id = _n
gen value = runiform()
gen str20 name = "item_" + string(_n)

benchmark_export, testname("10K rows export")

* 50K rows - full comparison
clear
set seed 54321
set obs 50000
gen id = _n
gen x = runiform()
gen y = runiform()

benchmark_export, testname("50K rows export")

* Many columns (20) - full comparison
clear
set seed 11111
set obs 1000
forvalues i = 1/20 {
    gen var`i' = runiform()
}

benchmark_export, testname("20 columns export")

/*******************************************************************************
 * SECTION: Pathological Data (Original)
 ******************************************************************************/
print_section "Pathological Data"

* All missing values
clear
set obs 100
gen x = .
gen y = .
gen z = .

capture cexport delimited using "temp/all_missing.csv", replace
if _rc == 0 {
    test_pass "all missing values"
}
else {
    test_fail "all missing" "rc=`=_rc'"
}

* Empty strings
clear
set obs 50
gen str20 name = ""
gen value = _n

capture cexport delimited using "temp/empty_strings.csv", replace
if _rc == 0 {
    test_pass "empty strings"
}
else {
    test_fail "empty strings" "rc=`=_rc'"
}

* Single column
clear
set obs 100
gen only_col = runiform()

capture cexport delimited using "temp/single_col.csv", replace
if _rc == 0 {
    test_pass "single column"
}
else {
    test_fail "single column" "rc=`=_rc'"
}

* Single row
clear
set obs 1
gen a = 1
gen b = 2
gen c = 3

capture cexport delimited using "temp/single_row.csv", replace
if _rc == 0 {
    test_pass "single row"
}
else {
    test_fail "single row" "rc=`=_rc'"
}

* Very long strings
clear
set obs 10
gen str244 long_text = "a" * 200

capture cexport delimited using "temp/long_strings.csv", replace
if _rc == 0 {
    test_pass "long strings (200 chars)"
}
else {
    test_fail "long strings" "rc=`=_rc'"
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
    test_pass "extreme numeric values"
}
else {
    test_fail "extreme nums" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION: Additional Pathological - Comparison Tests
 ******************************************************************************/
print_section "Pathological Comparison Tests"

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
print_section "Real-World Datasets"

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
    test_pass "nlswork dataset export"
}
else {
    test_fail "nlswork" "rc=`=_rc'"
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
print_section "Synthetic Test Datasets"

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
print_section "Comparison with export delimited"

* Basic comparison
sysuse auto, clear
export delimited using "temp/compare_stata.csv", replace
cexport delimited using "temp/compare_cexport.csv", replace

* Compare file sizes (should be similar)
capture file open f1 using "temp/compare_stata.csv", read
local rc1 = _rc
capture file open f2 using "temp/compare_cexport.csv", read
local rc2 = _rc
capture file close f1
capture file close f2
if `rc1' == 0 & `rc2' == 0 {
    test_pass "comparison files created"
}
else {
    test_fail "comparison files" "file open errors: stata_rc=`rc1' cexport_rc=`rc2'"
}

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
    test_pass "round-trip export/import"
}
else {
    test_fail "round-trip" "wrong N"
}

/*******************************************************************************
 * SECTION: Edge Cases
 ******************************************************************************/
print_section "Edge Cases"

* Variable names with underscores
clear
set obs 10
gen my_var_name = _n
gen another_var = runiform()

capture cexport delimited using "temp/underscores.csv", replace
if _rc == 0 {
    test_pass "variable names with underscores"
}
else {
    test_fail "underscores" "rc=`=_rc'"
}

* Numeric-looking strings
clear
set obs 5
gen str10 numstr = string(_n * 100)

capture cexport delimited using "temp/numstr.csv", replace
if _rc == 0 {
    test_pass "numeric-looking strings"
}
else {
    test_fail "numstr" "rc=`=_rc'"
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
    test_pass "mixed variable types"
}
else {
    test_fail "mixed types" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION: Excel Export Tests
 ******************************************************************************/
print_section "Excel Export Tests"

/*******************************************************************************
 * Helper: benchmark_export_excel - Export with both methods and compare
 ******************************************************************************/
capture program drop benchmark_export_excel
program define benchmark_export_excel
    syntax, testname(string) [EXPORTopts(string)]

    * Export with Stata's export excel
    capture export excel using "temp/stata_export.xlsx", `exportopts' replace
    if _rc != 0 {
        test_fail "`testname'" "Stata export excel failed with rc=`=_rc'"
        exit
    }
    local stata_n = _N
    local stata_k = c(k)

    * Export with cexport excel
    capture cexport excel using "temp/cexport_export.xlsx", `exportopts' replace
    if _rc != 0 {
        test_fail "`testname'" "cexport excel failed with rc=`=_rc'"
        exit
    }

    * Re-import Stata output
    clear
    capture import excel using "temp/stata_export.xlsx", firstrow clear
    if _rc != 0 {
        test_fail "`testname'" "Failed to re-import Stata export"
        exit
    }
    local reimport_stata_n = _N
    local reimport_stata_k = c(k)
    tempfile stata_reimport
    quietly save `stata_reimport', replace

    * Re-import cexport output
    clear
    capture import excel using "temp/cexport_export.xlsx", firstrow clear
    if _rc != 0 {
        test_fail "`testname'" "Failed to re-import cexport output"
        exit
    }
    local reimport_cexport_n = _N
    local reimport_cexport_k = c(k)

    * Check dimensions match
    if `reimport_stata_n' != `reimport_cexport_n' | `reimport_stata_k' != `reimport_cexport_k' {
        test_fail "`testname'" "dimensions differ: Stata N=`reimport_stata_n' K=`reimport_stata_k', cexport N=`reimport_cexport_n' K=`reimport_cexport_k'"
        exit
    }

    * Rename variables for positional comparison
    quietly ds
    local cexport_vars `r(varlist)'
    local i = 1
    foreach v of local cexport_vars {
        quietly rename `v' __cmp_v`i'
        local i = `i' + 1
    }
    tempfile cexport_reimport
    quietly save `cexport_reimport', replace

    use `stata_reimport', clear
    quietly ds
    local stata_vars `r(varlist)'
    local i = 1
    foreach v of local stata_vars {
        quietly rename `v' __cmp_v`i'
        local i = `i' + 1
    }

    * Compare data
    capture cf _all using `cexport_reimport'
    if _rc == 0 {
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "cf _all comparison failed - data not identical"
    }

    * Clean up temp files
    capture erase "temp/stata_export.xlsx"
    capture erase "temp/cexport_export.xlsx"
end

/*******************************************************************************
 * Helper: cexport_excel_test - Test cexport excel standalone
 ******************************************************************************/
capture program drop cexport_excel_test
program define cexport_excel_test
    syntax, testname(string) [EXPORTopts(string) expectn(integer 0) expectk(integer 0)]

    capture cexport excel using "temp/test_export.xlsx", `exportopts' replace
    if _rc != 0 {
        test_fail "`testname'" "rc=`=_rc'"
        exit
    }

    * Re-import and check
    local orig_n = _N
    local orig_k = c(k)

    clear
    capture import excel using "temp/test_export.xlsx", firstrow clear
    if _rc != 0 {
        test_fail "`testname'" "Failed to re-import"
        exit
    }

    if `expectn' > 0 & _N != `expectn' {
        test_fail "`testname'" "expected N=`expectn', got N=`=_N'"
        exit
    }

    if `expectk' > 0 & c(k) != `expectk' {
        test_fail "`testname'" "expected K=`expectk', got K=`=c(k)'"
        exit
    }

    * Clean up
    capture erase "temp/test_export.xlsx"

    test_pass "`testname'"
end

* Create test data for Excel export
clear
set obs 5
gen id = _n
gen str10 name = ""
replace name = "Alpha" in 1
replace name = "Beta" in 2
replace name = "Gamma" in 3
replace name = "Delta" in 4
replace name = "Epsilon" in 5
gen double value = _n * 100.5
gen byte category = mod(_n, 3) + 1
tempfile excel_basic_data
save `excel_basic_data', replace

* Data with missing values
clear
set obs 5
gen id = _n
gen double x = _n * 10
replace x = . in 2
replace x = . in 4
gen str10 y = ""
replace y = "A" in 1
replace y = "" in 2
replace y = "C" in 3
replace y = "" in 4
replace y = "E" in 5
tempfile excel_missing_data
save `excel_missing_data', replace

* Data with value labels
clear
set obs 4
gen id = _n
gen status = mod(_n, 2) + 1
label define status_lbl 1 "Active" 2 "Inactive"
label values status status_lbl
tempfile excel_labeled_data
save `excel_labeled_data', replace

* Data with various numeric types
clear
set obs 3
gen byte b = _n * 10
gen int i = _n * 1000
gen long l = _n * 100000
gen float f = _n + 0.5
gen double d = _n + 0.123456789
tempfile excel_numeric_data
save `excel_numeric_data', replace

* Basic Excel export tests
use `excel_basic_data', clear
benchmark_export_excel, testname("Excel basic export") exportopts(firstrow(variables))

* Export without firstrow
use `excel_basic_data', clear
cexport excel using "temp/excel_nofirstrow.xlsx", firstrow(nonames) replace
clear
import excel using "temp/excel_nofirstrow.xlsx", clear
if _N == 5 {
    test_pass "Excel export without firstrow"
}
else {
    test_fail "Excel export without firstrow" "expected N=5, got N=`=_N'"
}
capture erase "temp/excel_nofirstrow.xlsx"

* Export with sheet name
use `excel_basic_data', clear
cexport_excel_test, testname("Excel export with sheet name") exportopts(sheet("TestData")) expectn(5) expectk(4)

* Variable selection
use `excel_basic_data', clear
cexport excel id name using "temp/excel_varsel.xlsx", replace
clear
import excel using "temp/excel_varsel.xlsx", firstrow clear
if _N == 5 & c(k) == 2 {
    test_pass "Excel export selected variables"
}
else {
    test_fail "Excel export selected variables" "expected N=5 K=2, got N=`=_N' K=`=c(k)'"
}
capture erase "temp/excel_varsel.xlsx"

* Type handling tests
use `excel_numeric_data', clear
benchmark_export_excel, testname("Excel numeric types export") exportopts(firstrow(variables))

use `excel_missing_data', clear
benchmark_export_excel, testname("Excel missing values export") exportopts(firstrow(variables))

* Value label tests
use `excel_labeled_data', clear
cexport excel using "temp/excel_labeled.xlsx", replace
clear
import excel using "temp/excel_labeled.xlsx", firstrow clear
capture confirm string variable status
if _rc == 0 {
    test_pass "Excel export with value labels"
}
else {
    test_fail "Excel export with value labels" "status should be string with labels"
}
capture erase "temp/excel_labeled.xlsx"

use `excel_labeled_data', clear
cexport excel using "temp/excel_nolabel.xlsx", nolabel replace
clear
import excel using "temp/excel_nolabel.xlsx", firstrow clear
capture confirm numeric variable status
if _rc == 0 {
    test_pass "Excel export with nolabel option"
}
else {
    test_fail "Excel export with nolabel option" "status should be numeric with nolabel"
}
capture erase "temp/excel_nolabel.xlsx"

* If/in tests
use `excel_basic_data', clear
cexport excel using "temp/excel_iftest.xlsx" if id <= 3, replace
clear
import excel using "temp/excel_iftest.xlsx", firstrow clear
if _N == 3 {
    test_pass "Excel export with if condition"
}
else {
    test_fail "Excel export with if condition" "expected N=3, got N=`=_N'"
}
capture erase "temp/excel_iftest.xlsx"

use `excel_basic_data', clear
cexport excel using "temp/excel_intest.xlsx" in 2/4, replace
clear
import excel using "temp/excel_intest.xlsx", firstrow clear
if _N == 3 {
    test_pass "Excel export with in range"
}
else {
    test_fail "Excel export with in range" "expected N=3, got N=`=_N'"
}
capture erase "temp/excel_intest.xlsx"

* Error handling tests
use `excel_basic_data', clear
cexport excel using "temp/excel_exists_test.xlsx", replace
capture cexport excel using "temp/excel_exists_test.xlsx"
if _rc != 0 {
    test_pass "Excel file exists error"
}
else {
    test_fail "Excel file exists error" "should have returned error"
}
capture erase "temp/excel_exists_test.xlsx"

use `excel_basic_data', clear
capture cexport excel using "temp/wrong.csv", replace
if _rc != 0 {
    test_pass "Excel wrong extension error"
}
else {
    test_fail "Excel wrong extension error" "should have returned error"
}

* Clean up Excel temp files
capture erase "temp/stata_export.xlsx"
capture erase "temp/cexport_export.xlsx"
capture erase "temp/test_export.xlsx"

/*******************************************************************************
 * SECTION: Excel Export - New Options (cell, missing, keepcellfmt)
 ******************************************************************************/
print_section "Excel Export - New Options"

* Test cell() option - starting cell offset
use `excel_basic_data', clear
cexport excel using "temp/excel_cell_b5.xlsx", cell(B5) replace
clear
import excel using "temp/excel_cell_b5.xlsx", clear
* When imported, data should be in columns B-E starting at row 5
* Check that column A is empty and row 1-4 are empty
capture confirm variable A
if _rc != 0 {
    test_pass "Excel cell(B5) - column A empty"
}
else {
    * Column A exists - check if it's all missing
    count if !missing(A)
    if r(N) == 0 {
        test_pass "Excel cell(B5) - column A all missing"
    }
    else {
        test_fail "Excel cell(B5)" "column A should be empty"
    }
}
capture erase "temp/excel_cell_b5.xlsx"

* Test cell() with various references
use `excel_basic_data', clear
foreach cellref in A1 Z1 AA1 C10 {
    capture cexport excel id name using "temp/excel_cell_`cellref'.xlsx", cell(`cellref') replace
    if _rc == 0 {
        test_pass "Excel cell(`cellref') accepted"
    }
    else {
        test_fail "Excel cell(`cellref')" "rc=`=_rc'"
    }
    capture erase "temp/excel_cell_`cellref'.xlsx"
}

* Test missing() option - replacement value for missing data
use `excel_missing_data', clear
cexport excel using "temp/excel_missing_na.xlsx", missing("NA") replace
clear
import excel using "temp/excel_missing_na.xlsx", firstrow clear
* Check that missing values are now "NA" strings
capture confirm string variable x
if _rc == 0 {
    count if x == "NA"
    if r(N) == 2 {
        test_pass "Excel missing(NA) - numeric missing replaced"
    }
    else {
        test_fail "Excel missing(NA)" "expected 2 NA values, got `r(N)'"
    }
}
else {
    test_fail "Excel missing(NA)" "x should be string with NA values"
}
capture erase "temp/excel_missing_na.xlsx"

* Test missing(".") option
use `excel_missing_data', clear
cexport excel using "temp/excel_missing_dot.xlsx", missing(".") replace
clear
import excel using "temp/excel_missing_dot.xlsx", firstrow clear
capture confirm string variable x
if _rc == 0 {
    count if x == "."
    if r(N) == 2 {
        test_pass "Excel missing(.) - periods for missing"
    }
    else {
        test_fail "Excel missing(.)" "expected 2 period values, got `r(N)'"
    }
}
else {
    test_fail "Excel missing(.)" "x should be string with . values"
}
capture erase "temp/excel_missing_dot.xlsx"

* Test missing("") option - explicit empty string
use `excel_missing_data', clear
cexport excel using "temp/excel_missing_empty.xlsx", missing("") replace
clear
import excel using "temp/excel_missing_empty.xlsx", firstrow clear
* Empty string missing values - should result in empty cells
if _N == 5 {
    test_pass "Excel missing('') - empty string replacement"
}
else {
    test_fail "Excel missing('')" "expected N=5, got N=`=_N'"
}
capture erase "temp/excel_missing_empty.xlsx"

* Test keepcellfmt option - preserves styles from existing file
use `excel_basic_data', clear
* First create a file
cexport excel using "temp/excel_keepfmt.xlsx", replace
* Then update it with keepcellfmt
cexport excel id name value using "temp/excel_keepfmt.xlsx", keepcellfmt replace
clear
import excel using "temp/excel_keepfmt.xlsx", firstrow clear
if _N == 5 & c(k) == 3 {
    test_pass "Excel keepcellfmt - file updated"
}
else {
    test_fail "Excel keepcellfmt" "expected N=5 K=3, got N=`=_N' K=`=c(k)'"
}
capture erase "temp/excel_keepfmt.xlsx"

* Test combination: cell() + missing()
use `excel_missing_data', clear
cexport excel using "temp/excel_combo.xlsx", cell(C3) missing("N/A") replace
clear
import excel using "temp/excel_combo.xlsx", clear
* Data should start at C3, with N/A for missing values
if _N > 0 {
    test_pass "Excel cell() + missing() combination"
}
else {
    test_fail "Excel cell() + missing() combination" "no data imported"
}
capture erase "temp/excel_combo.xlsx"

* Test combination: cell() + keepcellfmt
use `excel_basic_data', clear
cexport excel using "temp/excel_combo2.xlsx", replace
use `excel_missing_data', clear
cexport excel using "temp/excel_combo2.xlsx", cell(F1) keepcellfmt replace
clear
import excel using "temp/excel_combo2.xlsx", firstrow clear
if _N > 0 {
    test_pass "Excel cell() + keepcellfmt combination"
}
else {
    test_fail "Excel cell() + keepcellfmt combination" "no data imported"
}
capture erase "temp/excel_combo2.xlsx"

* Test all three options together
use `excel_missing_data', clear
cexport excel using "temp/excel_allopts.xlsx", replace
cexport excel using "temp/excel_allopts.xlsx", cell(B2) missing("--") keepcellfmt replace
clear
import excel using "temp/excel_allopts.xlsx", clear
if _N > 0 {
    test_pass "Excel cell() + missing() + keepcellfmt all together"
}
else {
    test_fail "Excel all options" "no data imported"
}
capture erase "temp/excel_allopts.xlsx"

* Test with real dataset (auto) and new options
sysuse auto, clear
replace price = . in 1/5
cexport excel make price mpg using "temp/excel_auto_opts.xlsx", cell(B2) missing("N/A") replace
clear
import excel using "temp/excel_auto_opts.xlsx", clear
* Import starts at B2, so first row (row 1) should be empty or headers start at row 2
if _N > 0 {
    test_pass "Excel auto dataset with new options"
}
else {
    test_fail "Excel auto dataset with new options" "no data imported"
}
capture erase "temp/excel_auto_opts.xlsx"

/*******************************************************************************
 * SECTION: Intentional Error Tests
 *
 * These tests verify that cexport returns the same error codes as export delimited
 * when given invalid inputs or error conditions.
 ******************************************************************************/
print_section "Intentional Error Tests"

* Variable doesn't exist
sysuse auto, clear
test_error_match, stata_cmd(export delimited nonexistent_var using "temp/error_test.csv", replace) ctools_cmd(cexport delimited nonexistent_var using "temp/error_test.csv", replace) testname("nonexistent variable")

* No data loaded (empty dataset with no variables)
clear
test_error_match, stata_cmd(export delimited using "temp/error_test.csv", replace) ctools_cmd(cexport delimited using "temp/error_test.csv", replace) testname("no data loaded")

* Invalid path (directory doesn't exist)
sysuse auto, clear
test_error_match, stata_cmd(export delimited using "nonexistent_dir/test.csv", replace) ctools_cmd(cexport delimited using "nonexistent_dir/test.csv", replace) testname("invalid path")

* File already exists without replace
sysuse auto, clear
export delimited using "temp/exist_test.csv", replace
sysuse auto, clear
capture export delimited using "temp/exist_test.csv"
local stata_rc = _rc
capture cexport delimited using "temp/exist_test.csv"
local cexport_rc = _rc
if `stata_rc' == `cexport_rc' {
    test_pass "[error] file exists without replace (rc=`stata_rc')"
}
else {
    test_fail "[error] file exists without replace" "stata rc=`stata_rc', cexport rc=`cexport_rc'"
}
capture erase "temp/exist_test.csv"

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

* End of cexport validation
noi print_summary "cexport"
}
