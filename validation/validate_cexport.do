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
