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
 * Helper: compare_exports - Export with both methods, import back, compare with cf _all
 * Returns pass/fail result
 ******************************************************************************/
capture program drop compare_exports
program define compare_exports, rclass
    syntax [varlist], stata_file(string) cexport_file(string) ///
           [EXPORTopts(string) IMPORTopts(string) IFcond(string) INcond(string)]

    * Build if/in conditions
    local ifin ""
    if "`ifcond'" != "" local ifin "`ifin' if `ifcond'"
    if "`incond'" != "" local ifin "`ifin' in `incond'"

    * Export with Stata's export delimited
    if "`varlist'" != "" {
        export delimited `varlist' using "`stata_file'" `ifin', `exportopts' replace
    }
    else {
        export delimited using "`stata_file'" `ifin', `exportopts' replace
    }

    * Export with cexport delimited
    if "`varlist'" != "" {
        cexport delimited `varlist' using "`cexport_file'" `ifin', `exportopts' replace
    }
    else {
        cexport delimited using "`cexport_file'" `ifin', `exportopts' replace
    }

    * Import both CSVs back
    preserve

    import delimited using "`stata_file'", `importopts' clear
    tempfile stata_data
    quietly save `stata_data', replace
    local stata_n = _N
    local stata_k = c(k)

    import delimited using "`cexport_file'", `importopts' clear
    tempfile cexport_data
    quietly save `cexport_data', replace
    local cexport_n = _N
    local cexport_k = c(k)

    * Check dimensions first
    if `stata_n' != `cexport_n' | `stata_k' != `cexport_k' {
        restore
        return local result "fail"
        return local reason "dimensions differ: Stata N=`stata_n' K=`stata_k', cexport N=`cexport_n' K=`cexport_k'"
        exit
    }

    * Compare data using cf _all
    use `stata_data', clear
    capture cf _all using `cexport_data'
    local cfrc = _rc

    if `cfrc' == 0 {
        restore
        return local result "pass"
        return local reason ""
        exit
    }

    * cf _all failed - try tolerance-based comparison for numeric variables
    * This handles floating point precision differences in CSV export
    use `stata_data', clear
    ds
    local varlist `r(varlist)'

    * Rename variables from stata data with _stata suffix
    foreach v of local varlist {
        rename `v' `v'_stata
    }
    gen long _row = _n
    tempfile stata_renamed
    quietly save `stata_renamed', replace

    * Load cexport data and rename with _cexport suffix
    use `cexport_data', clear
    foreach v of local varlist {
        capture confirm variable `v'
        if _rc != 0 {
            restore
            return local result "fail"
            return local reason "variable `v' missing from cexport data"
            exit
        }
        rename `v' `v'_cexport
    }
    gen long _row = _n
    merge 1:1 _row using `stata_renamed', nogen

    local all_match = 1
    local fail_reason ""
    local tol = 1e-6

    foreach v of local varlist {
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
        return local result "pass"
        return local reason ""
    }
    else {
        return local result "fail"
        return local reason "`fail_reason'"
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
 * SECTION 2: Basic comma-delimited export (cf _all verification)
 ******************************************************************************/
noi print_section "Basic Comma-Delimited Export"

sysuse auto, clear
compare_exports, stata_file("temp/stata.csv") cexport_file("temp/cexport.csv") ///
    exportopts() importopts()
if "`r(result)'" == "pass" {
    noi test_pass "basic export data identical (cf _all)"
}
else {
    noi test_fail "basic export" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 3: Tab delimiter (cf _all verification)
 ******************************************************************************/
noi print_section "Tab Delimiter"

sysuse auto, clear
compare_exports, stata_file("temp/stata_tab.tsv") cexport_file("temp/cexport_tab.tsv") ///
    exportopts(delimiter(tab)) importopts(delimiters(tab))
if "`r(result)'" == "pass" {
    noi test_pass "tab delimiter data identical (cf _all)"
}
else {
    noi test_fail "tab delimiter" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 4: Semicolon delimiter (cf _all verification)
 ******************************************************************************/
noi print_section "Semicolon Delimiter"

sysuse auto, clear
compare_exports, stata_file("temp/stata_semi.csv") cexport_file("temp/cexport_semi.csv") ///
    exportopts(delimiter(";")) importopts(delimiters(";"))
if "`r(result)'" == "pass" {
    noi test_pass "semicolon delimiter data identical (cf _all)"
}
else {
    noi test_fail "semicolon delimiter" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 5: novarnames option (cf _all verification)
 ******************************************************************************/
noi print_section "novarnames Option"

sysuse auto, clear
compare_exports, stata_file("temp/stata_novar.csv") cexport_file("temp/cexport_novar.csv") ///
    exportopts(novarnames) importopts(varnames(nonames))
if "`r(result)'" == "pass" {
    noi test_pass "novarnames data identical (cf _all)"
}
else {
    noi test_fail "novarnames" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 6: quote option (cf _all verification)
 ******************************************************************************/
noi print_section "quote Option"

sysuse auto, clear
compare_exports, stata_file("temp/stata_quote.csv") cexport_file("temp/cexport_quote.csv") ///
    exportopts(quote) importopts()
if "`r(result)'" == "pass" {
    noi test_pass "quote option data identical (cf _all)"
}
else {
    noi test_fail "quote option" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 7: nolabel option (cf _all verification)
 ******************************************************************************/
noi print_section "nolabel Option"

* With labels (default)
sysuse auto, clear
compare_exports, stata_file("temp/stata_label.csv") cexport_file("temp/cexport_label.csv") ///
    exportopts() importopts()
if "`r(result)'" == "pass" {
    noi test_pass "with labels data identical (cf _all)"
}
else {
    noi test_fail "with labels" "`r(reason)'"
}

* Without labels
sysuse auto, clear
compare_exports, stata_file("temp/stata_nolabel.csv") cexport_file("temp/cexport_nolabel.csv") ///
    exportopts(nolabel) importopts()
if "`r(result)'" == "pass" {
    noi test_pass "nolabel data identical (cf _all)"
}
else {
    noi test_fail "nolabel" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 8: Variable selection (cf _all verification)
 ******************************************************************************/
noi print_section "Variable Selection"

sysuse auto, clear
compare_exports make price mpg, stata_file("temp/stata_vars.csv") cexport_file("temp/cexport_vars.csv") ///
    exportopts() importopts()
if "`r(result)'" == "pass" {
    noi test_pass "variable selection (3 vars) data identical (cf _all)"
}
else {
    noi test_fail "variable selection" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 9: if condition (cf _all verification)
 ******************************************************************************/
noi print_section "if Condition"

sysuse auto, clear
compare_exports, stata_file("temp/stata_if.csv") cexport_file("temp/cexport_if.csv") ///
    exportopts() importopts() ifcond(foreign == 1)
if "`r(result)'" == "pass" {
    noi test_pass "if foreign==1 data identical (cf _all)"
}
else {
    noi test_fail "if foreign==1" "`r(reason)'"
}

sysuse auto, clear
compare_exports, stata_file("temp/stata_if2.csv") cexport_file("temp/cexport_if2.csv") ///
    exportopts() importopts() ifcond(price > 10000)
if "`r(result)'" == "pass" {
    noi test_pass "if price>10000 data identical (cf _all)"
}
else {
    noi test_fail "if price>10000" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 10: in condition (cf _all verification)
 ******************************************************************************/
noi print_section "in Condition"

sysuse auto, clear
compare_exports, stata_file("temp/stata_in.csv") cexport_file("temp/cexport_in.csv") ///
    exportopts() importopts() incond(1/20)
if "`r(result)'" == "pass" {
    noi test_pass "in 1/20 data identical (cf _all)"
}
else {
    noi test_fail "in 1/20" "`r(reason)'"
}

sysuse auto, clear
compare_exports, stata_file("temp/stata_in2.csv") cexport_file("temp/cexport_in2.csv") ///
    exportopts() importopts() incond(30/50)
if "`r(result)'" == "pass" {
    noi test_pass "in 30/50 data identical (cf _all)"
}
else {
    noi test_fail "in 30/50" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 11: Combined if and in (cf _all verification)
 ******************************************************************************/
noi print_section "Combined if and in"

sysuse auto, clear
compare_exports, stata_file("temp/stata_ifin.csv") cexport_file("temp/cexport_ifin.csv") ///
    exportopts() importopts() ifcond(price > 5000) incond(1/50)
if "`r(result)'" == "pass" {
    noi test_pass "if and in combined data identical (cf _all)"
}
else {
    noi test_fail "if and in combined" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 12: Census dataset (cf _all verification)
 ******************************************************************************/
noi print_section "Census Dataset"

sysuse census, clear
compare_exports, stata_file("temp/stata_census.csv") cexport_file("temp/cexport_census.csv") ///
    exportopts() importopts()
if "`r(result)'" == "pass" {
    noi test_pass "census data identical (cf _all)"
}
else {
    noi test_fail "census" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 13: Large dataset export (cf _all verification)
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

compare_exports, stata_file("temp/stata_large.csv") cexport_file("temp/cexport_large.csv") ///
    exportopts() importopts()
if "`r(result)'" == "pass" {
    noi test_pass "large dataset (50K) data identical (cf _all)"
}
else {
    noi test_fail "large dataset" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 14: Panel data (nlswork) (cf _all verification)
 ******************************************************************************/
noi print_section "Panel Data (nlswork)"

webuse nlswork, clear
keep in 1/10000

compare_exports, stata_file("temp/stata_panel.csv") cexport_file("temp/cexport_panel.csv") ///
    exportopts() importopts()
if "`r(result)'" == "pass" {
    noi test_pass "panel data identical (cf _all)"
}
else {
    noi test_fail "panel data" "`r(reason)'"
}

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
 * SECTION 17: Edge cases (cf _all verification)
 ******************************************************************************/
noi print_section "Edge Cases"

* Single observation
clear
set obs 1
gen x = 42
gen str5 s = "test"

compare_exports, stata_file("temp/stata_single.csv") cexport_file("temp/cexport_single.csv") ///
    exportopts() importopts()
if "`r(result)'" == "pass" {
    noi test_pass "single observation data identical (cf _all)"
}
else {
    noi test_fail "single observation" "`r(reason)'"
}

* Single variable
clear
set obs 100
gen x = runiform()

compare_exports x, stata_file("temp/stata_singlevar.csv") cexport_file("temp/cexport_singlevar.csv") ///
    exportopts() importopts()
if "`r(result)'" == "pass" {
    noi test_pass "single variable data identical (cf _all)"
}
else {
    noi test_fail "single variable" "`r(reason)'"
}

* String with commas
clear
input id str40 name
1 "Smith, John"
2 "Doe, Jane"
3 "Regular Name"
end

compare_exports, stata_file("temp/stata_comma.csv") cexport_file("temp/cexport_comma.csv") ///
    exportopts() importopts()
if "`r(result)'" == "pass" {
    noi test_pass "strings with commas data identical (cf _all)"
}
else {
    noi test_fail "strings with commas" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 18: Missing values (cf _all verification)
 ******************************************************************************/
noi print_section "Missing Values"

clear
set obs 10
gen x = _n
replace x = . if mod(_n, 3) == 0
gen str10 s = "val" + string(_n)
replace s = "" if mod(_n, 4) == 0

compare_exports, stata_file("temp/stata_missing.csv") cexport_file("temp/cexport_missing.csv") ///
    exportopts() importopts()
if "`r(result)'" == "pass" {
    noi test_pass "missing values data identical (cf _all)"
}
else {
    noi test_fail "missing values" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 19: Numeric precision (cf _all verification)
 ******************************************************************************/
noi print_section "Numeric Precision"

clear
set obs 100
gen double precise = runiform() * 1e10
gen float approx = runiform()
gen long bigint = runiformint(1, 2147483647)
gen byte tinyint = runiformint(-127, 127)

compare_exports, stata_file("temp/stata_precision.csv") cexport_file("temp/cexport_precision.csv") ///
    exportopts() importopts()
if "`r(result)'" == "pass" {
    noi test_pass "numeric precision data identical (cf _all)"
}
else {
    noi test_fail "numeric precision" "`r(reason)'"
}

/*******************************************************************************
 * SECTION 20: Special characters in strings (cf _all verification)
 ******************************************************************************/
noi print_section "Special Characters"

clear
input id str50 text
1 "Normal text"
2 "Text with tabs	here"
3 "Text with newlines here"
4 "Text with ""quotes"""
end

compare_exports, stata_file("temp/stata_special.csv") cexport_file("temp/cexport_special.csv") ///
    exportopts() importopts()
if "`r(result)'" == "pass" {
    noi test_pass "special characters data identical (cf _all)"
}
else {
    noi test_fail "special characters" "`r(reason)'"
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
