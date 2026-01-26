/*******************************************************************************
 * validate_cexport_excel_opts.do
 * Validation tests for cexport excel new options: cell(), missing(), keepcellfmt
 ******************************************************************************/

clear all
set more off

* Add build directory to adopath so ctools commands are found
adopath ++ "build"

* Create test data with missing values
sysuse auto, clear
replace price = . in 1/3
replace mpg = . in 5/7

* Test 1: Basic export (baseline)
di as text "Test 1: Basic export"
cexport excel using "test_excel_basic.xlsx", replace
assert r(N) == _N
di as result "  PASSED"

* Test 2: Export starting at cell B5
di as text "Test 2: Export starting at cell(B5)"
cexport excel make price mpg using "test_excel_cell.xlsx", cell(B5) replace
assert r(N) == _N
di as result "  PASSED"

* Test 3: Export with missing value replacement
di as text "Test 3: Export with missing(NA)"
cexport excel make price mpg using "test_excel_missing.xlsx", missing("NA") replace
assert r(N) == _N
di as result "  PASSED"

* Test 4: Export with missing(".") replacement
di as text "Test 4: Export with missing(.)"
cexport excel make price mpg using "test_excel_missing_dot.xlsx", missing(".") replace
assert r(N) == _N
di as result "  PASSED"

* Test 5: Combination - cell offset + missing value
di as text "Test 5: Combination cell(C3) + missing(N/A)"
cexport excel make price mpg using "test_excel_combo.xlsx", cell(C3) missing("N/A") replace
assert r(N) == _N
di as result "  PASSED"

* Test 6: keepcellfmt (preserves styles from existing file)
di as text "Test 6: keepcellfmt option"
* First create a file
cexport excel make price using "test_excel_keepcellfmt.xlsx", replace
* Then update it with keepcellfmt
cexport excel make price mpg using "test_excel_keepcellfmt.xlsx", keepcellfmt replace
assert r(N) == _N
di as result "  PASSED"

* Test 7: Cell reference validation - various formats
di as text "Test 7: Various cell references"
cexport excel make price using "test_excel_cell_A1.xlsx", cell(A1) replace
cexport excel make price using "test_excel_cell_Z1.xlsx", cell(Z1) replace
cexport excel make price using "test_excel_cell_AA1.xlsx", cell(AA1) replace
cexport excel make price using "test_excel_cell_AB100.xlsx", cell(AB100) replace
di as result "  PASSED"

* Cleanup
foreach f in test_excel_basic test_excel_cell test_excel_missing test_excel_missing_dot ///
             test_excel_combo test_excel_keepcellfmt test_excel_cell_A1 test_excel_cell_Z1 ///
             test_excel_cell_AA1 test_excel_cell_AB100 {
    cap erase "`f'.xlsx"
}

di as text _n "All cexport excel option tests PASSED"
