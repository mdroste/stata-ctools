/*******************************************************************************
 * validate_cexport.do
 *
 * Comprehensive validation tests for cexport delimited vs export delimited
 * Tests all export options: delimiter, novarnames, quote, nolabel, if/in
 ******************************************************************************/

do "validate_setup.do"

quietly {

di as text ""
di as text "======================================================================"
di as text "              CEXPORT DELIMITED VALIDATION TEST SUITE"
di as text "======================================================================"

capture mkdir "temp"

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

export delimited using "temp/stata.csv", replace
cexport delimited using "temp/cexport.csv", replace

import delimited using "temp/stata.csv", clear
local stata_n = _N
local stata_k = c(k)

import delimited using "temp/cexport.csv", clear
local cexport_n = _N
local cexport_k = c(k)

if `stata_n' == `cexport_n' & `stata_k' == `cexport_k' {
    noi test_pass "basic export dimensions match"
}
else {
    noi test_fail "basic export" "dimensions differ"
}

/*******************************************************************************
 * SECTION 3: Tab delimiter
 ******************************************************************************/
noi print_section "Tab Delimiter"

sysuse auto, clear

export delimited using "temp/stata_tab.tsv", delimiter(tab) replace
cexport delimited using "temp/cexport_tab.tsv", delimiter(tab) replace

import delimited using "temp/stata_tab.tsv", delimiters(tab) clear
local stata_n = _N

import delimited using "temp/cexport_tab.tsv", delimiters(tab) clear
local cexport_n = _N

if `stata_n' == `cexport_n' {
    noi test_pass "tab delimiter N matches"
}
else {
    noi test_fail "tab delimiter" "N differs"
}

/*******************************************************************************
 * SECTION 4: Semicolon delimiter
 ******************************************************************************/
noi print_section "Semicolon Delimiter"

sysuse auto, clear

export delimited using "temp/stata_semi.csv", delimiter(";") replace
cexport delimited using "temp/cexport_semi.csv", delimiter(";") replace

import delimited using "temp/stata_semi.csv", delimiters(";") clear
local stata_n = _N

import delimited using "temp/cexport_semi.csv", delimiters(";") clear
local cexport_n = _N

if `stata_n' == `cexport_n' {
    noi test_pass "semicolon delimiter N matches"
}
else {
    noi test_fail "semicolon delimiter" "N differs"
}

/*******************************************************************************
 * SECTION 5: novarnames option
 ******************************************************************************/
noi print_section "novarnames Option"

sysuse auto, clear

export delimited using "temp/stata_novar.csv", novarnames replace
cexport delimited using "temp/cexport_novar.csv", novarnames replace

import delimited using "temp/stata_novar.csv", varnames(nonames) clear
local stata_n = _N

import delimited using "temp/cexport_novar.csv", varnames(nonames) clear
local cexport_n = _N

if `stata_n' == `cexport_n' {
    noi test_pass "novarnames N matches"
}
else {
    noi test_fail "novarnames" "N differs"
}

/*******************************************************************************
 * SECTION 6: quote option
 ******************************************************************************/
noi print_section "quote Option"

sysuse auto, clear

export delimited using "temp/stata_quote.csv", quote replace
cexport delimited using "temp/cexport_quote.csv", quote replace

import delimited using "temp/stata_quote.csv", clear
local stata_n = _N

import delimited using "temp/cexport_quote.csv", clear
local cexport_n = _N

if `stata_n' == `cexport_n' {
    noi test_pass "quote option N matches"
}
else {
    noi test_fail "quote option" "N differs"
}

/*******************************************************************************
 * SECTION 7: nolabel option
 ******************************************************************************/
noi print_section "nolabel Option"

sysuse auto, clear

* With labels (default)
export delimited using "temp/stata_label.csv", replace
cexport delimited using "temp/cexport_label.csv", replace

import delimited using "temp/stata_label.csv", clear
local stata_n = _N

import delimited using "temp/cexport_label.csv", clear
local cexport_n = _N

if `stata_n' == `cexport_n' {
    noi test_pass "with labels N matches"
}
else {
    noi test_fail "with labels" "N differs"
}

* Without labels
sysuse auto, clear

export delimited using "temp/stata_nolabel.csv", nolabel replace
cexport delimited using "temp/cexport_nolabel.csv", nolabel replace

import delimited using "temp/stata_nolabel.csv", clear
local stata_n = _N

import delimited using "temp/cexport_nolabel.csv", clear
local cexport_n = _N

if `stata_n' == `cexport_n' {
    noi test_pass "nolabel N matches"
}
else {
    noi test_fail "nolabel" "N differs"
}

/*******************************************************************************
 * SECTION 8: Variable selection
 ******************************************************************************/
noi print_section "Variable Selection"

sysuse auto, clear

export delimited make price mpg using "temp/stata_vars.csv", replace
cexport delimited make price mpg using "temp/cexport_vars.csv", replace

import delimited using "temp/stata_vars.csv", clear
local stata_k = c(k)

import delimited using "temp/cexport_vars.csv", clear
local cexport_k = c(k)

if `stata_k' == `cexport_k' & `stata_k' == 3 {
    noi test_pass "variable selection (3 vars)"
}
else {
    noi test_fail "variable selection" "K differs or not 3"
}

/*******************************************************************************
 * SECTION 9: if condition
 ******************************************************************************/
noi print_section "if Condition"

sysuse auto, clear

export delimited using "temp/stata_if.csv" if foreign == 1, replace
cexport delimited using "temp/cexport_if.csv" if foreign == 1, replace

import delimited using "temp/stata_if.csv", clear
local stata_n = _N

import delimited using "temp/cexport_if.csv", clear
local cexport_n = _N

if `stata_n' == `cexport_n' {
    noi test_pass "if condition N matches (foreign==1)"
}
else {
    noi test_fail "if condition" "N differs"
}

sysuse auto, clear

export delimited using "temp/stata_if2.csv" if price > 10000, replace
cexport delimited using "temp/cexport_if2.csv" if price > 10000, replace

import delimited using "temp/stata_if2.csv", clear
local stata_n = _N

import delimited using "temp/cexport_if2.csv", clear
local cexport_n = _N

if `stata_n' == `cexport_n' {
    noi test_pass "if condition N matches (price>10000)"
}
else {
    noi test_fail "if condition" "N differs"
}

/*******************************************************************************
 * SECTION 10: in condition
 ******************************************************************************/
noi print_section "in Condition"

sysuse auto, clear

export delimited using "temp/stata_in.csv" in 1/20, replace
cexport delimited using "temp/cexport_in.csv" in 1/20, replace

import delimited using "temp/stata_in.csv", clear
local stata_n = _N

import delimited using "temp/cexport_in.csv", clear
local cexport_n = _N

if `stata_n' == `cexport_n' & `stata_n' == 20 {
    noi test_pass "in 1/20 N=20"
}
else {
    noi test_fail "in condition" "N differs or not 20"
}

sysuse auto, clear

export delimited using "temp/stata_in2.csv" in 30/50, replace
cexport delimited using "temp/cexport_in2.csv" in 30/50, replace

import delimited using "temp/stata_in2.csv", clear
local stata_n = _N

import delimited using "temp/cexport_in2.csv", clear
local cexport_n = _N

if `stata_n' == `cexport_n' & `stata_n' == 21 {
    noi test_pass "in 30/50 N=21"
}
else {
    noi test_fail "in 30/50" "N differs or not 21"
}

/*******************************************************************************
 * SECTION 11: Combined if and in
 ******************************************************************************/
noi print_section "Combined if and in"

sysuse auto, clear

export delimited using "temp/stata_ifin.csv" if price > 5000 in 1/50, replace
cexport delimited using "temp/cexport_ifin.csv" if price > 5000 in 1/50, replace

import delimited using "temp/stata_ifin.csv", clear
local stata_n = _N

import delimited using "temp/cexport_ifin.csv", clear
local cexport_n = _N

if `stata_n' == `cexport_n' {
    noi test_pass "if and in combined N matches"
}
else {
    noi test_fail "if and in" "N differs"
}

/*******************************************************************************
 * SECTION 12: Census dataset
 ******************************************************************************/
noi print_section "Census Dataset"

sysuse census, clear

export delimited using "temp/stata_census.csv", replace
cexport delimited using "temp/cexport_census.csv", replace

import delimited using "temp/stata_census.csv", clear
local stata_n = _N
local stata_k = c(k)

import delimited using "temp/cexport_census.csv", clear
local cexport_n = _N
local cexport_k = c(k)

if `stata_n' == `cexport_n' & `stata_k' == `cexport_k' {
    noi test_pass "census dimensions match"
}
else {
    noi test_fail "census" "dimensions differ"
}

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

export delimited using "temp/stata_large.csv", replace
cexport delimited using "temp/cexport_large.csv", replace

import delimited using "temp/stata_large.csv", clear
local stata_n = _N

import delimited using "temp/cexport_large.csv", clear
local cexport_n = _N

if `stata_n' == `cexport_n' & `stata_n' == 50000 {
    noi test_pass "large dataset (50K) N matches"
}
else {
    noi test_fail "large dataset" "N differs or not 50K"
}

/*******************************************************************************
 * SECTION 14: Panel data (nlswork)
 ******************************************************************************/
noi print_section "Panel Data (nlswork)"

webuse nlswork, clear
keep in 1/10000

export delimited using "temp/stata_panel.csv", replace
cexport delimited using "temp/cexport_panel.csv", replace

import delimited using "temp/stata_panel.csv", clear
local stata_n = _N

import delimited using "temp/cexport_panel.csv", clear
local cexport_n = _N

if `stata_n' == `cexport_n' {
    noi test_pass "panel data N matches"
}
else {
    noi test_fail "panel data" "N differs"
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
 * SECTION 17: Edge cases
 ******************************************************************************/
noi print_section "Edge Cases"

* Single observation
clear
set obs 1
gen x = 42
gen str5 s = "test"

cexport delimited using "temp/single.csv", replace
import delimited using "temp/single.csv", clear
if _N == 1 {
    noi test_pass "single observation export"
}
else {
    noi test_fail "single observation" "N != 1"
}

* Single variable
clear
set obs 100
gen x = runiform()

export delimited x using "temp/singlevar_stata.csv", replace
cexport delimited x using "temp/singlevar_cexport.csv", replace

import delimited using "temp/singlevar_cexport.csv", clear
if c(k) == 1 {
    noi test_pass "single variable export"
}
else {
    noi test_fail "single variable" "K != 1"
}

* String with commas
clear
input id str40 name
1 "Smith, John"
2 "Doe, Jane"
3 "Regular Name"
end

export delimited using "temp/comma_stata.csv", replace
cexport delimited using "temp/comma_cexport.csv", replace

import delimited using "temp/comma_stata.csv", clear
local stata_n = _N

import delimited using "temp/comma_cexport.csv", clear
local cexport_n = _N

if `stata_n' == `cexport_n' & `stata_n' == 3 {
    noi test_pass "strings with commas"
}
else {
    noi test_fail "strings with commas" "N differs"
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
