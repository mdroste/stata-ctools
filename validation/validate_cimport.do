/*******************************************************************************
 * validate_cimport.do
 *
 * Comprehensive validation tests for cimport delimited vs import delimited
 * Tests all import options: delimiters, varnames, case, bindquotes, stripquotes
 ******************************************************************************/

do "validate_setup.do"

quietly {

di as text ""
di as text "======================================================================"
di as text "              CIMPORT DELIMITED VALIDATION TEST SUITE"
di as text "======================================================================"

capture mkdir "temp"

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

* Tab-delimited file
file open fh using "temp/tabfile.tsv", write replace
file write fh "id	name	score	group" _n
file write fh "1	John	85.5	X" _n
file write fh "2	Jane	92.3	Y" _n
file write fh "3	Bob	78.1	X" _n
file write fh "4	Alice	95.7	Z" _n
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
 * SECTION 2: Basic comma-delimited import
 ******************************************************************************/
noi print_section "Basic Comma-Delimited Import"

import delimited using "temp/basic.csv", clear
local stata_n = _N
local stata_k = c(k)

cimport delimited using "temp/basic.csv", clear
local cimport_n = _N
local cimport_k = c(k)

if `stata_n' == `cimport_n' & `stata_k' == `cimport_k' {
    noi test_pass "basic CSV dimensions match (N=`stata_n', K=`stata_k')"
}
else {
    noi test_fail "basic CSV" "N: `stata_n' vs `cimport_n', K: `stata_k' vs `cimport_k'"
}

/*******************************************************************************
 * SECTION 3: Tab-delimited import
 ******************************************************************************/
noi print_section "Tab-Delimited Import"

import delimited using "temp/tabfile.tsv", delimiters(tab) clear
local stata_n = _N
local stata_k = c(k)

cimport delimited using "temp/tabfile.tsv", delimiters(tab) clear
local cimport_n = _N
local cimport_k = c(k)

if `stata_n' == `cimport_n' & `stata_k' == `cimport_k' {
    noi test_pass "tab-delimited dimensions match"
}
else {
    noi test_fail "tab-delimited" "dimensions differ"
}

/*******************************************************************************
 * SECTION 4: Semicolon-delimited import
 ******************************************************************************/
noi print_section "Semicolon-Delimited Import"

import delimited using "temp/semicolon.csv", delimiters(";") clear
local stata_n = _N

cimport delimited using "temp/semicolon.csv", delimiters(";") clear
local cimport_n = _N

if `stata_n' == `cimport_n' {
    noi test_pass "semicolon-delimited N matches"
}
else {
    noi test_fail "semicolon-delimited" "N differs"
}

/*******************************************************************************
 * SECTION 5: varnames() option
 ******************************************************************************/
noi print_section "varnames() Option"

* varnames(1)
import delimited using "temp/basic.csv", varnames(1) clear
local stata_n = _N

cimport delimited using "temp/basic.csv", varnames(1) clear
local cimport_n = _N

if `stata_n' == `cimport_n' {
    noi test_pass "varnames(1) N matches"
}
else {
    noi test_fail "varnames(1)" "N differs"
}

* varnames(nonames)
import delimited using "temp/noheader.csv", varnames(nonames) clear
local stata_n = _N
local stata_k = c(k)

cimport delimited using "temp/noheader.csv", varnames(nonames) clear
local cimport_n = _N
local cimport_k = c(k)

if `stata_n' == `cimport_n' & `stata_k' == `cimport_k' {
    noi test_pass "varnames(nonames) dimensions match"
}
else {
    noi test_fail "varnames(nonames)" "dimensions differ"
}

/*******************************************************************************
 * SECTION 6: case() option
 ******************************************************************************/
noi print_section "case() Option"

* case(preserve)
capture cimport delimited using "temp/casetest.csv", case(preserve) clear
if _rc == 0 {
    noi test_pass "case(preserve) accepted"
}
else {
    noi test_fail "case(preserve)" "returned error `=_rc'"
}

* case(lower)
capture cimport delimited using "temp/casetest.csv", case(lower) clear
if _rc == 0 {
    ds
    local varlist `r(varlist)'
    local lower_ok = 1
    foreach v of local varlist {
        if "`v'" != lower("`v'") local lower_ok = 0
    }
    if `lower_ok' {
        noi test_pass "case(lower) produces lowercase vars"
    }
    else {
        noi test_fail "case(lower)" "vars not all lowercase"
    }
}
else {
    noi test_fail "case(lower)" "returned error `=_rc'"
}

* case(upper)
capture cimport delimited using "temp/casetest.csv", case(upper) clear
if _rc == 0 {
    ds
    local varlist `r(varlist)'
    local upper_ok = 1
    foreach v of local varlist {
        if "`v'" != upper("`v'") local upper_ok = 0
    }
    if `upper_ok' {
        noi test_pass "case(upper) produces uppercase vars"
    }
    else {
        noi test_fail "case(upper)" "vars not all uppercase"
    }
}
else {
    noi test_fail "case(upper)" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 7: Missing values handling
 ******************************************************************************/
noi print_section "Missing Values"

import delimited using "temp/missing.csv", clear
count if missing(x)
local stata_miss = r(N)

cimport delimited using "temp/missing.csv", clear
count if missing(x)
local cimport_miss = r(N)

if `stata_miss' == `cimport_miss' {
    noi test_pass "missing value count matches"
}
else {
    noi test_fail "missing values" "counts differ: `stata_miss' vs `cimport_miss'"
}

/*******************************************************************************
 * SECTION 8: Quoted fields
 ******************************************************************************/
noi print_section "Quoted Fields"

import delimited using "temp/quoted.csv", clear
local stata_n = _N

cimport delimited using "temp/quoted.csv", clear
local cimport_n = _N

if `stata_n' == `cimport_n' {
    noi test_pass "quoted fields N matches"
}
else {
    noi test_fail "quoted fields" "N differs"
}

/*******************************************************************************
 * SECTION 9: bindquotes() option
 ******************************************************************************/
noi print_section "bindquotes() Option"

capture cimport delimited using "temp/quoted.csv", bindquotes(strict) clear
if _rc == 0 {
    noi test_pass "bindquotes(strict) accepted"
}
else {
    noi test_fail "bindquotes(strict)" "returned error `=_rc'"
}

capture cimport delimited using "temp/quoted.csv", bindquotes(loose) clear
if _rc == 0 {
    noi test_pass "bindquotes(loose) accepted"
}
else {
    noi test_fail "bindquotes(loose)" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 10: stripquotes option
 ******************************************************************************/
noi print_section "stripquotes Option"

capture cimport delimited using "temp/quoted.csv", stripquotes clear
if _rc == 0 {
    noi test_pass "stripquotes option accepted"
}
else {
    noi test_fail "stripquotes" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 11: Large file import
 ******************************************************************************/
noi print_section "Large File Import"

import delimited using "temp/large.csv", clear
local stata_n = _N
local stata_k = c(k)

cimport delimited using "temp/large.csv", clear
local cimport_n = _N
local cimport_k = c(k)

if `stata_n' == `cimport_n' & `stata_k' == `cimport_k' {
    noi test_pass "large file (5000 rows) dimensions match"
}
else {
    noi test_fail "large file" "dimensions differ"
}

/*******************************************************************************
 * SECTION 12: verbose/fast options
 ******************************************************************************/
noi print_section "verbose/fast Options"

capture cimport delimited using "temp/basic.csv", verbose clear
if _rc == 0 {
    noi test_pass "verbose option accepted"
}
else {
    noi test_fail "verbose option" "returned error `=_rc'"
}

capture cimport delimited using "temp/basic.csv", fast clear
if _rc == 0 {
    noi test_pass "fast option accepted"
}
else {
    noi test_fail "fast option" "returned error `=_rc'"
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
 * SECTION 14: Value comparisons
 ******************************************************************************/
noi print_section "Value Comparisons"

import delimited using "temp/basic.csv", clear
local stata_id2 = id[2]
local stata_name2 = name[2]

cimport delimited using "temp/basic.csv", clear
local cimport_id2 = id[2]
local cimport_name2 = name[2]

if `stata_id2' == `cimport_id2' {
    noi test_pass "numeric values match"
}
else {
    noi test_fail "numeric values" "differ"
}

if "`stata_name2'" == "`cimport_name2'" {
    noi test_pass "string values match"
}
else {
    noi test_fail "string values" "differ"
}

/*******************************************************************************
 * SECTION 15: Round-trip tests
 ******************************************************************************/
noi print_section "Round-Trip Tests"

sysuse auto, clear
export delimited using "temp/auto_rt.csv", replace

import delimited using "temp/auto_rt.csv", clear
local stata_n = _N
local stata_k = c(k)

cimport delimited using "temp/auto_rt.csv", clear
local cimport_n = _N
local cimport_k = c(k)

if `stata_n' == `cimport_n' & `stata_k' == `cimport_k' {
    noi test_pass "auto dataset round-trip"
}
else {
    noi test_fail "round-trip" "dimensions differ"
}

sysuse census, clear
export delimited using "temp/census_rt.csv", replace

import delimited using "temp/census_rt.csv", clear
local stata_n = _N

cimport delimited using "temp/census_rt.csv", clear
local cimport_n = _N

if `stata_n' == `cimport_n' {
    noi test_pass "census dataset round-trip"
}
else {
    noi test_fail "census round-trip" "N differs"
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
