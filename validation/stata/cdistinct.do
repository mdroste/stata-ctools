gettoken __case 0 : 0
if "`__case'" == "" local __case cdistinct
if "`__case'" == "cdistinct" {
do "validation/support/validate_setup.do"
quietly {
clear
set obs 5000
set seed 37197
gen long order = _n
gen double x = mod(_n,13)
replace x = .a if mod(_n,17)==0
replace x = .b if mod(_n,19)==0
gen str8 s = "s"+string(mod(_n,7))
replace s = "" if mod(_n,23)==0
gen byte g = mod(_n,3)
foreach keys in "x" "s" "x s" {
    foreach missing in "" "missing" {
        foreach joint in "" "joint" {
            foreach selection in "" "if g==1" "in 17/901" "if 0" {
                * Freeze the row restriction before the reference sorts each
                * variable. distinct 1.2.1 re-evaluates in after those sorts.
                tempvar selected
                mark `selected' `selection'
                distinct `keys' if `selected', `missing' `joint'
                local n = r(N)
                local d = r(ndistinct)
                drop `selected'
                cdistinct `keys' `selection', `missing' `joint' threads(1)
                assert r(N)==`n' & r(ndistinct)==`d'
                assert order==_n
            }
        }
    }
}
test_pass "distinct oracle: numeric, string, joint, missing, if and in"
cdistinct x s, missing joint threads(4) verbose min(100) max(1) abbrev(3)
local n=r(N)
local d=r(ndistinct)
cdistinct x s, missing joint
assert r(N)==`n' & r(ndistinct)==`d'
test_pass "display filters do not affect results; thread parity"
by g, sort: cdistinct x s, joint
by g: distinct x s, joint
sort order
cdistinct
assert r(N)==_N
capture cdistinct x, threads(-1)
assert _rc==198
clear
set obs 0
gen x=.
cdistinct x
assert r(N)==0 & r(ndistinct)==0
test_pass "by prefix, default varlist, validation and empty data"

* Compare the rendered table with distinct, including name widths and joint mode.
clear
set obs 12
gen short = mod(_n,3)
gen very_long_variable_name = mod(_n,7)
foreach option in "" "abbrev(5)" "abbrev(32)" "joint" "minimum(4)" "maximum(0)" "verbose" {
    foreach command in distinct cdistinct {
        tempfile displaylog
        tempname displayhandle filehandle
        log using `displaylog', text replace name(`displayhandle')
        local displayoption `option'
        if "`command'" == "distinct" & "`option'" == "verbose" local displayoption
        noisily `command' short very_long_variable_name, `displayoption'
        log close `displayhandle'
        file open `filehandle' using `displaylog', read text
        local table
        local saw_header 0
        file read `filehandle' line
        while r(eof)==0 {
            if strpos(`"`line'"', "cdistinct:") == 1 {
                assert !`saw_header'
            }
            else if strpos(`"`line'"', "|") | strpos(`"`line'"', "+") | ///
                strpos(`"`line'"', "Observations") | strpos(`"`line'"', "distinct") | ///
                regexm(`"`line'"', "^ +[0-9]+ +[0-9]+ *$") {
                local table `"`table'`line';"'
            }
            if strpos(`"`line'"', "Observations") local saw_header 1
            file read `filehandle' line
        }
        file close `filehandle'
        if "`command'" == "distinct" local reference `"`table'"'
        else assert `"`table'"' == `"`reference'"'
    }
}
test_pass "table alignment, abbreviation, joint display, hidden rows, and timing placement"

* Exercise hash growth, the sorting fallback, and the ordered-input path.
clear
set obs 10000
gen long order = _n
gen double x = mod(_n*7919,10000)
gen str24 s = "key"+string(x,"%012.0f")
foreach keys in "x" "s" "x s" {
    foreach cardinality in 10000 3000 100 {
        replace x = mod(order*7919,`cardinality')
        replace s = "key"+string(x,"%012.0f")
        foreach threads in 1 4 {
            cdistinct `keys', joint threads(`threads')
            assert r(N)==10000 & r(ndistinct)==`cardinality'
            assert order==_n
        }
    }
}
replace x = cond(_n<5000,0,.a)
replace s = cond(_n<5000,"", "a")
cdistinct x s, joint missing
assert r(N)==10000 & r(ndistinct)==2
cdistinct x s, joint
assert r(N)==0 & r(ndistinct)==0
test_pass "hash growth, high-cardinality fallback, and sorted missing keys"
}
noisily print_summary "cdistinct"
global CTOOLS_CASE_COMPLETE "cdistinct"
}
