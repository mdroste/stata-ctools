* Reference tests require SSC rangejoin 1.1.3 and rangestat.
do "validation/validate_setup.do"
which rangejoin
quietly {
clear
set seed 62083
set obs 700
gen long uid=_n
gen double key=floor(runiform()*25)
replace key=.a if mod(_n,53)==0
gen double g=mod(_n,7)
replace g=.a if mod(_n,17)==0
replace g=.b if mod(_n,19)==0
gen str4 s="g"+string(mod(_n,3))
replace s="" if mod(_n,23)==0
gen str9 text="using"+string(_n)
gen byte flag=mod(_n,2)
label define flaglabel 0 "no" 1 "yes"
label values flag flaglabel
label variable text "Using text label"
format key %9.2f
tempfile u m ref
save `u'
clear
set obs 180
gen long mid=_n
gen double key=floor(runiform()*30)
replace key=.b if mod(_n,37)==0
gen double lo=key-1
gen double hi=key+2
replace lo=.a if mod(_n,11)==0
replace hi=.z if mod(_n,13)==0
replace lo=50 in 1
gen double g=mod(_n,7)
replace g=.a if mod(_n,17)==0
replace g=.b if mod(_n,19)==0
gen str4 s="g"+string(mod(_n,3))
replace s="" if mod(_n,23)==0
gen str9 text="master"+string(_n)
label variable text "Master text label"
save `m'
foreach groups in "" "g" "s" "g s" {
    local byopt
    if "`groups'" != "" local byopt by(`groups')
    foreach bounds in "lo hi" "-2 3" ". ." ". 2" "-1 ." "hi lo" {
        foreach opts in "" "keepusing(uid text flag) prefix(u_)" "all suffix(_other)" {
            use `m', clear
            rangejoin key `bounds' using `u', `byopt' `opts'
            save `ref', replace
            use `m', clear
            crangejoin key `bounds' using `u', `byopt' `opts' threads(1)
            cf _all using `ref'
            test_pass "rangejoin `groups' `bounds' `opts'"
        }
    }
}
* No key in master: numeric bounds are absolute.
use `m', clear
drop key
save `m', replace
rangejoin key 2 4 using `u', by(g s)
save `ref', replace
use `m', clear
crangejoin key 2 4 using `u', by(g s) threads(4) verbose
cf _all using `ref'
assert "`: variable label text_U'"=="Using text label"
assert "`: value label flag'"=="flaglabel"
assert "`: format key'"=="%9.2f"
test_pass "absolute bounds and variable/value labels"
* Zero matches retain every master row in order.
use `m', clear
crangejoin key 100 200 using `u'
assert _N==180 & mid==_n & missing(uid)
test_pass "unmatched rows preserved"
* Invalid calls restore the master exactly.
use `m', clear
gen text_U="collision"
save `m', replace
capture crangejoin key 0 2 using `u'
assert _rc==110
cf _all using `m'
test_pass "renaming conflict rolls back"
capture crangejoin key 5 2 using `u'
assert _rc==2000
cf _all using `m'
test_pass "all reversed intervals"
capture crangejoin absent 0 2 using `u'
assert _rc==111
cf _all using `m'
test_pass "missing using key rolls back"
* Parallel search and output, duplicate using keys retain input tie order.
clear
set obs 21000
gen long mid=_n
gen double key=mod(_n,25)
rangejoin key 0 0 using `u'
save `ref', replace
clear
set obs 21000
gen long mid=_n
gen double key=mod(_n,25)
crangejoin key 0 0 using `u', threads(4) verbose
cf _all using `ref'
test_pass "parallel join and stable tie ordering"
/*******************************************************************************
 * Rollback of failed joins
 *
 * A failed or interrupted crangejoin leaves the master as it was without a
 * snapshot of the dataset (Stata 16 or later); option snapshot preserves it
 * instead and option fast undoes nothing (see validation/rollback_helpers.do).
 * The hook $CTOOLS_crangejoin_fault fails the join:
 *   schema   with r(1) after the using variables, value labels, notes and
 *            characteristics were appended and the output rows added
 *   store    inside the plugin after half the output rows were written, with
 *            the memory error it simulates, r(920)
 *   post     with r(1) after the plugin wrote the whole join
 ******************************************************************************/
print_section "Rollback of Failed Joins"
capture confirm file "validation/rollback_helpers.do"
if !_rc do "validation/rollback_helpers.do"
else do "rollback_helpers.do"

clear
set obs 90
gen int grp = mod(_n, 7) + 1
gen str2 sgrp = "g" + string(mod(_n, 7) + 1)
gen double x = _n * 3
gen double val = _n / 4
gen str5 s = "m" + string(_n)
label define rb_grp 1 "first" 2 "second"
label values grp rb_grp
label data "Using data"
note: using dataset note
note val: using note on val
char _dta[rb_uchar] using characteristic
tempfile rb_using
save `rb_using'
clear
set seed 3141
set obs 400
gen long id = _n
gen int grp = mod(_n * 3, 8) + 1
gen str2 sgrp = "g" + string(grp)
gen double x = runiform() * 270
gen float val = rnormal()
gen str4 tag = "t" + string(mod(_n, 50))
label define rb_tag 1 "one"
label values grp rb_tag
label variable x "master key"
label data "crangejoin rollback master"
note: master dataset note
note x: master note on x
char _dta[rb_char] master characteristic
char val[rb_char] val characteristic
sort grp id
tempfile rb_copy rb_out
save `rb_copy'

foreach fault in schema store post {
    local rc = cond("`fault'" == "store", 920, 1)
    _ctools_rb_case, cmd(crangejoin x -10 10 using `rb_using', by(grp) keepusing(val s)) ///
        global(CTOOLS_crangejoin_fault) fault(`fault') rc(`rc') testname("[rollback] `fault'")
}
_ctools_rb_case, cmd(crangejoin x -5 5 using `rb_using', by(sgrp) all prefix(u_)) ///
    global(CTOOLS_crangejoin_fault) fault(post) testname("[rollback] string by(), all variables")

* Options snapshot and fast give the default's results; together they are refused
use `rb_copy', clear
crangejoin x -10 10 using `rb_using', by(grp)
save `rb_out', replace
foreach mode in snapshot fast {
    use `rb_copy', clear
    crangejoin x -10 10 using `rb_using', by(grp) `mode'
    capture cf _all using `rb_out', all
    if !_rc test_pass "[rollback] option `mode' gives the default results"
    else test_fail "[rollback] option `mode' gives the default results" "results differ"
}
use `rb_copy', clear
capture crangejoin x -10 10 using `rb_using', by(grp) snapshot fast
local rc = _rc
capture cf _all using `rb_copy', all
if `rc' == 198 & !_rc test_pass "[rollback] snapshot and fast together refused"
else test_fail "[rollback] snapshot and fast together refused" "rc=`rc', data changed=`=_rc!=0'"

}
print_summary "crangejoin"
global CTOOLS_COMPONENT_COMPLETE "crangejoin"
