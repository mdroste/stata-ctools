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
}
print_summary "crangejoin"
global CTOOLS_COMPONENT_COMPLETE "crangejoin"
