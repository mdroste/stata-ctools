gettoken __case 0 : 0
if "`__case'" == "astra_nov6_followup" {
do "validation/support/validate_setup.do"
set seed 615
clear
set obs 100
gen double x = _n
gen double z = rnormal()
gen double d = z+rnormal()
gen double y = x+d+rnormal()
gen g = mod(_n,5)
gen byte eligible = x<50.5
summarize x, meanonly
creghdfe y x if x<r(mean), absorb(g)
assert e(sample)==eligible
test_pass "ASTRA-17 OLS incoming r sample"
summarize x, meanonly
civreghdfe y x (d=z) if x<r(mean), absorb(g)
assert e(sample)==eligible
test_pass "ASTRA-17 IV incoming r sample"
summarize x, meanonly
crangestat (mean) y if x<r(mean), interval(x . .)
assert missing(y_mean)==!eligible
test_pass "ASTRA-17 range statistics incoming r sample"
clear
set obs 1
gen id=1
gen x=2
gen z=3
label define omitted 3 "omitted-variable label"
label values z omitted
tempfile usingdata
save `usingdata'
foreach noframes in 0 1 {
    clear
    set obs 1
    gen id=1
    global CTOOLS_cmerge_noframes `noframes'
    cmerge 1:1 id using `usingdata', keepusing(x)
    global CTOOLS_cmerge_noframes
    assert "`: label omitted 3'"=="omitted-variable label"
    test_pass "ASTRA-05 keepusing label definitions noframes=`noframes'"
}
clear
set obs 1
gen str3 s="a:b"
local delimiters ""
forvalues j=1/2046 {
    local delimiters "`delimiters' :"
    if `j'>=2044 {
        split s, parse(`delimiters') gen(ref)
        csplit s, parse(`delimiters') gen(got)
        assert ref1==got1 & ref2==got2
        drop ref* got*
        test_pass "ASTRA-42 `j' delimiter entries"
    }
}
clear
csample, count(0)
assert _csample_nobs==0 & _csample_kept==0 & _csample_ngroups==0
test_pass "ASTRA-45 zero-row count shortcut"
clear
set obs 100
gen g=mod(_n,5)
gen double x=rnormal()
gen double y=rpoisson(exp(.2*x))
gen double exposure=1
cppmlhdfe y x, absorb(g) d(saved_sum)
tempfile valid_ppml
save `valid_ppml'
foreach invalid in negative empty exposure {
    use `valid_ppml', clear
    cppmlhdfe y x, absorb(g)
    matrix previous_b=e(b)
    local qualifier ""
    local opts ""
    local expected_rc=198
    if "`invalid'"=="negative" replace y=-1 in 1
    if "`invalid'"=="empty" {
        local qualifier "if 0"
        local expected_rc=2000
    }
    if "`invalid'"=="exposure" {
        replace exposure=-1 in 1
        local opts "exposure(exposure)"
    }
    quietly datasignature
    local signature `"`r(datasignature)'"'
    capture cppmlhdfe y x `qualifier', absorb(g) `opts'
    assert _rc==`expected_rc'
    quietly datasignature
    assert `"`r(datasignature)'"'==`"`signature'"'
    assert mreldif(e(b),previous_b)==0
    test_pass "ASTRA-19 `invalid' exact error, data and previous estimates"
}
clear
set obs 2
gen str3 s="old"
cencode s, gen(n)
drop n
mata: __cencode_vl_names = 1
capture cencode s, gen(n)
assert _rc!=0
confirm new variable n
assert s=="old"
assert "`: label n 1'"=="old"
test_pass "ASTRA-39 label-backup failure aborts before mutation"
print_summary "astra_nov6_followup"
if $TESTS_FAILED exit 9
global CTOOLS_CASE_COMPLETE "astra_nov6_followup"
}
if "`__case'" == "astra_nov6_wrappers" {
do "validation/support/validate_setup.do"
set seed 611
clear
set obs 20
gen double x = _n
gen byte x_w = 99
label variable x_w "original output"
char x_w[audit] "keep me"
format x_w %8.0g
sort x
summarize x, meanonly
_ctools_undo record ASTRA_test_
if r(mean) == 10.5 test_pass "ASTRA-17 undo preserves incoming r()"
else test_fail "ASTRA-17 undo preserves incoming r()" "mean lost"
_ctools_undo clear ASTRA_test_
csort x, nosortedby
if "`: sortedby'" == "" & r(N) == 20 & "`r(cmd)'" == "csort" test_pass "ASTRA-02 shortcut results and nosortedby"
else test_fail "ASTRA-02 shortcut results and nosortedby" "marker or results incorrect"
foreach phase in outputs bounds store post {
    global CTOOLS_cwinsor_fault `phase'
    capture noisily cwinsor x, suffix(_w) replace
    local rc = _rc
    capture assert x_w == 99
    local bad = _rc
    if `rc' != 0 & !`bad' & "`: type x_w'" == "byte" & "`: variable label x_w'" == "original output" & "`: char x_w[audit]'" == "keep me" test_pass "ASTRA-34 `phase' restores existing target"
    else test_fail "ASTRA-34 `phase' restores existing target" "data or metadata changed"
    drop x_w
    gen byte x_w = 99
    label variable x_w "original output"
    char x_w[audit] "keep me"
}
macro drop CTOOLS_cwinsor_fault
foreach mode in "" snapshot fast {
    foreach op in < > <= == {
        summarize x, meanonly
        gen byte expected = x `op' r(mean)
        summarize x, meanonly
        capture noisily cwinsor x if x `op' r(mean), gen(got) `mode'
        local rc = _rc
        quietly count if expected
        if r(N) == 0 {
            if `rc' == 2000 test_pass "ASTRA-17 empty selection `op' `mode'"
            else test_fail "ASTRA-17 empty selection `op' `mode'" "wrong status"
        }
        else {
            capture assert !missing(got) == expected
            if !`rc' & !_rc test_pass "ASTRA-17 winsor selection `op' `mode'"
            else test_fail "ASTRA-17 winsor selection `op' `mode'" "wrong sample"
        }
        capture drop got
        drop expected
    }
}
tempfile base out
save `base'
foreach kind in spss sas sasxport5 sasxport8 dbase {
    foreach op in < > <= == {
        use `base', clear
        summarize x, meanonly
        gen byte expected = x `op' r(mean)
        quietly count if expected
        local want = r(N)
        summarize x, meanonly
        summarize x, meanonly
        capture noisily cexport `kind' x if x `op' r(mean) using "`out'", replace
        local rc = _rc
        if !`rc' & r(N) == `want' test_pass "ASTRA-10 `kind' selection `op'"
        else test_fail "ASTRA-10 `kind' selection `op'" "wrong exported count or status `rc'"
    }
}
print_summary "astra_nov6_wrappers"
if $TESTS_FAILED exit 9
global CTOOLS_CASE_COMPLETE "astra_nov6_wrappers"
}

if "`__case'" == "astra_nov6_contracts" {
do "validation/support/validate_setup.do"
set seed 612
clear
set obs 20
gen double x = _n
gen double y = x^2
gen double key = 21-_n
gen str8 payload = string(x)
tempfile original
save `original'
global CTOOLS_csort_fault post
capture noisily csort key, stream(1)
local rc = _rc
capture cf _all using `original', all
if `rc' == 1 & !_rc test_pass "ASTRA-01 streaming failure restores all columns"
else test_fail "ASTRA-01 streaming recovery" "status or data differs"
macro drop CTOOLS_csort_fault
csort key, stream(2) nosortedby
assert key == _n
assert x == 21-_n
assert payload == string(x)
test_pass "ASTRA-01 streaming mixed payload success after failure"
use `original', clear
foreach scale in 1 1e306 1e307 1e-300 {
    gen double w = `scale'
    cwinsor x [aw=w], gen(result) cuts(50 50)
    capture assert result == 10.5
    if !_rc test_pass "ASTRA-35 weighted median scale `scale'"
    else test_fail "ASTRA-35 weighted median scale `scale'" "expected 10.5"
    drop w result
}
foreach thread in -1 257 {
    capture csample 50, threads(`thread')
    if _rc == 198 test_pass "ASTRA-30 sample rejects threads `thread'"
    else test_fail "ASTRA-30 sample threads" "wrong status"
    capture crangestat (mean) x, interval(key -1 1) threads(`thread')
    if _rc == 198 test_pass "ASTRA-30 range rejects threads `thread'"
    else test_fail "ASTRA-30 range threads" "wrong status"
}
scalar _csample_kept = 999
csample 50 if 0, verbose
if _csample_kept == 0 & _csample_nobs == 0 & _csample_ngroups == 0 & _csample_threads == 0 test_pass "ASTRA-45 empty sample current counts"
else test_fail "ASTRA-45 empty sample" "stale counts"
foreach miss in . .a .z {
    crangestat (mean) result=x, interval(key `miss' `miss')
    assert result == 10.5
    drop result
}
test_pass "ASTRA-29 extended missing constant bounds"
crangestat (cov) x y, interval(key . .) verbose
assert !missing(cov_x)
test_pass "ASTRA-28 flexible-only verbose succeeds"
foreach tol in 0 -1 2 . {
    capture cqreg y x, tolerance(`tol')
    if _rc == 198 test_pass "ASTRA-21 rejects tolerance `tol'"
    else test_fail "ASTRA-21 tolerance `tol'" "wrong status"
}
capture cqreg y x, nopreprocess(-2)
if _rc == 198 test_pass "ASTRA-21 rejects negative preprocessing"
else test_fail "ASTRA-21 preprocessing" "wrong status"
* Subscript sizing is independent of the sampler's intentional RNG difference.
use `original', clear
bsample x[1] if x>5
local want = _N
use `original', clear
cbsample x[1] if x>5
if _N == `want' & _N == 6 test_pass "ASTRA-44 eligible-order sizing"
else test_fail "ASTRA-44 sizing" "wrong count"
* Force the strL/native path; counts describe categories before keep().
clear
set obs 2
gen id = _n+1
gen strL s = "using"
tempfile usingdata
save `usingdata'
foreach extra in "" "nogenerate" {
    clear
    set obs 2
    gen id = _n
    gen strL s = "master"
    cmerge 1:1 id using `usingdata', keep(match) `extra'
    if r(N)==1 & r(N_1)==1 & r(N_2)==1 & r(N_3)==1 & r(N_4)==0 & r(N_5)==0 & missing(r(time)) test_pass "ASTRA-04 fallback category counts `extra'"
    else test_fail "ASTRA-04 fallback `extra'" "wrong category counts"
}
print_summary "astra_nov6_contracts"
if $TESTS_FAILED exit 9
global CTOOLS_CASE_COMPLETE "astra_nov6_contracts"
}

if "`__case'" == "astra_nov6_bootstrap" {
do "validation/support/validate_setup.do"
set seed 613
clear
set obs 100
gen g = ceil(_n/10)
gen double x = rnormal()
gen double z = rnormal()
gen double d = z + rnormal()
gen double y = x + d + rnormal()
gen double yp = rpoisson(exp(.2*x+.1*d))
gen double __ctools_bsu1 = _n
gen str9 __ctools_bsu_user = "original"
gen double __ctools_wbv1 = -_n
gen str9 __ctools_wbv_user = "original"
label variable __ctools_bsu1 "retain label"
char __ctools_wbv1[audit] "retain char"
foreach cmd in creghdfe civreghdfe cppmlhdfe {
    local model "y x"
    if "`cmd'" == "civreghdfe" local model "y x (d=z)"
    if "`cmd'" == "cppmlhdfe" local model "yp x"
    foreach mode in pairs wild enum {
        local options "vce(bootstrap, reps(9) seed(71) nodots)"
        if "`mode'" == "wild" local options "vce(cluster g) wildbootstrap(reps(9) seed(71) noenum nodots)"
        if "`mode'" == "enum" local options "vce(cluster g) wildbootstrap(reps(1024) seed(71) nodots)"
        capture noisily `cmd' `model', absorb(g) `options'
        local rc = _rc
        capture assert __ctools_bsu1 == _n & __ctools_wbv1 == -_n & __ctools_bsu_user == "original" & __ctools_wbv_user == "original"
        if !`rc' & !_rc & "`: variable label __ctools_bsu1'" == "retain label" & "`: char __ctools_wbv1[audit]'" == "retain char" test_pass "ASTRA-16 `cmd' `mode' retains user variables"
        else test_fail "ASTRA-16 `cmd' `mode'" "rc=`rc' or user data changed"
    }
}
print_summary "astra_nov6_bootstrap"
if $TESTS_FAILED exit 9
global CTOOLS_CASE_COMPLETE "astra_nov6_bootstrap"
}

if "`__case'" == "astra_nov6_numerics" {
do "validation/support/validate_setup.do"
set seed 614
clear
set obs 30
gen double x = _n
gen double y = 2*x+mod(_n,3)
gen double q = cond(_n<=10, 2147483648, cond(_n<=20, 2^40, 2^50))
cbinscatter y x, xq(q) nograph linetype(none)
matrix bins=e(bindata)
assert bins[1,2]==2147483648 & bins[2,2]==2^40 & bins[3,2]==2^50
test_pass "ASTRA-22 exact external bin labels above INT_MAX"
* Covariate adjustment must still index bins by their internal dense IDs.
gen double z = mod(_n,7)
cbinscatter y x, xq(q) controls(z) nograph linetype(none)
matrix bins=e(bindata)
assert bins[1,2]==2147483648 & bins[2,2]==2^40 & bins[3,2]==2^50
test_pass "ASTRA-22 external labels with controls"
clear
set obs 5
gen byte treated = _n==5
gen double score = cond(_n==1,.1,cond(_n==2,.2,cond(_n==3,.8,cond(_n==4,.9,.85))))
gen double outcome = cond(_n==1,1e100,cond(_n==2,1e50,cond(_n==3,1,cond(_n==4,2,3))))
cpsmatch treated, pscore(score) outcome(outcome) radius caliper(.06)
assert _outcome[5]==1.5
test_pass "ASTRA-24 unrelated huge outcomes do not erase local mean"
clear
set obs 3
gen str8 s = "$12"
local ignore ""
forvalues j=1/255 {
    local ignore "`ignore'a"
}
local ignore "`ignore'$"
cdestring s, gen(n) ignore("`ignore'")
assert n==12
test_pass "ASTRA-36 complete ignore list beyond 255 bytes"
drop n
cdestring s if 0, gen(n) verbose
assert missing(n)
assert _cdestring_n_converted==0 & _cdestring_time_total==0
test_pass "ASTRA-37 empty verbose conversion"
* Export literal text safely, including a quote in the transport path.
gen byte code = 1
label define ASTRA 1 "O'Brien; sample"
label values code ASTRA
tempfile out
foreach kind in sasxport5 sasxport8 {
    local opts "vallabfile"
    if "`kind'"=="sasxport5" local opts "vallabfile(sascode)"
    cexport `kind' code using "`out'O'Brien.xpt", `opts' replace
    tempname fh
    file open `fh' using "`out'O'Brien.sas", read text
    local found_label 0
    local found_path 0
    file read `fh' line
    while !r(eof) {
        if strpos(`"`line'"', "O''Brien; sample") local found_label 1
        if strpos(`"`line'"', "O''Brien.xpt") local found_path 1
        file read `fh' line
    }
    file close `fh'
    assert `found_label' & `found_path'
    test_pass "ASTRA-11 SAS sidecar escaped `kind'"
    erase "`out'O'Brien.xpt"
    erase "`out'O'Brien.sas"
}
* The fit has raw K >= N, but only two independent columns including constant.
clear
set obs 10
gen double x = _n
gen double y = 2*x + mod(_n,3)
forvalues j=1/10 {
    gen double x`j' = x
}
qreg y x1-x10
matrix ref=e(b)
cqreg y x1-x10
assert abs(_b[x1]-ref[1,1])<1e-6
assert e(N)==10
test_pass "ASTRA-20 rank screen precedes sample-size rejection"
print_summary "astra_nov6_numerics"
if $TESTS_FAILED exit 9
global CTOOLS_CASE_COMPLETE "astra_nov6_numerics"
}
