* cinterpolate: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "cinterpolate"

if "`__case'" == "cinterpolate" {
* Reference comparisons for ctools linear interpolation.
do "validation/support/validate_setup.do"
quietly {
clear
set seed 40217
set obs 30000
gen long id = _n
gen double x = floor(runiform()*45)-10
gen double y = rnormal()*100
replace y = . if runiform()<.65
gen byte g = mod(_n,31)
replace g = .a if mod(_n,101)==0
replace g = .b if mod(_n,103)==0
gen str6 s = "g" + string(mod(_n,7))
replace s = "" if mod(_n,19)==0
replace x = .z if mod(_n,41)==0
replace y = .b if mod(_n,37)==0
foreach groups in "" "g" "s" "g s" {
    local byopt
    if "`groups'" != "" local byopt by(`groups')
    foreach ep in "" "epolate" {
        foreach sample in "" "if mod(id,3)!=0" "in 11/211" {
            sort id
            ipolate y x `sample', gen(ref) `byopt' `ep'
            * Restore original order so in refers to the same rows.
            sort id
            cinterpolate y x `sample', gen(got) `byopt' `ep' threads(1)
            assert_var_equal got ref $DEFAULT_SIGFIGS "ipolate `groups' `ep' `sample'"
            assert id == _n
            drop ref got
        }
    }
}
by g s, sort: ipolate y x, gen(ref) epolate
by g s, sort: cinterpolate y x, gen(got) epolate threads(4) verbose
assert_var_equal got ref $DEFAULT_SIGFIGS "by prefix and parallel groups"
drop ref got
* Known answers, duplicated x, all-missing and one-anchor groups.
clear
input double(x y g expected)
0 . 1 -10
1 0 1 0
2 . 1 10
3 15 1 20
3 25 1 20
4 . 1 30
0 . 2 .
1 5 2 5
2 . 2 .
1 . 3 .
2 . 3 .
. 9 1 .
end
cinterpolate y x, gen(got) by(g) epolate
assert_var_equal got expected $DEFAULT_SIGFIGS "known interpolation and endpoint extrapolation"
drop got
cinterpolate y x if 0, gen(got)
assert missing(got)
test_pass "empty selected sample"
capture cinterpolate y x, gen(got)
assert _rc == 110
test_pass "existing destination rejected"
capture by g, sort: cinterpolate y x, gen(bad) by(g)
assert _rc == 190
capture confirm variable bad
assert _rc == 111
test_pass "by option conflict leaves no output"
capture cinterpolate y x, gen(bad) threads(-1)
assert _rc == 198
test_pass "invalid thread count"
clear
set obs 0
gen double y=.
gen double x=.
cinterpolate y x, gen(z)
assert _N==0
test_pass "zero-observation data"
}
quietly {
* _ct_output runs a command and returns its output (r(text), lines joined
* by |) and return code (r(rc)). _ct_same_output runs a native and a ctools
* command on copies of the data in memory and passes when the output, return
* code, variable names and types, and (unless "nodata") values agree.
capture program drop _ct_output
program define _ct_output, rclass
    tempfile f
    tempname fh
    capture log close _ctout
    quietly log using `"`f'"', text replace name(_ctout)
    capture noisily `0'
    local rc = _rc
    quietly log close _ctout
    local text
    file open `fh' using `"`f'"', read text
    file read `fh' line
    while !r(eof) {
        if `"`macval(line)'"' != "" & substr(`"`macval(line)'"', 1, 5) != "r; t=" {
            local text `"`macval(text)'`macval(line)'|"'
        }
        file read `fh' line
    }
    file close `fh'
    return local text `"`macval(text)'"'
    return scalar rc = `rc'
end
capture program drop _ct_same_output
program define _ct_same_output
    args name native ctools nodata
    tempfile want
    preserve
    _ct_output `native'
    local want_text `"`r(text)'"'
    local want_rc = r(rc)
    local want_vars
    local want_types
    forvalues j = 1/`c(k)' {
        mata: st_local("v", st_varname(`j'))
        local want_vars `want_vars' `v'
        local want_types `want_types' `: type `v''
    }
    capture quietly save `"`want'"', emptyok
    restore, preserve
    _ct_output `ctools'
    local got_text `"`r(text)'"'
    local got_rc = r(rc)
    local got_vars
    local got_types
    forvalues j = 1/`c(k)' {
        mata: st_local("v", st_varname(`j'))
        local got_vars `got_vars' `v'
        local got_types `got_types' `: type `v''
    }
    local same_data = 1
    if "`nodata'" == "" & "`got_vars'" != "" {
        capture cf _all using `"`want'"'
        local same_data = !_rc
    }
    restore
    if "`nodata'" != "" {
        local want_vars
        local got_vars
        local want_types
        local got_types
    }
    if `"`macval(want_text)'"' == `"`macval(got_text)'"' & `want_rc' == `got_rc' & ///
        "`want_vars'" == "`got_vars'" & "`want_types'" == "`got_types'" & `same_data' {
        test_pass "`name'"
    }
    else {
        noisily di as error `"`name': ctools [`macval(got_text)'] rc `got_rc' vars `got_vars' (`got_types')"'
        noisily di as error `"`name': native [`macval(want_text)'] rc `want_rc' vars `want_vars' (`want_types') same data `same_data'"'
        test_fail "`name'" "output, return code, variables or data differ from the native command"
    }
end
* IPOL-01, IPOL-02: r(N) and the "(# missing value(s) generated)" message
clear
input double(x y g)
1 1 1
2 . 1
3 3 1
. 4 1
1 . 2
2 5 2
3 . 2
end
quietly summarize x
foreach pair in `"y x, gen(z)"' `"y x, gen(z) epolate"' `"y x, gen(z) by(g)"' ///
    `"y x if g == 1, gen(z)"' `"y x in 1/3, gen(z)"' `"y x if 0, gen(z)"' `"y x, gen(z) by(g) epolate"' {
    _ct_same_output "IPOL-02 message: `pair'" `"ipolate `pair'"' `"cinterpolate `pair'"'
    preserve
    ipolate `pair'
    local want = r(N)
    restore, preserve
    quietly summarize g
    cinterpolate `pair'
    local got = r(N)
    restore
    if `got' == `want' test_pass "IPOL-01 r(N): `pair'"
    else test_fail "IPOL-01 r(N): `pair'" "r(N)=`got', ipolate `want'"
}
preserve
keep in 1/3
replace x = . in 2
quietly summarize g
ipolate y x, gen(z)
local want = r(N)
cinterpolate y x, gen(z2)
if r(N) == `want' & `want' == 1 test_pass "IPOL-01 r(N) one missing value"
else test_fail "IPOL-01 r(N) one missing value" "r(N)=`r(N)', ipolate `want'"
restore
clear
set obs 0
gen double y = .
gen double x = .
cinterpolate y x, gen(z)
if r(N) == 0 test_pass "IPOL-01 r(N) with no observations"
else test_fail "IPOL-01 r(N) with no observations" "r(N)=`r(N)'"
}
* SSC differential oracle on unique x values, where its anchors are distinct.
* Double intermediates prevent the reference's default float storage from
* limiting the comparison; keep the suite's standard significant-figure rule.
quietly {
local oldtype = c(type)
set type double
clear
set seed 8722
set obs 3000
gen long id = _n
gen byte g = mod(_n, 7)
gen str4 s = "g" + string(g)
gen double x = _n + runiform()
gen double y = rnormal()*100
replace y = . if runiform() < .65
foreach groups in "" "g" "s" {
    local byopt
    if "`groups'" != "" local byopt by(`groups')
    foreach sample in "" "if mod(id,3)!=0" "in 101/2500" {
        sort id
        cipolate y x `sample', gen(ref) `byopt'
        sort id
        cinterpolate y x `sample', gen(got) `byopt' cubic
        assert_var_equal got ref $DEFAULT_SIGFIGS "SSC cubic `groups' `sample'"
        assert id == _n
        drop got ref
    }
}
set type `oldtype'
}
* Cubic polynomial oracle with irregular spacing, groups, and unsorted input.
quietly {
clear
set seed 87533
set obs 30000
gen long id = _n
gen byte g = mod(_n, 30)
gen str4 group = "g" + string(g)
gen double x = floor((_n-1)/30) + runiform()/2
gen double expected = 2 + x/10 - x^2/1000 + x^3/1000000 + g
gen double y = expected
replace y = . if mod(floor((_n-1)/30), 3) == 1
* Each group has two known points on each side only in the interior.
replace expected = . if missing(y) & (x < 2 | x > 998)
cinterpolate y x, gen(got) by(group) cubic threads(1)
assert_var_equal got expected $DEFAULT_SIGFIGS "cubic polynomial on irregular grouped x"
assert id == _n
cinterpolate y x, gen(parallel) by(group) cubic threads(4) verbose
assert got == parallel
bysort group: cinterpolate y x, gen(prefix) cubic
assert prefix == got
sort id
assert id == _n
test_pass "cubic threads, by prefix, and input order"
* Restriction must remove anchors as well as output observations.
replace expected = . if x < 300 | x > 700
replace expected = . if missing(y) & (x < 302 | x > 698)
cinterpolate y x if x >= 300 & x <= 700, gen(selected) by(group) cubic
assert_var_equal selected expected $DEFAULT_SIGFIGS "cubic selected anchors"
cinterpolate y x in 9001/21000, gen(ranged) by(group) cubic
assert_var_equal ranged expected $DEFAULT_SIGFIGS "cubic in selection"
* Duplicate anchors are averaged; repeated missing targets are both filled.
clear
input double(x y g expected)
-1 . 1 .
0 0 1 0
1 . 1 .
2 6 1 8
2 10 1 8
2 . 1 8
3 . 1 27
3 .a 1 27
4 64 1 64
5 . 1 .
6 216 1 216
7 . 1 .
. 10 1 .
.a 20 1 .
0 . 2 .
1 .a 2 .
0 7 3 7
1 . 3 .
0 0 .a 0
2 8 .a 8
3 . .a 27
4 64 .a 64
6 216 .a 216
3 . .b .
0 0 4 0
1 . 4 .
2 0 4 0
3 . 4 0
4 0 4 0
5 . 4 -1
6 0 4 0
7 . 4 .
8 16 4 16
end
cinterpolate y x, gen(got) by(g) cubic
assert_var_equal got expected $DEFAULT_SIGFIGS "cubic ties, sparse groups, and extended missing values"
* Translate and rescale x without changing the interpolant.
replace x = 1e12 + x/8 if x < .
cinterpolate y x, gen(shifted) by(g) cubic
assert_var_equal shifted expected $DEFAULT_SIGFIGS "cubic translated x"
replace x = (x - 1e12)*8e100 if x < .
cinterpolate y x, gen(scaled) by(g) cubic
assert_var_equal scaled expected $DEFAULT_SIGFIGS "cubic large x scale"
cinterpolate y x if 0, gen(empty) cubic
assert missing(empty)
capture cinterpolate y x, gen(bad) cubic epolate
assert _rc == 198
capture confirm variable bad
assert _rc == 111
test_pass "cubic rejects extrapolation without creating output"
clear
set obs 0
gen double x = .
gen double y = .
cinterpolate y x, gen(empty) cubic
assert _N == 0
test_pass "cubic empty data"
}
print_summary "cinterpolate"
global CTOOLS_COMPONENT_COMPLETE "cinterpolate"
    global CTOOLS_CASE_COMPLETE "cinterpolate"
    exit
}

di as error "Unknown cinterpolate case: `__case'"
exit 198
