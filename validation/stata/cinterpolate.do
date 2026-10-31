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
