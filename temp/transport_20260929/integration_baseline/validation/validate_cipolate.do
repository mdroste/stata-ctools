* Reference comparisons for ctools linear interpolation.
do "validation/validate_setup.do"
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
            cipolate y x `sample', gen(got) `byopt' `ep' threads(1)
            assert_var_equal got ref $DEFAULT_SIGFIGS "ipolate `groups' `ep' `sample'"
            assert id == _n
            drop ref got
        }
    }
}
by g s, sort: ipolate y x, gen(ref) epolate
by g s, sort: cipolate y x, gen(got) epolate threads(4) verbose
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
cipolate y x, gen(got) by(g) epolate
assert_var_equal got expected $DEFAULT_SIGFIGS "known interpolation and endpoint extrapolation"
drop got
cipolate y x if 0, gen(got)
assert missing(got)
test_pass "empty selected sample"
capture cipolate y x, gen(got)
assert _rc == 110
test_pass "existing destination rejected"
capture by g, sort: cipolate y x, gen(bad) by(g)
assert _rc == 190
capture confirm variable bad
assert _rc == 111
test_pass "by option conflict leaves no output"
capture cipolate y x, gen(bad) threads(-1)
assert _rc == 198
test_pass "invalid thread count"
clear
set obs 0
gen double y=.
gen double x=.
cipolate y x, gen(z)
assert _N==0
test_pass "zero-observation data"
}
print_summary "cipolate"
global CTOOLS_COMPONENT_COMPLETE "cipolate"
