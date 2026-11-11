* Run only against a test plugin whose host poll callback decrements
* __ctools_test_cancel_after and returns Break at zero. Production builds
* have no scalar-controlled cancellation hook.
local enabled : environment CTOOLS_TEST_CANCEL_PLUGIN
assert "`enabled'" == "1"
do "validation/support/validate_setup.do"
set seed 2626
set obs 1000
gen long id = _n
gen g = floor((_n-1)/20)
gen h = mod(floor((_n-1)/20) + (mod(_n,20)==0),50)
gen double x = rnormal()
gen double z = rnormal()
gen double d = .5*z+rnormal()
gen double y = .4*x+.3*d+rnormal()+g/10
gen double count = rpoisson(exp(.2*x+.1*d))
gen double w = 1+mod(_n,3)
label variable y "original outcome"
label data "cancellation rollback fixture"
tempfile original
save `original'
foreach mode in "generate(r_) sample(used)" "clear keepids keepvars(id)" {
    foreach weight in "" "[aw=w]" {
        foreach transform in sym kacz cim rand {
            foreach acceleration in cg none sd aitken hybrid {
                foreach checkpoint in 1 5 {
                    scalar __ctools_test_cancel_after = `checkpoint'
                    capture noisily chdfe y x `weight', absorb(g h) `mode' ///
                        transform(`transform') acceleration(`acceleration') ///
                        maxiterations(.) tolerance(1e-16) threads(4)
                    assert _rc == 1
                    cf _all using `original', all
                    local datalabel : data label
                    local ylabel : variable label y
                    assert "`datalabel'" == "cancellation rollback fixture"
                    assert "`ylabel'" == "original outcome"
                    confirm new variable r_y r_x used
                    test_pass "Break: chdfe `mode' `weight' `transform' `acceleration' at `checkpoint'"
                }
            }
        }
    }
}
foreach command in "creghdfe y x, absorb(g h)" ///
    "civreghdfe y x (d=z), absorb(g h)" ///
    "cppmlhdfe count x, absorb(g h)" {
    foreach checkpoint in 1 5 {
        scalar __ctools_test_cancel_after = `checkpoint'
        capture noisily `command'
        assert _rc == 1
        cf _all using `original', all
        test_pass "Break: `command' at `checkpoint'"
        scalar __ctools_test_cancel_after = -1
        quietly `command'
        assert e(N) > 0
        cf _all using `original', all
        test_pass "Fresh fit after Break: `command'"
    }
}
scalar __ctools_test_cancel_after = -1
quietly chdfe y x, absorb(g h) generate(r_)
confirm variable r_y r_x
test_pass "Fresh chdfe after Break"
print_summary "HDFE cancellation rollback"
