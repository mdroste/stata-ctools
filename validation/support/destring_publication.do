* Separate serial test plugin: publication_fault.h in cdestring_impl.c only.
local enabled : environment CTOOLS_TEST_PUBLICATION_PLUGIN
assert "`enabled'" == "1"
do "validation/support/validate_setup.do"
clear
set obs 40
gen str8 source = string(_n)
gen str8 second = string(2*_n)
label var source "original source"
label data "destring publication fixture"
tempfile original
save `original'
scalar __ctools_test_macro_after = 0
scalar __ctools_test_store_after = 0
foreach mode in generate replace {
    forvalues checkpoint = 1/9 {
        scalar __ctools_test_scalar_after = `checkpoint'
        if "`mode'"=="generate" capture noisily cdestring source second, generate(a b)
        else capture noisily cdestring source second, replace
        assert _rc == 459
        cf _all using `original', all
        confirm new variable a b
        local label : variable label source
        assert "`label'"=="original source"
        test_pass "destring `mode' scalar failure `checkpoint' leaves originals"
    }
    scalar __ctools_test_scalar_after = 0
    scalar __ctools_test_macro_after = 1
    if "`mode'"=="generate" capture noisily cdestring source second, generate(a b)
    else capture noisily cdestring source second, replace
    assert _rc == 459
    cf _all using `original', all
    confirm new variable a b
    test_pass "destring `mode' macro failure leaves originals"
}
scalar __ctools_test_macro_after = 0
forvalues checkpoint = 1/7 {
    scalar __ctools_test_scalar_after = `checkpoint'
    capture noisily cdestring source if 0, generate(a) verbose
    assert _rc == 459
    cf _all using `original', all
    confirm new variable a
    test_pass "destring empty scalar failure `checkpoint' publishes no name"
}
scalar __ctools_test_scalar_after = 0
cdestring source second, generate(a b)
assert a==_n & b==2*_n
test_pass "destring recovers after publication failures"
print_summary "destring publication rollback"
