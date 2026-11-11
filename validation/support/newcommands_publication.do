* Requires a separate serial plugin compiled with publication_fault.h included
* after the headers of csplit, cinterpolate and crangejoin. Native OpenMP tests
* use thread-safe mocks instead of reading these Stata injection scalars.
local enabled : environment CTOOLS_TEST_PUBLICATION_PLUGIN
assert "`enabled'" == "1"
do "validation/support/validate_setup.do"
clear
set obs 40
gen str8 source = "a:b:c"
gen double x = _n
gen double y = cond(mod(_n,2),.,2*_n)
label data "publication rollback fixture"
tempfile original
save `original'
scalar __ctools_test_macro_after = 0
foreach checkpoint in 1 60 120 {
    scalar __ctools_test_store_after = `checkpoint'
    capture noisily csplit source, generate(piece) parse(":")
    assert _rc == 459
    cf _all using `original', all
    confirm new variable piece1 piece2 piece3
    test_pass "split store failure `checkpoint' publishes no names"
    scalar __ctools_test_store_after = 0
    quietly csplit source, generate(piece) parse(":")
    assert piece1=="a" & piece2=="b" & piece3=="c"
    drop piece1 piece2 piece3
    test_pass "split recovers after store failure `checkpoint'"
}
scalar __ctools_test_macro_after = 1
capture noisily csplit source, generate(piece) parse(":")
assert _rc == 459
cf _all using `original', all
confirm new variable piece1 piece2 piece3
test_pass "split width publication failure publishes no names"
scalar __ctools_test_macro_after = 1
capture noisily cinterpolate y x, generate(interpolated) epolate
assert _rc == 459
cf _all using `original', all
confirm new variable interpolated
test_pass "interpolation final macro failure publishes no name"
scalar __ctools_test_macro_after = 0
quietly cinterpolate y x, generate(interpolated) epolate
assert interpolated == 2*x
drop interpolated
test_pass "interpolation recovers after publication failure"

clear
set obs 4
gen key=_n
gen value=10*_n
gen str5 tag="u"+string(_n)
tempfile usingdata master
save `usingdata'
clear
set obs 10
gen long id=_n
label data "original join master"
save `master'
foreach checkpoint in 1 60 120 {
    scalar __ctools_test_store_after = `checkpoint'
    capture noisily crangejoin key 1 3 using `usingdata'
    assert _rc == 459
    cf _all using `master', all
    assert _N==10 & id==_n
    confirm new variable key value tag
    local label : data label
    assert "`label'"=="original join master"
    test_pass "join store failure `checkpoint' restores master"
    scalar __ctools_test_store_after = 0
    quietly crangejoin key 1 3 using `usingdata'
    assert _N==30 & id==ceil(_n/3) & value==10*key
    test_pass "join recovers after store failure `checkpoint'"
    use `master', clear
}
scalar __ctools_test_macro_after = 1
capture noisily crangejoin key 1 3 using `usingdata'
assert _rc == 459
cf _all using `master', all
confirm new variable key value tag
test_pass "join count publication failure restores master"
print_summary "new command publication rollback"
