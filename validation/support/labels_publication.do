* Separate serial plugin with publication_fault.h in cencode/cdecode only.
local enabled : environment CTOOLS_TEST_PUBLICATION_PLUGIN
assert "`enabled'" == "1"
do "validation/support/validate_setup.do"
clear
set obs 40
gen str8 source=cond(mod(_n,2),"alpha","zz")
gen long coded=cond(mod(_n,2),-7,1000)
label define existing -7 "alpha" 3 "alpha" 9 "é" 1000 "b"
label values coded existing
label variable source "original text"
tempfile original
save `original'
mata: original_codes=J(0,1,.); original_text=J(0,1,""); got_codes=J(0,1,.); got_text=J(0,1,"")
mata: st_vlload("existing", original_codes, original_text)
scalar __ctools_test_macro_after=0
scalar __ctools_test_store_after=0
forvalues checkpoint=1/14 {
    scalar __ctools_test_scalar_after=`checkpoint'
    capture noisily cencode source, generate(encoded) label(existing)
    assert _rc==459
    cf _all using `original', all
    confirm new variable encoded
    mata: st_vlload("existing", got_codes, got_text)
    mata: assert(original_codes==got_codes & original_text==got_text)
    test_pass "encoder scalar failure `checkpoint' restores labels and publishes no output"
}
scalar __ctools_test_scalar_after=0
cencode source, generate(encoded) label(existing)
assert encoded==cond(mod(_n,2),-7,1001)
test_pass "encoder recovers and extends labels after failures"
foreach mode in generate replace {
    forvalues checkpoint=1/4 {
        use `original', clear
        scalar __ctools_test_scalar_after=`checkpoint'
        if "`mode'"=="generate" {
            cdecode coded, generate(decoded)
            assert r(N_plugin)==0
            assert decoded==cond(mod(_n,2),"alpha","b")
            assert coded==cond(mod(_n,2),-7,1000)
        }
        else {
            cdecode coded, replace
            assert r(N_plugin)==0
            assert coded==cond(mod(_n,2),"alpha","b")
        }
        test_pass "decoder `mode' scalar failure `checkpoint' falls back with complete correct output"
    }
}
scalar __ctools_test_scalar_after=0
use `original', clear
cdecode coded, generate(decoded)
assert r(N_plugin)==1
assert decoded==cond(mod(_n,2),"alpha","b")
test_pass "decoder recovers to plugin after publication failures"
print_summary "label publication rollback and fallback"
