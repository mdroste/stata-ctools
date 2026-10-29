/* Already ordered input without a Stata sort marker must retain every payload
   and skip stores. A final-row inversion must still match native stable sort.
   str2045 payloads force the column fallback; other cases use row staging. */
capture do "validation/validate_setup.do"
if _rc != 0 do "validate_setup.do"

capture program drop _csort_identity_case
program define _csort_identity_case
    syntax, LAYOUT(string) ALGorithm(string) THReads(integer)
    clear
    quietly set obs 4096
    generate double key = floor((_n - 1) / 16)
    quietly replace key = -0 in 1
    quietly replace key = . in 4094
    quietly replace key = .a in 4095
    quietly replace key = .z in 4096
    generate long row = _n
    generate double payload = sin(_n) * _n
    local keys "key"
    if "`layout'" == "mixed" {
        generate str8 text = string(floor((_n - 1) / 32), "%08.0f")
        local keys "text key"
    }
    if "`layout'" == "fallback" {
        generate str2045 text = "row:" + string(_n)
    }
    local marked : sortedby
    assert `"`marked'"' == ""
    quietly datasignature
    local before `"`r(datasignature)'"'
    quietly csort `keys', algorithm(`algorithm') threads(`threads') nostream nosortedby
    assert r(time_store) == 0
    assert r(stream) == 0
    assert row == _n
    quietly datasignature
    assert `"`r(datasignature)'"' == `"`before'"'

    /* Only the last adjacent pair is inverted. A stable native oracle checks
       both the full key order and row alignment; nosortedby prevents Stata's
       metadata sort from concealing an incorrect plugin order. */
    quietly replace key = -1 in 4096
    tempfile expected
    preserve
    quietly sort `keys', stable
    quietly save `expected'
    restore
    quietly csort `keys', algorithm(`algorithm') threads(`threads') nostream nosortedby
    assert r(time_store) > 0
    quietly cf _all using `expected'
    local marked : sortedby
    assert `"`marked'"' == ""
end

quietly {
    foreach threads in 1 8 {
        foreach layout in numeric mixed fallback {
            foreach algorithm in auto lsd msd timsort sample counting merge ips4o {
                capture noisily _csort_identity_case, layout(`layout') ///
                    algorithm(`algorithm') threads(`threads')
                if _rc {
                    test_fail "identity/inversion: `layout' `algorithm' threads(`threads')" "r(`=_rc')"
                }
                else test_pass "identity/inversion: `layout' `algorithm' threads(`threads')"
            }
        }
    }
}
print_summary "csort_identity"
global CTOOLS_COMPONENT_COMPLETE "csort_identity"
