/* Already ordered input without a Stata sort marker must retain every payload
   and skip stores. A final-row inversion must still match native stable sort.
   str2045 payloads force the column fallback; other cases use row staging. */
capture do "validation/support/validate_setup.do"
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

capture program drop _csort_strl_case
program define _csort_strl_case
    syntax, STATE(string) MODE(string) KEY(string) API(string)
    clear
    local n = cond("`state'" == "empty", 0, cond("`state'" == "single", 1, 4))
    quietly set obs `n'
    generate long key = _n
    generate strL text = "payload:" + string(_n)
    generate long row = _n
    if "`state'" == "inverted" quietly replace key = -1 in 4
    if "`state'" == "marked" quietly sort key, stable
    local marked : sortedby
    quietly datasignature
    local before `"`r(datasignature)'"'
    local options ""
    if "`mode'" == "standard" local options "nostream"
    if "`mode'" == "stream" local options "stream(1)"
    if "`api'" == "ado" {
        capture csort `key', `options' nosortedby
    }
    else {
        * Exercise C validation independently of the ado preflight.
        _ctools_load
        capture program ctools_plugin, plugin using("`__ctools_plugin'")
        if _rc != 0 & _rc != 110 exit 601
        local index = cond("`key'" == "key", 1, 2)
        capture plugin call ctools_plugin key text row, "csort `index' `options'"
    }
    local rc = _rc
    assert `rc' == 109
    assert row == _n & text == "payload:" + string(_n)
    quietly datasignature
    assert `"`r(datasignature)'"' == `"`before'"'
    local after : sortedby
    assert `"`after'"' == `"`marked'"'
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
    foreach state in empty single ordered inverted marked {
        foreach mode in auto standard stream {
            foreach key in key text {
                foreach api in ado plugin {
                    capture noisily _csort_strl_case, state(`state') mode(`mode') key(`key') api(`api')
                    if _rc test_fail "strL: `state' `mode' `key' `api'" "r(`=_rc')"
                    else test_pass "strL: `state' `mode' `key' `api'"
                }
            }
        }
    }
}
print_summary "csort_identity"
global CTOOLS_COMPONENT_COMPLETE "csort_identity"
