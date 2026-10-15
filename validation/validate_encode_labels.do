* Exact code and byte comparisons for the cencode label data handoff.
capture do "validation/validate_setup.do"
if _rc do "validate_setup.do"

capture program drop check_label_parity
program define check_label_parity
    syntax [, SUBSET(string) THREADS(integer 8)]
    capture drop got expected decoded_got decoded_expected
    capture label drop got expected
    cencode x `subset', gen(got) threads(`threads')
    encode x `subset', gen(expected)
    assert got == expected
    decode got, gen(decoded_got)
    decode expected, gen(decoded_expected)
    assert decoded_got == decoded_expected
    mata: st_vlload("got", gc=., gt=""); st_vlload("expected", ec=., et=""); assert(gc==ec); assert(gt==et)
    test_pass "exact codes, complete label maps and decoded bytes"
end

mata:
void encode_bad_fixture(string scalar path, real scalar kind)
{
    real scalar fh
    string scalar data, nl
    nl = char(10)
    data = "CTOOLS_LABELS_1" + nl + "1 3" + nl + "one"
    if (kind == 1) data = "bad magic" + nl
    if (kind == 2) data = data + "2 4" + nl + "two"
    if (kind == 3) data = data + "2 -1" + nl
    if (kind == 4) data = data + "2 2046" + nl
    if (kind == 5) data = data + "2 1.5" + nl
    if (kind == 6) data = data + "2" + nl
    if (kind == 7) data = data + "x 3" + nl + "two"
    if (kind == 8) data = data + "2 3" + nl + "two" + "extra"
    if (kind == 9) data = data + "2147483648 3" + nl + "two"
    fh = fopen(path, "w")
    fwrite(fh, data)
    fclose(fh)
}
end


quietly {
    clear
    set obs 270
    gen str2045 x = "prefix_" + char(mod(_n-1,255)+1) + "_suffix"
    replace x = char(96) + "local" + char(39) in 256
    replace x = char(36) + "CTOOLS_LABEL_SENTINEL" in 257
    replace x = char(34) + char(96) + char(34) + "compound" + char(34) + char(39) in 258
    replace x = "x" + char(13) + char(10) + "y" in 259
    replace x = " leading and trailing " in 260
    replace x = "é東京🚀" in 261
    mata: st_sstore(262, "x", invtokens(J(1,2045,"z"), ""))
    mata: st_sstore(263, "x", invtokens(J(1,681,"界"), ""))
    replace x = "" in 264
    mata: st_sstore(265, "x", invtokens(J(1,128,"a"), "") + char(10) + invtokens(J(1,128,"b"), ""))
    replace x = "end" + char(10) in 266
    replace x = char(13) + "begin" in 267
    replace x = "/* */ // ; mata: _error(9)" in 268
    replace x = "01 2045" + char(10) + "fake header" in 269
    replace x = "CTOOLS_LABELS_1" in 270
    global CTOOLS_LABEL_SENTINEL "must remain literal"
    foreach threads in 1 8 {
        check_label_parity, threads(`threads')
        check_label_parity, threads(`threads') subset("if mod(_n,3)==0 in 12/269")
    }
    assert "$CTOOLS_LABEL_SENTINEL" == "must remain literal"
    macro drop CTOOLS_LABEL_SENTINEL

    * Label counts at small sizes, former expression boundaries, and the limit.
    foreach categories in 1 2 127 128 129 1023 1024 1025 10000 65536 {
        clear
        set obs `categories'
        gen str20 x = "category" + string(_n,"%05.0f")
        check_label_parity
        assert got == _n
    }

    * Empty and entirely excluded samples create no labels.
    clear
    set obs 100
    gen str1 x = ""
    cencode x, gen(got)
    assert missing(got)
    mata: assert(!st_vlexists("got"))
    test_pass "all empty labels"
    drop got
    replace x = "x"
    cencode x if 0, gen(got)
    assert missing(got)
    mata: assert(!st_vlexists("got"))
    test_pass "empty selection"

    * Bad framing must fail before installing any part of the label map.
    * Ado-loaded Mata functions are private to their ado. Load the same source
    * into the test namespace to exercise the actual reader on malformed data.
    findfile cencode.ado
    local ado_source `"`r(fn)'"'
    program drop cencode
    quietly run `"`ado_source'"'
    tempfile bad
    forvalues kind=1/9 {
        capture erase `"`bad'"'
        mata: encode_bad_fixture(st_local("bad"), strtoreal(st_local("kind")))
        capture mata: _cencode_apply_labels(st_local("bad"), "bad_labels", 2)
        local reader_rc = _rc
        assert `reader_rc' == 610
        mata: assert(!st_vlexists("bad_labels"))
        test_pass "reject malformed label data before mutation"
    }
}
print_summary "encode_labels"
if $TESTS_FAILED exit 9
di "ENCODE_LABELS_COMPLETE"

global CTOOLS_COMPONENT_COMPLETE "encode_labels"
