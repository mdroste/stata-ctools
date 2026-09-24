/* Targeted regression tests for audit findings A18--A29. */
do validation/validate_setup.do

program p2_labels
    clear
    set obs 8
    gen str2045 source = ""
    global AUDIT_WORD changed
    replace source = char(36)+"AUDIT_WORD" in 1
    replace source = char(96)+"local_name"+char(39) in 2
    replace source = char(96)+char(34)+"quoted"+char(34)+char(39) in 3
    replace source = "a"+char(9)+"b"+char(10)+"c"+char(13)+"d" in 4
    replace source = "café 日本語" in 5
    replace source = char(34)+"plain quotes"+char(34) in 6
    mata: st_sstore(7,"source",invtokens(J(1,2045,"x"),""))
    cencode source, gen(code)
    decode code, gen(roundtrip)
    assert source == roundtrip
    macro drop AUDIT_WORD
    clear
    set obs 4
    gen str1 s = char(64+_n)
    label define signed -10 "A" 0 "B" 7 "C"
    capture encode s, gen(expected) label(signed) noextend
    local native_rc = _rc
    capture cencode s, gen(actual) label(signed) noextend
    assert _rc == `native_rc' & `native_rc' == 459
    capture drop expected actual
    drop in 4
    encode s, gen(expected) label(signed) noextend
    cencode s, gen(actual) label(signed) noextend
    assert actual == expected
end

program p2_clusters
    clear
    set seed 45021
    set obs 400
    gen double x = rnormal()
    gen double y = .7*x+rnormal()
    gen long g = ceil(_n/20)
    gen double fraction = g/100
    gen double large = 1e12+g
    gen str20 textgroup = string(g)
    cqreg y x, vce(cluster g)
    tempname original_b original_V
    matrix `original_b' = e(b)
    matrix `original_V' = e(V)
    foreach group in fraction large textgroup {
        cqreg y x, vce(cluster `group')
        assert e(N_clust) == 20
        assert_matrix_equal `original_b' e(b) $DEFAULT_SIGFIGS "cluster relabeling coefficients"
        assert_matrix_equal `original_V' e(V) $DEFAULT_SIGFIGS "cluster relabeling covariance"
    }
end

program p2_center
    clear
    set obs 30
    gen t = _n
    gen y = sin(t)
    gen x = cos(t)
    gen z = t^2
    tsset t
    capture civreghdfe y (x=z), robust bw(3) center
    assert _rc == 198
end

program p2_import
    clear
    set obs 2
    gen str5 name = cond(_n==1,"Alice","Bob")
    gen value = _n
    tempfile base
    local csv "`base' with spaces.csv"
    local xlsx "`base' with spaces.xlsx"
    export delimited using "`csv'", replace
    export excel using "`xlsx'", firstrow(variables) replace
    foreach syntax in positional using compound {
        if "`syntax'" == "positional" cimport delimited "`csv'", clear
        if "`syntax'" == "using" cimport delimited using "`csv'", clear
        if "`syntax'" == "compound" cimport delimited `"`csv'"', clear
        assert _N == 2 & name[1] == "Alice" & value[2] == 2
        if "`syntax'" == "positional" cimport excel "`xlsx'", firstrow clear
        if "`syntax'" == "using" cimport excel using "`xlsx'", firstrow clear
        if "`syntax'" == "compound" cimport excel `"`xlsx'"', firstrow clear
        assert _N == 2 & name[1] == "Alice" & value[2] == 2
    }
    foreach kind in delimited excel {
        capture cimport `kind' "`csv'" extra, clear
        assert _rc == 198
        assert _N == 2 & name[1] == "Alice"
    }
    foreach endian in le be {
        tempfile utf32
        tempname fh
        file open `fh' using "`utf32'", write binary replace
        if "`endian'" == "le" file write `fh' %1bu (255) %1bu (254) %1bu (0) %1bu (0)
        else file write `fh' %1bu (0) %1bu (0) %1bu (254) %1bu (255)
        foreach byte in 97 44 98 10 49 44 50 10 {
            if "`endian'" == "le" file write `fh' %1bu (`byte') %1bu (0) %1bu (0) %1bu (0)
            else file write `fh' %1bu (0) %1bu (0) %1bu (0) %1bu (`byte')
        }
        file close `fh'
        capture cimport delimited using "`utf32'", clear
        assert _rc != 0
        assert _N == 2 & name[1] == "Alice" & value[2] == 2
    }
    capture cimport delimited "`csv'", clear encoding(utf-32)
    assert _rc == 198
    assert _N == 2 & name[1] == "Alice"
    erase "`csv'"
    erase "`xlsx'"
end

program p2_sort
    clear
    set obs 600
    gen id = _n
    gen key = mod(600-_n,17)
    gen double nonkey = _n^2
    gen str20 text = "row"+string(_n)
    foreach qualifier in "if id<=2" "in 1/2" {
        capture csort key `qualifier'
        assert _rc == 101
        assert id == _n & nonkey == _n^2
    }
    tempfile source expected
    save `source'
    sort key, stable
    gen rank = _n
    save `expected'
    foreach algorithm in auto lsd msd merge timsort sample counting ips4o {
        use `source', clear
        csort key, algorithm(`algorithm') nosortedby
        gen rank = _n
        cf _all using `expected'
    }
    use `source', clear
    csort key, stream(1) nosortedby
    gen rank = _n
    cf _all using `expected'
end

program p2_sample_by
    foreach pct in 30 100 {
        clear
        set obs 100
        gen long id = _n
        gen byte g = ceil(_n/10)
        tempfile data expected
        save `data'
        sample `pct', by(g)
        collapse (count) n=id, by(g)
        save `expected'
        use `data', clear
        csample `pct', by(g)
        collapse (count) n=id, by(g)
        cf _all using `expected'
    }
end

program p2_io
    clear
    set obs 3
    gen strL text = "long string"
    gen id = _n
    cencode text, gen(encoded)
    decode encoded, gen(decoded)
    assert decoded == text
    drop encoded decoded
    label drop encoded
    capture csort id, stream(1) nosortedby
    assert _rc != 0
    assert text == "long string" & id == _n
    capture csort id, nosortedby
    assert _rc != 0
    assert text == "long string" & id == _n
    mata: st_sstore(., "text", J(3,1,invtokens(J(1,3000,"x"),"")))
    capture cencode text, gen(encoded)
    assert _rc != 0
    assert strlen(text) == 3000 & id == _n
end

program p2_identity
    ctools, version
    assert "`r(version)'" != ""
    assert "`r(version)'" == "`r(plugin_version)'"
    assert "`r(plugin_revision)'" != ""
end

foreach case in labels clusters center import sort sample_by io identity {
    capture noisily p2_`case'
    if _rc test_fail "A18-A29: `case'" "r(`=_rc')"
    else test_pass "A18-A29: `case'"
}
print_summary "P2 audit regressions"

* Reached only after the complete component script.
global CTOOLS_COMPONENT_COMPLETE "audit_p2"
