* Regressions for the eleven findings in the November 24 six-command audit.
args command

if "`command'" == "cdestring" {
    foreach mode in generate replace force {
        clear
        set obs 3
        gen str20 text = char(34) + "12" + char(34)
        replace text = char(34) + "-2.5" + char(34) in 2
        replace text = char(34) + ".a" + char(34) in 3
        local quote = char(34)
        if "`mode'" == "generate" {
            cdestring text, generate(value) ignore(`"`quote'"')
        }
        else {
            local force = cond("`mode'" == "force", "force", "")
            cdestring text, replace `force' ignore(`"`quote'"')
            rename text value
        }
        assert value == 12 in 1
        assert value == -2.5 in 2
        assert value == .a in 3
        test_pass "AUDIT-01 literal quote ignore: `mode'"
    }
    clear
    set obs 2
    gen str20 text = "12"
    foreach suffix in n_converted n_failed time_parse time_load time_convert time_store time_total threads_max openmp_enabled {
        gen str1 _cdestring_`suffix' = "x"
    }
    cdestring text, replace verbose
    assert text == 12
    assert scalar(_cdestring_n_converted) == 2 & scalar(_cdestring_n_failed) == 0
    test_pass "AUDIT-05 destring scalar names cannot be shadowed"
}

if "`command'" == "cencode" {
    foreach mode in generate replace extend {
        clear
        set obs 2
        gen str8 text = cond(_n == 1, "alpha", "beta")
        gen byte _cencode_n_unique = 0
        gen byte _cencode_n_new = 0
        foreach suffix in protocol fallback max_len time_parse time_load time_collect time_sort time_encode time_labels time_total threads_max openmp_enabled {
            gen str1 _cencode_`suffix' = "x"
        }
        if "`mode'" == "extend" {
            label define audit_label 7 "existing", replace
            cencode text, generate(value) label(audit_label) verbose
        }
        else if "`mode'" == "replace" {
            cencode text, replace verbose
            rename text value
        }
        else cencode text, generate(value) verbose
        decode value, generate(roundtrip)
        assert roundtrip == cond(_n == 1, "alpha", "beta")
        test_pass "AUDIT-05 encode retains labels despite scalar-name variables: `mode'"
    }
}

if "`command'" == "cdecode" {
    clear
    set obs 2
    gen byte value = _n
    label define audit_label 1 "alpha" 2 "beta", replace
    label values value audit_label
    gen str1 _cdecode_threads_max = "x"
    cdecode value, generate(text) verbose
    assert text == cond(_n == 1, "alpha", "beta")
    test_pass "AUDIT-05 decode diagnostic scalar cannot be shadowed"
}

if "`command'" == "csample" {
    clear
    set obs 20
    gen id = _n
    foreach suffix in kept dropped ngroups threads time_load time_sort time_sample time_store time_total {
        gen str1 _csample_`suffix' = "x"
    }
    csample 50
    assert _N == 10
    csample 50, verbose
    assert _N == 5
    test_pass "AUDIT-05 sample reporting succeeds with scalar-name variables"
}

if "`command'" == "crangestat" {
    foreach size in 7 8 10 63 64 1000 {
        foreach position in 1 `size' {
            clear
            set obs `size'
            gen key = 0
            gen double x = 1
            replace x = 1e16 in `position'
            crangestat (sum) total=x (skewness) skew=x (kurtosis) kurt=x, interval(key . .) excludeself
            assert total == `size' - 1 in `position'
            assert missing(skew) & missing(kurt) in `position'
            test_pass "AUDIT-02 excluded outlier: group `size', position `position'"
        }
    }
    foreach offset in 1.0000000000000002 .12345678901234567 123456789012345.67 {
        clear
        set obs 3
        gen double key = cond(_n == 1, 0, `offset')
        gen double x = _n
        gen double y = 2 * x
        crangestat (count) count=x (cov) x y, interval(key `offset' `offset')
        assert count == 2 in 1
        assert cov_x == 1 in 1
        test_pass "AUDIT-04 exact scalar interval survives C and Mata transport: `offset'"
    }
    * Scalar addition overflow becomes missing/unbounded, as with a bound
    * generated in Stata. Exercise small groups, cross-group and run paths.
    foreach size in 3 1001 {
        foreach groups in 1 8 {
            foreach offset in 1e307 -1e307 {
                foreach side in low high {
                    clear
                    set obs `=`size'*`groups''
                    gen g = floor((_n - 1) / `size')
                    gen double key = cond(mod(_n - 1, 3) == 0, -8e307, cond(mod(_n - 1, 3) == 1, 0, 8e307))
                    gen x = 1
                    gen double bound = key + `offset'
                    local lo = cond("`side'" == "low", "`offset'", ".")
                    local hi = cond("`side'" == "high", "`offset'", ".")
                    local lovar = cond("`side'" == "low", "bound", ".")
                    local hivar = cond("`side'" == "high", "bound", ".")
                    crangestat (count) expected=x, interval(key `lovar' `hivar') by(g) threads(1)
                    crangestat (count) actual=x, interval(key `lo' `hi') by(g) threads(4)
                    assert actual == expected
                    test_pass "AUDIT-07 overflow bounds: size `size', groups `groups', `side' `offset'"
                }
            }
        }
    }
    clear
    set obs 10
    gen key = _n
    gen x = _n
    foreach suffix in time_load time_sort time_groups time_compute time_store time_total threads_max openmp_enabled {
        gen str1 _crangestat_`suffix' = "x"
    }
    crangestat (count) count=x, interval(key . .) verbose
    assert count == 10
    test_pass "AUDIT-05 rangestat diagnostic scalars cannot be shadowed"
}

if "`command'" == "cbsample" {
    foreach strata in 0 1 {
        foreach filter in 0 1 {
            clear
            set obs 6
            gen c = 7 - _n
            gen s = 1
            gen double w = .
            local st = cond(`strata', "strata(s)", "")
            local selection = cond(`filter', "if c > 1", "")
            local expected = 1 + `filter'
            cbsample c[1] `selection', cluster(c) `st' weight(w)
            summarize w, meanonly
            assert r(sum) == `expected'
            test_pass "AUDIT-03 cluster-sorted expression: strata `strata', filter `filter'"
        }
    }
    clear
    set obs 20
    gen id = _n
    foreach suffix in idmax idrandom n_selected total_weight nclusters nstrata threads time_load time_sort time_sample time_store time_total {
        gen str1 _cbsample_`suffix' = "x"
    }
    cbsample 10, cluster(id) idcluster(newid) verbose
    assert _N == 10 & !missing(newid)
    test_pass "AUDIT-05 bootstrap diagnostic and id scalars cannot be shadowed"
    * cbsample keeps no copy of the data, as bsample does not. A failure
    * before the sample is drawn must still leave every row in place and no
    * idcluster() variable behind: here maxvar has no room for it.
    clear
    local oldmax = c(maxvar)
    set maxvar 2048
    set obs 12
    gen c = _n
    forvalues j = 1/2041 {
        gen byte v`j' = 1
    }
    capture cbsample 6, cluster(c) idcluster(newid)
    local rc = _rc
    assert `rc' == 900 & _N == 12 & c(k) == 2042
    assert c == _n
    capture confirm variable newid
    assert _rc == 111
    clear
    set maxvar `oldmax'
    test_pass "AUDIT-09 maxvar failure preserves bootstrap input"
}
