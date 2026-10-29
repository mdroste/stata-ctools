* Regression coverage for cross-group scheduling, tied windows, and sparse-key
* sliding scans. Reference comparisons retain the suite's seven-digit policy.
capture do "validation/support/validate_setup.do"
if _rc do "validate_setup.do"
quietly {
    foreach size in 63 64 499 500 999 1000 9999 10000 {
        foreach exclude in 0 1 {
            clear
            set obs `=8*`size''
            gen long id = _n
            gen long g = floor((_n-1)/`size')
            gen str8 gs = "g" + string(g)
            replace gs = "" if g==7
            replace g = .a if g==7
            gen double key = mod(_n-1,`size')
            if mod(`size',2) == 0 replace key = floor(key/4)
            gen double x = 100 + mod(id,29)/4
            replace x = . if mod(id,17)==0
            replace x = .b if mod(id,37)==0
            replace key = . if mod(id,131)==0
            replace key = .a if mod(id,137)==0
            * Original order intentionally differs from the sorted group order.
            gen double shuffle = mod(id*7919,8*`size'+1)
            sort shuffle id
            drop shuffle
            local option
            if `exclude' local option excludeself
            foreach stat in count sum mean sd min max median skewness kurtosis {
                rangestat (`stat') ref=x, interval(key -2 2) by(g gs) `option'
                crangestat (`stat') got=x, interval(key -2 2) by(g gs) `option' threads(8)
                capture assert missing(got)==missing(ref)
                local rc = _rc
                capture assert got==ref | abs(got-ref)/max(abs(got),abs(ref)) < $DEFAULT_TOL if !missing(ref)
                if `rc' | _rc test_fail "group windows `size' `exclude' `stat'" "reference values differ"
                else test_pass "group windows `size' `exclude' `stat'"
                drop got ref
            }
        }
    }

    * Unbounded, reversed, empty, and offset windows with many equal keys.
    foreach bounds in ". ." ". 0" "0 ." "5 2" "100 101" "2 4" {
        foreach exclude in 0 1 {
            clear
            set obs 1600
            gen long id = _n
            gen long g = floor((_n-1)/200)
            gen double key = floor(mod(_n-1,200)/20)
            gen double x = mod(_n,31)
            replace x = . if g==2
            replace x = .a if mod(id,19)==0
            replace key = . if mod(id,23)==0
            local option
            if `exclude' local option excludeself
            foreach stat in count mean min max median {
                rangestat (`stat') ref=x, interval(key `bounds') by(g) `option'
                crangestat (`stat') got=x, interval(key `bounds') by(g) `option' threads(4)
                capture assert missing(got)==missing(ref)
                local rc = _rc
                capture assert got==ref | abs(got-ref)/max(abs(got),abs(ref)) < $DEFAULT_TOL if !missing(ref)
                if `rc' | _rc test_fail "tied bounds `bounds' `exclude' `stat'" "reference values differ"
                else test_pass "tied bounds `bounds' `exclude' `stat'"
                drop got ref
            }
        }
    }
    * Both sides of the run-count dispatch boundary in one large group.
    foreach mode in 1 2 3 4 {
        clear
        set obs 2048
        gen double key = _n
        if `mode'==2 replace key = floor(_n/2)
        if `mode'==3 replace key = floor((_n-1)/2)
        if `mode'==4 replace key = floor(_n/4)
        gen double x = 100 + mod(_n,41)/4
        foreach stat in mean sd min max median {
            rangestat (`stat') ref=x, interval(key -2 2)
            crangestat (`stat') got=x, interval(key -2 2) threads(8)
            capture assert got==ref | abs(got-ref)/max(abs(got),abs(ref)) < $DEFAULT_TOL
            if _rc test_fail "run dispatch `mode' `stat'" "reference values differ"
            else test_pass "run dispatch `mode' `stat'"
            drop got ref
        }
    }
    * The additional cross-group cases must not send wide min/max requests
    * through the path whose per-group sparse-table buffer has 16 entries.
    foreach width in 16 17 {
        clear
        set obs 8000
        gen long g = floor((_n-1)/1000)
        gen double key = mod(_n-1,1000)
        local stats
        forvalues j=1/`width' {
            gen double x`j' = 100 + mod(_n*`j',47)
            local stats `stats' (min) got`j'=x`j'
        }
        rangestat (min) x1-x`width', interval(key -50 50) by(g)
        crangestat `stats', interval(key -50 50) by(g) threads(8)
        forvalues j=1/`width' {
            capture assert got`j'==x`j'_min
            if _rc test_fail "wide group windows `width' `j'" "reference values differ"
            else test_pass "wide group windows `width' `j'"
        }
    }
}
print_summary "group_windows"

global CTOOLS_COMPONENT_COMPLETE "group_windows"
