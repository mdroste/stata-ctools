* Combined target/by loading: index spaces, selection and exact grouping.
capture do "validation/validate_setup.do"
if _rc do "validate_setup.do"
quietly {
    foreach n in 1 2 3 1000 20000 {
        clear
        set seed 981723
        set obs `n'
        gen long id = _n
        gen double x = rnormal()*50
        gen double y = rnormal()*2
        gen double z = mod(_n,13)
        replace x = .a if mod(_n,17)==0
        replace y = .z if mod(_n,31)==0
        gen double g = mod(_n,7) + 1e-9*mod(_n,3)
        replace g = .b if mod(_n,23)==0
        gen str20 s = cond(mod(_n,3)==0,"",cond(mod(_n,2),"prefix_a","prefix_b"))
        gen double shuffle=runiform()
        sort shuffle
        drop shuffle
        foreach keys in none g s gs {
            local byopt
            if "`keys'" == "g" local byopt by(g)
            if "`keys'" == "s" local byopt by(s)
            if "`keys'" == "gs" local byopt by(g s)
            foreach threads in 1 8 {
                foreach mode in winsor trim {
                    local trimopt
                    if "`mode'" == "trim" local trimopt trim
                    foreach count in 1 3 {
                        local vars x
                        if `count' == 3 local vars x y z
                        local sample
                        if `n' > 3 local sample if mod(id,5)<2 in 2/`n'
                        noi di "CASE n=`n' keys=`keys' mode=`mode' vars=`vars' threads=`threads'"
                        cwinsor `vars' `sample', `byopt' `trimopt' cuts(5 95) suffix(_got) threads(`threads')
                        gstats winsor `vars' `sample', `byopt' `trimopt' cuts(5 95) suffix(_ref)
                        foreach v of local vars {
                            noi list id x x_got x_ref if missing(x_got)!=missing(x_ref), noobs
                            assert missing(`v'_got) == missing(`v'_ref)
                            assert `v'_got == `v'_ref | abs(`v'_got-`v'_ref)/max(abs(`v'_got),abs(`v'_ref)) < $DEFAULT_TOL if !missing(`v'_got)
                        }
                        drop *_got *_ref
                        test_pass "winsor combined transport matches gstats"
                    }
                }
            }
        }
    }
}
print_summary "winsor_transport"
if $TESTS_FAILED exit 9
di "WINSOR_TRANSPORT_COMPLETE"

global CTOOLS_COMPONENT_COMPLETE "winsor_transport"
