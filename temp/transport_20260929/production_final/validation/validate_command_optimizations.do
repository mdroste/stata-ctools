* Focused reference coverage for optimized kernels and wrapper fast paths.
do validation/validate_setup.do
quietly {
clear
set seed 43211
set obs 70001
gen long id=_n
gen double x=floor(runiform()*100)
gen double y=rnormal()
replace y=. if runiform()<.6
gen double g=floor(runiform()*100)
replace g=.a if mod(_n,97)==0
replace g=.z if mod(_n,101)==0
foreach groups in "g" "x" "y" "x y" {
    sort id
    ipolate y x, gen(ref) by(`groups') epolate
    sort id
    foreach threads in 1 4 12 {
        cipolate y x, gen(got) by(`groups') epolate threads(`threads')
        assert_var_equal got ref $DEFAULT_SIGFIGS "radix/alias by(`groups') threads(`threads')"
        assert id==_n
        drop got
    }
    drop ref
}
clear
set obs 40001
gen str2045 s=cond(mod(_n,3)==0, 150*"a::東京: "+" end  ", "  :aa::|--z;::end  ")
foreach trim in "" "notrim" {
    foreach limit in 1 5 200 {
        split s if mod(_n,7), gen(ref) parse("::" ":" "|" "--" ";" "東京") `trim' limit(`limit')
        local k=r(k_new)
        foreach threads in 1 4 12 {
            csplit s if mod(_n,7), gen(got) parse("::" ":" "|" "--" ";" "東京") `trim' limit(`limit') threads(`threads')
            assert r(k_new)==`k'
            forvalues j=1/`k' {
                assert got`j'==ref`j'
            }
            drop got*
            test_pass "packed tokens `trim' limit(`limit') threads(`threads')"
        }
        drop ref*
    }
}
* One master row expands across many output tiles; every row must retain ties.
clear
set obs 50001
gen long uid=_n
gen double key=mod(_n,31)
gen str16 text="using"+string(_n)
tempfile u m expected
save `u'
clear
set obs 1
gen double key=15
gen str12 master="one master"
save `m'
rangejoin key -20 20 using `u'
save `expected'
foreach threads in 1 4 12 {
    use `m', clear
    crangejoin key -20 20 using `u', threads(`threads')
    cf _all using `expected'
    test_pass "one-master tiled join threads(`threads')"
}
* Repeated restored master datasets retain metadata and stable using ties.
clear
set obs 40001
gen long mid=_n
gen double g=mod(_n,7)
replace g=.a if mod(_n,17)==0
replace g=.z if mod(_n,23)==0
gen double key=mod(_n,101)
label variable mid "Original master row"
format key %12.3f
save `u', replace
replace key=key+.5
save `m', replace
rangejoin key -1 1 using `u', by(g)
save `expected', replace
foreach threads in 1 4 12 {
    use `m', clear
    crangejoin key -1 1 using `u', by(g) threads(`threads')
    cf _all using `expected'
    assert "`: variable label mid'"=="Original master row"
    assert "`: format key'"=="%12.3f"
    test_pass "indexed extended-missing groups threads(`threads')"
}
}
print_summary "command_optimizations"
global CTOOLS_COMPONENT_COMPLETE "command_optimizations"
