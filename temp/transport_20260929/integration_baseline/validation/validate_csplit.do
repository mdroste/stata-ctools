* Reference comparisons for C string tokenization and result metadata.
do "validation/validate_setup.do"
quietly {
clear
input str100 s
"a b c"
"  d  e  "
"one"
""
" "
",a,,b,"
"a::b:c:::d"
" a , b , c "
"école 東京 café"
"x|y||z|"
"a,,"
end
local test = 0
foreach opts in "" `"parse(",")"' `"parse("::" ":")"' `"parse(":" "::")"' ///
    `"parse("|") notrim"' `"parse(",") notrim"' `"parse(" ")"' ///
    `"parse(" " ",")"' `"parse(",") limit(2)"' "limit(1)" {
    foreach sample in "" "if mod(_n,2)==0" "in 2/8" {
        split s `sample', gen(ref) `opts'
        local k = r(k_new)
        csplit s `sample', gen(got) `opts'
        assert r(k_new)==`k'
        forvalues j=1/`k' {
            assert_strvar_equal got`j' ref`j' "split `++test' field `j'"
            assert "`: type got`j''" == "`: type ref`j''"
        }
        drop ref* got*
    }
}
* Number conversion options retain native destring behavior.
clear
input str30 s
"1,2,3"
"4,.,6"
"7,x,9"
"10%,20%,30%"
end
foreach opts in `"parse(",") destring"' `"parse(",") destring force float"' ///
    `"parse(",") destring force percent"' `"parse(",") destring ignore("%" "x")"' {
    split s, gen(ref) `opts'
    csplit s, gen(got) `opts'
    forvalues j=1/3 {
        assert got`j'==ref`j'
        assert "`: type got`j''" == "`: type ref`j''"
    }
    test_pass "destring `opts'"
    drop ref* got*
}
* Repeated calls after a naming conflict cannot reuse stale cached strings.
gen got2=123
capture csplit s, gen(got) parse(",")
assert _rc==110
capture confirm variable got1
assert _rc==111
assert got2==123
test_pass "name conflict is atomic"
drop got2
csplit s, gen(got) parse(",")
assert got1=="1" in 1
drop got*
test_pass "recovery after error"
* Long strings use the documented native path, without truncation.
recast strL s
replace s = 2500*"x" + ",end" in 1
split s, gen(ref) parse(",")
csplit s, gen(got) parse(",") verbose
assert ref1==got1 & ref2==got2 & ref3==got3
test_pass "strL fallback preserves long fields"
clear
set obs 25000
gen str80 s=" first :: second :: 東京 ::last "
split s, gen(ref) parse("::")
csplit s, gen(got) parse("::") threads(4) verbose
forvalues j=1/4 {
    assert_strvar_equal got`j' ref`j' "parallel split field `j'"
}
capture csplit s if 0, gen(none)
assert _rc==2000
test_pass "empty selected sample"
capture csplit s, gen(none) force
assert _rc==198
test_pass "destring options require destring"
}
print_summary "csplit"
global CTOOLS_COMPONENT_COMPLETE "csplit"
