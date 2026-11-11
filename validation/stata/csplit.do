* csplit: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "csplit"

if "`__case'" == "csplit" {
* Reference comparisons for C string tokenization and result metadata.
do "validation/support/validate_setup.do"
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
* DAT-07: parse("") is the default space parsing, as in split
clear
input str20 name
"a b c"
"d  e"
end
split name, parse("") generate(ref)
local k = r(k_new)
capture csplit name, parse("") generate(got)
if _rc test_fail "DAT-07 empty parse()" "r(`=_rc')"
else {
    assert r(k_new) == `k'
    forvalues j = 1/`k' {
        assert_strvar_equal got`j' ref`j' "DAT-07 empty parse() field `j'"
    }
}
* DATA-10: destring_cmd names the final variables
clear
input str20 s
"1,2,3"
"4,5,6"
end
split s, parse(",") destring generate(r)
csplit s, parse(",") destring generate(c)
mata: st_local("same", strofreal(subinstr(st_global("r1[destring_cmd]"), "r1 r2 r3", "c1 c2 c3") == st_global("c1[destring_cmd]")))
if `same' test_pass "DATA-10 destring_cmd characteristic"
else test_fail "DATA-10 destring_cmd characteristic" "differs from split"
* _ct_output runs a command and returns its output (r(text), lines joined
* by |) and return code (r(rc)). _ct_same_output runs a native and a ctools
* command on copies of the data in memory and passes when the output, return
* code, variable names and types, and (unless "nodata") values agree.
capture program drop _ct_output
program define _ct_output, rclass
    tempfile f
    tempname fh
    capture log close _ctout
    quietly log using `"`f'"', text replace name(_ctout)
    capture noisily `0'
    local rc = _rc
    quietly log close _ctout
    local text
    file open `fh' using `"`f'"', read text
    file read `fh' line
    while !r(eof) {
        if `"`macval(line)'"' != "" & substr(`"`macval(line)'"', 1, 5) != "r; t=" {
            local text `"`macval(text)'`macval(line)'|"'
        }
        file read `fh' line
    }
    file close `fh'
    return local text `"`macval(text)'"'
    return scalar rc = `rc'
end
capture program drop _ct_same_output
program define _ct_same_output
    args name native ctools nodata
    tempfile want
    preserve
    _ct_output `native'
    local want_text `"`r(text)'"'
    local want_rc = r(rc)
    local want_vars
    local want_types
    forvalues j = 1/`c(k)' {
        mata: st_local("v", st_varname(`j'))
        local want_vars `want_vars' `v'
        local want_types `want_types' `: type `v''
    }
    capture quietly save `"`want'"', emptyok
    restore, preserve
    _ct_output `ctools'
    local got_text `"`r(text)'"'
    local got_rc = r(rc)
    local got_vars
    local got_types
    forvalues j = 1/`c(k)' {
        mata: st_local("v", st_varname(`j'))
        local got_vars `got_vars' `v'
        local got_types `got_types' `: type `v''
    }
    local same_data = 1
    if "`nodata'" == "" & "`got_vars'" != "" {
        capture cf _all using `"`want'"'
        local same_data = !_rc
    }
    restore
    if "`nodata'" != "" {
        local want_vars
        local got_vars
        local want_types
        local got_types
    }
    if `"`macval(want_text)'"' == `"`macval(got_text)'"' & `want_rc' == `got_rc' & ///
        "`want_vars'" == "`got_vars'" & "`want_types'" == "`got_types'" & `same_data' {
        test_pass "`name'"
    }
    else {
        noisily di as error `"`name': ctools [`macval(got_text)'] rc `got_rc' vars `got_vars' (`got_types')"'
        noisily di as error `"`name': native [`macval(want_text)'] rc `want_rc' vars `want_vars' (`want_types') same data `same_data'"'
        test_fail "`name'" "output, return code, variables or data differ from the native command"
    }
end
* SPLIT-01, SPLIT-02: split's report (variables created or born as string,
* the ds listing, destring's messages) and r(varlist)
clear
input str20 s str10 one
"1,2,a" "a b"
"3,4,b" "c d"
"5,,c" "e"
end
local i = 0
foreach pair in `"s, parse(,) destring"' `"s, parse(,) gen(t)"' `"s, parse(,) gen(u) destring force"' ///
    `"s, parse(,) gen(f) destring float force"' `"one, gen(q) limit(1)"' `"one, gen(q) destring"' ///
    `"s, parse(,) gen(v) destring ignore("a" "b")"' `"s, parse(,) gen(w) destring ignore("a,b")"' ///
    `"s, parse(,) gen(x) destring percent"' `"s, parse(,) gen(y) force"' `"s, parse(,) gen(y) float"' ///
    `"s, parse(,) gen(y) ignore("a")"' `"s, parse(,) gen(y) percent"' ///
    `"s, parse(,) gen(z) destring ignore(a,b)"' `"s, parse(,) gen(z) destring ignore("a", asbytes illegal)"' {
    local ++i
    _ct_same_output "SPLIT-01/02 output `i'" `"split `pair'"' `"csplit `pair'"'
}
split one, gen(ref) limit(1)
local want `"`r(varlist)'"'
local want_k = r(k_new)
csplit one, gen(got) limit(1)
if `"`r(varlist)'"' == subinstr(`"`want'"', "ref", "got", .) & r(k_new) == `want_k' {
    test_pass "SPLIT-02 r(varlist) as split"
}
else test_fail "SPLIT-02 r(varlist) as split" `"[`r(varlist)'] vs [`want']"'
}
print_summary "csplit"
global CTOOLS_COMPONENT_COMPLETE "csplit"
    global CTOOLS_CASE_COMPLETE "csplit"
    exit
}

di as error "Unknown csplit case: `__case'"
exit 198
