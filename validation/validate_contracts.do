/* Boundary regressions: invalid wire arguments must not mutate observations. */
do validation/validate_setup.do
clear
set obs 3
gen double result = 17
_ctools_load
program p1_contracts_plugin, plugin using("`__ctools_plugin'")
foreach args in "1 0 count=18446744073709551616" "1 0 percent=nan" "2147483648 0 count=1" "1 -1 count=1" "1 0 count=1 count=2" "1 0 seedhi=1 count=1" {
    capture plugin call p1_contracts_plugin result, "csample `args'"
    local rc = _rc
    if `rc' == 198 test_pass "rejected malformed sampling arguments"
    else test_fail "malformed sampling arguments" "expected 198, got `rc'"
    assert result == 17
}
foreach cmd in "csplit write" "crangejoin prepare 1" "crangejoin write" "cmerge execute" "cio column 1" {
    capture plugin call p1_contracts_plugin result, "`cmd'"
    local rc = _rc
    if `rc' == 198 test_pass "rejected stale command phase"
    else test_fail "stale command phase" "expected 198, got `rc'"
    assert result == 17
}
* Failure recovery and repeated cleanup remain safe in the real dispatcher.
foreach cmd in "cimport clear" "cio clear" "cmerge clear" "csplit clear" "crangejoin clear" {
    plugin call p1_contracts_plugin, "`cmd'"
    plugin call p1_contracts_plugin, "`cmd'"
    test_pass "repeated cleanup"
}
csample 100
assert _N == 3 & result == 17
test_pass "successful command after malformed calls"
program drop p1_contracts_plugin
print_summary "contracts"
global CTOOLS_COMPONENT_COMPLETE "contracts"
