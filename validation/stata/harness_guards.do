* Negative controls for shared comparators: deliberately incorrect candidates
* must fail, while identical results and matching reference errors must pass.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "comparison_guards"

if "`__case'" == "comparison_guards" {
do validation/support/validate_setup.do
do validation/support/coverage_helpers.do

program audit_check
    syntax, COMMAND(string asis) REJECT(integer)
    capture quietly `command'
    local rc = _rc
    local rejected = (`rc' != 0 | $TESTS_FAILED != 0)
    if `rejected' != `reject' {
        noisily di as error "Incorrect comparator verdict: `command' (rc=`rc', failures=$TESTS_FAILED)"
        exit 9
    }
    * Only these deliberately injected failures are consumed.
    reset_counters
    global AUDIT_GUARDS = $AUDIT_GUARDS + 1
end
global AUDIT_GUARDS = 0

* Byte comparisons must read through EOF, not accept an equal 10 MB prefix.
tempfile prefix_a prefix_b
mata:
chunk="x"
for (i=1;i<=20;i++) chunk=chunk+chunk
fa=fopen(st_local("prefix_a"),"w"); fb=fopen(st_local("prefix_b"),"w")
for (i=1;i<=10;i++) {
    fwrite(fa,chunk)
    fwrite(fb,chunk)
}
fwrite(fa,"a"); fwrite(fb,"b"); fclose(fa); fclose(fb)
assert(strlen(coverage_read_bytes(st_local("prefix_a")))==10*1048576+1)
end
audit_check, command(mata: assert(coverage_read_bytes("`prefix_a'")==coverage_read_bytes("`prefix_a'"))) reject(0)
audit_check, command(mata: assert(coverage_read_bytes("`prefix_a'")==coverage_read_bytes("`prefix_b'"))) reject(1)

audit_check, command(test_error_match, stata_cmd(error 198) ctools_cmd(error 198) testname("matching errors")) reject(0)
audit_check, command(test_error_match, stata_cmd(error 198) ctools_cmd(error 199) testname("different errors")) reject(1)
audit_check, command(test_error_match, stata_cmd(error 198) ctools_cmd(display 1) testname("error versus success")) reject(1)

* Extended missing codes are values: .a must not compare equal to .b or .
capture program drop cdestring
program cdestring
    syntax varlist [, GENerate(namelist) REPLACE]
    if "`replace'" != "" {
        destring `varlist', replace
        local output `varlist'
    }
    else {
        destring `varlist', generate(`generate')
        local output `generate'
    }
    if "$AUDIT_MUTATION" == "missing" {
        foreach v of local output {
            replace `v' = .b if `v' == .a
        }
    }
end
clear
set obs 3
gen str2 x = cond(_n == 1, ".a", "1")
foreach mutation in correct missing {
    global AUDIT_MUTATION `mutation'
    local bad = ("`mutation'" == "missing")
    audit_check, command(benchmark_destring x, generate(n) testname("generate missing code")) reject(`bad')
    audit_check, command(benchmark_destring x, replace testname("replace missing code")) reject(`bad')
}
program drop cdestring

capture program drop cencode
program cencode
    syntax varname, GENerate(name)
    if "$AUDIT_MUTATION" == "source" replace `varlist' = "z" in 2
    encode `varlist', generate(`generate')
    if "$AUDIT_MUTATION" == "missing" replace `generate' = .a if missing(`generate')
end
clear
set obs 3
gen str1 x = cond(_n == 1, "", "a")
foreach mutation in correct missing source {
    global AUDIT_MUTATION `mutation'
    local bad = ("`mutation'" != "correct")
    audit_check, command(benchmark_encode x, testname("encode missing code")) reject(`bad')
}
program drop cencode

* cf _all compares only variables in memory. Dropping a non-key output
* column must fail even if all remaining merge values/statuses are correct.
capture program drop cmerge
program cmerge
    merge `0'
    if "$AUDIT_MUTATION" == "dropvar" drop payload
end
clear
set obs 3
gen id = _n
gen payload = 10 * id
tempfile using_data
save "`using_data'"
drop payload
foreach mutation in correct dropvar {
    global AUDIT_MUTATION `mutation'
    local bad = ("`mutation'" == "dropvar")
    audit_check, command(benchmark_merge 1:1 id using "`using_data'", testname("merge dropped column")) reject(`bad')
}
program drop cmerge

tempfile ordered changed
save "`ordered'"
save "`changed'"
audit_check, command(assert_data_equal "`ordered'" "`changed'" "dataset identity") reject(0)
gsort -id
save "`changed'", replace
audit_check, command(assert_data_equal "`ordered'" "`changed'" "wrong row order") reject(1)
sort id
gen extra = 1
save "`changed'", replace
audit_check, command(assert_data_equal "`ordered'" "`changed'" "extra using variable") reject(1)

* Scalar dimension checks must still run when there are no observations.
clear
gen id = .
save "`ordered'", replace
save "`changed'", replace
audit_check, command(assert_data_equal "`ordered'" "`changed'" "empty dataset identity") reject(0)
gen extra = .
save "`changed'", replace
audit_check, command(assert_data_equal "`ordered'" "`changed'" "empty dataset extra variable") reject(1)

* Tiny score changes below 1e-7 and isolated missing scores used to bypass
* the propensity-score check. Incorrect treatment weights were never read.
capture program drop cpsmatch
program cpsmatch, rclass
    psmatch2 `0'
    if "$AUDIT_MUTATION" == "tiny" replace _pscore = 2 * _pscore in 1
    if "$AUDIT_MUTATION" == "missing" replace _pscore = . in 1
    if "$AUDIT_MUTATION" == "weight" replace _weight = 2 if _treated == 1
end
clear
set obs 8
gen byte treat = mod(_n, 2)
gen double p = _n * 1e-10
foreach mutation in correct tiny missing weight {
    global AUDIT_MUTATION `mutation'
    local bad = ("`mutation'" != "correct")
    audit_check, command(benchmark_psmatch2 treat, pscore(p) testname("matching mutation `mutation'")) reject(`bad')
}
program drop cpsmatch

* A no-op must not inherit the reference estimator's e() or the last graph.
program audit_noop, eclass
    exit
end
audit_noop
clear
set obs 8
gen x = _n
gen y = x + mod(_n, 3)
audit_check, command(coverage_estimate, reference(regress y x) candidate(regress y x) name("estimate identity")) reject(0)
audit_check, command(coverage_estimate, reference(regress y x) candidate(audit_noop) name("stale estimate")) reject(1)
audit_check, command(coverage_graph, reference(scatter y x) candidate(scatter y x) name("graph identity")) reject(0)
audit_check, command(coverage_graph, reference(scatter y x) candidate(audit_noop) name("stale graph")) reject(1)
audit_check, command(coverage_graph, reference(scatter y x) candidate(scatter x y) name("wrong graph coordinates")) reject(1)
program drop audit_noop

* Exercise the main e()-comparison helpers as well as coverage_estimate.
program audit_estimator, eclass
    if "$AUDIT_MUTATION" == "noop" exit
    $AUDIT_REFERENCE `0'
end
capture program drop cqreg
program cqreg, eclass
    audit_estimator `0'
end
capture program drop creghdfe
program creghdfe, eclass
    audit_estimator `0'
end
capture program drop civreghdfe
program civreghdfe, eclass
    audit_estimator `0'
end
capture program drop cppmlhdfe
program cppmlhdfe, eclass
    audit_estimator `0'
end
foreach pair in "qreg cqreg" "reghdfe creghdfe" "ivreghdfe civreghdfe" "ppmlhdfe cppmlhdfe" {
    gettoken ref candidate : pair
    global AUDIT_REFERENCE `ref'
    sysuse auto, clear
    foreach mutation in correct noop {
        global AUDIT_MUTATION `mutation'
        local bad = ("`mutation'" == "noop")
        if "`ref'" == "qreg" audit_check, command(benchmark_qreg price mpg weight, testname("qreg `mutation'")) reject(`bad')
        if "`ref'" == "reghdfe" audit_check, command(benchmark_reghdfe price mpg weight, absorb(foreign) testname("reghdfe `mutation'")) reject(`bad')
        if "`ref'" == "ivreghdfe" audit_check, command(benchmark_ivreghdfe price (mpg=weight), absorb(foreign) testname("IV `mutation'")) reject(`bad')
        if "`ref'" == "ppmlhdfe" audit_check, command(benchmark_ppmlhdfe price mpg weight, absorb(foreign) testname("PPML `mutation'")) reject(`bad')
    }
    program drop `candidate'
}
program drop audit_estimator
macro drop AUDIT_REFERENCE

* Summary/runner state must never describe skipped work as all tests passed.
test_skip "injected skip" "negative control"
capture print_summary "injected skip"
local rc = _rc
assert `rc' != 0
reset_counters
assert $TESTS_SKIPPED == 0

local guards = $AUDIT_GUARDS + 1
forvalues i = 1/`guards' {
    test_pass "comparator negative control `i'"
}
program drop audit_check
macro drop AUDIT_GUARDS AUDIT_MUTATION
print_summary "comparison_guards"
global CTOOLS_COMPONENT_COMPLETE "comparison_guards"
global CTOOLS_CASE_COMPLETE "comparison_guards"
exit
}
di as error "Unknown comparison_guards case: `__case'"
exit 198
