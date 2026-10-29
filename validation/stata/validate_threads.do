/* Thread settings: ctools, threads() session and saved settings, and how the
   dispatcher applies them. The developer's own session and saved settings are
   restored at the end. */
do validation/support/validate_setup.do
clear
set obs 2000
set seed 4412
gen double x = runiform()

* csort skips the plugin when x is already sorted, so reshuffle before each call.
capture program drop shuffle_x
program define shuffle_x
    quietly replace x = runiform()
end

local orig_session "$CTOOLS_THREADS"
global CTOOLS_THREADS
quietly ctools, environment_check
local orig_saved "`r(threads_saved)'"
local limit = r(threads_limit)
local cap = r(threads_cap)
local env_threads : env CTOOLS_NUM_THREADS

* Defaults: capped, never above the usable CPUs.
if `cap' == 16 & r(threads) >= 1 & r(threads) <= `limit' & r(cpus_available) <= r(cpus_detected) {
    test_pass "default is capped and within usable CPUs"
}
else test_fail "default thread policy" "threads=`r(threads)' limit=`limit' cap=`cap'"
if "`r(threads_source)'" != "session setting" & r(threads) == r(threads_default) {
    test_pass "no session setting after clearing"
}
else test_fail "cleared session setting" "source=`r(threads_source)'"

* Session setting applies to commands without threads().
quietly ctools, threads(2)
local expect = min(2, `limit')
if "$CTOOLS_THREADS" == "2" & r(threads) == `expect' & "`r(threads_source)'" == "session setting" {
    test_pass "session setting stored"
}
else test_fail "session setting" "global=$CTOOLS_THREADS threads=`r(threads)'"
shuffle_x
csort x
if r(threads_max) == `expect' test_pass "commands use the session setting"
else test_fail "session setting reaches commands" "threads_max=`r(threads_max)'"
shuffle_x
csort x, threads(1)
if r(threads_max) == 1 test_pass "threads() overrides the session setting"
else test_fail "threads() override" "threads_max=`r(threads_max)'"
shuffle_x
csort x, threads(256)
if r(threads_max) == `limit' test_pass "explicit request clamped to usable CPUs"
else test_fail "explicit clamp" "threads_max=`r(threads_max)' limit=`limit'"

* Invalid settings fail clearly and can be cleared.
global CTOOLS_THREADS abc
shuffle_x
capture csort x
if _rc == 198 test_pass "invalid session global rejected"
else test_fail "invalid session global" "rc=`=_rc'"
quietly ctools, environment_check
if "`r(threads_session)'" == "invalid" test_pass "invalid session global reported"
else test_fail "invalid session report" "session=`r(threads_session)'"
quietly ctools, threads(default)
if "$CTOOLS_THREADS" == "" & "`r(threads_source)'" != "session setting" test_pass "threads(default) clears the session"
else test_fail "threads(default)" "global=$CTOOLS_THREADS"
shuffle_x
csort x
test_pass "commands run after clearing"
foreach bad in 0 257 -1 1.5 abc "2 3" {
    capture ctools, threads(`bad')
    if _rc == 198 test_pass "rejected threads(`bad')"
    else test_fail "threads(`bad')" "rc=`=_rc'"
}
capture ctools, permanently
if _rc == 198 test_pass "permanently requires threads()"
else test_fail "permanently alone" "rc=`=_rc'"

* Saved setting: written outside the ado tree, applies in this session once
* the session setting is cleared, and is removed by threads(default).
if "`orig_saved'" == "invalid" {
    test_skip "saved setting" "existing settings file has an invalid threads= value"
}
else {
    quietly ctools, threads(3) permanently
    local file `"`r(threads_file)'"'
    capture confirm file `"`file'"'
    if !_rc & "`r(threads_saved)'" == "3" & "$CTOOLS_THREADS" == "3" test_pass "saved setting written"
    else test_fail "saved setting written" "saved=`r(threads_saved)' file=`file'"
    quietly ctools, threads(default)
    if "`env_threads'" == "" {
        if r(threads) == min(3, `limit') & "`r(threads_source)'" == "saved setting" {
            test_pass "saved setting is the default"
        }
        else test_fail "saved default" "threads=`r(threads)' source=`r(threads_source)'"
    }
    quietly ctools, threads(default) permanently
    if "`r(threads_saved)'" == "" & "$CTOOLS_THREADS" == "" test_pass "saved setting removed"
    else test_fail "saved setting removed" "saved=`r(threads_saved)'"
    if "`orig_saved'" != "" quietly ctools, threads(`orig_saved') permanently
}

global CTOOLS_THREADS `orig_session'
print_summary "threads"
global CTOOLS_COMPONENT_COMPLETE "threads"
