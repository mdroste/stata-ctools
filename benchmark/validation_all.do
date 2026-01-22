*! Master validation suite runner for all ctools commands

clear all
set more off
capture log close _all
adopath ++ "build"

* Configuration
local stop_on_failure = 0
local commands "csort cmerge cimport cexport creghdfe cqreg cbinscatter civreghdfe cencode"

* Initialize
log using "benchmark/validation_all.log", replace text

di ""
di as result "==============================================================================="
di as result "                    CTOOLS COMPREHENSIVE VALIDATION SUITE"
di as result "==============================================================================="
di ""
di as text "Date:    " as result c(current_date) " " c(current_time)
di as text "System:  " as result c(os) " " c(machine_type)
di as text "Stata:   " as result c(stata_version)
di ""

local start_time = clock(c(current_time), "hms")
local total_commands = 0
local failed_commands = 0
local failed_list ""

* Helper to count failures from log
capture program drop count_failures
program define count_failures, rclass
    args logfile
    local failures = 0
    local total = 0
    tempname fh
    capture file open `fh' using "`logfile'", read text
    if _rc != 0 {
        return scalar failures = -1
        return scalar total = 0
        exit
    }
    file read `fh' line
    while r(eof) == 0 {
        if strpos(`"`line'"', "[FAIL]") > 0 local ++failures
        if strpos(`"`line'"', "[PASS]") > 0 | strpos(`"`line'"', "[FAIL]") > 0 local ++total
        file read `fh' line
    }
    file close `fh'
    return scalar failures = `failures'
    return scalar total = `total'
end

* Run validation for each command
foreach cmd of local commands {
    local ++total_commands
    di ""
    di as text "{hline 70}"
    di as result "Running `cmd' validation..."
    di as text "{hline 70}"

    capture noisily do "benchmark/validation_`cmd'.do"
    local rc = _rc

    if `rc' != 0 {
        di as error "  `cmd' validation CRASHED with error `rc'"
        local ++failed_commands
        local failed_list "`failed_list' `cmd'"
        if `stop_on_failure' {
            di as error "Stopping due to failure."
            exit `rc'
        }
    }
    else {
        count_failures "benchmark/validation_`cmd'.log"
        if r(failures) > 0 {
            di as error "  `cmd': " r(failures) " of " r(total) " tests FAILED"
            local ++failed_commands
            local failed_list "`failed_list' `cmd'"
        }
        else {
            di as result "  `cmd': " r(total) " tests PASSED"
        }
    }
}

* Summary
local end_time = clock(c(current_time), "hms")
local elapsed = (`end_time' - `start_time') / 1000

di ""
di as result "==============================================================================="
di as result "                         VALIDATION SUITE SUMMARY"
di as result "==============================================================================="
di ""
di as text "Commands tested:   " as result `total_commands'
di as text "Commands passed:   " as result `=`total_commands' - `failed_commands''
di as text "Commands failed:   " as result `failed_commands'
di ""

if `failed_commands' == 0 {
    di as result "=============================================="
    di as result "   ALL VALIDATION TESTS PASSED SUCCESSFULLY"
    di as result "=============================================="
}
else {
    di as error "=============================================="
    di as error "   `failed_commands' COMMAND(S) HAD FAILURES"
    di as error "=============================================="
    di ""
    di as error "Failed commands:`failed_list'"
    di ""
    di as text "Check individual log files for details:"
    foreach cmd in `failed_list' {
        di as text "  - benchmark/validation_`cmd'.log"
    }
}

di ""
di as text "Total elapsed time: " as result %6.1f `elapsed' " seconds"
di ""

log close

if `failed_commands' > 0 exit 1
