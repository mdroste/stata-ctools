* Run the README speedup comparisons and report the timings.
*
* Run in Stata from the ctools folder (the one holding build/ and README.md):
*
*     do validation/performance/readme_benchmarks.do
*     do validation/performance/readme_benchmarks.do quick
*     do validation/performance/readme_benchmarks.do pairs(sort merge) sizes(1000000)
*
* Options follow the file name without a comma (do itself rejects a comma).
*     sizes(numbers)  observation counts; default 1000000 10000000 25000000, the
*                     README sizes. Each must be at least 1000 and a multiple of 200.
*     quick           smoke test with 10,000 observations.
*     pairs(names)    comparisons to run, named by the Stata command; default all:
*                     export import sort merge sample bsample encode decode
*                     destring distinct ipolate cipolate split rangejoin winsor
*                     rangestat binscatter reghdfe ivreghdfe ppmlhdfe hdfe
*                     psmatch2 qreg
*     processors(#)   Stata/MP processors for the Stata commands; default 4 (or
*                     all available, if fewer), as in the published run. ctools
*                     commands use their default thread count.
*     output(folder)  where results go; default temp/readme_benchmarks.
*
* The timed workloads are the profile_readme block of profiles.do, which
* produced the published numbers; this file adds checks and a report.
*
* Run time: a comparison times only the command itself. Making the test data
* for each size, reloading it before each command and checking each result are
* untimed steps that also take time. During the run, a PROFILE_STEP line after
* each step gives its seconds and the running total. The report at the end gives
* the seconds of every step, and the run time of each comparison, slowest first.
* Press Break to stop a run early; the report then covers the finished steps.
*
* Output, in the output folder:
*     report.txt     the tables printed at the end
*     speedups.csv   Stata and ctools seconds and the speedup, per command and size
*     timings.csv    timed seconds, one row per command, size and method
*     steps.csv      seconds of every step, timed or not
*     profile.log    Stata output from the run
* The speedup table is also left in memory.
*
* The full run needs several GB of memory and tens of GB of free disk space and
* can take an hour or more. A comparison whose reference command (an SSC
* package) is not installed is skipped and listed, with the ssc install
* commands that add it. The published figures ran each comparison in its own
* Stata session on an otherwise idle machine; see docs/benchmarks/readme_profile/.

* Messages print with noisily; quietly blocks keep the code out of the output.
quietly {
    local allpairs "export import sort merge sample bsample encode decode destring distinct ipolate cipolate split rangejoin winsor rangestat binscatter reghdfe ivreghdfe ppmlhdfe hdfe psmatch2 qreg"

    capture confirm file "validation/performance/profiles.do"
    if _rc {
        noisily di as error "Run this do-file from the ctools folder, for example:"
        noisily di as error `"    cd "path/to/stata-ctools""'
        noisily di as error "    do validation/performance/readme_benchmarks.do"
        exit 601
    }

    local 0 `", `0'"'
    syntax [, Sizes(string) QUICK PAIRS(string) PROCessors(integer 0) Output(string)]

    * Sizes
    if "`quick'" != "" {
        if `"`sizes'"' != "" {
            noisily di as error "specify sizes() or quick, not both"
            exit 198
        }
        local sizes 10000
    }
    if `"`sizes'"' == "" local sizes "1000000 10000000 25000000"
    local normalized
    local largest 0
    foreach N of local sizes {
        capture confirm integer number `N'
        local bad = _rc != 0
        if !`bad' local bad = `N' < 1000 | mod(`N', 200) != 0
        if `bad' {
            noisily di as error `"sizes(): `N' is not a whole number of at least 1000 that is a multiple of 200"'
            exit 198
        }
        local N = `N'
        local normalized `normalized' `N'
        local largest = max(`largest', `N')
    }
    local sizes : list uniq normalized

    * Comparisons
    if `"`pairs'"' == "" local pairs "`allpairs'"
    local unknown : list pairs - allpairs
    if `"`unknown'"' != "" {
        noisily di as error `"pairs(): unknown comparison `unknown'"'
        noisily di as error "choose from: `allpairs'"
        exit 198
    }
    local pairs : list uniq pairs

    * Processors
    local original_processors = c(processors)
    if `processors' == 0 local processors = min(4, c(processors_max))
    if !inrange(`processors', 1, c(processors_max)) {
        noisily di as error "processors() must be from 1 through `c(processors_max)'"
        exit 198
    }

    * Output folder
    if `"`output'"' == "" {
        capture mkdir "temp"
        local output "temp/readme_benchmarks"
    }
    capture mkdir `"`output'"'
    mata: st_local("exists", strofreal(direxists(st_local("output"))))
    if !`exists' {
        noisily di as error `"output(): cannot create folder `output'"'
        exit 693
    }

    * ctools must be found and its plugin must load.
    capture confirm file "build/ctools.ado"
    if !_rc adopath ++ "build"
    capture which ctools
    if _rc {
        noisily di as error "ctools was not found. Build it into build/ or install it, then rerun."
        exit 111
    }
    capture noisily ctools, version
    if _rc {
        local rc = _rc
        noisily di as error "The ctools plugin did not load, so the ctools commands cannot run."
        exit `rc'
    }

    * Reference commands from SSC, and the packages they call.
    local need_distinct   distinct
    local need_cipolate   cipolate
    local need_rangejoin  rangejoin rangestat
    local need_winsor     gstats
    local need_rangestat  rangestat
    local need_psmatch2   psmatch2
    local need_binscatter binscatter
    local need_reghdfe    reghdfe ftools
    local need_ivreghdfe  ivreghdfe ivreg2 ranktest reghdfe ftools
    local need_ppmlhdfe   ppmlhdfe reghdfe ftools
    local need_hdfe       hdfe
    local run
    local skipped
    local missing
    foreach pair of local pairs {
        local absent
        foreach ado of local need_`pair' {
            capture which `ado'
            if _rc local absent `absent' `ado'
        }
        if "`absent'" == "" local run `run' `pair'
        else {
            local skipped `skipped' `pair'
            local missing : list missing | absent
        }
    }
    * Older ctools releases installed a linear cipolate, which is no reference.
    local stale_cipolate
    if strpos(" `run' ", " cipolate ") {
        findfile cipolate.ado
        local stale_cipolate `"`r(fn)'"'
        tempname header
        file open `header' using `"`stale_cipolate'"', read text
        file read `header' line
        file close `header'
        if strpos(`"`macval(line)'"', "stata-ctools") {
            local drop cipolate
            local run : list run - drop
            local skipped `skipped' cipolate
        }
        else local stale_cipolate
    }
    if "`skipped'" != "" {
        noisily di as error _n "Skipping comparisons whose reference command is not available: `skipped'"
        if "`missing'" != "" {
            noisily di as text "To include them, install the missing packages and rerun:"
            foreach ado of local missing {
                local package = cond("`ado'" == "gstats", "gtools", "`ado'")
                noisily di as text "    ssc install `package'"
            }
        }
        if `"`stale_cipolate'"' != "" {
            noisily di as text `"cipolate: `stale_cipolate'"'
            noisily di as text "is an old ctools copy that does linear interpolation, not SSC cipolate."
            noisily di as text "Delete that file, then run: ssc install cipolate"
        }
    }
    if "`run'" == "" {
        noisily di as error "None of the selected comparisons can run."
        exit 111
    }

    noisily di as text _n "{hline 78}"
    noisily di as text "README benchmarks"
    noisily di as text "  Observations:  `sizes'"
    noisily di as text "  Comparisons:   `run'"
    noisily di as text "  Stata/MP processors for Stata commands: `processors'; ctools uses its default threads"
    noisily di as text `"  Output folder: `output'"'
    noisily di as text "  A PROFILE_STEP line after each step gives its seconds. Press Break to stop."
    if `largest' >= 10000000 {
        noisily di as text "  This needs several GB of memory and tens of GB of free disk space, and"
        noisily di as text "  can take an hour or more."
    }
    noisily di as text "{hline 78}" _n
}

* Timed run. profiles.do starts with clear all and leaves processors and
* graphics changed, so restore those afterward.
local original_graphics = c(graphics)
capture log close profile
capture noisily do "validation/performance/profiles.do" profile_readme `"`output'"' `"`sizes'"' `processors' `"`run'"'
local rc = _rc
set graphics `original_graphics'
capture set processors `original_processors'
if `rc' {
    * A stopped run can leave the profile's log and result files open.
    capture log close profile
    capture file close _all
    if `rc' == 1 noisily di as error _n "The run was stopped."
    else noisily di as error _n `"The run stopped with error `rc'; see `output'/profile.log."'
    noisily di as error "The report covers the steps that finished."
}

* Report
quietly {
    capture log close readme_report
    log using `"`output'/report.txt"', text replace name(readme_report)
    local count : word count `sizes'

    * Speedups, from timings.csv
    local have_speedups 0
    capture import delimited using `"`output'/timings.csv"', clear varnames(1) ///
        stringcols(1 3) numericcols(2 4 5 6) asdouble
    if !_rc & _N > 0 local have_speedups 1
    if `have_speedups' {
        generate long order = _n
        bysort command (order): replace order = order[1]
        reshape wide seconds rc check_rc, i(command observations) j(method) string
        foreach v in secondsref secondsc rcref rcc check_rcref check_rcc {
            capture confirm variable `v'
            if _rc generate double `v' = .
        }
        generate byte ok = rcref == 0 & rcc == 0 & check_rcref == 0 & check_rcc == 0 & ///
            !missing(secondsref, secondsc) & secondsref > 0 & secondsc > 0
        generate double speedup = secondsref / secondsc if ok
        generate str15 stata_command = command
        generate str19 ctools_command = "c" + command
        replace stata_command = "gstats winsor" if command == "winsor"
        replace ctools_command = "cinterpolate" if command == "ipolate"
        replace ctools_command = "cinterpolate, cubic" if command == "cipolate"
        replace ctools_command = "cpsmatch" if command == "psmatch2"
        sort order observations

        noisily di as text _n "Timed seconds; speedup = Stata time / ctools time"
        noisily di as text "{hline 92}"
        noisily di as text %-15s "Stata command" %-21s "ctools command" %14s "Observations" ///
            %12s "Stata (s)" %12s "ctools (s)" %18s "Speedup"
        noisily di as text "{hline 92}"
        forvalues i = 1/`=_N' {
            if ok[`i'] local result = strtrim(string(speedup[`i'], "%9.2f")) + "x"
            else if missing(rcref[`i']) | missing(rcc[`i']) local result "not run"
            else if rcref[`i'] == 1 | rcc[`i'] == 1 local result "stopped"
            else if rcref[`i'] local result "Stata r(`=rcref[`i']')"
            else if rcc[`i'] local result "ctools r(`=rcc[`i']')"
            else if check_rcref[`i'] | check_rcc[`i'] local result "check failed"
            else local result "too fast to time"
            noisily di as text %-15s stata_command[`i'] %-21s ctools_command[`i'] ///
                as result %14.0fc observations[`i'] ///
                %12.3f secondsref[`i'] %12.3f secondsc[`i'] %18s "`result'"
        }
        noisily di as text "{hline 92}"

        * Ranges are rounded outward as in the README: tenths below 2x, whole
        * numbers below 10x, then tens. A range needs every size.
        noisily di as text _n "Approximate speedup, as in the README table (observations: `sizes')"
        noisily di as text "{hline 54}"
        noisily di as text %-15s "Stata command" %-21s "ctools command" %18s "Approx. speedup"
        noisily di as text "{hline 54}"
        foreach pair of local run {
            count if command == "`pair'"
            if r(N) == 0 continue
            summarize speedup if command == "`pair'"
            if r(N) == `count' {
                local low = r(min)
                local high = r(max)
                if `low' < 2 local low = strtrim(string(floor(`low' * 10) / 10, "%9.1f"))
                else if `low' < 10 local low = floor(`low')
                else local low = floor(`low' / 10) * 10
                if `high' < 2 local high = strtrim(string(ceil(`high' * 10) / 10, "%9.1f"))
                else if `high' < 10 local high = ceil(`high')
                else local high = ceil(`high' / 10) * 10
                if "`low'" == "`high'" local range "`low'x"
                else local range "`low'-`high'x"
            }
            else local range "-"
            levelsof stata_command if command == "`pair'", local(stata_name) clean
            levelsof ctools_command if command == "`pair'", local(ctools_name) clean
            noisily di as text %-15s "`stata_name'" %-21s "`ctools_name'" as result %18s "`range'"
        }
        noisily di as text "{hline 54}"
        tempfile speedups
        save `speedups'
    }

    * Run time of every step, from steps.csv
    local have_steps 0
    capture import delimited using `"`output'/steps.csv"', clear varnames(1) ///
        stringcols(2 3 4) numericcols(1 5) asdouble
    if !_rc & _N > 0 local have_steps 1
    if `have_steps' {
        generate long order = _n
        bysort observations command method (order): replace order = order[1]
        reshape wide seconds, i(observations command method) j(step) string
        foreach v in secondssetup secondsload secondsrun secondscheck {
            capture confirm variable `v'
            if _rc generate double `v' = .
        }
        egen double total = rowtotal(secondssetup secondsload secondsrun secondscheck)
        generate str19 name = command
        replace name = "c" + command if method == "c"
        replace name = "gstats winsor" if command == "winsor" & method == "ref"
        replace name = "cinterpolate" if command == "ipolate" & method == "c"
        replace name = "cinterpolate, cubic" if command == "cipolate" & method == "c"
        replace name = "cpsmatch" if command == "psmatch2" & method == "c"
        replace name = "make test data" if command == "setup"
        sort order
        summarize total, meanonly
        local grand = r(sum)
        generate double share = 100 * total / `grand'

        noisily di as text _n "Run time of each step, in the order run (seconds)"
        noisily di as text "Command is the timed part; loading the data and checking the result are not timed."
        noisily di as text "{hline 92}"
        noisily di as text %-21s "Step" %14s "Observations" %11s "Load data" %11s "Command" ///
            %11s "Checks" %12s "Total" %12s "Share"
        noisily di as text "{hline 92}"
        forvalues i = 1/`=_N' {
            if command[`i'] == "setup" {
                noisily di as text %-21s name[`i'] as result %14.0fc observations[`i'] ///
                    %33s "" %12.3f total[`i'] %11.1f share[`i'] "%"
            }
            else {
                noisily di as text %-21s name[`i'] as result %14.0fc observations[`i'] ///
                    %11.3f secondsload[`i'] %11.3f secondsrun[`i'] %11.3f secondscheck[`i'] ///
                    %12.3f total[`i'] %11.1f share[`i'] "%"
            }
        }
        noisily di as text "{hline 92}"
        noisily di as text %-21s "Total" as result %59s "" %12.3f `grand'

        tempfile steps
        save `steps'
        generate double stata_total = total if method == "ref"
        generate double ctools_total = total if method == "c"
        collapse (sum) stata_total ctools_total total (min) order, by(command)
        generate str15 label = command
        replace label = "gstats winsor" if command == "winsor"
        replace label = "make test data" if command == "setup"
        generate double share = 100 * total / `grand'
        gsort -total order

        noisily di as text _n "Run time by comparison, slowest first (seconds, all sizes, untimed steps included)"
        noisily di as text "{hline 66}"
        noisily di as text %-21s "Comparison" %12s "Stata" %12s "ctools" %12s "Total" %9s "Share"
        noisily di as text "{hline 66}"
        forvalues i = 1/`=_N' {
            if command[`i'] == "setup" {
                noisily di as text %-21s label[`i'] as result %24s "" ///
                    %12.3f total[`i'] %8.1f share[`i'] "%"
            }
            else {
                noisily di as text %-21s label[`i'] as result %12.3f stata_total[`i'] ///
                    %12.3f ctools_total[`i'] %12.3f total[`i'] %8.1f share[`i'] "%"
            }
        }
        noisily di as text "{hline 66}"

        use `steps', clear
        collapse (sum) total, by(observations)
        noisily di as text _n "Run time by size (seconds)"
        forvalues i = 1/`=_N' {
            noisily di as text %14.0fc observations[`i'] " observations" as result %16.3f total[`i']
        }
        noisily di as text %-27s "Total" as result %16.3f `grand' as text "  (" ///
            strtrim(string(`grand' / 60, "%9.1f")) " minutes)"
    }

    if "`skipped'" != "" noisily di as error _n "Skipped (reference command not available): `skipped'"
    noisily di as text _n "Published figures, with the machine and settings behind them, are in"
    noisily di as text "docs/benchmarks/readme_profile/."
    log close readme_report

    if `have_speedups' {
        use `speedups', clear
        keep order stata_command ctools_command observations secondsref secondsc speedup
        rename (secondsref secondsc) (stata_seconds ctools_seconds)
        label variable stata_command "Stata command"
        label variable ctools_command "ctools command"
        label variable observations "Observations"
        label variable stata_seconds "Stata time (s)"
        label variable ctools_seconds "ctools time (s)"
        label variable speedup "Stata time / ctools time"
        export delimited stata_command ctools_command observations stata_seconds ///
            ctools_seconds speedup using `"`output'/speedups.csv"', replace
        drop order
    }
    else clear
    noisily di as text _n `"Files in `output': report.txt, speedups.csv, timings.csv, steps.csv, profile.log"'
}
if `rc' exit `rc'
