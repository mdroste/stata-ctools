*! version 1.1.0 20261016 github.com/mdroste/stata-ctools

program define ctools, rclass
    version 14.1

    syntax [, Version ENVironment_check Verbose UPDATE THReads(string) PERManently]

    if "`permanently'" != "" & `"`threads'"' == "" {
        di as error "option permanently requires threads()"
        exit 198
    }

    if "`update'" != "" {
        di as text ""
        di as text "Checking for ctools updates..."
        di as text ""

        capture noisily net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build") replace

        if _rc == 0 {
            di as text ""
            di as text "{bf:ctools has been updated successfully.}"
            di as text "You may need to restart Stata for all changes to take effect."
        }
        else if _rc == 631 {
            di as text ""
            di as text "{bf:ctools is already up to date.}"
        }
        else {
            di as error ""
            di as error "Update failed with error code `=_rc'."
            di as error "Please check your internet connection and try again."
        }
        exit
    }

    local __ado_version "1.1.0"
    if "`environment_check'" != "" | "`version'" != "" | `"`threads'"' != "" {
        * Load the plugin if not already loaded
        _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

        if `"`threads'"' != "" {
            * threads(#) sets the session setting; permanently also saves it
            * in a per-user file that later Stata sessions read.
            local __request = strtrim(`"`threads'"')
            if lower(`"`__request'"') == "default" local __n 0
            else {
                capture confirm integer number `__request'
                if _rc | `: word count `__request'' != 1 {
                    di as error "threads() must be an integer from 1 through 256, or default"
                    exit 198
                }
                local __n = real("`__request'")
                if !inrange(`__n', 1, 256) {
                    di as error "threads() must be an integer from 1 through 256, or default"
                    exit 198
                }
            }
            if "`permanently'" != "" {
                capture noisily plugin call ctools_plugin, "threads_save `__n'"
                if _rc exit _rc
            }
            if `__n' > 0 global CTOOLS_THREADS `__n'
            else global CTOOLS_THREADS
        }

        if "`environment_check'" != "" | "`version'" != "" {
            capture noisily plugin call ctools_plugin, "version"
            local plugin_rc = _rc
            if `plugin_rc' == 0 {
                di as text "ctools ado `__ado_version'; plugin `ctools_plugin_version' (`ctools_plugin_revision')"
                return local version "`__ado_version'"
                return local plugin_version "`ctools_plugin_version'"
                return local plugin_revision "`ctools_plugin_revision'"
            }
            else {
                di as error "ctools plugin identity could not be read."
                exit `plugin_rc'
            }
        }
        if "`environment_check'" == "" & `"`threads'"' == "" exit

        capture noisily plugin call ctools_plugin, "threads_info"
        if _rc exit _rc
        local __source "`ctools_threads_source'"
        if "`ctools_threads_source'" == "automatic" {
            local __source "automatic: usable CPUs, at most `ctools_threads_cap'"
        }
        if "`ctools_threads_session'" != "" & "`ctools_threads_session'" != "invalid" {
            local __source "session setting"
        }
        foreach __item in session saved env omp omp_limit {
            local __show_`__item' "`ctools_threads_`__item''"
            if "`__show_`__item''" == "" local __show_`__item' "(not set)"
        }
        local __openmp = cond("`ctools_threads_openmp'" == "1", "enabled", "not available")
        di as text ""
        di as text "{ul:Threads}"
        di as text "  Commands without threads() use" _col(40) as result "`ctools_threads_effective'" as text "  (`__source')"
        di as text "  CPUs detected / usable" _col(40) as result "`ctools_threads_detected' / `ctools_threads_available'"
        di as text "  Session (global CTOOLS_THREADS)" _col(40) as result "`__show_session'"
        di as text "  Saved setting" _col(40) as result "`__show_saved'"
        if `"`ctools_threads_file'"' != "" {
            di as text "  Settings file" _col(40) as result `"`ctools_threads_file'"'
        }
        di as text "  CTOOLS_NUM_THREADS" _col(40) as result "`__show_env'"
        di as text "  OMP_NUM_THREADS" _col(40) as result "`__show_omp'"
        di as text "  OMP_THREAD_LIMIT" _col(40) as result "`__show_omp_limit'"
        di as text "  OpenMP" _col(40) as result "`__openmp'"
        if "`ctools_threads_session'" == "invalid" {
            di as error "  global CTOOLS_THREADS is invalid; ctools commands will fail until"
            di as error "  you run {cmd:ctools, threads(default)} or set it to an integer from 1 through 256"
        }
        if `"`threads'"' != "" {
            if `__n' > `ctools_threads_limit' {
                di as text "  Note: only `ctools_threads_limit' CPUs are usable here, so ctools uses at most `ctools_threads_limit' threads."
            }
            if "`permanently'" != "" & "`ctools_threads_env'" != "" {
                di as text "  Note: CTOOLS_NUM_THREADS takes precedence over the saved setting in new sessions."
            }
        }

        return scalar threads = `ctools_threads_effective'
        return scalar threads_default = `ctools_threads_default'
        return scalar threads_cap = `ctools_threads_cap'
        return scalar threads_limit = `ctools_threads_limit'
        return scalar cpus_detected = `ctools_threads_detected'
        return scalar cpus_available = `ctools_threads_available'
        return local threads_source "`__source'"
        return local threads_session "`ctools_threads_session'"
        return local threads_saved "`ctools_threads_saved'"
        return local threads_file `"`ctools_threads_file'"'
        exit
    }

    di as text ""
    di as result " ██████╗████████╗ ██████╗  ██████╗ ██╗     ███████╗"
    di as result "██╔════╝╚══██╔══╝██╔═══██╗██╔═══██╗██║     ██╔════╝"
    di as result "██║        ██║   ██║   ██║██║   ██║██║     ███████╗"
    di as result "██║        ██║   ██║   ██║██║   ██║██║     ╚════██║"
    di as result "╚██████╗   ██║   ╚██████╔╝╚██████╔╝███████╗███████║"
    di as result " ╚═════╝   ╚═╝    ╚═════╝  ╚═════╝ ╚══════╝╚══════╝"
    di as text ""
    di as text "  {it:C-accelerated tools for Stata}{col 56}v1.1.0"
    di as text "{hline 60}"
    di as text ""
    di as text "{ul:Data Management}"
    di as text "  {help cimport:cimport}      Import delimited/Excel data"
    di as text "  {help cexport:cexport}      Export delimited/Excel data"
    di as text "  {help csort:csort}        Sort dataset"
    di as text "  {help cmerge:cmerge}       Merge datasets"
    di as text "  {help csample:csample}      Sample without replacement"
    di as text "  {help cbsample:cbsample}     Sample with replacement"
    di as text "  {help cencode:cencode}      String to labeled numeric"
    di as text "  {help cdecode:cdecode}      Labeled numeric to string"
    di as text "  {help cdestring:cdestring}    String to numeric"
    di as text "  {help cwinsor:cwinsor}      Winsorize variables"
    di as text "  {help crangestat:crangestat}   Range statistics"
    di as text ""
    di as text "  {help cipolate:cipolate}     Linear interpolation"
    di as text "  {help csplit:csplit}       Split strings into fields"
    di as text "  {help crangejoin:crangejoin}   Interval joins"
    di as text ""
    di as text "{ul:Estimation}"
    di as text "  {help creghdfe:creghdfe}     OLS with multi-way FE"
    di as text "  {help civreghdfe:civreghdfe}   2SLS/GMM with multi-way FE"
    di as text "  {help cqreg:cqreg}        Quantile regression"
    di as text "  {help cppmlhdfe:cppmlhdfe}    PPML with multi-way FE"
    di as text "  {help cpsmatch:cpsmatch}     Propensity score matching"
    di as text ""
    di as text "{ul:Visualization}"
    di as text "  {help cbinscatter:cbinscatter}  Binned scatter plots"
    di as text ""
    di as text "{ul:Utilities}"
    di as text "  {cmd:ctools, environment_check} - Check plugin status and threads"
    di as text "  {cmd:ctools, threads(#)}        - Set threads ({cmd:permanently} saves it)"
    di as text "  {cmd:ctools, update}            - Update to latest version"
    di as text ""
    di as text "For help: {cmd:help} {it:command}    Source: {browse github.com/mdroste/stata-ctools}"
    di as text "{hline 60}"
end
